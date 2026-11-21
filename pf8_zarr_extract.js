// Pf8 site-level extraction for poly-N tracts, run in a browser tab on
// https://pf8-release.cog.sanger.ac.uk (the zip is not reachable from scripted shells here).
// Reads Pf8.zarr.zip (1.98 TB) by HTTP Range requests: zip central directory -> per-chunk
// local headers -> blosc(lz4, byte shuffle) chunks. Only site arrays are touched
// (variants/POS, CHROM, FILTER_PASS, AN, AC, REF, ALT; ~330 MB in total), never calldata.
// Usage: paste into DevTools / javascript_tool; set window.REG = {chrom: [[start,end,tractId],...]}
// (sorted by start) before calling runScan(). Results in window.SCAN.recs:
// [tract, pos, pass, AN, REF, "ALT1,ALT2", "AC1,AC2"]. Run ONCE: a second concurrent run
// duplicates records (dedupe on tract|pos|REF|ALT if in doubt).

const ZIP = '/Pf8.zarr.zip';

async function zipIndex() {
  const size = +(await fetch(ZIP, {method: 'HEAD'})).headers.get('content-length');
  const tailN = 1 << 20;
  const b = new Uint8Array(await (await fetch(ZIP, {headers: {Range: `bytes=${size - tailN}-${size - 1}`}})).arrayBuffer());
  const dv = new DataView(b.buffer);
  let eocd = -1; for (let i = b.length - 22; i >= 0; i--) if (dv.getUint32(i, true) === 0x06054b50) { eocd = i; break; }
  let loc = -1; for (let i = eocd - 20; i >= Math.max(0, eocd - 100); i--) if (dv.getUint32(i, true) === 0x07064b50) { loc = i; break; }
  const z64 = Number(dv.getBigUint64(loc + 8, true));
  const zb = new DataView(await (await fetch(ZIP, {headers: {Range: `bytes=${z64}-${z64 + 56}`}})).arrayBuffer());
  const cdSize = Number(zb.getBigUint64(40, true)), cdOff = Number(zb.getBigUint64(48, true));
  const cd = new Uint8Array(await (await fetch(ZIP, {headers: {Range: `bytes=${cdOff}-${cdOff + cdSize - 1}`}})).arrayBuffer());
  const cv = new DataView(cd.buffer), td = new TextDecoder(), ents = new Map();
  let p = 0;
  while (p + 46 <= cd.length && cv.getUint32(p, true) === 0x02014b50) {
    const meth = cv.getUint16(p + 10, true); let csize = cv.getUint32(p + 20, true), usize = cv.getUint32(p + 24, true);
    const nl = cv.getUint16(p + 28, true), xl = cv.getUint16(p + 30, true), cl = cv.getUint16(p + 32, true);
    let off = cv.getUint32(p + 42, true); const name = td.decode(cd.subarray(p + 46, p + 46 + nl));
    let q = p + 46 + nl; const qe = q + xl;
    while (q < qe) { const id = cv.getUint16(q, true), sz = cv.getUint16(q + 2, true);
      if (id === 1) { let k = q + 4; if (usize === 0xFFFFFFFF) { usize = Number(cv.getBigUint64(k, true)); k += 8; }
        if (csize === 0xFFFFFFFF) { csize = Number(cv.getBigUint64(k, true)); k += 8; }
        if (off === 0xFFFFFFFF) { off = Number(cv.getBigUint64(k, true)); } }
      q += 4 + sz; }
    ents.set(name, {meth, csize, off}); p += 46 + nl + xl + cl;
  }
  return ents;
}

async function zread(E, name) {
  const e = E.get(name);
  const hb = new DataView(await (await fetch(ZIP, {headers: {Range: `bytes=${e.off}-${e.off + 29}`}})).arrayBuffer());
  const s = e.off + 30 + hb.getUint16(26, true) + hb.getUint16(28, true);
  let b = new Uint8Array(await (await fetch(ZIP, {headers: {Range: `bytes=${s}-${s + e.csize - 1}`}})).arrayBuffer());
  if (e.meth === 8) b = new Uint8Array(await new Response(new Blob([b]).stream().pipeThrough(new DecompressionStream('deflate-raw'))).arrayBuffer());
  return b;
}

function lz4(src, dst, dpos) {
  let s = 0; const n = src.length;
  while (s < n) { const tok = src[s++]; let lit = tok >> 4;
    if (lit === 15) { let b; do { b = src[s++]; lit += b; } while (b === 255); }
    dst.set(src.subarray(s, s + lit), dpos); s += lit; dpos += lit; if (s >= n) break;
    const off = src[s] | (src[s + 1] << 8); s += 2; let ml = tok & 15;
    if (ml === 15) { let b; do { b = src[s++]; ml += b; } while (b === 255); }
    ml += 4; let m = dpos - off; for (let i = 0; i < ml; i++) dst[dpos++] = dst[m++]; }
  return dpos;
}

function blosc(b) {
  const dv = new DataView(b.buffer, b.byteOffset, b.byteLength);
  const flags = b[2], ts = b[3], nbytes = dv.getUint32(4, true), bsize = dv.getUint32(8, true);
  const out = new Uint8Array(nbytes);
  if (flags & 2) { out.set(b.subarray(16, 16 + nbytes)); return out; }
  const nblocks = Math.ceil(nbytes / bsize), dontSplit = (flags & 0x10) >> 4, tmp = new Uint8Array(bsize);
  for (let i = 0; i < nblocks; i++) {
    const bstart = dv.getUint32(16 + 4 * i, true), leftover = (i === nblocks - 1) && (nbytes % bsize !== 0);
    const thisb = leftover ? nbytes - i * bsize : bsize, ns = (!dontSplit && !leftover && ts > 1) ? ts : 1, neb = thisb / ns;
    let p = bstart, dpos = 0; const target = (flags & 1) && ts > 1 ? tmp : out, base = target === tmp ? 0 : i * bsize;
    for (let j = 0; j < ns; j++) { const cb = dv.getUint32(p, true); p += 4;
      if (cb === neb) target.set(b.subarray(p, p + cb), base + dpos); else lz4(b.subarray(p, p + cb), target, base + dpos);
      p += cb; dpos += neb; }
    if (target === tmp) { const ne = thisb / ts, o = i * bsize;
      for (let e = 0; e < ne; e++) for (let k = 0; k < ts; k++) out[o + e * ts + k] = tmp[k * ne + e];
      for (let r = ne * ts; r < thisb; r++) out[o + r] = tmp[r]; }
  }
  return out;
}

function vlen(u8) {
  const dv = new DataView(u8.buffer, u8.byteOffset, u8.byteLength), td = new TextDecoder();
  const n = dv.getUint32(0, true), a = new Array(n); let p = 4;
  for (let i = 0; i < n; i++) { const l = dv.getUint32(p, true); p += 4; a[i] = td.decode(u8.subarray(p, p + l)); p += l; }
  return a;
}

async function runScan(PAD = 60) {
  const E = await zipIndex();
  const nChunks = Math.ceil(12493205 / 131072);   // 96 chunks of variants/*
  window.SCAN = {done: 0, recs: [], err: null, running: true};
  for (let ci = 0; ci < nChunks; ci++) {
    const pos = new Int32Array(blosc(await zread(E, 'variants/POS/' + ci)).buffer.slice(0));
    const chr = vlen(blosc(await zread(E, 'variants/CHROM/' + ci)));
    const hits = [];
    for (let i = 0; i < pos.length;) { const c = chr[i]; let k = i; while (k < pos.length && chr[k] === c) k++;
      const regs = window.REG[c];
      if (regs) { let r = 0; for (let v = i; v < k; v++) { const p = pos[v];
        while (r < regs.length && regs[r][1] < p) r++;
        for (let q = r; q < regs.length && regs[q][0] - PAD <= p; q++) if (p <= regs[q][1]) hits.push([v, regs[q][2]]); } }
      i = k; }
    if (hits.length) {
      const [fpb, anb, acb, refb, altb] = await Promise.all(['FILTER_PASS/' + ci, 'AN/' + ci, 'AC/' + ci + '.0', 'REF/' + ci, 'ALT/' + ci + '.0'].map(x => zread(E, 'variants/' + x)));
      const fp = blosc(fpb), an = new Int32Array(blosc(anb).buffer.slice(0)), ac = new Int32Array(blosc(acb).buffer.slice(0));
      const ref = vlen(blosc(refb)), alt = vlen(blosc(altb));
      for (const [v, ti] of hits) { const a = [], cnt = [];
        for (let x = 0; x < 6; x++) { const s = alt[v * 6 + x]; if (s) { a.push(s); cnt.push(ac[v * 6 + x]); } }
        window.SCAN.recs.push([ti, pos[v], fp[v], an[v], ref[v], a.join(','), cnt.join(',')]); }
    }
    window.SCAN.done = ci + 1;
  }
  window.SCAN.running = false;
}
// Per-tract metrics were then computed in the page exactly as in pop_length_variation.py
// (PASS sites overlapping [start-1, end]; '*' alleles ignored; length-changing alleles with AC>=1).
