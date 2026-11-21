// Check 1 — genome-wide slippage propensity, run in a browser tab on plasmodb.org (raw-files
// need the site session). For each release-71 genome: CDS intervals from the GFF; the rest of
// each nuclear contig (apicoplast / mitochondrion skipped) is "non-coding". Metrics per Mb of
// non-coding sequence: A/T homopolymer runs >= 10/15/20 bp; perfect period-2 repeats made of
// A/T only (e.g. (AT)n) >= 8 / 12 units; perfect period-3 A/T-only repeats (AAT, ATT, ...)
// >= 5 / 8 units; the same on CDS. Null: 10 Mb simulated from a first-order Markov chain fitted
// to that genome's non-coding dinucleotides (keeps base composition and A/T clustering).
// Usage: genomeStats('Pfalciparum3D7'); results in window.GS[org] (nc, co, sim).
window.GS = window.GS || {};
const CODE = {65: 1, 67: 2, 71: 3, 84: 4, 97: 1, 99: 2, 103: 3, 116: 4};
function scan(seq, n, res) {
  let i = 0;
  while (i < n) { const c = seq[i]; if (c === 0) { i++; continue; } let j = i + 1; while (j < n && seq[j] === c) j++;
    if (c === 1 || c === 4) res.hp[Math.min(j - i, 60)]++; i = j; }
  for (const k of [2, 3]) { let r = 0;
    for (let p = k; p <= n; p++) { const ok = p < n && seq[p] !== 0 && seq[p] === seq[p - k];
      if (ok) r++; else { if (r > 0) { const len = r + k, st = p - len; let homo = true;
          for (let q = st + 1; q < st + k; q++) if (seq[q] !== seq[st]) homo = false;
          if (!homo) { const units = Math.floor(len / k); const key = k === 2 ? 'd' : 't'; let at = true;
            for (let q = st; q < st + k; q++) if (seq[q] === 2 || seq[q] === 3) at = false;
            const nm = key + (at ? 'AT' : 'GC');
            if (units >= 5) res[nm + '5'] = (res[nm + '5'] || 0) + 1;
            if (units >= 8) res[nm + '8'] = (res[nm + '8'] || 0) + 1;
            if (units >= 12) res[nm + '12'] = (res[nm + '12'] || 0) + 1; } } r = 0; } } }
}
function newRes() { return {hp: new Array(61).fill(0), L: 0, base: [0, 0, 0, 0, 0]}; }
function markovSim(di, base, N) {
  const seq = new Uint8Array(N), cum = [];
  for (let a = 1; a <= 4; a++) { let s = 0; const row = []; for (let b = 1; b <= 4; b++) { s += di[a][b]; row.push(s); } cum[a] = row.map(x => x / s); }
  const tot = base[1] + base[2] + base[3] + base[4]; let r = Math.random() * tot, c = 1;
  for (let b = 1; b <= 4; b++) { r -= base[b]; if (r <= 0) { c = b; break; } }
  seq[0] = c; for (let i = 1; i < N; i++) { const u = Math.random(); const row = cum[seq[i - 1]]; let b = 0; while (b < 3 && u > row[b]) b++; seq[i] = b + 1; }
  return seq;
}
window.genomeStats = async function (org, prefix) {
  const st = window.GS[org] = {status: 'fetch'};
  const base = '/a/service/raw-files/release-71/' + org + '/', pfx = prefix || 'PlasmoDB';
  const gff = await (await fetch(base + 'gff/data/' + pfx + '-71_' + org + '.gff')).text();
  const cds = {};
  for (const line of gff.split('\n')) { if (!line || line[0] === '#') continue; const f = line.split('\t');
    if (f[2] !== 'CDS') continue; (cds[f[0]] = cds[f[0]] || []).push([+f[3] - 1, +f[4]]); }
  const txt = await (await fetch(base + 'fasta/data/' + pfx + '-71_' + org + '_Genome.fasta')).text();
  const nc = newRes(), co = newRes(), di = [null, [0, 0, 0, 0, 0], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0]];
  let ncontig = 0;
  for (const c of txt.split('>')) { if (!c) continue; const nl = c.indexOf('\n'); const head = c.slice(0, nl);
    const id = head.split(/\s|\|/)[0]; if (/API|MIT|apicoplast|mito/i.test(head)) continue;
    const s = c.slice(nl + 1).replace(/\s/g, ''), n = s.length; ncontig++;
    const mask = new Uint8Array(n); for (const [a, b] of (cds[id] || [])) mask.fill(1, a, Math.min(b, n));
    const A = new Uint8Array(n), B = new Uint8Array(n);
    for (let i = 0; i < n; i++) { const x = CODE[s.charCodeAt(i)] || 0; if (mask[i]) B[i] = x; else A[i] = x; }
    for (let i = 0; i < n; i++) { const x = A[i]; if (x) { nc.base[x]++; nc.L++; if (i > 0 && A[i - 1]) di[A[i - 1]][x]++; }
      const y = B[i]; if (y) { co.base[y]++; co.L++; } }
    scan(A, n, nc); scan(B, n, co); }
  const SIMN = 10000000, sim = newRes(), sq = markovSim(di, nc.base, SIMN); sim.L = SIMN; scan(sq, SIMN, sim);
  Object.assign(st, {status: 'done', ncontig, nc, co, sim});
};
