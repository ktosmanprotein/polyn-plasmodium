#!/usr/bin/env python3
"""Figure supplements, source data files and TIFF exports for the eLife submission.
Usage: python3 make_supplements_sourcedata.py <results_dir> <submission_dir>
Writes <submission_dir>/figures (TIFF, 300 dpi, main + supplements) and <submission_dir>/source_data (CSV).
"""
import csv, json, os, sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image

R, S = sys.argv[1], sys.argv[2]
FIG = os.path.join(S, 'figures')
SD = os.path.join(S, 'source_data')
os.makedirs(FIG, exist_ok=True)
os.makedirs(SD, exist_ok=True)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
BLUE, ORANGE, AQUA, YELLOW, GREY, LGREY, INK, MUTED = '#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#a9a8a2', '#e9e8e4', '#1a1a19', '#5c5b57'
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 7.5, 'axes.spines.top': False, 'axes.spines.right': False,
                     'axes.edgecolor': MUTED, 'axes.linewidth': 0.6, 'xtick.color': MUTED, 'ytick.color': MUTED,
                     'legend.frameon': False, 'pdf.fonttype': 42})


def J(p):
    return json.load(open(os.path.join(R, p)))


def T(p):
    return list(csv.DictReader(open(os.path.join(R, p)), delimiter='\t'))


def wcsv(name, header, rows):
    with open(os.path.join(SD, name), 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(header)
        w.writerows(rows)


def save_tiff(fig, name):
    p = os.path.join('/tmp', name + '.png')
    fig.savefig(p, dpi=400, bbox_inches='tight', facecolor='white')
    fig.savefig(os.path.join(FIG, name + '.pdf'), bbox_inches='tight', facecolor='white')
    plt.close(fig)
    to_tiff(p, os.path.join(FIG, name + '.tif'))
    os.remove(p)


def to_tiff(src, dst):
    im = Image.open(src).convert('RGB')          # 8-bit RGB
    im.save(dst, compression='tiff_lzw', dpi=(300, 300))
    return im.size


def title(ax, s):
    ax.set_title(s, loc='left', fontsize=8, fontweight='bold')


SPN = {'pfalciparum': 'P. falciparum', 'ppraefalciparum': 'P. praefalciparum', 'preichenowi': 'P. reichenowi',
       'pbillcollinsi': 'P. billcollinsi', 'pblacklocki': 'P. blacklocki', 'padleri': 'P. adleri', 'pgaboni': 'P. gaboni',
       'pgallinaceum': 'P. gallinaceum', 'prelictum': 'P. relictum', 'pberghei': 'P. berghei', 'pyoelii': 'P. yoelii',
       'pvivax': 'P. vivax', 'pknowlesi': 'P. knowlesi'}
GROUP = {k: ('Laverania' if k in ('pfalciparum', 'ppraefalciparum', 'preichenowi', 'pbillcollinsi', 'pblacklocki', 'padleri', 'pgaboni')
             else 'avian' if k in ('pgallinaceum', 'prelictum') else 'rodent' if k in ('pberghei', 'pyoelii') else 'primate') for k in SPN}
GC = {'Laverania': BLUE, 'avian': ORANGE, 'rodent': AQUA, 'primate': YELLOW}

# ---------------- main figures -> TIFF ----------------
sizes = {}
for i in range(1, 7):
    src = os.path.join(R, 'figures_v2', f'Figure{i}.png')
    sizes[f'Figure{i}'] = to_tiff(src, os.path.join(FIG, f'Figure{i}.tif'))
    os.system(f"cp '{os.path.join(R, 'figures_v2', f'Figure{i}.pdf')}' '{os.path.join(FIG, f'Figure{i}.pdf')}'")

# ---------------- Figure 1 supplement 1: full repeat spectra ----------------
spec = {r['species']: r for r in T('step9/repeat_length_spectrum.tsv')}
cols = ['ATn_3u', 'ATn_4u', 'ATn_5u', 'ATn_6u', 'ATn_7u', 'ATn_8to11u', 'ATn_12plus']
colt = ['AATtype_2u', 'AATtype_3u', 'AATtype_4u', 'AATtype_5to7u', 'AATtype_8plus']
fig, axs = plt.subplots(1, 2, figsize=(7.2, 2.8))
MK = {}
cnt = {}
for k in spec:
    g = GROUP[k]
    MK[k] = 'osD^v'[cnt.get(g, 0)]
    cnt[g] = cnt.get(g, 0) + 1
for k, r in spec.items():
    g = GROUP[k]
    axs[0].plot(range(len(cols)), [float(r[c]) for c in cols], '-', marker=MK[k], ms=3.5, lw=1.3, color=GC[g], label=SPN[k], alpha=0.9)
    axs[1].plot(range(len(colt)), [float(r[c]) for c in colt], '-', marker=MK[k], ms=3.5, lw=1.3, color=GC[g], alpha=0.9)
for ax, lab, c in [(axs[0], ['3', '4', '5', '6', '7', '8–11', '12+'], '(AT)n repeats per Mb'), (axs[1], ['2', '3', '4', '5–7', '8+'], 'A/T-only triplet repeats per Mb')]:
    ax.set_yscale('log')
    ax.set_xticks(range(len(lab)))
    ax.set_xticklabels(lab)
    ax.set_xlabel('Repeat length (units)')
    ax.set_ylabel(c)
title(axs[0], 'A  (AT)n microsatellites')
title(axs[1], 'B  A/T-only triplet repeats')
axs[0].legend(fontsize=6, ncol=1, loc='center left', bbox_to_anchor=(1.0, 0.5), handlelength=1.6, prop={'style': 'italic', 'size': 6})
fig.tight_layout()
save_tiff(fig, 'Figure1-figure supplement 1')
wcsv('Figure 1-source data 1.csv',
     ['species', 'group', 'pct_proteins_win20_N60', 'pct_proteins_win30_N50', 'pct_proteins_Nrun10', 'proteome_N_pct', 'noncoding_Mb', 'noncoding_AT',
      'ATn_ge8_per_Mb', 'ATn_ge12_per_Mb', 'AATtype_ge5_per_Mb', 'homopolymer_ge15_per_Mb', 'markov_ATn_ge8_per_Mb'] + cols + colt,
     [[k, GROUP[k], p['pct_win20_N60'], p['pct_win30_N50'], p['pct_Nrun10'], p['proteome_N_pct'], g['nc_Mb'], g['nc_AT'], g['dAT8_perMb'],
       g['dAT12_perMb'], g['tAT5_perMb'], g['hp15_perMb'], g['markov_dAT8']] + [spec.get(k, {}).get(c, '') for c in cols + colt]
      for k in SPN for p in [{r['species']: r for r in T('step8/absolute_polyN_prevalence.tsv')}[k]]
      for g in [{r['species']: r for r in T('step9/genome_repeats.tsv')}[k]]])

# ---------------- Figure 2 supplement 1: CA-geometry vs DSSP ----------------
ins = J('step8/insertion_site_dssp.json')
fig, ax = plt.subplots(figsize=(3.6, 2.4))
cats = ['element', 'boundary', 'loop']
labs = ['inside helix/strand', 'element edge', 'loop']
for j, (meth, col) in enumerate([('ca_geometry', GREY), ('dssp', ORANGE)]):
    for i, c in enumerate(cats):
        o, n = 100 * ins[meth][c]['obs'], 100 * ins[meth][c]['null']
        y = i + (j - 0.5) * 0.3
        ax.plot([n, o], [y, y], color=LGREY, lw=3, zorder=1)
        ax.scatter([n], [y], s=18, color='white', edgecolor=col, zorder=3)
        ax.scatter([o], [y], s=22, color=col, zorder=3)
ax.set_yticks(range(3))
ax.set_yticklabels(labs)
ax.invert_yaxis()
ax.set_xlabel('Sites (%)')
from matplotlib.lines import Line2D
ax.legend(handles=[Line2D([], [], marker='o', ls='', color=GREY, label='Cα geometry'), Line2D([], [], marker='o', ls='', color=ORANGE, label='DSSP'),
                   Line2D([], [], marker='o', ls='', mfc='white', mec=MUTED, label='random sites'), Line2D([], [], marker='o', ls='', color=MUTED, label='insertion sites')],
          fontsize=6, loc='center', bbox_to_anchor=(0.52, 0.3))
title(ax, f"Two secondary-structure assignments (n = {ins['n_used']})")
fig.tight_layout()
save_tiff(fig, 'Figure2-figure supplement 1')
s5 = J('step5/step5_results.json')['A_location_by_class']
locs = sorted({k for d in s5.values() for k in d})
wcsv('Figure 2-source data 1.csv', ['lcr_class'] + locs, [[c] + [s5[c].get(l, 0) for l in locs] for c in s5])
wcsv('Figure 2-source data 2.csv', ['method', 'site_class', 'observed_fraction', 'null_fraction', 'z', 'n_insertions'],
     [[m, c, ins[m][c]['obs'], ins[m][c]['null'], ins[m][c]['z'], ins['n_used']] for m in ('dssp', 'ca_geometry') for c in cats])

# ---------------- Figure 3 source data ----------------
seg = J('step5/segment_length_check.json')
ct = J('step5/compositional_tests.json')
wcsv('Figure 3-source data 1.csv', ['set', 'n', 'pct_longer', 'median_gain_aa', 'gain_over_region_len_median'],
     [[k, v.get('n'), v.get('pct_A_longer'), v.get('median_span_gain_aa'), v.get('gain_over_region_len_median')] for k, v in seg.items() if isinstance(v, dict) and 'n' in v])
rows = []
for p, d in ct['test2_seeding'].items():
    for side in ('opposite_Pf_Nrich', 'opposite_Pf_other_classes'):
        x = d[side]
        rows.append([p, side, x['n'], x['median_N_frac'], x['median_K_frac'], x['pct_N_run_ge3'], x['pct_K_run_ge3']])
wcsv('Figure 3-source data 2.csv', ['partner', 'set', 'n', 'median_N_frac', 'median_K_frac', 'pct_N_run_ge3', 'pct_K_run_ge3'], rows)
rows = []
for s, d in ct['test1_N_vs_K'].items():
    for part in ('whole_proteome', 'homologous_disordered_segments', 'segments_without_tracts', 'tracts_in_segments'):
        if part in d:
            rows.append([s, d['expected_N_over_K'], part, d[part]['N'], d[part]['K'], d[part]['N_over_K'], d[part]['log2_obs_over_expected']])
wcsv('Figure 3-source data 3.csv', ['species', 'expected_N_over_K', 'set', 'N', 'K', 'N_over_K', 'log2_obs_over_expected'], rows)

# ---------------- Figure 4 supplement 1: gains by branch under three rootings ----------------
anc = J('step8/ancestral.json')
fig, axs = plt.subplots(1, 2, figsize=(7.2, 3.0), sharey=True)
branches = ['Pf', 'n_FP', 'n_FPR', 'n_FPRB', 'CladeB', 'Lav', 'Pv', 'Pk', 'Pb', 'Py', 'Pgal', 'Prel']
blab = ['Pf terminal', 'Pf+Pprf', '+Pr', '+Pbilc', 'clade B stem', 'Laverania stem', 'P. vivax', 'P. knowlesi', 'P. berghei', 'P. yoelii', 'P. gallinaceum', 'P. relictum']
for ax, lab in zip(axs, ['strict', 'lenient']):
    for j, (tr, col) in enumerate([('T1_mammal', BLUE), ('T2_avian', ORANGE), ('T3_basal', AQUA)]):
        v = [anc[lab][tr]['expected_gains_by_branch'].get(b, 0) for b in branches]
        ax.barh(np.arange(len(branches)) + (j - 1) * 0.27, v, 0.26, color=col, label=tr.replace('_', ': '))
    ax.set_yticks(range(len(branches)))
    ax.set_yticklabels(blab, fontsize=6.5)
    ax.invert_yaxis()
    ax.set_xlabel('Expected tract gains')
    title(ax, f'{lab} presence calls')
axs[0].legend(fontsize=6, loc='lower right')
fig.tight_layout()
save_tiff(fig, 'Figure4-figure supplement 1')
dat = J('step6/pf_tract_dating.json')
datl = J('step6/pf_tract_dating_lenient.json')
cd = J('step8/codon_decay.json')
cw = J('step8/codon_within.json')
wcsv('Figure 4-source data 1.csv', ['origin', 'strict_n', 'lenient_n'], [[k, dat['origin_counts'][k], datl['origin_counts'].get(k, '')] for k in dat['origin_counts']])
wcsv('Figure 4-source data 2.csv', ['age_class', 'kt_from', 'kt_to', 'n', 'observed_AAT', 'se', 'predicted_born_pure', 'predicted_born_at_young_mean'],
     [[c['cls'], c['kt_from'], c['kt_to'], c['n'], c['observed'], c['se'], c['predicted_born_pure'], c['predicted_born_at_young_mean']] for c in cd['classes']])
wcsv('Figure 4-source data 3.csv', ['comparison', 'n_tracts', 'background_D', 'within_tract_D', 'ratio', 'same_length_nn', 'same_length_diff', 'frac_length_changed'],
     [[k, v['n_tracts'], v['background_D'], v['within_tract_D'], v['ratio_within_vs_background'], v['same_length_nn'], v['same_length_diff'], v['frac_length_changed']] for k, v in cw.items()])
wcsv('Figure 4-source data 4.csv', ['calls', 'rooting', 'loglik', 'gain_rate', 'loss_rate', 'branch', 'branch_length', 'expected_gains', 'expected_losses'],
     [[lab, tr, d['loglik'], d['gain_rate'], d['loss_rate'], b, d['branch_lengths'].get(b), d['expected_gains_by_branch'].get(b), d['expected_losses_by_branch'].get(b)]
      for lab in anc for tr, d in anc[lab].items() for b in d['branch_lengths']])

# ---------------- Figure 5 supplement 1: per-population frequencies ----------------
TT = J('step10/pf8_young_tract_segregation.json')
CC = J('step10/pf8_control_indels.json')
POPS = ['AF-W', 'AF-C', 'AF-E', 'AS-SE-W', 'AS-SE-E', 'AS-S-FE', 'OC-NG', 'SA']
sites = sorted(TT['sites'], key=lambda s: -s['f'])
M = np.array([[s['popf'].get(p) if s['popf'].get(p) is not None else np.nan for p in POPS] for s in sites], float)
fig, ax = plt.subplots(figsize=(3.6, 4.2))
cm5 = plt.get_cmap('Blues').copy()
cm5.set_bad(LGREY)
im = ax.imshow(np.ma.masked_invalid(M), aspect='auto', cmap=cm5, vmin=0, vmax=1, interpolation='nearest')
ax.set_xticks(range(len(POPS)))
ax.set_xticklabels(POPS, rotation=45, ha='right', fontsize=6.5)
ax.set_yticks([])
ax.set_ylabel(f'{len(sites)} young-tract sites (sorted by overall frequency)')
cb = fig.colorbar(im, ax=ax, fraction=0.05, pad=0.03)
cb.set_label('Frequency of non-3D7 length allele', fontsize=6.5)
title(ax, 'Young-tract length alleles by population')
fig.tight_layout()
save_tiff(fig, 'Figure5-figure supplement 1')
pop = J('step7/pop_length_variation.json')
lm = pop['length_matched_pct_any5pct_and_mean_het']
wcsv('Figure 5-source data 1.csv', ['length_bin'] + [f'{g} (n: %any5pct / mean het)' for g in lm['columns']], [[b] + lm[b] for b in lm if b != 'columns'])


def site_rows(S, kind):
    out = []
    for s in S:
        ld = s.get('ld_ALL') or {}
        out.append([kind, s.get('chrom'), s.get('pos'), s.get('age'), s.get('n'), s.get('f'), s.get('fst'), s.get('het_frac'), s.get('nsnp'),
                    ld.get('maf'), ld.get('max_r2'), ld.get('ctrl_median_max_r2')] + [s['popf'].get(p) for p in POPS])
    return out


wcsv('Figure 5-source data 2.csv', ['set', 'chrom', 'pos', 'age', 'n_called', 'freq_non3D7_length', 'fst', 'mixed_call_frac', 'n_flank_snps',
                                     'ld_maf', 'ld_max_r2', 'ld_ctrl_median_max_r2'] + [f'freq_{p}' for p in POPS],
     site_rows(TT['sites'], 'young_tract') + site_rows(CC['sites'], 'control_indel'))
os.system(f"cp '{os.path.join(R, 'step7', 'pop_length_variation_from_table.tsv')}' '{os.path.join(SD, 'Figure 5-source data 3.tsv')}'")

# ---------------- Figure 6 supplement 1: MMR genes ----------------
mm = J('step9/mmr_check.json')
SPP = ['ppraefalciparum', 'preichenowi', 'pbillcollinsi', 'pblacklocki', 'padleri', 'pgaboni', 'pvivax', 'pknowlesi', 'pberghei', 'pyoelii', 'pgallinaceum', 'prelictum']
genes = list(mm['genes'])
P = np.array([[mm['genes'][g]['species'][s].get('pident_percentile', np.nan) if mm['genes'][g]['species'][s].get('ortholog') else np.nan for s in SPP] for g in genes], float)
fig, ax = plt.subplots(figsize=(5.2, 3.2))
im = ax.imshow(P, cmap='Blues', vmin=0, vmax=100, aspect='auto')
for i in range(P.shape[0]):
    for j in range(P.shape[1]):
        ax.text(j, i, '–' if np.isnan(P[i, j]) else f'{P[i, j]:.0f}', ha='center', va='center', fontsize=5.5, color='white' if P[i, j] > 60 else INK)
ax.set_xticks(range(len(SPP)))
ax.set_xticklabels([SPN[s] for s in SPP], rotation=45, ha='right', fontsize=6, fontstyle='italic')
ax.set_yticks(range(len(genes)))
ax.set_yticklabels([mm['genes'][g]['name'] for g in genes], fontsize=6.5)
cb = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02)
cb.set_label('Genome-wide identity percentile vs P. falciparum', fontsize=6)
title(ax, 'Mismatch-repair and replication genes')
fig.tight_layout()
save_tiff(fig, 'Figure6-figure supplement 1')
k = J('step7/constraint_vs_age.json')
wcsv('Figure 6-source data 1.csv', ['gene_class', 'n', 'pct_essential_MIS', 'median_MIS', 'median_MFS', 'median_expr_pct', 'pb_n', 'pb_pct_essential', 'median_protein_len'],
     [[g] + [v.get(x) for x in ['n', 'pct_essential_MIS', 'median_MIS', 'median_MFS', 'median_expr_pct', 'pb_n', 'pb_pct_essential', 'median_protein_len']] for g, v in k['Q1_gene_classes'].items()])
wcsv('Figure 6-source data 2.csv', ['gene', 'name', 'species', 'ortholog', 'pident', 'pident_percentile', 'pf_coverage', 'motif_present'],
     [[g, mm['genes'][g]['name'], s, e.get('ortholog'), e.get('pident'), e.get('pident_percentile'), e.get('pf_cov'), e.get('motif')]
      for g in genes for s, e in mm['genes'][g]['species'].items()])
for f in sorted(os.listdir(FIG)):
    if f.endswith('.tif'):
        im = Image.open(os.path.join(FIG, f))
        print(f, im.size, f'{im.size[0] / 300 * 2.54:.1f} cm wide at 300 dpi')
print(sorted(os.listdir(SD)))
