#!/usr/bin/env python3
"""Main figures for the poly-N manuscript. Reads only saved result files under results/.
Usage: python3 make_figures.py <results_dir>   (writes <results_dir>/figures/Figure{1..6}.png/.pdf)
"""
import csv, json, os, sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

R = sys.argv[1]
OUT = os.path.join(R, 'figures')
os.makedirs(OUT, exist_ok=True)

# palette (dataviz reference palette): categorical slots + one-hue ordinal ramp
BLUE, ORANGE, AQUA, YELLOW = '#2a78d6', '#eb6834', '#1baf7a', '#eda100'
GREY, INK, MUTED = '#b9b8b2', '#0b0b0b', '#52514e'
AGE_RAMP = ['#86b6ef', '#5598e7', '#3987e5', '#2a78d6', '#1c5cab', '#104281']
GROUP_COL = {'Laverania': BLUE, 'avian': ORANGE, 'rodent': AQUA, 'primate': YELLOW}
GROUP_MK = {'Laverania': 'o', 'avian': 's', 'rodent': '^', 'primate': 'D'}
AGES = ['Pf only', 'Pf+Pprf', '+Pr', '+Pbilc', 'clade B', 'all Laverania']

plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 8, 'axes.titlesize': 9,
                     'axes.labelsize': 8, 'axes.edgecolor': MUTED, 'axes.labelcolor': INK,
                     'xtick.color': MUTED, 'ytick.color': MUTED, 'axes.spines.top': False,
                     'axes.spines.right': False, 'axes.grid': True, 'grid.color': '#e6e5e1',
                     'grid.linewidth': 0.5, 'axes.axisbelow': True, 'legend.frameon': False})


def J(p):
    return json.load(open(os.path.join(R, p)))


def T(p):
    return list(csv.DictReader(open(os.path.join(R, p)), delimiter='\t'))


def letter(ax, s):
    ax.text(-0.14, 1.06, s, transform=ax.transAxes, fontsize=11, fontweight='bold', va='bottom', color=INK)


def save(fig, name):
    fig.savefig(os.path.join(OUT, name + '.png'), dpi=300, bbox_inches='tight')
    fig.savefig(os.path.join(OUT, name + '.pdf'), bbox_inches='tight')
    plt.close(fig)


SP = [('pfalciparum', 'P. falciparum', 'Laverania'), ('ppraefalciparum', 'P. praefalciparum', 'Laverania'),
      ('preichenowi', 'P. reichenowi', 'Laverania'), ('pbillcollinsi', 'P. billcollinsi', 'Laverania'),
      ('pblacklocki', 'P. blacklocki', 'Laverania'), ('padleri', 'P. adleri', 'Laverania'),
      ('pgaboni', 'P. gaboni', 'Laverania'), ('pgallinaceum', 'P. gallinaceum', 'avian'),
      ('prelictum', 'P. relictum', 'avian'), ('pberghei', 'P. berghei', 'rodent'),
      ('pyoelii', 'P. yoelii', 'rodent'), ('pvivax', 'P. vivax', 'primate'), ('pknowlesi', 'P. knowlesi', 'primate')]

# ---------------- Figure 1: Laverania genomes are rich in short tandem repeats ----------------
prev = {r['species']: r for r in T('step8/absolute_polyN_prevalence.tsv')}
gen = {r['species']: r for r in T('step9/genome_repeats.tsv')}
spec = {r['species']: r for r in T('step9/repeat_length_spectrum.tsv')}
fig, axs = plt.subplots(1, 3, figsize=(7.2, 2.6), gridspec_kw={'width_ratios': [1.25, 1, 1]})
ax = axs[0]
names = [s[1] for s in SP][::-1]
vals = [float(prev[s[0]]['pct_win20_N60']) for s in SP][::-1]
cols = [GROUP_COL[s[2]] for s in SP][::-1]
ax.barh(range(len(SP)), vals, color=cols, height=0.7)
ax.set_yticks(range(len(SP)))
ax.set_yticklabels(names, fontstyle='italic', fontsize=7)
ax.set_xlabel('Proteins with a ≥60% N window (%)')
ax.grid(axis='y', visible=False)
for i, v in enumerate(vals):
    ax.text(v + 0.5, i, f'{v:.1f}', va='center', fontsize=6.5, color=MUTED)
ax.set_xlim(0, 33)
letter(ax, 'A')
ax = axs[1]
for s, nm, g in SP:
    x, y = float(gen[s]['nc_AT']), float(gen[s]['dAT8_perMb'])
    ax.scatter(x, y, s=28, color=GROUP_COL[g], marker=GROUP_MK[g], edgecolor='white', linewidth=0.6, zorder=3)
ax.set_yscale('log')
ax.set_xlabel('Non-coding AT content')
ax.set_ylabel('(AT)n ≥ 8 units per Mb (non-coding)')
for g, (x, y) in {'Laverania': (0.80, 3300), 'avian': (0.858, 300), 'rodent': (0.745, 230), 'primate': (0.635, 45)}.items():
    ax.text(x, y, g, color=GROUP_COL[g], fontsize=7, fontweight='bold')
letter(ax, 'B')
ax = axs[2]
lav = ['pfalciparum', 'preichenowi', 'padleri']
av = ['pgallinaceum', 'prelictum']
cols_u = ['ATn_3u', 'ATn_4u', 'ATn_5u', 'ATn_6u', 'ATn_7u', 'ATn_8to11u', 'ATn_12plus']
labs = ['3', '4', '5', '6', '7', '8–11', '≥12']
ratio = [np.mean([float(spec[s][c]) for s in lav]) / np.mean([float(spec[s][c]) for s in av]) for c in cols_u]
cols_t = ['AATtype_3u', 'AATtype_4u', 'AATtype_5to7u', 'AATtype_8plus']
ratio_t = [np.mean([float(spec[s][c]) for s in lav]) / np.mean([float(spec[s][c]) for s in av]) for c in cols_t]
ax.plot(range(len(labs)), ratio, '-o', color=BLUE, lw=2, ms=4, label='(AT)n')
ax.plot([1, 2, 3.5, 6], ratio_t, '-s', color=ORANGE, lw=2, ms=4, label='AAT/ATT-type')
ax.axhline(1, color=MUTED, lw=0.8, ls='--')
ax.set_yscale('log')
ax.set_yticks([1, 2, 4, 10])
ax.set_yticklabels(['1', '2', '4', '10'])
ax.set_xticks(range(len(labs)))
ax.set_xticklabels(labs, rotation=30)
ax.set_xlabel('Repeat length (units)')
ax.set_ylabel('Laverania / avian density')
ax.legend(loc='upper left', fontsize=7)
letter(ax, 'C')
fig.tight_layout()
save(fig, 'Figure1')

# ---------------- Figure 2: tracts sit in conserved disordered segments ----------------
s5 = J('step5/step5_results.json')['A_location_by_class']
ins = J('step8/insertion_site_dssp.json')
fig, axs = plt.subplots(1, 2, figsize=(7.2, 2.6), gridspec_kw={'width_ratios': [1.4, 1]})
ax = axs[0]
classes = ['N-rich', 'Q-rich', 'GS-rich', 'charged', 'other']
cat = [('homologous segment', ['internal_variable', 'terminal_variable'], BLUE),
       ('true insertion / extension', ['internal_insertion', 'terminal_extension'], ORANGE),
       ('folded core', ['core_shared', 'core_overlap'], GREY)]
left = np.zeros(len(classes))
for lab, keys, col in cat:
    v = []
    for c in classes:
        d = s5[c]
        tot = sum(d.get(k, 0) for _, ks, _ in cat for k in ks)
        v.append(100 * sum(d.get(k, 0) for k in keys) / tot)
    ax.barh(range(len(classes)), v, left=left, color=col, height=0.62, label=lab, edgecolor='white', linewidth=1)
    if lab == 'true insertion / extension':
        for i, x in enumerate(v):
            ax.text(100.4, i, f'{x:.1f}%', va='center', fontsize=6.5, color=ORANGE, fontweight='bold')
    left += v
ax.set_xlim(85, 102.2)
ax.set_xticks([86, 88, 90, 92, 94, 96, 98, 100])
ax.set_yticks(range(len(classes)))
ax.set_yticklabels(classes)
ax.invert_yaxis()
ax.set_xlabel('Placeable region–ortholog comparisons (%)  [axis starts at 85%]')
ax.grid(axis='y', visible=False)
ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.28), ncol=3, fontsize=6.5)
letter(ax, 'A')
ax = axs[1]
cats = ['element', 'boundary', 'loop']
labs = ['inside\nhelix/strand', 'element\nedge', 'loop']
x = np.arange(3)
for off, key, nm, col in [(-0.2, 'dssp', 'insertion sites', ORANGE), (0.2, 'dssp', 'random sites', GREY)]:
    vals = [100 * ins[key][c]['obs' if nm == 'insertion sites' else 'null'] for c in cats]
    ax.bar(x + off, vals, width=0.38, color=col, label=nm)
ax.set_xticks(x)
ax.set_xticklabels(labs, fontsize=7)
ax.set_ylabel('Sites (%) — DSSP')
ax.text(-0.2, 2, '0%', ha='center', fontsize=6.5, color=INK)
ax.text(0, 30, f"z = {ins['dssp']['element']['z']:.1f}".replace('-', '−'), ha='center', fontsize=7, color=INK)
ax.legend(fontsize=7, loc='upper left')
ax.grid(axis='x', visible=False)
letter(ax, 'B')
fig.tight_layout()
save(fig, 'Figure2')

# ---------------- Figure 3: Pf tracts are expansions from asparagine seeds ----------------
seg = J('step5/segment_length_check.json')
ct = J('step5/compositional_tests.json')
fig, axs = plt.subplots(1, 3, figsize=(7.2, 2.5))
ax = axs[0]
rows = [('Pf N-rich', 'Pf_Nrich', BLUE), ('Pf GS-rich', 'Pf_GS-rich', GREY), ('Pf charged', 'Pf_charged', GREY),
        ('Pf other', 'Pf_other', GREY), ('non-Pf N-rich', 'nonPf_Nrich', AQUA), ('non-Pf other', 'nonPf_other_classes', AQUA)]
v = [seg[k]['median_span_gain_aa'] for _, k, _ in rows]
ax.barh(range(len(rows)), v, color=[c for *_, c in rows], height=0.65)
for i, (_, k, _) in enumerate(rows):
    ax.text(v[i] + 1, i, f"{seg[k]['pct_A_longer']:.0f}% longer", va='center', fontsize=6.3, color=MUTED)
ax.set_yticks(range(len(rows)))
ax.set_yticklabels([r[0] for r in rows], fontsize=7)
ax.invert_yaxis()
ax.set_xlim(0, 75)
ax.set_xlabel('Median extra length of segment\nvs ortholog (aa)')
ax.grid(axis='y', visible=False)
letter(ax, 'A')
ax = axs[1]
seed = ct['test2_seeding']
partners = ['Pb', 'Pk', 'Pv']
x = np.arange(3)
a = [seed[p]['opposite_Pf_Nrich']['pct_N_run_ge3'] for p in partners]
b = [seed[p]['opposite_Pf_other_classes']['pct_N_run_ge3'] for p in partners]
ax.bar(x - 0.2, a, 0.38, color=BLUE, label='opposite Pf poly-N')
ax.bar(x + 0.2, b, 0.38, color=GREY, label='opposite other Pf LCRs')
ax.set_xticks(x)
ax.set_xticklabels(['P. berghei', 'P. knowlesi', 'P. vivax'], fontstyle='italic', fontsize=6.5, rotation=20)
ax.set_ylabel('Ortholog segments with\nan N run ≥ 3 (%)')
ax.legend(fontsize=6.5, loc='upper right')
ax.grid(axis='x', visible=False)
letter(ax, 'B')
ax = axs[2]
nk = ct['test1_N_vs_K']
sp4 = ['Pf', 'Pb', 'Pk', 'Pv']
x = np.arange(4)
bg = [nk[s]['segments_without_tracts']['log2_obs_over_expected'] for s in sp4]
tr = [nk[s]['tracts_in_segments']['log2_obs_over_expected'] for s in sp4]
ax.bar(x - 0.2, bg, 0.38, color=GREY, label='disordered background')
ax.bar(x + 0.2, tr, 0.38, color=BLUE, label='tracts')
ax.axhline(0, color=MUTED, lw=0.8)
ax.set_xticks(x)
ax.set_xticklabels(['P. falciparum', 'P. berghei', 'P. knowlesi', 'P. vivax'], fontstyle='italic', fontsize=6.5, rotation=20)
ax.set_ylabel('N:K vs mutation-bias\nexpectation (log2 O/E)')
ax.set_ylim(-0.3, 1.9)
ax.legend(fontsize=6.5, loc='upper right')
ax.grid(axis='x', visible=False)
letter(ax, 'C')
fig.tight_layout()
save(fig, 'Figure3')

# ---------------- Figure 4: dating, gains, codon decay ----------------
dat = J('step6/pf_tract_dating.json')
datl = J('step6/pf_tract_dating_lenient.json')
anc = J('step8/ancestral.json')
cd = J('step8/codon_decay.json')
cw = J('step8/codon_within.json')
fig, axs = plt.subplots(1, 3, figsize=(7.2, 2.6))
ax = axs[0]
keys = ['Pf only', 'Pf+Pprf', '+Pr', '+Pbilc', 'clade B (+Pblac)', 'all Laverania (+clade A)']
st = [dat['origin_counts'][k] for k in keys]
le = [datl['origin_counts'][k] for k in keys]
x = np.arange(6)
ax.bar(x - 0.2, st, 0.38, color=AGE_RAMP, label='strict')
ax.bar(x + 0.2, le, 0.38, color=AGE_RAMP, alpha=0.45, label='lenient', hatch='///', edgecolor='white')
ax.set_xticks(x)
ax.set_xticklabels(AGES, rotation=35, ha='right', fontsize=6.5)
ax.set_ylabel('Pf poly-N tracts (by origin)')
ax.legend(fontsize=6.5)
ax.grid(axis='x', visible=False)
letter(ax, 'A')
ax = axs[1]
cl = cd['classes']
obs = [c['observed'] for c in cl]
se = [c['se'] for c in cl]
pred = [c['predicted_born_pure'] for c in cl]
ax.errorbar(x, obs, yerr=[1.96 * s for s in se], fmt='o', color=BLUE, ms=4, capsize=2, label='observed')
ax.plot(x, pred, 's--', color=ORANGE, ms=3.5, lw=1.2, label='neutral prediction')
ax.axhline(cd['p_eq'], color=MUTED, lw=0.8, ls=':')
ax.text(-0.2, cd['p_eq'] - 0.004, 'genome equilibrium', fontsize=6, color=MUTED, ha='left', va='top')
ax.set_xticks(x)
ax.set_xticklabels(AGES, rotation=35, ha='right', fontsize=6.5)
ax.set_ylabel('AAT share of Asn codons')
ax.set_ylim(0.82, 1.01)
ax.legend(fontsize=6.3, loc='upper center', bbox_to_anchor=(0.5, 1.2), ncol=2)
letter(ax, 'B')
ax = axs[2]
labs = ['vs P. praefalciparum', 'vs P. reichenowi']
bgd = [cw['Pprf']['background_D'], cw['Pr']['background_D']]
wt = [cw['Pprf']['within_tract_D'], cw['Pr']['within_tract_D']]
x2 = np.arange(2)
ax.bar(x2 - 0.2, [100 * v for v in bgd], 0.38, color=GREY, label='genome N codons')
ax.bar(x2 + 0.2, [100 * v for v in wt], 0.38, color=BLUE, label='inside shared tracts')
for i in range(2):
    ax.text(i + 0.2, 100 * wt[i] + 0.2, f"{wt[i] / bgd[i]:.1f}×", ha='center', fontsize=7, color=INK)
ax.set_xticks(x2)
ax.set_xticklabels(labs, fontsize=6.5)
ax.set_ylabel('AAT↔AAC divergence (%)')
ax.legend(fontsize=6.5, loc='upper left')
ax.grid(axis='x', visible=False)
letter(ax, 'C')
fig.tight_layout()
save(fig, 'Figure4')

# ---------------- Figure 5: population turnover in Pf8 ----------------
pop = J('step7/pop_length_variation.json')
cc_t = J('step10/pf8_young_tract_segregation.json')['sites']
cc_c = J('step10/pf8_control_indels.json')['sites']
fig, axs = plt.subplots(1, 3, figsize=(7.2, 2.6), gridspec_kw={'width_ratios': [1.3, 1, 1]})
ax = axs[0]
lm = pop['length_matched_pct_any5pct_and_mean_het']
bins = ['<18 aa', '18-24', '25-34', '35-49', '>=50']
groups = lm['columns']
gcol = [AGE_RAMP[0], AGE_RAMP[2], AGE_RAMP[4], AGE_RAMP[5]]
x = np.arange(len(bins))
for j, g in enumerate(groups):
    v = [float(lm[b][j].split(': ')[1].split('%')[0]) for b in bins]
    ax.plot(x, v, '-o', color=gcol[j], ms=3.5, lw=1.6, label=g)
ax.set_xticks(x)
ax.set_xticklabels(['<18', '18–24', '25–34', '35–49', '≥50'])
ax.set_xlabel('Tract length (aa)')
ax.set_ylabel('Tracts with a length\nvariant ≥5% in Pf8 (%)')
ax.legend(fontsize=6.0, loc='upper center', bbox_to_anchor=(0.5, -0.3), ncol=2)
letter(ax, 'A')


def both(x):
    af = max(x['popf'].get(p) or 0 for p in ['AF-W', 'AF-C', 'AF-E'])
    a = max(x['popf'].get(p) or 0 for p in ['AS-SE-W', 'AS-SE-E'])
    return af >= 0.01 and a >= 0.01


ax = axs[1]
ft = [x['fst'] for x in cc_t if x['fst'] is not None]
fc = [x['fst'] for x in cc_c if x['fst'] is not None]
snp = J('step10/pf8_young_tract_segregation.json')['snp_fst'] + J('step10/pf8_control_indels.json')['snp_fst']
parts = ax.boxplot([ft, fc, snp], widths=0.55, patch_artist=True, showfliers=False,
                   medianprops=dict(color=INK, lw=1.2))
for p, c in zip(parts['boxes'], [BLUE, ORANGE, GREY]):
    p.set_facecolor(c)
    p.set_edgecolor('white')
ax.set_xticks([1, 2, 3])
ax.set_xticklabels(['young-tract\nvariants', 'control\nindels', 'flanking\nSNPs'], fontsize=6.5)
ax.set_ylabel('FST across 8 populations')
ax.grid(axis='x', visible=False)
letter(ax, 'B')
ax = axs[2]


def ldfrac(S):
    L = [s['ld_ALL'] for s in S if s.get('ld_ALL') and s['ld_ALL']['max_r2'] is not None and s['ld_ALL']['ctrl_median_max_r2'] is not None]
    return 100 * sum(l['max_r2'] < l['ctrl_median_max_r2'] for l in L) / len(L)


metrics = [('shared Africa\n+ SE Asia', 100 * np.mean([both(s) for s in cc_t]), 100 * np.mean([both(s) for s in cc_c])),
           ('weaker LD than\nmatched SNPs', ldfrac(cc_t), ldfrac(cc_c))]
x = np.arange(2)
ax.bar(x - 0.2, [m[1] for m in metrics], 0.38, color=BLUE, label='young-tract variants')
ax.bar(x + 0.2, [m[2] for m in metrics], 0.38, color=ORANGE, label='control indels')
ax.axhline(50, color=MUTED, lw=0.8, ls=':')
ax.set_xticks(x)
ax.set_xticklabels([m[0] for m in metrics], fontsize=6.5)
ax.set_ylabel('Sites (%)')
ax.set_ylim(0, 100)
ax.legend(fontsize=6.3, loc='upper right')
ax.grid(axis='x', visible=False)
letter(ax, 'C')
fig.tight_layout()
save(fig, 'Figure5')

# ---------------- Figure 6: no functional retention ----------------
k = J('step7/constraint_vs_age.json')
fig, axs = plt.subplots(1, 2, figsize=(5.2, 2.5))
ax = axs[0]
gc = k['Q1_gene_classes']
lab = [('no tract', 'no_tract'), ('young only', 'young_only'), ('mixed', 'mixed'), ('old only', 'old_only')]
v = [gc[g]['pct_essential_MIS'] for _, g in lab]
vb = [gc[g]['pb_pct_essential'] for _, g in lab]
x = np.arange(4)
ax.bar(x - 0.2, v, 0.38, color=BLUE, label='P. falciparum (piggyBac)')
ax.bar(x + 0.2, vb, 0.38, color=AQUA, label='P. berghei ortholog')
ax.set_xticks(x)
ax.set_xticklabels(['no\ntract', 'young\nonly', 'mixed', 'old\nonly'], fontsize=7)
ax.set_ylabel('Essential genes (%)')
ax.set_ylim(0, 70)
ax.legend(fontsize=6.3, loc='upper left')
ax.grid(axis='x', visible=False)
letter(ax, 'A')
ax = axs[1]
q3 = k['Q3_old_tracts_shared_outside_laverania']
ax.bar([0, 1], [q3['essential']['pct_shared_outside'], q3['dispensable']['pct_shared_outside']], 0.55, color=[BLUE, GREY])
ax.set_xticks([0, 1])
ax.set_xticklabels([f"essential\n(n = {q3['essential']['n_old']})", f"dispensable\n(n = {q3['dispensable']['n_old']})"], fontsize=7)
ax.set_ylabel('Old tracts also N-rich\noutside Laverania (%)')
ax.set_ylim(0, 60)
ax.grid(axis='x', visible=False)
letter(ax, 'B')
fig.tight_layout()
save(fig, 'Figure6')
print('done', os.listdir(OUT))
