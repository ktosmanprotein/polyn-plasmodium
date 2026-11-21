#!/usr/bin/env python3
"""Main figures, v2 (publication styling). Reads only saved result files under results/.
Usage: python3 make_figures_v2.py <results_dir>  -> <results_dir>/figures_v2/Figure{1..6}.png/.pdf
Design: one accent colour per figure carries the finding, grey carries context; every panel has a
one-line finding as its title; direct labels instead of legends; key statistic printed on the panel.
"""
import csv, json, os, sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle
import numpy as np

R = sys.argv[1]
OUT = os.path.join(R, 'figures_v2')
os.makedirs(OUT, exist_ok=True)

BLUE, ORANGE, AQUA, YELLOW = '#2a78d6', '#eb6834', '#1baf7a', '#eda100'
DEEP = '#104281'
GREY, LGREY, INK, MUTED = '#a9a8a2', '#e9e8e4', '#1a1a19', '#5c5b57'
AGE = ['#9ec5f4', '#6da7ec', '#3987e5', '#256abf', '#184f95', '#0d366b']
GCOL = {'Laverania': BLUE, 'avian': ORANGE, 'rodent': AQUA, 'primate': YELLOW}
GMK = {'Laverania': 'o', 'avian': 's', 'rodent': '^', 'primate': 'D'}

plt.rcParams.update({
    'font.family': 'DejaVu Sans', 'font.size': 7.5, 'axes.titlesize': 8, 'axes.labelsize': 7.5,
    'axes.edgecolor': MUTED, 'axes.linewidth': 0.6, 'axes.labelcolor': INK, 'text.color': INK,
    'xtick.color': MUTED, 'ytick.color': MUTED, 'xtick.major.width': 0.6, 'ytick.major.width': 0.6,
    'xtick.major.size': 2.5, 'ytick.major.size': 2.5, 'axes.spines.top': False, 'axes.spines.right': False,
    'axes.grid': False, 'legend.frameon': False, 'svg.fonttype': 'none', 'pdf.fonttype': 42})


def J(p):
    return json.load(open(os.path.join(R, p)))


def T(p):
    return list(csv.DictReader(open(os.path.join(R, p)), delimiter='\t'))


def panel(ax, letter, title):
    import textwrap
    fig = ax.figure
    w_in = ax.get_position().width * fig.get_size_inches()[0]
    lines = textwrap.wrap(title, width=max(14, int(w_in * 13.5)))
    top = 1.0 + 0.62 / ax.get_position().height / fig.get_size_inches()[1]
    ax.text(0.0, top, '\n'.join(lines), transform=ax.transAxes, fontsize=7.6, fontweight='bold', va='top', ha='left',
            color=INK, linespacing=1.15)
    ax.text(-0.04, top + 0.02, letter, transform=ax.transAxes, fontsize=12, fontweight='bold', va='top', ha='right')


def ygrid(ax):
    ax.yaxis.grid(True, color=LGREY, lw=0.6)
    ax.set_axisbelow(True)


def xgrid(ax):
    ax.xaxis.grid(True, color=LGREY, lw=0.6)
    ax.set_axisbelow(True)


def save(fig, name):
    fig.savefig(os.path.join(OUT, name + '.png'), dpi=400, bbox_inches='tight', facecolor='white')
    fig.savefig(os.path.join(OUT, name + '.pdf'), bbox_inches='tight', facecolor='white')
    plt.close(fig)


def pill(ax, x, y, s, color, trans=None, fs=7, ha='center'):
    ax.text(x, y, s, transform=trans or ax.transAxes, ha=ha, va='center', fontsize=fs, color='white', fontweight='bold',
            bbox=dict(boxstyle='round,pad=0.35,rounding_size=0.6', fc=color, ec='none'))


def dumbbell(ax, y, a, b, ca, cb, lw=2.2):
    ax.plot([a, b], [y, y], color=LGREY, lw=lw + 1.5, solid_capstyle='round', zorder=1)
    ax.scatter([a], [y], s=34, color=ca, zorder=3, edgecolor='white', linewidth=0.7)
    ax.scatter([b], [y], s=34, color=cb, zorder=3, edgecolor='white', linewidth=0.7)


SP = [('pfalciparum', 'P. falciparum', 'Laverania'), ('ppraefalciparum', 'P. praefalciparum', 'Laverania'),
      ('preichenowi', 'P. reichenowi', 'Laverania'), ('pbillcollinsi', 'P. billcollinsi', 'Laverania'),
      ('pblacklocki', 'P. blacklocki', 'Laverania'), ('padleri', 'P. adleri', 'Laverania'),
      ('pgaboni', 'P. gaboni', 'Laverania'), ('pgallinaceum', 'P. gallinaceum', 'avian'),
      ('prelictum', 'P. relictum', 'avian'), ('pberghei', 'P. berghei', 'rodent'),
      ('pyoelii', 'P. yoelii', 'rodent'), ('pvivax', 'P. vivax', 'primate'), ('pknowlesi', 'P. knowlesi', 'primate')]

# =========================== Figure 1 ===========================
prev = {r['species']: r for r in T('step8/absolute_polyN_prevalence.tsv')}
gen = {r['species']: r for r in T('step9/genome_repeats.tsv')}
spec = {r['species']: r for r in T('step9/repeat_length_spectrum.tsv')}
fig = plt.figure(figsize=(7.4, 2.9))
gs = fig.add_gridspec(1, 3, width_ratios=[1.35, 1, 1], wspace=0.55)
ax = fig.add_subplot(gs[0])
ys = np.arange(len(SP))[::-1].astype(float)
# group gaps
off = {'Laverania': 0, 'avian': -0.6, 'rodent': -1.2, 'primate': -1.8}
ys = np.array([ys[i] + off[s[2]] for i, s in enumerate(SP)])
for (key, nm, g), y in zip(SP, ys):
    v = float(prev[key]['pct_win20_N60'])
    ax.plot([0, v], [y, y], color=GCOL[g], lw=2.6, solid_capstyle='round', alpha=0.9)
    ax.scatter(v, y, s=26, color=GCOL[g], zorder=3, edgecolor='white', linewidth=0.6)
    ax.text(v + 1.0, y, f'{v:.1f}', va='center', fontsize=6, color=MUTED)
ax.set_yticks(ys)
ax.set_yticklabels([s[1] for s in SP], fontstyle='italic', fontsize=6.6)
ax.tick_params(axis='y', length=0)
ax.spines['left'].set_visible(False)
ax.set_xlim(0, 52)
ax.set_xticks([0, 10, 20, 30])
ax.set_xlabel('Proteins with a ≥60% asparagine window (%)')
xgrid(ax)
for tl, (k, n, g) in zip(ax.get_yticklabels(), SP):
    tl.set_color(GCOL[g])
for g in ['Laverania', 'avian', 'rodent', 'primate']:
    yy = [y for (k, n, gg), y in zip(SP, ys) if gg == g]
    ax.text(40, max(yy) if g == 'Laverania' else np.mean(yy), g, ha='left', va='center', fontsize=6.6, color=GCOL[g], fontweight='bold')
panel(ax, 'A', 'One in four Laverania proteins is poly-N')
ax = fig.add_subplot(gs[1])
ax.axhspan(0.5, 8, color=LGREY, zorder=0, lw=0)
ax.text(0.585, 3.2, 'Markov-model expectation', fontsize=5.8, color=MUTED, va='center')
for key, nm, g in SP:
    ax.scatter(float(gen[key]['nc_AT']), float(gen[key]['dAT8_perMb']), s=30, color=GCOL[g], marker=GMK[g],
               edgecolor='white', linewidth=0.6, zorder=3)
ax.set_yscale('log')
ax.set_ylim(0.6, 9000)
ax.set_xlim(0.58, 0.9)
lavy = np.mean([float(gen[k]['dAT8_perMb']) for k, _, g in SP if g == 'Laverania'])
avy = np.mean([float(gen[k]['dAT8_perMb']) for k, _, g in SP if g == 'avian'])
ax.annotate('', xy=(0.878, lavy), xytext=(0.878, avy), arrowprops=dict(arrowstyle='<->', color=INK, lw=0.8))
ax.text(0.872, np.sqrt(lavy * avy), '6.5×\nsame AT', ha='right', va='center', fontsize=6.4, fontweight='bold', color=INK)
for g, (x, y) in {'Laverania': (0.79, 4600), 'avian': (0.815, 190), 'rodent': (0.745, 75), 'primate': (0.615, 230)}.items():
    ax.text(x, y, g, color=GCOL[g], fontsize=6.6, fontweight='bold')
ax.set_xlabel('Non-coding AT content')
ax.set_ylabel('(AT)n repeats ≥ 8 units per Mb')
ygrid(ax)
panel(ax, 'B', 'Laverania DNA is rich in repeats')
ax = fig.add_subplot(gs[2])
lav = ['pfalciparum', 'preichenowi', 'padleri']
av = ['pgallinaceum', 'prelictum']
cu = ['ATn_3u', 'ATn_4u', 'ATn_5u', 'ATn_6u', 'ATn_7u', 'ATn_8to11u', 'ATn_12plus']
labs = ['3', '4', '5', '6', '7', '8–11', '12+']
ratio = [np.mean([float(spec[s][c]) for s in lav]) / np.mean([float(spec[s][c]) for s in av]) for c in cu]
ct = ['AATtype_3u', 'AATtype_4u', 'AATtype_5to7u', 'AATtype_8plus']
ratio_t = [np.mean([float(spec[s][c]) for s in lav]) / np.mean([float(spec[s][c]) for s in av]) for c in ct]
x = np.arange(len(labs))
ax.axvspan(-0.4, 1.4, color='#eef4fc', lw=0, zorder=0)
ax.text(0.5, 13.5, '≤ 8 bp', ha='center', va='top', fontsize=6, color=BLUE, fontweight='bold')
ax.plot(x, ratio, '-o', color=BLUE, lw=2, ms=4, mec='white', mew=0.6)
ax.plot([1, 2, 3.5, 6], ratio_t, '-s', color=ORANGE, lw=1.6, ms=3.6, mec='white', mew=0.6)
ax.text(6.3, ratio[-1], '(AT)n', color=BLUE, fontsize=6.4, fontweight='bold', va='center')
ax.text(6.3, ratio_t[-1], 'AAT-type', color=ORANGE, fontsize=6.4, fontweight='bold', va='center')
ax.axhline(1, color=MUTED, lw=0.7, ls=(0, (2, 2)))
ax.set_yscale('log')
ax.set_ylim(0.8, 15)
ax.set_yticks([1, 2, 4, 10])
ax.set_yticklabels(['1×', '2×', '4×', '10×'])
ax.set_xticks(x)
ax.set_xticklabels(labs, fontsize=6.2, rotation=35)
ax.set_xlim(-0.4, 8.0)
ax.set_xlabel('Repeat length (units)')
ax.set_ylabel('Laverania / avian density')
ygrid(ax)
panel(ax, 'C', 'Not an assembly artefact')
save(fig, 'Figure1')

# =========================== Figure 2 ===========================
s5 = J('step5/step5_results.json')['A_location_by_class']
ins = J('step8/insertion_site_dssp.json')
fig = plt.figure(figsize=(7.4, 2.7))
gs = fig.add_gridspec(1, 3, width_ratios=[0.9, 1.35, 1], wspace=0.6)
# schematic
ax = fig.add_subplot(gs[0])
ax.set_xlim(0, 10)
ax.set_ylim(0, 10)
ax.axis('off')


def prot(ax, y, segs, lab):
    for (a, b, kind) in segs:
        if kind == 'core':
            ax.add_patch(FancyBboxPatch((a, y - 0.45), b - a, 0.9, boxstyle='round,pad=0,rounding_size=0.25', fc=GREY, ec='none'))
        elif kind == 'dis':
            xs = np.linspace(a, b, 60)
            ax.plot(xs, y + 0.25 * np.sin(xs * 6), color=MUTED, lw=1.1)
        elif kind == 'N':
            xs = np.linspace(a, b, 60)
            ax.plot(xs, y + 0.25 * np.sin(xs * 6), color=BLUE, lw=2.4)
    ax.text(-0.2, y, lab, ha='right', va='center', fontsize=6.6, fontstyle='italic')


prot(ax, 8.0, [(0.3, 3.0, 'core'), (3.0, 7.2, 'dis'), (4.2, 6.6, 'N'), (7.2, 9.8, 'core')], 'Pf')
prot(ax, 6.2, [(0.3, 3.0, 'core'), (3.0, 4.6, 'dis'), (4.6, 7.2, 'core')], 'ortholog')
for xx in (3.0, ):
    ax.plot([xx, xx], [6.6, 7.6], color=INK, lw=0.6, ls=(0, (1.5, 1.5)))
ax.plot([7.2, 4.6], [7.6, 6.6], color=INK, lw=0.6, ls=(0, (1.5, 1.5)))
ax.text(5.4, 9.15, 'poly-N', color=BLUE, fontsize=6.6, fontweight='bold', ha='center')
ax.text(5.0, 4.7, 'Anchors come only from the\nconfident core (pLDDT ≥ 70);\ndisorder is never aligned\nto disorder.', ha='center', va='top', fontsize=6.0, color=MUTED, linespacing=1.3)
ax.text(1.65, 9.15, 'folded core', ha='center', fontsize=6, color=MUTED)
panel(ax, 'A', 'Structure-anchored placement')
ax = fig.add_subplot(gs[1])
classes = ['N-rich', 'Q-rich', 'GS-rich', 'charged', 'other']
cats = [('homologous segment', ['internal_variable', 'terminal_variable'], '#c9dcf3'),
        ('true insertion', ['internal_insertion', 'terminal_extension'], ORANGE),
        ('folded core', ['core_shared', 'core_overlap'], GREY)]
for i, c in enumerate(classes):
    d = s5[c]
    tot = sum(d.get(k, 0) for _, ks, _ in cats for k in ks)
    left = 0
    for lab, keys, col in cats:
        v = 100 * sum(d.get(k, 0) for k in keys) / tot
        ax.barh(i, v, left=left, color=col, height=0.66, edgecolor='white', linewidth=1)
        if lab == 'homologous segment':
            ax.text(left + v - 0.3, i, f'{v:.0f}%', ha='right', va='center', fontsize=6.3, color=DEEP, fontweight='bold')
        if lab == 'true insertion':
            ax.text(100.6, i, f'{v:.1f}%', ha='left', va='center', fontsize=6.3, color=ORANGE, fontweight='bold')
        left += v
ax.set_xlim(88, 104)
ax.set_xticks([88, 92, 96, 100])
ax.set_xticklabels(['88', '92', '96', '100%'])
ax.set_yticks(range(5))
ax.set_yticklabels(classes)
ax.tick_params(axis='y', length=0)
ax.spines['left'].set_visible(False)
ax.invert_yaxis()
ax.set_xlabel('Placeable comparisons (axis from 88%)')
panel(ax, 'B', '94–98% sit opposite a segment the ortholog has')
ax = fig.add_subplot(gs[2])
cats3 = ['element', 'boundary', 'loop']
lab3 = ['inside helix\nor strand', 'element\nedge', 'loop']
obs = [100 * ins['dssp'][c]['obs'] for c in cats3]
nul = [100 * ins['dssp'][c]['null'] for c in cats3]
for i in range(3):
    dumbbell(ax, i, nul[i], obs[i], GREY, ORANGE)
ax.set_yticks(range(3))
ax.set_yticklabels(lab3, fontsize=6.6)
ax.tick_params(axis='y', length=0)
ax.spines['left'].set_visible(False)
ax.invert_yaxis()
ax.set_xlim(-5, 95)
ax.set_xlabel('Sites (%), DSSP')
xgrid(ax)
ax.text(nul[0] + 3, -0.22, 'random sites', fontsize=6, color=MUTED)
ax.text(obs[0] + 1, 0.3, 'insertion sites: 0 of 195', fontsize=6, color=ORANGE, fontweight='bold', va='top')
pill(ax, 0.72, 0.12, 'z = −6.5', ORANGE)
panel(ax, 'C', 'Insertions avoid helices and strands')
save(fig, 'Figure2')

# =========================== Figure 3 ===========================
seg = J('step5/segment_length_check.json')
ct3 = J('step5/compositional_tests.json')
fig = plt.figure(figsize=(7.4, 2.7))
gs = fig.add_gridspec(1, 3, width_ratios=[1.1, 1, 1], wspace=0.62)
ax = fig.add_subplot(gs[0])
rows = [('Pf poly-N', 'Pf_Nrich', BLUE), ('Pf GS-rich', 'Pf_GS-rich', GREY), ('Pf charged', 'Pf_charged', GREY),
        ('Pf other', 'Pf_other', GREY), ('non-Pf poly-N', 'nonPf_Nrich', GREY), ('non-Pf other', 'nonPf_other_classes', GREY)]
for i, (lab, k, col) in enumerate(rows):
    v = seg[k]['median_span_gain_aa']
    ax.plot([0, v], [i, i], color=col, lw=2.6, solid_capstyle='round')
    ax.scatter(v, i, s=26, color=col, zorder=3, edgecolor='white', linewidth=0.6)
    ax.text(v + 2, i, f"+{v:.0f} aa", va='center', fontsize=6.2, color=col if col == BLUE else MUTED,
            fontweight='bold' if col == BLUE else 'normal')
ax.set_yticks(range(len(rows)))
ax.set_yticklabels([r[0] for r in rows], fontsize=6.6)
ax.tick_params(axis='y', length=0)
ax.spines['left'].set_visible(False)
ax.invert_yaxis()
ax.set_xlim(0, 62)
ax.set_xlabel('Median extra segment length vs ortholog')
xgrid(ax)
pill(ax, 0.62, 0.86, '86% longer', BLUE)
panel(ax, 'A', 'Pf tracts are added sequence')
ax = fig.add_subplot(gs[1])
seed = ct3['test2_seeding']
for i, (p, nm) in enumerate([('Pb', 'P. berghei'), ('Pk', 'P. knowlesi'), ('Pv', 'P. vivax')]):
    a = seed[p]['opposite_Pf_Nrich']['pct_N_run_ge3']
    b = seed[p]['opposite_Pf_other_classes']['pct_N_run_ge3']
    dumbbell(ax, i, b, a, GREY, BLUE)
    ax.text(max(a, b) + 3, i, f'{a:.0f} vs {b:.0f}%', va='center', fontsize=6.2, color=MUTED)
ax.set_yticks(range(3))
ax.set_yticklabels(['P. berghei', 'P. knowlesi', 'P. vivax'], fontstyle='italic', fontsize=6.6)
ax.tick_params(axis='y', length=0)
ax.spines['left'].set_visible(False)
ax.invert_yaxis()
ax.set_xlim(0, 75)
ax.set_ylim(2.5, -0.8)
ax.set_xlabel('Ortholog segments with an N run ≥ 3 (%)')
xgrid(ax)
ax.text(55, -0.45, 'opposite Pf poly-N', color=BLUE, fontsize=6, fontweight='bold', ha='center')
ax.text(36, 0.2, 'other Pf LCRs', color=MUTED, fontsize=5.8, ha='center', va='top')
panel(ax, 'B', 'They grow from asparagine seeds')
ax = fig.add_subplot(gs[2])
nk = ct3['test1_N_vs_K']
names = [('Pf', 'P. falciparum'), ('Pb', 'P. berghei'), ('Pk', 'P. knowlesi'), ('Pv', 'P. vivax')]
for i, (s, nm) in enumerate(names):
    bg = 2 ** nk[s]['segments_without_tracts']['log2_obs_over_expected']
    tr = 2 ** nk[s]['tracts_in_segments']['log2_obs_over_expected']
    dumbbell(ax, i, bg, tr, GREY, BLUE)
ax.axvline(1, color=MUTED, lw=0.7, ls=(0, (2, 2)))
ax.text(0.97, -0.38, 'mutation\nbias alone', fontsize=5.6, color=MUTED, ha='right', va='center')
ax.set_xscale('log', base=2)
ax.set_xticks([0.8, 1, 1.5, 2, 2.6])
ax.set_xticklabels(['0.8', '1', '1.5', '2', '2.6×'])
ax.set_yticks(range(4))
ax.set_yticklabels([n for _, n in names], fontstyle='italic', fontsize=6.6)
ax.tick_params(axis='y', length=0)
ax.spines['left'].set_visible(False)
ax.invert_yaxis()
ax.set_ylim(4.0, -0.6)
ax.set_xlabel('Asn:Lys relative to expectation')
xgrid(ax)
ax.text(2.45, -0.35, 'tracts', color=BLUE, fontsize=6, fontweight='bold', ha='center')
ax.text(1.33, 0.42, 'background', color=MUTED, fontsize=5.8, ha='center')
panel(ax, 'C', 'Asn over Lys is shared by AT-rich species')
save(fig, 'Figure3')

# =========================== Figure 4 ===========================
dat = J('step6/pf_tract_dating.json')
datl = J('step6/pf_tract_dating_lenient.json')
cd = J('step8/codon_decay.json')
cw = J('step8/codon_within.json')
AGES = ['Pf only', 'Pf+Pprf', '+Pr', '+Pbilc', 'clade B', 'all Laverania']
keys = ['Pf only', 'Pf+Pprf', '+Pr', '+Pbilc', 'clade B (+Pblac)', 'all Laverania (+clade A)']
fig = plt.figure(figsize=(7.4, 2.85))
gs = fig.add_gridspec(1, 3, width_ratios=[1.05, 1.15, 0.85], wspace=0.55)
ax = fig.add_subplot(gs[0])
st = [dat['origin_counts'][k] for k in keys]
le = [datl['origin_counts'][k] for k in keys]
x = np.arange(6)
ax.bar(x, le, 0.72, color='none', edgecolor=AGE, linewidth=0.9, linestyle=(0, (2, 1.5)))
ax.bar(x, st, 0.72, color=AGE)
for i in range(6):
    ax.text(i, max(st[i], le[i]) + 30, f'{st[i]:,}', ha='center', fontsize=6, color=INK)
ax.set_xticks(x)
ax.set_xticklabels(AGES, rotation=40, ha='right', fontsize=6.4)
ax.set_ylabel('Pf poly-N tracts')
ax.set_ylim(0, 2000)
ygrid(ax)
ax.text(5, 1880, 'dashed: lenient', ha='center', fontsize=5.8, color=MUTED)
ax.text(2.2, 1300, 'young  ←      → old', fontsize=6, color=MUTED, ha='center')
pill(ax, 0.40, 0.48, '63–81% predate Pf', DEEP)
panel(ax, 'A', 'Most tracts are Laverania-wide')
ax = fig.add_subplot(gs[1])
cl = cd['classes']
obs = np.array([c['observed'] for c in cl])
se = np.array([c['se'] for c in cl])
pred = np.array([c['predicted_born_pure'] for c in cl])
ax.fill_between(x, obs, pred, color=ORANGE, alpha=0.07, lw=0)
ax.plot(x, pred, '-', color=ORANGE, lw=1.6)
ax.scatter(x, pred, s=16, color=ORANGE, marker='s', zorder=3, edgecolor='white', linewidth=0.5)
ax.errorbar(x, obs, yerr=1.96 * se, fmt='none', ecolor=BLUE, elinewidth=1.0, capsize=0)
ax.scatter(x, obs, s=30, color=AGE, zorder=4, edgecolor=DEEP, linewidth=0.6)
ax.axhline(cd['p_eq'], color=MUTED, lw=0.7, ls=(0, (2, 2)))
ax.text(-0.2, cd['p_eq'] - 0.003, 'genome equilibrium', fontsize=5.6, color=MUTED, ha='left', va='top')
ax.text(0.2, 1.004, 'neutral point mutation', color=ORANGE, fontsize=6.2, fontweight='bold')
ax.text(2.0, 0.872, 'observed', color=BLUE, fontsize=6.2, fontweight='bold')
ax.text(3.2, 0.935, 'gap = faster\nthan neutral', color=ORANGE, fontsize=5.8, ha='center')
ax.set_xticks(x)
ax.set_xticklabels(AGES, rotation=40, ha='right', fontsize=6.4)
ax.set_ylabel('AAT share of Asn codons')
ax.set_ylim(0.83, 1.012)
ygrid(ax)
panel(ax, 'B', 'Purity decays faster than neutral')
ax = fig.add_subplot(gs[2])
for i, (k, nm) in enumerate([('Pprf', 'P. praefalciparum'), ('Pr', 'P. reichenowi')]):
    b, w = 100 * cw[k]['background_D'], 100 * cw[k]['within_tract_D']
    dumbbell(ax, i, b, w, GREY, BLUE)
    ax.text(w + 0.5, i, f'{w / b:.1f}×', va='center', fontsize=7, color=BLUE, fontweight='bold')
ax.set_yticks([0, 1])
ax.set_yticklabels(['Pf vs\nP. praefalciparum', 'Pf vs\nP. reichenowi'], fontsize=6.4)
ax.tick_params(axis='y', length=0)
ax.spines['left'].set_visible(False)
ax.invert_yaxis()
ax.set_ylim(1.6, -0.7)
ax.set_xlim(0, 11)
ax.set_xlabel('AAT↔AAC divergence (%)')
xgrid(ax)
ax.text(0.5, 0.3, 'genome', color=MUTED, fontsize=6, ha='center')
ax.text(2.6, -0.3, 'inside tracts', color=BLUE, fontsize=6, fontweight='bold')
panel(ax, 'C', '4–5× faster inside tracts')
save(fig, 'Figure4')

# =========================== Figure 5 ===========================
pop = J('step7/pop_length_variation.json')
TT = J('step10/pf8_young_tract_segregation.json')
CC = J('step10/pf8_control_indels.json')
fig = plt.figure(figsize=(7.4, 2.8))
gs = fig.add_gridspec(1, 3, width_ratios=[1.25, 1, 1], wspace=0.55)
ax = fig.add_subplot(gs[0])
lm = pop['length_matched_pct_any5pct_and_mean_het']
bins = ['<18 aa', '18-24', '25-34', '35-49', '>=50']
cols = [AGE[1], AGE[2], AGE[4], AGE[5]]
labels = ['Pf only', 'Pf+Pprf', 'middle ages', 'all Laverania']
x = np.arange(5)
for j in range(4):
    v = [float(lm[b][j].split(': ')[1].split('%')[0]) for b in bins]
    lw = 2.4 if j == 0 else 1.4
    ax.plot(x, v, '-o', color=cols[j], lw=lw, ms=3.6 if j == 0 else 2.8, mec='white', mew=0.5, zorder=4 - j)
    ax.text(4.15, [63, 55.5, 49, 42][j], labels[j],
            color=cols[j], fontsize=6.0, fontweight='bold' if j == 0 else 'normal', va='center')
ax.set_xticks(x)
ax.set_xticklabels(['<18', '18–24', '25–34', '35–49', '≥50'], fontsize=6.2)
ax.set_xlim(-0.25, 5.5)
ax.set_ylim(10, 80)
ax.set_xlabel('Tract length (aa)')
ax.set_ylabel('Tracts with a length allele ≥ 5% (%)')
ygrid(ax)
pill(ax, 0.32, 0.07, 'OR 3.1, p = 7×10⁻⁶', AGE[3])
panel(ax, 'A', 'Youngest tracts vary most in Pf8')


def both(s):
    af = max(s['popf'].get(p) or 0 for p in ['AF-W', 'AF-C', 'AF-E'])
    a = max(s['popf'].get(p) or 0 for p in ['AS-SE-W', 'AS-SE-E'])
    return af >= 0.01 and a >= 0.01


ax = fig.add_subplot(gs[1])
ft = [s['fst'] for s in TT['sites'] if s['fst'] is not None]
fc = [s['fst'] for s in CC['sites'] if s['fst'] is not None]
fs = TT['snp_fst'] + CC['snp_fst']
rng = np.random.default_rng(1)
for i, (d, col) in enumerate([(ft, BLUE), (fc, ORANGE), (fs, GREY)]):
    med = np.median(d)
    d = np.clip(np.array(d), -0.2, 0.9)
    vp = ax.violinplot([d], positions=[i], widths=0.8, showextrema=False)
    for b in vp['bodies']:
        b.set_facecolor(col)
        b.set_alpha(0.22)
        b.set_edgecolor('none')
    ax.scatter(i + rng.uniform(-0.12, 0.12, len(d)), d, s=2.2, color=col, alpha=0.55, lw=0)
    ax.plot([i - 0.22, i + 0.22], [med] * 2, color=INK, lw=1.4)
    ax.text(i, 0.97, f'{med:.3f}', ha='center', fontsize=6.2, color=col if col != GREY else MUTED, fontweight='bold')
ax.text(-0.45, 1.03, 'median', fontsize=5.6, color=MUTED)
ax.set_xticks([0, 1, 2])
ax.set_xticklabels(['young\ntracts', 'non-repeat\nindels', 'flanking\nSNPs'], fontsize=6.2)
ax.set_ylim(-0.25, 1.06)
ax.set_ylabel('FST across 8 populations')
ygrid(ax)
panel(ax, 'B', 'Less structured than other variants')
ax = fig.add_subplot(gs[2])


def ldfrac(S):
    L = [s['ld_ALL'] for s in S if s.get('ld_ALL') and s['ld_ALL']['max_r2'] is not None and s['ld_ALL']['ctrl_median_max_r2'] is not None]
    return 100 * sum(l['max_r2'] < l['ctrl_median_max_r2'] for l in L) / len(L)


def ladder(S):
    return 100 * np.mean([(s.get('nlen') or 0) >= 3 for s in S])


tl = [100 * np.mean([both(s) for s in TT['sites']]), ldfrac(TT['sites'])]
cl2 = [100 * np.mean([both(s) for s in CC['sites']]), ldfrac(CC['sites'])]
names = ['shared Africa +\nSoutheast Asia', 'weaker LD than\nmatched SNPs']
pv = ['p = 1×10⁻⁵', 'p = 6×10⁻⁴']
for i in range(2):
    dumbbell(ax, i, cl2[i], tl[i], ORANGE, BLUE)
    ax.text(tl[i] + 3, i, f'{tl[i]:.0f}%', va='center', fontsize=6.6, color=BLUE, fontweight='bold')
    ax.text(cl2[i] - 3, i, f'{cl2[i]:.0f}%', va='center', ha='right', fontsize=6.6, color=ORANGE, fontweight='bold')
    ax.text(57, i + 0.3, pv[i], fontsize=5.8, color=MUTED, ha='center')
ax.axvline(50, color=MUTED, lw=0.6, ls=(0, (2, 2)))
ax.set_yticks([0, 1])
ax.set_yticklabels(names, fontsize=6.4)
ax.tick_params(axis='y', length=0)
ax.spines['left'].set_visible(False)
ax.invert_yaxis()
ax.set_ylim(1.6, -0.75)
ax.set_xlim(15, 95)
ax.set_xlabel('Sites (%)')
xgrid(ax)
ax.text(78, -0.5, 'young tracts', color=BLUE, fontsize=6, fontweight='bold', ha='center')
ax.text(33, -0.5, 'control indels', color=ORANGE, fontsize=6, fontweight='bold', ha='center')
panel(ax, 'C', 'Length alleles arise again and again')
save(fig, 'Figure5')

# =========================== Figure 6 ===========================
k = J('step7/constraint_vs_age.json')
fig = plt.figure(figsize=(5.4, 2.6))
gs = fig.add_gridspec(1, 2, width_ratios=[1.3, 1], wspace=0.6)
ax = fig.add_subplot(gs[0])
gc = k['Q1_gene_classes']
lab = [('no tract', 'no_tract'), ('young tracts only', 'young_only'), ('young + old', 'mixed'), ('old tracts only', 'old_only')]
for i, (l, g) in enumerate(lab):
    a, b = gc[g]['pct_essential_MIS'], gc[g]['pb_pct_essential']
    ax.scatter(a, i, s=34, color=BLUE, zorder=3, edgecolor='white', linewidth=0.6)
    ax.scatter(b, i, s=34, color=AQUA, marker='D', zorder=3, edgecolor='white', linewidth=0.6)
    ax.plot([min(a, b), max(a, b)], [i, i], color=LGREY, lw=3.5, zorder=1, solid_capstyle='round')
ax.set_yticks(range(4))
ax.set_yticklabels([l for l, _ in lab], fontsize=6.6)
ax.tick_params(axis='y', length=0)
ax.spines['left'].set_visible(False)
ax.invert_yaxis()
ax.set_xlim(30, 60)
ax.set_ylim(4.3, -0.9)
ax.set_xlabel('Essential genes (%)')
xgrid(ax)
ax.text(38, -0.6, 'P. falciparum', color=BLUE, fontsize=6, fontweight='bold', fontstyle='italic', ha='center')
ax.text(53, -0.6, 'P. berghei ortholog', color=AQUA, fontsize=6, fontweight='bold', fontstyle='italic', ha='center')
pill(ax, 0.5, 0.06, 'old vs young: OR 0.84, p = 0.52', GREY)
panel(ax, 'A', 'Old tracts are not in essential genes')
ax = fig.add_subplot(gs[1])
q3 = k['Q3_old_tracts_shared_outside_laverania']
v = [q3['essential']['pct_shared_outside'], q3['dispensable']['pct_shared_outside']]
ax.bar([0, 1], v, 0.6, color=[GREY, BLUE])
for i in range(2):
    ax.text(i, v[i] + 1.2, f'{v[i]:.0f}%', ha='center', fontsize=7, fontweight='bold', color=[MUTED, BLUE][i])
ax.plot([0, 0, 1, 1], [52, 54, 54, 52], color=INK, lw=0.6)
ax.text(0.5, 55.5, 'p = 2×10⁻⁹', ha='center', fontsize=6.2)
ax.set_xticks([0, 1])
ax.set_xticklabels([f"essential\n(n = {q3['essential']['n_old']})", f"dispensable\n(n = {q3['dispensable']['n_old']})"], fontsize=6.4)
ax.set_ylim(0, 62)
ax.set_ylabel('Old tracts also poly-N\noutside Laverania (%)')
ygrid(ax)
panel(ax, 'B', 'Longest-lived in dispensable genes')
save(fig, 'Figure6')
print('done', sorted(os.listdir(OUT)))
