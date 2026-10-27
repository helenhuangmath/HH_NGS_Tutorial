#!/usr/bin/env python3
"""
Generate the example figures used in the tutorials (tutorials/figures/).

All data are SIMULATED. They only show what each analysis step produces
and how to lay the plots out in Nature style; they are not real results.

Style follows Nature figure guidelines:
  - sans-serif (Arial/Helvetica), 5-7 pt text
  - 0.5 pt axis lines, outward ticks, no top/right spines
  - single column 89 mm, double column 183 mm
  - bold lowercase panel labels (a, b, c ...)
  - colour-blind-safe palette (Okabe-Ito)

Usage:
    python3 scripts/make_example_figures.py
"""

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap, ListedColormap
from matplotlib.patches import Circle, Rectangle, FancyArrow
from scipy import stats

OUT = Path(__file__).resolve().parent.parent / "tutorials" / "figures"
OUT.mkdir(parents=True, exist_ok=True)

rng = np.random.default_rng(2023)

MM = 1 / 25.4
SINGLE = 89 * MM
DOUBLE = 183 * MM

# Okabe-Ito
OI = dict(black="#000000", orange="#E69F00", sky="#56B4E9", green="#009E73",
          yellow="#F0E442", blue="#0072B2", red="#D55E00", purple="#CC79A7",
          grey="#8C8C8C")
COND = {"TN": "#8C8C8C", "TMEM": OI["blue"], "TEX": OI["red"]}
UP, DOWN, NS = OI["red"], OI["blue"], "#C8C8C8"
RDBU = LinearSegmentedColormap.from_list(
    "rdbu", ["#2166AC", "#67A9CF", "#F7F7F7", "#EF8A62", "#B2182B"])

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 6,
    "axes.titlesize": 7,
    "axes.labelsize": 6,
    "xtick.labelsize": 5.5,
    "ytick.labelsize": 5.5,
    "legend.fontsize": 5.5,
    "axes.linewidth": 0.5,
    "xtick.major.width": 0.5,
    "ytick.major.width": 0.5,
    "xtick.major.size": 2,
    "ytick.major.size": 2,
    "xtick.direction": "out",
    "ytick.direction": "out",
    "axes.spines.top": False,
    "axes.spines.right": False,
    "lines.linewidth": 0.75,
    "legend.frameon": False,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.02,
    "pdf.fonttype": 42,
    "svg.fonttype": "none",
})


def label(ax, letter, x=-0.18, y=1.06):
    ax.text(x, y, letter, transform=ax.transAxes, fontsize=8,
            fontweight="bold", va="bottom", ha="left")


def save(fig, name):
    fig.savefig(OUT / f"{name}.png")
    fig.savefig(OUT / f"{name}.pdf")
    plt.close(fig)
    print("wrote", OUT / f"{name}.png")


def pstar(p):
    return "****" if p < 1e-4 else "***" if p < 1e-3 else "**" if p < 1e-2 \
        else "*" if p < 0.05 else "n.s."


def repel_labels(ax, x, y, labels, avoid_xy=None, fixed=(), fontsize=5, italic=True,
                 iters=400, pad=2.2, offset=(14, 8)):
    """Place text labels near (x, y) without overlapping each other, the points in
    `avoid_xy`, the artists in `fixed`, or the axes frame (ggrepel-like)."""
    fig = ax.figure
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    T, inv = ax.transData, ax.transData.inverted()
    pts = T.transform(np.c_[x, y])
    av = np.vstack([pts, T.transform(avoid_xy)]) if avoid_xy is not None else pts
    abb = ax.get_window_extent(r)
    cx = abb.x0 + abb.width / 2
    pos = pts.copy()
    pos[:, 0] += np.where(pts[:, 0] > cx, -offset[0], offset[0])
    pos[:, 1] += offset[1]
    texts = [ax.text(*inv.transform(q), lab, fontsize=fontsize, ha="center", va="center",
                     style="italic" if italic else "normal", zorder=5)
             for q, lab in zip(pos, labels)]
    size = np.array([[t.get_window_extent(r).width + 2 * pad,
                      t.get_window_extent(r).height + 2 * pad] for t in texts])
    fixed_bb = [a.get_window_extent(r) for a in fixed]
    n = len(texts)
    for _ in range(iters):
        d = np.zeros_like(pos)
        for i in range(n):                                   # label–label
            for j in range(i + 1, n):
                dx, dy = pos[i] - pos[j]
                ox = (size[i, 0] + size[j, 0]) / 2 - abs(dx)
                oy = (size[i, 1] + size[j, 1]) / 2 - abs(dy)
                if ox > 0 and oy > 0:
                    if ox < oy:
                        sgn = np.sign(dx) or 1; d[i, 0] += sgn * ox / 2; d[j, 0] -= sgn * ox / 2
                    else:
                        sgn = np.sign(dy) or 1; d[i, 1] += sgn * oy / 2; d[j, 1] -= sgn * oy / 2
            for bb in fixed_bb:                              # label–fixed text
                fx, fy = (bb.x0 + bb.x1) / 2, (bb.y0 + bb.y1) / 2
                dx, dy = pos[i, 0] - fx, pos[i, 1] - fy
                ox = (size[i, 0] + bb.width) / 2 - abs(dx)
                oy = (size[i, 1] + bb.height) / 2 - abs(dy)
                if ox > 0 and oy > 0:
                    d[i, 1] += (np.sign(dy) or -1) * oy
            dx, dy = pos[i, 0] - av[:, 0], pos[i, 1] - av[:, 1]   # label–points
            ox = size[i, 0] / 2 + 1.5 - np.abs(dx)
            oy = size[i, 1] / 2 + 1.5 - np.abs(dy)
            hit = (ox > 0) & (oy > 0)
            if hit.any():
                k = np.argmax(np.where(hit, np.minimum(ox, oy), -1))
                if ox[k] < oy[k]:
                    d[i, 0] += (np.sign(dx[k]) or 1) * ox[k] * 0.6
                else:
                    d[i, 1] += (np.sign(dy[k]) or 1) * oy[k] * 0.6
        pos += d
        pos[:, 0] = np.clip(pos[:, 0], abb.x0 + size[:, 0] / 2 + 2, abb.x1 - size[:, 0] / 2 - 2)
        pos[:, 1] = np.clip(pos[:, 1], abb.y0 + size[:, 1] / 2 + 2, abb.y1 - size[:, 1] / 2 - 2)
        if not np.abs(d).max() > 0.05:
            break
    for t, q, pnt, sz in zip(texts, pos, pts, size):
        t.set_position(inv.transform(q))
        near = np.clip(pnt, q - sz / 2 + pad, q + sz / 2 - pad)  # closest point on label box
        if np.hypot(*(near - pnt)) > 3:
            ax.plot(*inv.transform(np.vstack([pnt, near])).T, color="#555555", lw=0.3, zorder=4)
    return texts


SAMPLES = [f"{c}_Rep{r}" for c in COND for r in (1, 2)]


# --------------------------------------------------------------------------
# Fig. 1  Read QC: alignment summary, fragment sizes, FRiP
# --------------------------------------------------------------------------
def fig_qc():
    fig, axes = plt.subplots(1, 3, figsize=(DOUBLE, 1.9),
                             gridspec_kw=dict(wspace=0.45, width_ratios=[1.2, 1.1, 0.9]))

    ax = axes[0]
    total = rng.uniform(14, 26, len(SAMPLES))
    uniq = total * rng.uniform(0.70, 0.85, len(SAMPLES))
    dup = total * rng.uniform(0.05, 0.15, len(SAMPLES))
    unal = total - uniq - dup
    y = np.arange(len(SAMPLES))
    ax.barh(y, uniq, color=OI["blue"], height=0.7, label="Unique, deduplicated")
    ax.barh(y, dup, left=uniq, color=OI["sky"], height=0.7, label="PCR duplicates")
    ax.barh(y, unal, left=uniq + dup, color="#D9D9D9", height=0.7, label="Unaligned")
    ax.set_yticks(y, SAMPLES)
    ax.invert_yaxis()
    ax.set_xlabel("Read pairs (millions)")
    ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=3,
              handlelength=0.8, columnspacing=0.8)
    label(ax, "a", x=-0.45)

    ax = axes[1]
    x = np.arange(20, 800)
    def frag(mono, di, noise):
        d = (mono * stats.norm.pdf(x, 165, 25) + di * stats.norm.pdf(x, 330, 40)
             + 0.08 * stats.norm.pdf(x, 500, 60) + noise * stats.norm.pdf(x, 80, 25))
        return d / d.sum() * 100
    for name, args, col in [("H3K27ac", (1, 0.35, 0.10), OI["red"]),
                            ("H3K27me3", (1, 0.6, 0.05), OI["blue"]),
                            ("IgG", (0.6, 0.2, 0.6), "#8C8C8C")]:
        ax.plot(x, frag(*args), color=col, label=name)
    ax.set_xlabel("Fragment length (bp)")
    ax.set_ylabel("Fragments (%)")
    ax.set_xlim(0, 800)
    ax.legend(loc="upper right", handlelength=1.2)
    for pos, yv, txt in [(165, 1.06, "mono-"), (330, 0.37, "di-nucleosome")]:
        ax.text(pos + 12, yv, txt, ha="left", fontsize=5, color="#555555")
    label(ax, "b")

    ax = axes[2]
    frip = {c: rng.uniform(lo, hi, 2) for c, (lo, hi) in
            zip(COND, [(0.32, 0.4), (0.28, 0.36), (0.25, 0.33)])}
    frip["IgG"] = rng.uniform(0.01, 0.03, 2)
    for i, (c, v) in enumerate(frip.items()):
        col = COND.get(c, "#D9D9D9")
        ax.bar(i, v.mean(), color=col, width=0.65, alpha=0.85)
        ax.scatter(np.full(2, i) + [-0.12, 0.12], v, s=6, color="black", zorder=3, lw=0)
    ax.set_xticks(range(len(frip)), list(frip), rotation=0)
    ax.set_ylabel("FRiP")
    ax.set_ylim(0, 0.5)
    ax.axhline(0.1, ls="--", lw=0.5, color="#555555")
    ax.text(3.4, 0.11, "0.1", fontsize=5, ha="right", color="#555555")
    label(ax, "c", x=-0.3)
    save(fig, "fig1_qc")


# --------------------------------------------------------------------------
# Fig. 2  Peaks: counts, widths, replicate overlap
# --------------------------------------------------------------------------
def fig_peaks():
    fig, axes = plt.subplots(1, 3, figsize=(DOUBLE, 1.8),
                             gridspec_kw=dict(wspace=0.45))
    ax = axes[0]
    marks = ["H3K27ac", "H3K4me3", "H3K27me3", "CTCF"]
    base = np.array([42, 21, 18, 35]) * 1000
    w = 0.26
    for j, c in enumerate(COND):
        vals = base * rng.uniform(0.8, 1.2, 4)
        ax.bar(np.arange(4) + (j - 1) * w, vals / 1000, width=w, color=COND[c], label=c)
    ax.set_xticks(range(4), marks, rotation=30, ha="right")
    ax.set_ylabel("Peaks (×10$^3$)")
    ax.legend(ncol=3, loc="lower center", bbox_to_anchor=(0.5, 1.0),
              handlelength=0.8, columnspacing=0.8)
    label(ax, "a")

    ax = axes[1]
    bins = np.logspace(2, 5, 45)
    for m, (mu, sd), col in [("H3K4me3", (np.log(1200), 0.5), OI["green"]),
                             ("H3K27ac", (np.log(1800), 0.7), OI["red"]),
                             ("H3K27me3", (np.log(6000), 1.0), OI["blue"])]:
        ax.hist(rng.lognormal(mu, sd, 20000), bins=bins, histtype="step",
                color=col, lw=0.8, density=True, label=m)
    ax.set_xscale("log")
    ax.set_xlabel("Peak width (bp)")
    ax.set_ylabel("Density")
    ax.legend(loc="upper right")
    label(ax, "b")

    ax = axes[2]
    ax.set_aspect("equal")
    ax.add_patch(Circle((0.38, 0.5), 0.3, fc=OI["blue"], alpha=0.35, ec=OI["blue"], lw=0.6))
    ax.add_patch(Circle((0.66, 0.5), 0.27, fc=OI["red"], alpha=0.35, ec=OI["red"], lw=0.6))
    for (x, yv, t) in [(0.22, 0.5, "4,812"), (0.52, 0.5, "31,406"), (0.81, 0.5, "2,957")]:
        ax.text(x, yv, t, ha="center", va="center", fontsize=6)
    ax.text(0.30, 0.86, "Rep1", ha="center", color=OI["blue"])
    ax.text(0.74, 0.83, "Rep2", ha="center", color=OI["red"])
    ax.text(0.52, 0.08, "Consensus = overlap (≥2 reps)", ha="center", fontsize=5)
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    ax.axis("off")
    label(ax, "c", x=0.0)
    save(fig, "fig2_peaks")


# --------------------------------------------------------------------------
# Fig. 3  deepTools-style profile + heatmap
# --------------------------------------------------------------------------
def fig_profile():
    n, nb = 1500, 80
    pos = np.linspace(-1, 1, nb)  # kb
    strength = np.sort(rng.gamma(2, 1.5, n))[::-1]
    width = rng.uniform(0.15, 0.4, n)
    base = strength[:, None] * np.exp(-(pos[None, :] ** 2) / (2 * width[:, None] ** 2))
    scale = {"TN": 0.7, "TMEM": 1.0, "TEX": 1.35}
    mats = {c: base * s + rng.gamma(1, 0.25, (n, nb)) for c, s in scale.items()}

    fig = plt.figure(figsize=(SINGLE * 1.35, 3.2))
    gs = fig.add_gridspec(2, 4, height_ratios=[1, 3.2], width_ratios=[1, 1, 1, 0.07],
                          hspace=0.08, wspace=0.22)
    vmax = 10
    for j, c in enumerate(COND):
        axp = fig.add_subplot(gs[0, j])
        axp.plot(pos, mats[c].mean(0), color=COND[c], lw=1)
        axp.set_ylim(0, 6)
        axp.set_xticks([])
        axp.set_title(c, color=COND[c], fontsize=7, pad=2)
        if j == 0:
            axp.set_ylabel("Mean CPM")
            label(axp, "a", x=-0.55, y=1.1)
        else:
            axp.set_yticklabels([])
        axh = fig.add_subplot(gs[1, j])
        im = axh.imshow(mats[c], aspect="auto", cmap="Reds", vmin=0, vmax=vmax,
                        interpolation="nearest", extent=[-1, 1, n, 0])
        axh.set_xticks([-1, 0, 1], ["−1", "0", "+1"])
        if j == 1:
            axh.set_xlabel("Distance from peak centre (kb)")
        axh.set_yticks([])
        if j == 0:
            axh.set_ylabel(f"H3K27ac peaks (n = {n:,})")
            label(axh, "b", x=-0.55, y=1.0)
        for s in axh.spines.values():
            s.set_visible(True)
    cax = fig.add_subplot(gs[1, 3])
    cb = fig.colorbar(im, cax=cax)
    cb.set_label("CPM")
    cb.outline.set_linewidth(0.5)
    save(fig, "fig3_profile_heatmap")


# --------------------------------------------------------------------------
# Fig. 4  Differential analysis: PCA, correlation, MA, volcano
# --------------------------------------------------------------------------
def simulate_de(n=20000):
    base = rng.lognormal(4, 1.6, n)
    lfc = rng.normal(0, 0.25, n)
    de = rng.random(n) < 0.08
    lfc[de] += rng.choice([-1, 1], de.sum()) * rng.gamma(2.5, 0.7, de.sum())
    se = 0.15 + 2.2 / np.sqrt(base)
    z = lfc / se
    p = 2 * stats.norm.sf(np.abs(z))
    padj = np.minimum(1, p * n / (stats.rankdata(p)))
    return base, lfc, padj


def fig_diff():
    fig, axes = plt.subplots(2, 2, figsize=(DOUBLE * 0.78, 3.7),
                             gridspec_kw=dict(wspace=0.75, hspace=0.6))
    axes = axes.ravel()
    # a PCA
    ax = axes[0]
    centers = {"TN": (-14, 2), "TMEM": (3, -6), "TEX": (12, 5)}
    for c, (cx, cy) in centers.items():
        pts = np.array([cx, cy]) + np.array([[-1.6, 1.1], [1.4, -1.0]]) + rng.normal(0, 0.4, (2, 2))
        ax.scatter(pts[:, 0], pts[:, 1], s=14, color=COND[c], label=c, lw=0.3, ec="white")
    ax.axhline(0, color="#D9D9D9", lw=0.4, zorder=0)
    ax.axvline(0, color="#D9D9D9", lw=0.4, zorder=0)
    ax.set_xlabel("PC1 (61.4%)")
    ax.set_ylabel("PC2 (18.2%)")
    ax.set_xlim(-20, 20); ax.set_ylim(-10, 10)
    ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=3,
              handletextpad=0.1, columnspacing=0.6)
    label(ax, "a", x=-0.25)

    # b correlation
    ax = axes[1]
    grp = np.repeat([0, 1, 2], 2)
    r = np.where(grp[:, None] == grp[None, :], 0.97, 0.86)
    r[(grp[:, None] - grp[None, :]) ** 2 == 4] = 0.78
    r = r + rng.normal(0, 0.008, r.shape)
    r = (r + r.T) / 2
    np.fill_diagonal(r, 1)
    im = ax.imshow(r, cmap="Blues", vmin=0.75, vmax=1)
    ax.set_xticks(range(6), SAMPLES, rotation=90)
    ax.set_yticks(range(6), SAMPLES)
    ax.tick_params(length=0)
    for s in ax.spines.values():
        s.set_visible(False)
    cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03)
    cb.set_label("Pearson r")
    cb.outline.set_linewidth(0.5)
    label(ax, "b", x=-0.55)

    base, lfc, padj = simulate_de()
    sig = (padj < 0.05) & (np.abs(lfc) > np.log2(1.5))
    col = np.where(sig & (lfc > 0), UP, np.where(sig & (lfc < 0), DOWN, NS))
    order = np.argsort(sig)

    # c MA
    ax = axes[2]
    ax.scatter(np.log10(base[order]), lfc[order], c=col[order], s=0.6, lw=0, rasterized=True)
    ax.axhline(0, color="black", lw=0.4)
    ax.set_xlabel("log$_{10}$ mean normalized count")
    ax.set_ylabel("log$_2$ fold change (TEX/TN)")
    ax.set_ylim(-7, 7)
    label(ax, "c", x=-0.25)

    # d volcano
    ax = axes[3]
    ax.scatter(lfc[order], -np.log10(padj[order]), c=col[order], s=0.6, lw=0, rasterized=True)
    ax.axvline(np.log2(1.5), ls="--", color="#555555", lw=0.4)
    ax.axvline(-np.log2(1.5), ls="--", color="#555555", lw=0.4)
    ax.axhline(-np.log10(0.05), ls="--", color="#555555", lw=0.4)
    ax.set_xlabel("log$_2$ fold change (TEX/TN)")
    ax.set_ylabel("−log$_{10}$ adjusted P")
    ax.set_xlim(-7, 7)
    ax.set_ylim(-5, (-np.log10(padj)).max() * 1.25)
    nu, nd = (sig & (lfc > 0)).sum(), (sig & (lfc < 0)).sum()
    t_up = ax.text(0.98, 0.98, f"Up {nu:,}", color=UP, transform=ax.transAxes, ha="right", va="top")
    t_dn = ax.text(0.02, 0.98, f"Down {nd:,}", color=DOWN, transform=ax.transAxes, ha="left", va="top")
    ord_p = np.argsort(padj)
    top_up = [i for i in ord_p if lfc[i] > 0][:5]
    top_dn = [i for i in ord_p if lfc[i] < 0][:5]
    idx = top_up + top_dn
    names = ["Tox", "Pdcd1", "Havcr2", "Entpd1", "Cd160", "Tcf7", "Sell", "Il7r", "Ccr7", "Lef1"]
    repel_labels(ax, lfc[idx], -np.log10(padj[idx]), names,
                 avoid_xy=np.c_[lfc[sig], -np.log10(padj[sig])], fixed=(t_up, t_dn))
    label(ax, "d", x=-0.25)
    save(fig, "fig4_differential")


# --------------------------------------------------------------------------
# Fig. 5  Annotation + boxplot with Wilcoxon tests
# --------------------------------------------------------------------------
def fig_anno():
    fig, axes = plt.subplots(1, 3, figsize=(DOUBLE, 1.9),
                             gridspec_kw=dict(wspace=0.5, width_ratios=[1.35, 1, 0.9]))
    ax = axes[0]
    cats = ["Promoter", "5′ UTR", "3′ UTR", "Exon", "Intron", "Intergenic", "TTS"]
    cols = ["#0072B2", "#56B4E9", "#009E73", "#F0E442", "#E69F00", "#D55E00", "#CC79A7"]
    sets = {"All peaks": [28, 1, 2, 4, 33, 30, 2],
            "Up in TEX": [14, 1, 2, 3, 41, 37, 2],
            "Down in TEX": [36, 2, 2, 4, 28, 26, 2]}
    for i, (k, v) in enumerate(sets.items()):
        v = np.array(v) / sum(v) * 100
        left = 0
        for val, c, cat in zip(v, cols, cats):
            ax.barh(i, val, left=left, color=c, height=0.65, label=cat if i == 0 else None,
                    ec="white", lw=0.3)
            left += val
    ax.set_yticks(range(3), list(sets))
    ax.invert_yaxis()
    ax.set_xlabel("Peaks (%)")
    ax.set_xlim(0, 100)
    ax.legend(loc="lower center", bbox_to_anchor=(0.45, 1.0), ncol=4,
              handlelength=0.8, columnspacing=0.6)
    label(ax, "a", x=-0.38)

    ax = axes[1]
    d = np.concatenate([rng.normal(0, 0.8, 9000), rng.normal(0, 25, 21000)])
    d = d[np.abs(d) < 100]
    ax.hist(d, bins=np.linspace(-100, 100, 81), color=OI["blue"], alpha=0.85)
    ax.axvspan(-2, 0.5, color=OI["red"], alpha=0.2, lw=0)
    ax.set_xlabel("Distance to TSS (kb)")
    ax.set_ylabel("Peaks")
    ax.text(3, ax.get_ylim()[1] * 0.92, "promoter\n(−2 kb, +0.5 kb)", fontsize=5,
            color=OI["red"], va="top")
    label(ax, "b")

    ax = axes[2]
    vals = [rng.normal(m, 1.1, 1500) for m in (5.0, 5.6, 6.4)]
    bp = ax.boxplot(vals, widths=0.6, patch_artist=True, showfliers=False,
                    medianprops=dict(color="black", lw=0.8),
                    whiskerprops=dict(lw=0.5), capprops=dict(lw=0.5),
                    boxprops=dict(lw=0.5))
    for patch, c in zip(bp["boxes"], COND):
        patch.set_facecolor(COND[c])
        patch.set_alpha(0.85)
    ax.scatter([1, 2, 3], [v.mean() for v in vals], marker="^", s=8, color="white",
               zorder=3, lw=0.3, ec="black")
    ax.axhline(np.median(vals[0]), ls="--", lw=0.4, color="#555555")
    ax.set_xticks([1, 2, 3], list(COND))
    ax.set_ylabel("log$_2$(normalized count + 1)")
    ytop = max(np.percentile(v, 99) for v in vals)
    for k, (i, j) in enumerate([(0, 1), (1, 2), (0, 2)]):
        p = stats.mannwhitneyu(vals[i], vals[j]).pvalue
        h = ytop + 0.4 + k * 0.75
        ax.plot([i + 1, i + 1, j + 1, j + 1], [h - 0.15, h, h, h - 0.15], lw=0.5, color="black")
        ax.text((i + j) / 2 + 1, h + 0.02, pstar(p), ha="center", va="bottom", fontsize=5.5)
    ax.set_ylim(top=ytop + 2.6)
    label(ax, "c", x=-0.4)
    save(fig, "fig5_annotation_boxplot")


# --------------------------------------------------------------------------
# Fig. 6  Motifs & pathway enrichment
# --------------------------------------------------------------------------
def fig_motif():
    fig, axes = plt.subplots(1, 3, figsize=(DOUBLE, 2.4),
                             gridspec_kw=dict(wspace=1.25, width_ratios=[0.8, 1.2, 1.1]))
    ax = axes[0]
    motifs = ["BATF", "AP-1", "JunB", "NFAT", "IRF4", "RUNX", "TCF7", "LEF1", "ETS1", "KLF2"]
    comps = ["TEX up", "TEX down"]
    fc = np.array([[3.1, 1.0], [2.8, 1.1], [2.6, 1.0], [2.1, 1.2], [1.9, 1.1],
                   [1.3, 1.6], [1.0, 2.7], [1.0, 2.4], [1.1, 1.9], [1.0, 1.8]])
    lp = fc * rng.uniform(9, 14, fc.shape) - 8
    lp = np.clip(lp, 1, None)
    yy = np.arange(len(motifs))
    cmap = plt.get_cmap("viridis")
    norm = mpl.colors.Normalize(0, lp.max())
    for j in range(2):
        ax.scatter(np.full(len(motifs), j), yy, s=fc[:, j] ** 2 * 6, c=lp[:, j],
                   cmap=cmap, norm=norm, lw=0.3, ec="black")
    ax.set_xticks([0, 1], comps, rotation=40, ha="right")
    ax.set_yticks(yy, motifs)
    ax.invert_yaxis()
    ax.set_xlim(-0.6, 1.6)
    sm = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
    cb = fig.colorbar(sm, ax=ax, fraction=0.08, pad=0.08, aspect=12)
    cb.set_label("−log$_{10}$ P")
    cb.outline.set_linewidth(0.5)
    for s in (1.5, 2.5, 3.0):
        ax.scatter([], [], s=s ** 2 * 6, color="#BBBBBB", ec="black", lw=0.3, label=f"{s}")
    ax.legend(title="Fold\nenrichment", loc="upper left", bbox_to_anchor=(1.55, 0.42),
              labelspacing=1.0, title_fontsize=5.5, borderaxespad=0)
    label(ax, "a", x=-0.9, y=1.02)

    ax = axes[1]
    n = 400
    lpv = rng.exponential(3, n)
    fcd = (1 + rng.exponential(0.25, n)) * rng.choice([-1, 1], n)
    hits = {"BATF": (35, 3.1), "AP-1": (31, 2.8), "NFAT": (18, 2.1), "TCF7": (29, -2.7),
            "LEF1": (24, -2.4), "ETS1": (15, -1.9)}
    for k, (a, b) in hits.items():
        lpv = np.append(lpv, a); fcd = np.append(fcd, b)
    c = np.where(fcd > 0, UP, DOWN)
    ax.scatter(lpv, fcd, s=3, c=c, alpha=0.6, lw=0)
    offs = {"BATF": (-2, 7), "AP-1": (-6, -9), "NFAT": (0, 7), "TCF7": (-2, -8),
            "LEF1": (-6, 7), "ETS1": (0, -8)}
    for k, (a, b) in hits.items():
        ax.annotate(k, (a, b), xytext=offs[k], textcoords="offset points",
                    fontsize=5, ha="center", va="center")
    ax.axhline(0, color="black", lw=0.4)
    ax.set_xlabel("−log$_{10}$ P")
    ax.set_ylabel("Signed fold enrichment")
    ax.set_ylim(-4.2, 4.2)
    ax.text(0.02, 0.98, "enriched in TEX-up peaks", transform=ax.transAxes, ha="left",
            va="top", color=UP, fontsize=5)
    ax.text(0.02, 0.02, "enriched in TEX-down peaks", transform=ax.transAxes, ha="left",
            va="bottom", color=DOWN, fontsize=5)
    label(ax, "b", x=-0.25, y=1.02)

    ax = axes[2]
    terms = ["T cell activation", "Regulation of cytokine production",
             "Response to interferon-γ", "Leukocyte differentiation",
             "PD-1 signalling", "Positive regulation of apoptosis",
             "Antigen processing via MHC I", "Cell–cell adhesion"]
    gr = np.sort(rng.uniform(0.04, 0.16, len(terms)))[::-1]
    cnt = (gr * 120).astype(int)
    padj = np.sort(10 ** -rng.uniform(2, 9, len(terms)))
    yy = np.arange(len(terms))
    sc = ax.scatter(gr, yy, s=cnt * 1.6, c=-np.log10(padj), cmap="magma_r", lw=0.3, ec="black")
    ax.set_yticks(yy, terms)
    ax.invert_yaxis()
    ax.set_xlabel("Gene ratio")
    ax.set_xlim(0.02, 0.18)
    cb = fig.colorbar(sc, ax=ax, fraction=0.07, pad=0.04, aspect=14)
    cb.set_label("−log$_{10}$ adj. P")
    cb.outline.set_linewidth(0.5)
    label(ax, "c", x=-1.3, y=1.02)
    save(fig, "fig6_motif_enrichment")


# --------------------------------------------------------------------------
# Fig. 7  RNA-seq k-means heatmap (ComplexHeatmap layout)
# --------------------------------------------------------------------------
def fig_rna_heatmap():
    k = 7
    sizes = rng.integers(80, 260, k)
    patterns = np.array([[1.2, -0.4, -0.8], [-1.0, 1.0, 0.0], [-0.8, -0.4, 1.2],
                         [1.0, 0.2, -1.2], [-1.2, 0.4, 0.8], [0.6, -1.2, 0.6],
                         [0.0, 1.1, -1.1]])
    rows, cl = [], []
    for i, (s, pat) in enumerate(zip(sizes, patterns)):
        m = np.repeat(pat, 4)[None, :] + rng.normal(0, 0.35, (s, 12))
        rows.append(m); cl += [i] * s
    Z = np.vstack(rows)
    Z = (Z - Z.mean(1, keepdims=True)) / Z.std(1, keepdims=True)
    cl = np.array(cl)
    tn, tmem, tex = Z[:, :4].mean(1), Z[:, 4:8].mean(1), Z[:, 8:].mean(1)
    diff = np.c_[tmem - tn, tex - tn, tex - tmem]
    de = np.where(np.abs(diff) < 0.9, 0, np.sign(diff))

    fig = plt.figure(figsize=(SINGLE, 3.4))
    gs = fig.add_gridspec(2, 4, height_ratios=[0.05, 1], width_ratios=[0.06, 1, 0.08, 0.32],
                          hspace=0.02, wspace=0.06, right=0.78)
    clcols = plt.get_cmap("BuPu")(np.linspace(0.3, 0.95, k))
    axc = fig.add_subplot(gs[1, 0])
    axc.imshow(cl[:, None], aspect="auto", cmap=ListedColormap(clcols), interpolation="nearest")
    axc.set_xticks([]); axc.set_yticks([])
    bounds = np.r_[0, np.cumsum(sizes)]
    for i in range(k):
        axc.text(-0.9, (bounds[i] + bounds[i + 1]) / 2, f"C{i + 1}", ha="right", va="center",
                 fontsize=5.5)
    for s in axc.spines.values():
        s.set_visible(False)

    axh = fig.add_subplot(gs[1, 1])
    im = axh.imshow(Z, aspect="auto", cmap=RDBU, vmin=-2, vmax=2, interpolation="nearest")
    for b in bounds[1:-1]:
        axh.axhline(b - 0.5, color="white", lw=1.2)
    axh.set_xticks([]); axh.set_yticks([])
    for s in axh.spines.values():
        s.set_visible(False)
    axt = fig.add_subplot(gs[0, 1])
    axt.imshow(np.repeat([0, 1, 2], 4)[None, :], aspect="auto",
               cmap=ListedColormap(list(COND.values())))
    axt.set_xticks([1.5, 5.5, 9.5], list(COND))
    axt.xaxis.tick_top(); axt.tick_params(length=0, pad=1)
    axt.set_yticks([])
    for s in axt.spines.values():
        s.set_visible(False)

    axd = fig.add_subplot(gs[1, 3])
    axd.imshow(de, aspect="auto", cmap=ListedColormap([DOWN, "white", UP]), vmin=-1, vmax=1,
               interpolation="nearest")
    axd.set_xticks([0, 1, 2], ["TMEM vs TN", "TEX vs TN", "TEX vs TMEM"], rotation=90)
    axd.set_yticks([])
    for b in bounds[1:-1]:
        axd.axhline(b - 0.5, color="white", lw=1.2)
    for s in axd.spines.values():
        s.set_visible(False)

    cax = fig.add_axes([0.83, 0.55, 0.025, 0.2])
    cb = fig.colorbar(im, cax=cax)
    cb.set_label("Row z-score")
    cb.outline.set_linewidth(0.5)
    lax = fig.add_axes([0.81, 0.3, 0.15, 0.15]); lax.axis("off")
    for c, t in [(UP, "Up"), (DOWN, "Down")]:
        lax.bar(0, 0, color=c, label=t)
    lax.legend(title="DE call", loc="center left", handlelength=0.8, title_fontsize=5.5)
    fig.text(0.5, 0.965, f"DEGs (FDR < 0.05, FC > 1.5; n = {len(Z):,})", ha="center", fontsize=7)
    save(fig, "fig7_rna_kmeans_heatmap")


# --------------------------------------------------------------------------
# Fig. 8  Genome-browser style tracks (Gviz layout)
# --------------------------------------------------------------------------
def fig_tracks():
    start, end = 1_000_000, 1_060_000
    x = np.linspace(start, end, 1200)
    peaks = [(1_008_000, 1.0), (1_021_000, 0.6), (1_034_500, 1.4), (1_052_000, 0.8)]
    gain = {"TN": [1, 0.4, 0.2, 1], "TMEM": [1, 0.8, 0.7, 1], "TEX": [1, 1.4, 1.6, 1]}
    tracks = [("H3K27ac " + c, COND[c], gain[c]) for c in COND] + [("IgG", "#8C8C8C", [0.05] * 4)]

    fig, axes = plt.subplots(len(tracks) + 2, 1, figsize=(DOUBLE * 0.75, 2.6), sharex=True,
                             gridspec_kw=dict(height_ratios=[1] * len(tracks) + [0.3, 0.6],
                                              hspace=0.15))
    ymax = 20
    for ax, (name, col, g) in zip(axes, tracks):
        y = rng.gamma(1, 0.25, x.size)
        for (p, h), gg in zip(peaks, g):
            y += 8 * h * gg * np.exp(-((x - p) / 1500) ** 2)
        ax.fill_between(x, 0, y, color=col, lw=0)
        ax.set_ylim(0, ymax)
        ax.set_yticks([])
        ax.spines["left"].set_visible(False)
        ax.spines["bottom"].set_visible(False)
        ax.tick_params(bottom=False)
        ax.text(start, ymax * 0.95, name, fontsize=5.5, va="top", color=col)
        ax.text(end, ymax * 0.95, f"[0–{ymax}]", fontsize=5, va="top", ha="right", color="#555555")
    ax = axes[-2]
    for p, _ in peaks[1:3]:
        ax.add_patch(Rectangle((p - 2500, 0.2), 5000, 0.6, color=OI["red"], lw=0))
    ax.set_ylim(0, 1); ax.set_yticks([]); ax.axis("off")
    ax.text(start, 0.5, "TEX-up DARs", fontsize=5.5, va="center", color=OI["red"])
    ax = axes[-1]
    ax.axis("off")
    ax.set_ylim(0, 1)
    ax.plot([1_006_000, 1_048_000], [0.5, 0.5], color="black", lw=0.6)
    for a, b in [(1_006_000, 1_008_500), (1_015_000, 1_016_000), (1_026_000, 1_027_500),
                 (1_040_000, 1_041_000), (1_046_000, 1_048_000)]:
        ax.add_patch(Rectangle((a, 0.3), b - a, 0.4, color="black", lw=0))
    for xx in np.arange(1_010_000, 1_046_000, 5000):
        ax.annotate("", (xx + 900, 0.5), (xx, 0.5),
                    arrowprops=dict(arrowstyle="-|>", lw=0.4, color="black", mutation_scale=4))
    ax.text(1_027_000, -0.05, "Gene1", style="italic", ha="center", va="top", fontsize=5.5)
    ax.plot([end - 10_000, end], [0.95, 0.95], color="black", lw=0.6)
    ax.text(end - 5_000, 1.0, "10 kb", ha="center", va="bottom", fontsize=5)
    axes[0].set_title("chr1:1,000,000–1,060,000", fontsize=6, loc="left", pad=2)
    save(fig, "fig8_genome_tracks")


# --------------------------------------------------------------------------
# Fig. 9  SON TSA-seq
# --------------------------------------------------------------------------
def fig_tsa():
    fig, axes = plt.subplots(1, 2, figsize=(DOUBLE, 1.8),
                             gridspec_kw=dict(width_ratios=[3, 1], wspace=0.3))
    ax = axes[0]
    pos = np.linspace(0, 50, 2500)
    sig = np.convolve(rng.normal(0, 1, pos.size), np.ones(60) / 60, mode="same") * 7
    sig += 1.4 * np.sin(pos / 3.1) + 0.8 * np.sin(pos / 1.3)
    ax.fill_between(pos, 0, sig, where=sig > 0, color=OI["red"], lw=0)
    ax.fill_between(pos, 0, sig, where=sig <= 0, color=OI["blue"], lw=0)
    ax.axhline(0, color="black", lw=0.4)
    ax.set_xlim(0, 50)
    ax.set_xlabel("chr2 position (Mb)")
    ax.set_ylabel("SON TSA-seq\nlog$_2$(pull-down/IgG)")
    ax.set_title("20-kb sliding window, 100-bp step", loc="right", fontsize=5.5, pad=2)
    label(ax, "a", x=-0.07)

    ax = axes[1]
    data = [np.log2(rng.lognormal(0.8 + 0.32 * d, 1.1, 600) + 1) for d in range(10)]
    bp = ax.boxplot(data, widths=0.6, patch_artist=True, showfliers=False,
                    medianprops=dict(color="black", lw=0.6), whiskerprops=dict(lw=0.4),
                    capprops=dict(lw=0.4), boxprops=dict(lw=0.4))
    cm = LinearSegmentedColormap.from_list("tsa", [OI["blue"], "#F7F7F7", OI["red"]])
    for k, b_ in enumerate(bp["boxes"]):
        b_.set_facecolor(cm(k / 9))
    ax.set_xticks([1, 5, 10], ["1\n(far)", "5", "10\n(near)"])
    ax.set_xlabel("SON TSA-seq decile")
    ax.set_ylabel("log$_2$(TPM + 1)")
    label(ax, "b", x=-0.35)
    save(fig, "fig9_tsa_seq")


# --------------------------------------------------------------------------
# Fig. 12  Hi-C: contact map, insulation, distance decay
# --------------------------------------------------------------------------
def fig_hic():
    fig = plt.figure(figsize=(DOUBLE, 2.5))
    gs = fig.add_gridspec(2, 4, height_ratios=[3, 1], width_ratios=[1.1, 0.04, 0.32, 1],
                          hspace=0.08, wspace=0.12)
    n = 120
    i, j = np.indices((n, n))
    d = np.abs(i - j) + 1
    M = 1 / d ** 1.1
    tads = [0, 25, 55, 80, 120]
    for a, b in zip(tads[:-1], tads[1:]):
        M[a:b, a:b] *= 3
    M += rng.gamma(1, 0.002, M.shape)
    M = (M + M.T) / 2
    ax = fig.add_subplot(gs[0, 0])
    im = ax.imshow(np.log10(M), cmap="YlOrRd", vmin=-2.6, vmax=0.4, interpolation="nearest",
                   aspect="auto")
    for a, b in zip(tads[:-1], tads[1:]):
        ax.plot([a, b, b], [a, a, b], color=OI["blue"], lw=0.6)
    ax.set_xticks([]); ax.set_yticks([0, 60, 119], ["10", "13", "16"])
    ax.set_ylabel("chr6 (Mb)")
    label(ax, "a", x=-0.18)
    cax = fig.add_subplot(gs[0, 1])
    cb = fig.colorbar(im, cax=cax); cb.set_label("log$_{10}$ ICE contacts"); cb.outline.set_linewidth(0.5)

    ax = fig.add_subplot(gs[1, 0])
    w = 5
    ins = np.array([np.log2(M[max(0, k - w):k, k + 1:k + 1 + w].mean() + 1e-9) if 0 < k < n - 1 else np.nan
                    for k in range(n)])
    ins = ins - np.nanmean(ins)
    ax.plot(np.arange(n), ins, color="black", lw=0.6)
    for b in tads[1:-1]:
        ax.axvline(b, color=OI["blue"], lw=0.5, ls="--")
    ax.set_xlim(0, n - 1)
    ax.set_xticks([0, 60, 119], ["10", "13", "16"])
    ax.set_xlabel("chr6 (Mb)")
    ax.set_ylabel("Insulation")
    ax.set_yticks([])

    ax = fig.add_subplot(gs[:, 3])
    s = np.logspace(4, 7.5, 60)
    for name, slope, col in [("TN", -1.05, COND["TN"]), ("TEX", -0.85, COND["TEX"])]:
        p = (s / 1e4) ** slope * 1e-2 * np.exp(rng.normal(0, 0.03, s.size))
        ax.loglog(s, p, color=col, label=name, lw=0.8)
    ax.loglog([1e5, 1e6], [2e-3, 2e-4], color="#555555", lw=0.5, ls="--")
    ax.text(3.5e5, 8e-4, "slope −1", fontsize=5, color="#555555", rotation=-28)
    ax.set_xlabel("Genomic distance s (bp)")
    ax.set_ylabel("Contact probability P(s)")
    ax.legend(loc="upper right")
    label(ax, "b", x=-0.28)
    save(fig, "fig12_hic")


# --------------------------------------------------------------------------
# Shared single-cell simulation
# --------------------------------------------------------------------------
CELLTYPES = ["Naive", "Tcm", "Tem", "Tpex", "Tex", "Cycling"]
CT_COL = dict(zip(CELLTYPES, ["#8C8C8C", OI["sky"], OI["blue"], OI["green"], OI["red"], OI["purple"]]))
CT_CENTER = {"Naive": (-7, 4), "Tcm": (-4, 7.5), "Tem": (1, 7), "Tpex": (0, -1),
             "Tex": (6, -3), "Cycling": (6.5, 4.5)}


def sim_cells(n=4000):
    w = np.array([0.18, 0.14, 0.2, 0.16, 0.24, 0.08])
    ct = rng.choice(len(CELLTYPES), n, p=w)
    xy = np.array([CT_CENTER[CELLTYPES[c]] for c in ct]) + rng.normal(0, 1.1, (n, 2))
    return ct, xy


def umap_ax(ax, xy, colors, title=None, s=1.2):
    perm = rng.permutation(len(xy))
    ax.scatter(xy[perm, 0], xy[perm, 1], s=s, c=np.asarray(colors)[perm], lw=0, rasterized=True)
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_xlabel("UMAP1", labelpad=1); ax.set_ylabel("UMAP2", labelpad=1)
    if title:
        ax.set_title(title, fontsize=6.5, pad=2)


# --------------------------------------------------------------------------
# Fig. 13  scRNA-seq
# --------------------------------------------------------------------------
def fig_scrna():
    ct, xy = sim_cells()
    n = len(ct)
    fig = plt.figure(figsize=(DOUBLE, 4.3))
    outer = fig.add_gridspec(2, 1, height_ratios=[1, 1.15], hspace=0.45)
    top = outer[0].subgridspec(1, 4, width_ratios=[1.7, 1.3, 0.75, 0.75], wspace=0.35)
    bot = outer[1].subgridspec(1, 2, width_ratios=[1.45, 0.9], wspace=0.55)

    # a QC violins
    qc = top[0].subgridspec(1, 3, wspace=0.7)
    nfeat = np.r_[rng.normal(2600, 650, n), rng.normal(600, 150, 300)]
    ncount = np.exp(np.log(nfeat) * 1.15 + rng.normal(0.6, 0.15, nfeat.size))
    mt = np.r_[rng.gamma(2.5, 1.4, n), rng.gamma(6, 3, 300)]
    for k, (v, nm, thr, logy) in enumerate([(nfeat, "Genes per cell", (200, 6000), False),
                                           (ncount, "UMIs per cell", None, True),
                                           (mt, "% mitochondrial", (10,), False)]):
        ax = fig.add_subplot(qc[k])
        vals = np.log10(v) if logy else v
        vp = ax.violinplot(vals, widths=0.8, showextrema=False)
        for b_ in vp["bodies"]:
            b_.set_facecolor("#BBBBBB"); b_.set_alpha(1); b_.set_edgecolor("none")
        ax.scatter(1 + rng.uniform(-0.25, 0.25, 600), rng.choice(vals, 600), s=0.3, c="black",
                   lw=0, alpha=0.5, rasterized=True)
        for t in (thr or ()):
            ax.axhline(t, color=OI["red"], lw=0.5, ls="--")
        ax.set_xticks([])
        ax.set_title(nm, fontsize=5.5, pad=2)
        if logy:
            ax.set_yticks([3, 4], ["10$^3$", "10$^4$"])
        if k == 0:
            label(ax, "a", x=-1.1)

    # b UMAP by cell type
    ax = fig.add_subplot(top[1])
    umap_ax(ax, xy, [CT_COL[CELLTYPES[c]] for c in ct], "Annotated clusters")
    for c, nm in enumerate(CELLTYPES):
        cx, cy = xy[ct == c].mean(0)
        ax.text(cx, cy, nm, fontsize=5.5, ha="center", va="center", fontweight="bold")
    label(ax, "b", x=-0.12)

    # c feature plots
    expr = {"Tcf7": np.isin(ct, [0, 1, 3]) * rng.gamma(3, 0.6, n),
            "Havcr2": np.isin(ct, [4]) * rng.gamma(3, 0.7, n)}
    cmap = LinearSegmentedColormap.from_list("feat", ["#E5E5E5", "#6A51A3", "#3F007D"])
    for k, (g, e) in enumerate(expr.items()):
        ax = fig.add_subplot(top[2 + k])
        e = e * (rng.random(n) > 0.25)
        o = np.argsort(e)
        ax.scatter(xy[o, 0], xy[o, 1], s=0.8, c=e[o], cmap=cmap, vmin=0, vmax=np.percentile(e, 99),
                   lw=0, rasterized=True)
        ax.set_xticks([]); ax.set_yticks([])
        ax.set_title(g, fontsize=6.5, style="italic", pad=2)
        if k == 0:
            label(ax, "c", x=-0.2)

    # d dot plot of markers
    ax = fig.add_subplot(bot[0])
    genes = ["Sell", "Ccr7", "Il7r", "Gzmb", "Cx3cr1", "Tcf7", "Slamf6", "Havcr2", "Tox", "Entpd1",
             "Mki67", "Top2a"]
    owner = [0, 0, 1, 2, 2, 3, 3, 4, 4, 4, 5, 5]
    pct = np.full((len(CELLTYPES), len(genes)), 0.0)
    avg = np.full_like(pct, 0.0)
    for gi, o in enumerate(owner):
        pct[:, gi] = rng.uniform(3, 25, len(CELLTYPES))
        pct[o, gi] = rng.uniform(65, 95)
        if o in (0, 1):
            pct[1 - o if o in (0, 1) else o, gi] = rng.uniform(40, 60)
        avg[:, gi] = (pct[:, gi] - pct[:, gi].mean()) / pct[:, gi].std()
    for r in range(len(CELLTYPES)):
        sc_ = ax.scatter(np.arange(len(genes)), np.full(len(genes), r), s=pct[r] * 0.45,
                         c=avg[r], cmap=RDBU, vmin=-1.5, vmax=1.5, lw=0.2, ec="#555555")
    ax.set_xticks(range(len(genes)), genes, rotation=45, ha="right", style="italic")
    ax.set_yticks(range(len(CELLTYPES)), CELLTYPES)
    ax.set_xlim(-0.6, len(genes) - 0.4); ax.set_ylim(len(CELLTYPES) - 0.5, -0.5)
    cax = ax.inset_axes([1.03, 0.55, 0.025, 0.45])
    cb = fig.colorbar(sc_, cax=cax)
    cb.set_label("Scaled mean\nexpression", fontsize=5); cb.outline.set_linewidth(0.5)
    for v in (25, 50, 75):
        ax.scatter([], [], s=v * 0.45, color="#BBBBBB", ec="#555555", lw=0.2, label=f"{v}%")
    ax.legend(title="% expressing", loc="upper left", bbox_to_anchor=(1.01, 0.42),
              labelspacing=0.7, title_fontsize=5.5, borderaxespad=0, handletextpad=0.2)
    label(ax, "d", x=-0.12)

    # e composition per condition
    ax = fig.add_subplot(bot[1])
    comp = {"Day 8": [0.30, 0.15, 0.30, 0.12, 0.05, 0.08],
            "Day 15": [0.15, 0.12, 0.22, 0.18, 0.23, 0.10],
            "Day 30": [0.08, 0.10, 0.12, 0.20, 0.42, 0.08]}
    left = np.zeros(len(comp))
    yy = np.arange(len(comp))
    for c, nm in enumerate(CELLTYPES):
        v = np.array([comp[k][c] for k in comp]) * 100
        ax.barh(yy, v, left=left, color=CT_COL[nm], height=0.65, label=nm, ec="white", lw=0.3)
        left += v
    ax.set_yticks(yy, list(comp)); ax.invert_yaxis()
    ax.set_xlim(0, 100); ax.set_xlabel("Cells (%)")
    ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=3, handlelength=0.8,
              columnspacing=0.6)
    label(ax, "e", x=-0.3, y=1.02)
    save(fig, "fig13_scrnaseq")


# --------------------------------------------------------------------------
# Fig. 14  scATAC-seq: QC, clusters, chromVAR, coverage tracks
# --------------------------------------------------------------------------
def fig_scatac():
    ct, xy = sim_cells(3500)
    n = len(ct)
    fig = plt.figure(figsize=(DOUBLE, 4.6))
    outer = fig.add_gridspec(2, 1, height_ratios=[1, 1.25], hspace=0.42)
    top = outer[0].subgridspec(1, 3, wspace=0.42)
    bot = outer[1].subgridspec(1, 3, width_ratios=[0.95, 0.85, 1.25], wspace=0.5)

    # a fragment size
    ax = fig.add_subplot(top[0])
    x = np.arange(10, 700)
    dens = (2.2 * stats.norm.pdf(x, 60, 25) + 1.0 * stats.norm.pdf(x, 190, 30)
            + 0.35 * stats.norm.pdf(x, 380, 45) + 0.12 * stats.norm.pdf(x, 560, 60))
    dens += 0.4 * np.exp(-x / 120) / 120 + 0.002 * np.sin(x / 10.4 * 2 * np.pi) * np.exp(-x / 150)
    dens /= dens.sum()
    ax.fill_between(x, 0, dens * 100, color=OI["blue"], lw=0, alpha=0.85)
    for xp, t in [(60, "NFR"), (190, "mono"), (380, "di"), (560, "tri")]:
        ax.text(xp, dens[xp - 10] * 100 + 0.04, t, ha="center", fontsize=5, color="#555555")
    ax.set_xlabel("Fragment length (bp)"); ax.set_ylabel("Fragments (%)")
    ax.set_xlim(0, 700); ax.set_ylim(0, dens.max() * 115)
    label(ax, "a", x=-0.3)

    # b TSS enrichment vs fragments
    ax = fig.add_subplot(top[1])
    good = np.c_[rng.normal(3.9, 0.35, 5000), rng.normal(13, 3.5, 5000)]
    bad = np.c_[rng.normal(2.8, 0.45, 2500), rng.gamma(2, 1.5, 2500)]
    allp = np.vstack([good, bad])
    kde = stats.gaussian_kde(allp.T)(allp.T)
    o = np.argsort(kde)
    ax.scatter(allp[o, 0], allp[o, 1], c=kde[o], cmap="viridis", s=0.6, lw=0, rasterized=True)
    ax.axvline(np.log10(1000), color=OI["red"], lw=0.5, ls="--")
    ax.axhline(4, color=OI["red"], lw=0.5, ls="--")
    ax.set_xlabel("log$_{10}$ unique fragments"); ax.set_ylabel("TSS enrichment")
    ax.set_ylim(0, 25)
    ax.text(0.98, 0.97, f"pass: {((allp[:,0] > 3) & (allp[:,1] > 4)).sum():,} cells",
            transform=ax.transAxes, ha="right", va="top", fontsize=5)
    label(ax, "b", x=-0.3)

    # c UMAP
    ax = fig.add_subplot(top[2])
    umap_ax(ax, xy, [CT_COL[CELLTYPES[c]] for c in ct], "LSI → UMAP (label transfer)")
    for c, nm in enumerate(CELLTYPES):
        cx, cy = xy[ct == c].mean(0)
        ax.text(cx, cy, nm, fontsize=5.5, ha="center", va="center", fontweight="bold")
    label(ax, "c", x=-0.12)

    # d chromVAR heatmap
    ax = fig.add_subplot(bot[0])
    motifs = ["TCF7", "LEF1", "KLF2", "ETS1", "RUNX3", "TBX21", "EOMES", "NR4A1", "BATF",
              "JUNB", "IRF4", "E2F1"]
    prof = np.array([[2.0, 1.4, 0.6, 1.2, -0.6, -1.4],    # TCF7
                     [1.8, 1.5, 0.2, 1.0, -0.8, -1.2],    # LEF1
                     [1.6, 1.0, 0.4, -0.6, -1.4, -0.6],   # KLF2
                     [1.0, 1.2, 0.6, 0.0, -0.8, -0.6],    # ETS1
                     [-0.4, 0.2, 1.2, 0.4, 0.8, 0.2],     # RUNX3
                     [-1.2, -0.2, 1.8, -0.2, 0.8, 0.4],   # TBX21
                     [-1.0, 0.0, 1.4, 0.2, 1.0, 0.2],     # EOMES
                     [-1.4, -1.0, -0.2, 0.6, 1.8, 0.0],   # NR4A1
                     [-1.6, -1.0, 0.4, 0.8, 2.0, 0.6],    # BATF
                     [-1.4, -0.8, 0.6, 0.6, 1.8, 0.4],    # JUNB
                     [-1.2, -0.6, 0.2, 0.4, 1.6, 0.8],    # IRF4
                     [-0.4, -0.4, -0.2, 0.0, 0.2, 2.0]])  # E2F1
    prof += rng.normal(0, 0.15, prof.shape)
    im = ax.imshow(prof, cmap=RDBU, vmin=-2, vmax=2, aspect="auto")
    ax.set_xticks(range(len(CELLTYPES)), CELLTYPES, rotation=45, ha="right")
    ax.set_yticks(range(len(motifs)), motifs)
    ax.tick_params(length=0)
    for s_ in ax.spines.values():
        s_.set_visible(False)
    cb = fig.colorbar(im, ax=ax, fraction=0.06, pad=0.04, aspect=15)
    cb.set_label("chromVAR deviation\n(mean z-score)"); cb.outline.set_linewidth(0.5)
    label(ax, "d", x=-0.45)

    # e UMAP coloured by chromVAR TCF7
    ax = fig.add_subplot(bot[1])
    dev = np.array([prof[0, c] for c in ct]) + rng.normal(0, 0.8, n)
    o = np.argsort(dev)
    sc_ = ax.scatter(xy[o, 0], xy[o, 1], c=dev[o], cmap=RDBU, vmin=-3, vmax=3, s=0.9, lw=0,
                     rasterized=True)
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_xlabel("UMAP1", labelpad=1); ax.set_ylabel("UMAP2", labelpad=1)
    ax.set_title("TCF7 motif deviation", fontsize=6.5, pad=2)
    cb = fig.colorbar(sc_, ax=ax, fraction=0.05, pad=0.03, aspect=15)
    cb.set_label("z-score"); cb.outline.set_linewidth(0.5)
    label(ax, "e", x=-0.15)

    # f coverage tracks per cell type
    tg = bot[2].subgridspec(len(CELLTYPES) + 2, 1, height_ratios=[1] * len(CELLTYPES) + [0.35, 0.6],
                            hspace=0.12)
    xs = np.linspace(0, 40, 800)
    peaks = [5, 11, 18, 24, 31, 36]
    acc = {"Naive":   [1.0, 0.9, 0.2, 1.0, 0.1, 0.2], "Tcm": [0.9, 0.8, 0.3, 0.9, 0.2, 0.2],
           "Tem":     [0.6, 0.4, 0.5, 0.7, 0.3, 0.3], "Tpex": [0.8, 0.9, 0.4, 0.9, 0.5, 0.2],
           "Tex":     [0.2, 0.2, 1.0, 0.5, 1.0, 0.6], "Cycling": [0.4, 0.3, 0.4, 0.5, 0.4, 0.3]}
    for k, nm in enumerate(CELLTYPES):
        ax = fig.add_subplot(tg[k])
        y = rng.gamma(1, 0.03, xs.size)
        for p, a in zip(peaks, acc[nm]):
            y += a * np.exp(-((xs - p) / 0.45) ** 2)
        ax.fill_between(xs, 0, y, color=CT_COL[nm], lw=0)
        ax.set_ylim(0, 1.3); ax.set_xlim(0, 40)
        ax.axis("off")
        ax.text(-0.5, 0.4, nm, ha="right", va="center", fontsize=5.5, color=CT_COL[nm])
        if k == 0:
            ax.text(40, 1.25, "[0–1.3]", ha="right", va="top", fontsize=4.5, color="#555555")
            label(ax, "f", x=-0.28, y=1.3)
    ax = fig.add_subplot(tg[-2]); ax.axis("off"); ax.set_xlim(0, 40); ax.set_ylim(0, 1)
    for p in peaks:
        ax.add_patch(Rectangle((p - 0.7, 0.2), 1.4, 0.6, color="#555555", lw=0))
    ax.text(-0.5, 0.5, "Peaks", ha="right", va="center", fontsize=5.5)
    ax = fig.add_subplot(tg[-1]); ax.axis("off"); ax.set_xlim(0, 40); ax.set_ylim(0, 1)
    ax.plot([4, 30], [0.55, 0.55], color="black", lw=0.5)
    for a, b in [(4, 5.5), (12, 12.8), (19, 19.6), (28.5, 30)]:
        ax.add_patch(Rectangle((a, 0.35), b - a, 0.4, color="black", lw=0))
    ax.text(17, 0.0, "Tcf7", style="italic", ha="center", va="top", fontsize=5.5)
    ax.plot([34, 39], [0.55, 0.55], color="black", lw=0.6)
    ax.text(36.5, 0.65, "5 kb", ha="center", va="bottom", fontsize=5)
    save(fig, "fig14_scatac")


# --------------------------------------------------------------------------
# Fig. 15  Spatial transcriptomics integration
# --------------------------------------------------------------------------
def fig_spatial():
    ncl, per = 6, 350
    ang = np.linspace(0, 2 * np.pi, ncl, endpoint=False)
    centers = np.c_[7 * np.cos(ang), 7 * np.sin(ang)] + rng.normal(0, 0.8, (ncl, 2))
    batch_shift = np.array([[0, 0], [5, -4]])
    pts, cls, bat = [], [], []
    for c in range(ncl):
        for b in range(2):
            pts.append(centers[c] + batch_shift[b] + rng.normal(0, 0.9, (per // 2, 2)))
            cls += [c] * (per // 2); bat += [b] * (per // 2)
    pts, cls, bat = np.vstack(pts), np.array(cls), np.array(bat)
    after = pts - batch_shift[bat] + rng.normal(0, 0.3, pts.shape)
    pal = plt.get_cmap("tab10")
    bcol = np.array([OI["blue"], OI["orange"]])

    fig = plt.figure(figsize=(DOUBLE, 3.9))
    gs = fig.add_gridspec(2, 3, hspace=0.45, wspace=0.3)
    ax = fig.add_subplot(gs[0, 0])
    umap_ax(ax, pts, bcol[bat], "Before Harmony", s=0.8)
    for b, nm in enumerate(["WT", "KO"]):
        ax.scatter([], [], s=6, color=bcol[b], label=nm)
    ylo, yhi = ax.get_ylim(); ax.set_ylim(ylo, yhi + 0.3 * (yhi - ylo))
    ax.legend(loc="upper left", ncol=2, handletextpad=0.1, borderaxespad=0.1)
    label(ax, "a", x=-0.1)
    ax = fig.add_subplot(gs[0, 1])
    umap_ax(ax, after, pal(cls), "After Harmony", s=0.8)
    for c in range(ncl):
        cx, cy = after[cls == c].mean(0)
        ax.text(cx, cy, str(c), fontsize=6, ha="center", va="center", fontweight="bold")
    label(ax, "b", x=-0.1)

    ax = fig.add_subplot(gs[0, 2])
    secs = ["WT_S1", "WT_S2", "KO_S1", "KO_S2"]
    base = np.array([0.22, 0.18, 0.2, 0.15, 0.13, 0.12])
    left = np.zeros(4)
    for c in range(ncl):
        v = np.array([base[c] * (1 + (0.6 if (c == 4 and s.startswith("KO")) else 0)
                                 + rng.normal(0, 0.05)) for s in secs])
        ax.barh(np.arange(4), v, left=left, color=pal(c), height=0.65, ec="white", lw=0.3,
                label=str(c))
        left += v
    ax.set_xlim(0, left.max() * 1.02)
    ax.set_yticks(range(4), secs); ax.invert_yaxis()
    ax.set_xticks([])
    ax.set_xlabel("Spot / cell fraction per cluster")
    ax.legend(title="Cluster", loc="center left", bbox_to_anchor=(1.0, 0.5), handlelength=0.8,
              title_fontsize=5.5)
    label(ax, "c", x=-0.3)

    # tissue sections
    def tissue(seed_shift):
        g = np.mgrid[-1:1:0.028, -1:1:0.028].reshape(2, -1).T
        g = g[(g[:, 0] ** 2 / 1.0 + g[:, 1] ** 2 / 0.6) < 1]
        r = np.sqrt(g[:, 0] ** 2 + (g[:, 1] / 0.77) ** 2)
        band = np.digitize(r + 0.05 * np.sin(6 * np.arctan2(g[:, 1], g[:, 0]) + seed_shift),
                           [0.2, 0.38, 0.55, 0.7, 0.85])
        return g, band
    for k, (title, sh) in enumerate([("WT_S1 — clusters", 0.0), ("KO_S1 — clusters", 1.0)]):
        ax = fig.add_subplot(gs[1, k])
        g, band = tissue(sh)
        if k == 1:
            band = np.where((band == 3) & (g[:, 0] > 0), 4, band)
        ax.scatter(g[:, 0], g[:, 1], c=pal(band), s=1.3, lw=0, rasterized=True)
        ax.set_aspect("equal"); ax.axis("off")
        ax.set_title(title, fontsize=6.5, pad=2)
        if k == 0:
            ax.plot([0.55, 0.95], [-0.85, -0.85], color="black", lw=0.8)
            ax.text(0.75, -0.82, "500 µm", ha="center", va="bottom", fontsize=5)
        label(ax, "de"[k], x=0.0)
    ax = fig.add_subplot(gs[1, 2])
    g, band = tissue(1.0)
    expr = np.where(band == 4, rng.gamma(4, 0.5, len(g)), rng.gamma(1, 0.2, len(g)))
    expr[(band == 3) & (g[:, 0] > 0)] = rng.gamma(4, 0.5, ((band == 3) & (g[:, 0] > 0)).sum())
    sc_ = ax.scatter(g[:, 0], g[:, 1], c=expr, cmap="magma_r", s=1.3, lw=0, rasterized=True,
                     vmax=np.percentile(expr, 99))
    ax.set_aspect("equal"); ax.axis("off")
    ax.set_title("KO_S1 — Gene X", fontsize=6.5, pad=2)
    cb = fig.colorbar(sc_, ax=ax, fraction=0.04, pad=0.02); cb.set_label("log-normalized")
    cb.outline.set_linewidth(0.5)
    label(ax, "f", x=0.0)
    save(fig, "fig15_spatial")


# --------------------------------------------------------------------------
# Fig. 10  Motif analysis results: logos, centrality, footprints
# --------------------------------------------------------------------------
from matplotlib.textpath import TextPath
from matplotlib.patches import PathPatch
from matplotlib.font_manager import FontProperties
import matplotlib.transforms as mtrans

BASE_COL = {"A": "#109648", "C": "#255C99", "G": "#F7B32B", "T": "#D62839"}
_FP = FontProperties(family="Arial", weight="bold")
_GLYPH = {b: TextPath((0, 0), b, size=1, prop=_FP) for b in "ACGT"}
IUPAC = {"A": "A", "C": "C", "G": "G", "T": "T", "S": "CG", "W": "AT", "R": "AG",
         "Y": "CT", "K": "GT", "M": "AC", "N": "ACGT"}


def pwm_from_consensus(cons, strength=0.86):
    rows = []
    for ch in cons:
        allowed = IUPAC[ch]
        p = np.full(4, (1 - strength) / (4 - len(allowed)) if len(allowed) < 4 else 0.25)
        for b in allowed:
            p["ACGT".index(b)] = strength / len(allowed) if len(allowed) < 4 else 0.25
        p = np.clip(p + rng.normal(0, 0.02, 4), 0.005, None)
        rows.append(p / p.sum())
    return np.array(rows)


def draw_logo(ax, pwm):
    ic = 2 + (pwm * np.log2(pwm)).sum(1)          # information content (bits)
    for i, (p, h) in enumerate(zip(pwm, ic)):
        y = 0
        for k in np.argsort(p):
            hh = p[k] * h
            if hh < 0.01:
                continue
            b = "ACGT"[k]
            g = _GLYPH[b]; bb = g.get_extents()
            t = (mtrans.Affine2D().translate(-bb.x0, -bb.y0)
                 .scale(0.92 / bb.width, hh / bb.height).translate(i + 0.04, y))
            ax.add_patch(PathPatch(t.transform_path(g), fc=BASE_COL[b], lw=0))
            y += hh
    ax.set_xlim(0, len(pwm)); ax.set_ylim(0, 2)
    ax.set_xticks([]); ax.set_yticks([0, 2]); ax.tick_params(labelsize=4.5, length=1.5)
    ax.spines["bottom"].set_visible(False)


def fig_motif_results():
    fig = plt.figure(figsize=(DOUBLE, 4.3))
    outer = fig.add_gridspec(2, 1, height_ratios=[1.05, 1], hspace=0.55)
    top = outer[0].subgridspec(1, 2, width_ratios=[1.55, 1], wspace=0.3)
    bot = outer[1].subgridspec(1, 3, width_ratios=[1, 1, 1.1], wspace=0.55)

    # a  HOMER de novo results table with logos
    motifs = [("ATGASTCAT", "BATF::JUN (bZIP)", "1e-41", 24.8, 6.9),
              ("TTTCCA",    "NFATC1 (RHD)",     "1e-22", 31.5, 17.2),
              ("GAAASTGAAA", "IRF4 (IRF)",      "1e-18", 12.1, 4.6),
              ("CTTTGATG",  "TCF7 (HMG)",       "1e-27", 19.7, 7.8)]
    tbl = top[0].subgridspec(len(motifs) + 1, 5, width_ratios=[0.12, 1.5, 1.25, 0.55, 0.85],
                             hspace=0.35, wspace=0.08)
    heads = ["", "De novo motif", "Best match", "P", "% target / bg"]
    for j, h in enumerate(heads):
        axh = fig.add_subplot(tbl[0, j]); axh.axis("off")
        axh.text(0 if j == 1 else 0.5, 0.15, h, ha="left" if j == 1 else "center",
                 va="bottom", fontsize=6, fontweight="bold")
        if j == 0:
            label(axh, "a", x=-3.0, y=0.6)
    for i, (cons, match, p, tg, bgp) in enumerate(motifs, start=1):
        ax = fig.add_subplot(tbl[i, 1]); draw_logo(ax, pwm_from_consensus(cons))
        ax.set_ylabel("bits", fontsize=4.5, labelpad=0.5)
        ax = fig.add_subplot(tbl[i, 2]); ax.axis("off")
        ax.text(0.03, 0.5, match, va="center", fontsize=5.5)
        ax = fig.add_subplot(tbl[i, 3]); ax.axis("off")
        b, e = p.split("e")
        ax.text(0.5, 0.5, f"10$^{{{int(e)}}}$", ha="center", va="center", fontsize=5.5)
        ax = fig.add_subplot(tbl[i, 4])
        ax.barh([1, 0], [tg, bgp], color=[UP if i < 4 else DOWN, "#C8C8C8"], height=0.75)
        ax.set_xlim(0, 35); ax.set_yticks([])
        ax.spines["left"].set_visible(False)
        ax.text(tg + 0.8, 1, f"{tg:.0f}", va="center", fontsize=4.5)
        ax.text(bgp + 0.8, 0, f"{bgp:.0f}", va="center", fontsize=4.5, color="#555555")
        if i < len(motifs):
            ax.set_xticks([])
            ax.spines["bottom"].set_visible(False)
        else:
            ax.set_xticks([0, 15, 30]); ax.tick_params(labelsize=4.5)
            ax.set_xlabel("% regions", fontsize=5, labelpad=1)

    # b  motif centrality around peak summits (annotatePeaks.pl -m -hist)
    ax = fig.add_subplot(top[1])
    x = np.arange(-500, 501, 10)
    for name, h, col in [("TEX-up peaks", 0.055, UP), ("Union background", 0.012, "#8C8C8C")]:
        y = 1000 * (0.008 + h * np.exp(-(x / 70) ** 2) + rng.normal(0, 0.0012, x.size))
        ax.plot(x, y, color=col, label=name)
    ax.set_xlabel("Distance from peak summit (bp)")
    ax.set_ylabel("BATF motifs per kb per peak")
    ax.set_xlim(-500, 500)
    ax.legend(loc="upper right", handlelength=1.2)
    label(ax, "b", x=-0.25)

    # c  % of target regions with motif, up vs down regions
    ax = fig.add_subplot(bot[0])
    names = ["BATF", "AP-1", "NFAT", "IRF4", "RUNX", "ETS1", "LEF1", "TCF7"]
    up = np.array([25, 28, 31, 12, 20, 18, 9, 8])
    dn = np.array([7, 9, 17, 5, 22, 27, 18, 20])
    bgv = np.array([7, 9, 17, 5, 19, 19, 8, 8])
    yy = np.arange(len(names))
    ax.barh(yy - 0.27, up, height=0.27, color=UP, label="Up regions")
    ax.barh(yy, dn, height=0.27, color=DOWN, label="Down regions")
    ax.barh(yy + 0.27, bgv, height=0.27, color="#C8C8C8", label="Background")
    ax.set_yticks(yy, names); ax.invert_yaxis()
    ax.set_xlabel("Regions with motif (%)")
    ax.legend(loc="lower center", bbox_to_anchor=(0.45, 1.0), ncol=3, handlelength=0.8,
              columnspacing=0.6, handletextpad=0.3)
    label(ax, "c", x=-0.35, y=1.12)

    # d  HINT footprint: bias-corrected ATAC cleavage at BATF sites
    ax = fig.add_subplot(bot[1])
    x = np.arange(-100, 101)
    def prof(acc, depth):
        flank = acc * (0.35 + 0.65 * np.exp(-(x / 60) ** 2))
        dip = 1 - depth * np.exp(-(x / 7) ** 2)
        return flank * dip + rng.normal(0, 0.004, x.size)
    ax.plot(x, prof(0.30, 0.55), color="#0072B2", label="TMEM")
    ax.plot(x, prof(0.42, 0.70), color=UP, label="TEX")
    ax.axvspan(-5, 5, color="#E0E0E0", lw=0, zorder=0)
    ax.set_xlabel("Distance from BATF motif (bp)")
    ax.set_ylabel("Corrected Tn5 cleavage (a.u.)")
    ax.set_xlim(-100, 100)
    ax.legend(loc="upper right", handlelength=1.2)
    ax.text(0, 0.02, "motif", ha="center", fontsize=5, color="#555555",
            transform=mtrans.blended_transform_factory(ax.transData, ax.transAxes))
    label(ax, "d", x=-0.3)

    # e  HINT differential TF activity
    ax = fig.add_subplot(bot[2])
    n = 450
    act = rng.normal(0, 0.04, n)
    pv = np.clip(rng.uniform(0, 1, n) ** 0.7, 1e-6, 1)
    hits = {"BATF": (0.32, 1e-9), "JUNB": (0.27, 1e-7), "IRF4": (0.21, 1e-5), "NFATC1": (0.17, 1e-4),
            "TCF7": (-0.29, 1e-8), "LEF1": (-0.24, 1e-6), "KLF2": (-0.18, 1e-4)}
    act = np.r_[act, [v[0] for v in hits.values()]]
    pv = np.r_[pv, [v[1] for v in hits.values()]]
    sig = pv < 1e-3
    col = np.where(sig & (act > 0), UP, np.where(sig & (act < 0), DOWN, NS))
    ax.scatter(act, -np.log10(pv), s=3, c=col, lw=0)
    offs = {"BATF": (-3, 4), "JUNB": (-3, 4), "IRF4": (-3, 4), "NFATC1": (-4, 4),
            "TCF7": (3, 4), "LEF1": (3, 4), "KLF2": (3, 4)}
    for k, (a, b) in hits.items():
        ax.annotate(k, (a, -np.log10(b)), xytext=offs[k], textcoords="offset points",
                    fontsize=5, ha="right" if a > 0 else "left")
    ax.axhline(3, ls="--", lw=0.4, color="#555555")
    ax.axvline(0, lw=0.4, color="black")
    ax.set_xlim(-0.42, 0.42)
    ax.set_xlabel("TF activity score (TEX − TMEM)")
    ax.set_ylabel("−log$_{10}$ P")
    label(ax, "e", x=-0.3)
    save(fig, "fig10_motif_results")


# --------------------------------------------------------------------------
# Fig. 11  Pooled CRISPR screen: guide distributions, QC, ranks, sgRNA level
# --------------------------------------------------------------------------
def fig_crispr():
    n_genes, g_per = 2000, 4
    n_ntc = 400
    ng = n_genes * g_per + n_ntc
    samples = ["Plasmid", "Day 0", "Spleen R1", "Spleen R2", "Tumour R1", "Tumour R2"]
    disp = [0.35, 0.45, 0.9, 0.95, 1.5, 1.6]                  # in vivo bottleneck widens spread
    base = rng.normal(0, 0.35, ng)

    # gene effects (log2FC tumour vs day 0); a few strong hits
    eff = rng.normal(0, 0.45, n_genes) * (rng.random(n_genes) < 0.35)
    pos_hits = {"Dhx37": 3.0, "Pdcd1": 2.6, "Cblb": 2.4, "Ptpn2": 2.2, "Socs1": 2.0, "Rasa2": 1.8}
    neg_hits = {"Myb": -2.9, "Tbx21": -2.5, "Cd28": -2.2, "Il2rb": -2.0, "Bcl2": -1.8}
    names = [f"Gene{i}" for i in range(n_genes)]
    for k, (g, v) in enumerate({**pos_hits, **neg_hits}.items()):
        names[k] = g; eff[k] = v
    guide_eff = np.r_[np.repeat(eff, g_per) * rng.uniform(0.5, 1.2, n_genes * g_per),
                      np.zeros(n_ntc)]
    guide_lfc_r1 = guide_eff + rng.normal(0, 0.35, ng)
    guide_lfc_r2 = guide_eff + rng.normal(0, 0.35, ng)

    fig = plt.figure(figsize=(DOUBLE, 5.6))
    gs = fig.add_gridspec(2, 1, height_ratios=[1, 2.2], hspace=0.32)
    top = gs[0].subgridspec(1, 3, wspace=0.45)
    bot = gs[1].subgridspec(1, 3, width_ratios=[1.4, 1.05, 1.35], wspace=0.6)

    # a  sgRNA count distributions
    ax = fig.add_subplot(top[0, 0])
    cols = ["#000000", "#8C8C8C", OI["sky"], OI["blue"], OI["orange"], OI["red"]]
    grid = np.linspace(0, 16, 300)
    ginis, zeros = [], []
    for nm, dsp, c in zip(samples, disp, cols):
        logc = 9 + base + rng.normal(0, dsp, ng) - (0.8 * dsp ** 2) * (rng.random(ng) < 0.15)
        cnt = np.floor(2 ** logc * (rng.random(ng) > 0.02 * dsp ** 2))
        zeros.append(100 * (cnt == 0).mean())
        x = np.sort(cnt); ginis.append((2 * np.arange(1, ng + 1) - ng - 1) @ x / (ng * x.sum()))
        kde = stats.gaussian_kde(np.log2(cnt + 1))
        ax.plot(grid, kde(grid), color=c, lw=0.7, label=nm)
    ax.set_xlabel("log$_2$(normalized sgRNA count + 1)")
    ax.set_ylabel("Density")
    ax.set_xlim(2, 14)
    ax.legend(loc="upper left", handlelength=1, labelspacing=0.25, fontsize=5)
    label(ax, "a", x=-0.28)

    # b  library QC: Gini index and zero-count guides
    ax = fig.add_subplot(top[0, 1])
    yy = np.arange(len(samples))
    ax.barh(yy, ginis, color=cols, height=0.65)
    ax.axvline(0.2, ls="--", lw=0.4, color="#555555")
    ax.set_yticks(yy, samples); ax.invert_yaxis()
    ax.set_xlabel("Gini index")
    ax.set_xlim(0, max(ginis) * 1.45)
    for yv, gval, z in zip(yy, ginis, zeros):
        ax.text(gval + 0.01, yv, f"{z:.1f}% zero", va="center", fontsize=5, color="#555555")
    label(ax, "b", x=-0.42)

    # c  replicate concordance (gene-level LFC)
    ax = fig.add_subplot(top[0, 2])
    g1 = np.median(guide_lfc_r1[: n_genes * g_per].reshape(n_genes, g_per), 1)
    g2 = np.median(guide_lfc_r2[: n_genes * g_per].reshape(n_genes, g_per), 1)
    r_val = stats.pearsonr(g1, g2)[0]
    ax.scatter(g1, g2, s=1.5, color="#BBBBBB", lw=0, rasterized=True)
    hi = np.arange(len(pos_hits) + len(neg_hits))
    ax.scatter(g1[hi], g2[hi], s=5, c=[UP] * len(pos_hits) + [DOWN] * len(neg_hits), lw=0)
    lim = (-4, 4)
    ax.plot(lim, lim, ls="--", lw=0.4, color="#555555")
    ax.set_xlim(lim); ax.set_ylim(lim)
    ax.set_xlabel("Gene log$_2$FC, replicate 1")
    ax.set_ylabel("Gene log$_2$FC, replicate 2")
    ax.text(0.04, 0.96, f"Pearson r = {r_val:.2f}", transform=ax.transAxes, va="top", fontsize=5.5)
    label(ax, "c", x=-0.3)

    # d  gene rank plot (tall panel; rank 1 = most enriched, on the left)
    ax = fig.add_subplot(bot[0, 0])
    gl = (g1 + g2) / 2
    order = np.argsort(-gl)
    rank = np.empty(n_genes, int); rank[order] = np.arange(1, n_genes + 1)
    colr = np.where(np.isin(np.arange(n_genes), hi[: len(pos_hits)]), UP,
                    np.where(np.isin(np.arange(n_genes), hi[len(pos_hits):]), DOWN, "#BBBBBB"))
    ax.scatter(rank, gl, s=2, c=colr, lw=0, rasterized=True)
    ax.axhline(0, color="black", lw=0.4)
    ax.axhspan(-0.5, 0.5, color="#F0F0F0", lw=0, zorder=0)
    ax.set_xlim(-120, n_genes + 120)
    ax.set_xticks([0, 500, 1000, 1500, 2000])
    ax.set_ylim(-4, 4)
    ax.set_xlabel("Gene rank")
    ax.set_ylabel("Median sgRNA log$_2$FC (tumour / day 0)")
    t_up = ax.text(0.97, 0.99, "enriched\n(KO advantage)", transform=ax.transAxes,
                   ha="right", va="top", fontsize=5, color=UP)
    t_dn = ax.text(0.03, 0.01, "depleted\n(KO disadvantage)", transform=ax.transAxes,
                   ha="left", va="bottom", fontsize=5, color=DOWN)
    avoid = np.abs(gl) > 0.4
    hp, hn = hi[: len(pos_hits)], hi[len(pos_hits):]
    t_pos = repel_labels(ax, rank[hp], gl[hp], list(pos_hits), offset=(14, 0),
                         fixed=(t_up, t_dn), avoid_xy=np.c_[rank[avoid], gl[avoid]])
    repel_labels(ax, rank[hn], gl[hn], list(neg_hits), offset=(14, 0),
                 fixed=(t_up, t_dn, *t_pos), avoid_xy=np.c_[rank[avoid], gl[avoid]])
    label(ax, "d", x=-0.32, y=1.02)

    # e  sgRNA-level plot for top hits vs non-targeting controls
    ax = fig.add_subplot(bot[0, 1])
    show = ["Dhx37", "Pdcd1", "Cblb", "Myb", "Tbx21", "Cd28"]
    lfc_g = (guide_lfc_r1 + guide_lfc_r2) / 2
    ntc = lfc_g[n_genes * g_per:]
    yy = np.arange(len(show) + 1)
    lo, hi_q = np.percentile(ntc, [2.5, 97.5])
    ax.axvspan(lo, hi_q, color="#E8E8E8", lw=0, zorder=0)
    ax.axvline(0, color="black", lw=0.4)
    for k, g in enumerate(show):
        gi = names.index(g)
        v = lfc_g[gi * g_per:(gi + 1) * g_per]
        c = UP if g in pos_hits else DOWN
        ax.scatter(v, np.full(g_per, k) + rng.uniform(-0.12, 0.12, g_per), s=9, color=c,
                   lw=0.3, ec="white", zorder=3)
        ax.plot([np.median(v)] * 2, [k - 0.3, k + 0.3], color="black", lw=0.8)
    vp = ax.violinplot(ntc, positions=[len(show)], vert=False, widths=0.8, showextrema=False)
    for b_ in vp["bodies"]:
        b_.set_facecolor("#BBBBBB"); b_.set_alpha(0.9); b_.set_edgecolor("none")
    ax.set_yticks(yy, show + ["Non-targeting"])
    for tl in ax.get_yticklabels()[:-1]:
        tl.set_style("italic")
    ax.invert_yaxis()
    ax.set_xlabel("sgRNA log$_2$FC (tumour / day 0)")
    ax.set_xlim(-4.5, 4.5)
    ax.set_xticks([-4, -2, 0, 2, 4])
    label(ax, "e", x=-0.55, y=1.02)

    # f  sgRNAs of example genes against the distribution of all sgRNAs:
    #    per gene, the library distribution (top) and that gene's guides (bottom strip)
    dist = ["Pdcd1", "Dhx37", "Ptpn2", "Non-targeting", "Tbx21", "Myb"]
    sub = bot[0, 2].subgridspec(len(dist), 1, hspace=0.35)
    bins = np.linspace(-4.5, 4.5, 91)
    for k, g in enumerate(dist):
        pair = sub[k, 0].subgridspec(2, 1, height_ratios=[1.6, 1], hspace=0.08)
        ax_h = fig.add_subplot(pair[0, 0])
        ax_r = fig.add_subplot(pair[1, 0], sharex=ax_h)
        ax_h.hist(lfc_g, bins=bins, color="#C8C8C8", lw=0, density=True)
        if g == "Non-targeting":
            v, c = ntc[:g_per], "#000000"
        else:
            gi = names.index(g)
            v, c = lfc_g[gi * g_per:(gi + 1) * g_per], UP if g in pos_hits else DOWN
        ax_r.vlines(v, 0, 1, color=c, lw=0.9)
        ax_r.set_ylim(0, 1)
        ax_h.set_xlim(-4.5, 4.5)
        for a_ in (ax_h, ax_r):
            a_.set_yticks([])
            a_.spines["left"].set_visible(False)
        ax_h.tick_params(axis="x", labelbottom=False, length=0)
        ax_h.spines["bottom"].set_color("#BBBBBB")
        ax_r.spines["bottom"].set_visible(False)
        ax_r.set_facecolor("#F4F4F4")
        ax_h.text(-0.03, 0.0, g, transform=ax_h.transAxes, ha="right", va="center", fontsize=5.5,
                  style="normal" if g == "Non-targeting" else "italic")
        if k < len(dist) - 1:
            ax_r.tick_params(axis="x", labelbottom=False, length=0)
        else:
            ax_r.spines["bottom"].set_visible(True)
            ax_r.set_xlabel("sgRNA log$_2$FC (tumour / day 0)")
        if k == 0:
            label(ax_h, "f", x=-0.42, y=1.15)
    save(fig, "fig11_crispr_screen")


if __name__ == "__main__":
    fig_qc()
    fig_peaks()
    fig_profile()
    fig_diff()
    fig_anno()
    fig_motif()
    fig_rna_heatmap()
    fig_tracks()
    fig_tsa()
    fig_motif_results()
    fig_crispr()
    fig_hic()
    fig_scrna()
    fig_scatac()
    fig_spatial()
