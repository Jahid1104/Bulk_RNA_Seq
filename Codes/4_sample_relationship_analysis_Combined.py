"""
Sample Relationship Analysis - Combined (PE + SE) gene counts
=================================================================
Recreates Figure 2 style panels from Li et al. (2025, Czech J. Genet. Plant
Breed.) using the combined PE+SE gene count table produced by
`tximport_combined.R` (New_Name/Results/tximport/Combined_gene_counts.csv):

    (A) PCA of all expression samples, colored by 1F/1L/2F/2L group,
        with a 95% confidence ellipse drawn around each group
    (B) Sample-to-sample Pearson correlation heat map
    (C) Violin plot of the gene expression distribution (log2(count+1)) per sample

Sample naming / groups
-----------------------
Column names are the Salmon sample names produced by tximport, e.g.
"1F_T0_24h_R1", "2L_T1_72h_R3", etc. Per Metadata.csv:
    - leading digit  -> genotype: 1 = CLN1466EA, 2 = NC123S
    - trailing letter -> tissue:   F = Flower,    L = Leaf
So the leading "<digit><letter>" combination (1F, 1L, 2F, 2L) defines four
genotype x tissue groups, and every plot below is colored by that group
(instead of just Flower/Leaf).

All figures are exported at 1000 dpi (same as the notebook) as PNG + PDF
into `output_dir` below - no Colab upload/download needed, everything
reads/writes local files.

STYLE / FONT CONTROL
---------------------
Every font used across every panel (titles, axis labels, tick labels,
legends, PCA sample labels, panel tags a/b/c) is driven from the single
STYLE dictionary in section 0 below. Change size / weight ("bold" or
"normal") / style ("italic" or "normal") there and every panel picks it up
automatically - no need to hunt through the plotting code.
"""
# -----------------------------
# Check and install required packages
# -----------------------------
import importlib.util
import subprocess
import sys

_required = {
    "pandas": "pandas",
    "numpy": "numpy",
    "matplotlib": "matplotlib",
    "seaborn": "seaborn",
    "sklearn": "scikit-learn",
}

_missing = [
    pip_name
    for mod_name, pip_name in _required.items()
    if importlib.util.find_spec(mod_name) is None
]

if _missing:
    print(f"Installing missing packages: {_missing}")
    subprocess.check_call([
        sys.executable, "-m", "pip", "install", *_missing
    ])

import re
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.transforms as transforms
from matplotlib.patches import Ellipse
import seaborn as sns
from sklearn.decomposition import PCA

sns.set_style("white")
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['svg.fonttype'] = 'none'


# =====================================================================
# 0. STYLE / FONT SETTINGS  <- EDIT ANYTHING HERE
# =====================================================================
# size   -> font size in points
# weight -> "normal" or "bold"
# style  -> "normal" or "italic"
#
# Every ax.set_title / set_xlabel / tick_params / legend call below reads
# from this dict, so changing a value here changes it everywhere that
# element is used (standalone panels AND the combined Figure2_Combined_All).
STYLE = {
    "panel_title":   {"size": 13, "weight": "bold",   "style": "normal"},  # "(A) Principal Component Analysis" etc.
    "axis_label":    {"size": 12, "weight": "normal",  "style": "normal"},  # PC1/PC2, "log2(count+1)"
    "tick_label":    {"size": 9,  "weight": "normal",  "style": "normal"},  # numeric / sample-name tick labels
    "sample_label":  {"size": 6,  "weight": "normal",  "style": "normal"},  # sample-name text next to PCA points
    "legend_title":  {"size": 10, "weight": "normal",  "style": "normal"},
    "legend_text":   {"size": 9,  "weight": "normal",  "style": "normal"},
    "heatmap_annot": {"size": 6,  "weight": "normal",  "style": "normal"},  # numbers inside the correlation heatmap
    "panel_tag":     {"size": 14, "weight": "bold",    "style": "normal"},  # (A) (B) (C) tags on the combined figure
}


def _font_kwargs(key):
    """Return matplotlib fontdict-style kwargs (size/weight/style) for STYLE[key]."""
    s = STYLE[key]
    return {"fontsize": s["size"], "fontweight": s["weight"], "fontstyle": s["style"]}


def apply_tick_style(ax, axis="both"):
    """Apply the tick_label style to an axes' tick labels (works after rotation too)."""
    s = STYLE["tick_label"]
    if axis in ("both", "x"):
        for lbl in ax.get_xticklabels():
            lbl.set_fontsize(s["size"])
            lbl.set_fontweight(s["weight"])
            lbl.set_fontstyle(s["style"])
    if axis in ("both", "y"):
        for lbl in ax.get_yticklabels():
            lbl.set_fontsize(s["size"])
            lbl.set_fontweight(s["weight"])
            lbl.set_fontstyle(s["style"])


def style_legend(legend):
    """Apply legend_title / legend_text style to an existing legend object."""
    if legend is None:
        return
    t = legend.get_title()
    t.set_fontsize(STYLE["legend_title"]["size"])
    t.set_fontweight(STYLE["legend_title"]["weight"])
    t.set_fontstyle(STYLE["legend_title"]["style"])
    for txt in legend.get_texts():
        txt.set_fontsize(STYLE["legend_text"]["size"])
        txt.set_fontweight(STYLE["legend_text"]["weight"])
        txt.set_fontstyle(STYLE["legend_text"]["style"])


def confidence_ellipse(x, y, ax, n_std=2.4477, **kwargs):
    """
    Draw a 95% confidence ellipse for the 2D points (x, y) onto ax.

    n_std = 2.4477 corresponds to sqrt(chi2.ppf(0.95, df=2)), i.e. the
    scaling that makes this a 95% confidence ellipse for a bivariate
    normal (this is the standard "95% CI ellipse" used for PCA score
    plots, matching e.g. ggplot2's stat_ellipse(level = 0.95, type="norm")).
    Needs at least 3 points per group; groups with fewer points are skipped.
    """
    if len(x) < 3:
        return None

    cov = np.cov(x, y)
    pearson = cov[0, 1] / np.sqrt(cov[0, 0] * cov[1, 1])
    ell_radius_x = np.sqrt(1 + pearson)
    ell_radius_y = np.sqrt(1 - pearson)

    ellipse = Ellipse((0, 0), width=ell_radius_x * 2, height=ell_radius_y * 2, **kwargs)

    scale_x = np.sqrt(cov[0, 0]) * n_std
    mean_x = np.mean(x)
    scale_y = np.sqrt(cov[1, 1]) * n_std
    mean_y = np.mean(y)

    transf = (
        transforms.Affine2D()
        .rotate_deg(45)
        .scale(scale_x, scale_y)
        .translate(mean_x, mean_y)
    )
    ellipse.set_transform(transf + ax.transData)
    return ax.add_patch(ellipse)


def plot_pca_panel(ax, pca_df, palette, group_order, explained,
                    point_size=55, label_points=True, draw_ellipses=True):
    """Shared PCA scatter + 95% ellipse plotting, used by both the standalone
    Panel A figure and the combined Figure2_Combined_All panel A."""
    all_groups = group_order + [g for g in pca_df['Group'].unique() if g not in group_order]

    for grp in all_groups:
        sub = pca_df[pca_df['Group'] == grp]
        if sub.empty:
            continue
        color = palette.get(grp, '#888888')

        ax.scatter(sub['PC1'], sub['PC2'], s=point_size, alpha=0.85,
                   label=grp, color=color, edgecolor='black', linewidth=0.4)

        if draw_ellipses:
            confidence_ellipse(
                sub['PC1'].values, sub['PC2'].values, ax,
                facecolor=color, edgecolor=color, alpha=0.15, linewidth=1.2, linestyle='--'
            )

    if label_points:
        for sname, row in pca_df.iterrows():
            ax.annotate(sname, (row['PC1'], row['PC2']), xytext=(3, 3), textcoords='offset points',
                        fontsize=STYLE["sample_label"]["size"],
                        fontweight=STYLE["sample_label"]["weight"],
                        fontstyle=STYLE["sample_label"]["style"])

    ax.set_xlabel(f"PC1 ({explained[0]:.1f}%)", **_font_kwargs("axis_label"))
    ax.set_ylabel(f"PC2 ({explained[1]:.1f}%)", **_font_kwargs("axis_label"))
    apply_tick_style(ax)


# -----------------------------
# 1. User paths (edit these)
# -----------------------------
# Gene count CSV produced by tximport_PE.R (write.csv(txi_pe$counts, ...))
csv_filename = "New_Name/Results/3_tximport/gene_counts.csv"

# Where the figures for this run get saved
output_dir = "New_Name/Results/4_Sample_Relationships"
os.makedirs(output_dir, exist_ok=True)


# -----------------------------
# 2. Load data and assign sample groups
# First column is treated as the gene ID and set as the index.
# -----------------------------
df = pd.read_csv(csv_filename, index_col=0)
print("Shape (genes x samples):", df.shape)
print(df.head())

# Keep only numeric expression columns (safety check)
df = df.apply(pd.to_numeric, errors='coerce')
df = df.dropna(how='all')
print("After cleaning, shape:", df.shape)


# Assign group (1F / 1L / 2F / 2L) based on the leading "<digit><letter>"
# of each sample name, e.g. "2L_T1_72h_R3" -> "2L"
def get_group(sample_name):
    s = str(sample_name).strip().upper()
    m = re.match(r'^(\d+)([FL])', s)
    if m:
        return f"{m.group(1)}{m.group(2)}"
    return "Unknown"


sample_groups = pd.Series({s: get_group(s) for s in df.columns}, name='Group')
print(sample_groups.value_counts())

if (sample_groups == 'Unknown').any():
    print("\nWARNING: some samples could not be classified into a 1F/1L/2F/2L group:")
    print(sample_groups[sample_groups == 'Unknown'])


# -----------------------------
# 3. Filter and transform expression data
# - Remove genes with zero expression across all samples.
# - Apply a log2(count + 1) transform (same convention as log2(FPKM+1)
#   used in the paper's Figure 2C).
# -----------------------------
# Remove genes that are all-zero across samples
expr_raw = df.loc[(df.sum(axis=1) > 0)]
print("Genes retained after removing all-zero rows:", expr_raw.shape[0])

# log2(count + 1) transform
expr_log = np.log2(expr_raw + 1)
print(expr_log.head())


# -----------------------------
# Shared color palette for the four genotype x tissue groups
# Flower = warm colors, Leaf = cool/green colors; genotype 1 vs 2 = shade
# -----------------------------
palette = {
    '1F': '#E69F00',   # genotype 1, Flower - orange
    '2F': '#D55E00',   # genotype 2, Flower - vermillion
    '1L': '#009E73',   # genotype 1, Leaf   - teal green
    '2L': '#004D40',   # genotype 2, Leaf   - dark green
    'Unknown': '#888888',
}
group_order = ['1F', '1L', '2F', '2L']


# -----------------------------
# 4. Panel A - PCA of all expression samples (colored by 1F/1L/2F/2L,
#    with a 95% confidence ellipse per group)
# -----------------------------
# Samples as rows, genes as columns for PCA
X = expr_log.T.values
sample_names = expr_log.columns.tolist()
groups = sample_groups.loc[sample_names].values

# Center (mean-subtract) genes - matches typical prcomp()-style PCA on expression data
X_centered = X - X.mean(axis=0)

pca = PCA(n_components=min(10, X_centered.shape[0] - 1))
pca_scores = pca.fit_transform(X_centered)
explained = pca.explained_variance_ratio_ * 100

pca_df = pd.DataFrame(pca_scores[:, :2], columns=['PC1', 'PC2'], index=sample_names)
pca_df['Group'] = groups

fig, ax = plt.subplots(figsize=(12, 10))
plot_pca_panel(ax, pca_df, palette, group_order, explained, point_size=55)
ax.set_title("Principal Component Analysis", **_font_kwargs("panel_title"))
legend = ax.legend(title="Group", frameon=False)
style_legend(legend)
sns.despine()
plt.tight_layout()
plt.savefig(os.path.join(output_dir, "Figure2A_PCA_Combined.png"), dpi=1000, bbox_inches='tight')
plt.savefig(os.path.join(output_dir, "Figure2A_PCA_Combined.pdf"), bbox_inches='tight')
plt.show()


# -----------------------------
# 5. Panel B - Sample correlation heat map
# -----------------------------
corr_matrix = expr_log.corr(method='pearson')


# Order samples: group blocks in the order 1F, 1L, 2F, 2L, and within each
# group, natural-sort the full sample name (so "..._R2" sorts after "..._R1"
# and "10" sorts after "2", not before it - same idea as sort_samples() in
# the R tximport scripts).
def natural_key(sample_name):
    s = str(sample_name).strip()
    grp = sample_groups[sample_name]
    grp_rank = group_order.index(grp) if grp in group_order else len(group_order)
    # pad every run of digits in the full name for natural sorting
    padded = re.sub(r'\d+', lambda m: m.group().zfill(10), s)
    return (grp_rank, padded)


order = sorted(sample_groups.index, key=natural_key)
corr_matrix = corr_matrix.loc[order, order]

n_samples = corr_matrix.shape[0]
annotate = n_samples <= 20   # auto-hide numeric annotations if too many samples to read

fig, ax = plt.subplots(figsize=(max(8, n_samples * 0.25), max(10, n_samples * 0.25)))
sns.heatmap(
    corr_matrix,
    cmap="RdBu_r",
    vmin=corr_matrix.values.min(), vmax=1.0,
    square=True,
    annot=annotate, fmt=".2f" if annotate else None,
    annot_kws={
        "size": STYLE["heatmap_annot"]["size"],
        "weight": STYLE["heatmap_annot"]["weight"],
        "style": STYLE["heatmap_annot"]["style"],
    } if annotate else None,
    cbar_kws={
        'label': 'Value',
        'shrink': 0.4,     # smaller = shorter colorbar (fraction of ax height); try 0.3-0.5
        'aspect': 25,       # bigger = thinner/narrower colorbar; try 20-40
        'pad': 0.02         # space between heatmap and colorbar
    },
    linewidths=0.2, linecolor='white',
    ax=ax
)
ax.set_title("Sample correlation", **_font_kwargs("panel_title"))
ax.set_xlabel("")
ax.set_ylabel("")
plt.xticks(rotation=90)
plt.yticks(rotation=0)
apply_tick_style(ax)
plt.tight_layout()
plt.savefig(os.path.join(output_dir, "Figure2B_CorrelationHeatmap_Combined.png"), dpi=1000, bbox_inches='tight')
plt.savefig(os.path.join(output_dir, "Figure2B_CorrelationHeatmap_Combined.pdf"), bbox_inches='tight')
plt.show()


# -----------------------------
# 6. Panel C - Violin plot of expression distribution per sample
# -----------------------------
plot_df = expr_log.reset_index().melt(id_vars=expr_log.index.name or 'index',
                                       var_name='Sample', value_name='log2(count+1)')
plot_df['Group'] = plot_df['Sample'].map(sample_groups)

# keep the same sample order as the correlation heatmap
sample_order = order

fig, ax = plt.subplots(figsize=(max(10, n_samples * 0.3), 6))
sns.violinplot(
    data=plot_df, x='Sample', y='log2(count+1)', order=sample_order,
    hue='Group', hue_order=group_order, dodge=False, palette=palette, cut=0, linewidth=0.5,
    ax=ax
)
ax.set_title("Expression distribution", **_font_kwargs("panel_title"))
ax.set_xlabel("")
ax.set_ylabel(r"$\log_2(\mathrm{count}+1)$", **_font_kwargs("axis_label"))
plt.xticks(rotation=90)
apply_tick_style(ax)
legend = ax.legend(frameon=False, loc='upper right')
style_legend(legend)
sns.despine()
plt.tight_layout()
plt.savefig(os.path.join(output_dir, "Figure2C_ViolinPlot_Combined.png"), dpi=1000, bbox_inches='tight')
plt.savefig(os.path.join(output_dir, "Figure2C_ViolinPlot_Combined.pdf"), bbox_inches='tight')
plt.show()


# -----------------------------
# 7. Combined Figure 2 (A + B + C panels together)
# -----------------------------
fig = plt.figure(figsize=(22, 16))
gs = fig.add_gridspec(2, 30, height_ratios=[1.1, 1.0],
                       hspace=0.2, wspace=0.0)

# Panel A - row 1, first third
ax1 = fig.add_subplot(gs[0, 0:10])
plot_pca_panel(ax1, pca_df, palette, group_order, explained, point_size=45)
ax1.set_title("(A) Principal Component Analysis", loc='left',
              fontsize=STYLE["panel_tag"]["size"], fontweight=STYLE["panel_tag"]["weight"],
              fontstyle=STYLE["panel_tag"]["style"])
legend1 = ax1.legend(frameon=False)
style_legend(legend1)
sns.despine(ax=ax1)


# Panel B - row 1, remaining columns
ax2 = fig.add_subplot(gs[0, 12:30])
sns.heatmap(
    corr_matrix,
    cmap="RdBu_r",
    vmin=corr_matrix.values.min(), vmax=1.0,
    square=False,
    annot=annotate, fmt=".2f" if annotate else None,
    annot_kws={
        "size": STYLE["heatmap_annot"]["size"],
        "weight": STYLE["heatmap_annot"]["weight"],
        "style": STYLE["heatmap_annot"]["style"],
    } if annotate else None,
    cbar_kws={
        'label': 'Value',
        'shrink': 0.4,
        'aspect': 25,
        'pad': 0.01
    },
    linewidths=0.2, linecolor='white',
    ax=ax2
)
# no set_aspect / set_anchor here - let the heatmap fill the full gridspec cell,
# so its height matches Panel A's height (same row) automatically
ax2.tick_params(axis='x', labelrotation=90)
ax2.tick_params(axis='y', labelrotation=0)
apply_tick_style(ax2)
ax2.set_title("(B) Sample correlation", loc='left',
              fontsize=STYLE["panel_tag"]["size"], fontweight=STYLE["panel_tag"]["weight"],
              fontstyle=STYLE["panel_tag"]["style"])


# Panel C - row 2, full width
ax3 = fig.add_subplot(gs[1, :])
sns.violinplot(data=plot_df, x='Sample', y='log2(count+1)', order=sample_order,
               hue='Group', hue_order=group_order, dodge=False, palette=palette, cut=0,
               linewidth=0.4, ax=ax3)
ax3.set_title("(C) Expression distribution", loc='left',
              fontsize=STYLE["panel_tag"]["size"], fontweight=STYLE["panel_tag"]["weight"],
              fontstyle=STYLE["panel_tag"]["style"])
ax3.set_xlabel("")
ax3.set_ylabel(r"$\log_2(\mathrm{count}+1)$", **_font_kwargs("axis_label"))
ax3.tick_params(axis='x', labelrotation=90)
apply_tick_style(ax3)
legend3 = ax3.legend(frameon=False, loc='upper right')
style_legend(legend3)
sns.despine(ax=ax3)

plt.savefig(os.path.join(output_dir, "Figure2_Combined_All.png"), dpi=1000, bbox_inches='tight')
plt.savefig(os.path.join(output_dir, "Figure2_Combined_All.pdf"), bbox_inches='tight')
plt.show()

print(f"\nDone. All Combined figures saved to: {output_dir}")