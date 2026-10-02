##Bulk RNA Seq Data Analysis Workflow by Md Jahid Hasan Jone##

# Sample relationship analysis from the tximport gene count table.
# One combined figure (Figure 1):
#   (A) PCA of all samples, colored by group, with a 95% confidence ellipse per group
#   (B) Sample-to-sample Pearson correlation heat map
#   (C) Violin plot of log2(count + 1) per sample
#
# STYLE / FONT CONTROL
# Every font (titles, axis labels, tick labels, legends, PCA sample labels,
# panel tags) is set in the STYLE dictionary below. Change size / weight
# ("bold" or "normal") / style ("italic" or "normal") there.

import re
import os
import numpy as np
import pandas as pd
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
# STYLE / FONT SETTINGS  <- EDIT ANYTHING HERE
# =====================================================================
STYLE = {
    "panel_title":   {"size": 13, "weight": "bold",   "style": "normal"},
    "axis_label":    {"size": 12, "weight": "normal",  "style": "normal"},  # PC1/PC2, "log2(count+1)"
    "tick_label":    {"size": 9,  "weight": "normal",  "style": "normal"},  # numeric / sample-name tick labels
    "sample_label":  {"size": 6,  "weight": "normal",  "style": "normal"},  # sample-name text next to PCA points
    "legend_title":  {"size": 10, "weight": "normal",  "style": "normal"},
    "legend_text":   {"size": 9,  "weight": "normal",  "style": "normal"},
    "heatmap_annot": {"size": 6,  "weight": "normal",  "style": "normal"},  # numbers inside the correlation heatmap
    "panel_tag":     {"size": 14, "weight": "bold",    "style": "normal"},  # (A) (B) (C) panel titles
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
    n_std = 2.4477 = sqrt(chi2.ppf(0.95, df=2)).
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
    """PCA scatter + 95% ellipse per group."""
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
# Paths
# -----------------------------
csv_filename = "/.../.../Bulk_RNA_seq/5_Results/3_tximport/gene_counts.csv"
output_dir = "/.../.../Bulk_RNA_seq/5_Results/4_Sample_Relationship_Analysis"

os.makedirs(output_dir, exist_ok=True)


# -----------------------------
# Sample groups (edit these two lines to match your sample names)
# -----------------------------
# The group is taken from the sample name: split the name at group_sep and
# join the first group_fields pieces.
#   group_fields =  1   "2L_T1_72h_R3"   -> "2L"          (default)
#   group_fields =  2   "2L_T1_72h_R3"   -> "2L_T1"
#   group_fields = -1   "2L_T1_72h_R3"   -> "2L_T1_72h"   (drops only the last piece, usually the replicate)
group_sep = "_"
group_fields = 1

# Optional: your own colors, e.g. {'Control': '#1f77b4', 'Heat': '#d62728'}.
# Leave empty to get colors automatically. Groups missing from this dict also get automatic colors.
custom_palette = {}


# -----------------------------
# Load data (first column = gene ID)
# -----------------------------
df = pd.read_csv(csv_filename, index_col=0)
print("Shape (genes x samples):", df.shape)

df = df.apply(pd.to_numeric, errors='coerce')
df = df.dropna(how='all')
print("After cleaning, shape:", df.shape)


def get_group(sample_name):
    pieces = str(sample_name).strip().split(group_sep)
    group = group_sep.join(pieces[:group_fields])
    return group if group else "Unknown"


sample_groups = pd.Series({s: get_group(s) for s in df.columns}, name='Group')
print(sample_groups.value_counts())

if sample_groups.nunique() == len(sample_groups):
    print("\nWARNING: every sample is in its own group. Check group_sep and group_fields.")
if sample_groups.nunique() == 1:
    print("\nWARNING: all samples are in one group. Check group_sep and group_fields.")


# -----------------------------
# Filter and transform
# -----------------------------
# Remove genes with zero counts in all samples, then log2(count + 1)
expr_raw = df.loc[(df.sum(axis=1) > 0)]
print("Genes retained after removing all-zero rows:", expr_raw.shape[0])

expr_log = np.log2(expr_raw + 1)


# -----------------------------
# Group order and colors
# -----------------------------
def natural_pad(text):
    # pad every run of digits so "10" sorts after "2"
    return re.sub(r'\d+', lambda m: m.group().zfill(10), str(text).strip())


group_order = sorted(sample_groups.unique(), key=natural_pad)

# One color per group: tab10 for up to 10 groups, evenly spaced hues for more
auto_colors = sns.color_palette("tab10" if len(group_order) <= 10 else "husl", len(group_order))
palette = {g: custom_palette.get(g, c) for g, c in zip(group_order, auto_colors)}


# -----------------------------
# PCA (Panel A data)
# -----------------------------
# Samples as rows, genes as columns
X = expr_log.T.values
sample_names = expr_log.columns.tolist()
groups = sample_groups.loc[sample_names].values

# Center (mean-subtract) genes
X_centered = X - X.mean(axis=0)

pca = PCA(n_components=min(10, X_centered.shape[0] - 1))
pca_scores = pca.fit_transform(X_centered)
explained = pca.explained_variance_ratio_ * 100

pca_df = pd.DataFrame(pca_scores[:, :2], columns=['PC1', 'PC2'], index=sample_names)
pca_df['Group'] = groups


# -----------------------------
# Correlation matrix (Panel B data)
# -----------------------------
corr_matrix = expr_log.corr(method='pearson')


# Order samples: group blocks in group_order, and within each group a natural
# sort of the full sample name
def natural_key(sample_name):
    grp_rank = group_order.index(sample_groups[sample_name])
    return (grp_rank, natural_pad(sample_name))


order = sorted(sample_groups.index, key=natural_key)
corr_matrix = corr_matrix.loc[order, order]

n_samples = corr_matrix.shape[0]
annotate = n_samples <= 20   # hide numbers in the heat map if there are too many samples


# -----------------------------
# Violin plot data (Panel C data)
# -----------------------------
plot_df = expr_log.reset_index().melt(id_vars=expr_log.index.name or 'index',
                                       var_name='Sample', value_name='log2(count+1)')
plot_df['Group'] = plot_df['Sample'].map(sample_groups)

sample_order = order


# -----------------------------
# Figure 1: A + B + C combined
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

plt.savefig(os.path.join(output_dir, "Figure1_Combined.png"), dpi=1000, bbox_inches='tight')
plt.savefig(os.path.join(output_dir, "Figure1_Combined.pdf"), bbox_inches='tight')
plt.close()

print(f"\nDone. Figure 1 saved to: {output_dir}")