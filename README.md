# Bulk RNA-Seq Workflow

Scripts for a bulk RNA-seq analysis, from raw reads to a gene count matrix and downstream analysis. Reads are trimmed with fastp and quantified with Salmon against the transcriptome (cDNA fasta file). Transcript estimates are summarized to genes with tximport. In this repository, we used tomato heat stress tolerance as an example.

The workflow runs in three places:

| Where | What runs there |
| --- | --- |
| Linux HPC (LSF, `bsub`) | Steps 1 and 2: `1_Bulk_RNASeq_fastp.sh` and `2_Salmon.sh` |
| Your computer, RStudio | Every R script: Steps 3.1, 3.2 and 5 to 8 |
| Your computer, VS Code | Every Python script: Step 4 and the three helper scripts |

The HPC storage is mounted on your computer, so the HPC jobs and your local RStudio and VS Code sessions all read and write the same project folder. Nothing is uploaded or downloaded. Only the path to the folder differs (see [Before you run anything](#before-you-run-anything)). The two shell scripts run paired-end reads by default. The changes for single-end reads are listed in [Single-end reads](#single-end-reads).

Author: Md Jahid Hasan Jone

## Workflow

| Step | Script | What it does |
| --- | --- | --- |
| 1 | `1_Bulk_RNASeq_fastp.sh` | Trims and filters raw reads with fastp; writes an HTML and JSON report per sample |
| 1 (helper) | `fastp_summary_from_json_files.py` | Reads the fastp `.json` files in `5_Results/1_fastp/` and writes `fastp_summary.xlsx` there, one row per sample (read counts, Q20/Q30, GC content, filtered reads, duplication rate, adapter-trimmed reads) |
| 1 (helper) | `Fastp_HTML_to_PowerPoint.py` | Opens the fastp `.html` files in `5_Results/1_fastp/` in a headless browser and writes `fastp_charts.pptx` there, one slide per sample with all fastp charts (paired-end reports) |
| 2 | `2_Salmon.sh` | Builds the Salmon index (if missing) and quantifies each sample |
| 2 (helper) | `collect_salmon_mapping_rate.py` | Reads `aux_info/meta_info.json` in each sample folder of `5_Results/2_Salmon/` and writes `salmon_mapping_rates.csv` there with the percent of reads mapped and the library type per sample |
| 3.1 | `3.1_make_tx2gene.R` | Builds the transcript-to-gene table from the GFF |
| 3.2 | `3.2_tximport_combined.R` | Imports Salmon output with tximport and writes gene count and TPM tables |
| 4 | `4_Sample_Relationship_Analysis.py` | Combined figure: PCA, sample correlation heat map and expression violin plot |
| 5 | `5_DESeq2_Analysis.R` | Differential expression with DESeq2 for a list of comparisons, with volcano, MA, heat map, bar and Venn figures |
| 6.1 | `6.1_GeneOntology.R` | GO enrichment (topGO) of the up- and down-regulated genes of each comparison, with a functional annotation figure and a functional enrichment figure |
| 6.2 | `6.2_KEGG_Analysis.R` | KEGG pathway annotation and enrichment of the same DEGs, with the same two figures |
| 7 | `7_TF_Heatmap.R` | Heat map of the log2 fold change of a list of selected genes (for example transcription factors) across selected comparisons, with significance stars |
| 8 | `8_WGCNA.R` | WGCNA co-expression network analysis: modules, module-trait heat maps, hub genes and a hub gene annotation table |

## Requirements

### Linux HPC (Steps 1 and 2)

- An LSF cluster (`bsub`) with conda
- A conda environment with `fastp` and `salmon`:

```bash
conda create -n gatk_rnaseq -c bioconda -c conda-forge fastp salmon
```

Nothing else is needed on the HPC. R and Python are not used there.

### Python (VS Code, on your computer)

Python 3.9 or newer, with the Python extension in VS Code. Step 4 needs:

```bash
pip install pandas numpy matplotlib seaborn scikit-learn
```

The helper scripts also run in VS Code. `collect_salmon_mapping_rate.py` needs nothing extra. `fastp_summary_from_json_files.py` needs `openpyxl`, and `Fastp_HTML_to_PowerPoint.py` needs `python-pptx`, `playwright` and `pillow`, plus a one-time browser download:

```bash
pip install openpyxl python-pptx playwright pillow
playwright install chromium
```

On a Linux computer, if Chromium does not start, also run `playwright install-deps chromium`.

### R (RStudio, on your computer)

Install R and RStudio, then install everything the R scripts use in one go. `BiocManager::install()` handles both Bioconductor and CRAN packages:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

bioc_pkgs <- c("tximport", "rtracklayer", "DESeq2", "topGO", "clusterProfiler",
               "KEGGREST", "impute", "preprocessCore", "GO.db", "AnnotationDbi")

cran_pkgs <- c("dplyr", "tidyr", "readr", "stringr", "purrr", "tibble", "jsonlite",
               "ggplot2", "ggrepel", "ggtext", "ggVennDiagram", "patchwork", "ragg",
               "pheatmap", "RColorBrewer", "extrafont", "WGCNA", "igraph", "magick")

BiocManager::install(c(bioc_pkgs, cran_pkgs), ask = FALSE, update = FALSE)
```

| Packages | Used in |
| --- | --- |
| `rtracklayer`, `dplyr`, `readr` | 3.1 |
| `tximport`, `readr`, `stringr` | 3.2 |
| `DESeq2`, `ggplot2`, `ggrepel`, `pheatmap`, `RColorBrewer`, `ggVennDiagram`, `extrafont`, `patchwork` | 5 |
| `topGO`, `ggplot2`, `dplyr`, `tidyr`, `stringr`, `readr`, `purrr`, `tibble`, `ggtext` | 6.1 |
| `clusterProfiler`, `KEGGREST`, `jsonlite`, `ragg`, plus the CRAN packages of 6.1 | 6.2 (needs an internet connection the first time, to download the KEGG pathway tables) |
| `pheatmap`, `RColorBrewer`, `extrafont` | 7 |
| `WGCNA`, `DESeq2`, `impute`, `preprocessCore`, `ggplot2`, `pheatmap`, `RColorBrewer`, `igraph`, `magick` | 8 |

The scripts for steps 5 to 8 also install any missing package themselves, so the block above is optional. On Windows and Mac, `magick` installs ImageMagick with it. On Linux, install the system libraries first (on Ubuntu, `libmagick++-dev`), and the same applies to the system libraries that `igraph` needs.

## Folder structure

Create one project folder named `Bulk_RNA_seq` (the placeholder paths in the scripts use this spelling, and Linux is case-sensitive). Everything (scripts, references, reads, results, sample information) lives inside it. The same folder layout is used on the HPC and on your computer.

```
Bulk_RNA_seq/
├── 1_Codes/
│   ├── 1_Bulk_RNASeq_fastp.sh
│   ├── 2_Salmon.sh
│   ├── 3.1_make_tx2gene.R
│   ├── 3.2_tximport_combined.R
│   ├── 4_Sample_Relationship_Analysis.py
│   ├── 5_DESeq2_Analysis.R
│   ├── 6.1_GeneOntology.R
│   ├── 6.2_KEGG_Analysis.R
│   ├── 7_TF_Heatmap.R
│   ├── 8_WGCNA.R
│   ├── fastp_summary_from_json_files.py
│   ├── Fastp_HTML_to_PowerPoint.py
│   └── collect_salmon_mapping_rate.py
├── 2_References/
│   ├── ITAG4.0_cDNA.fasta
│   ├── gene_annotation.gff
│   ├── ITAG4.0_goterms.txt
│   ├── go-basic.obo
│   ├── query.ko.txt
│   ├── salmon_tmt_index/      (created by Step 2)
│   └── tx2gene.csv            (created by Step 3.1)
├── 3_Raw_Reads/
├── 4_Trimmed_Reads/
├── 5_Results/
│   ├── 1_fastp/
│   ├── 2_Salmon/
│   ├── 3_tximport/
│   ├── 4_Sample_Relationship_Analysis/
│   ├── 5_DESeq2/
│   ├── 6.1_GO/
│   ├── 6.2_KEGG/
│   ├── 7_TF/
│   │   └── TF_Genes.csv
│   └── 8_WGCNA/
└── Metadata.csv
```

The project folder sits on the HPC storage that is mounted on your computer. The HPC sees it at its cluster path (for example `/share/group_name/user_name/Bulk_RNA_seq`) and your computer sees the same folder through the mount (for example `Z:/Bulk_RNA_seq` on Windows). Steps 1 and 2 only use `1_Codes/1_Bulk_RNASeq_fastp.sh`, `1_Codes/2_Salmon.sh`, `2_References/ITAG4.0_cDNA.fasta`, `3_Raw_Reads/`, `4_Trimmed_Reads/` and `5_Results/`. Put `Metadata.csv` and `5_Results/7_TF/TF_Genes.csv` in the project folder yourself.

To create all the folders at once, run this once from the directory where you want the project, in an HPC terminal (on Windows, if you create the folders from your computer instead, use Git Bash or make them by hand). Add `Metadata.csv` and `5_Results/7_TF/TF_Genes.csv` yourself:

```bash
mkdir -p Bulk_RNA_seq/{1_Codes,2_References,3_Raw_Reads,4_Trimmed_Reads} \
         Bulk_RNA_seq/5_Results/{1_fastp,2_Salmon,3_tximport,4_Sample_Relationship_Analysis,5_DESeq2,6.1_GO,6.2_KEGG,7_TF,8_WGCNA}
```

## Before you run anything

The project folder is the same on both sides, but its path is not: the HPC uses the cluster path and your computer uses the mounted path. The shell scripts run on the HPC and the R and Python scripts run on your computer, so they get different paths even though they sit in the same `1_Codes` folder. Change the following:

**File paths.** Every path starts with the placeholder `/.../.../Bulk_RNA_seq/`. Replace `/.../.../` with the folder that holds your `Bulk_RNA_seq` project folder, as described next.

**HPC paths.** Run this from `1_Codes` on the HPC (replace `/your/hpc/path` with your own path):

```bash
sed -i 's#/\.\.\./\.\.\./Bulk_RNA_seq#/your/hpc/path/Bulk_RNA_seq#g' *.sh
```

**Local paths.** Open the `1_Codes` folder in VS Code, press `Ctrl+Shift+H` (Replace in Files), search for `/.../.../Bulk_RNA_seq`, and replace it with the mounted path of the project folder, for example `Z:/Bulk_RNA_seq` on Windows or `/Volumes/your_mount/Bulk_RNA_seq` on Mac. Use forward slashes on Windows; R and Python both accept them. Set "files to include" to `*.R, *.py`. On Mac, Linux or Git Bash, the `sed` command above also works if you use `*.R *.py` instead of `*.sh`.

**Sample name in `3.1_make_tx2gene.R`.** The last check reads `5_Results/2_Salmon/.../quant.sf`. Replace the `...` with the name of any sample folder.

**Conda location (HPC only).** In `1_Bulk_RNASeq_fastp.sh` and `2_Salmon.sh`, change the path in `conda activate /.../.../usrapps/group/gatk_rnaseq` to the path of your conda environment that has `fastp` and `salmon` (`conda env list` shows the path). The scripts are written for an environment named `gatk_rnaseq`.

## Name your sequencing files correctly

Each paired-end sample has two files, and only the final read indicator differs:

```
SampleName_1.fq.gz    SampleName_2.fq.gz
```

The scripts find samples by looking for `*_1.fq.gz` and then expect a matching `_2.fq.gz`. Two files share a prefix only when they are Read 1 and Read 2 of the same biological sample.

Correct:

```
DP_1F_R1_1.fq.gz    DP_1F_R1_2.fq.gz
```

Incorrect (the `_1` and `_2` files come from different samples, so they will be treated as a pair):

```
DP_1F_1.fq.gz    ← sample A
DP_1F_2.fq.gz    ← sample B
```

Including genotype, treatment and replicate in the prefix helps (for example `Tomato_Heat_R1_1.fq.gz`). Use the same prefix for that sample in `Metadata.csv`.

## Step 1: Trim reads (fastp)

Run this on the HPC.

```bash
cd 1_Codes
bsub < 1_Bulk_RNASeq_fastp.sh
bjobs            # check job status
```

Reads `3_Raw_Reads/`. Writes:

- `4_Trimmed_Reads/TMD_<sample>_1.fq.gz` and `TMD_<sample>_2.fq.gz`
- `5_Results/1_fastp/<sample>.html` and `<sample>.json`

## Step 2: Quantify transcripts (Salmon)

Run this on the HPC.

```bash
bsub < 2_Salmon.sh
bjobs
```

Builds the index in `2_References/salmon_tmt_index` if it does not exist yet, then quantifies every `TMD_*_1.fq.gz` pair in `4_Trimmed_Reads/`. Writes one folder per sample in `5_Results/2_Salmon/`, each containing `quant.sf`.

Set `transcriptome` in the script to your reference transcriptome fasta.

## Step 3.1: Transcript-to-gene table

The GFF (`gene_annotation.gff`) and the transcriptome fasta must come from the same ITAG4.0 release, otherwise transcript IDs in the Salmon output will not match this table.

Run `3.1_make_tx2gene.R` in RStudio. It keeps the `mRNA` features from the GFF and saves `2_References/tx2gene.csv` with two columns, `TXNAME` and `GENEID`. The `mRNA:` and `gene:` prefixes are removed if present.

The script ends with a check. Point it at any sample's `quant.sf`; the printed fraction should be close to 1. A low value means the transcript IDs in the GFF and the fasta do not match.

## Step 3.2: Count matrix (tximport)

Run `3.2_tximport_combined.R` in RStudio. It reads every `quant.sf` under `5_Results/2_Salmon/`, sorts the samples in natural order (`2F` before `10F`), and summarizes to genes with `countsFromAbundance = "lengthScaledTPM"`.

Written to `5_Results/3_tximport/`:

- `txi.rds`: the full tximport object, saved for reference (Step 5 reads `gene_counts.csv`, not this file)
- `gene_counts.csv`
- `TPM.csv`

Each sample folder in `5_Results/2_Salmon/` must have a unique name, because the folder name becomes the column name.

## Step 4: Sample relationships

Run `4_Sample_Relationship_Analysis.py` in VS Code (open the file and use Run Python File, or run this in the VS Code terminal):

```bash
python 4_Sample_Relationship_Analysis.py
```

Reads `5_Results/3_tximport/gene_counts.csv`, removes genes with zero counts in every sample, and applies log2(count + 1). Writes one combined figure, `Figure1_Combined.png` (1000 dpi) and `Figure1_Combined.pdf`, to `5_Results/4_Sample_Relationship_Analysis/`:

- (A) PCA with a 95% confidence ellipse per group
- (B) Pearson correlation between samples
- (C) log2(count + 1) distribution per sample (violin plot)

**Groups.** The group is taken from each sample name: the name is split at `group_sep` and the first `group_fields` pieces are joined. Set both near the top of the script.

| Setting | `2L_T1_72h_R3` becomes |
| --- | --- |
| `group_fields = 1` (default) | `2L` |
| `group_fields = 2` | `2L_T1` |
| `group_fields = -1` (drops the last piece, usually the replicate) | `2L_T1_72h` |

The script prints the group sizes and warns if every sample ends up in its own group or all samples end up in one. Check that printout before using the figure.

**Colors.** Group colors are assigned automatically (10 or fewer groups use a standard 10-color palette, more groups get evenly spaced hues). To use your own colors, fill in `custom_palette` in the script.

**Fonts.** Sizes, weights and styles for all text are set in the `STYLE` dictionary at the top of the script.

The figure is saved at 1000 dpi, and with many samples the PNG is large. Lower the `dpi` in the PNG `plt.savefig` call if it runs out of memory (the PDF is vector and has no `dpi`).

## Step 5: Differential expression (DESeq2)

Run `5_DESeq2_Analysis.R` in RStudio. It runs DESeq2 separately for every comparison in the `comparisons` list, then makes tables and figures for each comparison and for all comparisons together.

**Experimental setup.** The script is written for this experiment:

| Metadata column | Values |
| --- | --- |
| `Genotype` | `1` = CLN1466EA, `2` = NC123S |
| `Tissue` | `F` = Flower, `L` = Leaf |
| `Temperature` | `T0` = control, `T1` = heat |
| `Time` | `24h`, `48h`, `72h` |

Sample names follow `Genotype+Tissue_Temperature_Time_Replicate`, for example `2F_T0_48h_R1` (genotype 2, Flower, T0, 48h, replicate 1).

**Input files.**

- `gene_counts.csv` from Step 3.2 (gene IDs in the first column, samples in the other columns)
- `Metadata.csv`: one row per sample, with the columns `Sample_ID`, `Genotype`, `Tissue`, `Temperature` and `Time`. Every `Sample_ID` must match a column name in `gene_counts.csv` exactly.

There are two ways to give the script each file. Use one per file and comment out the other:

1. **Direct path (default).** Replace `/.../.../` in `counts_file` and `meta_file` with your project path.
2. **`file.choose()`.** A window opens so you can pick the file. Remove the `#` from the `message(...)` and `file.choose()` lines for that file and put a `#` in front of the direct-path line. This works in RStudio.

Set `output_dir` to your `5_Results/5_DESeq2` folder. The script creates the `Tables` and `Figures` subfolders it needs.

**Parameters** (section 3 of the script): adjusted p-value cutoff (`padj_cutoff = 0.05`), log2 fold-change cutoff (`lfc_cutoff = 2`), the minimum total reads for a gene to be kept in a comparison (`min_gene_count = 10`), and how many genes are labeled on volcano plots and shown in heat maps.

**Fold-change shrinkage.** For each comparison the script runs `results()` and then `lfcShrink(type = "normal")`. The `Up` and `Down` calls, the volcano and MA plots and the files read by Steps 6 and 7 all use the shrunken `log2FoldChange`, so the `lfc_cutoff` is applied to the shrunken value. `padj` comes from the Wald test before shrinkage.

**Editing the comparisons.** Each comparison is one line in the `comparisons` list (section 5):

```r
list(name = "Flower vs Leaf", group_col = "Tissue", group1 = "F", group2 = "L", filter = NULL)
```

| Field | Meaning |
| --- | --- |
| `name` | Label used in file names, plots and tables. Must be unique and a valid file name. |
| `group_col` | Metadata column that holds the two groups |
| `group1` | Group compared against `group2`. A positive log2 fold change means higher in `group1`. |
| `group2` | Reference group |
| `filter` | `NULL` to use all samples, or a restriction such as `list(Tissue = "F", Time = "24h")` |

Each comparison is its own DESeq2 model (`design = ~ group_col`) on the samples that pass the filter, so factors that are in neither `group_col` nor `filter` are pooled. To add a comparison, copy a line and change the fields; to remove one, delete the line or put a `#` in front of it. Every line needs a comma at the end except the last one.

After changing the list, check these too:

- `venn_comparisons` (section 10) and `final_comparisons` (near the end) must use the comparison names exactly as written in `name`. The Venn diagram has 6 colors, so it handles up to 6 comparisons.
- The supplementary MA and volcano figures are 5 x 5 grids (25 comparisons, panels a to y). For more than 25 comparisons, increase `ncol` and `nrow` in the two `wrap_plots()` calls.

The current list has 25 comparisons. For a different experiment, also edit the metadata column names used in the PCA and heat map in section 4 of the script (`Genotype`, `Tissue`, `Temperature`, `Time`).

**Outputs** in `5_Results/5_DESeq2/`:

- `Tables/Comparisons/<name>_DESeq2_results.csv`: full results for each comparison, with a `Regulation` column (`Up`, `Down` or `NS`)
- `Tables/`: top genes per comparison, the top 10 genes overall for qPCR validation, the up/down gene counts per comparison, and legend files that match panel letters to comparison names
- `Figures/Volcano_Plots/`, `Figures/MA_Plots/` and `Figures/Heatmaps/`: one figure per comparison (heat maps are skipped when a comparison has fewer than 2 DEGs). `Figures/Heatmaps/` also holds `Heatmap_TopVariableGenes_AllSamples.tiff`, the heat map of the most variable genes across all samples
- `Figures/`: overall sample PCA, bar plot of up/down gene counts, Venn diagram, and the combined MA and volcano figures

All figures are TIFF files saved at 1000 dpi.

**Fonts.** Figures use Times New Roman through `extrafont`. The first time you use `extrafont` on a computer, run `font_import()` once (it takes a few minutes). The script picks the font device by itself (`"win"` on Windows, `"pdf"` on Mac and Linux), so nothing needs changing there. If Times New Roman is not found on your computer, set `FONT <- "serif"`.

## Step 6: GO and KEGG enrichment

Both scripts read the DESeq2 result tables written by Step 5 (`5_Results/5_DESeq2/Tables/Comparisons/<name>_DESeq2_results.csv`). A gene is a DEG when `padj` is below `PADJ_CUTOFF` (0.05) and its absolute log2 fold change is above `LFC_CUTOFF` (2), set in section 2 of each script. These are the same cutoffs Step 5 uses for Up and Down genes (`padj_cutoff = 0.05`, `lfc_cutoff = 2`), so all three scripts call the same DEGs. If you change a cutoff in one script, change it in the others too.

**Input files.**

| File | Used by | What it is |
| --- | --- | --- |
| `2_References/ITAG4.0_goterms.txt` | 6.1 | One gene per line: gene ID, a tab, then comma-separated GO IDs (for example `Solyc01g005000.2` then `GO:0016831,GO:0019752`). Genes without GO terms can have the ID only. |
| `2_References/go-basic.obo` | 6.1 | GO ontology file, used for term names and for the molecular function / cellular component / biological process split |
| `2_References/query.ko.txt` | 6.2 | KAAS output that links each gene ID to a K number |
| `5_Results/5_DESeq2/Tables/Comparisons/` | 6.1, 6.2 | DESeq2 result tables from Step 5 |

As in Step 5, each input file can be given by a direct path (default) or picked in a window with `file.choose()`. Use one per file and comment out the other. For `deseq_dir`, the `file.choose()` option asks you to pick any one DESeq2 result csv and uses the folder it is in.

**Editing the comparisons.** Both scripts take the comparisons from a `comparison_names` vector in section 3 (6.1) or section 2 (6.2):

```r
comparison_names <- c(
  "CLN1466EA vs NC123S",
  "Flower vs Leaf"
)
```

Each name must match `name` in the `comparisons` list of `5_DESeq2_Analysis.R` exactly, because the script reads `<name>_DESeq2_results.csv`. Letters a, b, c and so on are given in the order listed, and they label the figure panels. To add a comparison, add a line; to remove one, delete the line or put a `#` in front of it. Every line needs a comma at the end except the last one. The current list has 8 comparisons. Figure 1 has one row per comparison and Figure 2 one column per comparison, so after changing the number, adjust `ANNOTATION_HEIGHT` and `ENRICHMENT_WIDTH` in 6.1, or `FIG1_HEIGHT` and `FIG2_WIDTH` in 6.2. Use the same list in both scripts so the panel letters agree.

### Step 6.1: Gene Ontology (`6.1_GeneOntology.R`)

Run it in RStudio. For each comparison, the up- and down-regulated genes are tested separately for each GO category (MF, CC, BP) with topGO (`weight01` algorithm, Fisher test, minimum node size 5). The background is all genes in the DESeq2 table that have a GO annotation. p-values are adjusted with Benjamini-Hochberg and called significant at FDR 0.05.

Written to `5_Results/6.1_GO/`:

- `Individual_Results/`: the background, up and down gene lists and the GO result table for each comparison, direction and category
- `GO_All_Results.csv` and `GO_Significant_FDR_0.05.csv`: results for all comparisons combined
- `Functional_Annotation_Top10.csv` and `Functional_Enrichment_Top30.csv`: the terms shown in the two figures
- `Functional_Annotation_Figure_Data.csv`, `Functional_Enrichment_Figure_Data.csv` and `DEG_Summary_All_Comparisons.csv`: plotting data and DEG counts
- `Comparison_Key.csv`: which letter is which comparison
- `Figures/Figure_Functional_Annotation_Top10.tiff` (and `.pdf`): top 10 GO terms per category, with up (red) and down (steel blue) gene counts for each comparison
- `Figures/Figure_Functional_Enrichment_Top30.tiff` (and `.pdf`): top 30 GO terms (10 per category), bar length = number of DEGs, color = -log10(FDR)
- `sessionInfo.txt`

Thresholds, the number of terms shown, fonts and figure sizes are set in section 2. The TIFF files are saved at 1000 dpi.

### Step 6.2: KEGG pathways (`6.2_KEGG_Analysis.R`)

Run it in RStudio. The script attaches K numbers from `query.ko.txt` to the DEGs, maps them to KEGG pathways, and keeps the five main KEGG categories (Metabolism, Genetic Information Processing, Environmental Information Processing, Cellular Processes, Organismal Systems). It downloads the pathway tables from KEGG on the first run and saves them in the output folder as `kegg_category_map.rds` and `ko_pathway_map.rds`, so later runs reuse those two copies. The `enrichKEGG()` step still contacts KEGG on every run; without an internet connection it is skipped with a warning and the rest of the script still runs.

Gene IDs in `query.ko.txt` and in the DESeq2 tables can differ by a trailing transcript number (`Solyc10g079470.3.1` and `Solyc10g079470.3`). The script tries a few ways of trimming the IDs, keeps the one that matches the DESeq2 IDs best, and prints the match rate. Check that number the first time you run it.

Figure 2 uses fold enrichment, calculated directly as (DEGs in the pathway / all DEGs) / (background genes in the pathway / all background genes). A statistical `enrichKEGG()` test is also run and saved, but the figure does not depend on it, so every selected pathway is always shown.

Written to `5_Results/6.2_KEGG/`:

- `01_KO_annotation/`: gene to K number table
- `02_Annotated_DESeq2/`: each DESeq2 table with K numbers and KEGG categories added
- `03_DEG_tables/`: all, up and down DEGs for each comparison
- `04_KEGG_enrichment/`: fold enrichment for all comparisons and the `enrichKEGG()` result for each comparison
- `Pathway_Ranking_AllComparisons.csv`, `Top_Pathways_Figure1.csv` and `Top_Pathways_Figure2.csv`
- `kegg_category_map.rds` and `ko_pathway_map.rds`: the saved KEGG tables described above
- `Comparison_Key.csv`: which letter is which comparison
- `Figure1_Functional_Annotation.tiff`: top 30 pathways with up and down DEG counts
- `Figure2_Functional_Enrichment.tiff`: top 20 pathways by fold enrichment

The number of pathways shown (`TOP_N_ANNOTATION`, `TOP_N_ENRICHMENT`), colors, fonts and figure sizes are set in section 2.

## Step 7: Heat map of selected genes (`7_TF_Heatmap.R`)

Run it in RStudio. It takes a list of genes, such as transcription factors, and draws their log2 fold change across the comparisons you select.

**Input files.**

- A gene list csv with the columns `Gene Name` and `Gene ID` (default `5_Results/7_TF/TF_Genes.csv`). The row order of this file is the row order of the heat map.
- The DESeq2 result tables from Step 5 in `5_Results/5_DESeq2/Tables/Comparisons/`.

Each input can be given by a direct path (default) or picked in a window with `file.choose()` (gene list only), as in Step 5.

**Editing the comparisons.** The columns come from the `selected_comparisons` vector in section 3. Each name must match `name` in the `comparisons` list of `5_DESeq2_Analysis.R` exactly. Letters a, b, c and so on are given in the order listed. To add or remove a comparison, add or delete a line (every line needs a comma at the end except the last one). The current order is the same as `final_comparisons` in Step 5, so the panel letters agree with the combined volcano figure of Step 5. They do not agree with the letters of Step 6, which follow the order of `comparison_names` in those scripts.

**How genes are matched.** A gene is looked up by its `Gene ID` in the `Gene` column of each DESeq2 table. If there is no exact match, the script tries again with the version suffix removed from both sides (`Solyc01g060400.1` and `Solyc01g060400`). The script prints how many genes were found in each comparison. Genes that are not found show as grey cells.

**Heat map.** Colors show log2 fold change on a symmetric blue-white-red scale, so white is no change. Stars show the adjusted p-value (`*` below 0.05, `**` below 0.01, `***` below 0.001). Rows and columns are not clustered. The figure is 6.27 x 8 inches, which fits about 25 genes and 8 comparisons; change the width and height in section 7 for other sizes.

**Fonts.** The figure uses Times New Roman through `extrafont`. The first time you use `extrafont` on a computer, run `font_import()` once in the R console (it takes a few minutes); this script does not run it for you, and you can skip it if you already did it for Step 5. The script picks the font device by itself (`"win"` on Windows, `"pdf"` on Mac and Linux). If the font is not found, set `FONT <- "serif"`.

Written to `5_Results/7_TF/`:

- `TF_Heatmap.tiff` (1000 dpi)
- `TF_Heatmap_Column_Legend.csv`: which letter is which comparison

## Step 8: WGCNA (`8_WGCNA.R`)

Run it in RStudio. With 5000 genes and 1000 dpi figures it can take a while.

**Input files.**

- `gene_counts.csv` from Step 3.2
- `Metadata.csv`: `Sample_ID` in the first column, and the columns `Genotype`, `Tissue` (`F` or `L`), `Temperature` and `Time`. The sample IDs must match the column names of the count file. The trait section (section 10) is written for these columns, so edit it for a different experiment.
- `2_References/gene_annotation.gff`, used only for the hub gene annotation table. If the file is not found, everything else still runs.

**What it does.** Genes with fewer than 1 CPM in at least 20% of the samples (and at least 3 samples) are removed, the counts are transformed with a variance stabilizing transformation, and the 5000 most variable genes are used to build a signed network (bicor correlation, minimum module size 30, modules with eigengenes closer than 0.25 are merged). The soft-thresholding power is the lowest one with a scale-free fit R2 of at least 0.80, or the best one if none reaches it. These settings are in section 3.

**Panels.** Each panel is saved as its own TIFF in `Figures/`:

| Panel | Content |
| --- | --- |
| A | Soft threshold: scale independence and mean connectivity |
| B | Gene dendrogram with modules and per-treatment-group correlation rows |
| C | Eigengene dendrogram and adjacency heat map |
| D | TOM heat map of a random subsample of genes |
| E | Module-trait heat map (Genotype, Tissue, Temperature, Time) |
| F | Module-group heat map, one column per treatment group |
| G | Module membership against gene significance for Genotype within Flower and within Leaf |
| H | Expression heat map and eigengene bar plot of the focal modules |
| I | Co-expression network of the top 5 hub genes of each module in `panel_I_modules`, plus a 2x2 combined image |

The focal modules for G and H are the modules most correlated with Genotype within Flower and within Leaf. Module colors depend on your data, so run the script once, check `Module_Sizes.csv`, and then set `panel_I_modules` in section 23. The 2x2 image is written for 4 modules.

**Other outputs** in `5_Results/8_WGCNA/`:

- Tables: `Filtered_Counts.csv`, `VST_Expression_Matrix.csv`, `WGCNA_Selected_Genes.csv`, `WGCNA_Trait_Matrix.csv`, `Soft_Thresholding_Results.csv`, `Module_Sizes.csv`, `Module_Eigengenes.csv`, `Module_Trait_Correlations.csv`, `Module_Group_Correlations.csv`, `Module_Group_Pvalues.csv`, `Combined_Experimental_Trait_Matrix.csv`, `Metadata_used_for_WGCNA.csv`, `Gene_Module_Membership_KME.csv`, `Gene_Trait_Significance.csv` and `Complete_Gene_WGCNA_Information.csv` (with their p-value files)
- `Hub_Genes/Hub_Genes_<module>.csv`: all genes with |kME| of at least `KME_THRESHOLD` (0.80) in each module
- `Top_5_Hub_Genes_Per_Module.csv`, `Table_1_Top_5_Hub_Genes.csv` and `Hub_Gene_Annotation_Table.csv`: the hub genes shown in Panel I, with the annotation from the GFF (the parsed GFF annotation is also saved as `ITAG4.0_Gene_Annotation.csv`)
- `WGCNA_Network.rds`, `Module_Eigengenes.rds`, `WGCNA_Expression_Matrix.rds` and `Complete_WGCNA_Analysis.RData`
- `01_Sample_Clustering.pdf`, `02_Sample_Correlation_Heatmap.pdf` and `07_Module_Size.pdf`

## Single-end reads

Both shell scripts (run on the HPC) default to paired-end. Their loops look for `*_1.fq.gz` files and skip any sample that has no matching `_2.fq.gz`, so single-end files are skipped unless you make the changes below. The comment above the `fastp` command ("remove both lines for single end") only covers the `-I` and `-O` lines, not the loop. For single-end data:

**`1_Bulk_RNASeq_fastp.sh`**

1. In the `for` loop and in `basename`, change `_1.fq.gz` to your single-end file suffix (for example `*.fq.gz` and `.fq.gz`).
2. Set `filename1=${sample}.fq.gz` (same suffix).
3. Delete the `filename2` line and remove `&& [ -f $INPUTDIR/$filename2 ]` from the `if` test.
4. Delete the `-I` and `-O` lines from the `fastp` command.

**`2_Salmon.sh`**

1. In the `for` loop and in `basename`, change `_1.fq.gz` to your single-end suffix (for example `TMD_*.fq.gz` and `.fq.gz`).
2. Set `filename1=TMD_${sample}.fq.gz`.
3. Delete the `filename2` line and remove `&& [ -f $INPUTDIR/$filename2 ]` from the `if` test.
4. In `salmon quant`, replace the `-1` and `-2` lines with `-r $INPUTDIR/$filename1 \`.

The R scripts need no change for single-end data.

## License

MIT. See `LICENSE`.
