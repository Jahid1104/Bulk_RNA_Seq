# Bulk RNA-Seq Workflow

Scripts for a bulk RNA-seq analysis of tomato, from raw reads to a gene count matrix and downstream analysis. Reads are trimmed with fastp and quantified with Salmon against the ITAG4.0 transcriptome. Transcript estimates are summarized to genes with tximport.

The trimming and quantification scripts are LSF job scripts (`bsub`) and run paired-end reads by default. The changes for single-end reads are listed in [Single-end reads](#single-end-reads).

Author: Md Jahid Hasan Jone

## Workflow

| Step | Script | What it does |
| --- | --- | --- |
| 1 | `1_Bulk_RNASeq_fastp.sh` | Trims and filters raw reads with fastp; writes an HTML and JSON report per sample |
| 2 | `2_Salmon.sh` | Builds the Salmon index (if missing) and quantifies each sample |
| 3.1 | `3.1_make_tx2gene.R` | Builds the transcript-to-gene table from the GFF |
| 3.2 | `3.2_tximport_combined.R` | Imports Salmon output with tximport and writes gene count and TPM tables |
| 4 | `4_Sample_Relationship_Analysis.py` | Combined figure: PCA, sample correlation heat map and expression violin plot |
| 5 | `5_DESeq2_Analysis.R` | Differential expression with DESeq2 for a list of comparisons, with volcano, MA, heat map, bar and Venn figures |
| 6.1 | `6.1_GeneOntology.R` | GO enrichment (topGO) of the up- and down-regulated genes of each comparison, with a functional annotation figure and a functional enrichment figure |
| 6.2 | `6.2_KEGG_Analysis.R` | KEGG pathway annotation and enrichment of the same DEGs, with the same two figures |
| 7 | `7_TF_Heatmap.R` | Transcription factor heatmap |
| 8 | `8_WGCNA.R` | WGCNA co-expression analysis |

Helper scripts: `fastp_summary_from_json_files.ipynb`, `Fastp_HTML_to_PowerPoint.ipynb`, `collect_salmon_mapping_rate.py`, `submit_R.sh`.

Steps 7 and 8 and the helper scripts are being added to this repository one at a time as they are finalized. Steps 1 to 6.2 are documented in full below.

## Requirements

- An LSF cluster with `bsub`, and conda
- A conda environment containing `fastp` and `salmon`
- Python 3 with `pandas`, `numpy`, `matplotlib`, `seaborn` and `scikit-learn`
- R with these packages: `rtracklayer`, `dplyr`, `readr`, `stringr`, `tximport`
- For step 5, also: `DESeq2` (Bioconductor), `ggplot2`, `ggrepel`, `pheatmap`, `RColorBrewer`, `ggVennDiagram`, `extrafont` and `patchwork`. The script installs any that are missing.
- For step 6.1, also: `topGO` (Bioconductor), `ggplot2`, `dplyr`, `tidyr`, `stringr`, `readr`, `purrr`, `tibble` and `ggtext`. The script installs any that are missing.
- For step 6.2, also: `clusterProfiler` and `KEGGREST` (Bioconductor), `jsonlite`, `ragg` and the packages used in step 6.1 except `topGO`. The script installs any that are missing, and needs an internet connection the first time it runs (it downloads the KEGG pathway tables).

## Folder structure

Create one project folder named `Bulk_RNA_Seq`. Everything (scripts, references, reads, results, sample information) lives inside it.

```
Bulk_RNA_Seq/
├── 1_Codes/
│   ├── 1_Bulk_RNA_Seq_fastp.sh
│   ├── 2_Salmon.sh
│   ├── 3.1_make_tx2gene.R
│   ├── 3.2_tximport_combined.R
│   ├── 4_sample_relationship_analysis.py
│   ├── 5_DESeq2_Analysis.R
│   ├── 6.1_GeneOntology.R
│   ├── 6.2_KEGG_Analysis.R
│   ├── 7_TF_Heatmap.R
│   ├── 8_WGCNA.R
│   ├── fastp_summary_from_json_files.ipynb
│   ├── Fastp_HTML_to_PowerPoint.ipynb
│   ├── collect_salmon_mapping_rate.py
│   └── submit_R.sh
├── 2_References/
│   ├── gene_annotation.gff
│   ├── ITAG4.0_cDNA.fasta
│   ├── ITAG4.0_goterms.txt
│   ├── go-basic.obo
│   └── query.ko.txt
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
│   └── 8_WGCNA/
└── Metadata.csv
```

To create all the folders at once, run this from the directory where you want the project (add `Metadata.csv` yourself):

```bash
mkdir -p Bulk_RNA_Seq/{1_Codes,2_References,3_Raw_Reads,4_Trimmed_Reads} \
         Bulk_RNA_Seq/5_Results/{1_fastp,2_Salmon,3_tximport,4_Sample_Relationship_Analysis,5_DESeq2,6.1_GO,6.2_KEGG,7_TF,8_WGCNA}
```

## Before you run anything

The scripts contain placeholder paths written as `/.../.../`. Replace them with your own paths in every script:

- `INPUTDIR`, `OUTPUTDIR`, `RESULTDIR`, `seq_path`
- the conda environment path in `conda activate ...`
- the GFF, tx2gene and `quant.sf` paths in the R scripts
- `counts_file`, `meta_file` and `output_dir` in `5_DESeq2_Analysis.R` (or use `file.choose()` for the two input files, see Step 5)
- `go_file`, `obo_file`, `deseq_dir` and `output_dir` in `6.1_GeneOntology.R`, and `ko_file`, `deseq_dir` and `output_dir` in `6.2_KEGG_Analysis.R` (or use `file.choose()` for the input files, see Step 6)

The GFF (`gene_annotation.gff`) and the transcriptome fasta must come from the same ITAG4.0 release. If they don't, transcript IDs in Salmon output won't match the tx2gene table.

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

```bash
cd 1_Codes
bsub < 1_Bulk_RNA_Seq_fastp.sh
bjobs            # check job status
```

Reads `3_Raw_Reads/`. Writes:

- `4_Trimmed_Reads/TMD_<sample>_1.fq.gz` and `TMD_<sample>_2.fq.gz`
- `5_Results/1_fastp/<sample>.html` and `<sample>.json`

## Step 2: Quantify transcripts (Salmon)

```bash
bsub < 2_Salmon.sh
bjobs
```

Builds the index in `2_References/salmon_tmt_index` if it does not exist yet, then quantifies every `TMD_*_1.fq.gz` pair in `4_Trimmed_Reads/`. Writes one folder per sample in `5_Results/2_Salmon/`, each containing `quant.sf`.

Set `transcriptome` in the script to your reference transcriptome fasta.

## Step 3.1: Transcript-to-gene table

Run `3.1_make_tx2gene.R` in RStudio or with `Rscript`. It keeps the `mRNA` features from the GFF and saves `2_References/tx2gene.csv` with two columns, `TXNAME` and `GENEID`. The `mRNA:` and `gene:` prefixes are removed if present.

The script ends with a check. Point it at any sample's `quant.sf`; the printed fraction should be close to 1. A low value means the transcript IDs in the GFF and the fasta do not match.

## Step 3.2: Count matrix (tximport)

Run `3.2_tximport_combined.R`. It reads every `quant.sf` under `5_Results/2_Salmon/`, sorts the samples in natural order (`2F` before `10F`), and summarizes to genes with `countsFromAbundance = "lengthScaledTPM"`.

Written to `5_Results/3_tximport/`:

- `txi.rds`: the full tximport object, for DESeq2 later
- `gene_counts.csv`
- `TPM.csv`

Each sample folder in `5_Results/2_Salmon/` must have a unique name, because the folder name becomes the column name.

## Step 4: Sample relationships

```bash
python 4_sample_relationship_analysis.py
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

The figure is saved at 1000 dpi, and with many samples the PNG is large. Lower the `dpi` in the two `plt.savefig` calls if it runs out of memory.

## Step 5: Differential expression (DESeq2)

Run `5_DESeq2_Analysis.R` in R (RStudio or `Rscript`). It runs DESeq2 separately for every comparison in the `comparisons` list, then makes tables and figures for each comparison and for all comparisons together.

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
2. **`file.choose()`.** A window opens so you can pick the file. Remove the `#` from the `message(...)` and `file.choose()` lines for that file and put a `#` in front of the direct-path line. This needs an interactive R session such as RStudio.

Set `output_dir` to your `5_Results/5_DESeq2` folder. The script creates the `Tables` and `Figures` subfolders it needs.

**Parameters** (section 3 of the script): adjusted p-value cutoff (`padj_cutoff = 0.05`), log2 fold-change cutoff (`lfc_cutoff = 2`), the minimum total reads for a gene to be kept in a comparison (`min_gene_count = 10`), and how many genes are labeled on volcano plots and shown in heat maps.

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
- `Figures/Volcano_Plots/`, `Figures/MA_Plots/` and `Figures/Heatmaps/`: one figure per comparison (heat maps are skipped when a comparison has fewer than 2 DEGs)
- `Figures/`: overall sample PCA, heat map of the most variable genes, bar plot of up/down gene counts, Venn diagram, and the combined MA and volcano figures

All figures are TIFF files saved at 1000 dpi.

**Fonts.** Figures use Times New Roman through `extrafont`. The first time you use `extrafont` on a computer, run `font_import()` once (it takes a few minutes). The script loads fonts with `loadfonts(device = "win")`, which is for Windows; on Mac or Linux use `device = "pdf"` and set `FONT` to a font you have (for example `"serif"`).

## Step 6: GO and KEGG enrichment

Both scripts read the DESeq2 result tables written by Step 5 (`5_Results/5_DESeq2/Tables/Comparisons/<name>_DESeq2_results.csv`). A gene is a DEG when `padj` is below `PADJ_CUTOFF` (0.05) and its absolute log2 fold change is above `LFC_CUTOFF` (1), set in section 2 of each script. Step 5 calls Up and Down genes with `lfc_cutoff = 2`, so set `LFC_CUTOFF` to 2 in Step 6 if you want the same DEGs.

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

Run it in R (RStudio or `Rscript`). For each comparison, the up- and down-regulated genes are tested separately for each GO category (MF, CC, BP) with topGO (`weight01` algorithm, Fisher test, minimum node size 5). The background is all genes in the DESeq2 table that have a GO annotation. p-values are adjusted with Benjamini-Hochberg and called significant at FDR 0.05.

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

Run it in R. The script attaches K numbers from `query.ko.txt` to the DEGs, maps them to KEGG pathways, and keeps the five main KEGG categories (Metabolism, Genetic Information Processing, Environmental Information Processing, Cellular Processes, Organismal Systems). It downloads the pathway tables from KEGG on the first run and saves them in the output folder, so later runs work from the saved copies.

Gene IDs in `query.ko.txt` and in the DESeq2 tables can differ by a trailing transcript number (`Solyc10g079470.3.1` and `Solyc10g079470.3`). The script tries a few ways of trimming the IDs, keeps the one that matches the DESeq2 IDs best, and prints the match rate. Check that number the first time you run it.

Figure 2 uses fold enrichment, calculated directly as (DEGs in the pathway / all DEGs) / (background genes in the pathway / all background genes). A statistical `enrichKEGG()` test is also run and saved, but the figure does not depend on it, so every selected pathway is always shown.

Written to `5_Results/6.2_KEGG/`:

- `01_KO_annotation/`: gene to K number table
- `02_Annotated_DESeq2/`: each DESeq2 table with K numbers and KEGG categories added
- `03_DEG_tables/`: all, up and down DEGs for each comparison
- `04_KEGG_enrichment/`: fold enrichment for all comparisons and the `enrichKEGG()` result for each comparison
- `Pathway_Ranking_AllComparisons.csv`, `Top_Pathways_Figure1.csv` and `Top_Pathways_Figure2.csv`
- `Comparison_Key.csv`: which letter is which comparison
- `Figure1_Functional_Annotation.tiff`: top 30 pathways with up and down DEG counts
- `Figure2_Functional_Enrichment.tiff`: top 20 pathways by fold enrichment

The number of pathways shown (`TOP_N_ANNOTATION`, `TOP_N_ENRICHMENT`), colors, fonts and figure sizes are set in section 2.

## Single-end reads

Both shell scripts default to paired-end. For single-end data:

**`1_Bulk_RNA_Seq_fastp.sh`**

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
