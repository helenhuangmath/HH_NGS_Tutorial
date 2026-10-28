# NGS Analysis Tutorial

Useful functions to NGS analysis: CUT&RUN, ChIP-seq, ATAC-seq, RNA-seq, single-cell, spatial and related assays, from raw reads to publication figures.

Every part gives parameterised commands, example outputs, Nature-style example figures directly after the relevant functions, and a list of references and resources. Figures are shown inline in each part; click a figure label to open its vector PDF.

---

## Workflow

```text
FASTQ ─► QC / trimming ─► alignment ─► signal tracks + peaks or gene counts
      ─► differential analysis ─► annotation, enrichment, motifs ─► figures
```

Parts 0–8 follow this path for CUT&RUN, ChIP-seq and ATAC-seq; parts 9–15 adapt it to other assays; parts 16–18 are reusable R figure recipes.

## Parts

Open a part for functions and commands, or a figure link for its PDF.

| # | Part | Covers | Example figure |
|---|---|---|:-:|
| 0 | [Setup](tutorials/00_setup.md) | LSF `bsub`, project layout, reference files, modules, conda, tool installs | — |
| 1 | [Data transfer](tutorials/01_data_transfer.md) | BaseSpace, `bcl2fastq`, SRA, Box, Dropbox, GEO upload, down-sampling | — |
| 2 | [QC & trimming](tutorials/02_qc_trimming.md) | FastQC/MultiQC, Trimmomatic, cutadapt, insert size, strandedness, FRiP | [Fig. 2 — Library QC (PDF)](tutorials/figures/part02_qc.pdf) |
| 3 | [CUT&RUN / ChIP alignment](tutorials/03_alignment_cutrun_chip.md) | bowtie2, filtering, Picard dedup, bigWig, spike-in, parallel jobs | — |
| 4 | [Peak calling](tutorials/04_peak_calling.md) | MACS2 broad/IgG, filtering, consensus & union peaks, overlaps | [Fig. 4 — Peaks (PDF)](tutorials/figures/part04_peaks.pdf) |
| 5 | [Differential analysis](tutorials/05_differential_analysis.md) | `multicov`/featureCounts, DESeq2, PCA, correlation, MA, volcano, box plots | [Fig. 5 — Differential (PDF)](tutorials/figures/part05_differential.pdf) |
| 6 | [Annotation & enrichment](tutorials/06_annotation_enrichment.md) | HOMER annotation, GREAT, REVIGO, biomaRt, clusterProfiler, Enrichr, Venn | [Fig. 6 — Annotation (PDF)](tutorials/figures/part06_annotation_boxplot.pdf) |
| 7 | [Motif analysis](tutorials/07_motif_analysis.md) | HOMER, MEME AME/FIMO, logos, motif centrality, HINT footprinting, Cytoscape | [Fig. 7.1 — Motifs & GO (PDF)](tutorials/figures/part07_motif_enrichment.pdf) · [Fig. 7.2 — Motif results (PDF)](tutorials/figures/part07_motif_results.pdf) |
| 8 | [Visualization](tutorials/08_visualization.md) | deepTools heatmaps, pooled bigWigs, UCSC tracks, Gviz | [Fig. 8.1 — Profile heatmap (PDF)](tutorials/figures/part08_profile_heatmap.pdf) · [Fig. 8.2 — Tracks (PDF)](tutorials/figures/part08_genome_tracks.pdf) |
| 9 | [RNA-seq](tutorials/09_rnaseq.md) | STAR, featureCounts, count cleaning, ComplexHeatmap k-means, GSEA (GSEA-CLI, fgsea, running-ES plots), SQuIRE | [Fig. 9 — RNA k-means (PDF)](tutorials/figures/part09_rna_kmeans_heatmap.pdf) |
| 10 | [SON TSA-seq](tutorials/10_tsa_seq.md) | sliding-window normalisation, speckle proximity vs expression | [Fig. 10 — TSA-seq (PDF)](tutorials/figures/part10_tsa_seq.pdf) |
| 11 | [Hi-C / HiChIP](tutorials/11_hic.md) | HiC-Pro, HiCExplorer (P(s), TADs, insulation, track plots), FitHiChIP, HOMER | [Fig. 11 — Hi-C (PDF)](tutorials/figures/part11_hic.pdf) |
| 12 | [CRISPR screens](tutorials/12_crispr_screens.md) | MAGeCK count/test, merging per-sample count tables, library QC, rank plot, sgRNA plots, RIGER, CRISPRO | [Fig. 12 — CRISPR screen (PDF)](tutorials/figures/part12_crispr_screen.pdf) |
| 13 | [scRNA-seq](tutorials/13_scrnaseq.md) | Cell Ranger, QC, doublets, Harmony, clustering, markers, annotation, pseudobulk DE, SCENIC | [Fig. 13 — scRNA-seq (PDF)](tutorials/figures/part13_scrnaseq.pdf) |
| 14 | [scATAC-seq](tutorials/14_scatacseq.md) | Signac QC, LSI, label transfer, cell-type peaks, chromVAR, footprints, coverage tracks, ArchR | [Fig. 14 — scATAC-seq (PDF)](tutorials/figures/part14_scatac.pdf) |
| 15 | [Spatial transcriptomics](tutorials/15_spatial.md) | Scanpy + Harmony integration, clusters on tissue, Squidpy neighbourhoods | [Fig. 15 — Spatial (PDF)](tutorials/figures/part15_spatial.pdf) |
| 16 | [PCA plots](tutorials/16_pca_plot.md) | rlog/VST input, colour by condition + shape by guide, PC1/2 and PC2/3 | [Fig. 16 — PCA (PDF)](tutorials/figures/part16_pca.pdf) |
| 17 | [ComplexHeatmap](tutorials/17_complexheatmap.md) | row z-scores, column split by condition, labelled genes (`anno_mark`) | [Fig. 17 — Heatmap (PDF)](tutorials/figures/part17_complexheatmap.pdf) |
| 18 | [Gviz highlighted tracks](tutorials/18_gviz_tracks.md) | bigWig tracks per condition, TxDb gene models, highlighted DERs and motif sites, loop over regions | [Fig. 18 — Tracks (PDF)](tutorials/figures/part18_gviz_highlight.pdf) |

## Figure style (Nature guidelines)

| Element | Setting |
|---|---|
| Width | 89 mm (single column) or 183 mm (double column) |
| Font | Arial/Helvetica, 5–7 pt; bold lowercase panel letters (**a**, **b**, …) |
| Lines | 0.25–0.5 pt axes, outward ticks, no top/right frame |
| Colour | colour-blind-safe Okabe–Ito palette; red–blue diverging for z-scores |
| Export | vector PDF with embedded fonts (`pdf.fonttype 42` / `cairo_pdf`); dense scatter/heatmap layers rasterized at 300 dpi |
| Legends | state *n*, the statistical test, the error-bar definition and the thresholds |

R users: `theme_nature()` and `save_fig()` are defined in [part 5](tutorials/05_differential_analysis.md#53-plot-functions-nature-style). Python users: copy the `rcParams` block at the top of `scripts/make_example_figures.py`.

Regenerate the figures of parts 2–15 (PDF, plus the PNG shown on each page):

```bash
python3 scripts/make_example_figures.py      # needs numpy, scipy, matplotlib
```

---

## Common pitfalls covered in this tutorial

| Where | Pitfall | Recommendation |
|---|---|---|
| MACS2 on paired-end CUT&RUN (part 4) | `-f BAM` ignores mate information | `-f BAMPE` |
| Spatial / single-cell integration (parts 13, 15) | `neighbors()` ran on uncorrected PCs; Harmony block skipped normalisation | `use_rep="X_pca_harmony"`; normalise once before PCA |
| HOMER motif plots (part 7) | `log10(exp(LogP))` underflows to −∞ | `LogP / log(10)` |
| rGREAT (part 6) | `adv_oneDistance` is in kb (5000 = 5 Mb) | documented; default 1000 kb |
| DESeq2 peaks (part 5) | factor levels follow file order | list the reference condition first in the sample sheet |
| bowtie2 logging (part 3) | `2>` after a pipe captures samtools, not bowtie2, stderr | redirect on the bowtie2 call |
| bigWig via bedGraph (part 3) | unsorted bedGraph passed to `bedGraphToBigWig` | sort first |
| Pooling bigWigs (part 8) | summing replicates changes the scale | `bigwigAverage` / `wiggletools mean` |
| scATAC barcode renaming (part 14) | global substitution edits the whole SAM line | edit the `CB:Z:` tag only |
| featureCounts ≥ 2.0.2 (parts 5, 9) | `-p` alone counts reads | add `--countReadPairs` |
| Credentials (part 1) | tokens/passwords written into scripts | keep them in `~/.secrets` (chmod 600) |
| scRNA-seq DE (part 13) | cells treated as replicates inflate significance | pseudobulk per sample + DESeq2 |
| scATAC LSI (part 14) | LSI component 1 tracks sequencing depth | use components 2–30 |
| CRISPR hits (part 12) | a single extreme guide | require several concordant guides outside the non-targeting range |
| MAGeCK input (part 12) | pre-normalised, rounded counts passed to `mageck test` | give raw counts; let MAGeCK normalise (`--norm-method control`) |
| GSEA plots (part 9) | running ES drawn by joining only the hit positions from the GSEA table | rebuild the full curve (`es_from_gsea()`), which drops between hits |
