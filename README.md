# NGS Analysis Tutorial

A practical guide to NGS analysis: CUT&RUN, ChIP-seq, ATAC-seq, RNA-seq, single-cell, spatial and related assays, from raw reads to publication figures.

Every chapter gives parameterised commands, example outputs, Nature-style example figures and a list of references and resources.

> **About the figures.** All figures in this tutorial are generated from **simulated data** by [`scripts/make_example_figures.py`](scripts/make_example_figures.py). They show what each step produces and how to lay it out for publication; they are not results from any experiment.

---

## Workflow

```text
             ┌──────────────┐   ┌────────────┐   ┌───────────────────┐   ┌──────────────┐
 BaseSpace ─►│ 1 Data in/out│──►│ 2 QC / trim│──►│ 3 Align → BAM → bw│──►│ 4 Peak calls │
 SRA / GEO   └──────────────┘   └────────────┘   └───────────────────┘   └──────┬───────┘
                                                                                │ union peaks
         ┌──────────────────┐   ┌────────────────────┐   ┌──────────────────────▼──────┐
         │ 8 Tracks/heatmaps│◄──│ 7 Motifs / networks│◄──│ 5 Counts → DESeq2 → plots   │
         └──────────────────┘   └────────────────────┘   │ 6 Annotation / GO / GREAT   │
                                                         └─────────────────────────────┘
 RNA-seq ─► 9 STAR → featureCounts → DESeq2 → clustered heatmap → GSEA
 Nuclear organisation ─► 10 SON TSA-seq · 11 Hi-C / HiChIP
 Screens              ─► 12 pooled CRISPR (MAGeCK)
 Single-cell          ─► 13 scRNA-seq · 14 scATAC-seq (chromVAR, tracks) · 15 spatial
```

## Chapters

Click a figure to open its chapter.

| # | Chapter | Covers | Example figure |
|---|---|---|:-:|
| 0 | [Setup](tutorials/00_setup.md) | LSF `bsub`, project layout, reference files, modules, conda, tool installs | — |
| 1 | [Data transfer](tutorials/01_data_transfer.md) | BaseSpace, `bcl2fastq`, SRA, Box, Dropbox, GEO upload, down-sampling | — |
| 2 | [QC & trimming](tutorials/02_qc_trimming.md) | FastQC/MultiQC, Trimmomatic, cutadapt, insert size, strandedness, FRiP | <a href="tutorials/02_qc_trimming.md"><img src="tutorials/figures/fig1_qc.png" width="260" alt="Fig. 1"></a><br><sub><b>Fig. 1</b> Library QC</sub> |
| 3 | [CUT&RUN / ChIP alignment](tutorials/03_alignment_cutrun_chip.md) | bowtie2, filtering, Picard dedup, bigWig, spike-in, parallel jobs | — |
| 4 | [Peak calling](tutorials/04_peak_calling.md) | MACS2 broad/IgG, filtering, consensus & union peaks, overlaps | <a href="tutorials/04_peak_calling.md"><img src="tutorials/figures/fig2_peaks.png" width="260" alt="Fig. 2"></a><br><sub><b>Fig. 2</b> Peaks</sub> |
| 5 | [Differential analysis](tutorials/05_differential_analysis.md) | `multicov`/featureCounts, DESeq2, PCA, correlation, MA, volcano, box plots | <a href="tutorials/05_differential_analysis.md"><img src="tutorials/figures/fig4_differential.png" width="260" alt="Fig. 4"></a><br><sub><b>Fig. 4</b> Differential</sub> |
| 6 | [Annotation & enrichment](tutorials/06_annotation_enrichment.md) | HOMER annotation, GREAT, REVIGO, biomaRt, clusterProfiler, Enrichr, Venn | <a href="tutorials/06_annotation_enrichment.md"><img src="tutorials/figures/fig5_annotation_boxplot.png" width="260" alt="Fig. 5"></a><br><sub><b>Fig. 5</b> Annotation</sub> |
| 7 | [Motif analysis](tutorials/07_motif_analysis.md) | HOMER, MEME AME/FIMO, logos, motif centrality, HINT footprinting, Cytoscape | <a href="tutorials/07_motif_analysis.md"><img src="tutorials/figures/fig6_motif_enrichment.png" width="260" alt="Fig. 6"></a><br><sub><b>Fig. 6</b> Motifs & GO</sub><br><a href="tutorials/07_motif_analysis.md"><img src="tutorials/figures/fig10_motif_results.png" width="260" alt="Fig. 10"></a><br><sub><b>Fig. 10</b> Motif results</sub> |
| 8 | [Visualization](tutorials/08_visualization.md) | deepTools heatmaps, pooled bigWigs, UCSC tracks, Gviz | <a href="tutorials/08_visualization.md"><img src="tutorials/figures/fig3_profile_heatmap.png" width="260" alt="Fig. 3"></a><br><sub><b>Fig. 3</b> Profile heatmap</sub><br><a href="tutorials/08_visualization.md"><img src="tutorials/figures/fig8_genome_tracks.png" width="260" alt="Fig. 8"></a><br><sub><b>Fig. 8</b> Tracks</sub> |
| 9 | [RNA-seq](tutorials/09_rnaseq.md) | STAR, featureCounts, count cleaning, ComplexHeatmap k-means, GSEA (GSEA-CLI, fgsea, running-ES plots), SQuIRE | <a href="tutorials/09_rnaseq.md"><img src="tutorials/figures/fig7_rna_kmeans_heatmap.png" width="260" alt="Fig. 7"></a><br><sub><b>Fig. 7</b> RNA k-means</sub> |
| 10 | [SON TSA-seq](tutorials/10_tsa_seq.md) | sliding-window normalisation, speckle proximity vs expression | <a href="tutorials/10_tsa_seq.md"><img src="tutorials/figures/fig9_tsa_seq.png" width="260" alt="Fig. 9"></a><br><sub><b>Fig. 9</b> TSA-seq</sub> |
| 11 | [Hi-C / HiChIP](tutorials/11_hic.md) | HiC-Pro, HiCExplorer (P(s), TADs, insulation, track plots), FitHiChIP, HOMER | <a href="tutorials/11_hic.md"><img src="tutorials/figures/fig12_hic.png" width="260" alt="Fig. 12"></a><br><sub><b>Fig. 12</b> Hi-C</sub> |
| 12 | [CRISPR screens](tutorials/12_crispr_screens.md) | MAGeCK count/test, merging per-sample count tables, library QC, rank plot, sgRNA plots, RIGER, CRISPRO | <a href="tutorials/12_crispr_screens.md"><img src="tutorials/figures/fig11_crispr_screen.png" width="260" alt="Fig. 11"></a><br><sub><b>Fig. 11</b> CRISPR screen</sub> |
| 13 | [scRNA-seq](tutorials/13_scrnaseq.md) | Cell Ranger, QC, doublets, Harmony, clustering, markers, annotation, pseudobulk DE, SCENIC | <a href="tutorials/13_scrnaseq.md"><img src="tutorials/figures/fig13_scrnaseq.png" width="260" alt="Fig. 13"></a><br><sub><b>Fig. 13</b> scRNA-seq</sub> |
| 14 | [scATAC-seq](tutorials/14_scatacseq.md) | Signac QC, LSI, label transfer, cell-type peaks, chromVAR, footprints, coverage tracks, ArchR | <a href="tutorials/14_scatacseq.md"><img src="tutorials/figures/fig14_scatac.png" width="260" alt="Fig. 14"></a><br><sub><b>Fig. 14</b> scATAC-seq</sub> |
| 15 | [Spatial transcriptomics](tutorials/15_spatial.md) | Scanpy + Harmony integration, clusters on tissue, Squidpy neighbourhoods | <a href="tutorials/15_spatial.md"><img src="tutorials/figures/fig15_spatial.png" width="260" alt="Fig. 15"></a><br><sub><b>Fig. 15</b> Spatial</sub> |

## Figure style (Nature guidelines)

| Element | Setting |
|---|---|
| Width | 89 mm (single column) or 183 mm (double column) |
| Font | Arial/Helvetica, 5–7 pt; bold lowercase panel letters (**a**, **b**, …) |
| Lines | 0.25–0.5 pt axes, outward ticks, no top/right frame |
| Colour | colour-blind-safe Okabe–Ito palette; red–blue diverging for z-scores |
| Export | vector PDF with embedded fonts (`pdf.fonttype 42` / `cairo_pdf`); dense scatter/heatmap layers rasterized at 300 dpi |
| Legends | state *n*, the statistical test, the error-bar definition and the thresholds |

R users: `theme_nature()` and `save_fig()` are defined in [chapter 5](tutorials/05_differential_analysis.md#53-plot-functions-nature-style). Python users: copy the `rcParams` block at the top of `scripts/make_example_figures.py`.

Regenerate every figure (PNG + PDF):

```bash
python3 scripts/make_example_figures.py      # needs numpy, scipy, matplotlib
```

---

## Common pitfalls covered in this tutorial

| Where | Pitfall | Recommendation |
|---|---|---|
| MACS2 on paired-end CUT&RUN (ch. 4) | `-f BAM` ignores mate information | `-f BAMPE` |
| Spatial / single-cell integration (ch. 13, 15) | `neighbors()` ran on uncorrected PCs; Harmony block skipped normalisation | `use_rep="X_pca_harmony"`; normalise once before PCA |
| HOMER motif plots (ch. 7) | `log10(exp(LogP))` underflows to −∞ | `LogP / log(10)` |
| rGREAT (ch. 6) | `adv_oneDistance` is in kb (5000 = 5 Mb) | documented; default 1000 kb |
| DESeq2 peaks (ch. 5) | factor levels follow file order | list the reference condition first in the sample sheet |
| bowtie2 logging (ch. 3) | `2>` after a pipe captures samtools, not bowtie2, stderr | redirect on the bowtie2 call |
| bigWig via bedGraph (ch. 3) | unsorted bedGraph passed to `bedGraphToBigWig` | sort first |
| Pooling bigWigs (ch. 8) | summing replicates changes the scale | `bigwigAverage` / `wiggletools mean` |
| scATAC barcode renaming (ch. 14) | global substitution edits the whole SAM line | edit the `CB:Z:` tag only |
| featureCounts ≥ 2.0.2 (ch. 5, 9) | `-p` alone counts reads | add `--countReadPairs` |
| Credentials (ch. 1) | tokens/passwords written into scripts | keep them in `~/.secrets` (chmod 600) |
| scRNA-seq DE (ch. 13) | cells treated as replicates inflate significance | pseudobulk per sample + DESeq2 |
| scATAC LSI (ch. 14) | LSI component 1 tracks sequencing depth | use components 2–30 |
| CRISPR hits (ch. 12) | a single extreme guide | require several concordant guides outside the non-targeting range |
| MAGeCK input (ch. 12) | pre-normalised, rounded counts passed to `mageck test` | give raw counts; let MAGeCK normalise (`--norm-method control`) |
| GSEA plots (ch. 9) | running ES drawn by joining only the hit positions from the GSEA table | rebuild the full curve (`es_from_gsea()`), which drops between hits |
