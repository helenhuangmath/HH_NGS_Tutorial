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

Open a part for functions and commands; click a figure to enlarge it.

| # | Part | Covers | Example figure |
|---|---|---|:-:|
| 0 | [Setup](tutorials/00_setup.md) | LSF `bsub`, project layout, reference files, modules, conda, tool installs | — |
| 1 | [Data transfer](tutorials/01_data_transfer.md) | BaseSpace, `bcl2fastq`, SRA, Box, Dropbox, GEO upload, down-sampling | — |
| 2 | [QC & trimming](tutorials/02_qc_trimming.md) | FastQC/MultiQC, Trimmomatic, cutadapt, insert size, strandedness, FRiP | <img src="tutorials/figures/overview/part02_qc.png" alt="Fig. 2" width="220"> |
| 3 | [CUT&RUN / ChIP alignment](tutorials/03_alignment_cutrun_chip.md) | bowtie2, filtering, Picard dedup, bigWig, spike-in, parallel jobs | — |
| 4 | [Peak calling](tutorials/04_peak_calling.md) | MACS2 broad/IgG, filtering, consensus & union peaks, overlaps | <img src="tutorials/figures/overview/part04_peaks.png" alt="Fig. 4" width="220"> |
| 5 | [Differential analysis](tutorials/05_differential_analysis.md) | `multicov`/featureCounts, DESeq2, PCA, correlation, MA, volcano, box plots | <img src="tutorials/figures/overview/part05_differential.png" alt="Fig. 5" width="220"> |
| 6 | [Annotation & enrichment](tutorials/06_annotation_enrichment.md) | HOMER annotation, GREAT, REVIGO, biomaRt, clusterProfiler, Enrichr, Venn | <img src="tutorials/figures/overview/part06_annotation.png" alt="Fig. 6" width="220"> |
| 7 | [Motif analysis](tutorials/07_motif_analysis.md) | HOMER, MEME AME/FIMO, logos, motif centrality, HINT footprinting, Cytoscape | <img src="tutorials/figures/overview/part07_motif_enrichment.png" alt="Fig. 7" width="220"> <img src="tutorials/figures/overview/part07_motif_results.png" alt="Fig. 7" width="220"> |
| 8 | [Visualization](tutorials/08_visualization.md) | deepTools heatmaps, pooled bigWigs, UCSC tracks, Gviz | <img src="tutorials/figures/overview/part08_profile_heatmap.png" alt="Fig. 8" width="220"> <img src="tutorials/figures/part08_genome_tracks.png" alt="Fig. 8" width="220"> |
| 9 | [RNA-seq](tutorials/09_rnaseq.md) | STAR, featureCounts, count cleaning, ComplexHeatmap k-means, GSEA (GSEA-CLI, fgsea, running-ES plots), SQuIRE | <img src="tutorials/figures/part09_kmeans_heatmap.png" alt="Fig. 9" width="140"> |
| 10 | [SON TSA-seq](tutorials/10_tsa_seq.md) | sliding-window normalisation, speckle proximity vs expression | <img src="tutorials/figures/overview/part10_tsa_seq.png" alt="Fig. 10" width="220"> |
| 11 | [Hi-C / HiChIP](tutorials/11_hic.md) | HiC-Pro, HiCExplorer (P(s), TADs, insulation, track plots), FitHiChIP, HOMER | <img src="tutorials/figures/overview/part11_hic.png" alt="Fig. 11" width="220"> |
| 12 | [CRISPR screens](tutorials/12_crispr_screens.md) | MAGeCK count/test, merging per-sample count tables, library QC, rank plot, sgRNA plots, RIGER, CRISPRO | <img src="tutorials/figures/overview/part12_crispr_screen.png" alt="Fig. 12" width="220"> |
| 13 | [scRNA-seq](tutorials/13_scrnaseq.md) | Cell Ranger, QC, doublets, Harmony, clustering, markers, annotation, pseudobulk DE, SCENIC | <img src="tutorials/figures/overview/part13_scrnaseq.png" alt="Fig. 13" width="220"> |
| 14 | [scATAC-seq](tutorials/14_scatacseq.md) | Signac QC, LSI, label transfer, cell-type peaks, chromVAR, footprints, coverage tracks, ArchR | <img src="tutorials/figures/overview/part14_scatac.png" alt="Fig. 14" width="220"> |
| 15 | [Spatial transcriptomics](tutorials/15_spatial.md) | Scanpy + Harmony integration, clusters on tissue, Squidpy neighbourhoods | <img src="tutorials/figures/overview/part15_spatial.png" alt="Fig. 15" width="220"> |
| 16 | [PCA plots](tutorials/16_pca_plot.md) | rlog/VST input, colour by condition + shape by guide, PC1/2 and PC2/3 | <img src="tutorials/figures/part16_pca.png" alt="Fig. 16" width="220"> |
| 17 | [ComplexHeatmap](tutorials/17_complexheatmap.md) | row z-scores, column split by condition, labelled genes (`anno_mark`) | <img src="tutorials/figures/part17_complexheatmap.png" alt="Fig. 17" width="90"> |
| 18 | [Gviz highlighted tracks](tutorials/18_gviz_tracks.md) | bigWig tracks per condition, TxDb gene models, highlighted DERs and motif sites, loop over regions | <img src="tutorials/figures/part18_gviz_highlight.png" alt="Fig. 18" width="220"> |

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

Regenerate the figures of parts 2–15 (PDF, plus the PNGs shown here and on each page):

```bash
python3 scripts/make_example_figures.py      # needs numpy, scipy, matplotlib
```
