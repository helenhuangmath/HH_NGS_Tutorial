# Part 18 · Genome tracks with highlighted regions (Gviz)

[← ComplexHeatmap](17_complexheatmap.md) · [Index](../README.md)

Builds on the basic Gviz figure in [part 8.4](08_visualization.md#84-publication-track-figures-in-r-gviz): one track per condition, peaks, gene models from a TxDb, and **highlight boxes** over differential regions (grey) and motif sites (red), looped over a table of regions.

---

## 18.1 Inputs

```r
library(Gviz); library(rtracklayer); library(data.table); library(dplyr)
library(TxDb.Mmusculus.UCSC.mm10.knownGene); library(org.Mm.eg.db)

gen <- "mm10"
bw  <- c(cond1 = "Wig/pool_cond1.rpm.ave.bw",
         cond2 = "Wig/pool_cond2.rpm.ave.bw",
         cond3 = "Wig/pool_cond3.rpm.ave.bw")
cols <- c(cond1 = "#b8b8b8", cond2 = "#2c7fb8", cond3 = "#a1dab4")

peaks   <- import("Peak/peak_Union.bed")
DER     <- fread("deseq2_AllRes_All3Compares_Sig_padj0.05_FC1.5.txt")   # Chr Start End …
motif   <- fread("fimo_TCF7.Merge_wDER_padj0.05_FC1.5_noNA.txt")        # Chr motif_start motif_stop …
regions <- fread("regions_toplot.txt")                                  # Gene Chr Start End ymax
```

---

## 18.2 Tracks

```r
axis_track <- GenomeAxisTrack(scale = 0.2, labelPos = "above")      # 20%-width scale bar
peak_track <- AnnotationTrack(peaks, genome = gen, name = "", fill = "grey20", col = "grey20")

## gene models: longest transcript, Entrez IDs → symbols
gene_track <- GeneRegionTrack(TxDb.Mmusculus.UCSC.mm10.knownGene, genome = gen, name = "",
                              collapseTranscripts = "longest", transcriptAnnotation = "symbol",
                              fill = "black", col = "black", protein_coding = "black",
                              utr5 = "black", utr3 = "black", fontcolor.group = "black")
has_id <- !is.na(gene(gene_track))
symbol(gene_track)[has_id] <- mapIds(org.Mm.eg.db, keys = sub("\\.\\d+$", "", gene(gene_track)[has_id]),
                                     keytype = "ENTREZID", column = "SYMBOL")

signal_tracks <- function(chr) lapply(names(bw), function(s)
  DataTrack(range = bw[[s]], genome = gen, chromosome = chr, name = s, type = "histogram",
            fill.histogram = cols[[s]], col.histogram = cols[[s]],
            baseline = 0, col.baseline = cols[[s]], col.axis = "black", fontcolor.title = "black"))

## DERs (grey) and motif sites (red) inside the window
highlight_sites <- function(chr, from, to) bind_rows(
  DER   |> filter(Chr == chr, Start >= from, End <= to) |>
    transmute(start = Start, end = End, fill = "#e8e8e8", col = "black"),
  motif |> filter(Chr == chr, motif_start >= from, motif_stop <= to) |>
    transmute(start = motif_start, end = motif_stop, fill = "red", col = "red")) |> distinct()
```

---

## 18.3 Plot every region

```r
for (i in seq_len(nrow(regions))) {
  r   <- regions[i]
  sig <- signal_tracks(r$Chr)
  hl  <- highlight_sites(r$Chr, r$Start, r$End)
  if (nrow(hl))
    sig <- list(HighlightTrack(sig, chromosome = r$Chr, start = hl$start, end = hl$end,
                               fill = hl$fill, col = hl$col, alpha = 0.6, lwd = 0.7,
                               inBackground = FALSE))

  pdf(sprintf("p_track_%s_%s_%d_%d.pdf", r$Gene, r$Chr, r$Start, r$End), width = 9, height = 5)
  plotTracks(c(list(axis_track), sig, list(peak_track, gene_track)),
             chromosome = r$Chr, from = r$Start, to = r$End, ylim = c(0, r$ymax),
             background.title = "white", col.border.title = NA, title.width = 1.6)
  dev.off()
}
```

<img src="figures/part18_gviz_highlight.png" alt="Fig. 18" width="640">

**[Fig. 18](figures/part18_gviz_highlight.pdf) | Signal at the *Prf1* locus with a highlighted differential region.** Pooled, RPM-normalized tracks for each condition on a shared y axis; grey box, differential region (DESeq2 padj < 0.05, |fold change| > 1.5); red line, TCF7 motif (FIMO); black bar, peak. Scale bar, 3 kb. *Example output.*

Tips:
- Use the same `ylim` for every condition (here the `ymax` column), otherwise samples cannot be compared; part 8.4 shows how to compute it from the bigWigs.
- To show all transcripts, drop `collapseTranscripts = "longest"` and add `stacking = "squish"`.
- Gene models from a GTF instead of a TxDb package: `GeneRegionTrack(makeTxDbFromGFF("gencode.vM25.annotation.gtf"), …)`.
- `sizes = c(0.6, 5, 0.3, 1)` in `plotTracks()` sets the relative height of each track (axis, highlight block, peaks, genes).

---

## References and resources

- Hahne, F. & Ivanek, R. Visualizing genomic data using Gviz and Bioconductor. *Methods Mol. Biol.* **1418**, 335–351 (2016). — <https://bioconductor.org/packages/Gviz>
- Gviz vignette — <https://bioconductor.org/packages/release/bioc/vignettes/Gviz/inst/doc/Gviz.html>
- Using Gviz (worked examples) — <https://davetang.org/muse/2013/10/03/using-gviz/>; <https://www.sthda.com/english/wiki/gviz-visualize-genomic-data>

---

[← ComplexHeatmap](17_complexheatmap.md) · [Index](../README.md)
