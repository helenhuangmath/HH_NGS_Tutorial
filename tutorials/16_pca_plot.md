# Part 16 · PCA plots of samples (ggplot2)

[← Spatial transcriptomics](15_spatial.md) · [Index](../README.md) · [Next: ComplexHeatmap →](17_complexheatmap.md)

A PCA plot with colour = condition and shape = guide (or replicate group), drawn for PC1/PC2 and PC2/PC3. Use it on **rlog- or VST-transformed** counts from DESeq2 ([part 5](05_differential_analysis.md#52-deseq2-on-peak-counts)) so a few very strong peaks or genes do not dominate the PCs.

---

## 16.1 Input and sample sheet

```r
library(data.table); library(dplyr); library(ggplot2)
source("Code/plot_funcs.R")                  # theme_nature(), save_fig() from part 5.3

comps <- c("cond1", "cond2", "cond3")

## rlog counts: one row per peak/gene, one column per sample (+ annotation columns)
rlog_counts <- fread("NormRC_rlog_byDESeq_wAnno.txt")
mat <- rlog_counts |> select(contains(comps)) |> as.matrix()
rownames(mat) <- rlog_counts$GeneInfo

## one row per sample, in the column order of mat
samples <- data.frame(
  sample = colnames(mat),
  group  = rep(c("cond1", "cond2", "cond3"), times = c(2, 4, 4)),
  guide  = c("ctrl", "ctrl", "g1", "g1", "g2", "g2", "g1", "g1", "g2", "g2"))

group_cols   <- c(cond1 = "#b8b8b8", cond2 = "#2c7fb8", cond3 = "#a1dab4")
guide_shapes <- c(ctrl = 22, g1 = 21, g2 = 24)        # filled square, circle, triangle
```

Without DESeq2 output, `rlog(dds)` or `vst(dds)` gives the same matrix: `mat <- assay(vst(dds))`.

---

## 16.2 Plot function

```r
plot_pca <- function(mat, samples, name, pcs = c(1, 2), ntop = 5000) {
  top <- head(order(apply(mat, 1, var), decreasing = TRUE), ntop)   # most variable rows
  pca <- prcomp(t(mat[top, ]))
  ve  <- round(100 * summary(pca)$importance[2, ], 1)
  x <- paste0("PC", pcs[1]); y <- paste0("PC", pcs[2])

  p <- ggplot(cbind(samples, pca$x), aes(.data[[x]], .data[[y]], fill = group, shape = guide)) +
    geom_hline(yintercept = 0, colour = "grey85", linewidth = 0.2) +
    geom_vline(xintercept = 0, colour = "grey85", linewidth = 0.2) +
    geom_point(size = 3, stroke = 0.3) +
    scale_fill_manual(values = group_cols) +
    scale_shape_manual(values = guide_shapes) +
    guides(fill = guide_legend(override.aes = list(shape = 21))) +   # show fill colours in the legend
    labs(x = sprintf("%s (%s%%)", x, ve[pcs[1]]),
         y = sprintf("%s (%s%%)", y, ve[pcs[2]]), fill = NULL, shape = NULL) +
    theme_nature()
  save_fig(p, sprintf("plot_PCA_%s_%s%s.pdf", name, x, y), w = 70, h = 50)
  p
}

plot_pca(mat, samples, "ATAC_rlog", pcs = c(1, 2))
plot_pca(mat, samples, "ATAC_rlog", pcs = c(2, 3))
```

<img src="figures/part16_pca.png" alt="Fig. 16" width="420">

**[Fig. 16](figures/part16_pca.pdf) | PCA of ATAC-seq samples.** PC1 against PC2 of rlog-transformed counts; colour, condition; shape, guide. Both guides of each condition cluster together and away from the cond1 control. *Example output.*

Tips:
- Replicates should cluster together. If samples split by sequencing batch instead, add the batch to the DESeq2 design (`~ batch + Condition`) and remove it with `limma::removeBatchEffect()` before PCA (for plotting only).
- To find which peaks or genes drive a PC, look at its loadings (`pca$rotation`), see [PCA by hand](05_differential_analysis.md#pca-by-hand-loadings).
- DESeq2 shortcut: `plotPCA(vst(dds), intgroup = c("group", "guide"), returnData = TRUE)` returns the coordinates for the same ggplot.

---

## References and resources

- Ringnér, M. What is principal component analysis? *Nat. Biotechnol.* **26**, 303–304 (2008).
- Love, M. I., Huber, W. & Anders, S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biol.* **15**, 550 (2014). — <https://bioconductor.org/packages/DESeq2>
- ggplot2 — <https://ggplot2.tidyverse.org>

---

[← Spatial transcriptomics](15_spatial.md) · [Index](../README.md) · [Next: ComplexHeatmap →](17_complexheatmap.md)
