# Part 17 · Annotated heatmaps with ComplexHeatmap

[← PCA plots](16_pca_plot.md) · [Index](../README.md) · [Next: Gviz highlighted tracks →](18_gviz_tracks.md)

**ComplexHeatmap** draws publication-quality heatmaps with row and column annotations, split columns and labels for selected genes. This example plots row z-scores of differential peaks that carry a TCF7 motif, splits columns by condition and labels genes of interest. For a k-means clustered heatmap of RNA-seq DEGs see [part 9.3](09_rnaseq.md#93-k-means-clustered-heatmap-of-degs-complexheatmap).

---

## 17.1 Plot function

```r
library(data.table); library(dplyr); library(ComplexHeatmap); library(circlize)

plot_heatmap <- function(mat, group, group_cols, name, label_at = NULL, row_km = NULL,
                         cols = c("#000000", "#5DBAE9"), width = 4.2, height = 8) {
  z <- t(scale(t(mat)))                                        # row z-score

  top <- HeatmapAnnotation(
    group = anno_block(gp = gpar(fill = group_cols[levels(group)], col = NA),
                       labels = levels(group), labels_gp = gpar(col = "white", fontsize = 12)),
    show_annotation_name = FALSE)
  right <- if (length(label_at))
    rowAnnotation(gene = anno_mark(at = label_at, labels = rownames(mat)[label_at],
                                   labels_gp = gpar(fontsize = 12, fontface = "italic")))

  set.seed(100)                                                # reproducible row_km
  ht <- Heatmap(z, name = "Z-score", col = colorRamp2(c(-3, 3), cols),
                row_km = row_km, cluster_rows = TRUE, cluster_columns = FALSE,
                show_row_dend = FALSE, show_row_names = FALSE, show_column_names = FALSE,
                column_split = group, column_title = NULL, column_gap = unit(0, "mm"),
                top_annotation = top, right_annotation = right, use_raster = TRUE,
                heatmap_legend_param = list(direction = "horizontal", at = c(-3, 0, 3),
                                            legend_width = unit(3, "inch")))

  pdf(paste0("p_ComplexHeatmap_", name, ".pdf"), width = width, height = height)
  draw(ht, heatmap_legend_side = "bottom")
  dev.off()
}
```

---

## 17.2 Example: differential peaks with a TCF7 motif

```r
comps      <- c("cond1", "cond2", "cond3")
group_cols <- c(cond1 = "#b8b8b8", cond2 = "#2c7fb8", cond3 = "#a1dab4")

NormRC <- fread("NormRC_wAnno.txt")                            # peakName, Gene.Name, <samples>
fimo   <- fread("fimo_TCF7.Merge_wDER_padj0.05_FC1.5_noNA.txt")  # DERs with a TCF7 motif (part 7)
genes_to_label <- fread("genes_to_label.csv")$Gene             # one column named Gene

dt  <- NormRC |> filter(peakName %in% fimo$peakName)
mat <- dt |> select(contains(comps)) |> as.matrix()
rownames(mat) <- dt$Gene.Name

## sample columns look like cond2.g1.Rep1 → group cond2
group <- factor(sub("[.].*", "", colnames(mat)), levels = names(group_cols))

plot_heatmap(mat, group, group_cols, name = "DER_padj0.05_FC1.5_TCF7",
             label_at = which(dt$Gene.Name %in% genes_to_label))
```

<img src="figures/part17_complexheatmap.png" alt="Fig. 17" width="260">

**[Fig. 17](figures/part17_complexheatmap.pdf) | Differential peaks with a TCF7 motif.** Row z-scores of normalized counts at differential peaks (DESeq2 padj < 0.05, |fold change| > 1.5) carrying a TCF7 motif (FIMO), columns split by condition; selected nearby genes are labelled. *Example output.*

Tips:
- `anno_mark()` labels only the rows you pick and draws connector lines, so labels never overlap even on thousands of rows.
- `colorRamp2(c(-3, 3), …)` caps the colour scale at ±3; use `c(-2, 0, 2)` with three colours for a diverging red–blue scale.
- Add `row_km = 4` to split rows into k-means clusters; `row_order(draw(ht))` returns the rows in each cluster.
- Several heatmaps side by side (e.g. ATAC + RNA for the same genes): `draw(ht1 + ht2)`; rows follow the first heatmap.

---

## References and resources

- Gu, Z., Eils, R. & Schlesner, M. Complex heatmaps reveal patterns and correlations in multidimensional genomic data. *Bioinformatics* **32**, 2847–2849 (2016).
- Gu, Z. Complex heatmap visualization. *iMeta* **1**, e43 (2022).
- ComplexHeatmap reference book — <https://jokergoo.github.io/ComplexHeatmap-reference/book/>
- Vignette — <https://www.bioconductor.org/packages/release/bioc/vignettes/ComplexHeatmap/inst/doc/complex_heatmap.html>; frequently asked questions — <https://www.bioconductor.org/packages/release/bioc/vignettes/ComplexHeatmap/inst/doc/most_probably_asked_questions.html>

---

[← PCA plots](16_pca_plot.md) · [Index](../README.md) · [Next: Gviz highlighted tracks →](18_gviz_tracks.md)
