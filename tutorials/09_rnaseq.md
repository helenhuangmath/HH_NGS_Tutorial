# Part 9 · Bulk RNA-seq: STAR → counts → clustered heatmap → GSEA

[← Visualization](08_visualization.md) · [Index](../README.md) · [Next: TSA-seq →](10_tsa_seq.md)

Differential expression uses the same DESeq2 code as part 5, with genes in place of peaks.

---

## 9.1 Alignment with STAR

```bash
#!/bin/bash
## align_rna.sh <SAMPLE>
set -euo pipefail
module load STAR/2.7.1a samtools-1.9

sample=$1
datapath=/path/to/Project/Fastq
ref=/path/to/indexes/mm10/STAR_Gencode
BLACKLIST=/path/to/indexes/mm10/mm10-blacklist.clean.bed
cd /path/to/Project/RNA

STAR --runThreadN 8 --genomeDir ${ref} \
     --readFilesIn ${datapath}/${sample}_R1.fastq.gz ${datapath}/${sample}_R2.fastq.gz \
     --readFilesCommand zcat \
     --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 42000000000 \
     --quantMode GeneCounts \
     --outFileNamePrefix ${sample}_

samtools index ${sample}_Aligned.sortedByCoord.out.bam

bamCoverage -b ${sample}_Aligned.sortedByCoord.out.bam -o ${sample}.cpm.bw \
  -bl ${BLACKLIST} -p 8 --normalizeUsing CPM --binSize 5 --smoothLength 10 \
  --ignoreForNormalization chrM --skipNAs --minMappingQuality 5
# stranded tracks: add --filterRNAstrand forward / reverse
```

For single-end reads, give one FASTQ to `--readFilesIn`. Build the index once with `STAR --runMode genomeGenerate --genomeDir … --genomeFastaFiles genome.fa --sjdbGTFfile gencode.gtf --sjdbOverhang 100`.

**Key outputs**

| File | Content |
|---|---|
| `*_Log.final.out` | uniquely mapped %, multi-mapped %, splice junctions |
| `*_ReadsPerGene.out.tab` | counts: col 2 unstranded, col 3 forward, col 4 reverse |
| `*_Aligned.sortedByCoord.out.bam` | coordinate-sorted BAM |

Check strandedness with RSeQC (part 2.4) and pick the matching column.

---

## 9.2 Gene counts with featureCounts

```bash
featureCounts -T 8 -p --countReadPairs -s 2 \
  -a gencode.v42.annotation.gtf -t exon -g gene_id \
  -o RawRC_RNA_featureCounts.txt *_Aligned.sortedByCoord.out.bam
```

### Clean the count table (add gene names, drop all-zero genes)

```r
## clean_featureCounts.R  — usage: Rscript … RawRC_RNA_featureCounts.txt
suppressPackageStartupMessages({ library(data.table); library(dplyr) })

fname <- commandArgs(trailingOnly = TRUE)[1]

## BED made from the GTF: chr start end gene_id gene_name strand
gencode <- fread("/path/to/index/hg38/gencode.v42.annotation.gtf.bed",
                 col.names = c("Chr", "Start", "End", "Ensembl", "Gene.Name", "Strand"))

dt <- fread(fname, skip = 1)                                    # line 1 is the command
colnames(dt) <- sub("_Aligned.sortedByCoord.out.bam$", "", basename(colnames(dt)))

dt <- dt |>
  select(-Chr, -Start, -End, -Strand, -Length) |>
  inner_join(gencode, by = c("Geneid" = "Ensembl")) |>
  select(Chr, Start, End, Geneid, Gene.Name, Strand, everything())
fwrite(dt, sub(".txt$", ".clean.txt", fname), sep = "\t")

expressed <- apply(dt[, -(1:6)], 1, max) > 0
fwrite(dt[expressed, ], sub(".txt$", ".clean_noZero.txt", fname), sep = "\t")
message(sum(expressed), " / ", nrow(dt), " genes with ≥1 count")
```

Make the GTF → BED lookup once:

```bash
awk 'BEGIN{OFS="\t"} $3=="gene" {
       match($0, /gene_id "[^"]+"/);   id=substr($0, RSTART+9,  RLENGTH-10);
       match($0, /gene_name "[^"]+"/); nm=substr($0, RSTART+11, RLENGTH-12);
       print $1, $4-1, $5, id, nm, $7 }' gencode.v42.annotation.gtf > gencode.v42.annotation.gtf.bed
```

---

## 9.3 k-means clustered heatmap of DEGs (ComplexHeatmap)

Goal: z-scored expression of all DEGs (any pairwise contrast), split into *k* clusters, with annotation columns showing the direction of each contrast.

```r
library(data.table); library(dplyr); library(ComplexHeatmap); library(circlize); library(RColorBrewer)

conds <- c(TN = "#8C8C8C", TMEM = "#0072B2", TEX = "#D55E00")
k     <- 7
set.seed(1)

## log-normalized counts of DEGs (rows = genes, columns = samples ordered by condition)
expr <- fread("NormRC_RNA_DEG_any.txt") |> as.data.frame()
mat  <- log2(as.matrix(expr[, -1]) + 1); rownames(mat) <- expr$Gene.Name
z    <- t(scale(t(mat)))                                             # row z-score

## choose k: elbow / silhouette (factoextra::fviz_nbclust(z, kmeans, method = "wss"))
km  <- kmeans(z, centers = k, nstart = 50, iter.max = 100)
## order clusters by the condition with the highest mean z (makes the figure readable)
peak_cond <- apply(km$centers, 1, \(v) which.max(tapply(v, sub("_.*", "", colnames(z)), mean)))
cl <- factor(paste0("C", match(km$cluster, order(peak_cond))), levels = paste0("C", 1:k))
fwrite(data.frame(Gene.Name = rownames(z), Cluster = cl), "DEG_clusters_km7.txt", sep = "\t")

## DE direction per contrast (−1 / 0 / 1) from DESeq2 tables
de <- fread("DEG_direction_3comps.txt")              # Gene.Name TMEMvsTN TEXvsTN TEXvsTMEM
de <- de[match(rownames(z), de$Gene.Name)]
dir_col <- c("-1" = "#0072B2", "0" = "white", "1" = "#D55E00")

col_fun  <- colorRamp2(c(-2, 0, 2), c("#2166AC", "#F7F7F7", "#B2182B"))
cond_vec <- factor(sub("_.*", "", colnames(z)), levels = names(conds))

top  <- HeatmapAnnotation(Condition = cond_vec, col = list(Condition = conds),
                          simple_anno_size = unit(2, "mm"), annotation_name_gp = gpar(fontsize = 6))
left <- rowAnnotation(Cluster = cl, col = list(Cluster = setNames(brewer.pal(k + 1, "BuPu")[-1], levels(cl))),
                      simple_anno_size = unit(2, "mm"), show_annotation_name = FALSE)
right <- rowAnnotation(`TMEM vs TN` = as.character(de$TMEMvsTN),
                       `TEX vs TN`  = as.character(de$TEXvsTN),
                       `TEX vs TMEM`= as.character(de$TEXvsTMEM),
                       col = list(`TMEM vs TN` = dir_col, `TEX vs TN` = dir_col, `TEX vs TMEM` = dir_col),
                       simple_anno_size = unit(2.5, "mm"), annotation_name_gp = gpar(fontsize = 6),
                       show_legend = c(TRUE, FALSE, FALSE),
                       annotation_legend_param = list(title = "DE call", at = c("1", "-1"),
                                                      labels = c("Up", "Down")))

ht <- Heatmap(z, name = "Row z-score", col = col_fun,
              cluster_rows = FALSE, cluster_columns = FALSE,
              row_split = cl, row_gap = unit(0.6, "mm"), row_title_gp = gpar(fontsize = 6),
              show_row_names = FALSE, show_column_names = FALSE,
              top_annotation = top, left_annotation = left, right_annotation = right,
              use_raster = TRUE, raster_quality = 3,
              column_title = sprintf("DEGs (FDR < 0.05, FC > 1.5; n = %s)", format(nrow(z), big.mark = ",")),
              column_title_gp = gpar(fontsize = 7),
              heatmap_legend_param = list(title_gp = gpar(fontsize = 6), labels_gp = gpar(fontsize = 6)))

pdf(paste0("plot_heat_DEG_km", k, ".pdf"), width = 89 / 25.4, height = 4.2)
draw(ht, merge_legend = TRUE, heatmap_legend_side = "right")
dev.off()
```

![Fig. 9.1](figures/part09_kmeans_heatmap.png)

**[Fig. 9.1](figures/part09_kmeans_heatmap.pdf) | k-means clustered heatmap of DEGs.** Row-scaled expression (z-score) of DEGs from any pairwise comparison (DESeq2 FDR < 0.05, fold change > 1.5) across four replicates per condition, split into seven k-means clusters (C1–C7, left bar). Right, direction of change in each pairwise contrast (orange, up; blue, down; white, not significant). *Simulated data.*

> `use_raster = TRUE` keeps PDFs with thousands of rows small. Journals accept rasterized heatmap cells inside a vector figure.

For a heatmap split by condition with selected genes labelled, see [part 17](17_complexheatmap.md).

---

## 9.4 Gene set enrichment analysis (GSEA)

GSEA tests whether the genes of a set (e.g. exhausted T-cell subset signatures) sit at the top or bottom of a list of **all** genes ranked by the DE statistic. Unlike over-representation analysis (part 6), no DEG cut-off is used, so do not pre-filter the ranked list.

### Ranked list from DESeq2

```r
library(data.table); library(dplyr)

res <- fread("DESeq2_KO_vs_Ctrl.txt")             # one row per gene: gene_name, stat, ...
rnk <- res |>
  filter(!is.na(stat), gene_name != "") |>
  group_by(gene_name) |> slice_max(abs(stat), n = 1, with_ties = FALSE) |> ungroup() |>
  arrange(desc(stat)) |> select(gene_name, stat)
fwrite(rnk, "KO_vs_Ctrl.rnk", sep = "\t", col.names = FALSE)
```

Rank by the Wald `stat` (or `sign(log2FoldChange) * -log10(pvalue)`). Shrunken log₂ fold changes alone ignore the variance and give a poorer ranking.

### Option A: GSEA preranked (Broad GSEA-CLI)

```bash
# gene sets: GMT file, one set per line (name, description, genes …), e.g. Tex_subsets.gmt
gsea-cli.sh GSEAPreranked -rnk KO_vs_Ctrl.rnk -gmx Tex_subsets.gmt \
  -collapse No_Collapse -scoring_scheme weighted -nperm 1000 -rnd_seed 149 \
  -set_min 15 -set_max 500 -make_sets true -create_svgs false -zip_report false \
  -rpt_label KO_vs_Ctrl -out gsea_out
# → gsea_out/KO_vs_Ctrl.GseaPreranked.<id>/
#     gsea_report_for_na_pos_<id>.tsv, gsea_report_for_na_neg_<id>.tsv   (NES, FDR q-val per set)
#     <SET_NAME>.tsv   (one row per gene of the set: RANK IN GENE LIST, RANK METRIC SCORE, RUNNING ES)
```

GSEA ≤ 4.0 writes the per-set tables with an `.xls` extension; they are still tab-delimited text (`read.delim` reads both).

### Option B: fgsea in R

```r
library(fgsea); library(tibble)

ranks    <- deframe(rnk)                          # named vector, decreasing
pathways <- gmtPathways("Tex_subsets.gmt")        # or split(msigdbr(...)$gene_symbol, ...$gs_name)
fg <- fgseaMultilevel(pathways, ranks, minSize = 15, maxSize = 500, eps = 0) |> arrange(padj)
fwrite(fg, "fgsea_KO_vs_Ctrl.txt", sep = "\t", sep2 = c("", ",", ""))
```

### Running-enrichment-score plot

Both options feed the same plot. GSEA's per-set table lists only the genes in the set, so `es_from_gsea()` rebuilds the full curve. Between two hits the score falls by 1/(*N* − *N*ₕ) per gene. Joining only the hit positions would miss these drops and overstate the curve.

```r
library(ggplot2)
# theme_nature() and save_fig() are defined in part 5

## full running-ES curve from a GSEA <SET_NAME>.tsv; n_genes = length of the ranked list
es_from_gsea <- function(tab, n_genes) {
  r  <- tab$RANK.IN.GENE.LIST                     # 0-based positions of the hits
  es <- tab$RUNNING.ES
  d  <- 1 / (n_genes - length(r))                 # decrement per non-hit gene
  jump <- es - (c(0, head(es, -1)) - (r - c(-1, head(r, -1)) - 1) * d)
  step <- rep(-d, n_genes); step[r + 1] <- jump
  hit  <- seq_len(n_genes) %in% (r + 1)
  metric <- rep(NA_real_, n_genes); metric[r + 1] <- tab$RANK.METRIC.SCORE
  data.frame(rank = seq_len(n_genes) - 1, es = cumsum(step), hit = hit, metric = metric)
}

## same curve computed directly from a ranked vector (weighted GSEA, p = 1)
es_from_ranks <- function(ranks, genes) {
  ranks <- sort(ranks, decreasing = TRUE)
  hit <- names(ranks) %in% genes
  inc <- abs(ranks) * hit
  es  <- cumsum(inc / sum(inc) - (!hit) / sum(!hit))
  data.frame(rank = seq_along(ranks) - 1, es = es, hit = hit, metric = unname(ranks))
}

plot_gsea <- function(curve, title, colour = "#D55E00", label = NULL,
                      groups = c("KO", "Ctrl")) {
  lo <- min(0, curve$es); hi <- max(0, curve$es); span <- hi - lo
  tick <- c(lo - 0.22 * span, lo - 0.08 * span)   # hit ticks in a strip below the curve
  hits <- curve[curve$hit, ]
  ggplot(curve, aes(rank, es)) +
    geom_hline(yintercept = 0, linewidth = 0.25, colour = "grey50") +
    geom_line(colour = colour, linewidth = 0.5) +
    geom_segment(data = hits, aes(x = rank, xend = rank, y = tick[1], yend = tick[2],
                                  colour = metric > 0), linewidth = 0.15) +
    scale_colour_manual(values = c(`TRUE` = "#D55E00", `FALSE` = "#0072B2"), guide = "none") +
    annotate("text", x = 0, y = tick[1], label = groups[1], hjust = 0, vjust = 1.6, size = 1.8) +
    annotate("text", x = max(curve$rank), y = tick[1], label = groups[2], hjust = 1, vjust = 1.6,
             size = 1.8) +
    { if (!is.null(label)) annotate("text", x = Inf, y = Inf, label = label, hjust = 1.05,
                                    vjust = 1.4, size = 2) } +
    scale_y_continuous(expand = expansion(mult = c(0.12, 0.05))) +
    labs(title = title, x = "Rank in ordered gene list", y = "Running enrichment score") +
    theme_nature() + theme(aspect.ratio = 0.8)
}

## Option A input: one plot per gene set, NES and FDR from the GSEA summary tables
run_dir <- list.files("gsea_out", pattern = "^KO_vs_Ctrl\\.GseaPreranked", full.names = TRUE)[1]
summ <- list.files(run_dir, pattern = "^gsea_report_for_na_(pos|neg).*\\.(tsv|xls)$",
                   full.names = TRUE) |> lapply(read.delim) |> bind_rows()
n_genes <- nrow(rnk)
set_cols <- c(TEXPROG1 = "#EBA352", TEXPROG2 = "#EBA352", TEXINT = "#DC8D9C", TEXTERM = "#D95959")

for (s in names(set_cols)) {
  f   <- list.files(run_dir, pattern = paste0("^", s, "\\.(tsv|xls)$"), full.names = TRUE)
  tab <- read.delim(f)
  st  <- summ[summ$NAME == s, ]
  p <- plot_gsea(es_from_gsea(tab, n_genes), s, set_cols[[s]],
                 label = sprintf("NES = %.2f\nFDR = %.2g", st$NES, st$FDR.q.val))
  save_fig(p, paste0("p_GSEA_KO_vs_Ctrl_", s, ".pdf"), w = 60, h = 55)
}

## Option B input
# p <- plot_gsea(es_from_ranks(ranks, pathways[["TEXTERM"]]), "TEXTERM", "#D95959",
#                label = with(fg[fg$pathway == "TEXTERM", ], sprintf("NES = %.2f\nFDR = %.2g", NES, padj)))
```

> Ticks are coloured by the sign of the ranking metric: orange = higher in the first group (left end of the list), blue = higher in the second. A set is enriched in the first group when the curve peaks early and most of its ticks fall in the orange part of the strip. Report NES and FDR, not the nominal *P* value. With many gene sets, FDR < 0.25 is the GSEA convention for exploratory work; use a stricter cut-off for claims.

---

## 9.5 Transposable-element expression (SQuIRE)

```bash
conda create -n squire -c bioconda -c conda-forge -c r squire   # see the SQuIRE README for pinned versions
conda activate squire

# Preparation (once per genome)
bsub -M 140000 -o log_fetch.o -e log_fetch.e \
  squire Fetch -b mm10 -o ~/tools/squire/mm10 -f -c -r -g -p 12 -v -x
squire Clean -r ~/tools/squire/mm10/mm10_rmsk.txt -b mm10 -o ~/tools/squire/mm10/cleaned

# Quantification (per sample)
bsub -n 16 -M 65536 -o log_map.o -e log_map.e \
  squire Map -1 ${s}_R1.fastq.gz -2 ${s}_R2.fastq.gz -o squire_map/${s} \
             -f ~/tools/squire/mm10 -r 150 -n ${s} -b mm10 -p 16 -v
bsub -n 16 -M 65536 -o log_cnt.o -e log_cnt.e \
  squire Count -m squire_map/${s} -c ~/tools/squire/mm10/cleaned -o squire_count/${s} \
               -t squire_count/tmp -f ~/tools/squire/mm10 -r 150 -n ${s} -b mm10 -p 16 -s 2 -v
```

`-r` = read length; `-s` = strandedness (0 unstranded, 1 forward, 2 reverse). Output: locus- and subfamily-level TE counts for `squire Call` (DESeq2).

---

## References and resources

- Dobin, A. *et al.* STAR: ultrafast universal RNA-seq aligner. *Bioinformatics* **29**, 15–21 (2013). — <https://github.com/alexdobin/STAR>
- Liao, Y., Smyth, G. K. & Shi, W. featureCounts. *Bioinformatics* **30**, 923–930 (2014).
- Love, M. I., Huber, W. & Anders, S. DESeq2. *Genome Biol.* **15**, 550 (2014).
- Subramanian, A. *et al.* Gene set enrichment analysis: a knowledge-based approach for interpreting genome-wide expression profiles. *Proc. Natl Acad. Sci. USA* **102**, 15545–15550 (2005). — <https://www.gsea-msigdb.org/gsea/>
- Korotkevich, G. *et al.* Fast gene set enrichment analysis. *bioRxiv* 060012 (2021). (fgsea) — <https://bioconductor.org/packages/fgsea>
- Gu, Z., Eils, R. & Schlesner, M. Complex heatmaps reveal patterns and correlations in multidimensional genomic data. *Bioinformatics* **32**, 2847–2849 (2016). — book: <https://jokergoo.github.io/ComplexHeatmap-reference/book/>
- Yang, W. R. *et al.* SQuIRE reveals locus-specific regulation of interspersed repeat expression. *Nucleic Acids Res.* **47**, e27 (2019). — <https://github.com/wyang17/SQuIRE>
- Frankish, A. *et al.* GENCODE reference annotation for the human and mouse genomes. *Nucleic Acids Res.* **47**, D766–D773 (2019).
- RNA-seq workflow (Bioconductor) — <https://bioconductor.org/packages/rnaseqGene>

---

[← Visualization](08_visualization.md) · [Index](../README.md) · [Next: TSA-seq →](10_tsa_seq.md)
