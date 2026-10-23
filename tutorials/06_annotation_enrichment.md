# 6 · Peak annotation and functional enrichment

[← Differential analysis](05_differential_analysis.md) · [Index](../README.md) · [Next: Motif analysis →](07_motif_analysis.md)

---

## 6.1 HOMER `annotatePeaks.pl`

Input: a BED with a **unique ID in column 4** and a **strand in column 6**:

```text
chr1  1000512  1002240  peak1  .  +
```

```bash
export PATH=$PATH:/path/to/tools/homer/bin      # or: module load homer-v4.6
species=hg19                                       # or a custom genome, e.g. mm10_GencodeV24 (chapter 0)
p=peak_Union_H3K27ac

annotatePeaks.pl ${p}.bed6.bed ${species} 1> ${p}.homer_anno.txt 2> log_homer_${p}.txt
Rscript Code/clean_homer.R ${p}.homer_anno.txt          # → ${p}.homer_anno.simp.txt
```

`clean_homer.R` simplifies the HOMER table and defines the promoter as **−2 kb to +0.5 kb from the TSS**:

```r
## clean_homer.R  — usage: Rscript clean_homer.R <homer_anno.txt>
suppressPackageStartupMessages({ library(data.table); library(dplyr); library(stringr) })

annolevels <- c("Promoter", "TTS", "5' UTR", "3' UTR", "Exon", "Intron", "Intergenic", "Non-Coding")
f  <- commandArgs(trailingOnly = TRUE)[1]
dt <- fread(f); colnames(dt)[1] <- "peakN"          # column 1 header is "PeakID (cmd=...)"

simp <- data.frame(dt) |>                            # data.frame() → "Gene Name" becomes Gene.Name
  mutate(width = End - Start,
         Anno  = str_to_title(str_replace(Annotation, "\\s\\(.*", "")),
         Anno  = recode(Anno, "Promoter-Tss" = "Promoter", "5' Utr" = "5' UTR",
                        "3' Utr" = "3' UTR", "Tts" = "TTS"),
         Anno  = ifelse(Distance.to.TSS > -2000 & Distance.to.TSS < 500, "Promoter", Anno),
         Anno  = factor(Anno, levels = annolevels)) |>
  arrange(Chr, Start) |>
  select(Chr, Start, End, peakN, width, Gene.Name, Strand, Distance.to.TSS, Anno) |>
  distinct()

print(table(simp$Anno))
fwrite(simp, sub("homer_anno.txt$", "homer_anno.simp.txt", f), sep = "\t")
```

Example output (illustrative):

```text
  Promoter        TTS     5' UTR     3' UTR       Exon     Intron Intergenic Non-Coding
     11542        812        402        790       1630      13551      12233        276
```

Merge the annotation onto the DESeq2 table:

```r
anno <- fread("peak_Union_H3K27ac.homer_anno.simp.txt")
res  <- fread("deseq2_H3K27ac.txt") |> left_join(anno |> select(peakN, Gene.Name, Distance.to.TSS, Anno))
```

### Annotation distribution plot

```r
anno_cols <- c(Promoter = "#0072B2", `5' UTR` = "#56B4E9", `3' UTR` = "#009E73", Exon = "#F0E442",
               Intron = "#E69F00", Intergenic = "#D55E00", TTS = "#CC79A7", `Non-Coding` = "grey70")
df <- bind_rows(All = res, Up = filter(res, Sig.Dir == 1), Down = filter(res, Sig.Dir == -1), .id = "Set") |>
  count(Set, Anno) |> group_by(Set) |> mutate(pct = 100 * n / sum(n))
ggplot(df, aes(pct, Set, fill = Anno)) + geom_col(width = 0.65, colour = "white", linewidth = 0.2) +
  scale_fill_manual(values = anno_cols) + labs(x = "Peaks (%)", y = NULL, fill = NULL) + theme_nature()
```

(Stacked bars make it easier than pie charts to compare peak sets side by side.)

---

## 6.2 GREAT (regulatory-domain gene assignment)

GREAT assigns each region to the nearby genes whose regulatory domains it falls in. This works better than nearest-TSS for distal enhancers.

```r
library(data.table); library(GenomicRanges); library(rGREAT)

peaks <- fread("peak_Union_H3K27ac.bed.txt", select = 1:3, col.names = c("chr", "start", "end"))
gr    <- makeGRangesFromDataFrame(peaks)

## online GREAT (submitGreatJob). adv_oneDistance is in **kb**: 1000 = 1 Mb (GREAT default)
job <- submitGreatJob(gr, species = "hg19", rule = "oneClosest", adv_oneDistance = 1000)

assoc <- getRegionGeneAssociations(job)    # older rGREAT: plotRegionGeneAssociationGraphs(job)
fwrite(as.data.frame(assoc), "peak_Union_H3K27ac_wAnno_rGreat.txt", sep = "\t")

tbl <- getEnrichmentTables(job)            # GO MF / BP / CC
bp  <- tbl[["GO Biological Process"]]
head(bp[order(bp$Binom_Adjp_BH), c("name", "Binom_Fold_Enrichment", "Binom_Adjp_BH")])
```

> **Units:** `adv_oneDistance` is in kilobases, so `5000` means 5 Mb, not 5 kb.
>
> rGREAT ≥ 2.0 can also run **locally** without the web server: `great(gr, "GO:BP", "txdb:hg19")`. Use this for genomes or gene sets the server does not support.

GO bar plot per direction (the logic of `plotgreat()`):

```r
top <- head(bp[order(bp$Binom_Adjp_BH), ], 10)
ggplot(top, aes(-log10(Binom_Adjp_BH), reorder(name, -Binom_Adjp_BH))) +
  geom_col(fill = "#0072B2", width = 0.7) +
  geom_vline(xintercept = -log10(0.05), linetype = 2, linewidth = 0.2) +
  labs(x = expression(-log[10]~adjusted~italic(P)), y = NULL) + theme_nature()
```

### Reduce redundant GO terms (REVIGO)

GREAT and clusterProfiler often return dozens of near-identical GO terms. REVIGO groups them by semantic similarity. Input is two columns, `GO-ID  P-value`:

```text
GO:0042110  1.2e-12
GO:0050863  3.4e-11
```

**Offline in R (rrvgo)**

```r
library(rrvgo)
go  <- fread("GREAT_GO_BP_Up.txt")                 # columns ID, Binom_Adjp_BH
sim <- calculateSimMatrix(go$ID, orgdb = "org.Mm.eg.db", ont = "BP", method = "Rel")
red <- reduceSimMatrix(sim, scores = setNames(-log10(go$Binom_Adjp_BH), go$ID),
                       threshold = 0.7, orgdb = "org.Mm.eg.db")
fwrite(red, "GREAT_GO_BP_Up_reduced.csv")
treemapPlot(red)
```

**REVIGO web service (Python)**

REVIGO offers a job API for scripted submission (see *revigo.irb.hr → FAQ → programmatic access*; check the current endpoint names there before use):

```python
#!/usr/bin/env python3
## usage: python revigo_auto.py GO_pvalues.txt   →  GO_pvalues_revigo.csv
import sys, time, requests

infile = sys.argv[1]
payload = {"cutoff": "0.7", "valueType": "pvalue", "speciesTaxon": "10090",   # 10090 mouse, 9606 human
           "measure": "SIMREL", "goList": open(infile).read()}
job = requests.post("http://revigo.irb.hr/StartJob", data=payload).json()["jobid"]

while requests.post("http://revigo.irb.hr/QueryJob",
                    data={"jobid": job, "type": "jstatus"}).json()["running"]:
    time.sleep(1)

# namespace 1 = biological process, 2 = cellular component, 3 = molecular function
table = requests.post("http://revigo.irb.hr/QueryJob",
                      data={"jobid": job, "type": "table", "namespace": "1"}).text
open(infile.replace(".txt", "_revigo.csv"), "w").write(table)
```

---

## 6.3 Gene ID conversion (biomaRt)

```r
library(biomaRt)
mart_mm <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "mmusculus_gene_ensembl",
                   host = "https://oct2022.archive.ensembl.org")     # pin an archive for reproducibility
attrs <- c("ensembl_gene_id", "external_gene_name", "chromosome_name",
           "start_position", "end_position", "strand", "band", "description")
anno <- getBM(attributes = attrs, filters = "external_gene_name",
              values = c("Tox", "Tcf7"), mart = mart_mm)
# Ensembl IDs → symbols: filters = "ensembl_gene_id", values = sub("\\..*", "", ids)  (drop version suffix)
```

---

## 6.4 Over-representation analysis (clusterProfiler)

```r
library(clusterProfiler); library(msigdbr); library(dplyr); library(ggplot2); library(stringr)

run_enricher <- function(genes, term2gene, oname, universe = NULL) {
  er <- enricher(genes, TERM2GENE = term2gene, universe = universe,
                 pAdjustMethod = "BH", pvalueCutoff = 0.1, qvalueCutoff = 1)
  if (is.null(er) || nrow(as.data.frame(er)) == 0) { message("no enriched terms"); return(invisible(NULL)) }
  saveRDS(er, paste0("Enricher_", oname, ".rds"))
  fwrite(as.data.frame(er), paste0("Enricher_", oname, ".csv"))
  p <- dotplot(er, showCategory = 15) +
       scale_y_discrete(labels = \(x) str_wrap(x, 45)) + theme_nature()
  ggsave(paste0("p_Dot_Enricher_", oname, ".pdf"), p, width = 120, height = 90, units = "mm")
  er
}

## reference gene sets: MSigDB C5 (GO) for mouse
## msigdbr ≥ 10: use collection = "C5"; older versions: category = "C5"
C5 <- msigdbr(species = "Mus musculus", category = "C5") |> select(gs_name, gene_symbol)

genes <- c("Tox", "Klrg1", "Entpd1", "Pdcd1", "Clnk", "Gzma", "Zeb2", "Cx3cr1",
           "Fasl", "Prf1", "Ccr5", "Havcr2", "Ly6c2", "Ifng", "Gzmk", "Nr4a1")
er <- run_enricher(genes, C5, "MsigDB_C5_TestGeneSet")
```

> `enricher()` returns `NULL` or an empty result when nothing passes the cut-offs, which the function checks for. Supply `universe =` (all genes with a peak or detected in RNA-seq) to avoid inflated *P* values.

### Enrichr web API (`enrichR` package)

```r
library(enrichR)
setEnrichrSite("Enrichr")
dbs <- c("GO_Biological_Process_2023", "KEGG_2019_Mouse")
enriched <- enrichr(c("Runx1", "Gfi1", "Gfi1b", "Spi1", "Gata1", "Kdr"), dbs)
head(enriched[["GO_Biological_Process_2023"]][, c("Term", "Overlap", "Adjusted.P.value")])
plotEnrich(enriched[[1]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
```

Use `listEnrichrDbs()` to see the current library names; they change with each Enrichr release.

---

## 6.5 Overlap of gene lists (Venn / Euler)

```r
library(VennDiagram); library(eulerr)

sets <- list(Int = list1$Gene.Name, Prog1 = list2$Gene.Name,
             Prog2 = list3$Gene.Name, Term = list4$Gene.Name)
cols <- c("#E69F00", "#56B4E9", "#009E73", "#CC79A7")

## table of which genes fall in which intersection
ol <- calculate.overlap(sets)
gene_ol <- rbindlist(lapply(ol, \(g) data.table(Gene.Name = g)), idcol = "OL_Group")
fwrite(gene_ol, "OL_Venn_Groups.txt", sep = "\t")

## area-proportional Euler diagram (ellipses handle 4+ sets better than circles)
fit <- euler(sets, shape = "ellipse")
pdf("p_Venn_euler.pdf", width = 3, height = 2.2)
plot(fit, fills = list(fill = cols, alpha = 0.5), edges = FALSE,
     quantities = list(cex = 0.6), labels = list(cex = 0.7))
dev.off()
```

Other options: `VennDiagram::venn.diagram()` (classic), `ggVennDiagram()` (fill by count) and `ggvenn()` (quick). For more than four sets, use an UpSet plot (`ComplexHeatmap::UpSet`).

---

## Example output

![Fig. 5](figures/fig5_annotation_boxplot.png)

**Fig. 5 | Genomic annotation of H3K27ac regions.** **a**, HOMER annotation of all union peaks and of regions gained (Up) or lost (Down) in TEX; promoter = −2 kb to +0.5 kb from the TSS. **b**, Distance from each peak to the nearest TSS; the shaded band marks the promoter window. **c**, see chapter 5. *Simulated data.*

![Fig. 6c](figures/fig6_motif_enrichment.png)

**Fig. 6 | (panel c)** Over-representation analysis (clusterProfiler `enricher`, MSigDB C5) of genes near TEX-gained regions. Dot size = genes in the term; colour = −log₁₀ BH-adjusted *P*; x-axis = gene ratio. Panels a and b are described in chapter 7. *Simulated data.*

---

## References and resources

- Heinz, S. *et al.* *Mol. Cell* **38**, 576–589 (2010). (HOMER `annotatePeaks.pl`) — <http://homer.ucsd.edu/homer/ngs/annotation.html>
- McLean, C. Y. *et al.* GREAT improves functional interpretation of cis-regulatory regions. *Nat. Biotechnol.* **28**, 495–501 (2010). — <http://great.stanford.edu>
- Gu, Z. & Hübschmann, D. rGREAT: an R/Bioconductor package for functional enrichment on genomic regions. *Bioinformatics* **39**, btac745 (2023).
- Supek, F., Bošnjak, M., Škunca, N. & Šmuc, T. REVIGO summarizes and visualizes long lists of Gene Ontology terms. *PLoS ONE* **6**, e21800 (2011). — <http://revigo.irb.hr>
- Durinck, S. *et al.* Mapping identifiers for the integration of genomic datasets with the R/Bioconductor package biomaRt. *Nat. Protoc.* **4**, 1184–1191 (2009).
- Wu, T. *et al.* clusterProfiler 4.0: a universal enrichment tool for interpreting omics data. *Innovation* **2**, 100141 (2021). — book: <https://yulab-smu.top/biomedical-knowledge-mining-book/>
- Liberzon, A. *et al.* The Molecular Signatures Database hallmark gene set collection. *Cell Syst.* **1**, 417–425 (2015). — <https://www.gsea-msigdb.org/gsea/msigdb>
- Kuleshov, M. V. *et al.* Enrichr: a comprehensive gene set enrichment analysis web server 2016 update. *Nucleic Acids Res.* **44**, W90–W97 (2016). — <https://maayanlab.cloud/Enrichr/>
- Chen, H. & Boutros, P. C. VennDiagram: a package for the generation of highly-customizable Venn and Euler diagrams in R. *BMC Bioinformatics* **12**, 35 (2011); eulerr — <https://jolars.github.io/eulerr/>

---

[← Differential analysis](05_differential_analysis.md) · [Index](../README.md) · [Next: Motif analysis →](07_motif_analysis.md)
