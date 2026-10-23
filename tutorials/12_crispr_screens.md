# 12 · Pooled CRISPR screens

[← Hi-C](11_hic.md) · [Index](../README.md) · [Next: scRNA-seq →](13_scrnaseq.md)

A pooled knockout screen compares sgRNA abundance between a reference (plasmid or day-0 cells) and the selected population (e.g. tumour-infiltrating T cells). Genes whose guides are **enriched** confer an advantage when knocked out; genes whose guides are **depleted** are required. The workflow:

1. **Count** guides per sample and check library QC (mapping rate, zero-count guides, Gini index).
2. **Test** each gene (MAGeCK RRA or MLE; RIGER as an alternative).
3. **Plot**: guide-count distributions, QC, replicate concordance, gene rank plot, sgRNA-level plot of top hits.

## 12.1 MAGeCK: count and test

```bash
conda create -n mageck -c bioconda -c conda-forge mageck
conda activate mageck

# library.csv: sgRNA_id,sequence,gene   ;   ntc.txt: one non-targeting sgRNA id per line
mageck count -l library.csv -n Screen \
  --sample-label Plasmid,Day0,Spleen_R1,Spleen_R2,Tumour_R1,Tumour_R2 \
  --fastq Plasmid.fq.gz Day0.fq.gz Spleen_R1.fq.gz Spleen_R2.fq.gz Tumour_R1.fq.gz Tumour_R2.fq.gz \
  --norm-method control --control-sgrna ntc.txt
# → Screen.count.txt, Screen.count_normalized.txt, Screen.countsummary.txt

# robust rank aggregation (RRA): tumour vs day 0, replicates pooled
mageck test -k Screen.count.txt -t Tumour_R1,Tumour_R2 -c Day0 -n Tumour_vs_Day0 \
  --norm-method control --control-sgrna ntc.txt
# → Tumour_vs_Day0.gene_summary.txt, Tumour_vs_Day0.sgrna_summary.txt
```

Library QC targets (`Screen.countsummary.txt`):

| Metric | Plasmid / day 0 | In vivo samples |
|---|---|---|
| Reads mapped to library | > 65 % | > 65 % |
| Zero-count sgRNAs | < 1 % | a few % expected after bottlenecks |
| Gini index | < 0.1–0.2 | rises with selection; very high (> 0.5) means too few cells per guide |
| Coverage | ≥ 300–500 reads per guide | as high as cell numbers allow |

## 12.2 Plots in R

```r
library(data.table); library(dplyr); library(tidyr); library(ggplot2); library(ggrepel)
# theme_nature() is defined in chapter 5

## a. guide-count distribution per sample
cnt <- fread("Screen.count_normalized.txt") |>
  pivot_longer(-c(sgRNA, Gene), names_to = "Sample", values_to = "Count")
ggplot(cnt, aes(log2(Count + 1), colour = Sample)) +
  geom_density(linewidth = 0.4) +
  labs(x = expression(log[2]*"(normalized sgRNA count + 1)"), y = "Density", colour = NULL) +
  theme_nature()

## b. library QC: Gini index and zero-count guides
qc <- fread("Screen.countsummary.txt")
ggplot(qc, aes(GiniIndex, reorder(Label, GiniIndex))) +
  geom_col(fill = "grey60", width = 0.65) +
  geom_text(aes(label = sprintf("%.1f%% zero", 100 * Zerocounts / TotalsgRNAs)),
            hjust = -0.1, size = 1.8) +
  geom_vline(xintercept = 0.2, linetype = 2, linewidth = 0.2) +
  labs(x = "Gini index", y = NULL) + theme_nature()

## c. replicate concordance (gene-level median log2FC per replicate)
raw <- fread("Screen.count_normalized.txt")
lfc <- raw |>
  mutate(R1 = log2((Tumour_R1 + 0.5) / (Day0 + 0.5)),
         R2 = log2((Tumour_R2 + 0.5) / (Day0 + 0.5))) |>
  group_by(Gene) |> summarise(R1 = median(R1), R2 = median(R2))
r <- cor(lfc$R1, lfc$R2)
ggplot(lfc, aes(R1, R2)) + geom_point(size = 0.3, colour = "grey60") +
  geom_abline(linetype = 2, linewidth = 0.2) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.5, size = 2,
           label = sprintf("Pearson r = %.2f", r)) +
  labs(x = expression("Gene"~log[2]*"FC, replicate 1"), y = expression("Gene"~log[2]*"FC, replicate 2")) +
  theme_nature()

## d. gene rank plot (MAGeCK RRA): enriched genes on the left, depleted on the right
gs <- fread("Tumour_vs_Day0.gene_summary.txt") |>
  transmute(Gene = id, LFC = `pos|lfc`, FDR = pmin(`pos|fdr`, `neg|fdr`)) |>
  filter(Gene != "NonTargeting") |>
  arrange(desc(LFC)) |> mutate(Rank = row_number(),
                         Hit = case_when(FDR < 0.05 & LFC > 0 ~ "Enriched",
                                         FDR < 0.05 & LFC < 0 ~ "Depleted", TRUE ~ "NS"))
lab <- bind_rows(slice_min(gs, LFC, n = 5), slice_max(gs, LFC, n = 6)) |> filter(Hit != "NS")
ggplot(gs, aes(Rank, LFC, colour = Hit)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -0.5, ymax = 0.5, fill = "grey95") +
  geom_hline(yintercept = 0, linewidth = 0.25) +
  geom_point(size = 0.4) +
  geom_text_repel(data = lab, aes(label = Gene), size = 2, fontface = "italic", colour = "black",
                  segment.size = 0.2, min.segment.length = 0, box.padding = 0.3,
                  nudge_x = ifelse(lab$LFC > 0, 150, -150), direction = "y", max.overlaps = Inf) +
  scale_colour_manual(values = c(Enriched = "#D55E00", Depleted = "#0072B2", NS = "grey75")) +
  labs(x = "Gene rank", y = expression("Median sgRNA"~log[2]*"FC")) +
  theme_nature() + theme(legend.position = "none")

## e. sgRNA-level plot for top hits vs non-targeting controls
sg  <- fread("Tumour_vs_Day0.sgrna_summary.txt")
ntc <- sg |> filter(Gene == "NonTargeting")
show <- c("Dhx37", "Pdcd1", "Cblb", "Myb", "Tbx21", "Cd28")
d <- sg |> filter(Gene %in% show) |> mutate(Gene = factor(Gene, levels = rev(show)))
ggplot(d, aes(LFC, Gene)) +
  annotate("rect", xmin = quantile(ntc$LFC, 0.025), xmax = quantile(ntc$LFC, 0.975),
           ymin = -Inf, ymax = Inf, fill = "grey92") +
  geom_vline(xintercept = 0, linewidth = 0.25) +
  geom_jitter(aes(colour = LFC > 0), height = 0.12, size = 1) +
  stat_summary(fun = median, geom = "crossbar", width = 0.5, linewidth = 0.2) +
  scale_colour_manual(values = c(`TRUE` = "#D55E00", `FALSE` = "#0072B2"), guide = "none") +
  labs(x = expression("sgRNA"~log[2]*"FC"), y = NULL) +
  theme_nature() + theme(axis.text.y = element_text(face = "italic"))
```

> A real hit has **several concordant guides** clearly outside the non-targeting range (grey band). A single extreme guide is usually an off-target or bottleneck artefact. With strong bottlenecks (in vivo screens), use more replicates or animals and consider `mageck mle` with a replicate-aware design matrix.

## 12.3 RIGER (gene ranking from guide-level scores)

Input (tab-separated, this exact header order):

```text
Construct   GeneSymbol   NormalizedScore   Construct Rank   HairpinWeight
sgTox_1     Tox          3.21              1                1
```

```bash
rigerj=/path/to/tools/rigerj/target/rigerj-2.0.2-assembly.jar
java -jar ${rigerj} -inputFile Input_ForRIGER_d15_spleen.txt \
     -scoringMethod SecondBestRank -outputFile RIGER_d15_spleen_SBR
# -scoringMethod: KSbyScore | WtSum | SecondBestRank
```

Output columns: `Gene Rank`, `Gene Name`, `Score`, `p-value`, `p-value rank`, `Hairpin ranks`. For negative selection, rank guides in the opposite direction.

## 12.4 CRISPRO (map scores to protein domains/structure)

```bash
source activate crispro
crispro --scores Donovan2017_inputCRISPRO.csv \
        -i Donovan2017_IDs.txt \
        -a annotations/CRISPRO.GRCh37.ensembl90.20180901.csv.gz \
        -o out_CRISPRO
conda deactivate
```

---

## Example output

![Fig. 11](figures/fig11_crispr_screen.png)

**Fig. 11 | Pooled in vivo CRISPR knockout screen.** **a**, Distribution of normalized sgRNA counts per sample; selection in vivo widens the distribution compared with the plasmid and day-0 library. **b**, Library evenness (Gini index; dashed line, 0.2) and percentage of zero-count sgRNAs per sample. **c**, Concordance of gene-level log₂ fold changes (tumour versus day 0) between biological replicates; coloured points are significant hits. **d**, Genes ranked by median sgRNA log₂ fold change, from most enriched (left) to most depleted (right); significantly enriched (orange) and depleted (blue) genes (MAGeCK RRA, FDR < 0.05) are labelled; grey band, |log₂FC| < 0.5. **e**, Individual sgRNA log₂ fold changes for top hits (dots; bar, median) against the distribution of non-targeting controls (violin; grey band, central 95%). *Simulated data for illustration.*

---

## References and resources

- Li, W. *et al.* MAGeCK enables robust identification of essential genes from genome-scale CRISPR/Cas9 knockout screens. *Genome Biol.* **15**, 554 (2014). — <https://sourceforge.net/p/mageck/wiki/Home/>
- Li, W. *et al.* Quality control, modeling, and visualization of CRISPR screens with MAGeCK-VISPR. *Genome Biol.* **16**, 281 (2015).
- Wang, B. *et al.* Integrative analysis of pooled CRISPR genetic screens using MAGeCKFlute. *Nat. Protoc.* **14**, 756–780 (2019). — <https://bioconductor.org/packages/MAGeCKFlute>
- Luo, B. *et al.* Highly parallel identification of essential genes in cancer cells. *Proc. Natl Acad. Sci. USA* **105**, 20380–20385 (2008). (RIGER)
- He, W. *et al.* De novo identification of essential protein domains from CRISPR-Cas9 tiling-sgRNA knockout screens. *Nat. Commun.* **10**, 4541 (2019). (CRISPRO) — <https://gitlab.com/bauerlab/crispro>
- Hart, T. & Moffat, J. BAGEL: a computational framework for identifying essential genes from pooled library screens. *BMC Bioinformatics* **17**, 164 (2016).
- Dong, M. B. *et al.* Systematic immunotherapy target discovery using genome-scale in vivo CRISPR screens in CD8 T cells. *Cell* **178**, 1189–1204 (2019). (example of an in vivo T-cell screen)

---

[← Hi-C](11_hic.md) · [Index](../README.md) · [Next: scRNA-seq →](13_scrnaseq.md)
