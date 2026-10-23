# 2 · Read QC and trimming

[← Data transfer](01_data_transfer.md) · [Index](../README.md) · [Next: CUT&RUN / ChIP alignment →](03_alignment_cutrun_chip.md)

| Step | Tool | Input | Output |
|---|---|---|---|
| Per-file QC | FastQC | `*.fastq.gz` | `QC/*_fastqc.html` |
| Summary report | MultiQC | FastQC / aligner / Picard logs | `QC/MultiQC/multiqc_fastq.html` |
| Adapter trimming (if needed) | Trimmomatic or cutadapt | FASTQ | trimmed FASTQ |
| Fragment size | Picard `CollectInsertSizeMetrics` | BAM | histogram PDF + table |
| Library strandedness (RNA-seq) | RSeQC `infer_experiment.py` | BAM + BED12 | strand fractions |

---

## 2.1 FastQC + MultiQC

`run_qc.sh` takes the working and data directories as arguments, so the same script works for any project:

```bash
#!/bin/bash
## run_qc.sh  — usage: bash run_qc.sh <WORKDIR> <DATADIR>
module load FastQC-0.11.2

WORKDIR=$1
DATADIR=$2
cd ${WORKDIR} && mkdir -p QC

find ${DATADIR}/ -name '*.fastq*' | xargs fastqc -t 8 -o QC
multiqc -d QC -o QC/MultiQC -i multiqc_fastq     # -i = report title / file name

rm -f QC/*_fastqc.zip                            # keep the HTML only
```

```bash
bsub -n 8 -M 32000 -e log_qc.e -o log_qc.o bash Code/run_qc.sh $PWD/Data /path/to/Project/Fastq
```

**What to look at**

| FastQC module | CUT&RUN / ChIP expectation | Action if it fails |
|---|---|---|
| Per-base quality | Phred > 28 across most of the read | trim low-quality 3′ ends |
| Adapter content | < 5 % at read end | trim adapters (§2.2) |
| Sequence duplication | high in CUT&RUN (low input) | expected; removed after alignment |
| Overrepresented sequences | adapter dimers | trim |

> Run MultiQC again after alignment (`multiqc -d BAM QC`). It also parses bowtie2 logs, `samtools flagstat` and Picard metrics, so all alignment statistics end up in one table.

---

## 2.2 Adapter trimming

Only trim when MultiQC shows adapter read-through. CUT&RUN fragments are often shorter than the read length, so check every library. (Henikoff's CUT&RUN protocol skips trimming because bowtie2 `--local` soft-clips adapters.)

**Trimmomatic (paired-end)**

```bash
TRIM=/path/to/tools/Trimmomatic-0.38
java -jar ${TRIM}/trimmomatic-0.38.jar PE -threads 8 -phred33 \
  ${id}_R1.fastq.gz ${id}_R2.fastq.gz \
  -baseout Fastq_Trim/${id}_trim.fastq.gz \
  ILLUMINACLIP:${TRIM}/adapters/TruSeq3-PE.fa:2:15:4:4:true \
  LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:25
# → ${id}_trim_1P.fastq.gz / _2P (paired) and _1U / _2U (orphans)
```

> For single-end ChIP-seq use `SE` with `TruSeq3-SE.fa`.

**cutadapt (explicit adapter sequence)**

```bash
conda activate cutadapt
# Illumina TruSeq adapter on both reads, quality-trim at Q20, drop reads < 20 nt
cutadapt -j 8 -q 20 -m 20 \
  -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
  -o ${id}_R1.trim.fastq.gz -p ${id}_R2.trim.fastq.gz \
  ${id}_R1.fastq.gz ${id}_R2.fastq.gz > ${id}.cutadapt.log
```

---

## 2.3 Fragment-size distribution

```bash
java -jar picard.jar CollectInsertSizeMetrics \
  I=BAM/${id}.map.dedup.bam \
  O=QC/${id}.insert_size_metrics.txt \
  H=QC/${id}.insert_size_histogram.pdf \
  M=0.5
```

Good CUT&RUN libraries for histone marks show a **nucleosomal ladder**: a mono-nucleosome peak near 150–200 bp and a smaller di-nucleosome peak near 300–350 bp (**Fig. 1b**). In IgG controls, an extra population below 120 bp usually comes from untethered MNase cutting.

---

## 2.4 Library strandedness (RNA-seq)

```bash
infer_experiment.py -r /path/to/indexes/mm10/gencode/anno_mm10_gencode.bed12 -i sample.bam
```

Example output:

```text
This is PairEnd Data
Fraction of reads failed to determine: 0.0102
Fraction of reads explained by "1++,1--,2+-,2-+": 0.0153
Fraction of reads explained by "1+-,1-+,2++,2--": 0.9745
```

Here ~97 % of reads fit `1+-,1-+` → **reverse-stranded** (e.g. Illumina TruSeq Stranded). Use `featureCounts -s 2` or STAR's `ReadsPerGene` column 4.

| Pattern dominating | Library | featureCounts `-s` |
|---|---|---|
| ~50 % / 50 % | unstranded | 0 |
| `1++,1--,2+-,2-+` | stranded (forward) | 1 |
| `1+-,1-+,2++,2--` | reverse-stranded | 2 |

---

## Example output

![Fig. 1](figures/fig1_qc.png)

**Fig. 1 | Library-level QC for a CUT&RUN experiment.** **a**, Read pairs per library split into unique deduplicated, PCR-duplicate and unaligned reads (from bowtie2 logs and Picard `MarkDuplicates` metrics). **b**, Fragment-length distribution from paired-end alignments; histone libraries show mono- and di-nucleosome peaks, and IgG shows extra sub-nucleosomal fragments. **c**, Fraction of reads in peaks (FRiP) per condition (bars, mean; dots, replicates); the dashed line marks the ENCODE guideline of 0.1. *Simulated data for illustration* (`scripts/make_example_figures.py`).

Compute FRiP for one sample:

```bash
total=$(samtools view -c -F 0x904 BAM/${id}.map.dedup.bam)
inpk=$(samtools view -c -F 0x904 -L Peak/${id}.final.bed BAM/${id}.map.dedup.bam)
echo -e "${id}\t$(echo "scale=4; ${inpk}/${total}" | bc)"
```

---

## References and resources

- Andrews, S. FastQC: a quality control tool for high throughput sequence data — <https://www.bioinformatics.babraham.ac.uk/projects/fastqc/>
- Ewels, P. *et al.* MultiQC: summarize analysis results for multiple tools and samples in a single report. *Bioinformatics* **32**, 3047–3048 (2016). — <https://multiqc.info>
- Bolger, A. M., Lohse, M. & Usadel, B. Trimmomatic: a flexible trimmer for Illumina sequence data. *Bioinformatics* **30**, 2114–2120 (2014).
- Martin, M. Cutadapt removes adapter sequences from high-throughput sequencing reads. *EMBnet.journal* **17**, 10–12 (2011). — <https://cutadapt.readthedocs.io>
- Picard tools — <https://broadinstitute.github.io/picard/>
- Wang, L., Wang, S. & Li, W. RSeQC: quality control of RNA-seq experiments. *Bioinformatics* **28**, 2184–2185 (2012). — <https://rseqc.sourceforge.net>
- Landt, S. G. *et al.* ChIP-seq guidelines and practices of the ENCODE and modENCODE consortia. *Genome Res.* **22**, 1813–1831 (2012). (FRiP and library QC)

---

[← Data transfer](01_data_transfer.md) · [Index](../README.md) · [Next: CUT&RUN / ChIP alignment →](03_alignment_cutrun_chip.md)
