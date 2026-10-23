# 3 · CUT&RUN and ChIP-seq alignment → BAM → bigWig

[← QC & trimming](02_qc_trimming.md) · [Index](../README.md) · [Next: Peak calling →](04_peak_calling.md)

```text
FASTQ (R1, R2)
   │  bowtie2 --local --very-sensitive-local -I 10 -X 700
   ▼
sample.bam ──► samtools view -f 2 (proper pairs) ──► sort ──► sample.map.sort.bam
   │                                                         │ Picard MarkDuplicates (remove)
   │                                                         ▼
   │                                                sample.map.dedup.bam (+ .bai)
   ├─► samtools flagstat                                     │
   │                                                         ├─► bamCoverage (CPM, blacklist) ──► sample.cpm.bw
   │                                                         └─► MACS2 (chapter 4)
```

---

## 3.1 Per-sample alignment script

`align_cutrun.sh <SAMPLE> <SPECIES> <WORKDIR> <DATADIR>`

```bash
#!/bin/bash
## CUT&RUN alignment (paired-end), Henikoff-lab bowtie2 settings
set -euo pipefail

module load FastQC-0.11.2 picard/1.96 bowtie2/2.3.4.1
module use.own /path/to/lab/privatemodules
module load samtools-1.9          # samtools 1.1 has known bugs; use ≥ 1.9

Sample_ID=$1
SPECIES=$2                        # hg19 | hg38 | mm10 …
WORKDIR=$3
DATADIR=$4
THREADS=8

LIBRARIES=/path/to/indexes
GNM_IDX=${LIBRARIES}/${SPECIES}
BOWTIE_IDX=${GNM_IDX}/bowtie2/${SPECIES}
BLACKLIST=${GNM_IDX}/${SPECIES}-blacklist.clean.bed

cd ${WORKDIR}
mkdir -p BAM Wig Peak ReadCount Log

fq1=${DATADIR}/${Sample_ID}_R1.fastq.gz
fq2=${DATADIR}/${Sample_ID}_R2.fastq.gz

## 1. Align ---------------------------------------------------------------
# --local --very-sensitive-local : soft-clip adapters, high sensitivity
# --no-mixed --no-discordant      : keep concordant pairs only
# -I 10 -X 700                    : fragment length 10–700 bp
bowtie2 -p ${THREADS} --local --very-sensitive-local --no-unal --no-mixed \
        --no-discordant --phred33 -I 10 -X 700 \
        -x ${BOWTIE_IDX} -1 ${fq1} -2 ${fq2} \
        2> BAM/${Sample_ID}.bowtie2.log \
  | samtools view -b -h -@ 2 -o BAM/${Sample_ID}.bam -

## 2. Keep properly paired reads, sort ------------------------------------
samtools view -b -f 2 -@ 2 BAM/${Sample_ID}.bam \
  | samtools sort -@ 4 -O BAM -o BAM/${Sample_ID}.map.sort.bam -
# add  -q 10  to the view call to also drop multi-mapping reads (MAPQ < 10)

## 3. Remove PCR duplicates -----------------------------------------------
picard MarkDuplicates \
  I=BAM/${Sample_ID}.map.sort.bam \
  O=BAM/${Sample_ID}.map.dedup.bam \
  M=BAM/${Sample_ID}.dedupMetrics.txt \
  VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=true \
  2> BAM/${Sample_ID}.dedup.log
samtools index BAM/${Sample_ID}.map.dedup.bam

## 4. Alignment statistics ------------------------------------------------
samtools flagstat BAM/${Sample_ID}.bam           > BAM/${Sample_ID}.bam_flagstat
samtools flagstat BAM/${Sample_ID}.map.dedup.bam > BAM/${Sample_ID}.map.dedup.bam_flagstat

## 5. CPM-normalised bigWig ----------------------------------------------
bamCoverage -b BAM/${Sample_ID}.map.dedup.bam -o Wig/${Sample_ID}.cpm.bw \
  --binSize 5 --smoothLength 25 --normalizeUsing CPM \
  --extendReads --blackListFileName ${BLACKLIST} \
  --ignoreForNormalization chrM --skipNAs -p ${THREADS}

## 6. Free disk space (after checking the logs) ---------------------------
# rm BAM/${Sample_ID}.bam BAM/${Sample_ID}.map.sort.bam
```

> If you trimmed reads (chapter 2), point `fq1`/`fq2` at the `*_trim_1P` / `*_trim_2P` files.

### Example output

`BAM/<sample>.bowtie2.log`

```text
15234567 reads; of these:
  15234567 (100.00%) were paired; of these:
    1203511 (7.90%) aligned concordantly 0 times
    11634892 (76.37%) aligned concordantly exactly 1 time
    2396164 (15.73%) aligned concordantly >1 times
92.10% overall alignment rate
```

`BAM/<sample>.map.dedup.bam_flagstat`

```text
24120314 + 0 in total (QC-passed reads + QC-failed reads)
0 + 0 secondary
0 + 0 supplementary
0 + 0 duplicates
24120314 + 0 mapped (100.00% : N/A)
24120314 + 0 paired in sequencing
24120314 + 0 properly paired (100.00% : N/A)
```

*(Values are illustrative.)* As a rule of thumb, CUT&RUN libraries need roughly 5–10 M unique fragments for histone marks; transcription factors can need fewer.

---

## 3.2 Run all samples in parallel

```bash
#!/bin/bash
## submit_all.sh
codepath=/path/to/Project/Code
workdir=/path/to/Project/Data
datadir=/path/to/Project/Fastq
species=hg19
cd ${workdir} && mkdir -p Log

# QC once
bsub -n 8 -M 32000 -e Log/qc.e -o Log/qc.o bash ${codepath}/run_qc.sh ${workdir} ${datadir}

# one alignment job per sample (skip the header line)
for sample in $(awk 'NR>1 {print $1}' sample_list_all.txt); do
  bsub -n 8 -M 64000 -R "span[hosts=1]" \
       -J align_${sample} -e Log/${sample}.e -o Log/${sample}.o \
       bash ${codepath}/align_cutrun.sh ${sample} ${species} ${workdir} ${datadir}
done
```

To re-run only failed samples, list them one per line in `sample_list_rerun.txt` and loop over `$(cat sample_list_rerun.txt)` instead.

Collect a summary table once all jobs have finished:

```bash
printf "sample\ttotal_pairs\toverall_rate\tdedup_reads\tdup_rate\n" > QC/align_summary.tsv
for s in $(awk 'NR>1{print $1}' sample_list_all.txt); do
  tot=$(head -1 BAM/${s}.bowtie2.log | cut -d' ' -f1)
  rate=$(grep overall BAM/${s}.bowtie2.log | cut -d' ' -f1)
  dedup=$(head -1 BAM/${s}.map.dedup.bam_flagstat | cut -d' ' -f1)
  # find PERCENT_DUPLICATION by name (its column number differs between Picard versions)
  dup=$(awk -F'\t' '/^LIBRARY/{for(i=1;i<=NF;i++) if($i=="PERCENT_DUPLICATION") c=i; getline; print $c}' \
        BAM/${s}.dedupMetrics.txt)
  printf "%s\t%s\t%s\t%s\t%s\n" $s $tot $rate $dedup $dup >> QC/align_summary.tsv
done
column -t QC/align_summary.tsv
```

---

## 3.3 Two ways to make a bigWig

**A. deepTools `bamCoverage`** (recommended; used above)

```bash
bamCoverage -b ${sample}.bam -o ${sample}.cpm.bw \
  -bl ${BLACKLIST} -p 4 --normalizeUsing CPM \
  --binSize 5 --smoothLength 10 --extendReads \
  --ignoreForNormalization chrM --ignoreDuplicates --skipNAs --minMappingQuality 5
```

| Option | Why |
|---|---|
| `--normalizeUsing CPM` | counts per million mapped reads; compares libraries of different depth (RPKM/BPM/RPGC also available) |
| `--extendReads` (`-e`) | paired-end: fills the whole fragment, not just the read |
| `-bl` | zero out ENCODE blacklist regions |
| `--binSize / --smoothLength` | resolution vs. file size; 5–10 bp bins for TFs, 25–50 bp for broad marks |

**B. bedtools → bedGraph → bigWig** (UCSC tools)

```bash
GSIZE=${GNM_IDX}/${SPECIES}.genomesize

# paired-end: one interval per fragment (name-sorted BAM required)
samtools sort -n -@ 4 -o ${sample}.nsort.bam ${sample}.map.dedup.bam
bedtools bamtobed -bedpe -i ${sample}.nsort.bam \
  | awk 'BEGIN{OFS="\t"} $1==$4 && $6>$2 {print $1,$2,$6,$6-$2}' \
  | grep -v -e random -e chrUn -e chrM -e _ \
  | sort -k1,1 -k2,2n > ${sample}.frag.bed

# remove fragments overlapping the blacklist
bedtools intersect -v -a ${sample}.frag.bed -b ${BLACKLIST} > ${sample}.frag.deblack.bed

scale=$(wc -l < ${sample}.frag.deblack.bed | awk '{print 1000000/$1}')   # RPM
bedtools genomecov -bg -scale ${scale} -i ${sample}.frag.deblack.bed -g ${GSIZE} \
  | sort -k1,1 -k2,2n > ${sample}.rpm.bedgraph
bedGraphToBigWig ${sample}.rpm.bedgraph ${GSIZE} ${sample}.rpm.bw
```

> `bedGraphToBigWig` requires a bedGraph sorted with `sort -k1,1 -k2,2n`. Always convert the **sorted** file.
>
> `bedtools subtract` trims a fragment where it overlaps the blacklist. `bedtools intersect -v` drops it entirely, which is usually what you want.

---

## 3.4 Spike-in (E. coli / yeast) normalisation

When global levels may change (e.g. inhibitor treatments), normalise to spike-in reads rather than library size:

```bash
# align the same reads to the spike-in genome (repeat-masked to avoid cross-mapping)
bowtie2 -p 8 --local --very-sensitive-local --no-unal --no-mixed --no-discordant \
        --no-overlap --no-dovetail -I 10 -X 700 \
        -x ${LIBRARIES}/Ecoli/Ecoli -1 ${fq1} -2 ${fq2} 2> BAM/${id}.spike.log \
  | samtools view -b -f 2 -o BAM/${id}.spike.bam -

spike=$(samtools view -c -f 64 BAM/${id}.spike.bam)       # read-1 count = fragments
scale=$(echo "10000 / ${spike}" | bc -l)                  # constant 10,000 is arbitrary

bamCoverage -b BAM/${id}.map.dedup.bam -o Wig/${id}.spikein.bw \
  --scaleFactor ${scale} --binSize 10 --extendReads -bl ${BLACKLIST} -p 8
```

The Henikoff-lab `spike_in_calibration.csh` does the same calculation on fragment BED files and can also restrict fragment sizes (e.g. 100–300 bp for nucleosomal signal).

---

## 3.5 Single-end ChIP-seq variant

Differences from CUT&RUN: single-end reads, end-to-end alignment, and an input control.

```bash
bowtie2 -p 8 --very-sensitive -x ${BOWTIE_IDX} -U ${id}.fastq.gz 2> BAM/${id}.bowtie2.log \
  | samtools view -b -q 10 - | samtools sort -@ 4 -o BAM/${id}.sort.bam -
picard MarkDuplicates I=BAM/${id}.sort.bam O=BAM/${id}.dedup.bam \
  M=BAM/${id}.dedupMetrics.txt REMOVE_DUPLICATES=true
samtools index BAM/${id}.dedup.bam
bamCoverage -b BAM/${id}.dedup.bam -o Wig/${id}.cpm.bw --normalizeUsing CPM \
  --extendReads 200 -bl ${BLACKLIST} -p 8      # extend single-end reads to ~fragment size
```

> Older samtools (0.1.x) used `samtools sort in.bam prefix`. With samtools ≥ 1.3, use `samtools sort -o out.bam in.bam`.

---

## References and resources

- Skene, P. J. & Henikoff, S. An efficient targeted nuclease strategy for high-resolution mapping of DNA binding sites. *eLife* **6**, e21856 (2017). (CUT&RUN)
- Meers, M. P., Bryson, T. D., Henikoff, J. G. & Henikoff, S. Improved CUT&RUN chromatin profiling tools. *eLife* **8**, e46314 (2019). (spike-in calibration)
- Langmead, B. & Salzberg, S. L. Fast gapped-read alignment with Bowtie 2. *Nat. Methods* **9**, 357–359 (2012). — <https://bowtie-bio.sourceforge.net/bowtie2/>
- Ramírez, F. *et al.* deepTools2: a next generation web server for deep-sequencing data analysis. *Nucleic Acids Res.* **44**, W160–W165 (2016). — <https://deeptools.readthedocs.io>
- Quinlan, A. R. & Hall, I. M. BEDTools: a flexible suite of utilities for comparing genomic features. *Bioinformatics* **26**, 841–842 (2010). — <https://bedtools.readthedocs.io>
- Picard `MarkDuplicates` — <https://broadinstitute.github.io/picard/>
- CUT&Tag data processing and analysis tutorial (Henikoff lab; applies to CUT&RUN) — <https://yezhengstat.github.io/CUTTag_tutorial/>

---

[← QC & trimming](02_qc_trimming.md) · [Index](../README.md) · [Next: Peak calling →](04_peak_calling.md)
