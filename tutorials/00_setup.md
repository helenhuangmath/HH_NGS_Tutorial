# Part 0 · Setup: HPC, project layout and software

[← Index](../README.md) · [Next: Data transfer →](01_data_transfer.md)

All pipelines in this tutorial were written for an **HPC cluster running the IBM LSF scheduler**. Commands that start with `bsub` submit a job; everything else can run interactively on a compute node.

---

## 0.1 Submitting jobs with LSF

```bash
# 8 cores, 64 GB RAM, all cores on one node, separate log files per job (%J = job ID)
bsub -n 8 -M 64000 -R "rusage[mem=64000] span[hosts=1]" \
     -e log_%J.e -o log_%J.o \
     bash my_script.sh arg1 arg2

bjobs                 # list running jobs
bjobs -l <jobID>      # details / memory usage
bkill <jobID>         # cancel
```

> **Tip** – match the thread count of the tool (`-p`, `-T`, `--threads`) to `bsub -n`. Asking for 8 cores and running 16 threads slows every job on the node.

---

## 0.2 Recommended project layout

```text
Project/
├── Code/                    # scripts (version-controlled)
├── Data/
│   ├── Fastq/               # raw reads (read-only)
│   ├── QC/                  # FastQC / MultiQC
│   ├── BAM/                 # alignments, flagstat, dedup metrics
│   ├── Wig/                 # bigWig tracks
│   ├── Peak/                # MACS2 output, consensus peaks
│   ├── ReadCount/           # count matrices
│   ├── Log/                 # LSF logs
│   └── sample_list_*.txt    # sample sheets
└── Results/                 # DESeq2, figures, tables
```

Create the sample sheet from the FASTQ names (one line per sample, `-` = no control):

```bash
DATADIR=/path/to/Project/Fastq

ls ${DATADIR}/*fastq* | xargs -n 1 basename \
  | sed -E 's/_R[12](_001)?\.fastq(\.gz)?//' | sort -u \
  | awk 'BEGIN{OFS="\t"; print "sample","control"} {print $1,"-"}' \
  > Data/sample_list_all.txt

# subsets by target (keep the header line)
grep -e sample -e IGG -e H3K27ac Data/sample_list_all.txt > Data/sample_list_H3K27ac.txt
```

Example `sample_list_all.txt`:

```text
sample                   control
DMSO_H3K27ac_Rep1        -
DMSO_H3K27ac_Rep2        -
cLTP_H3K27ac_Rep1        -
DMSO_IGG_Rep1            -
```

---

## 0.3 Reference files

Every pipeline expects the same layout under one `indexes/` folder so that only `SPECIES` changes between projects:

```text
indexes/
└── hg19/                               (also mm10, mm9, hg38 …)
    ├── bowtie2/hg19.*.bt2              # bowtie2 index
    ├── STAR_Gencode/                   # STAR index (RNA-seq)
    ├── hg19.fa                         # genome FASTA
    ├── hg19.genomesize                 # chrom<TAB>length (UCSC chrom.sizes)
    ├── hg19-blacklist.clean.bed        # ENCODE blacklist v2
    └── gencode.v42.annotation.gtf(.bed)
```

| Variable | Typical value |
|---|---|
| Effective genome size (MACS2 `-g`) | human `2.7e9` (or `hs`), mouse `1.87e9` (or `mm`) |
| Blacklist | ENCODE blacklist v2 (Amemiya *et al.*, *Sci. Rep.* 2019) |

---

## 0.4 Modules and private R libraries

```bash
# shared (group) module directory
module use.own /path/to/lab/privatemodules
module avail
module load samtools-1.9 bowtie2/2.3.4.1 picard/1.96 FastQC-0.11.2
module load R/4.0.2 hdf5/1.8.21
```

Point R at a shared package library (run inside R):

```r
.libPaths()                                        # check current paths
.libPaths("/path/to/lab/R_libs")   # use shared library
# or prepend a personal library
.libPaths(c("~/R/userLibrary", .libPaths()))
```

Install a **specific** package version (needed when a newer release breaks old code):

```r
install.packages("remotes")
remotes::install_version("emmeans",  version = "1.4.5")
remotes::install_version("spatstat", version = "1.64-1")   # Seurat v3 dependency
remotes::install_version("Seurat",   version = "3.2.3")
```

> To run several R versions with RStudio, see Posit's guides on installing multiple R versions and selecting one in RStudio Workbench.

---

## 0.5 Conda environments

Install Miniconda into a shared project folder (home directories are small on most clusters):

```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh -b -u -p /path/to/tools/miniconda3
source /path/to/tools/miniconda3/etc/profile.d/conda.sh
```

One environment per tool family avoids dependency clashes:

```bash
conda create -n cutadapt   -c bioconda -c conda-forge cutadapt
conda create -n ngs        -c bioconda -c conda-forge \
      macs2 deeptools bedtools samtools subread multiqc homer meme
conda create -n hicexplorer -c bioconda -c conda-forge hicexplorer
conda create -n r42        -c conda-forge -c bioconda r-base=4.2 \
      bioconductor-complexheatmap bioconductor-deseq2 r-seurat
conda activate ngs
```

> R ≥ 4.2 is recommended for ComplexHeatmap; older R versions produced errors with recent releases.

---

## 0.6 Tool-specific installs

<details>
<summary><b>HOMER + custom GENCODE genome</b></summary>

HOMER's bundled annotation can be outdated. Load a GENCODE genome so annotation uses current gene models:

```bash
export PATH=$PATH:/path/to/tools/homer/bin
ref=/path/to/indexes/mm10/Gencode_vM24

# all genes; -gid = report gene_id instead of transcript_id
loadGenome.pl -name mm10_GencodeV24 -org null \
  -fasta ${ref}/GRCm38.primary_assembly.genome.fa \
  -gtf   ${ref}/gencode.vM24.annotation.gtf -gid

# optional: drop processed pseudogenes / Gm models, or keep protein-coding only
grep -v -e processed_pseudogene -e 'gene_name "Gm' ${ref}/gencode.vM24.primary_assembly.annotation.gtf \
  > ${ref}/gencode.vM24.nopseudo.gtf
grep protein_coding ${ref}/gencode.vM24.primary_assembly.annotation.gtf \
  > ${ref}/gencode.vM24.codinggene.gtf
loadGenome.pl -name mm10_GencodeV24_codinggene -org null \
  -fasta ${ref}/GRCm38.primary_assembly.genome.fa \
  -gtf   ${ref}/gencode.vM24.codinggene.gtf -gid
```
Then call `annotatePeaks.pl peaks.bed mm10_GencodeV24`.
</details>

<details>
<summary><b>HiC-Pro (v3, Python 3)</b></summary>

```bash
git clone https://github.com/nservant/HiC-Pro.git && cd HiC-Pro
pip install pysam bx-python numpy scipy iced
nano config-install.txt        # set PREFIX, scheduler (LSF), paths to bowtie2/samtools
make configure && make install
export PATH=$PATH:/path/to/tools/HiC-Pro/bin
HiC-Pro --help
```
Older HiC-Pro (≤ 2.11) required Python 2.7; v3 runs on Python ≥ 3.7.
</details>

<details>
<summary><b>FitHiChIP</b></summary>

Prerequisites: Python 3 (`networkx`, `hic-straw`, `cooler`), R ≥ 3.4 (`optparse ggplot2 data.table splines fdrtool parallel tools plyr dplyr GenomicRanges edgeR`), bedtools ≥ 2.26, samtools/htslib ≥ 1.6, bowtie2, HiC-Pro, MACS2.

```bash
pip install networkx hic-straw cooler
git clone https://github.com/ay-lab/FitHiChIP.git
bash FitHiChIP/FitHiChIP_HiCPro.sh -C /path/to/configfile
```
</details>

<details>
<summary><b>HiCExplorer</b></summary>

```bash
conda create -n hicexplorer -c bioconda -c conda-forge hicexplorer
conda activate hicexplorer
```
(Older releases required Python 3.6 – a dedicated environment avoids conflicts.)
</details>

<details>
<summary><b>HINT (RGT) footprinting</b></summary>

```bash
pip install RGT                     # creates ~/rgtdata
cd ~/rgtdata
python setupGenomicData.py --mm10   # also --hg19, --hg38 …
python setupLogoData.py --all       # motif logos
```
</details>

<details>
<summary><b>Seurat / SCENIC</b></summary>

```r
# Seurat (current)
install.packages("Seurat")
# Seurat v3.2.3 (legacy code)
remotes::install_version("Seurat", version = "3.2.3")

# SCENIC
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("AUCell", "RcisTarget", "GENIE3",
                       "zoo", "mixtools", "DT", "NMF", "ComplexHeatmap",
                       "R2HTML", "Rtsne", "doMC", "doRNG"))
remotes::install_github("aertslab/SCopeLoomR", build_vignettes = TRUE)
remotes::install_github("aertslab/SCENIC")
```
</details>

<details>
<summary><b>diffReps, RIGER (RigerJ)</b></summary>

```bash
# diffReps (Perl)
cpan Statistics::TTest Math::CDF Parallel::ForkManager Time::HiRes
git clone https://github.com/shenlab-sinai/diffreps.git && cd diffreps
perl Makefile.PL PREFIX=$HOME && make && make test && make install

# RigerJ (Java + Maven); skip javadoc to avoid build errors
export JAVA_HOME=/path/to/tools/jdk-12.0.2
export PATH=/path/to/tools/maven-3.6.1/bin:$JAVA_HOME/bin:$PATH
cd rigerj && mvn package -Dmaven.javadoc.skip=true
```
</details>

---

## References and resources

- IBM Spectrum LSF documentation (`bsub`, `bjobs`, resource strings) — <https://www.ibm.com/docs/en/spectrum-lsf>
- Grüning, B. *et al.* Bioconda: sustainable and comprehensive software distribution for the life sciences. *Nat. Methods* **15**, 475–476 (2018). — <https://bioconda.github.io>
- Huber, W. *et al.* Orchestrating high-throughput genomic analysis with Bioconductor. *Nat. Methods* **12**, 115–121 (2015). — <https://bioconductor.org>
- Amemiya, H. M., Kundaje, A. & Boyle, A. P. The ENCODE blacklist: identification of problematic regions of the genome. *Sci. Rep.* **9**, 9354 (2019). — <https://github.com/Boyle-Lab/Blacklist>
- Heinz, S. *et al.* Simple combinations of lineage-determining transcription factors prime cis-regulatory elements required for macrophage and B cell identities. *Mol. Cell* **38**, 576–589 (2010). (HOMER) — <http://homer.ucsd.edu/homer/>
- GENCODE annotation — <https://www.gencodegenes.org>; UCSC chrom.sizes and genome FASTA — <https://hgdownload.soe.ucsc.edu/downloads.html>

---

[← Index](../README.md) · [Next: Data transfer →](01_data_transfer.md)
