# Part 1 · Getting data in and out

[← Setup](00_setup.md) · [Index](../README.md) · [Next: QC & trimming →](02_qc_trimming.md)

> ⚠️ **Never write tokens or passwords into scripts.** Keep them in a file only you can read (`chmod 600 ~/.secrets`) and `source` it, or use each tool's own login command. All credentials below are placeholders.

---

## 1.1 Illumina BaseSpace (`bs` CLI)

```bash
# one-time install
mkdir -p $HOME/bin
wget "https://launch.basespace.illumina.com/CLI/latest/amd64-linux/bs" -O $HOME/bin/bs
chmod u+x $HOME/bin/bs

# one-time login; stores a config named "mywork" in ~/.basespace/
bs auth -c mywork

bs -c mywork list projects
bs -c mywork list runs
```

Download either **demultiplexed FASTQs (project)** or the **raw run folder (BCL)**:

```bash
OUT=/path/to/Project/Data

# FASTQs from a project
bsub -n 2 -e log_bs.e -o log_bs.o \
  bs -c mywork download project --id <PROJECT_ID> -o ${OUT}/Fastq

# whole run folder (for bcl2fastq)
bsub -n 2 -e log_bs.e -o log_bs.o \
  bs -c mywork download run --id <RUN_ID> -o ${OUT}/0_BaseSpace
```

BaseSpace splits FASTQs by lane; combine lanes per read (gzip files can be concatenated directly):

```bash
for s in $(ls Fastq/*_L001_R1_001.fastq.gz | xargs -n1 basename | sed 's/_L001_R1_001.fastq.gz//'); do
  cat Fastq/${s}_L00*_R1_001.fastq.gz > Fastq/${s}_R1.fastq.gz
  cat Fastq/${s}_L00*_R2_001.fastq.gz > Fastq/${s}_R2.fastq.gz
done
```

## 1.2 Demultiplex BCL → FASTQ (`bcl2fastq`)

```bash
module load bcl2fastq2/v2.20.0.422
DATADIR=/path/to/Project/Data

bsub -n 8 -M 32000 -e log_bcl2fq.e -o log_bcl2fq.o \
  bcl2fastq \
    --runfolder-dir ${DATADIR}/0_BaseSpace \
    --output-dir    ${DATADIR}/Fastq \
    --sample-sheet  ${DATADIR}/SampleSheet.csv \
    --no-lane-splitting -p 8
```

`SampleSheet.csv` holds the index sequences (and optional adapter sequences for trimming). `--no-lane-splitting` gives one FASTQ per sample and read, so the lane-merging step is not needed.

## 1.3 Public data from SRA / GEO

```bash
# sra-tools ≥ 2.10: prefetch + fasterq-dump is faster than fastq-dump
for srr in SRR4205310 SRR4205311 SRR4205315; do
  bsub -n 8 -e log_${srr}.e -o log_${srr}.o \
    "prefetch ${srr} && fasterq-dump --split-3 -e 8 ${srr} && gzip ${srr}*.fastq"
done
```

`--split-3` writes `_1.fastq`/`_2.fastq` for paired-end runs and a single file for unpaired reads.

### ENCODE metadata (R)

```r
library(ENCODExplorer)   # may be missing from recent Bioconductor releases; the ENCODE portal REST API is an alternative
data(encode_df, package = "ENCODExplorer")

res <- queryEncode(df = encode_df, organism = "Homo sapiens",
                   biosample_name = "MCF-7", file_format = "fastq", fixed = FALSE)
fz  <- fuzzySearch(searchTerm = c("Bach1"), filterVector = "target", database = encode_df)
design <- createDesign(input = fz, df = encode_df, fileFormat = "fastq")
downloadEncode(res, df = encode_df, format = "fastq")
```

---

## 1.4 Sync with Box (`lftp`)

Set an external FTP password in Box first (*Account Settings → Box FTP*; this password is **not** your institutional single-sign-on password).

```bash
# HPC  → Box   (mirror -R = reverse/upload)
lftp -p 990 -e "open <box_login_email>@ftp.box.com; mirror -R /path/to/Project/Results Box_folder/Results; exit"

# Box  → HPC
lftp -p 990 -e "open <box_login_email>@ftp.box.com; mirror Box_folder/RawData /path/to/Project/RawData; exit"
```

## 1.5 Sync with Dropbox

```bash
# Option A: dbxcli (lab module)
module use.own /path/to/lab/privatemodules
module load dropbox_cli            # or install dbxcli from GitHub releases
dbxcli account                                  # first run: authorise
dbxcli get /Lab/Project/Wig/sample.bw ./sample.bw
dbxcli put ./pool_CTCF_Arm.bw /Lab/Project/Wig/pool_CTCF_Arm.bw

# Option B: Dropbox-Uploader (folders)
git clone https://github.com/andreafabrizi/Dropbox-Uploader.git
chmod +x Dropbox-Uploader/dropbox_uploader.sh
./Dropbox-Uploader/dropbox_uploader.sh            # first run: set up app token
./Dropbox-Uploader/dropbox_uploader.sh upload  ./Results  /Lab/Project/Results
./Dropbox-Uploader/dropbox_uploader.sh download /Lab/Project/Raw ./Raw
```

## 1.6 Submit to GEO (FTP)

GEO gives every submitter a personal upload folder and password on the *Transfer files* page. Put the metadata spreadsheet, raw FASTQs and processed files into one local folder, then:

```bash
# credentials from https://www.ncbi.nlm.nih.gov/geo/info/submissionftp.html
source ~/.secrets      # defines GEO_FTP_PASS and GEO_UPLOAD_DIR (e.g. uploads/<email>_<code>)

ncftpput -F -R -z -u geoftp -p "${GEO_FTP_PASS}" \
  ftp-private.ncbi.nlm.nih.gov "./${GEO_UPLOAD_DIR}" /path/to/Project/Data_ToGEO
```

`-R` uploads recursively and `-z` resumes interrupted transfers. Notify GEO through the submission form when the upload has finished.

---

## 1.7 Down-sample a BAM (depth matching)

```bash
id=ChIP_sample
# keep ~10 M random reads (single-end)
cat <(samtools view -H ${id}.bam) \
    <(samtools view ${id}.bam | shuf -n 10000000 --random-source=<(yes)) \
  | samtools sort -@ 4 -T ${id}_tmp -o ${id}_rand10M.bam -

# paired-end: subsample by fraction so mates stay together (seed 42, fraction 0.25)
samtools view -b -s 42.25 ${id}.bam > ${id}_sub25.bam
```

> `shuf` samples individual lines, so it splits read pairs. Use `samtools view -s SEED.FRACTION` for paired-end data. To work out the fraction, divide the target read count by `samtools view -c ${id}.bam`.

---

## References and resources

- Illumina BaseSpace CLI — <https://developer.basespace.illumina.com/docs/content/documentation/cli/cli-overview>
- Illumina bcl2fastq2 conversion software guide (and its successor, BCL Convert) — <https://support.illumina.com>
- SRA Toolkit (`prefetch`, `fasterq-dump`) — <https://github.com/ncbi/sra-tools>
- NCBI GEO submission guide and FTP upload — <https://www.ncbi.nlm.nih.gov/geo/info/submission.html>
- ENCODE portal and REST API — <https://www.encodeproject.org/help/rest-api/>
- Danecek, P. *et al.* Twelve years of SAMtools and BCFtools. *GigaScience* **10**, giab008 (2021). — <http://www.htslib.org>

---

[← Setup](00_setup.md) · [Index](../README.md) · [Next: QC & trimming →](02_qc_trimming.md)
