#!/bin/bash
set -euo pipefail

BASE_DIR="${BASE_DIR:-$HOME/rnaseq_project}"
FASTQ_DIR="${FASTQ_DIR:-$BASE_DIR/fastq}"
GENOME_DIR="${GENOME_DIR:-$BASE_DIR/genome}"
OUT_DIR="${OUT_DIR:-$BASE_DIR/dual_rnaseq_out}"

HOST_GENOME="${HOST_GENOME:-$GENOME_DIR/host.fna}"
PATHOGEN_GENOME="${PATHOGEN_GENOME:-$GENOME_DIR/pathogen.fna}"
HOST_IDX="${HOST_IDX:-$GENOME_DIR/host_idx}"          # hisat2 index prefix

TRIM_JAR="${TRIM_JAR:-/opt/Trimmomatic-0.39/trimmomatic-0.39.jar}"
ADAPTERS="${ADAPTERS:-/opt/Trimmomatic-0.39/adapters/TruSeq3-PE-2.fa}"

THREADS="${THREADS:-8}"
ORDER="${1:-${ORDER:-host_first}}"   # host_first | pathogen_first

QC_DIR="$OUT_DIR/qc"
TRIM_DIR="$OUT_DIR/trimmed"
UNMAPPED_DIR="$OUT_DIR/unmapped"
HOST_DIR="$OUT_DIR/host_aligned"
PATHOGEN_DIR="$OUT_DIR/pathogen_aligned"
LOG_DIR="$OUT_DIR/logs"
SUMMARY="$OUT_DIR/mapping_summary.tsv"

case "$ORDER" in
    host_first)     FIRST_ORG="host";     SECOND_ORG="pathogen" ;;
    pathogen_first) FIRST_ORG="pathogen"; SECOND_ORG="host" ;;
    *) echo "ERROR: ORDER must be 'host_first' or 'pathogen_first' (got: $ORDER)" >&2; exit 1 ;;
esac

mkdir -p "$QC_DIR" "$TRIM_DIR" "$UNMAPPED_DIR" "$HOST_DIR" "$PATHOGEN_DIR" "$LOG_DIR"

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }

for tool in fastqc java bwa hisat2 hisat2-build samtools bedtools; do
    command -v "$tool" >/dev/null 2>&1 || { echo "ERROR: required tool '$tool' not found in PATH" >&2; exit 1; }
done
[[ -f "$TRIM_JAR" ]]  || { echo "ERROR: Trimmomatic jar not found at $TRIM_JAR" >&2; exit 1; }
[[ -f "$ADAPTERS" ]]  || { echo "ERROR: adapter file not found at $ADAPTERS" >&2; exit 1; }
[[ -f "$HOST_GENOME" ]]     || { echo "ERROR: host genome not found at $HOST_GENOME" >&2; exit 1; }
[[ -f "$PATHOGEN_GENOME" ]] || { echo "ERROR: pathogen genome not found at $PATHOGEN_GENOME" >&2; exit 1; }

log "=== Dual RNA-seq pipeline starting (order: $ORDER) ==="

log "Step 1: FastQC on raw reads"
fastqc -o "$QC_DIR" -t "$THREADS" "$FASTQ_DIR"/*.fastq.gz

log "Step 2: Trimming reads with Trimmomatic"
for R1 in "$FASTQ_DIR"/*_1.fastq.gz; do
    R2="${R1/_1.fastq.gz/_2.fastq.gz}"
    SAMPLE="$(basename "$R1" _1.fastq.gz)"

    if [[ -s "$TRIM_DIR/${SAMPLE}_1.fq.gz" && -s "$TRIM_DIR/${SAMPLE}_2.fq.gz" ]]; then
        log "  $SAMPLE: trimmed files already exist, skipping"
        continue
    fi

    java -jar "$TRIM_JAR" PE -threads "$THREADS" \
        "$R1" "$R2" \
        "$TRIM_DIR/${SAMPLE}_1.fq.gz" "$TRIM_DIR/${SAMPLE}_1.unpaired.fq.gz" \
        "$TRIM_DIR/${SAMPLE}_2.fq.gz" "$TRIM_DIR/${SAMPLE}_2.unpaired.fq.gz" \
        ILLUMINACLIP:"$ADAPTERS":2:30:10 LEADING:3 TRAILING:3 \
        SLIDINGWINDOW:4:15 MINLEN:36 \
        2> "$LOG_DIR/${SAMPLE}_trimmomatic.log"
done

log "Step 3: Building genome indices (skipped if already present)"
if [[ ! -f "${HOST_IDX}.1.ht2" ]]; then
    log "  Building hisat2 index for host genome"
    hisat2-build -p "$THREADS" "$HOST_GENOME" "$HOST_IDX" > "$LOG_DIR/hisat2_build.log" 2>&1
else
    log "  Host hisat2 index already exists, skipping"
fi

if [[ ! -f "${PATHOGEN_GENOME}.bwt" ]]; then
    log "  Building bwa index for pathogen genome"
    bwa index "$PATHOGEN_GENOME" > "$LOG_DIR/bwa_index.log" 2>&1
else
    log "  Pathogen bwa index already exists, skipping"
fi

align_host() {
    local R1="$1" R2="$2" SAMPLE="$3" OUTDIR="$4"
    hisat2 -x "$HOST_IDX" -1 "$R1" -2 "$R2" -p "$THREADS" \
        --summary-file "$LOG_DIR/${SAMPLE}_host_hisat2.log" \
        | samtools sort -@ "$THREADS" -o "$OUTDIR/${SAMPLE}.bam" -
    samtools index "$OUTDIR/${SAMPLE}.bam"
}

align_pathogen() {
    local R1="$1" R2="$2" SAMPLE="$3" OUTDIR="$4"
    bwa mem -t "$THREADS" "$PATHOGEN_GENOME" "$R1" "$R2" \
        2> "$LOG_DIR/${SAMPLE}_pathogen_bwa.log" \
        | samtools sort -@ "$THREADS" -o "$OUTDIR/${SAMPLE}.bam" -
    samtools index "$OUTDIR/${SAMPLE}.bam"
}

extract_unmapped() {
    local BAM="$1" SAMPLE="$2" PREFIX="$3"
    samtools view -b -f 12 -F 256 -@ "$THREADS" "$BAM" > "$UNMAPPED_DIR/${SAMPLE}_${PREFIX}_unmapped.bam"
    bedtools bamtofastq \
        -i "$UNMAPPED_DIR/${SAMPLE}_${PREFIX}_unmapped.bam" \
        -fq  "$UNMAPPED_DIR/${SAMPLE}_${PREFIX}_1.fq" \
        -fq2 "$UNMAPPED_DIR/${SAMPLE}_${PREFIX}_2.fq"
}

log "Step 4: Aligning samples ($FIRST_ORG first, then $SECOND_ORG)"
echo -e "sample\t${FIRST_ORG}_mapped_reads\t${SECOND_ORG}_mapped_reads" > "$SUMMARY"

for R1 in "$TRIM_DIR"/*_1.fq.gz; do
    R2="${R1/_1.fq.gz/_2.fq.gz}"
    SAMPLE="$(basename "$R1" _1.fq.gz)"
    log "  Sample: $SAMPLE"

    if [[ "$FIRST_ORG" == "host" ]]; then
        FIRST_DIR="$HOST_DIR"
        align_host "$R1" "$R2" "$SAMPLE" "$FIRST_DIR"
    else
        FIRST_DIR="$PATHOGEN_DIR"
        align_pathogen "$R1" "$R2" "$SAMPLE" "$FIRST_DIR"
    fi

    extract_unmapped "$FIRST_DIR/${SAMPLE}.bam" "$SAMPLE" "$FIRST_ORG"

    U1="$UNMAPPED_DIR/${SAMPLE}_${FIRST_ORG}_1.fq"
    U2="$UNMAPPED_DIR/${SAMPLE}_${FIRST_ORG}_2.fq"

    if [[ "$SECOND_ORG" == "host" ]]; then
        SECOND_DIR="$HOST_DIR"
        align_host "$U1" "$U2" "$SAMPLE" "$SECOND_DIR"
    else
        SECOND_DIR="$PATHOGEN_DIR"
        align_pathogen "$U1" "$U2" "$SAMPLE" "$SECOND_DIR"
    fi

    FIRST_MAPPED=$(samtools view -c -F 260 "$FIRST_DIR/${SAMPLE}.bam")
    SECOND_MAPPED=$(samtools view -c -F 260 "$SECOND_DIR/${SAMPLE}.bam")
    echo -e "${SAMPLE}\t${FIRST_MAPPED}\t${SECOND_MAPPED}" >> "$SUMMARY"
done

log "=== Dual RNA-seq pipeline finished ==="
log "Host BAMs:      $HOST_DIR"
log "Pathogen BAMs:  $PATHOGEN_DIR"
log "Mapping summary: $SUMMARY"
