import sys
from subprocess import call

SAMPLE_PAIRS = [
    ("C-10_S42", "I-10_S14"),
    ("C-9_S41",  "I-9_S13"),
    ("C-19_S51", "I-19_S23"),
    ("C-20_S52", "I-20_S24"),
]

THREADS = 24

def run(cmd, name):
    if call(cmd, shell=True) != 0:
        sys.exit(f"{name} failed")

run(f"seqkit grep -r -n -p '^chr([1-9]|1[0-9]|2[0-2]|[XY])$' hg38.fa"
    f" | bowtie2-build --threads {THREADS} - hg38_filtered",
    f"seqkit grep | bowtie2-build")

def preprocessing(fastq1, fastq2, prefix):
    T = THREADS
    run(
        f"trim_galore --paired --gzip -o trimmed/ {fastq1} {fastq2}",
        f"trim_galore ({prefix})"
    )
    # bowtie2 -> SAM-to-BAM -> name sort -> fixmate -> pos sort -> markdup
    run(
        f"bowtie2 -p {T}  --local  --very-sensitive-local -N 1 -L 16 -i S,1,0.25 -D 25 -R 4  --score-min G,10,6 --mp 3,1 "# --local
        f" -x /data/budkina/reference/filtered/hg38_filtered "
        f" -1 trimmed/{prefix}_R1_001_val_1.fq.gz "        # f" -1 "
        f" -2 trimmed/{prefix}_R2_001_val_2.fq.gz "        # f" -2 "
        f" 2> {prefix}.align.report"
        f" | samtools view -@ {T} -bS -"
        f" | samtools fixmate -@ {T} -m - -"
        f" | samtools sort -@ {T} -"
        f" | samtools markdup -@ {T} -r -s - {prefix}.dedup.bam"
        f" 2> {prefix}.dedup.report",
        f"preprocessing pipeline ({prefix})"
    )
    run(
        f"samtools sort -@ {T} {prefix}.dedup.bam -o {prefix}.sorted_for_index.bam"
        f" && samtools index -@ {T} {prefix}.sorted_for_index.bam"
        f" && bamCoverage -b {prefix}.sorted_for_index.bam --normalizeUsing RPKM --extendReads -o {prefix}.bw",
        f"sort-index-coverage ({prefix})"
    )

for C_ID, I_ID in SAMPLE_PAIRS:
    C_R1_FASTQ = f"{C_ID}_R1_001.fastq.gz"
    C_R2_FASTQ = f"{C_ID}_R2_001.fastq.gz"
    I_R1_FASTQ = f"{I_ID}_R1_001.fastq.gz"
    I_R2_FASTQ = f"{I_ID}_R2_001.fastq.gz"
    preprocessing(C_R1_FASTQ, C_R2_FASTQ, C_ID)
    preprocessing(I_R1_FASTQ, I_R2_FASTQ, I_ID)

# run(
#     f"""samtools view -h -f 2 {prefix}.dedup.bam | awk 'substr($1,1,1)=="@" || ($9>=300 || $9<=-300)'| samtools view -b -o {prefix}.frag300.bam - """,
#     f"filter pipeline ({prefix})"
# )