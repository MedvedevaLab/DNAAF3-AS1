import sys
from subprocess import call

C_ID = "C-10"
I_ID = "I-10"
C_R1_FASTQ = "C-10_S42_R1_001.fastq.gz"
C_R2_FASTQ = "C-10_S42_R2_001.fastq.gz"
I_R1_FASTQ = "I-10_S14_R1_001.fastq.gz"
I_R2_FASTQ = "I-10_S14_R2_001.fastq.gz"
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
        f"fastp --in1 {fastq1} --in2 {fastq2} "
        f" --trim_poly_g --adapter_fasta adapters.fa "
        f"--length_required 36 -w {T} --poly_g_min_len 5 "
        f"--out1 {prefix}_1.fastq.gz --out2 {prefix}_2.fastq.gz", 
        f"fastp ({prefix})"
    )

    # bowtie2 -> SAM-to-BAM -> name sort -> fixmate -> pos sort -> markdup
    run(
        f"bowtie2 -p {T} --local"
        f" -x hg38_filtered"
        f" -1 {prefix}_1.fastq.gz"
        f" -2 {prefix}_2.fastq.gz"
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

preprocessing(C_R1_FASTQ, C_R2_FASTQ, C_ID)
preprocessing(I_R1_FASTQ, I_R2_FASTQ, I_ID)

