import sys
from subprocess import call, check_output

# Fraction of reads in peaks (FRiP) for each ChIP experiment.
#
# Uses the deduplicated, indexed BAMs produced by chipseq_processing.py and the
# broad peaks called by peak_calling.py (macs3 --broad, written to {C_ID}/).
# FRiP = reads overlapping called peaks / total mapped reads.

# ChIP samples and the experiment label each maps to (see the cp block in
# peak_calling.py: C-9 -> ASO_01, C-10 -> ASO_03, C-19/C-20 -> NCA controls).
SAMPLES = [
    ("C-9_S41",  "ASO_01"),
    ("C-10_S42", "ASO_03"),
    ("C-19_S51", "NCA_01"),
    ("C-20_S52", "NCA_02"),
]

BLACKLIST = "/data/budkina/himorna_chipseq_data/new_processing/hg38-blacklist.v3.bed"


def run(cmd, name):
    if call(cmd, shell=True) != 0:
        sys.exit(f"{name} failed")


def count(cmd, name):
    try:
        return int(check_output(cmd, shell=True).strip())
    except Exception as exc:                      # noqa: BLE001
        sys.exit(f"{name} failed: {exc}")


rows = []
for c_id, label in SAMPLES:
    bam       = f"{c_id}.sorted_for_index.bam"
    peaks     = f"{c_id}/{c_id}_peaks.broadPeak"
    peaks_bed = f"{c_id}.peaks.bed"
    # macs3 broadPeak is BED-like; take the first three columns, drop blacklist
    # regions, and merge so overlapping peaks are not double-counted.
    run(f"cut -f1-3 {peaks} "
        f"| bedtools intersect -v -a - -b {BLACKLIST} "
        f"| sort -k1,1 -k2,2n "
        f"| bedtools merge -i - > {peaks_bed}",
        f"prepare peaks ({c_id})")
    # count only primary, mapped, non-duplicate reads for both totals
    filt = "-F 0x404 -q 1"
    total   = count(f"samtools view -c {filt} {bam}", f"total reads ({c_id})")
    in_peak = count(f"samtools view -c {filt} -L {peaks_bed} {bam}",
                    f"reads in peaks ({c_id})")
    frip = in_peak / total if total else 0.0
    rows.append((label, total, in_peak, frip))
    print(f"{label}\t{c_id}\ttotal={total}\tin_peak={in_peak}\tFRIP={frip:.4f}")

with open("/data/budkina/ctd/frip_stat.csv", "w") as fh:
    fh.write("Sample\tTotal_reads\tReads_in_peaks\tFRIP\n")
    for label, total, in_peak, frip in rows:
        fh.write(f"{label}\t{total}\t{in_peak}\t{frip:.6f}\n")
