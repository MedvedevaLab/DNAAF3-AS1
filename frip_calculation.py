import sys
from subprocess import call, check_output

# Fraction of reads in peaks (FRiP) for each ChIP experiment.
SAMPLES = [
    ("C-9_S41",  "ASO_01"),
    ("C-10_S42", "ASO_03"),
    ("C-19_S51", "NCA_01"),
    ("C-20_S52", "NCA_02"),
]

BLACKLIST = "GRCh38-blacklist.v3.bed"

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
    run(f"cut -f1-3 {peaks} "
        f"| bedtools intersect -v -a - -b {BLACKLIST} "
        f"| sort -k1,1 -k2,2n "
        f"| bedtools merge -i - > {peaks_bed}",
        f"prepare peaks ({c_id})")
    filt = "-F 0x404 -q 1"
    total   = count(f"samtools view -c {filt} {bam}", f"total reads ({c_id})")
    in_peak = count(f"samtools view -c {filt} -L {peaks_bed} {bam}",
                    f"reads in peaks ({c_id})")
    frip = in_peak / total if total else 0.0
    rows.append((label, total, in_peak, frip))
    print(f"{label}\t{c_id}\ttotal={total}\tin_peak={in_peak}\tFRIP={frip:.4f}")

with open("frip_stat.csv", "w") as fh:
    fh.write("Sample\tTotal_reads\tReads_in_peaks\tFRIP\n")
    for label, total, in_peak, frip in rows:
        fh.write(f"{label}\t{total}\t{in_peak}\t{frip:.6f}\n")
