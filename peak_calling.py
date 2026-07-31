import sys
from subprocess import call

BED_REGIONS_HIGH = "ctd_high.bed"
SPP_PATH = "/data/budkina/himorna_chipseq_data/new_processing/phantompeakqualtools/run_spp.R"
BLACKLIST = "/data/budkina/himorna_chipseq_data/new_processing/hg38-blacklist.v3.bed"

SAMPLE_PAIRS = [
    ("C-10_S42", "I-10_S14"),
    ("C-9_S41",  "I-9_S13"),
    ("C-19_S51", "I-19_S23"),
    ("C-20_S52", "I-20_S24")
]

def run(cmd, name):
    if call(cmd, shell=True) != 0:
        sys.exit(f"{name} failed")

for C_ID, I_ID in SAMPLE_PAIRS:
    run(f"macs3 callpeak --broad -t {C_ID}.sorted_for_index.bam -c {I_ID}.sorted_for_index.bam"
        f" -f BAM -g hs --broad-cutoff 0.1 -n {C_ID} --outdir {C_ID}",
        "macs3 callpeak broad")
    run(f"Rscript {SPP_PATH} -rf -c={C_ID}.sorted_for_index.bam"
        f" -savp={C_ID}.spp.pdf -savd={C_ID}.spp.Rdata -out={C_ID}.spp.out",
        "phantomtools")

cp '/data/budkina/himorna_chipseq_data/test/C-9_S41.bw' 'ASO_G0267577_01 signal.bw'
cp '/data/budkina/himorna_chipseq_data/test/C-10_S42.bw' 'ASO_G0267577_03 signal.bw'
cp '/data/budkina/himorna_chipseq_data/test/I-9_S13.bw' 'ASO_G0267577_01 input.bw'
cp '/data/budkina/himorna_chipseq_data/test/I-10_S14.bw' 'ASO_G0267577_03 input.bw'
cp '/data/budkina/himorna_chipseq_data/test/C-19_S51.bw' 'Control1 signal.bw'
cp '/data/budkina/himorna_chipseq_data/test/C-20_S52.bw' 'Control2 signal.bw'
cp '/data/budkina/himorna_chipseq_data/test/I-19_S23.bw' 'Control1 input.bw'
cp '/data/budkina/himorna_chipseq_data/test/I-20_S24.bw' 'Control2 input.bw'

run(f"computeMatrix scale-regions --metagene -m 2000 -S 'ASO_G0267577_01 signal.bw' 'ASO_G0267577_03 signal.bw' 'ASO_G0267577_01 input.bw'  'ASO_G0267577_03 input.bw' "
f"-R {BED_REGIONS_HIGH} "
f"-a 2000 -b 2000 "
f"--blackListFileName {BLACKLIST} -o matrix1.mtx ",
"computeMatrix")

run(f"plotProfile   -m matrix1.mtx --plotTitle 'DNAAF3-AS1 KD, H3K36me3' -out Figure2A1.png --perGroup --plotHeight 12 --plotWidth 14 "
f"--colors '#E41A1C' '#377EB8' '#4DAF4A' '#984EA3'  --yMin 1 --yMax 12",
"plotProfile")

run(f"computeMatrix scale-regions --metagene -m 2000 -S 'Control1 signal.bw' 'Control2 signal.bw' 'Control1 input.bw'  'Control2 input.bw' "
f"-R {BED_REGIONS_HIGH} "
f"-a 2000 -b 2000 "
f"--blackListFileName {BLACKLIST} -o matrix2.mtx ",
"computeMatrix")

run(f"plotProfile   -m matrix2.mtx --plotTitle 'DNAAF3-AS1 Control, H3K36me3' -out Figure2A2.png --perGroup --plotHeight 12 --plotWidth 14 "
f"--colors '#E41A1C' '#377EB8' '#4DAF4A' '#984EA3'  --legendLocation upper-right --yMin 1 --yMax 12",
"plotProfile")

# --- Clusterisation of metagene profiles for each bw coverage file ---
# For every bigwig, build a per-region metagene matrix and group the regions
# into KMEANS_CLUSTERS clusters by the shape of their profile. k-means is run
# once (via plotProfile --kmeans) to get the region-to-cluster assignment, then
# each cluster is plotted in its own plotProfile call so every metagene plot
# gets its own, independently auto-scaled y-axis.
KMEANS_CLUSTERS = 8

BW_FILES = [ ("ASO_G0267577_01 signal.bw", "ASO_G0267577_01 input.bw"),
    ("ASO_G0267577_03 signal.bw", "ASO_G0267577_03 input.bw"),
    ("Control1 signal.bw", "Control1 input.bw"),
    ("Control2 signal.bw", "Control2 input.bw")]

def split_bed_by_cluster(sorted_bed, prefix):
    """Split a deepTools --outFileSortedRegions BED into one BED6 per cluster.
    The last column of the sorted-regions file holds the deepTools group label
    (e.g. 'cluster_1'). Returns an ordered list of (cluster_label, bed_path).
    """
    clusters = {}
    with open(sorted_bed) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            cluster = fields[-1]
            clusters.setdefault(cluster, []).append("\t".join(fields[:6]))
    out = []
    for cluster in sorted(clusters):
        bed_path = f"{prefix}.{cluster}.bed"
        with open(bed_path, "w") as fh:
            fh.write("\n".join(clusters[cluster]) + "\n")
        out.append((cluster, bed_path))
    return out

for bw_s,bw_i  in BW_FILES:
    # a filesystem-safe prefix for the output files of this bw
    prefix = bw_s.replace(".bw", "").replace(" ", "_")
    run(f"computeMatrix scale-regions --metagene -m 2000 -S '{bw_s}' "
        f"-R {BED_REGIONS_HIGH} "
        f"-a 2000 -b 2000 "
        f"--blackListFileName {BLACKLIST} -o '{prefix}.cluster.mtx' ",
        f"computeMatrix ({bw_s})")
    # Run k-means once, only to obtain the region-to-cluster assignment.
    run(f"plotProfile -m '{prefix}.cluster.mtx' "
        f"--kmeans {KMEANS_CLUSTERS} "
        f"--outFileSortedRegions '{prefix}.clusters.bed' "
        f"-out '{prefix}.clusters.png'",
        f"plotProfile kmeans ('{bw_s}' '{bw_i}')")
    # Plot each cluster on its own, so each metagene plot has its own y-scale.
    for cluster, cluster_bed in split_bed_by_cluster(f"{prefix}.clusters.bed", prefix):
        run(f"computeMatrix scale-regions --metagene -m 2000 -S '{bw_s}' '{bw_i}' "
            f"-R '{cluster_bed}' "
            f"-a 2000 -b 2000 "
            f"--blackListFileName {BLACKLIST} -o '{prefix}.{cluster}.mtx' ",
            f"computeMatrix ('{bw_s}' '{bw_i}' {cluster})")
        run(f"plotProfile -m '{prefix}.{cluster}.mtx' "
            f"--plotTitle '{cluster}' --perGroup "
            f"-out '{prefix}.{cluster}.png'",
            f"plotProfile ('{bw_s}' '{bw_i}' {cluster})")
