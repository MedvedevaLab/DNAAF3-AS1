import sys
from subprocess import call

BED_REGIONS_HIGH = "ctd_high.bed"
SPP_PATH = "/data/budkina/himorna_chipseq_data/new_processing/phantompeakqualtools/run_spp.R"
BLACKLIST = "/data/budkina/himorna_chipseq_data/new_processing/hg38-blacklist.v3.bed"

SAMPLE_PAIRS = [
    ("C-10_S42", "I-10_S14"),
    ("C-9_S41",  "I-9_S13")
    ("C-19_S51", "I-19_S23"),
    ("C-20_S52", "I-20_S24"),
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

# cp 'C-9_S41.bw' 'ASO_G0267577_01 signal.bw'
# cp 'C-10_S42.bw' 'ASO_G0267577_03 signal.bw'
# cp 'I-9_S13.bw' 'ASO_G0267577_01 input.bw'
# cp 'I-10_S14.bw' 'ASO_G0267577_03 input.bw'
# cp 'C-19_S51.bw' 'Control1 signal.bw'
# cp 'C-20_S52.bw' 'Control2 signal.bw'
# cp 'I-19_S23.bw' 'Control1 input.bw'
# cp 'I-20_S24.bw' 'Control2 input.bw'

run(f"computeMatrix scale-regions --metagene -m 2000 -S 'ASO_G0267577_01 signal.bw' 'ASO_G0267577_03 signal.bw' 'ASO_G0267577_01 input.bw'  'ASO_G0267577_03 input.bw' "
f"-R {BED_REGIONS_HIGH} "
f"-a 2000 -b 2000 "
f"--blackListFileName {BLACKLIST} -o matrix1.mtx ",
"computeMatrix")

run(f"plotProfile   -m matrix1.mtx --plotTitle 'DNAAF3-AS1 KD, H3K36me3' -out Figure2A1.png --perGroup --plotHeight 12 --plotWidth 14 "
f"--colors '#E41A1C' '#377EB8' '#4DAF4A' '#984EA3'  --yMin 1 --yMax 15",
"plotProfile")

run(f"computeMatrix scale-regions --metagene -m 2000 -S 'Control1 signal.bw' 'Control2 signal.bw' 'Control1 input.bw'  'Control2 input.bw' "
f"-R {BED_REGIONS_HIGH} "
f"-a 2000 -b 2000 "
f"--blackListFileName {BLACKLIST} -o matrix2.mtx ",
"computeMatrix")

run(f"plotProfile   -m matrix2.mtx --plotTitle 'DNAAF3-AS1 Control, H3K36me3' -out Figure2A2.png --perGroup --plotHeight 12 --plotWidth 14 "
f"--colors '#E41A1C' '#377EB8' '#4DAF4A' '#984EA3'  --legendLocation upper-right --yMin 1 --yMax 15",
"plotProfile")
