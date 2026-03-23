import sys
from subprocess import call

BED_REGIONS = "ctd_high.bed ctd_medium.bed ctd_low.bed ctd_zero.bed"
BED_REGIONS_HIGH = "ctd_high.bed"
SPP_PATH = "/data/budkina/himorna_chipseq_data/new_processing/phantompeakqualtools/run_spp.R"
BLACKLIST = "/data/budkina/himorna_chipseq_data/new_processing/hg38-blacklist.v3.bed"

C_IDS = ["C-9_", "C-10", "C-19", "C-20"]
I_IDS = ["I-9_", "I-10", "I-19", "I-20"]

def run(cmd, name):
    if call(cmd, shell=True) != 0:
        sys.exit(f"{name} failed")

for C_ID, I_ID in zip(C_IDS, I_IDS):
    run(f"java --add-modules=jdk.incubator.vector -Xmx8G -jar /omnipeak/build/libs/omnipeak-1.3.build.jar analyze"
        f" -t {C_ID}.sorted_for_index.bam -c {I_ID}.sorted_for_index.bam --cs hg38.chrom.sizes -p {C_ID}.peak",
        "omnipeak")

    run(f"Rscript {SPP_PATH} -rf -c={C_ID}.sorted_for_index.bam"
        f" -savp={C_ID}.spp.pdf -savd={C_ID}.spp.Rdata -out={C_ID}.spp.out",
        "phantomtools")

run(f"computeMatrix scale-regions --metagene -m 2000 -S 'ASO_G0267577_01 signal.bw' 'ASO_G0267577_03 signal.bw' 'ASO_G0267577_01 input.bw'  'ASO_G0267577_03 input.bw' "
f"-R {BED_REGIONS_HIGH} "
f"-a 2000 -b 2000 "
f"--blackListFileName {BLACKLIST} -o matrix1.mtx ",
"computeMatrix")

run(f"plotProfile   -m matrix1.mtx --plotTitle 'DNAAF3-AS1 KD, H3K36me3' -out Figure2A1.png --perGroup --plotHeight 12 --plotWidth 14 "
f"--colors '#E41A1C' '#377EB8' '#4DAF4A' '#984EA3'  --yMin 1 --yMax 9",
"plotProfile")

run(f"computeMatrix scale-regions --metagene -m 2000 -S 'Control1 signal.bw' 'Control2 signal.bw' 'Control1 input.bw'  'Control2 input.bw' "
f"-R {BED_REGIONS_HIGH} "
f"-a 2000 -b 2000 "
f"--blackListFileName {BLACKLIST} -o matrix2.mtx ",
"computeMatrix")

run(f"plotProfile   -m matrix2.mtx --plotTitle 'DNAAF3-AS1 Control, H3K36me3' -out Figure2A2.png --perGroup --plotHeight 12 --plotWidth 14 "
f"--colors '#E41A1C' '#377EB8' '#4DAF4A' '#984EA3'  --legendLocation upper-right --yMin 1 --yMax 9",
"plotProfile")
