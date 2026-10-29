# DNAAF3-AS1

Analysis code accompanying **Knockdown of the long non-coding RNA DNAAF3-AS1 activates the p53 pathway and remodels the H3K36me3 landscape**.

**Authors:** Anna Budkina, Anatoliy Zubritskiy, Daria Marakulina, Andrey Buyan, Ivan Kulakovsky, Chi Wai Yip, Jude Aneke, Jay Shin, and Yulia A. Medvedeva.

## Overview

This study investigates the role of the long non-coding RNA **DNAAF3-AS1** in chromatin regulation. It integrates H3K36me3 ChIP-seq following DNAAF3-AS1 knockdown in human dermal fibroblasts with FANTOM6 CAGE, HiMoRNA associations, RNA-chromatin interaction data, and single-cell RNA-seq from adult human testis.

The analyses examine changes in H3K36me3 distribution, differential gene expression and promoter usage, p53 pathway activity, and DNAAF3-AS1-associated gene modules during spermatogenesis.

DNAAF3-AS1 is also referred to as **CTD-2587H24.5**, **AC010327.4**, and **ENSG00000267577** in the source datasets. The two knockdown reagents are **ASO_G0267577_01** (ASO_01) and **ASO_G0267577_03** (ASO_03).

## Repository contents

| Script | Purpose |
| --- | --- |
| `chipseq_processing.py` | Reference filtering, paired-end read trimming, alignment, duplicate removal, BAM indexing, and RPKM-normalized bigWig generation. |
| `peak_calling.py` | Broad peak calling, strand cross-correlation QC, metagene profiles, and k-means clustering of H3K36me3 profiles. |
| `frip_calculation.py` | Fraction of reads in peaks (FRiP), excluding blacklisted regions. |
| `DiffBind.R` | Differential H3K36me3 binding analysis and export of peak tables and BED files. |
| `genes_for_metagene.R` | Selection of gene sets from FANTOM6 control expression and export of BED files for metagene analysis. |
| `gsea.R` | GO biological process GSEA using FANTOM6 differential expression results. |
| `decoupler.R` | PROGENy pathway activity inference with decoupleR, including p53 activity across FANTOM6 perturbations. |
| `write_bed6_promoters.R` | Export of promoter tracks for MDM2, NR2F2, FASTK, and CFDP1 from precomputed differential promoter usage results. |
| `rna_chrom_contacts.R` | Integration and clustering of DNAAF3-AS1 RNA-chromatin contact regions. |
| `genometricorr.R` | Spatial association analysis of differential H3K36me3 peaks and RNA-chromatin contact regions. |
| `testis.R` | Seurat preprocessing, clustering, and manual cell-type annotation of adult human testis single-cell data. |
| `cell_profiler.R` | Visualization of cell morphology measurements. |
| `Figure1.R`–`Figure4.R` | Main analysis figures; mapping to manuscript numbering is given below. |
| `Supplementary_Figure1.R`–`Supplementary_Figure8.R` | Supplementary figure assembly. |

### Main figure mapping

Main figure script names retain earlier numbering. In the accompanying manuscript, Figure 1 is the study schematic.

| Script | Manuscript figure | Content |
| --- | --- | --- |
| `Figure1.R` | Figure 2 | DNAAF3-AS1 locus, tissue expression, differential expression, pathway activity, and GO enrichment. |
| `Figure2.R` | Figure 3 | Differential promoter usage and coding-region promoter changes. |
| `Figure3.R` | Figure 4 | HiMoRNA associations, differential H3K36me3 peaks, metagene profiles, and spatial association. |
| `Figure4.R` | Figure 5 | Testis cell populations, DNAAF3-AS1/DNAAF3 expression, and module-score correlations. |

Output filenames follow the script numbering. Supplementary scripts cover the corresponding supplementary figures, although panel order may differ from the manuscript.

## Data and reference resources

| Resource | Use |
| --- | --- |
| [BioProject PRJNA1551676](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1551676) | Raw H3K36me3 ChIP-seq reads; accession reported in the manuscript. |
| [FANTOM6](https://fantom.gsc.riken.jp/6/datafiles/Core_FANTOM6/) | Human dermal fibroblast CAGE expression and differential expression following ASO perturbations. |
| [GEO GSE112013](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE112013) | Adult human testis single-cell expression matrix. |
| [HiMoRNA](https://himorna.fbras.ru/) | Associations between DNAAF3-AS1 expression and histone modification peaks. The manuscript uses the database updated on 24 April 2024. |

### ChIP-seq sample identifiers

| Condition | ChIP prefix | Matched input prefix |
| --- | --- | --- |
| DNAAF3-AS1 ASO_01 | `C-9_S41` | `I-9_S13` |
| DNAAF3-AS1 ASO_03 | `C-10_S42` | `I-10_S14` |
| Non-targeting control 1 | `C-19_S51` | `I-19_S23` |
| Non-targeting control 2 | `C-20_S52` | `I-20_S24` |
