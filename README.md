# ATAC-seq_TCGA

# Cross validating specific noncoding mutation regions

The ATAC-seq (assay for transposase-accessible chromatin using sequencing) BigWig files for TCGA breast, lung, prostate, and bladder cancer samples were sourced from the Corces study (1) (TCGA https://gdc.cancer.gov/about-data/publications/ATACseq-AWG). These files, quantifying the normalized chromatin accessibility at 100 bp resolution, were prepared (1) based on the Hg38 genome assembly. 

Download the bigwig files using tar -zxvf file_name.tgz --strip-components 8
Convert bigwig files to wig files conda install bioconda::ucsc-bigwigtowig

Variant coordinates, originally mapped to Hg19, were converted to Hg38 using the online Lift Genome Annotation tool from the UCSC Genome Browser to match the ATAC data. 

For breast cancer, a total of 80 samples were analysed, including 18 basal subtype samples and 60 non-basal subtype samples as defined by PAM50 (prediction analysis of microarray 50) classification from the TCGA breast cancer study (2). 

RNA-seq expression data for 10 matched basal subtype samples and 40 matched non-basal subtype samples was obtained from the International Cancer Genome Consortium (ICGC) to further compare the expression between samples with high and low ATAC-seq signal. 

For lung cancers, the two cohorts annotated in TCGA were used as subtypes: LUAD (Lung adenocarcinoma, 19 samples) and LUSC (Lung squamous cell carcinoma, 11 samples). Prostate cancer and bladder cancer analyses included 52 and 20 samples, respectively.

# Genome-wide regulatory differences for ATAC-seq signal 

Calculate 





1.	Corces, M. R., Granja, J. M., Shams, S., Louie, B. H., Seoane, J. A., Zhou, W., Silva, T. C., Groeneveld, C., Wong, C. K., Cho, S. W., Satpathy, A. T., Mumbach, M. R., Hoadley, K. A., Robertson, A. G., Sheffield, N. C., Felau, I., Castro, M. A. A., Berman, B. P., Staudt, L. M., Zenklusen, J. C., … Chang, H. Y. The chromatin accessibility landscape of primary human cancers. Science 362, eaav1898 (2018). https://doi.org/10.1126/science.aav1898.
2.	The Cancer Genome Atlas Network. Comprehensive molecular portraits of human breast tumours.Nature 490,   61–70 (2012). https://doi.org/10.1038/nature11412
