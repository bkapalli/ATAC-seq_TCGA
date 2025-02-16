# ATAC-seq_TCGA

## Cross-Validating Specific Noncoding Mutation Regions

The ATAC-seq (Assay for Transposase-Accessible Chromatin using sequencing) BigWig files for TCGA breast, lung, prostate, and bladder cancer samples were sourced from the Corces study [(1)](https://gdc.cancer.gov/about-data/publications/ATACseq-AWG). These files quantify normalized chromatin accessibility at 100 bp resolution and were prepared using the Hg38 genome assembly.

To download the BigWig files:
```
tar -zxvf file_name.tgz --strip-components 8
```
To convert BigWig files to Wig files:
```
conda install bioconda::ucsc-bigwigtowig
```

Variant coordinates, originally mapped to Hg19, were converted to Hg38 using the UCSC Genome Browser's Lift Genome Annotation tool to match the ATAC-seq data.

For breast cancer, a total of 80 samples were analyzed, including 18 basal subtype samples and 60 non-basal subtype samples, classified using PAM50 (Prediction Analysis of Microarray 50) from the TCGA breast cancer study [(2)](https://doi.org/10.1038/nature11412).

RNA-seq expression data for 10 matched basal subtype samples and 40 matched non-basal subtype samples was obtained from the International Cancer Genome Consortium (ICGC). PCAWG was used to decrypt the barcodes, allowing comparison of gene expression between samples with high and low ATAC-seq signal.

For lung cancer, two TCGA cohorts were analyzed: LUAD (Lung Adenocarcinoma, 19 samples) and LUSC (Lung Squamous Cell Carcinoma, 11 samples). Prostate cancer and bladder cancer analyses included 52 and 20 samples, respectively.

### Script Descriptions
- `gene.py` → Converts Wig files to count matrix files for specific gene mutation regions.
- `mean_signal.py` → Generates mean ATAC-seq signal plots highlighting mutation regions.
- `expressionplot.py` → Creates scatter plots showing correlation between gene expression and ATAC-seq signals.
- `foldchangeplot.py` → Produces line plots showing ATAC-seq fold changes across different cancer subtypes.

## Genome-Wide Regulatory Differences in ATAC-seq Signal

Genomic signal data was processed based on ATAC-seq base pair resolution to ensure a score for each 100 bp region. The data was categorized by genomic region subtypes, including 

The genomic signal data was further classified by functional genomic subtypes:
- **CX** → Coding nonsynonymous
- **CR** → Splice variant
- **CS** → Coding synonymous
- **F5** → 5' flank
- **U5** → 5' UTR
- **I5** → 5' end intron
- **IX** → Other introns
- **U3** → 3' UTR
- **F3** → 3' flank
- **R** → Noncoding RNA
- **J** → Intergenic

Two output files were generated—one for mutation data and one for coverage.

### Script Descriptions
- `mean.py` → Calculates the mean ATAC-seq signal and maps it to the NCBI gene list [(NCBI RefSeq)](https://www.ncbi.nlm.nih.gov/refseq/).
- `merge.py` → Merges genomic signal and ATAC-seq data.
- `subtype_analysis_plot1.py` → Visualizes genomic region distributions relative to mean ATAC-seq signal. A curve above the diagonal (y = x) suggests mutation enrichment, while alignment with the diagonal indicates a random distribution.
- `subtype_analysis_plot2.py` → Examines the distribution of promoter, intronic, and coding regions relative to mean ATAC-seq signal.


## References
1. Corces, M. R., et al. *The chromatin accessibility landscape of primary human cancers.* Science 362, eaav1898 (2018). [https://doi.org/10.1126/science.aav1898](https://doi.org/10.1126/science.aav1898)
2. The Cancer Genome Atlas Network. *Comprehensive molecular portraits of human breast tumors.* Nature 490, 61–70 (2012). [https://doi.org/10.1038/nature11412](https://doi.org/10.1038/nature11412)
