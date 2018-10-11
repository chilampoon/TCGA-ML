# Differential Expression Analysis

Before applying the machine learning algorithms, here we first did a differential exression analysis for feature selection using this TCGA RNA-seq dataset.

|Data               | File             |
|-------------------|------------------|
|Tximport object    |SalmonTXI_TCGA.rds|
|Patient information| tcga-cdr.csv     |
|Sample information | sample.info.tsv  |

| Scripts         |Step |
|-----------------|-----|
|data.process.Rmd |1,2,3|
|DE1.Rmd          |4-1  |
|allGenes.DE1.Rmd |4-2  |
|Single.DE.Rmd    |4-3  |
|GSEA.Rmd         |5    |


### Main steps

#### 1. Filt samples
- Discard patients with multiple metastatic sites;
- Keep primary tumor samples only.
 
#### 2. Split the dataset according to biotypes
- Use biomaRt package to extract protein-coding genes;
- Confirmed that the tximport object was generated using ENSEMBL v75 hg19 annnotation GTF file (tximport is to summarize abundance from transcript level to gene level).
 
#### 3. Generate gene-level counts-from-abundance
- Scale the count matrix for afterward limma-voom implementation

#### 4. Filter genes
Since different gene filtering ways cause different DE results, and some types don't have metastasis samples or just 1-5 samples. Also I set up the DE metastatic sample threshold is 10, therefore here using the `keep <- rowSums(cpm(cfa) > 1) > 5` where the 5 is 10/2.

#### 5. Limma-voom 
##### 5-1. For protein-coding genes pancancer subset
##### 5-2. For all genes pancancer set
##### 5-3. Implement limma-voom in cancer types with more than 10 metastasis samples for coding gene subsets
- Fitler out samples with large proportion of lowly-expressed genes (outliers);
- Remove lowly-expressed genes using functions in DESeq2 & edgeR;
- Normalize the counts using TMM method;
- Voom and limma:
  - design1: ~ Metastasis
  - design2: ~ Metastasis + type (adding the confounding factor cancer type)
- Draw boxplots with top several DEGs for each cancer type and heatmap using 30 DEGs (q-vals and logFC threshold see the file).

#### 6. Functional enrichment analysis
 - Conducted pre-ranked GSEA using gene lists ranked by the t-statistics from the results of DE analysis;
 - Evaluation: GOseq -- a bit more different terms

---
### Conclusions:

1. The results DEA of single type are strange, there is no DEGs in BRCA or BLCA etc. tissues that have enough sample numbers, so then we can look into what happens here -- check any bugs? draw the boxplot of some top genes ranked by logFC.

2. From the heatmaps, the DEGs (qval<0.05 & FC > 2) of design 1 in PCG subset and ALL gene subset can classify the metastatic and primary tumors better. each can try.

3. Do the classification for 
   - single cancer type 
   - pan cancer

One can use the cancer census genes to check if it's a cancer oncogene or suppressor genes

---
### To do:

1. Redo the single type DEA, adding PCA analysis.
2. Extract BRCA & HNSC subsets, first spilt into training and validation parts, then to redo the differential expression analysis **on training part**, then to try the machine learning pipeline!


