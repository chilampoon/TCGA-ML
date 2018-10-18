# Machine Learning Section

To do：

1. Firstly spilt out small datasets (those with more than 10 metastasis samples) -- done
2. Use the whole small dataset to do PCA and differential expression analysis using different formula, note the results.
- The PCA... I don't know how to scale the RNA-seq data and how to plot?
- If correct the gender and race, the q-vals are even higher in BRCA
- Though the protein coding genes are not differentially expressed in BRCA, the non-coding RNAs do. 

2-1. Do the all gene set for those 11 cancer types -> gold bar plot -- done

2-2. Choose some genes and plot some barplot in BRCA -- done

---
## Try the glmnet cross-validation first:

Since there is another extra TCGA dataset, this 6k dataset can be used as a whole training part. 

1. Do differential expression analysis using filtered all gene set:
    - Design 1: ~ Metastasis
    - Design 2: ~ Metastasis + Type
    
  After DEA, there are several available DEG sets:
    1) q<0.1 & FC > 1.5 in design 1 only
    2) q<0.1 & FC > 1.5 in design 2 only
    3) the overlapped subset between design 1 & 2 (q<0.1 & FC > 1.5 || q < 0.1 only)

Get different gene esID list as features.

2. Use the feature sets to do glmnet first, finding the best combination data
  - 5-fold cross validation
  - reduce the gene numbers or not
