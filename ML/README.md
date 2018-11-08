# Machine Learning Section


## Try the glmnet cross-validation first:

Since there is another extra TCGA dataset, this 6k dataset can be used as a whole training part. 

1. Do differential expression analysis using filtered all gene set for
    - Design 1: ~ Metastasis
    - Design 2: ~ Metastasis + Type
    
  After DEA, there are several available DEG sets:
  
    1. DEGs with q<0.1 & FC > 1.5 in design 1 only (2656)

    2. DEGs with q<0.1 & FC > 2 in design 1 only (556)

    3. DEGs with q<0.1 in design 2 only (1256)

    4. DEGs with q<0.1 & FC > 1.2 in design 2 only (380)

    5. Overlapped subset between design 1 & 2 (q<0.1 only) (993)

    6. Overlapped subset between design 1 & 2 (q<0.1 & FC > 1.2) (203)

 Get different gene esID list as features.

2. Use the feature sets to do glmnet first, finding the best combination data
    - 5-fold cross validation
    - reduce the gene numbers or not

##### To be completed:
1. Select feature set 2, 4, 5 to test the performance in the 
    - whole dataset excluding training data (all are non-metastatic samples)
    - every single cancer type 

2. Select top genes with high coefficients in lasso regression in results from 1, then do ROC tests

#### Results:

The performance of these 6 feature sets are not quite well, so I will try to do feature engineering and then try other ML algorithms like random forest which can take tissues into consideration.

---

## Feature engineering

### GO terms as features

Since the enrichment analysis is often followed by differential expression analyses, we can use the top GO terms list, sum up the count values of gene in those gene lists.


