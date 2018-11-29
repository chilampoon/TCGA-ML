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

3. Select feature set 2, 4, 5 to test the performance in the 
    - whole dataset excluding training data (all are non-metastatic samples)
    - every single cancer type 
    

#### Results:

The performance of these 6 feature sets are not quite well, so I will try to do feature engineering and then try other ML algorithms like random forest which can take tissues into consideration.



# To Do:

1. ~~Redo the preprocessing and differential expression analyses using the `~ Metastasis + Type` design matrix~~ (done)
2.0 Try to utilize all `y=0` samples (rotation?)
2.1 Add tissue type dummy variables into glmnet matrix
2. Implement *glmnet* to the undersampled whole dataset and every individual tumor type **using log2(TPM + 1)**
3. Plot the performances
4. GOseq analysis using esIDs-GOterms `gene2Cat`
5. Plot the TPM desity/distribution plots
6. Sum up the TPM values in GOterms, then redo the *glmnet* and apply to *random forest* using new features
7. Try *xgboost*, *SVM* and other algorithms and other feature engineering methods!




