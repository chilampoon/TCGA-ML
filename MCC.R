# MCC (Matthew's Correlation Coefficient)
# Correlation coefficient between the observed and predicted binary classifications
sum1 <- TP+FP
sum2 <- TP+FN
sum3 <- TN+FP
sum4 <- TN+FN
# as.double to avoid overflow error on large products
denom <- as.double(sum1)*sum2*sum3*sum4
if (any(sum1==0, sum2==0, sum3==0, sum4==0)) {
  denom <- 1
}
results[['MCC']] <- ((TP*TN)-(FP*FN)) / sqrt(denom)