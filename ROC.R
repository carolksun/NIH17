library(ROCR)
library(nnet)

dfadj <- read.csv("test.csv", header = T)
dfun <- read.csv("testun.csv", header = T)

dfun <- read.csv("Adjusted_median_all_unadjusted_deleted.csv", header = T)
dfadj <- read.csv("Adjusted_median_all_adjusted_deleted.csv", header = T)

dfun <- read.csv("Adjusted_median_8_unadjusted.csv", header = T)
dfadj <- read.csv("Adjusted_median_8_adjusted.csv", header = T)

dfun <- read.csv("Adjusted_median_5_unadjusted.csv", header = T)
dfadj <- read.csv("Adjusted_median_5_adjusted.csv", header = T)

dfun <- read.csv("Ratios_unadjusted.csv", header = T)
dfadj <- read.csv("Ratios_adjusted.csv", header = T)

#create model for prediciton
mymodel <- multinom(Diagnosis~.,data=dfun)
pred <- predict(mymodel, dfun, type = "prob")

prednew <- prediction(pred, dfun$Diagnosis)
rocun <- performance(prednew,"tpr","fpr")

#Area Under Curve (AUC)
aucun <- performance(prednew,"auc")
aucun <- unlist(slot(performance(prednew, "auc"), "y.values"))

#create model for prediciton
mymodeladj <- multinom(Diagnosis~.,data=dfadj)
predadj <- predict(mymodeladj, dfadj, type = "prob")

prednewadj <- prediction(predadj, dfadj$Diagnosis)
rocadj <- performance(prednewadj,"tpr","fpr")

#Area Under Curve (AUC)
aucadj <- performance(prednewadj,"auc")
aucadj <- unlist(slot(performance(prednewadj, "auc"), "y.values"))

png(filename="5Proteins.png")
par(family="serif")
plot(rocun, col= "black",lwd=3)
par(new=TRUE)
plot(rocadj, col= "red",main="SD=0.8",lwd=3, lty=2)
abline(a=0,b=1, col = "grey71")
legend("bottomright", legend=c(paste("Adjusted (AUC = ", round(aucadj,3),")", sep = ""), paste("Unadjusted (AUC = ", round(aucun,3),")", sep = "")), col=c("red", "black"),lty=2,lwd=3)
dev.off()