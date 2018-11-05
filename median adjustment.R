#import data
df <- read.csv("TNCalldiagnosis.csv", header = TRUE)
dftraining <- read.csv("TNCtraining.csv", header = TRUE)

varfraction <- matrix(nrow=ncol(df)-2, ncol = 1)
samplefraction <- matrix(nrow=ncol(df)-2, ncol = 1)
for(i in 1:1){
  #split data for training
  df.split.training<-split(dftraining,dftraining[[i+1]])
  rsatraining <- df.split.training[[2]][,c("Diagnosis", colnames(df)[i+1], "Protein")]
  
  #deals with missing data and 2 groups of alleles
  if(nlevels(rsatraining[[2]]) == 4){ 
    rsatraining <- df.split.training[[2]][,c("Diagnosis", colnames(df)[i+1], "Protein")]
    rsbtraining <- df.split.training[[3]][,c("Diagnosis", colnames(df)[i+1], "Protein")]
    rsctraining <- df.split.training[[4]][,c("Diagnosis", colnames(df)[i+1], "Protein")]  
    rstraining <- rbind(rsatraining, rsbtraining, rsctraining)
  }
  if (nlevels(rsatraining[[2]]) == 2){ 
    rsatraining <- df.split.training[[1]][,c("Diagnosis", colnames(df)[i+1], "Protein")]
    rsbtraining <- df.split.training[[2]][,c("Diagnosis", colnames(df)[i+1], "Protein")]
    rstraining <- rbind(rsatraining, rsbtraining)
  }
  if (nlevels(rsatraining[[2]]) == 3){
    rsatraining <- df.split.training[[1]][,c("Diagnosis", colnames(df)[i+1], "Protein")]
    rsbtraining <- df.split.training[[2]][,c("Diagnosis", colnames(df)[i+1], "Protein")]
    rsctraining <- df.split.training[[3]][,c("Diagnosis", colnames(df)[i+1], "Protein")]
    rstraining <- rbind(rsatraining, rsbtraining, rsctraining)
  }
  
  #orders groups from major to minor
  if(nrow(rsatraining) < nrow(rsctraining)){
    fctstraining <- c(as.character(rsctraining[1,2]), as.character(rsbtraining[1,2]), as.character(rsatraining[1,2]))
  }else{
    fctstraining <- c(as.character(rsatraining[1,2]), as.character(rsbtraining[1,2]), as.character(rsctraining[1,2]))
  }
  
  #same process for all patients
  df.split<-split(df,df[[i+1]])
  rsa <- df.split[[2]][,c("Diagnosis", colnames(df)[i+1], "Protein")]
  if(nlevels(rsa[[2]]) == 4){ 
    rsa <- df.split[[2]][,c("Diagnosis", colnames(df)[i+1], "Protein")]
    rsb <- df.split[[3]][,c("Diagnosis", colnames(df)[i+1], "Protein")]
    rsc <- df.split[[4]][,c("Diagnosis", colnames(df)[i+1], "Protein")]  
    rs <- rbind(rsa, rsb, rsc)
  }
  if (nlevels(rsa[[2]]) == 2){ 
    rsa <- df.split[[1]][,c("Diagnosis", colnames(df)[i+1], "Protein")]
    rsb <- df.split[[2]][,c("Diagnosis", colnames(df)[i+1], "Protein")]
    rs <- rbind(rsa, rsb)
  }
  if (nlevels(rsa[[2]]) == 3){
    rsa <- df.split[[1]][,c("Diagnosis", colnames(df)[i+1], "Protein")]
    rsb <- df.split[[2]][,c("Diagnosis", colnames(df)[i+1], "Protein")]
    rsc <- df.split[[3]][,c("Diagnosis", colnames(df)[i+1], "Protein")]
    rs <- rbind(rsa, rsb, rsc)
  }
  e1 = as.numeric(rs$Protein)
  s1 = as.numeric(rs$rs16936516)
  lm1 = lm(e1 ~ s1)
  tidy(lm1)
  
  newdata <- subset(rs, Protein >= 604.85 & Protein <= 725.9, 
                    select=c(Protein))
  var(log(newdata$Protein))/var(log(rs$Protein))
  
  rsanew <- rsa
  rsbnew <- rsb
  rscnew <- rsc
  rsanew$Protein <- log(rsa$Protein) + (median(log(rsbtraining$Protein)) - median(log(rsatraining$Protein)))
  rsbnew$Protein <- log(rsb$Protein)
  rscnew$Protein <- log(rsc$Protein) + (median(log(rsbtraining$Protein)) - median(log(rsctraining$Protein)))
  rsnew <- rbind(rsanew,rsbnew,rscnew)
  
  #reorder data
  rn <- as.numeric(row.names(rsnew))
  new1<-data.frame(rsnew, rn)
  orderedadj = new1[do.call(order, new1["rn"]), ]
  
  rn <- as.numeric(row.names(rs))
  new1<-data.frame(rs, rn)
  orderedorg = new1[do.call(order, new1["rn"]), ]
  
  varfraction[i,1]<- var(log(rs$Protein))/var(rsnew$Protein)
  if(nrow(rsa)>nrow(rsc)){
    samplefraction[i,1]<- nrow(rsc)/nrow(rsa)
  }
  if(nrow(rsc)>nrow(rsa)){
    samplefraction[i,1]<- nrow(rsa)/nrow(rsc)
  }
}

median(rsa$Protein)
median(rsc$Protein)
#write table
#write.table(table, "C:\\Users\\sunck\\Documents\\plink_win32\\HDMS\\TLR4table.txt", sep = "\t")
write.table(orderedadj, "C:\\Users\\sunck\\Documents\\plink_win32\\HDMS\\MAPK14newadj.txt", sep = "\t")
test <- orderedorg
newdata <- subset(test, Protein >=606.3  & Protein <= 734.3, 
                  select=c(Diagnosis, colnames(df)[i+1],rn))