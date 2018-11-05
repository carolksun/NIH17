simulationalleles <- runif(200,0,1)
alleles <- simulationalleles
allelefactor <- alleles
diagnosiscontrol <- matrix(0,nrow=100,ncol=1)
diagnosiscase <- matrix(1,nrow=100,ncol=1)
diagnosis <- rbind(diagnosiscontrol,diagnosiscase)

for(i in 1:200){
  if(simulationalleles[i]<0.49){
    allelefactor[i] <- "AA"
  }
  if(simulationalleles[i]>0.91){
    allelefactor[i] <- "aa"
  }
  if(simulationalleles[i]>0.49 && simulationalleles[i]<0.91){
    allelefactor[i] <- "Aa"
  }
}

for(i in 1:200){
  if(allelefactor[i]=="AA"){
    alleles[i] <- "0"
  }
  if(allelefactor[i]=="aa"){
    alleles[i] <- "2"
  }
  if(allelefactor[i]=="Aa"){
    alleles[i] <- "1"
  }
}

Protein <- alleles
allelenumeric <- as.numeric(alleles)
for(i in 1:100){
  Protein[i] <- 0+1*(allelenumeric[i]-1)+rnorm(1,0,1)
}

for(i in 101:200){
  Protein[i] <- 1+1*(allelenumeric[i]-1)+rnorm(1,0,1)
}
alleles <- as.factor(allelefactor)
df <- cbind(diagnosis,allelefactor,Protein)
colnames(df)<- c("Diagnosis", "SNP", "Protein")
df <- as.data.frame(df)
df$Protein <- as.numeric(as.character(df$Protein))
dfun <- df
dfun$Protein <- exp(dfun$Protein)


df.split<-split(df,df[[2]])
aa <- df.split[[1]][,c("Diagnosis", "SNP", "Protein")]
Aa <- df.split[[2]][,c("Diagnosis", "SNP", "Protein")]
AA <- df.split[[3]][,c("Diagnosis", "SNP", "Protein")]

rsafake <- aa
rsbfake <- Aa
rscfake <- AA
rsafake$Protein <- rsafake$Protein + (median(rsbfake$Protein)-median(rsafake$Protein))
rscfake$Protein <- rscfake$Protein + (median(rsbfake$Protein)-median(rscfake$Protein))
logrs <- rbind(rsafake,rsbfake,rscfake)

dat <- logrs[sample(nrow(logrs), 200), ]
rn <- rownames(dat)
dfadj <- dat[order(as.numeric(rn)), ]

dfadj$Protein <- exp(dfadj$Protein)

