#import data
df <- read.csv("MMP7alldiagnosis.csv", header = TRUE)
dftraining <- read.csv("MMP7training.csv", header = TRUE)

#set table
table <- matrix(nrow = ncol(df)-2, ncol = 5)
colnames(table) = c("Variance (ln)", "Less than Original?", "Variances(ln) Statistically Different?", "Mann-Whitney Signed Rank Test (ln)","Significant P-value? (<0.01)")

for(i in 1:4){
#split data for training
df.split.training<-split(dftraining,dftraining[[i+1]])
rsatraining <- df.split.training[[1]][,c("Diagnosis", colnames(df)[i+1], "Protein")]

#deals with missing data and 2 groups of alleles
if(nlevels(rsatraining[[2]]) == 4){ 
  rsatraining <- df.split.training[[2]][,c("Diagnosis", colnames(df)[i+1], "Protein")]
  rsbtraining <- df.split.training[[3]][,c("Diagnosis", colnames(df)[i+1], "Protein")]
  rsctraining <- df.split.training[[4]][,c("Diagnosis", colnames(df)[i+1], "Protein")]  
  rstraining <- rbind(rsatraining, rsbtraining, rsctraining)
}
if (nlevels(rstraining[[2]]) == 2){ 
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
rsa <- df.split[[1]][,c("Diagnosis", colnames(df)[i+1], "Protein")]
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

#slope calculation
elog = log(as.numeric(rstraining$Protein))
slog = as.numeric(factor(rstraining[[2]],fctstraining))
logmlm = mblm(elog ~ slog, repeated = F)
logline = tidy(logmlm)

#actual adjustment
rsaadj <- rsa
rsbadj <- rsb
rscadj <- rsc
rsaadj$Protein <- log(rsa$Protein) - logline[2,2]
rscadj$Protein <- log(rsc$Protein) + logline[2,2]
rsbadj$Protein <- log(rsb$Protein)
rsadj <- rbind(rsaadj,rsbadj,rscadj)

#reorder data
rn <- as.numeric(row.names(rsadj))
new1<-data.frame(rsadj, rn)
orderedadj = new1[do.call(order, new1["rn"]), ]

rn <- as.numeric(row.names(rs))
new1<-data.frame(rs, rn)
orderedorg = new1[do.call(order, new1["rn"]), ]

#compile table
table[i,1] <- format(var(rsadj$Protein), scientific = FALSE, digits = 7)

if (var(rsadj$Protein) < var(log(orderedorg$Protein))){
  table[i,2] <- "Yes"
}
if (var(rsadj$Protein) > var(log(orderedorg$Protein))){
  table[i,2] <- "No"
}

table[i,3] <- round(glance(var.test(orderedadj$Protein,log(orderedorg$Protein), alternative ="less"))$p.value,3)

WMW <- wilcox.test(orderedadj$Protein,log(orderedorg$Protein),paired=TRUE)
table[i,4] <- format(glance(WMW)$p.value, scientific = TRUE, digits = 4)

if (glance(WMW)$p.value< 0.01){
  table[i,5] <- "Yes"
}
if (glance(WMW)$p.value> 0.01){
  table[i,5] <- "No"
}

}
#write table
write.table(table, "C:\\Users\\sunck\\Documents\\plink_win32\\HDMS\\TLR4table.txt", sep = "\t")
write.table(orderedadj, "C:\\Users\\sunck\\Documents\\plink_win32\\HDMS\\TNFRSF17adj.txt", sep = "\t")
