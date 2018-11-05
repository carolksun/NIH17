library(dplyr)

Blood <- read.csv("Blood.csv", header = FALSE)
Nerve <- read.csv("Nerve.csv", header = FALSE)
CD48 <- read.table("chr1.txt", header = TRUE)

TLR4 <- read.table("chr9.txt", header = TRUE)

PRTN3 <- read.table("chr19.txt", header = TRUE)
FLT4 <- read.table("chr5.txt", header = TRUE)
CRK <- read.table("chr17.txt", header = TRUE)
MAPK14 <- read.table("chr6.txt", header = TRUE)

TNFRSF17 <- read.table("chr16.txt", header = TRUE)
PDCD1LG2 <- read.table("chr9.txt", header = TRUE)
TNC <- read.table("chr9.txt", header = TRUE)
MMP7 <- read.table("chr11.txt", header = TRUE)
colnames(Blood) <- c("SNP")
colnames(Nerve) <- c("SNP")
matchingMAPK14 <- inner_join(Blood, MAPK14)

nerveMMP7 <- inner_join(Nerve, MMP7)

mydata <- rbind(matchingcd48, matchingCRK, matchingFLT4, matchingMAPK14, matchingPRTN3, matchingtlr4)
mydata1 <- rbind(nerveCSF, nerveMMP7, nervePDCD1LG2, nervePLA2G7, nerveTNC, nerveTNFRSF17)
write.table(mydata1, "C:\\Users\\sunck\\Documents\\nervematching.txt", sep="\t")
