library(dplyr)
library(data.table)

#import Peter data
Peter1 <- read.csv("Peter1.csv", header=FALSE)
Peter2 <- read.csv("Peter2.csv", header=FALSE)
Peter3 <- read.csv("Peter3.csv", header=FALSE)
Peter4 <- read.csv("Peter4.csv", header=FALSE)

#import online SNPs
CSNPs <- read.csv("SNPs_C.csv", header=FALSE)
MSNPs <- read.csv("SNPs_M.csv", header=FALSE)
combined <- read.csv("Combined.csv", header=FALSE)
proteins <- read.csv("proteins.csv", header=FALSE, stringsAsFactors = FALSE)

#delete column V2
matching1$V2 <- NULL
matching2$V2 <- NULL
matching3$V2 <- NULL
matching4$V2 <- NULL
MSNPs$V2 <- NULL
combined$V2 <- NULL

#Matching using inner_join
matching1 <- inner_join(combined, Peter1)
matching2 <- inner_join(combined, Peter2)
matching3 <- inner_join(combined, Peter3)
matching4 <- inner_join(combined, Peter4)
matchingprotein1 <- inner_join(proteins, Peter1)
matchingprotein2 <- inner_join(proteins, Peter2)
matchingprotein3 <- inner_join(proteins, Peter3)
matchingprotein4 <- inner_join(proteins, Peter4)

#Matching using fintersect
fmatching1 <- fintersect(setDT(combined), setDT(Peter1))
fmatching2 <- fintersect(setDT(combined), setDT(Peter2))
fmatching3 <- fintersect(setDT(combined), setDT(Peter3))
fmatching4 <- fintersect(setDT(combined), setDT(Peter4))

#combine all matching
total <- rbind(matching1,matching2,matching3,matching4)

#check indexes
index <- which(Peter1$V1=="rs12501028")

#for loop to check indexes
for(i in 1:nrow(total)){
  if(length(which(Peter1$V1==total[i,1]))==1 | length(which(Peter2$V1==total[i,1]))==1 | 
     length(which(Peter3$V1==total[i,1]))==1 | length(which(Peter4$V1==total[i,1]))==1){
    print(TRUE)
  }else{
    print(FALSE)    
  }
}
#importing Meghana's matches
Mmatches <- read.csv("Mmatches.csv", header=FALSE, stringsAsFactors = FALSE)
#checking Meghana's matches
for(i in 1:nrow(Mmatches)){
  if(length(which(Peter1$V1==Mmatches[i,1]))==1 | length(which(Peter2$V1==Mmatches[i,1]))==1 | 
     length(which(Peter3$V1==Mmatches[i,1]))==1 | length(which(Peter4$V1==Mmatches[i,1]))==1){
    print(TRUE)
  }else{
    print(FALSE)    
  }
}