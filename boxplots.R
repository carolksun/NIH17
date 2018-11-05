#Creating boxplots and plotting Median-Based Linear Model
#Carol Sun
#July 26th, 2017

#------------------------------------------------------------------------------
#load in packages
#------------------------------------------------------------------------------

library(ggplot2)
library(broom)
library(gridExtra)
library(mblm)

#------------------------------------------------------------------------------
#load in alleles and Protein Levels
#------------------------------------------------------------------------------

df <- read.csv("PDCDprotein.csv", header = TRUE)
df <- read.csv("TNFR.csv", header = TRUE)
df <- read.csv("MAPK14.csv", header = TRUE)
df <- read.csv("CRK.csv", header = TRUE)
df <- read.csv("PLA.csv", header = TRUE)
df <- read.csv("FLT4.csv", header = TRUE)

df <- read.csv("PDCDlittlepnerve.csv", header = TRUE)
df <- read.csv("PDCDeven.csv", header = TRUE)
dfnew <- read.csv("PDCDlittlepall.csv", header = TRUE)
dfnew <- read.csv("PDCDevenall.csv", header = TRUE)
pvalues <- read.csv("pvalues.csv", header = TRUE)

colnames(df) <- c("Diagnosis", colnames(df)[i+1],"Protein")
#------------------------------------------------------------------------------
#update title center and plot font
#------------------------------------------------------------------------------

theme_update(plot.title = element_text(hjust = 0.5))
theme_update(text=element_text(size=16, family = "serif"))

#------------------------------------------------------------------------------
#For loop to make plots and save lines
#------------------------------------------------------------------------------

for(i in 19:19){
  #------------------------------------------------------------------------------
  #split data frame based on alleles
  #------------------------------------------------------------------------------
  
  df.split<-split(df,df[[i+1]])
  rsa <- df.split[[1]][,c("Diagnosis", colnames(df)[i+1], "Protein")]
  #------------------------------------------------------------------------------
  #if/else statement to fix problem when allele is missing and replaced by 0 0
  #------------------------------------------------------------------------------
  
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

  #------------------------------------------------------------------------------
  #if/else statement to create list to organize factors in order of major -> minor
  #------------------------------------------------------------------------------    
  
  if(nrow(rsa) < nrow(rsc)){
    fcts <- c(as.character(rsc[1,2]), as.character(rsb[1,2]), as.character(rsa[1,2]))
  }else{
    fcts <- c(as.character(rsa[1,2]), as.character(rsb[1,2]), as.character(rsc[1,2]))
  }
  #------------------------------------------------------------------------------
  #Create table of MS/HD and alleles
  #------------------------------------------------------------------------------   
  
  counta1 <- length(which(rsc$Diagnosis == "MS"))
  counta2 <- length(which(rsc$Diagnosis == "non-MS")) 
  countb1 <- length(which(rsb$Diagnosis == "MS"))
  countb2 <- length(which(rsb$Diagnosis == "non-MS"))
  countc1 <- length(which(rsa$Diagnosis == "MS"))
  countc2 <- length(which(rsa$Diagnosis == "non-MS"))
  
  table <-matrix(c(counta1,counta2, countb1, countb2, countc1, countc2),ncol=2,byrow=TRUE)
  colnames(table) <- c("MS", "non-MS")
  rownames(table) <- fcts
  table <- as.table(table)
  
  tt <- ttheme_minimal(
    core=list(bg_params = list(fill = blues9[1:4], col=NA),
            fg_params=list(fontface=3)),
    colhead=list(fg_params=list(col="navyblue", fontface=4L)))
  tbl <- tableGrob(table, rows=rownames(table), cols = colnames(table), theme=tt)
  
  #------------------------------------------------------------------------------   
  #Generate MBLM (Median-Based Linear Model) for ln and normal values
  #------------------------------------------------------------------------------   
  
  elog = log(as.numeric(rs$Protein))
  slog = as.numeric(factor(rs[[2]],fcts))
  logmlm = mblm(elog ~ slog, repeated = F)
  logline = tidy(logmlm)
      
  efake = log(as.numeric(rsadj$Protein))
  sfake = as.numeric(factor(rsadj[[2]],fcts))
  logfake = mblm(efake ~ sfake, repeated = F)
  logfake = tidy(logfake)
  #------------------------------------------------------------------------------   
  #save separated and combined box plots
  #------------------------------------------------------------------------------   
  rsafake <- rsa
  rsbfake <- rsb
  rscfake <- rsc
  rsafake$Protein <- log(rsa$Protein) + logline[2,2]
  rscfake$Protein <- log(rsc$Protein) - logline[2,2]
  rsbfake$Protein <- log(rsb$Protein)
  logrs <- rbind(rsafake,rsbfake,rscfake)
  logrs$Protein <- exp(logrs$Protein)
  pdf("C:\\Users\\sunck\\Documents\\plink_win32\\HDMS\\newgraphs.pdf")
  print(ggplot(data = rs, aes(factor(rs[[2]],fcts), log(Protein))) 
          + geom_boxplot(aes(fill=Diagnosis)) + ylab("TNC Expression (ln(RFU))") + xlab("Alleles")
          + ggtitle(paste(colnames(df)[i+1], "(p=", pvalues[which(pvalues$SNP == colnames(df)[i+1]),2], ")")) 
          + geom_abline(slope = logline[2,2], intercept = logline[1,2])
          + annotate(geom="text", x=3, y=5.5, label=paste("slope = ", round(logline[2,2], digits = 5)), colour="navy",
             size=6.5, family="serif", fontface="plain")
          )
  
  print(ggplot(data = rs, aes(factor(rs[[2]],fcts), log(Protein))) 
          + geom_boxplot() + ylab("TNC Expression (ln(RFU)") + xlab("Alleles")
          + ggtitle(paste(colnames(df)[i+1], "(p=", pvalues[which(pvalues$SNP == colnames(df)[i+1]),2], ")")) 
          + geom_abline(slope = logline[2,2], intercept = logline[1,2])
          + annotate(geom="text", x=3, y=5.5, label=paste("slope = ", round(logline[2,2], digits = 5)), colour="navy",
             size=6.5, family="serif", fontface="plain")
          )
  
  print(ggplot(data = rsadj, aes(factor(rsadj[[2]],fcts), Protein)) 
            + geom_boxplot() + ylab("TNC Expression (ln(RFU)") + xlab("Alleles")
            + ggtitle(paste(colnames(df)[i+1], "(p=", pvalues[which(pvalues$SNP == colnames(df)[i+1]),2], ")")) 
            + annotate(geom="text", x=3, y=5.5, label=paste("slope = 0"), colour="navy",
                       size=6.5, family="serif", fontface="plain")
  )
  
print(ggplot(data = rsadj, aes(factor(rsadj[[2]],fcts), Protein)) 
    + geom_boxplot(aes(fill=Diagnosis)) + ylab("TNC Expression (ln(RFU)") + xlab("Alleles")
    + ggtitle(paste(colnames(df)[i+1], "(p=", pvalues[which(pvalues$SNP == colnames(df)[i+1]),2], ")")) 
    + annotate(geom="text", x=3, y=5.5, label=paste("slope = 0"), colour="navy",
               size=6.5, family="serif", fontface="plain")
  )
dev.off()
  #------------------------------------------------------------------------------   
  #Export to pdf
  #------------------------------------------------------------------------------ 
  
  pdf("C:\\Users\\sunck\\Documents\\graphs.pdf")
  grid.arrange(plot1)
  grid.arrange(plot2)
  grid.arrange(plot3)
  grid.arrange(plot4)
  grid.arrange(tbl)

  dev.off()
}

print(tbl)
