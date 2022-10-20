args <- commandArgs(trailingOnly = TRUE)
 codingfile <- args[1]
 snpresults <- args[2]
 outname <-args[3]
 snp <- args[4] 
 
   codingfile="results_filtered/forest_plots/biov_biov_eur.coding"
  snpresults="results_filtered/forest_plots/biov_biov_eur.rs78201023"
   outname="results_filtered/forest_plots/biov_biov_eur.rs78201023.formatted"
   snp="rs78201023"

library(data.table)

codingdata <- fread(codingfile,stringsAsFactors=F,header=T,nrow=2,data.table=F)
snpdata <- fread(snpresults,stringsAsFactors=F,header=T,nrow=2,data.table=F, sep=c(" ","\t"))

rfdata <- as.data.frame(matrix(nrow=1,ncol=7))
names(rfdata ) <- c("SNP","A1","A2","MAF","Effect","N","P") #reformatted dataset

rfdata$SNP <- snp
rfdata$study <- codingdata$study[1]
rfdata$ancestry <- codingdata$ancestry[1]

#Only do if SNP found.. otherwise export NAs
if(dim(snpdata)[1] >0)
{

#Reformat odds ratios here. 
if(codingdata$EFFECT[1] %in% c("OR", "orbeta"))
{
 rfdata$Effect <-  log(snpdata[,codingdata$EFFECT[1]])
} else {
 rfdata$Effect <-  snpdata[,codingdata$EFFECT[1]]
}

#If it is numeric in the coding file, don't use the data in the SNP file. Instead give N in the coding data version. If it is a column name, however, then use the SNP file.
if(is.numeric(codingdata$WEIGHT))
{ 
 rfdata$N <- codingdata$WEIGHT
} else {
rfdata$N <- snpdata[,codingdata$WEIGHT[1]]
}


rfdata$SNP <- snpdata[,codingdata$SNP[1]]
rfdata$A1 <- snpdata[,codingdata$A1[1]]
rfdata$A2 <- snpdata[,codingdata$A2[1]]
rfdata$MAF <- snpdata[,codingdata$MAF[1]] #These are reformatted in the meta-analysis step
rfdata$P <- snpdata[,codingdata$P[1]] #These are reformatted in the meta-analysis step

}

write.table(subset(rfdata,select=c(study,ancestry,SNP,A1,A2,MAF,Effect,P,N)),file=outname,quote=F,row.names=F)
