args <- commandArgs(trailingOnly = TRUE)
 input_snp <- args[1]
 descriptor_file <-args[2]
 outfile <- args[3]

#Example code
# input_snp="rs762897"
# descriptor_file="eur_descriptors_july2_2021.csv"
# outfile="forest_plots/rs762897"

#Load metafor library and plyr (for mapped values)
library(metafor)
library(plyr)
library(data.table)

#Load study info. This should include how columns are named (to harmonize), study abbreviations, rescaling factors, study characteristics (for meta regression)
 descriptors <- read.csv(descriptor_file,header=T,stringsAsFactors=F)

#Make a matrix to store the study info
 input_data <- as.data.frame(matrix(nrow=dim(descriptors)[1], ncol=12))
 names(input_data ) <- c("study","SNP","CHR","BP","A1","A2","Freq","Info","Effect","SE","P","Neff")
 
#Load info from each study
#We will use the 'descriptors' sheet to pull in the correct column names to harmonize the data sources
for (i in 1:dim(descriptors))
{
 #Load SNP results
  tempdat <- read.table(paste('tophits_meta/',descriptors[i,]$File,".tophits",sep=""),stringsAsFactors=F,header=T,check.names = F)
 #Subset to just the plotting SNP
  tempdat$SNP <- tempdat[,descriptors[i,]$SNP]
 #Define a study variable for indexing
  tempdat$study <- descriptors[i,]$study
 if(input_snp %in% tempdat$SNP)
 {
  #print(input_snp %in% tempdat$SNP)
  tempdat <- subset(tempdat,SNP==input_snp) #subset data to just the relevant SNP
  
  tempdat$CHRN <- tempdat[,descriptors[i,]$CHR]
  tempdat$BPN <- tempdat[,descriptors[i,]$BP]
  tempdat$A1N <- tempdat[,descriptors[i,]$A1]
  tempdat$A2N <- tempdat[,descriptors[i,]$A2]
  tempdat$FreqN <- tempdat[,descriptors[i,]$Freq]
  tempdat$InfoN <- tempdat[,descriptors[i,]$SNP] #Get the info at some point
  tempdat$EffectN <- tempdat[,descriptors[i,]$Effect]
  
  #If the effect is named OR, transform to log OR!
   if(descriptors[i,]$Effect %in% c("OR","orbeta"))
   {
    tempdat$EffectN <- log(tempdat$EffectN)
   }
  
  tempdat$PN <- tempdat[,descriptors[i,]$P]
  
  #If the SE is NA, estimate it using P and Beta
   if(is.na(descriptors[i,]$SE))
   {
     Zv <- qnorm(tempdat$P,lower.tail=F)
     tempdat$SEN <- tempdat[,descriptors[i,]$Effect] / Zv
   } else {
    tempdat$SEN <- tempdat[,descriptors[i,]$SE]
   }
   
 #In case R is stupid, this turns TRUE back to T ALLELE!
  # tempdat$A1 <- as.factor(revalue(tempdat$A1, c("TRUE" = "T", "FALSE" = "F")))
  # tempdat$A2 <- as.factor(revalue(tempdat$A2, c("TRUE" = "T", "FALSE" = "F")))
  # tempdat$A1<- as.character(tempdat$A1)
  # tempdat$A2  <- as.character(tempdat$A2)

  #If N is specified, use that. If it's a column name, use that.
  if(!is.na(as.numeric(descriptors[i,]$N)))
   {
    tempdat$Neff <- descriptors[i,]$N
   } else {
    tempdat$Neff <- tempdat[,descriptors[i,]$N]
    }
   
  #Assign the harmonized data into the results matrix
  input_data[i,] <- subset(tempdat,select=c(study,SNP,CHRN,BPN,A1N,A2N,FreqN,InfoN,EffectN,SEN,PN,Neff))
 } else { 
   print (paste("Input snp ", input_snp, "is not in", tempdat$study))
   }
 #Remove the data in case the next file fails to load and it defaults to the existing file
 rm(tempdat)
}

#Only take data with SNP values

input_data <- subset(input_data, !is.na(Effect))
##Align alleles

#Set a reference for pivoting (Just choose row 1 in the data)
pivot=c(input_data[1,"A1"],input_data[1,"A2"])

#This is not coded for ambiguous alleles! That requires extra code to align based on frequencies #I can steal the code I already made when I compared fingen to esbb, in the ehr folder

input_data$Effect2 <- input_data$Effect
input_data$Freq2  <- input_data$Freq

#For every row, check allele alignment
for ( i in 1:dim(input_data)[1])
{
 betaval=input_data[i,]$Effect
 freqval=input_data[i,]$Freq
 snp=c(input_data[i,"A1"],input_data[i,"A2"])
 flipped=mapvalues(snp, c("A", "C", "G", "T"),  c("T", "G", "C", "A"))

 if (pivot[1] == snp[1] & pivot[2] == snp[2])
 {
  #No changes need to be made, alleles already aligned
     print(paste("No changes for", input_data[i,]$study,sep=" "))
 } else if (pivot[1] == snp[2] & pivot[2] == snp[1]) 
 {
  #A1 and A2 need to be flipped, change sign of beta value
   betaval=-1*betaval
   freqval=1-freqval
   print(paste("Direction flip for", input_data[i,]$study,sep=" "))
 } else if ( !(flipped[1] %in% pivot) |  !(flipped[2] %in% pivot)) 
 {
  print(paste("study is unflippable", input_data[i,]$study,sep=" "))
  #Check if allele might be flippable, if not, give NA value 
  betaval=NA
 } else if (flipped[1] == snp[1] &  flipped[2] == snp[2]) 
 {
    print(paste("Wrong strand but correct flip for ", input_data[i,]$study,sep=" "))
   #It's on the wrong strand but nothing actually needs to be done
 } else if (flipped[1] == snp[2] &  flipped[2] == snp[1]) 
 {
  #It's on the wrong strand and needs to be flipped
   print(paste("Direction and allele flip for", input_data[i,]$study,sep=" "))
   betaval=-1*betaval
   freqval=1-freqval
 } else {
 betaval=NA
 }
 
 #Now reassign adjusted beta
 input_data[i,]$Effect2 <- betaval
 input_data[i,]$Freq2 <- freqval

}

#Merge the input data and descriptors sheet for meta regression
input_data_desc <- merge(descriptors,input_data,by="study",all.x=T,suffixes=c("_annotat","")) #possibly remove the all.x at some point, useful for book keeping but otherwise annoying

#Rescale effect sizes for quantitative scores
input_data_desc$Effect3  <- input_data_desc$Effect2
input_data_desc$SE3 <- input_data_desc$SE

input_data_desc[which(input_data_desc$Included.As == "Quantitative"),]$Effect3 <- input_data_desc[which(input_data_desc$Included.As == "Quantitative"),]$Effect2 * as.numeric(input_data_desc[which(input_data_desc$Included.As == "Quantitative"),]$ptsd_scale) * 100
input_data_desc[which(input_data_desc$Included.As == "Quantitative"),]$SE3 <- input_data_desc[which(input_data_desc$Included.As == "Quantitative"),]$SE * as.numeric(input_data_desc[which(input_data_desc$Included.As == "Quantitative"),]$ptsd_scale) * 100




input_data_desc <- input_data_desc[order(as.numeric(input_data_desc$Neff),decreasing=c(T)),]

#Which studys are NOT in the original input_data_desca 
print("Not included (missing input_data_desca):")
print(descriptors[-which(descriptors$study %in% input_data$study),]$study)

print("N Studies analyzed:")
print(dim(input_data)[1]) 


#These will match if the number of inputs in the starting directory 
#matches the number of inputs in the study description file

input_data_desc_qt <- subset(input_data_desc,Included.As == "Quantitative")
input_data_desc_cc <- subset(input_data_desc,Included.As != "Quantitative")


meta_qt <- rma(yi=input_data_desc_qt$Effect3, sei=input_data_desc_qt$SE3, slab=input_data_desc_qt$study,method="FE" )
meta_cc <- rma(yi=input_data_desc_cc$Effect3, sei=input_data_desc_cc$SE3, slab=input_data_desc_cc$study ,method="FE")

#this should be it for SSW
#meta_zscore <- rma(yi=input_data_desc$Effect3/input_data_desc$SE3*sqrt(as.numeric(input_data_desc$Neff)), sei=1/sqrt(as.numeric(input_data_desc$Neff)), slab=input_data_desc$study ,method="FE")

pdf(paste(outfile,'.pdf',sep=''),7,9)
plot(meta_qt)
plot(meta_cc)

dev.off()

write.table(input_data_desc, paste(outfile,'.txt',sep=''),row.names=F,quote=F)


