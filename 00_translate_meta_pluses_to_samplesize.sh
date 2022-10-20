#Given +/- and ? in meta results, and an ordered list with sample sizes that they correspond to, get sample size counts
library(data.table)
#d1 <- fread('zcat eur_ptsdcasecontrol_pcs_v5_jan4_2022.temp1.tbl',data.table=F)

#Make sure that this is sorted by the data entry order into meta analysis
#samplecounts <- fread('grotzinger_case_control_counts.csv',data.table=F)
d1 <- fread('zcat eur_ptsdcasecontrolnomvp_pcs_v4_aug3_2021.temp1.tbl ',data.table=F)
samplecounts <- fread('grotzinger_case_control_counts_nomvp.csv',data.table=F)
samplecounts$neff <- samplecounts$nef
samplecounts$ncontrol <- samplecounts$ncon

samplecounts <- samplecounts[which(samplecounts$study12 != ""),]

#convert +/- to equal 1, ? to equal 0
d1$metal1 <- sapply(d1$Direction,gsub,pattern="0",replacement="1 ",fixed=TRUE) #sometimes effects have sign NULL, coded as 0. Since the study was entered, the indicator is 1
d1$metal1 <- sapply(d1$metal1,gsub,pattern="+",replacement="1 ",fixed=TRUE)
d1$metal1 <- sapply(d1$metal1,gsub,pattern="-",replacement="1 ",fixed=TRUE)
d1$metal1 <- sapply(d1$metal1,gsub,pattern="?",replacement="0 ",fixed=TRUE)


#Splitfun 

splitfun <- function(x)
{
 as.numeric(unlist(strsplit(x,split=" ")))
}

#Now for each row, multiply the split by the matrix

splitmat <- function(x)
{
 x <- splitfun(x)

 casecount <- x %*% samplecounts$ncase

 controlcount <- x %*% samplecounts$ncontrol
 neffcount <- x %*% samplecounts$neff
 return(c(casecount,controlcount,neffcount))
}

d1$nchar <- sapply(d1$metal1,nchar) #all lines should have the same number of attempted studies. If this is off, the string replacemnet is wrong

casecounts <- t(sapply(d1$metal1,splitmat))

#for each row, transpose the metal +/- colmn. 
#Do a strsplit

# take hte product of sample size and these 1s, this returns sample size

d1$ncase <- as.vector(casecounts[,1])
d1$ncontrol <- casecounts[,2]
d1$neff <- casecounts[,3]
d1$ntot <- d1$ncase + d1$ncontrol

 totalN=463244
 totalN=342624

 #totalN=51252
 
 percentN=0.8
 d1_exp <- subset(d1,neff >= totalN*percentN & Freq1 >= 0.01 & Freq1 <= 0.99)

#write.table(d1_exp[,-which(names(d1_exp) %in% c("nchar","metal1"))],file='eur_ptsdcasecontrol_pcs_v5_jan4_2022.temp1.tbl.neffX',quote=F,row.names=F)
write.table(d1_exp[,-which(names(d1_exp) %in% c("nchar","metal1"))],file='eur_ptsdcasecontrolnomvp_pcs_v5_jan4_2022.tbl.neff',quote=F,row.names=F)

#write.table(d1_exp[,-which(names(d1_exp) %in% c("nchar","metal1"))],file='aam_ptsdcasecontrol_pcs_v5_jan4_2022.temp1.tbl.neff',quote=F,row.names=F)


test<- gsub(x=gsub(x=gsub(x=metalcheck,"+","1 ",fixed=TRUE),"-","1 ",fixed=TRUE),"?","0 ",fixed=TRUE)

as.numeric(unlist(strsplit(test,split=" ")))

