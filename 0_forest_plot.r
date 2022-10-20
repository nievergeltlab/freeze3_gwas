args <- commandArgs(trailingOnly = TRUE)
 snpresults <- args[1]
 descriptor <-args[2]
  outfile <- args[3]
 ancestry <- args[4]


# snpresults="results_filtered/forest_plots/rs6802567.allstudies"
# descriptor="f3_studies_forest_descriptors.txt"
# outfile="results_filtered/forest_plots/rs6802567"
# ancestry="eur"


#Load metafor library and plyr (for mapped values)
library(metafor)
library(plyr)
library(data.table)

#study info
descriptors <- fread(descriptor,data.table=F)

#Subset to ancestry if need be
if(ancestry !="all")
{
 descriptors <- subset(descriptors,ancestry==ancestry)
}

#Load data (should be in described format of file name followed by 
 d1 <- fread(paste(snpresults,sep=""),data.table=F)


dat0 <- subset(d1,study %in% descriptors$study & !is.na(P))

#all studies get rescaled based on Z score
dat0$Zscore <- abs(qnorm(dat0$P/2,lower.tail=F)) * sign(dat0$Effect)
sed=3.28
#use Rescaling code..
 dat0$B = dat0$Zscore * sed / sqrt(2* dat0$N* dat0$MAF*(1- dat0$MAF))
 dat0$SE <-  dat0$B /  dat0$Zscore 
 
 dat0$B2 <- NA
##Align alleles

#Set a reference for pivoting (Just choose row 1 in the data)
pivot=c(dat0[1,"A1"],dat0[1,"A2"])

#this is not coded for C G alleles!!

#For every row, check allele alignment
#No code for ambiguous alleles?

for ( i in 1:dim(dat0)[1])
{
 betaval=dat0[i,]$B
 snp=c(dat0[i,"A1"],dat0[i,"A2"])
 flipped=mapvalues(snp, c("A", "C", "G", "T"),  c("T", "G", "C", "A"))

 if (pivot[1] == snp[1] & pivot[2] == snp[2])
 {
  #No changes need to be made, alleles already aligned
     print(paste("No changes for", dat0[i,]$Study,sep=" "))
 } else if (pivot[1] == snp[2] & pivot[2] == snp[1]) 
 {
  #A1 and A2 need to be flipped, change sign of beta value
   betaval=-1*betaval
   print(paste("Direction flip for", dat0[i,]$Study,sep=" "))
 } else if ( !(flipped[1] %in% pivot) |  !(flipped[2] %in% pivot)) 
 {
  print(paste("Study is unflippable", dat0[i,]$Study,sep=" "))
  #Check if allele might be flippable, if not, give NA value 
  betaval=NA
 } else if (flipped[1] == snp[1] &  flipped[2] == snp[2]) 
 {
    print(paste("Wrong strand but correct flip for ", dat0[i,]$Study,sep=" "))
   #It's on the wrong strand but nothing actually needs to be done
 } else if (flipped[1] == snp[2] &  flipped[2] == snp[1]) 
 {
  #It's on the wrong strand and needs to be flipped
   print(paste("Direction and allele flip for", dat0[i,]$Study,sep=" "))
   betaval=-1*betaval
 } else {
 betaval=NA
 }
 
 #Now reassign adjusted beta
 dat0[i,]$B2 <- betaval

}



#Merge in study descriptors for the plot

descriptors$Study = descriptors$filename2              
dat <- merge(dat0,descriptors,by="study",suffixes=c("","_descriptor"))



dat <- dat[order(dat$Type,dat$Abbr_meta,decreasing=c(F,T)),]

### a little helper function to add Q-test, I^2, and tau^2 estimate info

#These will match if the number of inputs in the starting directory 
#matches the number of inputs in the study description file


results <- rma(yi=dat$B2, sei=dat$SE, method="FE",data=dat,slab=dat$Abbr_meta)

#analysis within MVP, EHR, 2.5
results_ehr   <- rma(yi=dat$B2, sei=dat$SE, method="FE",data=dat,slab=dat$Abbr_meta,subset=(Type=="BEHR"))
results_25  <- rma(yi=dat$B2, sei=dat$SE, method="FE",data=dat,slab=dat$Abbr_meta,subset=(Type=="AF2.5"))
results_mvp  <- rma(yi=dat$B2, sei=dat$SE, method="FE",data=dat,slab=dat$Abbr_meta,subset=(Type=="MVP"))
sink(file=paste(outfile,'.log',sep='_'))

print(results)

#Which studys are NOT in the original data 
print("Not included (missing data):")
print(descriptors[-which(descriptors$study %in% dat$study),]$study)

print("N Studies analyzed:")
print(dim(dat)[1])


#List the variant name and coded allele on the plot somewhere
paste(dat0$SNP[1], " (",dat0[1,"A1"],")",sep="")
#save(meta_res_qt,file=paste(outfile,"_fp",'.R',sep=''))

#Need dimenssion of data for plotting purposes
nstudies=dim(dat)[1]

#need n groups
ngroups=3
#Determine, given ordered data, on which rows the subsets would stop. add blank spots to allow for meta analysis of subgroups


pdf(paste(outfile,'.pdf',sep=''),6,11)

plotstartpoint=results$k+ngroups*2*2  # #The limits of the graph include extra spaces at the top of the plot, plus spaces for the headings, see forest.rma in https://cran.r-project.org/web/packages/metafor/metafor.pdf
#the second *2 is for double spacing

rowpoints <- c((plotstartpoint-4):(plotstartpoint-4-results_25$k+1), (plotstartpoint-4-results_25$k-4):(plotstartpoint-4-results_25$k-results_ehr$k-3),(plotstartpoint-4-results_25$k-results_ehr$k-8))

forest(results,header="Study",rows= rowpoints,decreasing=TRUE,ylim=c(-1,plotstartpoint) ,cex=0.75,xlim=c(-6, 6),col="blue",border="blue",mlab="Freeze 3 Meta" )


#notice how I cut out an empty space for the first meta... now for the second...
#Add text for subgroups as well

text(-6.2, c(plotstartpoint-3,(plotstartpoint-4-results_25$k-3),(plotstartpoint-4-results_25$k-results_ehr$k-7)), pos=4, c("Freeze 2.5",
                               "EHR",
                               "MVP"))
                               
addpoly(results_25, mlab="F2.5 Meta", row= (plotstartpoint-4-results_25$k-1),col="blue",border="blue" )
addpoly(results_ehr,  mlab="EHR Meta",  row=(plotstartpoint-4-results_25$k-results_ehr$k-5),col="blue",border="blue" )
#addpoly(results_mvp , row= (results$k+ngroups-results_25$k-results_ehr$k)-3)
dev.off()

pdf(paste(outfile,'_diagnostic.pdf',sep=''),8,8)
plot(results)
dev.off()



# #Write output file for metasoft
# dout <- unlist(matrix(c(dat$B, dat$SE),ncol=2,byrow=F))

# write.table(t(c(outfile,t(dout + 0.0000001))),paste("results_cat/", outfile,"_",anc,'.msin',sep=''),quote=F,row.names=F,col.names=F)
# write.table(paste(dat$abbr,sep="" ),paste("results_cat/", outfile,"_",anc,'.msinstudynames',sep=''),quote=F,row.names=F,col.names=F) # Label by Study Number and abbr
# write.table(dat$caroline_order2,paste("results_cat/", outfile,"_",anc,'.msinstudyorder',sep=''),quote=F,row.names=F,col.names=F) #Order by number of cases

# print(warnings())


# sum(dat$B * 1/dat$SE^2) / sum(1/dat$SE^2)

