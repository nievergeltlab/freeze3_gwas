library(data.table)
library(janitor)
library(plyr)
 unlist_split <- function(x, ...)
  {
	toret <- unlist(strsplit(x, ...) )
	return(t(toret))
  }
  
 tdata="psy3"

 dm1a <- fread(paste('data/p2_feen.csv',sep=''),data.table=F,na.strings=c(NA,"-9"))
 dm1b <- fread(paste('data/p2_dcsr.csv',sep=''),data.table=F,na.strings=c(NA,"-9"))
 dm1c1 <- fread(paste('data/p2_teic_Kaufman__updated_20190220.csv',sep=''),data.table=F,na.strings=c(NA,"-9")) #This is 39_TEIC1_Study-level data form_PGC-PTSD.docx
 dm1c2 <- fread(paste('data/p2_teic_ls_sn.csv',sep=''),data.table=F,na.strings=c(NA,"-9")) #These subjects have no PTSD QT measure!! #this is 39_TEIC_Study-level data form_PGC-PTSD.docx
 
 #dm1c2 
 dm1a <- remove_empty(dm1a, which = c("cols"))
 dm1b <- remove_empty(dm1b, which = c("cols"))
 dm1c1 <- remove_empty(dm1c1, which = c("cols"))
 dm1c2 <- remove_empty(dm1c2, which = c("cols"))
 
 dm1a$Current_PTSD_Continuous_harmonized <- dm1a$Current_PTSD_Continuous / 51
 dm1b$Current_PTSD_Continuous_harmonized <- dm1b$Current_PTSD_Continuous / 136
 dm1c1$Current_PTSD_Continuous_harmonized <-  dm1c1$Current_PTSD_Continuous / 80 
 dm1c2$Current_PTSD_Continuous_harmonized <- 0#everyone is a control..., impute to 0
 
dm1a$LT_Count_harmonized <- dm1a$LT_Count / 10
dm1b$LT_Count_harmonized <- dm1b$LT_Count / 11
dm1c1$LT_Count_harmonized <- dm1c1$LT_Count / 20
dm1c2$LT_Count_harmonized <- dm1c2$LT_Count / 12
   
 
 
 dm1 <- rbind.fill(dm1a,dm1b,dm1c1,dm1c2)
 
 
 ##harmonize to 0-1 based on theoretical ranges
 
 #split up teic kaufman etc based on substudy
 
 fam <- fread(paste('/home/pgcdac/DWFV2CJb8Piv_0116_pgc_data/pts/wave3/v1/',tdata,'/qc/pts_',tdata,'_mix_am-qc.fam',sep=''),data.table=F)
 names(fam) <- c("FID","IID","M","F","G","P")

 for (ancgroup in c("eur","aam"))
 {
  pcs <-  read.table(paste('/home/pgcdac/DWFV2CJb8Piv_0116_pgc_data/pts/wave3/v1/',tdata,'/covariates/pts_',tdata,'_mix_am-qc-',ancgroup,'_pca.menv.mds_cov',sep=''),header=T,stringsAsFactors=F)
  
  dma <- merge(pcs,fam,by=c("FID","IID"),suffixes=c("_pc","")) 
  dm1x <- merge(dm1,dma,by=c("FID","IID"),suffixes=c("_fam","")) 
  
  
  dm1_exp <- subset(dm1x,!is.na(Case),select=c(FID,IID,Case))
 
  print(ancgroup)
 print(  table(dm1_exp$Case))
  
  print(table(is.na(dm1_exp$Case)))
    print(table(dm1_exp$Case))
     counts <- table(dm1_exp$Case)
  4/((1/counts[1]) + (1/counts[2]))

  
  write.table(dm1_exp,file=paste("pheno/p2_",tdata,"_feen","_",ancgroup,".pheno",sep=""),quote=F,row.names=F)
    
  dm_cov <- dm1x  

  covexp4 <- subset(dm_cov,select=c("FID","IID","C1","C2","C3","C4","C5")) #mrs doesnt need sex!
  
  write.table(covexp4,file=paste("pheno/p2_",tdata,"_",ancgroup,"_feen","_pcs.cov",sep=""),quote=F,row.names=F)
 }  
  
  