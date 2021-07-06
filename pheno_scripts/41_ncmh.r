library(data.table)
library(janitor)
 unlist_split <- function(x, ...)
  {
	toret <- unlist(strsplit(x, ...) )
	return(t(toret))
  }
  
tdata="psy3"

 dm1 <- fread(paste('data/p3_NCMH_20190320.csv',sep=''),data.table=F,na.strings=c(NA,"-9"))
 
dm1 <- remove_empty(dm1, which = c("cols"))

 fam <- fread(paste('/home/pgcdac/DWFV2CJb8Piv_0116_pgc_data/pts/wave3/v1/',tdata,'/qc/pts_',tdata,'_mix_am-qc.fam',sep=''),data.table=F)
 names(fam) <- c("FID","IID","M","F","G","P")

 dm1$Current_PTSD_Continuous <- dm1$PTSD2E_1+ dm1$PTSD2E_2+ dm1$PTSD2E_3+ dm1$PTSD2E_4	+ dm1$PTSD2E_5	+ dm1$PTSD2E_6	+ dm1$PTSD2E_7	+ dm1$PTSD2E_8	+ dm1$PTSD2E_9	+ dm1$PTSD2E_10
 

 for (ancgroup in c("eur", "aam")) #,"lat","aam"))
 {
  pcs <-  read.table(paste('/home/pgcdac/DWFV2CJb8Piv_0116_pgc_data/pts/wave3/v1/',tdata,'/covariates/pts_',tdata,'_mix_am-qc-',ancgroup,'_pca.menv.mds_cov',sep=''),header=T,stringsAsFactors=F)
  
  dma <- merge(pcs,fam,by=c("FID","IID"),suffixes=c("_pc","")) 
  dm1x <- merge(dm1,dma,by=c("IID"),suffixes=c("_fam",""))  #iid only 
  
  
  dm1_exp <- subset(dm1x,!is.na(Current_PTSD_Continuous),select=c(FID,IID,Current_PTSD_Continuous))
 
  table(is.na(dm1_exp$Current_PTSD_Continuous))
  
  
  write.table(dm1_exp,file=paste("pheno/p2_",tdata,"_ncmh","_",ancgroup,".pheno",sep=""),quote=F,row.names=F)
    
  dm_cov <- dm1x  

  covexp4 <- subset(dm_cov,select=c("FID","IID","C1","C2","C3","C4","C5")) #mrs doesnt need sex!
  
  write.table(covexp4,file=paste("pheno/p2_",tdata,"_",ancgroup,"_ncmh","_pcs.cov",sep=""),quote=F,row.names=F)
 }  
  