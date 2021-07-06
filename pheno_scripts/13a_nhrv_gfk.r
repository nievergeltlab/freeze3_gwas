library(data.table)
library(janitor)

 unlist_split <- function(x, ...)
  {
	toret <- unlist(strsplit(x, ...) )
	return(t(toret))
  }
  
tdata="nhrv"

 dm1 <- fread(paste('data/p3_NHRV_GFK_13.csv',sep=''),data.table=F,na.strings=c(NA,"#N/A","-9"))
 
dm1 <- remove_empty(dm1, which = c("cols"))

 fam <- fread(paste('/home/pgcdac/DWFV2CJb8Piv_0116_pgc_data/pts/wave3/v1/',tdata,'/qc/pts_',tdata,'_mix_am-qc.fam',sep=''),data.table=F)
 names(fam) <- c("FID","IID","M","F","G","P")
 
 dm1$Lifetime_PTSD_Continuous <- dm1$pcl_total_lifetime


 for (ancgroup in c("eur"))
 {
  pcs <-  read.table(paste('/home/pgcdac/DWFV2CJb8Piv_0116_pgc_data/pts/wave3/v1/',tdata,'/covariates/pts_',tdata,'_mix_am-qc-',ancgroup,'_pca.menv.mds_cov',sep=''),header=T,stringsAsFactors=F)
  
  dma <- merge(pcs,fam,by=c("FID","IID"),suffixes=c("_pc","")) 
  dm1x <- merge(dm1,dma,by=c("IID"),suffixes=c("","_fam"))  
  
  
  dm1_exp <- subset(dm1x,!is.na(Lifetime_PTSD_Continuous),select=c(FID,IID,Lifetime_PTSD_Continuous))
 
  table(is.na(dm1_exp$Lifetime_PTSD_Continuous))
  
  #now 604
  write.table(dm1_exp,file=paste("pheno/p2_",tdata,"_gfk","_",ancgroup,".pheno",sep=""),quote=F,row.names=F)
    
  dm_cov <- dm1x  

  covexp4 <- subset(dm_cov,select=c("FID","IID","C1","C2","C3","C4","C5")) #mrs doesnt need sex!
  
  write.table(covexp4,file=paste("pheno/p2_",tdata,"_",ancgroup,"_gfk","_pcs.cov",sep=""),quote=F,row.names=F)
 }  
  