library(data.table)
library(janitor)
 unlist_split <- function(x, ...)
  {
	toret <- unlist(strsplit(x, ...) )
	return(t(toret))
  }

tdata='mrsc'

 dm1 <- fread(paste('data/p2_mrsc_CAPS-IV.csv',sep=''),data.table=F,na.strings=c(NA,"-9"))
 dm1  <- subset(dm1,IID !=2 )
 dm1 <- remove_empty(dm1, which = c("cols"))
 
 fam <- fread(paste('/home/pgcdac/DWFV2CJb8Piv_0116_pgc_data/pts/wave3/v1/',tdata,'/qc/pts_',tdata,'_mix_am-qc.fam',sep=''),data.table=F)
 names(fam) <- c("FID","IID","M","F","G","P")
 
 for (ancgroup in c("eur","lat","aam"))
 {
  pcs <-  read.table(paste('/home/pgcdac/DWFV2CJb8Piv_0116_pgc_data/pts/wave3/v1/',tdata,'/covariates/pts_',tdata,'_mix_am-qc-',ancgroup,'_pca.menv.mds_cov',sep=''),header=T,stringsAsFactors=F)
  
  dma <- merge(pcs,fam,by=c("FID","IID"),suffixes=c("_pc","")) 
  dm1x <- merge(dm1,dma,by=c("FID","IID"),suffixes=c("_fam","")) 
  
  
  dm1_exp <- subset(dm1x,select=c(FID,IID,Current_PTSD_Continuous))
  
  table(is.na(dm1_exp$Current_PTSD_Continuous))
  
  write.table(dm1_exp,file=paste("pheno/p2_",tdata,"_mrsc","_",ancgroup,".pheno",sep=""),quote=F,row.names=F)
    
  dm_cov <- dm1x  

  covexp4 <- subset(dm_cov,select=c("FID","IID","C1","C2","C3","C4","C5")) #mrs doesnt need sex!
  
  write.table(covexp4,file=paste("pheno/p2_",tdata,"_",ancgroup,"_mrsc","_pcs.cov",sep=""),quote=F,row.names=F)
     
  }
