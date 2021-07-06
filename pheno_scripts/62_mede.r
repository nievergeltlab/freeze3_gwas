library(data.table)
library(janitor)
 unlist_split <- function(x, ...)
  {
	toret <- unlist(strsplit(x, ...) )
	return(t(toret))
  }
  
tdata="mede"

 dm1 <- fread(paste('data/medeIndividual_data_CVBPGCPTSD_Your_Study_final_.csv',sep=''),data.table=F,na.strings=c(NA,"-9"),colClasses=c("IID"="character"))

 #THIS IS A HUGE PAIN IN THE ASS WITH THE LEADING ZEROS... just take the one already encoded in the FAM FILE!
 
#Dont use CurrentDX for STARRS, its not great
 #dm1 <- fread(paste('data/p2_',tdata,'.csv',sep=''),data.table=F,na.strings=c(NA,"-9"))
 
dm1 <- remove_empty(dm1, which = c("cols"))
dm1$IID_pc <- dm1$IID



 fam <- fread(paste('/home/pgcdac/DWFV2CJb8Piv_0116_pgc_data/pts/wave3/v1/',tdata,'/qc/pts_',tdata,'_mix_am-qc1.fam',sep=''),data.table=F)
 names(fam) <- c("FID","IID","M","F","G","P")
 
 for (ancgroup in c("lat","pue"))
 {
  pcs <-  read.table(paste('/home/pgcdac/DWFV2CJb8Piv_0116_pgc_data/pts/wave3/v1/',tdata,'/covariates/pts_',tdata,'_mix_am-qc-',ancgroup,'_pca.menv.mds_cov',sep=''),header=T,stringsAsFactors=F)
  
  dma <- merge(pcs,fam,by=c("FID"),suffixes=c("_pc","")) 
  dm1x <- dma
  dm1x$Lifetime_PTSD_Dx <- dm1x$P
  #dm1x <- merge(dm1,dma,by=c("IID_pc"),suffixes=c("_fam","")) #only has IID
  
  dm1x$Case <- dm1x$P
  dm1_exp <- subset(dm1x,!is.na(Case),select=c(FID,IID,Case))
 
  print(table(dm1_exp$Case))
  
  
  write.table(dm1_exp,file=paste("pheno/p2_",tdata,"_mede","_",ancgroup,".pheno",sep=""),quote=F,row.names=F)
    
  dm_cov <- dm1x  

  covexp4 <- subset(dm_cov,select=c("FID","IID","C1","C2","C3","C4","C5")) #mrs doesnt need sex!
  
  write.table(covexp4,file=paste("pheno/p2_",tdata,"_",ancgroup,"_mede","_pcs.cov",sep=""),quote=F,row.names=F)
 }  
  