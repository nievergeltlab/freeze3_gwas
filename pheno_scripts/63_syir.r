library(data.table)
library(janitor)
 unlist_split <- function(x, ...)
  {
	toret <- unlist(strsplit(x, ...) )
	return(t(toret))
  }
  
tdata="syir"
 dm1 <- fread(paste('data/P3_SYIR_63_final.csv',sep=''),data.table=F,na.strings=c(NA,"-9"))

 dm1 <- remove_empty(dm1, which = c("cols"))

 fam <- fread(paste('/home/pgcdac/DWFV2CJb8Piv_0116_pgc_data/pts/wave3/v1/',tdata,'/qc/pts_',tdata,'_mix_am-qc1.fam',sep=''),data.table=F)
 names(fam) <- c("FID","IID","M","F","G","P")

 dm1$Current_PTSD_Continuous <- dm1$cj_memories+ dm1$cj_dreams	+dm1$cj_reliving	+ dm1$cj_upset	+dm1$cj_reaction+	dm1$cj_thinking	+
                                dm1$cj_activities	+dm1$cj_remembering+	dm1$cj_lossinterest+	dm1$cj_distant	+dm1$cj_emotionally	+
                                dm1$cj_future+	dm1$cj_asleep	+dm1$cj_irritable+	dm1$cj_concentrating + dm1$cj_alert + 	dm1$cj_startled

 dm2 <- dm1[!duplicated(dm1$IID),]
 
 for (ancgroup in c("oth"))
 {
  pcs <-  read.table(paste('/home/pgcdac/DWFV2CJb8Piv_0116_pgc_data/pts/wave3/v1/',tdata,'/covariates/pts_',tdata,'_mix_am-qc-',ancgroup,'_pca.menv.mds_cov',sep=''),header=T,stringsAsFactors=F)
  
  dma <- merge(pcs,fam,by=c("FID","IID"),suffixes=c("_pc","")) 
  dm1x <- merge(dm2,dma,by=c("IID"),suffixes=c("_fam","")) 
  
  
  dm1_exp <- subset(dm1x,!is.na(Current_PTSD_Continuous),select=c(FID,IID,Current_PTSD_Continuous))
 
  print(table(is.na(dm1_exp$Current_PTSD_Continuous)))
  
  
  write.table(dm1_exp,file=paste("pheno/p2_",tdata,"_syir","_",ancgroup,".pheno",sep=""),quote=F,row.names=F)
    
  dm_cov <- dm1x  

  covexp4 <- subset(dm_cov,select=c("FID","IID","C1","C2","C3","C4","C5")) #mrs doesnt need sex!
  
  write.table(covexp4,file=paste("pheno/p2_",tdata,"_",ancgroup,"_syir","_pcs.cov",sep=""),quote=F,row.names=F)
 }  
  