library(data.table)
library(janitor)
library(plyr)
 unlist_split <- function(x, ...)
  {
	toret <- unlist(strsplit(x, ...) )
	return(t(toret))
  }
  
 tdata="pts1"


 ##harmonize to 0-1 based on theoretical ranges
 
 #split up teic kaufman etc based on substudy
 /home/pgcdac/DWFV2CJb8Piv_0116_pgc_data/pts/wave3/v1/pts1/qc/pts_pts1_mix_am-qc1.fam
 /home/pgcdac/DWFV2CJb8Piv_0116_pgc_data/pts/wave3/v1/pts1/qc/pts_pts1_mix_am-qc1.fam
 fam <- fread(paste('/home/pgcdac/DWFV2CJb8Piv_0116_pgc_data/pts/wave3/v1/',tdata,'/qc/pts_',tdata,'_mix_am-qc1.fam',sep=''),data.table=F)
 names(fam) <- c("FID","IID","M","F","G","P")
 for (ancgroup in c("eur","aam"))
 {
  pcs <-  read.table(paste('/home/pgcdac/DWFV2CJb8Piv_0116_pgc_data/pts/wave3/v1/',tdata,'/covariates/pts_',tdata,'_mix_am-qc1-',ancgroup,'_pca.menv.mds_cov',sep=''),header=T,stringsAsFactors=F)

  dma <- merge(pcs,fam,by=c("FID","IID"),suffixes=c("_pc","")) 
  #dm1x <- merge(dm1,dma,by=c("FID","IID"),suffixes=c("_fam","")) 
  dm1x <- dma
  dm1x$Case <- dm1x$P # Current_PTSD_Dx
  dm1_exp <- subset(dm1x,!is.na(Case),select=c(FID,IID,Case))
 
  print(table(is.na(dm1_exp$Case)))
  
  print(table(is.na(dm1_exp$Case)))
    print(table(dm1_exp$Case))
     counts <- table(dm1_exp$Case)
  4/((1/counts[1]) + (1/counts[2]))

  
  write.table(dm1_exp,file=paste("pheno/p2_",tdata,"_pts1","_",ancgroup,".pheno",sep=""),quote=F,row.names=F)
    
  dm_cov <- dm1x  

  covexp4 <- subset(dm_cov,select=c("FID","IID","C1","C2","C3","C4","C5")) #mrs doesnt need sex!
  
  write.table(covexp4,file=paste("pheno/p2_",tdata,"_",ancgroup,"_pts1","_pcs.cov",sep=""),quote=F,row.names=F)
 }  
  