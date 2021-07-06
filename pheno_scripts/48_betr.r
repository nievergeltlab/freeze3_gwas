library(data.table)
library(janitor)
library(plyr)

 unlist_split <- function(x, ...)
  {
	toret <- unlist(strsplit(x, ...) )
	return(t(toret))
  }
  
tdata="betr"

 dm1a <- fread(paste('data/p2_betr_control_updated_20190220.csv',sep=''),data.table=F,na.strings=c(NA,"-9"))
 dm1b <- fread(paste('data/p2_betr_updated_20190424.csv',sep=''),data.table=F,na.strings=c(NA,"-9"))
  dm1 <- remove_empty(dm1a, which = c("cols"))
 dm1 <- remove_empty(dm1b, which = c("cols"))

 dm1 <- rbind.fill(dm1a,dm1b)
 

 #dm1a$Current_PTSD_Continuous <- 0 #Set better controls to 0
 

 fam <- fread(paste('/home/pgcdac/DWFV2CJb8Piv_0116_pgc_data/pts/wave3/v1/',tdata,'/qc/pts_',tdata,'_mix_am-qc.fam',sep=''),data.table=F)
 names(fam) <- c("FID","IID","M","F","G","P")
 for (ancgroup in c("eur"))
 {
  pcs <-  read.table(paste('/home/pgcdac/DWFV2CJb8Piv_0116_pgc_data/pts/wave3/v1/',tdata,'/covariates/pts_',tdata,'_mix_am-qc-',ancgroup,'_pca.menv.mds_cov',sep=''),header=T,stringsAsFactors=F)
  
  dma <- merge(pcs,fam,by=c("FID","IID"),suffixes=c("_pc","")) 
  dm1x <- merge(dm1,dma,by=c("FID","IID"),suffixes=c("_fam","")) 
  
   dm1x$Case <- dm1x$P  # Current_PTSD_Dx
 
  dm1_exp <- subset(dm1x,!is.na(Case),select=c(FID,IID,Case))
 
  print(table(is.na(dm1_exp$Case)))
  #Neff 
  counts <- table(dm1_exp$Case)
  4/((1/counts[1]) + (1/counts[2]))

  
  write.table(dm1_exp,file=paste("pheno/p2_",tdata,"_betr","_",ancgroup,".pheno",sep=""),quote=F,row.names=F)
    
  dm_cov <- dm1x  

  covexp4 <- subset(dm_cov,select=c("FID","IID","C1","C2","C3","C4","C5")) #mrs doesnt need sex!
  
  write.table(covexp4,file=paste("pheno/p2_",tdata,"_",ancgroup,"_betr","_pcs.cov",sep=""),quote=F,row.names=F)
 }  
    
  