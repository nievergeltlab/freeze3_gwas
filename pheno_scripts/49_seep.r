library(data.table)
library(janitor)
 unlist_split <- function(x, ...)
  {
	toret <- unlist(strsplit(x, ...) )
	return(t(toret))
  }
  
tdata="psy4"

 dm1 <- fread(paste('data/p2_seep.csv',sep=''),data.table=F,na.strings=c(NA,"-9"))
 
 
 dm1$PTSD_whichmax <- apply(dm1[,c("Current_PTSD_Dx","Lifetime_PTSD_Dx")],1,which.max)

 takecol <- function(x)
 {
  column_to_pick <- x[3]
  return(x[column_to_pick])
  print(x[column_to_pick])
 }
 dm1$PTSD_dx <- apply(dm1[,c("Current_PTSD_Dx","Lifetime_PTSD_Dx","PTSD_whichmax")],1,takecol)


dm1 <- remove_empty(dm1, which = c("cols"))

#There are people with ONLY current measure but no lifetime and (vice versa!!). I will set them to being PTSD people
#dont uise it, not as predictive as case status
# currentcases <- which( is.na(dm1$Lifetime_PTSD_Continuous) & !is.na(dm1$Current_PTSD_Continuous) )
# ltcases <- which( !is.na(dm1$Lifetime_PTSD_Continuous) & is.na(dm1$Current_PTSD_Continuous) )
# controls <- which(
                 # (dm1$Lifetime_PTSD_Dx == 1 & is.na(dm1$Current_PTSD_Dx)) | 
                 # ( is.na(dm1$Lifetime_PTSD_Dx) & dm1$Current_PTSD_Dx == 1  ) |
                 # (dm1$Lifetime_PTSD_Dx == 1 & dm1$Current_PTSD_Dx == 1)    
                 # )
# dm1$PTSD_Continuous <- NA
# dm1[currentcases,]$PTSD_Continuous   <-  dm1[currentcases,]$Current_PTSD_Continuous           
# dm1[ltcases,]$PTSD_Continuous   <-  dm1[ltcases,]$Lifetime_PTSD_Continuous           
# dm1[controls,]$PTSD_Continuous   <-  0    
                 



 fam <- fread(paste('/home/pgcdac/DWFV2CJb8Piv_0116_pgc_data/pts/wave3/v1/',tdata,'/qc/pts_',tdata,'_mix_am-qc.fam',sep=''),data.table=F)
 names(fam) <- c("FID","IID","M","F","G","P")

 dm1$Case <- dm1$PTSD_dx

 
 for (ancgroup in c("eur"))
 {
  pcs <-  read.table(paste('/home/pgcdac/DWFV2CJb8Piv_0116_pgc_data/pts/wave3/v1/',tdata,'/covariates/pts_',tdata,'_mix_am-qc-',ancgroup,'_pca.menv.mds_cov',sep=''),header=T,stringsAsFactors=F)
  
  dma <- merge(pcs,fam,by=c("FID","IID"),suffixes=c("_pc","")) 
  dm1x <- merge(dm1,dma,by=c("FID","IID"),suffixes=c("_fam","")) 
  
  
  dm1_exp <- subset(dm1x,!is.na(Case),select=c(FID,IID,Case))
 
  print(table(dm1_exp$Case))
  
  
  write.table(dm1_exp,file=paste("pheno/p2_",tdata,"_psy4","_",ancgroup,".pheno",sep=""),quote=F,row.names=F)
    
  dm_cov <- dm1x  

  covexp4 <- subset(dm_cov,select=c("FID","IID","C1","C2","C3","C4","C5")) #mrs doesnt need sex!
  
  write.table(covexp4,file=paste("pheno/p2_",tdata,"_",ancgroup,"_psy4","_pcs.cov",sep=""),quote=F,row.names=F)
 }  
  