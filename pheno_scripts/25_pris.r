library(data.table)
library(janitor)
 unlist_split <- function(x, ...)
  {
	toret <- unlist(strsplit(x, ...) )
	return(t(toret))
  }
  
tdata="pris"

 dm1 <- fread(paste('data/P3_Item_level_PRISMO.csv',sep=''),data.table=F,na.strings=c(NA,"-9"))
 
 dm1$IID <- dm1$GeneticsID
 
dm1 <- remove_empty(dm1, which = c("cols"))


 wm2 <- function(x,...)
 {
  ff <- which.max(x)
  if(length(ff) > 0)
  {
   return(ff)
  } 
   if(length(ff) == 0)
  {
   return(NA)
  }  
 }
  
 dm1$Current_PTSD_Continuous_whichmax <- apply(dm1[,c("SRIP.TotalA",	"SRIP.TotalB","SRIP.TotalC"	,"SRIP.TotalD"	,"SRIP.TotalE")],1,wm2)

 takecol <- function(x)
 {
  column_to_pick <- as.numeric(x[length(x)])
  if(!is.na(column_to_pick))
  { 
   return(x[column_to_pick])
  } 
   if(is.na(column_to_pick))
  { 
   return(NA)
  } 
   
 }
 
 dm1$Current_PTSD_Continuous <- apply(dm1[,c("SRIP.TotalA",	"SRIP.TotalB","SRIP.TotalC"	,"SRIP.TotalD"	,"SRIP.TotalE","Current_PTSD_Continuous_whichmax")],1,takecol)
 dm1$IID <- dm1$FID
 dm1$FID <- NULL #FID and IID are swapped!


 fam <- fread(paste('/home/pgcdac/DWFV2CJb8Piv_0116_pgc_data/pts/wave3/v1/',tdata,'/qc/pts_',tdata,'_mix_am-qc.fam',sep=''),data.table=F)
 names(fam) <- c("FID","IID","M","F","G","P")
 
 for (ancgroup in c("eur")) # ,"lat","aam"))
 {
  pcs <-  read.table(paste('/home/pgcdac/DWFV2CJb8Piv_0116_pgc_data/pts/wave3/v1/',tdata,'/covariates/pts_',tdata,'_mix_am-qc-',ancgroup,'_pca.menv.mds_cov',sep=''),header=T,stringsAsFactors=F)
  
  dma <- merge(pcs,fam,by=c("FID","IID"),suffixes=c("_pc","")) 
  
  #Only have IID in this file..
  dm1x <- merge(dm1,dma,by=c("IID"),suffixes=c("_fam","")) 
  
  
  dm1_exp <- subset(dm1x,!is.na(Current_PTSD_Continuous),select=c(FID,IID,Current_PTSD_Continuous))
 
  print(table(is.na(dm1_exp$Current_PTSD_Continuous)))
  
  
  write.table(dm1_exp,file=paste("pheno/p2_",tdata,"_pris","_",ancgroup,".pheno",sep=""),quote=F,row.names=F)
    
  dm_cov <- dm1x  

  covexp4 <- subset(dm_cov,select=c("FID","IID","C1","C2","C3","C4","C5")) #mrs doesnt need sex!
  
  write.table(covexp4,file=paste("pheno/p2_",tdata,"_",ancgroup,"_pris","_pcs.cov",sep=""),quote=F,row.names=F)
 }  
   
  