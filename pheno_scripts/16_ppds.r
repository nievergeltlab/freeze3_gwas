library(data.table)
library(janitor)
 unlist_split <- function(x, ...)
  {
	toret <- unlist(strsplit(x, ...) )
	return(t(toret))
  }
  
tdata="nss1"

 dm1 <- fread(paste('data/p2_ppds_updated_20191019_PCL-IV and PCL-6.csv',sep=''),data.table=F,na.strings=c(NA,"-9"))
 #WTF is with the duplicates??
 #dm1[duplicated(dm1$IID),]
 #No idea why they exist but they seem to just be exact duplicates so lets just toss them...
 dm1 <- dm1[!duplicated(dm1$IID),]
 
 
 dm1 <- remove_empty(dm1, which = c("cols"))
 
 
which.max2 <- function(x)
{
 wm <- which.max(x)

 if(length(wm) > 0)
 {
  return(wm)
 }
 else if (length(wm == 0)){

 return(NA)} else  { return (NA)}
 
 
 }
 
  takecol <- function(x)
 {
  if(is.na(x[4]))
  { return (NA)}
  
  column_to_pick <- x[4]
  return(x[column_to_pick])
 }
 
 
 dm1$Current_PTSD_Continuous_T1_recode <- 2.72*(dm1$Current_PTSD_Continuous_T1 -6) + 17 #rescale measure to pcl range. Revised to account for the fact that 5 should be NA!!
 dm1$Current_PTSD_Continuous_whichmax <- apply(dm1[,c("Current_PTSD_Continuous_T2","Current_PTSD_Continuous_T3","Current_PTSD_Continuous_T1_recode")],1,which.max2)
 dm1$Current_PTSD_Dx_whichmax <- apply(dm1[,c("Current_PTSD_Dx_T2","Current_PTSD_Dx_T3","Current_PTSD_Dx_T1")],1,which.max)
 

 dm1$Current_PTSD_Continuous <- apply(dm1[,c("Current_PTSD_Continuous_T2","Current_PTSD_Continuous_T3","Current_PTSD_Continuous_T1_recode","Current_PTSD_Continuous_whichmax")],1,takecol)
 #dm1$Current_PTSD_Dx <- apply(dm1[,c("Current_PTSD_Dx_T2","Current_PTSD_Dx_T3","Current_PTSD_Dx_T1","Current_PTSD_Dx_whichmax")],1,takecol)

 
 fam <- fread(paste('/home/pgcdac/DWFV2CJb8Piv_0116_pgc_data/pts/wave3/v1/',tdata,'/qc/pts_',tdata,'_mix_am-qc.fam',sep=''),data.table=F)
 names(fam) <- c("FID","IID","M","F","G","P")
 for (ancgroup in c("eur","lat","aam"))
 {
  pcs <-  read.table(paste('/home/pgcdac/DWFV2CJb8Piv_0116_pgc_data/pts/wave3/v1/',tdata,'/covariates/pts_',tdata,'_mix_am-qc-',ancgroup,'_pca.menv.mds_cov',sep=''),header=T,stringsAsFactors=F)
  
  dma <- merge(pcs,fam,by=c("FID","IID"),suffixes=c("_pc","")) 
  dm1x <- merge(dm1,dma,by=c("FID","IID"),suffixes=c("_fam","")) 
  
  
  dm1_exp <- subset(dm1x,!is.na(Current_PTSD_Continuous),select=c(FID,IID,Current_PTSD_Continuous))
 
  table(is.na(dm1_exp$Current_PTSD_Continuous))
  
  write.table(dm1_exp,file=paste("pheno/p2_",tdata,"_ppds","_",ancgroup,".pheno",sep=""),quote=F,row.names=F)
    
  dm_cov <- dm1x  

  covexp4 <- subset(dm_cov,select=c("FID","IID","C1","C2","C3","C4","C5")) #mrs doesnt need sex!
  
  write.table(covexp4,file=paste("pheno/p2_",tdata,"_",ancgroup,"_ppds","_pcs.cov",sep=""),quote=F,row.names=F)
 }  
  
  