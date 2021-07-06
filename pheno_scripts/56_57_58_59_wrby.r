library(data.table)
library(janitor)
library(plyr)
 unlist_split <- function(x, ...)
  {
	toret <- unlist(strsplit(x, ...) )
	return(t(toret))
  }
  
 tdata="wrby"

 dm1a <- fread(paste('data/p2_wang_updated_20180910.csv',sep=''),data.table=F,na.strings=c(NA,"-9"))
 dm1b <- fread(paste('data/p2_yehu.csv',sep=''),data.table=F,na.strings=c(NA,"-9"))
 dm1c <- fread(paste('data/p2_vris_update_2019313.csv',sep=''),data.table=F,na.strings=c(NA,"-9")) 
 exposed <- which(dm1c$LT_Count >=1)
 unexposed <- which(dm1c$LT_Count ==0)
 dm1c$Exposure <- NA
 dm1c[exposed,]$Exposure <- 1
 dm1c[unexposed,]$Exposure <- 0
 
 dm1d <- fread(paste('data/baker_geno.csv',sep=''),data.table=F,na.strings=c(NA,"-9")) 
 exposed <- which(dm1d$LT_Count >=1)
 unexposed <- which(dm1d$LT_Count ==0)
 dm1d$Exposure <- NA
 dm1d[exposed,]$Exposure <- 1
 dm1d[unexposed,]$Exposure <- 0
 
  famx <- fread(paste('/home/pgcdac/DWFV2CJb8Piv_0116_pgc_data/pts/wave3/v1/',tdata,'/qc/pts_',tdata,'_mix_am-qc.fam',sep=''),data.table=F) #some stupid merging stuff I have to do, the original file does not have FID
 names(famx) <- c("FID","IID","M","F","G","P")
 famx <- famx[,1:2]
 dm1d <- merge(dm1d,famx,by=c("IID"))
 
 # dm1d1$LT_Count <- dm1d1$LT_Count_T0
 # #Need to just manually calculate best PTSD score for the cohort 1-4, 11-3 subjects..
 # mrs1caps <- fread("data/CAPS_critA_revised_9-17-12.csv",data.table=F)[,1:20]
 # mrs1capsw <- reshape(mrs1caps,timevar='visit',idvar='studyid',sep="_V",direction="wide")

 # mrs1capsw$Current_PTSD_Continuous <- apply(mrs1capsw[,c("CAPStots_V0","CAPStots_V2","CAPStots_V3")],1,which.max)
 
 
 # mrs2caps <- fread('data/MRSII_CAPS_FINAL_v2.csv',data.table=F)[,1:20]

 # mrs2capsw <- reshape(mrs1caps,timevar='visit',idvar='studyid',sep="_V",direction="wide")
 # mrs2capsw$Current_PTSD_Continuous <- apply(mrs2capsw[,c("CAPStots_V0","CAPStots_V2","CAPStots_V3")],1,which.max)
 
 
  
 
 # takecol <- function(x)
 # {
  # column_to_pick <- x[4]
  # return(x[column_to_pick])
 # }
 # dm1$Current_PTSD_Continuous <- apply(dm1[,c("Current_PTSD_Continuous_T2","Current_PTSD_Continuous_T3","Current_PTSD_Continuous_T1_recode","Current_PTSD_Continuous_whichmax")],1,takecol)
 # dm1$Current_PTSD_Dx <- apply(dm1[,c("Current_PTSD_Dx_T2","Current_PTSD_Dx_T3","Current_PTSD_Dx_T1","Current_PTSD_Dx_whichmax")],1,takecol)


 #dm1c2 
 dm1a <- remove_empty(dm1a, which = c("cols"))
 dm1b <- remove_empty(dm1b, which = c("cols"))
 dm1c <- remove_empty(dm1c, which = c("cols"))
 dm1d <- remove_empty(dm1d, which = c("cols"))
 
 # dm1a$Current_PTSD_Continuous_harmonized <- dm1a$Current_PTSD_Continuous / 51
 # dm1b$Current_PTSD_Continuous_harmonized <- dm1b$Current_PTSD_Continuous / 136
 # dm1c1$Current_PTSD_Continuous_harmonized <-  dm1c1$Current_PTSD_Continuous / 80 
 # dm1c2$Current_PTSD_Continuous_harmonized <- 0#everyone is a control..., impute to 0
 
# dm1a$LT_Count_harmonized <- dm1a$LT_Count / 10
# dm1b$LT_Count_harmonized <- dm1b$LT_Count / 11
# dm1c1$LT_Count_harmonized <- dm1c1$LT_Count / 20
# dm1c2$LT_Count_harmonized <- dm1c2$LT_Count / 12
   
 dm1a$study <- "wang"
 dm1b$study <- "yehu"
 dm1c$study <- "vris"
 dm1d$study <- "bake"
 
 dm1 <- rbind.fill(dm1a,dm1b,dm1c,dm1d)
 
 
 ##harmonize to 0-1 based on theoretical ranges
 
 #split up teic kaufman etc based on substudy
 
 fam <- fread(paste('/home/pgcdac/DWFV2CJb8Piv_0116_pgc_data/pts/wave3/v1/',tdata,'/qc/pts_',tdata,'_mix_am-qc.fam',sep=''),data.table=F)
 names(fam) <- c("FID","IID","M","F","G","P")
 for (ancgroup in c("eur","aam"))
 {
  pcs <-  read.table(paste('/home/pgcdac/DWFV2CJb8Piv_0116_pgc_data/pts/wave3/v1/',tdata,'/covariates/pts_',tdata,'_mix_am-qc-',ancgroup,'_pca.menv.mds_cov',sep=''),header=T,stringsAsFactors=F)
  
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

  
  write.table(dm1_exp,file=paste("pheno/p2_",tdata,"_wrby","_",ancgroup,".pheno",sep=""),quote=F,row.names=F)
    
  dm_cov <- dm1x  

  covexp4 <- subset(dm_cov,select=c("FID","IID","C1","C2","C3","C4","C5")) #mrs doesnt need sex!
  
  write.table(covexp4,file=paste("pheno/p2_",tdata,"_",ancgroup,"_wrby","_pcs.cov",sep=""),quote=F,row.names=F)
 }  
  