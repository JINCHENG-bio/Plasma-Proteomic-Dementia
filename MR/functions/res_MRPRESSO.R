res_MRPRESSO <- function(dat=dat, NbD=1000,SignifThreshold = 0.05) {
      if(nrow(dat)<4){
        resMRPRESSO <- data.frame(check.names = FALSE,
                              Exposure=c(NA,NA), "MR Analysis"=NA , "Causal Estimate"=NA,
                              Sd=NA, "T-stat"=NA ,"P-value"=NA,"Outlier"=NA,"globalP"=NA
        )
        resMRPRESSO[1,1:3]="dat<4SNP"
        return(resMRPRESSO)
      }
      # 
      res_presso <- run_mr_presso(dat,NbDistribution = NbD,SignifThreshold = SignifThreshold)
      resMRPRESSO=res_presso[[1]][["Main MR results"]]
      if(is.null(res_presso[[1]][[2]][["Outlier Test"]])){
        resMRPRESSO$Outlier[2]=NA
        resMRPRESSO$globalP[2]=res_presso[[1]][[2]][[1]][["Pvalue"]]
        return(resMRPRESSO)
      } 
      aaa=res_presso[[1]][["MR-PRESSO results"]][["Distortion Test"]][["Outliers Indices"]]
      outliers <- dat$SNP[aaa]
      resMRPRESSO$Outlier[2]=list(outliers)
      resMRPRESSO$globalP[2]=res_presso[[1]][[2]][[1]][["Pvalue"]]
      resMRPRESSO$distortionP[2]=res_presso[[1]]$`MR-PRESSO results`$`Distortion Test`$Pvalue
      return(resMRPRESSO)
    }