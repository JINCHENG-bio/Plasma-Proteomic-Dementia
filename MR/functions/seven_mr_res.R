 seven_mr_res <- function(dat){
         ###MR
        library(TwoSampleMR)
         res  <- TwoSampleMR::mr(dat)
        # check the least SNP counts
        #print(paste0(id,"_SNP数_",res$nsnp[1]))
        
        ##增加Contamination mixture和cML-MA-BIC方法    
        ##Contamination mixture：（1）robust；（2）异构数据中查找因果机制；（3）找出有效的IV和无效的IV
        library(MRcML)
        library(MendelianRandomization)
        dat <- subset(dat,mr_keep)    ###二次确保满足mr_keep的snp保留
        
        MRinputObject <- mr_input(bx = dat$beta.exposure,
                              bxse = dat$se.exposure,
                              by = dat$beta.outcome,
                              byse = dat$se.outcome)
                              
                              
         conmix <- tryCatch(
           mr_conmix(MRinputObject, CIStep = 0.01),
           error = function(e) {
           mr_conmix(MRinputObject, CIStep = -0.01)
          }
        )
        
        CM <- data.frame(id.exposure=c(dat$id.exposure[1]),
                 id.outcome=c(dat$id.outcome[1]),
                 outcome=c(dat$outcome[1]),
                 exposure=c(dat$exposure[1]),
                 method=c("Contamination mixture"),
                 nsnp=c(conmix@SNPs),
                 b=c(conmix@Estimate),
                 se=c(0),
                 low = c(min(conmix@CIRange)),
                 up = c(max(conmix@CIRange)),
                 pval=c(conmix@Pvalue),
                 nvalidsnp=c(length(conmix@Valid)),
                 ValidSNP=c(paste0(conmix@Valid,collapse = ","))
                 )

        ###cML-MA-BIC：不依赖于InSIDE假设，可以控制相关和不相关的多效性
        cML_result = MRcML::mr_cML(dat$beta.exposure,
                            dat$beta.outcome,
                            dat$se.exposure,
                            dat$se.outcome,
                            n =dat$samplesize.exposure[1],
                            random_start = 100,
                            random_seed = 1)

        cML_MA_BIC <- data.frame(id.exposure=c(dat$id.exposure[1]),
                   id.outcome=c(dat$id.outcome[1]),
                   outcome=c(dat$outcome[1]),
                   exposure=c(dat$exposure[1]),
                   method=c("cML-MA-BIC"),
                   nsnp=c(nrow(dat)-length(cML_result$BIC_invalid)),
                   b=c(cML_result$MA_BIC_theta[1]),
                   se=c(cML_result$MA_BIC_se[1]),
                   pval=c(cML_result$MA_BIC_p[1]),
                   snps=c(ifelse(length(cML_result$BIC_invalid) ==0,"", cML_result$BIC_invalid))
                    )#cML_MA_BIC

        res <- bind_rows(res,CM[,c(1:8,11)],cML_MA_BIC[,1:9])
        seven_mr_res <- list(res,CM,cML_MA_BIC)
        return(seven_mr_res)
    }