 seven_mr_res <- function(dat){
         ###MR
        library(TwoSampleMR)
         res  <- TwoSampleMR::mr(dat)
        # check the least SNP counts
        #print(paste0(id,"_SNP��_",res$nsnp[1]))
        
        ##����Contamination mixture��cML-MA-BIC����    
        ##Contamination mixture����1��robust����2���칹�����в���������ƣ���3���ҳ���Ч��IV����Ч��IV
        library(MRcML)
        library(MendelianRandomization)
        dat <- subset(dat,mr_keep)    ###����ȷ������mr_keep��snp����
        
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

        ###cML-MA-BIC����������InSIDE���裬���Կ�����غͲ���صĶ�Ч��
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