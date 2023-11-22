

get_cox <- function(df,olinkrange,outpath){
    library(tidyverse)
    library("survival")
    library("survminer")
    olink <- names(df)[olinkrange]
    
    ###原始的status中 1=患病；2=没有患病，为了运行cox回归需要设置1=没有患病，2=患病
    
    df <- df %>%
    mutate(allcause = case_when(
    allcause == 1 ~ 2,
    allcause == 2 ~ 1,
    TRUE ~ allcause  # 如果不满足以上条件，保持不变
    ),
    AD = case_when(
    AD == 1 ~ 2,
    AD == 2 ~ 1,
    TRUE ~ AD  # 如果不满足以上条件，保持不变
    ),
    VD = case_when(
    VD == 1 ~ 2,
    VD == 2 ~ 1,
    TRUE ~ VD  # 如果不满足以上条件，保持不变
    )
    )
    
    ##此时1=没有患病，2=患病
    allcause <- df  
    AD <- df %>% dplyr::filter(allcause == 1 | (VD == 1 & AD == 2)) 
    VD <- df %>% dplyr::filter(allcause == 1 | (AD == 1 & VD == 2)) 
    
    ######{model1: age, sex, education level, and APOE ε4 status}
    allcauseout1 <- data.frame()
    for (i in olink) {
      expr = allcause[,i]
      cox = coxph(Surv(time_year, allcause) ~ expr+
                    age+sex+edu_g+apoe_e4carrier
                    ,allcause)
      coxsummary = summary(cox)
      allcauseout1=rbind(allcauseout1,cbind(protein=i,
                              HR_allcause=round(coxsummary$coefficients[1,"exp(coef)"],2), 
                              z_allcause=round(coxsummary$coefficients[1,"z"],2), 
                              "95% CI_allcause"=paste(round(coxsummary$conf.int[1,3], 2), round(coxsummary$conf.int[1,4],2), sep = "-"),
                              pvalue_allcause=coxsummary$coefficients[1,"Pr(>|z|)"])) 
      }
    allcauseout1[,"fdr_allcause"] <- p.adjust(allcauseout1$pvalue_allcause,method = "fdr")
    allcauseout1[,"bonf_allcause"] <- p.adjust(allcauseout1$pvalue_allcause,method = "bonferroni")

  ##AD
  ADout1 <- data.frame()
  for (i in olink) {
      expr = AD[,i]
      cox = coxph(Surv(time_year, AD) ~ expr+
                  age+sex+edu_g+apoe_e4carrier
                  ,AD)
      coxsummary = summary(cox)
      ADout1=rbind(ADout1,cbind(protein=i,
                            HR_AD=round(coxsummary$coefficients[1,"exp(coef)"],2), 
                            z_AD=round(coxsummary$coefficients[1,"z"],2), 
                            "95% CI_AD"=paste(round(coxsummary$conf.int[1,3], 2), round(coxsummary$conf.int[1,4],2), sep = "-"),
                            pvalue_AD=coxsummary$coefficients[1,"Pr(>|z|)"])) 
    }
  ADout1[,"fdr_AD"] <- p.adjust(ADout1$pvalue_AD,method = "fdr")
  ADout1[,"bonf_AD"] <- p.adjust(ADout1$pvalue_AD,method = "bonferroni")
  
  ##VD
  VDout1 <- data.frame()
  for (i in olink) {
    expr = VD[,i]
    cox = coxph(Surv(time_year, VD) ~ expr+
                age+sex+edu_g+apoe_e4carrier
                ,VD)
    coxsummary = summary(cox)
    VDout1=rbind(VDout1,cbind(protein=i,
                          HR_VD=round(coxsummary$coefficients[1,"exp(coef)"],2), 
                          z_VD=round(coxsummary$coefficients[1,"z"],2), 
                          "95% CI_VD"=paste(round(coxsummary$conf.int[1,3], 2), round(coxsummary$conf.int[1,4],2), sep = "-"),
                          pvalue_VD=coxsummary$coefficients[1,"Pr(>|z|)"])) 
    }
  VDout1[,"fdr_VD"] <- p.adjust(VDout1$pvalue_VD,method = "fdr")
  VDout1[,"bonf_VD"] <- p.adjust(VDout1$pvalue_VD,method = "bonferroni")
    
    ######{model2: model1 + ethnicity, TDI, BMI, smoking, eGFR, hypertension (yes/no), diabetes mellitus (yes/no), and CVD (yes/no)}
    ##allcause
    allcauseout2 <- data.frame()
    for (i in olink) {
        expr = allcause[,i]
        cox = coxph(Surv(time_year, allcause) ~ expr+
                    age+TDI+sex+BMI+smoking+drinking+ethnicity_g+
                    edu_g+diabetes+hypertension+eGFRcr+apoe_e4carrier+CVD
                    ,allcause)
        coxsummary = summary(cox)
        allcauseout2=rbind(allcauseout2,cbind(protein=i,
                              HR_allcause=round(coxsummary$coefficients[1,"exp(coef)"],2), 
                              z_allcause=round(coxsummary$coefficients[1,"z"],2), 
                              "95% CI_allcause"=paste(round(coxsummary$conf.int[1,3], 2), round(coxsummary$conf.int[1,4],2), sep = "-"),
                              pvalue_allcause=coxsummary$coefficients[1,"Pr(>|z|)"])) 
      }
    allcauseout2[,"fdr_allcause"] <- p.adjust(allcauseout2$pvalue_allcause,method = "fdr")
    allcauseout2[,"bonf_allcause"] <- p.adjust(allcauseout2$pvalue_allcause,method = "bonferroni")

    ##AD
    ADout2 <- data.frame()
    for (i in olink) {
        expr = AD[,i]
        cox = coxph(Surv(time_year, AD) ~ expr+
                    age+TDI+sex+BMI+smoking+drinking+ethnicity_g+
                    edu_g+diabetes+hypertension+eGFRcr+apoe_e4carrier+CVD
                    ,AD)
      coxsummary = summary(cox)
      ADout2=rbind(ADout2,cbind(protein=i,
                            HR_AD=round(coxsummary$coefficients[1,"exp(coef)"],2), 
                            z_AD=round(coxsummary$coefficients[1,"z"],2), 
                            "95% CI_AD"=paste(round(coxsummary$conf.int[1,3], 2), round(coxsummary$conf.int[1,4],2), sep = "-"),
                            pvalue_AD=coxsummary$coefficients[1,"Pr(>|z|)"])) 
      }
    ADout2[,"fdr_AD"] <- p.adjust(ADout2$pvalue_AD,method = "fdr")
    ADout2[,"bonf_AD"] <- p.adjust(ADout2$pvalue_AD,method = "bonferroni")

    ##VD
    VDout2 <- data.frame()
    for (i in olink) {
      expr = VD[,i]
      cox = coxph(Surv(time_year, VD) ~ expr+
                  age+TDI+sex+BMI+smoking+drinking+ethnicity_g+
                  edu_g+diabetes+hypertension+eGFRcr+apoe_e4carrier+CVD
                  ,VD)
      coxsummary = summary(cox)
      VDout2=rbind(VDout2,cbind(protein=i,
                          HR_VD=round(coxsummary$coefficients[1,"exp(coef)"],2), 
                          z_VD=round(coxsummary$coefficients[1,"z"],2), 
                          "95% CI_VD"=paste(round(coxsummary$conf.int[1,3], 2), round(coxsummary$conf.int[1,4],2), sep = "-"),
                          pvalue_VD=coxsummary$coefficients[1,"Pr(>|z|)"])) 
      }
    VDout2[,"fdr_VD"] <- p.adjust(VDout2$pvalue_VD,method = "fdr")
    VDout2[,"bonf_VD"] <- p.adjust(VDout2$pvalue_VD,method = "bonferroni")
    
    
    
    ######{model3: model2 + baseline cognition(fluid intelligence)}
    allcause <- df %>% dplyr::filter(!is.na(fluid_intelligence))
    AD <- allcause %>% dplyr::filter(allcause == 1 | (VD == 1 & AD == 2)) 
    VD <- allcause %>% dplyr::filter(allcause == 1 | (AD == 1 & VD == 2)) 
    ##allcause
    allcauseout3 <- data.frame()
    for (i in olink) {
        expr = allcause[,i]
        cox = coxph(Surv(time_year, allcause) ~ expr+
                      age+TDI+sex+BMI+smoking+drinking+ethnicity_g+
                      edu_g+diabetes+hypertension+eGFRcr+apoe_e4carrier+CVD+fluid_intelligence
                      ,allcause)
        coxsummary = summary(cox)
        allcauseout3=rbind(allcauseout3,cbind(protein=i,
                              HR_allcause=round(coxsummary$coefficients[1,"exp(coef)"],2), 
                              z_allcause=round(coxsummary$coefficients[1,"z"],2), 
                              "95% CI_allcause"=paste(round(coxsummary$conf.int[1,3], 2), round(coxsummary$conf.int[1,4],2), sep = "-"),
                              pvalue_allcause=coxsummary$coefficients[1,"Pr(>|z|)"])) 
      }
    allcauseout3[,"fdr_allcause"] <- p.adjust(allcauseout3$pvalue_allcause,method = "fdr")
    allcauseout3[,"bonf_allcause"] <- p.adjust(allcauseout3$pvalue_allcause,method = "bonferroni")

    ##AD
    ADout3 <- data.frame()
    for (i in olink) {
        expr = AD[,i]
        cox = coxph(Surv(time_year, AD) ~ expr+
                      age+TDI+sex+BMI+smoking+drinking+ethnicity_g+
                      edu_g+diabetes+hypertension+eGFRcr+apoe_e4carrier+CVD+fluid_intelligence
                      ,AD)
        coxsummary = summary(cox)
        ADout3=rbind(ADout3,cbind(protein=i,
                            HR_AD=round(coxsummary$coefficients[1,"exp(coef)"],2), 
                            z_AD=round(coxsummary$coefficients[1,"z"],2), 
                            "95% CI_AD"=paste(round(coxsummary$conf.int[1,3], 2), round(coxsummary$conf.int[1,4],2), sep = "-"),
                            pvalue_AD=coxsummary$coefficients[1,"Pr(>|z|)"])) 
      }
    ADout3[,"fdr_AD"] <- p.adjust(ADout3$pvalue_AD,method = "fdr")
    ADout3[,"bonf_AD"] <- p.adjust(ADout3$pvalue_AD,method = "bonferroni")

    ##VD
    VDout3 <- data.frame()
    for (i in olink) {
      expr = VD[,i]
      cox = coxph(Surv(time_year, VD) ~ expr+
                      age+TDI+sex+BMI+smoking+drinking+ethnicity_g+
                      edu_g+diabetes+hypertension+eGFRcr+apoe_e4carrier+CVD+fluid_intelligence
                      ,VD)
      coxsummary = summary(cox)
      VDout3=rbind(VDout3,cbind(protein=i,
                          HR_VD=round(coxsummary$coefficients[1,"exp(coef)"],2), 
                          z_VD=round(coxsummary$coefficients[1,"z"],2), 
                          "95% CI_VD"=paste(round(coxsummary$conf.int[1,3], 2), round(coxsummary$conf.int[1,4],2), sep = "-"),
                          pvalue_VD=coxsummary$coefficients[1,"Pr(>|z|)"])) 
      }
    VDout3[,"fdr_VD"] <- p.adjust(VDout3$pvalue_VD,method = "fdr")
    VDout3[,"bonf_VD"] <- p.adjust(VDout3$pvalue_VD,method = "bonferroni")

    #######{collect result and output}
    
    outTab1 <- allcauseout1 %>%  left_join(.,ADout1,by=c("protein"="protein")) %>% 
      left_join(.,VDout1,by=c("protein"="protein"))
  
    outTab2 <- allcauseout2 %>%  left_join(.,ADout2,by=c("protein"="protein")) %>% 
      left_join(.,VDout2,by=c("protein"="protein"))
    
    outTab3 <- allcauseout3 %>%  left_join(.,ADout3,by=c("protein"="protein")) %>% 
      left_join(.,VDout3,by=c("protein"="protein"))
    
    writexl::write_xlsx(x = list("model1"=outTab1,"model2"=outTab2,"model3"=outTab3),path=outpath)
    
  }
  
  
  