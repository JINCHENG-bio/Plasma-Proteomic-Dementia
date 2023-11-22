library("tidyverse")
library("data.table")
library("stringr")
library("openxlsx")
library("grDevices")
library(readxl)
library(dplyr)
library(MRcML)
library(MendelianRandomization)
library(TwoSampleMR)
library("doParallel")      #����doParallel������֮��ע�����
library("foreach")         #����foreach��


file_list <- list.files("/mnt/data/lijincheng/UKB/Olink/0627/UKB_74pGWAS/74pGWAS/", pattern = "\\.rds$", full.names = F)
file_list <- file_list[file_list != "annotation.rds"]
annotation <- readRDS("/mnt/data/lijincheng/UKB/Olink/0627/UKB_74pGWAS/74pGWAS/annotation.rds")
out <- "/mnt/data/lijincheng/UKB/Olink/0627/final_MR1009/results/"
setwd(out)

pGWAS <- function(i){
  id <- sub("(.*)\\.rds", "\\1", file_list[i])
  out_each <- paste0(out,id,"/")
  if (   (!dir.exists(out_each) | length(dir(out_each)) == 0 |  length(list.files(out_each, pattern = "\\.xlsx$", full.names = TRUE)) == 0 )   ){ ##��ѡ��ĵ��ף������ļ��в����ڻ����ļ��д��ڵ���Ϊ��
  dir.create(id)
  df <- readRDS(paste0("/mnt/data/lijincheng/UKB/Olink/0627/UKB_74pGWAS/74pGWAS/",file_list[i])) %>% dplyr::filter(pval < 5e-8)
  df <- left_join(df,annotation,by=c("ID"="ID"))
  df <- df[order(df$pval),]


  blood_exp_dat <- TwoSampleMR::format_data(
          dat = df,
          type = "exposure",
          snps = NULL,
          header = TRUE,
          phenotype_col = "phe",
          snp_col = "rsid",
          beta_col = "BETA",
          se_col = "SE",
          eaf_col = 'A1FREQ',
          effect_allele_col = "ALLELE1",
          other_allele_col = "ALLELE0",
          pval_col = "pval",
          units_col = FALSE,
          ncase_col = FALSE,
          ncontrol_col = FALSE,
          samplesize_col = "N",
          gene_col = FALSE,
          id_col = FALSE,
          min_pval = 1e-200,
          z_col = FALSE,
          info_col = FALSE,
          chr_col = FALSE,
          pos_col = FALSE,
          log_pval = FALSE
    )
  #### makesure all SNP input in the 1000g panel
  #rsid <- ieugwasr::afl2_rsid(blood_exp_dat$SNP, reference = "1000g")
  #save(rsid,file="/mnt/data/lijincheng/UKB/Olink/0627/final_MR/UVMR/1000gSNP.Rdata")
  #blood_exp_dat <- blood_exp_dat %>% dplyr::filter(SNP %in% rsid$rsid)
  
    repeat{
      try({
        ###clumb
        blood_exp_dat <- TwoSampleMR::clump_data(dat = blood_exp_dat,
                                                 clump_r2 = 0.001,
                                                 clump_kb = 10000, ##10Mb
                                                 pop = "EUR")
      })
      if(exists("blood_exp_dat")) break
      Sys.sleep(2)
    }
    
    repeat{
      try({
      AD_out_dat <-  TwoSampleMR::extract_outcome_data(snps =blood_exp_dat$SNP, outcomes =   'ieu-b-2', #2013lambert  ##'ieu-b-2', ###kunkle
                                                       maf_threshold = 0.01)
      })
      if(exists("AD_out_dat")) break
      Sys.sleep(2)
    }
                                       
    
     ###-----{2}harmonise-----
    dat <- TwoSampleMR::harmonise_data(exposure_dat = blood_exp_dat,
                          outcome_dat = AD_out_dat)
    
    dat$r.exposure <- TwoSampleMR::get_r_from_bsen(b = dat$beta.exposure,
                                  dat$se.exposure,
                                 dat$samplesize.exposure)
  
  ###�����ۺϿ��Ƕ���MR�����Ľ�����������ﱣ������INSIDE�����SNP
    dat <- subset(dat,mr_keep)    ###������mr_keep�����snp��������Ϊ��һ��������
    
    ###dat$r.outcome Ĭ����get_r_from_pn����
    ### LOADȱ��EAF��Ϣ�����ã�get_r_from_bsen����R2
    #dat$r.outcome <- get_r_from_bsen(b = dat$beta.outcome,
    #                                 dat$se.outcome,
    #                                dat$samplesize.outcome)
  
  ###------{3}��Ҫ��װ�����Ķ�ȡ----- 
    source("/mnt/data/lijincheng/UKB/Olink/0627/final_MR/function/res_MRPRESSO.R")
    source("/mnt/data/lijincheng/UKB/Olink/0627/final_MR/function/seven_mr_res.R")
    ###-----{4}��Ҫ��������------
    ###-----{4.1}MRPRESSO�ж��Ƿ����outlier����������distortionʧ��������������޳���Ⱥֵ�ٽ��к�������------
    ###-----Ϊ��ȷ��������SNPˮƽ�������ԣ�����һ����ͨ��������ֱ������MRPRESSO����֮������Ⱥֵoutlier֮��
    ###�ٽ��к����ķ���

    set.seed(123)
    res_presso_step1 = res_MRPRESSO(dat = dat,NbD = 5000,SignifThreshold = 0.05)
    outlier <- res_presso_step1$Outlier[2][[1]]
        
    library(tidyverse)
    dat_input <- dat
    while(!is.na(outlier[1])){
      ###�ж�ֱ��û��outlier֮������ѭ��
      dat_input <- dat_input %>% filter(!SNP %in% outlier )
      set.seed(i)
      tmp <- res_MRPRESSO(dat = dat_input,NbD = 5000,SignifThreshold = 0.05)
      res_presso_step1$globalP  <- as.character(res_presso_step1$globalP)
      tmp$globalP <- as.character(tmp$globalP)
      
      res_presso_step1$distortionP  <- as.character(res_presso_step1$distortionP)
      #tmp$distortionP <- as.character(tmp$distortionP)
      #res_presso_step1 <- bind_rows(res_presso_step1,tmp)
      outlier <- tmp$Outlier[2][[1]]
    }
  
  
   ###-----{4.2}������Ⱥֵ������½���7��MR����------
    seven_mr_res <- seven_mr_res(dat_input)
    res <- seven_mr_res[[1]]
    CM <- seven_mr_res[[2]]
    cML_MA_BIC <- seven_mr_res[[3]]    
        
    
   ###-----{5}main_result------
    results <- TwoSampleMR::generate_odds_ratios(res)
    results[which(results$method =="Contamination mixture"),c("or_lci95","or_uci95")] <- c(exp(CM[1,"low"]),exp(CM[1,"up"]))
    # OR
    results$estimate <- paste0(
        format(round(results$or, 3), nsmall = 2), " (", 
        format(round(results$or_lci95, 3), nsmall = 2), "-",
        format(round(results$or_uci95, 3), nsmall = 2), ")")
    resdata <- dat_input
    R2 <- (dat_input$r.exposure)^2
    N <- dat_input$samplesize.exposure
    resdata$F_statistic <- (R2/(1- R2))*(N-2)
    # Assumption 1 and 3
    names(resdata)
    Assumption13 <- subset(resdata,mr_keep==TRUE,
                           select = c("SNP","pval.exposure",
                                      "pval.outcome", "F_statistic",
                                      "mr_keep"))
                                      
                                      
    # -----{6}--Sensitive_analysis------
    res_hete <- TwoSampleMR::mr_heterogeneity(dat_input)
    res_plei <- TwoSampleMR::mr_pleiotropy_test(dat_input)
    res_leaveone <- mr_leaveoneout(dat_input)  # 
    # 
    #set.seed(123)
    #res_presso = res_MRPRESSO(dat = dat_input,NbD = 2500,SignifThreshold = 0.05)
    res_presso = res_presso_step1
    ##Steiger fltering
    res_dir  <-  directionality_test(dat_input)  ##dat$r.outcome Ĭ����get_r_from_pn����
    
    #### Isq statistics
    I2.exposure <- TwoSampleMR::Isq(dat_input$beta.exposure,dat_input$se.exposure)
    I2.outcome <- TwoSampleMR::Isq(dat_input$beta.outcome,dat_input$se.outcome)
    res_isq = data.frame(I2.exposure=I2.exposure,I2.outcome=I2.outcome)
    
    # [["MR-PRESSO results"]][["Global Test"]][["Pvalue"]]
    #res_presso[[1]][[2]][[1]][["Pvalue"]]
    #sink(paste0(out_each,id,"_PRESSO.txt"),append=FALSE,split = FALSE) 
    #print(res_presso)
    #sink()
    #print(res_presso)
    
    # ------{7}Merge_results--------
    # Export
    writexl::write_xlsx(x = list(
    "step1_MRPRESSO" = res_presso_step1,
    "main"=results,
    "Assumption13"=Assumption13,
    "pleiotropy"=res_plei,
    "heterogeneity"=res_hete,
    "final_MRPRESSO"=res_presso,
    "leaveone"=res_leaveone,
    "ReverseMR_Steiger_fltering"=res_dir,
    "Isq_statistic" = res_isq,
    "Conmix"=CM,
    "cML_MA_BIC"=cML_MA_BIC,
    "dat_input"=dat_input),
    path =  paste0(out_each,id,".xlsx"))
    
    # -----{8}--Sensitive_analysis_out------
    p1 <- mr_scatter_plot(res, dat_input)
    #p1[[1]]
    pdf(paste0(out_each,id,"_scatter.pdf"))
    print(p1[[1]])
    dev.off()
    
    p11 <- mr_scatter_plot(res[which(res$method != "Simple mode"),], dat_input)
    #p1[[1]]
    pdf(paste0(out_each,id,"_scatter1.pdf"))
    print(p11[[1]])
    dev.off()
    
    res_single <- mr_singlesnp(dat_input)
    p2 <- mr_forest_plot(res_single)
    pdf(paste0(out_each,id,"_forest.pdf"))
    print(p2[[1]])
    dev.off()
    # 
    p3 <- mr_funnel_plot(res_single)
    pdf(paste0(out_each,id,"_funnel.pdf"))
    print(p3[[1]])
    dev.off()
    # 
    res_loo <- mr_leaveoneout(dat_input)
    pdf(paste0(out_each,id,"_leave_one_out.pdf"))
    print(mr_leaveoneout_plot(res_loo))
    dev.off()
    
    
    }
   }
  
  
  
  

system.time({
  cl<- makeCluster(10)      
  registerDoParallel(cl)       #���н���ע��
  mydata1 <- foreach(
              i=1:length(file_list),   .errorhandling="pass",       #����ȴ�����Ĳ���
              .packages = c("tidyverse", "data.table","readxl","dplyr","TwoSampleMR","stringr","openxlsx","grDevices",
                            "MRcML","MendelianRandomization") 
              #������̹�����ϵͳ����
  ) %dopar% pGWAS(i)
  stopCluster(cl)
})













for(i in seq_along(file_list)){
  id <- sub("(.*)\\.rds", "\\1", file_list[i])
  out_each <- paste0(out,id,"/")
  print(paste0("----------------running ", i," ", id ," ! -------------------------------"))
  if (   (!dir.exists(out_each) | length(dir(out_each)) == 0 |  length(list.files(out_each, pattern = "\\.xlsx$", full.names = TRUE)) == 0 )   ){ ##��ѡ��ĵ��ף������ļ��в����ڻ����ļ��д��ڵ���Ϊ��
  dir.create(id)
  df <- readRDS(paste0("/mnt/data/lijincheng/UKB/Olink/0627/UKB_74pGWAS/74pGWAS/",file_list[i])) %>% dplyr::filter(pval < 5e-8)
  df <- left_join(df,annotation,by=c("ID"="ID"))
  df <- df[order(df$pval),]


  blood_exp_dat <- TwoSampleMR::format_data(
          dat = df,
          type = "exposure",
          snps = NULL,
          header = TRUE,
          phenotype_col = "phe",
          snp_col = "rsid",
          beta_col = "BETA",
          se_col = "SE",
          eaf_col = 'A1FREQ',
          effect_allele_col = "ALLELE1",
          other_allele_col = "ALLELE0",
          pval_col = "pval",
          units_col = FALSE,
          ncase_col = FALSE,
          ncontrol_col = FALSE,
          samplesize_col = "N",
          gene_col = FALSE,
          id_col = FALSE,
          min_pval = 1e-200,
          z_col = FALSE,
          info_col = FALSE,
          chr_col = FALSE,
          pos_col = FALSE,
          log_pval = FALSE
    )
  #### makesure all SNP input in the 1000g panel
  #rsid <- ieugwasr::afl2_rsid(blood_exp_dat$SNP, reference = "1000g")
  #save(rsid,file="/mnt/data/lijincheng/UKB/Olink/0627/final_MR/UVMR/1000gSNP.Rdata")
  #blood_exp_dat <- blood_exp_dat %>% dplyr::filter(SNP %in% rsid$rsid)
  
    repeat{
      try({
        ###clumb
        blood_exp_dat <- TwoSampleMR::clump_data(dat = blood_exp_dat,
                                                 clump_r2 = 0.001,
                                                 clump_kb = 10000, ##10Mb
                                                 pop = "EUR")
      })
      if(exists("blood_exp_dat")) break
      Sys.sleep(2)
    }
    
    repeat{
      try({
      AD_out_dat <-  TwoSampleMR::extract_outcome_data(snps =blood_exp_dat$SNP, outcomes =   'ieu-b-2', #2013lambert  ##'ieu-b-2', ###kunkle
                                                       maf_threshold = 0.01)
      })
      if(exists("AD_out_dat")) break
      Sys.sleep(2)
    }
                                       
    
     tryCatch({ # ���Խ���һЩ���ܻ�����Ĳ���
    
     ###-----{2}harmonise-----
    dat <- TwoSampleMR::harmonise_data(exposure_dat = blood_exp_dat,
                          outcome_dat = AD_out_dat)
    
    dat$r.exposure <- TwoSampleMR::get_r_from_bsen(b = dat$beta.exposure,
                                  dat$se.exposure,
                                 dat$samplesize.exposure)
  
  ###�����ۺϿ��Ƕ���MR�����Ľ�����������ﱣ������INSIDE�����SNP
    dat <- subset(dat,mr_keep)    ###������mr_keep�����snp��������Ϊ��һ��������
    
    ###dat$r.outcome Ĭ����get_r_from_pn����
    ### LOADȱ��EAF��Ϣ�����ã�get_r_from_bsen����R2
    #dat$r.outcome <- get_r_from_bsen(b = dat$beta.outcome,
    #                                 dat$se.outcome,
    #                                dat$samplesize.outcome)
  
  ###------{3}��Ҫ��װ�����Ķ�ȡ----- 
    source("/mnt/data/lijincheng/UKB/Olink/0627/final_MR/function/res_MRPRESSO.R")
    source("/mnt/data/lijincheng/UKB/Olink/0627/final_MR/function/seven_mr_res.R")
    ###-----{4}��Ҫ��������------
    ###-----{4.1}MRPRESSO�ж��Ƿ����outlier����������distortionʧ��������������޳���Ⱥֵ�ٽ��к�������------
    ###-----Ϊ��ȷ��������SNPˮƽ�������ԣ�����һ����ͨ��������ֱ������MRPRESSO����֮������Ⱥֵoutlier֮��
    ###�ٽ��к����ķ���

    set.seed(123)
    res_presso_step1 = res_MRPRESSO(dat = dat,NbD = 5000,SignifThreshold = 0.05)
    outlier <- res_presso_step1$Outlier[2][[1]]
        
    library(tidyverse)
    dat_input <- dat
    while(!is.na(outlier[1])){
      ###�ж�ֱ��û��outlier֮������ѭ��
      dat_input <- dat_input %>% filter(!SNP %in% outlier )
      set.seed(i)
      tmp <- res_MRPRESSO(dat = dat_input,NbD = 5000,SignifThreshold = 0.05)
      res_presso_step1$globalP  <- as.character(res_presso_step1$globalP)
      tmp$globalP <- as.character(tmp$globalP)
      
      res_presso_step1$distortionP  <- as.character(res_presso_step1$distortionP)
      #tmp$distortionP <- as.character(tmp$distortionP)
      #res_presso_step1 <- bind_rows(res_presso_step1,tmp)
      outlier <- tmp$Outlier[2][[1]]
    }
  
  
   ###-----{4.2}������Ⱥֵ������½���7��MR����------
    seven_mr_res <- seven_mr_res(dat_input)
    res <- seven_mr_res[[1]]
    CM <- seven_mr_res[[2]]
    cML_MA_BIC <- seven_mr_res[[3]]    
        
    
   ###-----{5}main_result------
    results <- TwoSampleMR::generate_odds_ratios(res)
    results[which(results$method =="Contamination mixture"),c("or_lci95","or_uci95")] <- c(exp(CM[1,"low"]),exp(CM[1,"up"]))
    # OR
    results$estimate <- paste0(
        format(round(results$or, 3), nsmall = 2), " (", 
        format(round(results$or_lci95, 3), nsmall = 2), "-",
        format(round(results$or_uci95, 3), nsmall = 2), ")")
    resdata <- dat_input
    R2 <- (dat_input$r.exposure)^2
    N <- dat_input$samplesize.exposure
    resdata$F_statistic <- (R2/(1- R2))*(N-2)
    # Assumption 1 and 3
    names(resdata)
    Assumption13 <- subset(resdata,mr_keep==TRUE,
                           select = c("SNP","pval.exposure",
                                      "pval.outcome", "F_statistic",
                                      "mr_keep"))
                                      
                                      
    # -----{6}--Sensitive_analysis------
    res_hete <- TwoSampleMR::mr_heterogeneity(dat_input)
    res_plei <- TwoSampleMR::mr_pleiotropy_test(dat_input)
    res_leaveone <- mr_leaveoneout(dat_input)  # 
    # 
    #set.seed(123)
    #res_presso = res_MRPRESSO(dat = dat_input,NbD = 2500,SignifThreshold = 0.05)
    res_presso = res_presso_step1
    ##Steiger fltering
    res_dir  <-  directionality_test(dat_input)  ##dat$r.outcome Ĭ����get_r_from_pn����
    
    
    #### Isq statistics
    I2.exposure <- TwoSampleMR::Isq(dat_input$beta.exposure,dat_input$se.exposure)
    I2.outcome <- TwoSampleMR::Isq(dat_input$beta.outcome,dat_input$se.outcome)
    res_isq = data.frame(I2.exposure=I2.exposure,I2.outcome=I2.outcome)
    
    # [["MR-PRESSO results"]][["Global Test"]][["Pvalue"]]
    #res_presso[[1]][[2]][[1]][["Pvalue"]]
    #sink(paste0(out_each,id,"_PRESSO.txt"),append=FALSE,split = FALSE) 
    #print(res_presso)
    #sink()
    #print(res_presso)
    
    # ------{7}Merge_results--------
    # Export
    writexl::write_xlsx(x = list(
    "step1_MRPRESSO" = res_presso_step1,
    "main"=results,
    "Assumption13"=Assumption13,
    "pleiotropy"=res_plei,
    "heterogeneity"=res_hete,
    "final_MRPRESSO"=res_presso,
    "leaveone"=res_leaveone,
    "ReverseMR_Steiger_fltering"=res_dir,
    "Isq_statistic" = res_isq,
    "Conmix"=CM,
    "cML_MA_BIC"=cML_MA_BIC,
    "dat_input"=dat_input),
    path =  paste0(out_each,id,".xlsx"))
    
    # -----{8}--Sensitive_analysis_out------
    p1 <- mr_scatter_plot(res, dat_input)
    #p1[[1]]
    pdf(paste0(out_each,id,"_scatter.pdf"))
    print(p1[[1]])
    dev.off()
    
    p11 <- mr_scatter_plot(res[which(res$method != "Simple mode"),], dat_input)
    #p1[[1]]
    pdf(paste0(out_each,id,"_scatter1.pdf"))
    print(p11[[1]])
    dev.off()
    
    res_single <- mr_singlesnp(dat_input)
    p2 <- mr_forest_plot(res_single)
    pdf(paste0(out_each,id,"_forest.pdf"))
    print(p2[[1]])
    dev.off()
    # 
    p3 <- mr_funnel_plot(res_single)
    pdf(paste0(out_each,id,"_funnel.pdf"))
    print(p3[[1]])
    dev.off()
    # 
    res_loo <- mr_leaveoneout(dat_input)
    pdf(paste0(out_each,id,"_leave_one_out.pdf"))
    print(mr_leaveoneout_plot(res_loo))
    dev.off()
    
    print(paste0("----------------finished ", i," ", id ," ! -------------------------------"))
    
        }, error = function(e) {
    # ��������ʱ�Ĵ���
    cat("Error in iteration", i, ": ", id, " get wrong", "\n")
  })
    
    }
   }












