
#########################--------------{01}imputation 01 covariates------------------ ###################


dat_imp <- dat_input %>% dplyr::select(c(1:6,18,8,9:17,19:21))
library(lattice)
library(MASS)
library(nnet)
library(mice)
library(foreign)
library(VIM)
library(car)

nimp <- 100 * nic(dat_imp) / nrow(dat_imp)   ## 5.565778%, M=6

summary(dat_imp)

###multiple imputation
meth<-c(rep(c(""),2),"pmm","pmm",rep(c(""),6),"pmm","","pmm",
        rep(c(""),7))

###############################################################################
# MI using the  PMM for the continuous variables
################################################################################
imp.parm <- mice(dat_imp,m=6,method=meth,maxit=20,seed=123)

fit=with(imp.parm,glm(allcause~age+sex+TDI+BMI+smoking+drinking+ethnicity+edu+diabetes+
                        hypertension+eGFRcr+eGFRcr_cys+
                        apoe_e4carrier+CVD
                      ,family = "binomial"))
fit
summary(fit)

pdf(file = "/mnt/data/lijincheng/UKB/Olink/0627/impu1.pdf")
densityplot(imp.parm)
dev.off()

pdf(file = "/mnt/data/lijincheng/UKB/Olink/0627/impu2.pdf")
stripplot(imp.parm,  TDI+BMI+eGFRcr+eGFRcr_cys~ .imp, 
          col=c("grey",mdc(2)),pch=c(1,20))
dev.off()

pdf(file = "/mnt/data/lijincheng/UKB/Olink/0627/impu3.pdf")
bwplot(imp.parm)
dev.off()


###according to the outcome of fit choose AIC  minimum imputation dataset

completeData <- complete(imp.parm,1)
completeData$eid <- row.names(dat_imp)
completeData<-completeData[,c(ncol(completeData),1:(ncol(completeData)-1))]
write.csv(completeData,"/mnt/data/lijincheng/UKB/Olink/0627/final_input_covriates.csv",row.names = F)
save(dat_imp,file="/mnt/data/lijincheng/UKB/Olink/0627/covariates_before_inputation.Rdata")

#########################--------------{02}imputation 02 olink------------------ ###################

olink_input  <- olink_input  %>%  filter( eid %in% include_id) %>%
  dplyr::select(-`NPM1;Nucleophosmin`,-`PCOLCE;Procollagen C-endopeptidase enhancer 1`,-`CTSS;Cathepsin S`)
row.names(olink_input ) <- olink_input [,1]
olink_input  <- olink_input[,-c(1,2)]


###KNN?岹????????
library(caret)
set.seed(123)
###
preProcValues <- preProcess(olink_input,
                            method = c("knnImpute"),
                            k = 10,
                            knnSummary = mean)
set.seed(123)
impute_olink <- predict(preProcValues,olink_input)

##keep center and scale data
final_impute_olink <- scale(impute_olink, center = TRUE, scale = TRUE )

#########################--------------{03} calculate time ------------------ ###################

dat_time <- dat %>% dplyr::select(c("f.eid","indate","allcause_d_date","AD_date","VD_date",
                                    "FTD_date","death_date0","lost_date")) %>% 
  dplyr::filter(f.eid %in% include_id) ##48368-1

dat_time <- dat_time %>% mutate(end_date0 = as.IDate(ifelse(is.na(allcause_d_date) & is.na(death_date0) & is.na(lost_date), NA,
                                                            as.IDate(pmin(allcause_d_date,AD_date,VD_date,FTD_date,death_date0,lost_date,na.rm = TRUE)))))

#max(dat_time$end_date0,na.rm = T) ###"2021-11-08"
#as.Date("2021-12-31")

dat_time <- dat_time %>% mutate(end_date_final = as.IDate(ifelse(!is.na(end_date0),end_date0,
                                                                 as.Date("2021-12-31")))) %>% 
  mutate(time_day =as.numeric(difftime(as.Date(end_date_final), as.Date(indate), units = "days")),
         time_week =as.numeric(difftime(as.Date(end_date_final), as.Date(indate), units = "weeks"))) %>% 
  mutate(time_year = round(time_day/365,2))

dat_time<- as.data.frame(dat_time)
row.names(dat_time) <- dat_time[,1]

save(dat_time,file = "/mnt/data/lijincheng/UKB/Olink/0627/dat_time.Rdata")

final_input <- cbind(cbind(completeData,dat_time[,c(2,11)]) ,final_impute_olink)

write.csv(final_input ,"/mnt/data/lijincheng/UKB/Olink/0627/01final_input.csv",row.names=F)

save.image(file="/mnt/data/lijincheng/UKB/Olink/0627/01final_all_imputation.Rdata")



#########################--------------{03}cox regression------------------ ###################
###################--------------------{01 collect data}-----------------------------------

protein <- colnames(final_input)[24:1483]
split_vector <- strsplit(protein, ";") 
result <- sapply(split_vector, function(x) x[1])
colnames(final_input)[24:1483] <- result
protein_reflection <- data.frame(fullname=protein ,shortname=result)


load(file="/mnt/data/lijincheng/UKB/Olink/0627/cognitive.Rdata")
cognitive$f.eid <- as.character(cognitive$f.eid)
total_input <- final_input %>% mutate(within_two_year = ifelse(time_day < 730, 1,0)) %>%
  left_join(.,cognitive,by=c("eid" = "f.eid")) %>%
  dplyr::select(c(1:21,23,1484:1493,24:1483))

names(total_input)[24:32] <- c("numeric_memory","fluid_intelligence","reaction_time","prospective_memory","visual_memory_round1","visual_memory_round2","visual_memory_round3","TMT_A","TMT_B")

total_input <- as.data.frame(total_input)
total_input <-total_input[,-c(31,32)]

load("/mnt/data/lijincheng/UKB/Olink/0627/stroke.Rdata")
id_stroke <- unique(long_stroke[which(long_stroke$stroke_group == 1),"f.eid"]) 

total_input <- total_input %>%  mutate(stroke = ifelse(eid %in% id_stroke,1,0 )) 
table(total_input$stroke)
##   0     1
## 47552   815 

total_input  <- total_input  %>% dplyr::select(c(1:30,1491,31:1490))

total_input <- total_input %>% mutate(edu_g = ifelse(edu == "-7","low",
                                                     ifelse(edu == "-3","missing","high")),
                                      ethnicity_g = ifelse(ethnicity == "White", "White",
                                                           ifelse( ethnicity == "no answer", "missing" ,"not White")))

total_input$edu_g <- factor(total_input$edu_g ,levels=c("high","low","missing"),labels=c("high","low","missing"))
total_input$ethnicity_g <- factor(total_input$ethnicity_g ,levels=c("White","not White","missing"),labels =c("White","not White","missing"))


total_input  <- total_input  %>% dplyr::select(c(1:8,1493,9,1492,10:1491))


total_input$allcause <- as.numeric(as.character(total_input$allcause))
total_input$allcause[which(total_input$allcause == 0)] <- 2  ###1:?в???2??û??
total_input$CVD <- factor(total_input$CVD,levels=c(0,1),labels = c(0,1))
total_input$sex <- relevel(total_input$sex,ref = "M")
total_input$smoking <- relevel(total_input$smoking,ref = "Never")
total_input$drinking <- relevel(total_input$drinking,ref = "Never")
total_input$ethnicity_g <- relevel(total_input$ethnicity_g,ref = "White")
total_input$edu_g <- relevel(total_input$edu_g ,ref = "low")
total_input$diabetes <- relevel(total_input$diabetes ,ref = "F")
total_input$hypertension <- relevel(total_input$hypertension ,ref = "F")
total_input$apoe_e4carrier <- relevel(total_input$apoe_e4carrier ,ref = "0")


total_input$AD <- as.numeric(as.character(total_input$AD))
total_input$AD[which(total_input$AD == 0)] <- 2  ###1:?в???2??û??

total_input$VD <- as.numeric(as.character(total_input$VD))
total_input$VD[which(total_input$VD == 0)] <- 2  ###1:?в???2??û??

total_input <- total_input %>% mutate(time_year = time_day / 365) %>%  dplyr::select(c(1:24,1494,25:1493))  


save(total_input,file="/mnt/data/lijincheng/UKB/Olink/0627/total_input.Rdata")

load(file="/mnt/data/lijincheng/UKB/Olink/0627/total_input.Rdata")
library("survival")
library("survminer")
source("/mnt/data/lijincheng/UKB/Olink/0627/get_cox.R")

get_cox(df = total_input,
        olinkrange=c(35:1494),
        outpath="/mnt/data/lijincheng/UKB/Olink/0627/cox_main.xlsx")




################------------------{02 sensitivity }-----------------------------------
#####-----------------------------{get totally independent proteins}-----------------------------------------------
#################------------{ exclude (1)dementia within two years; (2)those with history of stroke; (3)add fluid intelligence as a covariate}
#####1
load(file="/mnt/data/lijincheng/UKB/Olink/0627/total_input.Rdata")
library("survival")
library("survminer")
library(tidyverse)
source("/mnt/data/lijincheng/UKB/Olink/0627/get_cox.R")
se_input1 <- total_input %>% dplyr::filter(!(within_two_year ==1 & allcause ==1)) 
get_cox(df = se_input1,
        olinkrange=c(35:1494),
        outpath="/mnt/data/lijincheng/UKB/Olink/0627/cox_exclude_withintwo.xlsx")

#####2
load(file="/mnt/data/lijincheng/UKB/Olink/0627/total_input.Rdata")
library("survival")
library("survminer")
library(tidyverse)
source("/mnt/data/lijincheng/UKB/Olink/0627/get_cox.R")
se_input2_nonstroke <- total_input %>% filter(stroke == 0) 
get_cox(df = se_input2_nonstroke,
        olinkrange=c(35:1494),
        outpath="/mnt/data/lijincheng/UKB/Olink/0627/cox_exclude_stroke.xlsx")

#####3 in cox_main model3


###############----------------------------{03 subgroup }-------------------------------------
###############----------------------------{(1) age: <65 and >=65}---------------------------

load(file="/mnt/data/lijincheng/UKB/Olink/0627/total_input.Rdata")
library("survival")
library("survminer")
library(tidyverse)
source("/mnt/data/lijincheng/UKB/Olink/0627/get_cox.R")
age_less_65 <- total_input %>% filter(age < 65)
age_more_65 <- total_input %>% filter(age == 65 | age > 65)

get_cox(df = age_less_65,
        olinkrange=c(35:1494),
        outpath="/mnt/data/lijincheng/UKB/Olink/0627/cox_age_less_65.xlsx")

get_cox(df = age_more_65,
        olinkrange=c(35:1494),
        outpath="/mnt/data/lijincheng/UKB/Olink/0627/cox_age_more_65.xlsx")

###############----------------------------{(2) sex: male and female}---------------------------
load(file="/mnt/data/lijincheng/UKB/Olink/0627/total_input.Rdata")
library("survival")
library("survminer")
library(tidyverse)
source("/mnt/data/lijincheng/UKB/Olink/0627/get_cox.R")
male <- total_input %>% filter(sex == "M")
female <- total_input %>% filter(sex == "F")
get_cox(df = male,
        olinkrange=c(35:1494),
        outpath="/mnt/data/lijincheng/UKB/Olink/0627/cox_male.xlsx")

get_cox(df = female,
        olinkrange=c(34:1493),
        outpath="/mnt/data/lijincheng/UKB/Olink/0627/cox_female.xlsx")
###############----------------------------{(3) APOE carrier: 0 and >= 1}---------------------------
load(file="/mnt/data/lijincheng/UKB/Olink/0627/total_input.Rdata")
library("survival")
library("survminer")
library(tidyverse)
source("/mnt/data/lijincheng/UKB/Olink/0627/get_cox.R")
apoe_carrier <- total_input %>% filter(apoe_e4carrier == ">=1")
apoe_noncarrier <- total_input %>% filter(apoe_e4carrier == "0")
get_cox(df = apoe_carrier,
        olinkrange=c(35:1494),
        outpath="/mnt/data/lijincheng/UKB/Olink/0627/cox_apoe_carrier.xlsx")

get_cox(df = apoe_noncarrier,
        olinkrange=c(35:1494),
        outpath="/mnt/data/lijincheng/UKB/Olink/0627/cox_apoe_noncarrier.xlsx")

################---------------------------{04 association with cognitive function}------------------------------------
load(file="/mnt/data/lijincheng/UKB/Olink/0627/total_input.Rdata")
##numeric_memory????<0??ֵʵ??ΪNA
total_input$numeric_memory[which(total_input$numeric_memory < 0)] = NA

## prospective_memory ???¸?ֵscores
total_input$prospective_memory[which(total_input$prospective_memory == 0)]<-"Instruction not recalled,either skipped or incorrect"
total_input$prospective_memory[which(total_input$prospective_memory == 1)]<-"Correct recall on first attempt"
total_input$prospective_memory[which(total_input$prospective_memory == 2)]<-"Correct recall on second attempt"

total_input$prospective_memory[which(total_input$prospective_memory == "Instruction not recalled,either skipped or incorrect")]<- 0
total_input$prospective_memory[which(total_input$prospective_memory == "Correct recall on first attempt")]<- 2
total_input$prospective_memory[which(total_input$prospective_memory == "Correct recall on second attempt")]<-1


#####
##lm_reg <- total_input %>% gather(key = "olink",value = "epr", -c(names(total_input)[1:34]))  

library(stats)
library(broom)
library(purrr)
library(betareg)
glm_mcog <- total_input %>% 
  gather(key = "olink",value = "value", -c(names(total_input)[1:34]))  

###cognitive function regression
#####----------------------------------{ 01 Fluid intelligence }-------------------------------

fluid_beta_model2 <- glm_mcog %>%
  mutate(trans_fluid = (fluid_intelligence - 0) / (13 - 0)) %>%
  group_by(olink) %>%
  nest() %>%
  mutate(model_results = map(data, ~ {
    value_sd <- sd(.$value)
    .$value_scaled <- .$value / value_sd
    .$trans_fluid <- (.$trans_fluid * (sum(!is.na(.$fluid_intelligence)) - 1) + 0.5) / sum(!is.na(.$fluid_intelligence))
    tidy(betareg(trans_fluid ~ value_scaled +
                   age + TDI + sex + BMI + smoking + drinking + ethnicity_g +
                   edu_g + diabetes + hypertension + eGFRcr + apoe_e4carrier + CVD,
                 .))
  })) %>%
  unnest(model_results) %>%
  mutate(model = "fluid_beta_model2") %>%
  dplyr::filter(term == "value_scaled")

fluid_beta_model2$fdr = p.adjust(fluid_beta_model2$p.value, method = "fdr")
fluid_beta_model2$bonf = p.adjust(fluid_beta_model2$p.value, method = "bonferroni")
fluid_beta_model2 <- fluid_beta_model2 %>% mutate(cog= "fluid intellegence")
##########---------------------------{ 02 Numeric memory }-----------------------------------------------

numeric_beta_model2 <- glm_mcog %>%
  mutate(trans_numeric = (numeric_memory - 2) / (12 - 0)) %>%
  group_by(olink) %>%
  nest() %>%
  mutate(model_results = map(data, ~ {
    value_sd <- sd(.$value)
    .$value_scaled <- .$value / value_sd
    .$trans_numeric <- (.$trans_numeric * (sum(!is.na(.$numeric_memory)) - 1) + 0.5) / sum(!is.na(.$numeric_memory))
    tidy(betareg(trans_numeric ~ value_scaled +
                   age + TDI + sex + BMI + smoking + drinking + ethnicity_g +
                   edu_g + diabetes + hypertension + eGFRcr + apoe_e4carrier + CVD,
                 .))
  })) %>%
  unnest(model_results) %>%
  mutate(model = "numeric_beta_model2") %>%
  dplyr::filter(term == "value_scaled")

numeric_beta_model2$fdr = p.adjust(numeric_beta_model2$p.value, method = "fdr")
numeric_beta_model2$bonf = p.adjust(numeric_beta_model2$p.value, method = "bonferroni")
numeric_beta_model2 <- numeric_beta_model2 %>% mutate(cog= "numeric memory")
##########---------------------------{ 03 Prospective memory }-----------------------------------------------

prospective_logistic_model2<-  glm_mcog %>%
  mutate(prospective = ifelse(prospective_memory ==2,
                              1,
                              ifelse(prospective_memory < 2, 0,
                                     NA))) %>% 
  group_by(olink) %>%
  nest() %>%
  mutate(model_results = map(data, ~ {
    value_sd <- sd(.$value)
    .$value_scaled <- .$value / value_sd
    tidy(glm(prospective ~ value_scaled +
               age + TDI + sex + BMI + smoking + drinking + ethnicity_g +
               edu_g + diabetes + hypertension + eGFRcr + apoe_e4carrier + CVD,
             .,family = binomial(link = "logit")))
  })) %>%
  unnest(model_results) %>%
  mutate(model = "prospective_binomial_model2") %>%
  dplyr::filter(term == "value_scaled")

prospective_logistic_model2$fdr = p.adjust(prospective_logistic_model2$p.value, method = "fdr")
prospective_logistic_model2$bonf = p.adjust(prospective_logistic_model2$p.value, method = "bonferroni")

prospective_logistic_model2 <- prospective_logistic_model2 %>% 
  mutate(cog= "prospective memory")  
##########---------------------------{ 04 Visual memory/visual_memory_round2 }-----------------------------------------------

visual_possion_model<-  glm_mcog %>%
  group_by(olink) %>%
  nest() %>%
  mutate(model_results = map(data, ~ {
    value_sd <- sd(.$value)
    .$value_scaled <- .$value / value_sd
    tidy(glm(visual_memory_round2 ~ value_scaled +
               age + TDI + sex + BMI + smoking + drinking + ethnicity_g +
               edu_g + diabetes + hypertension + eGFRcr + apoe_e4carrier + CVD,
             .,family = poisson(link = "log")))
  })) %>%
  unnest(model_results) %>%
  mutate(model = "visual_possion_model2") %>%
  dplyr::filter(term == "value_scaled")
  
  visual_possion_model2$fdr = p.adjust(visual_possion_model2$p.value, method = "fdr")
visual_possion_model2$bonf = p.adjust(visual_possion_model2$p.value, method = "bonferroni")

visual_possion_model2 <- visual_possion_model2 %>% 
  mutate(cog= "visual memory")

###########---------------------------{ 05 Reaction time/reaction_time }-----------------------------------------------
reaction_poisson_model2 <- glm_mcog%>% 
  group_by(olink) %>%
  nest() %>%
  mutate(model_results = map(data, ~ {
    value_sd <- sd(.$value)
    .$value_scaled <- .$value / value_sd
    tidy(glm(reaction_time ~ value_scaled +
               age + TDI + sex + BMI + smoking + drinking + ethnicity_g +
               edu_g + diabetes + hypertension + eGFRcr + apoe_e4carrier + CVD,
             .,family = poisson(link = "log")))
  })) %>%
  unnest(model_results) %>%
  mutate(model = "reaction_poisson_model2") %>%
  dplyr::filter(term == "value_scaled") 

reaction_possion_model2$fdr = p.adjust(reaction_possion_model2$p.value, method = "fdr")
reaction_possion_model2$bonf = p.adjust(reaction_possion_model2$p.value, method = "bonferroni")
reaction_possion_model2 <- reaction_possion_model2 %>%
  mutate(cog= "reaction time")

olink_reg_cog <- bind_rows(fluid_beta_model2[,c("Network","estimate","std.error","statistic","p.value","bonf","cog")],
                            numeric_beta_model2[,c("Network","estimate","std.error","statistic","p.value","bonf","cog")],
                            prospective_logistic_model2[,c("Network","estimate","std.error","statistic","p.value","bonf","cog")],
                            visual_possion_model2[,c("Network","estimate","std.error","statistic","p.value","bonf","cog")],
                            reaction_possion_model2[,c("Network","estimate","std.error","statistic","p.value","bonf","cog")])

write.csv(olink_reg_cog,"/mnt/data/lijincheng/UKB/Olink/0627/olink_reg_cog.csv",row.names=F)

#########-----------------------------------------{prediction}---------------------------------------------------------------------
###run by Python

########-----------------------------------------{MRI regression}--------------------------------------------------------------------
setwd("/mnt/data/lijincheng/UKB/Olink/0627/MRI/")
load(file="/mnt/data/lijincheng/UKB/Olink/0627/total_input.Rdata")
library(tidyverse)
library(data.table)
library(vroom)
library(openxlsx)
MRI <- vroom("/mnt/data/lijincheng/UKB/Olink/0627/ukb673193.tab")

MRI1<- read.xlsx("/mnt/data/lijincheng/UKB/Olink/0627/MRI/UKB_MRIidlist.xlsx",sheet=1) %>% mutate(eid = paste0("f.",feid,".2.0")) ##56
MRI2<- read.xlsx("/mnt/data/lijincheng/UKB/Olink/0627/MRI/UKB_MRIidlist.xlsx",sheet=2) %>% mutate(eid = paste0("f.",feid,".2.0")) ##121
MRI3<- read.xlsx("/mnt/data/lijincheng/UKB/Olink/0627/MRI/UKB_MRIidlist.xlsx",sheet=3) %>% mutate(eid = paste0("f.",feid,".2.0")) ##14
MRI4<- read.xlsx("/mnt/data/lijincheng/UKB/Olink/0627/MRI/UKB_MRIidlist.xlsx",sheet=4) %>% mutate(eid = paste0("f.",feid,".2.0")) ##66
MRI5<- read.xlsx("/mnt/data/lijincheng/UKB/Olink/0627/MRI/UKB_MRIidlist.xlsx",sheet=5) %>% mutate(eid = paste0("f.",feid,".2.0")) ##48
MRI6<- read.xlsx("/mnt/data/lijincheng/UKB/Olink/0627/MRI/UKB_MRIidlist.xlsx",sheet=6) %>% mutate(eid = paste0("f.",feid,".2.0")) ##48

MRI_input <- MRI %>% dplyr::select(c("f.eid",MRI1$eid,MRI2$eid,MRI3$eid,MRI4$eid,MRI5$eid,MRI6$eid,"f.25009.2.0")) ###f.25009.2.0 	Volume of brain, grey+white matter (normalised for head size)
MRI_input <- as.data.frame(MRI_input)
MRI_input$f.eid <- as.character(MRI_input$f.eid)
##olink=35-1494 ; MRI=1495-1848
total_MRI <- dplyr::left_join(total_input,MRI_input,by=c("eid"="f.eid")) 
load(file="/mnt/data/lijincheng/UKB/Olink/0627/MRI/gaptime.Rdata")
names(UKB_v0v2_gap)[4:5] <- c("gap_day","gap_year")
UKB_v0v2_gap$f.eid <- as.character(UKB_v0v2_gap$f.eid)
total_MRI <- dplyr::left_join(total_MRI, UKB_v0v2_gap[,c(1,4,5)],by=c("eid"="f.eid"))
save(total_MRI,file="/mnt/data/lijincheng/UKB/Olink/0627/MRI/total_MRI.Rdata")

library(stats)
library(broom)
library(purrr)
response_vars <- c(names(total_MRI)[1495:1847]) # ??Ӧ??��??????
predictors <- data.frame(rawname =  c(names(total_MRI)[35:1494]) ,
                         runname = gsub("-", "_",  c(names(total_MRI)[35:1494]) ))
predictor_vars <- predictors$runname  # ?Ա?��??????
names(total_MRI)[35:1494] <- predictors$runname

regression_results <- data.frame()

for (response_var in response_vars) {
  for (predictor_var in predictor_vars) {
    total_MRI[[response_var]] <- total_MRI[[response_var]]
    total_MRI[[predictor_var]] <- total_MRI[[predictor_var]] / sd(total_MRI[[predictor_var]])
    
    # ???嵱ǰ?ع?ģ?͵? formula
    lm_formula <- as.formula(paste0(response_var, "~", predictor_var, "+", 
                                    "age+TDI+sex+BMI+smoking+drinking+ethnicity_g+edu_g+diabetes+hypertension+eGFRcr+apoe_e4carrier+CVD+gap_year"))
    
    # ִ?лع?????
    lm_result <- lm(lm_formula, data = total_MRI)
    tidy_result <- tidy(lm_result)
    tidy_result$ResponseVar <- response_var
    tidy_result$PredictorVar <- predictor_var
    
    # ????ǰ??Ӧ??��?????лع????????���һ?????ݿ򣬲????ӵ????????б???
    regression_results <- bind_rows( regression_results,tidy_result[2,])
  }
}
write.csv(regression_results,"/mnt/data/lijincheng/UKB/Olink/0627/final_olink_MRI_regression1011.csv",row.names=F)



#### module MRI regression
library(openxlsx)
library(tidyverse)
library(stats)
library(broom)
library(purrr)
setwd("/mnt/data/lijincheng/UKB/Olink/0627/MRI/")
load(file="/mnt/data/lijincheng/UKB/Olink/0627/MRI/total_MRI.Rdata")
allprotein_module <- read.xlsx("/mnt/data/lijincheng/UKB/Olink/0627/final_olink_WGCNA_module.xlsx",sheet=2)

module_MRI_reg <- function(predictors,response_vars,inputdata){
  regres <- data.frame()
  for (response_var in response_vars) {
    # ????һ???յ??б���???ڴ洢??ǰ??Ӧ??��?????лع?????
    # ѭ???????Ա?��
    for (predictor_var in predictors) {
      # ???嵱ǰ?ع?ģ?͵? formula
      lm_formula <- as.formula(paste0(response_var, "~", predictor_var, "+", 
                                      "age+TDI+sex+BMI+smoking+drinking+ethnicity_g+edu_g+diabetes+hypertension+eGFRcr+apoe_e4carrier+CVD+gap_year"))
      
      inputdata[[response_var]] <- inputdata[[response_var]] ## inputdata$f.25009.2.0  ??????????????  inputdata[[response_var]] <- inputdata[[response_var]] / inputdata$f.25009.2.0  ## ????????????
      # ?????Ա?��?ı?׼??
      predictor_sd <- sd(inputdata[[predictor_var]])
      # ???Ա?��??ֵ?????????ı?׼??
      inputdata[[predictor_var]] <- inputdata[[predictor_var]] / predictor_sd
      # ִ?лع?????
      lm_result <- lm(lm_formula, data = inputdata)
      tidy_result <- tidy(lm_result)
      tidy_result$ResponseVar <- response_var
      tidy_result$PredictorVar <- predictor_var
      # ????ǰ??Ӧ??��?????лع????????���һ?????ݿ򣬲????ӵ????????б???
      regres <- bind_rows(regres, tidy_result[2,])
    }
  }
  return(regres)
}

allprotein_MRI_module <- dplyr::left_join(allprotein_module,total_MRI,by= c("sample"="eid"))
allprotein_predictors <- names(allprotein_module)[2:13]
response_vars <- c(names(total_MRI)[1495:1847])
allprotein_reg <- module_MRI_reg(allprotein_predictors,response_vars,allprotein_MRI_module)

write.csv(allprotein_reg,"/mnt/data/lijincheng/UKB/Olink/0627/final_MRI_module_linearreg.csv",row.names=F)




###supplement module cox reg in model2
library("survival")
library("survminer")
library(tidyverse)
library(openxlsx)
get_module_cox <- function(input,module,disease){
  ######{model2: model1 + ethnicity, TDI, BMI, smoking, eGFR, hypertension (yes/no), diabetes mellitus (yes/no), and CVD (yes/no)}
  out <- data.frame()
  if(disease == "allcause"){
    for (i in module) {
      expr = input[,i]
      cox = coxph(Surv(time_year, allcause) ~ expr+
                    age+TDI+sex+BMI+smoking+drinking+ethnicity_g+
                    edu_g+diabetes+hypertension+eGFRcr+apoe_e4carrier+CVD
                  ,input)
      coxsummary = summary(cox)
      out=rbind(out,cbind(protein=i,
                          HR=round(coxsummary$coefficients[1,"exp(coef)"],2), 
                          z=round(coxsummary$coefficients[1,"z"],2), 
                          "95% CI"=paste(round(coxsummary$conf.int[1,3], 2), round(coxsummary$conf.int[1,4],2), sep = "-"),
                          pvalue=coxsummary$coefficients[1,"Pr(>|z|)"])) 
    }
    
  }else if (disease == "AD"){
    for (i in module) {
      expr = input[,i]
      cox = coxph(Surv(time_year, AD) ~ expr+
                    age+TDI+sex+BMI+smoking+drinking+ethnicity_g+
                    edu_g+diabetes+hypertension+eGFRcr+apoe_e4carrier+CVD
                  ,input)
      coxsummary = summary(cox)
      out=rbind(out,cbind(protein=i,
                          HR=round(coxsummary$coefficients[1,"exp(coef)"],2), 
                          z=round(coxsummary$coefficients[1,"z"],2), 
                          "95% CI"=paste(round(coxsummary$conf.int[1,3], 2), round(coxsummary$conf.int[1,4],2), sep = "-"),
                          pvalue=coxsummary$coefficients[1,"Pr(>|z|)"])) 
    }
    
  }else if (disease == "VD"){
    for (i in module) {
      expr = input[,i]
      cox = coxph(Surv(time_year, VD) ~ expr+
                    age+TDI+sex+BMI+smoking+drinking+ethnicity_g+
                    edu_g+diabetes+hypertension+eGFRcr+apoe_e4carrier+CVD
                  ,input)
      coxsummary = summary(cox)
      out=rbind(out,cbind(protein=i,
                          HR=round(coxsummary$coefficients[1,"exp(coef)"],2), 
                          z=round(coxsummary$coefficients[1,"z"],2), 
                          "95% CI"=paste(round(coxsummary$conf.int[1,3], 2), round(coxsummary$conf.int[1,4],2), sep = "-"),
                          pvalue=coxsummary$coefficients[1,"Pr(>|z|)"])) 
    }
  }
  
  
  out[,"fdr"] <- p.adjust(out$pvalue,method = "fdr")
  out[,"bonf"] <- p.adjust(out$pvalue,method = "bonferroni")
  
  return(out)
}

allprotein_module <- read.xlsx("/mnt/data/lijincheng/UKB/Olink/0627/final_olink_WGCNA_module.xlsx",sheet=2)

allprotein_module[,c(2:13)] <- allprotein_module[,c(2:13)] / sapply(allprotein_module[,c(2:13)], sd)


load(file="/mnt/data/lijincheng/UKB/Olink/0627/total_input.Rdata")

allprotein_module_input <- dplyr::left_join(allprotein_module,total_input,by=c("sample"="eid")) 
###ԭʼ??status?? 1=??????2=û?л?????Ϊ??????cox?ع???Ҫ????1=û?л?????2=????
allprotein_module_input  <-  allprotein_module_input %>%
  mutate(allcause = case_when(
    allcause == 1 ~ 2,
    allcause == 2 ~ 1,
    TRUE ~ allcause  # ???????????????????????ֲ???
  ),
  AD = case_when(
    AD == 1 ~ 2,
    AD == 2 ~ 1,
    TRUE ~ AD  # ???????????????????????ֲ???
  ),
  VD = case_when(
    VD == 1 ~ 2,
    VD == 2 ~ 1,
    TRUE ~ VD  # ???????????????????????ֲ???
  ))


##??ʱ1=û?л?????2=????
allcause_module_input <- allprotein_module_input  
AD_module_input <- allprotein_module_input %>% dplyr::filter(allcause == 1 | (VD == 1 & AD == 2)) 
VD_module_input <- allprotein_module_input %>% dplyr::filter(allcause == 1 | (AD == 1 & VD == 2)) 


modulename <- names(allprotein_module)[2:13]

allcause_module_cox <-  get_module_cox(allcause_module_input,modulename,"allcause")
AD_module_cox <- get_module_cox(AD_module_input,modulename,"AD")
VD_module_cox <- get_module_cox(VD_module_input,modulename,"VD")

writexl::write_xlsx(x = list("allcause_module_cox"=allcause_module_cox,"AD_module_cox"=AD_module_cox,"VD_module_cox"=VD_module_cox),path="/mnt/data/lijincheng/UKB/Olink/0627/final_module_coxall.xlsx")





