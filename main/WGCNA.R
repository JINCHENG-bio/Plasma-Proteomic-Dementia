#####WGCNA
rm(list=ls())
load(file = "E://UKB//Olink//0627//total_input.Rdata")
library(tidyverse)
setwd("E://UKB//Olink//0627//final_WGCNA1004/")

library(openxlsx)
cox <- read.xlsx("E://UKB//Olink//0627//cox_res1004//cox_main.xlsx",sheet = 2)

cox[,c(2,3,5,6,7,8,9,11,12,13,14,15,17,18,19)] <- apply(cox[,c(2,3,5,6,7,8,9,11,12,13,14,15,17,18,19)],2,as.numeric)

FS <- cox[order(cox$pvalue_AD,decreasing=F), ]

#top_50 <- head(FS, 50) %>% pull("protein")


allcause_protein <- cox %>% filter(bonf_allcause < 0.05) %>% pull(protein)
AD_protein <- cox %>% filter(bonf_AD < 0.05) %>% pull(protein)
VD_protein <- cox %>% filter(bonf_VD < 0.05) %>% pull(protein)

selected_protein <- unique(c(allcause_protein,AD_protein,VD_protein))
length(unique(c(allcause_protein,AD_protein,VD_protein)))
allprotein <- cox$protein
####
####WGCNA
library(WGCNA)
options(stringsAsFactors = FALSE);
library(reshape2)
library(stringr)
### 1 allcause   
row.names(total_input) <- total_input[,1]
df <- total_input[,allprotein]

###检查离群值
gsg = goodSamplesGenes(df, verbose = 3);
gsg$allOK
###TRUE
if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}


sampleTree = hclust(dist(df), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
pdf(file = "sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)

dev.off()



powers = c(c(1:10), seq(from = 12, to=30, by=2))
sft = pickSoftThreshold(df, powerVector=powers, 
                        networkType="unsigned", verbose=5,
)
cex1 = 0.9;

plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
abline(h=0.90,col="red")  #查看位于0.9以上的点，可以改变高度值

plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5],cex=cex1, labels=powers,col="red")


sft$powerEstimate  #自动推荐的结果为6

net = blockwiseModules(df, power = 6, deepSplit = 4,##最敏感
                       TOMType = "signed", minModuleSize = 10,
                       reassignThreshold = 0, mergeCutHeight = 0.05,
                       useBranchEigennodeDissim = T,
                       corType = "bicor",
                       randomSeed = 123,
                       TOMDenom = "mean",
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "allprotein",
                       verbose = 3)
# power = 6是刚才选择的软阈值
#minModuleSize：模块中最少的基因数
#mergeCutHeight ：模块合并阈值，阈值越大，模块越少（重要）
#saveTOMs = TRUE,saveTOMFileBase = "allprotein"保存TOM矩阵，名字为"allprotein"
#net$colors 包含模块分配，net$MEs 包含模块的模块特征基因。

table(net$colors)


sizeGrWindow(12, 9)
mergedColors = labels2colors(net$colors)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

pdf(file = "module.pdf", width = 12, height = 9);
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

#colnames(df)[moduleColors=="turquoise"]#返回属于棕色模块的基因ID
#annot = read.csv(file = "GeneAnnotation.csv");
#dim(sigDEG)
#names(sigDEG)
#probes = colnames(expinput) # 匹配信息
#probes2annot = match(probes, sigDEG$mirna);
#sum(is.na(probes2annot)) # 检测是否有没有匹配上的ID号，正常来说为0，即全匹配上了。
#输出必要的信息：
# 重新计算带有颜色标签的模块
MEs0 = moduleEigengenes(df, mergedColors)$eigengenes
MEs = orderMEs(MEs0)
geneModuleMembership = as.data.frame(cor(df, MEs, use = "p",method = "spearman"));
moduleColors <- mergedColors
olinkInfo0 = data.frame(olink = colnames(df),
                       moduleColor = moduleColors,
                       geneModuleMembership)

MEs$sample <- row.names(MEs)

MEs <- MEs[,c(13,1:12)]

openxlsx::write.xlsx(x=list("olinkInfo0" = olinkInfo0,
                         "Module"=MEs),file = "final_olink_WGCNA_module.xlsx")
save.image(file="WGCNA.Rdata")
#load(file="WGCNA/WGCNA.Rdata")


#####----------------{模块内关联网络，根据membership得到hub蛋白}---------------------------------
rm(list = ls())
setwd("E://UKB//Olink//0627//final_WGCNA1004/")
library(tidyverse)
library(openxlsx)
library(data.table)
load(file = "E://UKB//Olink//0627//total_input.Rdata")
module <- openxlsx::read.xlsx("E://UKB//Olink//0627//final_WGCNA1004/final_olink_WGCNA_module.xlsx",sheet = 1)
network <- openxlsx::read.xlsx("E://UKB//Olink//0627//final_WGCNA1004/final_olink_WGCNA_module.xlsx",sheet = 2)
olinkref <- fread("E://UKB//Olink//0627//final_WGCNA1004/protein_code.tsv")
panel <- fread("E://UKB//Olink//0627//final_WGCNA1004/olink_assay.dat.txt")

protein <- fread("E://UKB//Olink//0627//final_WGCNA1004/protein_code.tsv") %>% 
  separate(.,col="meaning",into = c("name","full name"),sep = ";")

panel <- left_join(panel,protein[,c(2,3)],by=c("Assay"="name"))

cox_main <- read.xlsx("E://UKB//Olink//0627//cox_res1004//cox_main.xlsx",sheet=1)

full <- left_join(cox_main[,c(1,2)],panel,by=c("protein"="Assay")) %>% 
  dplyr::select(c(1,3,4,5))
write.csv(full,"E://UKB//Olink//0627//final_WGCNA1004/full.csv",row.names = F)

############----------------------{将模块与认知关联}----------------------
##numeric_memory含有<0的值实际为NA
total_input$numeric_memory[which(total_input$numeric_memory < 0)] = NA

## prospective_memory 重新赋值scores
total_input$prospective_memory[which(total_input$prospective_memory == 0)]<-"Instruction not recalled,either skipped or incorrect"
total_input$prospective_memory[which(total_input$prospective_memory == 1)]<-"Correct recall on first attempt"
total_input$prospective_memory[which(total_input$prospective_memory == 2)]<-"Correct recall on second attempt"

total_input$prospective_memory[which(total_input$prospective_memory == "Instruction not recalled,either skipped or incorrect")]<- 0
total_input$prospective_memory[which(total_input$prospective_memory == "Correct recall on first attempt")]<- 2
total_input$prospective_memory[which(total_input$prospective_memory == "Correct recall on second attempt")]<-1

mcog <- left_join(total_input[,c(1:34)],network,by=c("eid"="sample"))

###矫正协变量，线性回归
library(stats)
library(broom)
library(purrr)
library(betareg)
glm_mcog <- mcog%>% 
  gather(key = "Network",value = "value", -c(names(mcog)[1:34]))  

glm_mcog$prospective_memory <- as.numeric(glm_mcog$prospective_memory)
###beta回归 - 有边界的认知测试：fluid intelligence；Numeric memory；Prospective Memory 
fluid_beta_model2 <- glm_mcog %>%
  mutate(trans_fluid = (fluid_intelligence - 0) / (13 - 0)) %>%
  group_by(Network) %>%
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

numeric_beta_model2 <- glm_mcog %>%
  mutate(trans_numeric = (numeric_memory - 2) / (12 - 0)) %>%
  group_by(Network) %>%
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


###广义线性回归 - 全正数的认知测试：
##prospective memory(二项式分布)
##Visual memory(泊松分布) ; 
##Reaction time(泊松分布)
prospective_logistic_model2 <-  glm_mcog %>%
  mutate(prospective = ifelse(prospective_memory ==2,
                              1,
                              ifelse(prospective_memory < 2, 0,
                                     NA))) %>% 
  group_by(Network) %>%
  nest() %>%
  mutate(model_results = map(data, ~ {
    # 计算自变量 value 的标准差
    value_sd <- sd(.$value)
    # 将 value 除以标准差，得到新的变量 value_scaled
    .$value_scaled <- .$value / value_sd
    # 进行线性回归，并返回结果
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


visual_possion_model2 <- glm_mcog %>%
  group_by(Network) %>%
  nest() %>%
  mutate(model_results = map(data, ~ {
    # 计算自变量 value 的标准差
    value_sd <- sd(.$value)
    # 将 value 除以标准差，得到新的变量 value_scaled
    .$value_scaled <- .$value / value_sd
    # 进行线性回归，并返回结果
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

reaction_possion_model2 <- glm_mcog %>%
  group_by(Network) %>%
  nest() %>%
  mutate(model_results = map(data, ~ {
    # 计算自变量 value 的标准差
    value_sd <- sd(.$value)
    # 将 value 除以标准差，得到新的变量 value_scaled
    .$value_scaled <- .$value / value_sd
    # 进行线性回归，并返回结果
    tidy(glm(reaction_time ~ value_scaled +
              age + TDI + sex + BMI + smoking + drinking + ethnicity_g +
              edu_g + diabetes + hypertension + eGFRcr + apoe_e4carrier + CVD,
             .,family = poisson(link = "log")))
  })) %>%
  unnest(model_results) %>%
  mutate(model = "reaction_possion_model2") %>%
  dplyr::filter(term == "value_scaled")

reaction_possion_model2$fdr = p.adjust(reaction_possion_model2$p.value, method = "fdr")
reaction_possion_model2$bonf = p.adjust(reaction_possion_model2$p.value, method = "bonferroni")
reaction_possion_model2 <- reaction_possion_model2 %>%
  mutate(cog= "reaction time")




#############

module_reg_cog <- bind_rows(fluid_beta_model2[,c("Network","estimate","std.error","statistic","p.value","bonf","cog")],
                           numeric_beta_model2[,c("Network","estimate","std.error","statistic","p.value","bonf","cog")],
                           prospective_logistic_model2[,c("Network","estimate","std.error","statistic","p.value","bonf","cog")],
                           visual_possion_model2[,c("Network","estimate","std.error","statistic","p.value","bonf","cog")],
                           reaction_possion_model2[,c("Network","estimate","std.error","statistic","p.value","bonf","cog")])

module_reg_cog <- module_reg_cog %>% 
  mutate(signif =ifelse(bonf < 0.001 , "***",
                        ifelse(bonf < 0.01, "**",
                               ifelse(bonf < 0.05, "*",""))))

module_ref <- as.data.frame(table(module$moduleColor))
module_ref$ME <- paste("ME", module_ref$Var1, sep = "")
module_ref <- module_ref[order(module_ref$Freq,decreasing = T),]
module_ref <- module_ref %>% 
  mutate(network = paste0("network", 1:12))

module_reg_cog <- left_join(module_reg_cog,module_ref[,c(3,4)],by=c("Network"="ME")) %>% 
                  mutate(module = case_when(
                    network == "network1" ~ "M0",
                    network == "network2" ~ "M1",
                    network == "network3" ~ "M2",
                    network == "network4" ~ "M3",
                    network == "network5" ~ "M4",
                    network == "network6" ~ "M5",
                    network == "network7" ~ "M6",
                    network == "network8" ~ "M7",
                    network == "network9" ~ "M8",
                    network == "network10" ~ "M9",
                    network == "network11" ~ "M10",
                    network == "network12" ~ "M11",
                    TRUE ~ NA_character_
                  ))

#module_reg_cog$network <- factor(module_reg_cog$network,levels = rev(c(paste0("network", 1:11))),
#                                labels = rev(c(paste0("network", 1:11))))
module_reg_cog$cog <- factor(module_reg_cog$cog,
                            levels = c("fluid intellegence","numeric memory",
                                       "prospective memory",
                                       "visual memory","reaction time"),
                            labels = c("fluid intellegence","numeric memory",
                                       "prospective memory",
                                       "visual memory","reaction time"))

module_reg_cog <- module_reg_cog %>% 
  mutate(type = ifelse( cog == "fluid intellegence" | cog == "numeric memory",
                        "Beta regression",
                        ifelse( cog == "prospective memory" , "binomial",
                                "poisson")))

library(ggplot2)
library(ggsci)
# Plot 
p_module_cog <- ggplot(aes(x=cog, y=module , fill=estimate), data=module_reg_cog)
fig_module_cog <- p_module_cog + geom_tile() + scale_fill_gradient2(low="#2b8cbe", mid="white", high="#e34a33") + 
  #   geom_text(aes(label=stars, color=value), size=8) + scale_colour_gradient(low="grey30", high="white", guide="none") +
  geom_text(aes(label=signif), color="black", size=10) + 
  labs(y=NULL, x=NULL, fill="") + 
  scale_y_discrete(limits = c("M11","M10","M9","M8","M7","M6","M5","M4","M3","M2","M1","M0"))+
  ggtitle("")+
  theme_minimal()  + 
  theme(axis.text.x=element_text(angle = 0, hjust = 0,size = 12),
        axis.text.y=element_text(angle = 0, hjust = 0,size = 15))+
  facet_grid(~ type, scales='free_x', space="free_x")

fig_module_cog
###
write.csv(module_reg_cog,"allprotein_modeule_cogreg.csv",row.names = F)
save.image(file = "WGCNA_cog_reg.Rdata")

###########-------------------{蛋白模块的cox回归}-------------------------
library(survival)
library(dplyr)
library(tidyr)
library(purrr)
library(broom)
allcause_module_cox <- glm_mcog %>%
  group_by(Network) %>%
  nest() %>%
  mutate(model_results = map(data, ~ {
    # 创建生存对象
    surv_obj <- Surv(.$time_year, .$allcause)
    # 计算自变量 value 的标准差
    value_sd <- sd(.$value)
    # 将 value 除以标准差，得到新的变量 value_scaled
    .$value_scaled <- .$value / value_sd
    # 进行Cox回归，并返回结果
    tidy(coxph(surv_obj ~ value_scaled +
                 age + TDI + sex + BMI + smoking + drinking + ethnicity_g +
                 edu_g + diabetes + hypertension + eGFRcr + apoe_e4carrier + CVD,
               data = .),exponentiate = TRUE,conf.int=T)
  })) %>%
  unnest(model_results) %>%
  mutate(model = "allcause_cox_model2") %>%
  filter(term == "value_scaled") 

allcause_module_cox$bonf <- p.adjust(allcause_module_cox$p.value, method = "bonferroni")
  
AD_module_cox <- glm_mcog %>%
  dplyr::filter(allcause == 2 | (VD == 2 & AD == 1)) %>% 
  group_by(Network) %>%
  nest() %>%
  mutate(model_results = map(data, ~ {
    # 创建生存对象
    surv_obj <- Surv(.$time_year, .$AD)
    # 计算自变量 value 的标准差
    value_sd <- sd(.$value)
    # 将 value 除以标准差，得到新的变量 value_scaled
    .$value_scaled <- .$value / value_sd
    # 进行Cox回归，并返回结果
    tidy(coxph(surv_obj ~ value_scaled +
                 age + TDI + sex + BMI + smoking + drinking + ethnicity_g +
                 edu_g + diabetes + hypertension + eGFRcr + apoe_e4carrier + CVD,
               data = .),exponentiate = TRUE,conf.int=T)
  })) %>%
  unnest(model_results) %>%
  mutate(model = "AD_cox_model2") %>%
  filter(term == "value_scaled")
AD_module_cox$bonf <- p.adjust(AD_module_cox$p.value, method = "bonferroni")


VD_module_cox <- glm_mcog %>%
  dplyr::filter(allcause == 2 | (AD == 2 & VD == 1)) %>% 
  group_by(Network) %>%
  nest() %>%
  mutate(model_results = map(data, ~ {
    # 创建生存对象
    surv_obj <- Surv(.$time_year, .$VD)
    # 计算自变量 value 的标准差
    value_sd <- sd(.$value)
    # 将 value 除以标准差，得到新的变量 value_scaled
    .$value_scaled <- .$value / value_sd
    # 进行Cox回归，并返回结果
    tidy(coxph(surv_obj ~ value_scaled +
                 age + TDI + sex + BMI + smoking + drinking + ethnicity_g +
                 edu_g + diabetes + hypertension + eGFRcr + apoe_e4carrier + CVD,
               data = .),exponentiate = TRUE,conf.int=T)
  })) %>%
  unnest(model_results) %>%
  mutate(model = "VD_cox_model2") %>%
  filter(term == "value_scaled") 
VD_module_cox$bonf <- p.adjust(VD_module_cox$p.value, method = "bonferroni")


allcause_module_cox <- allcause_module_cox[,-2] %>% left_join(.,module_ref[,c(3,4)],by=c("Network"="ME"))
allcause_module_cox <- allcause_module_cox %>% dplyr::select(c(1,12,2:11))

AD_module_cox <- AD_module_cox[,-2] %>% left_join(.,module_ref[,c(3,4)],by=c("Network"="ME"))
AD_module_cox <- AD_module_cox %>% dplyr::select(c(1,12,2:11))

VD_module_cox <- VD_module_cox[,-2] %>% left_join(.,module_ref[,c(3,4)],by=c("Network"="ME"))
VD_module_cox <- VD_module_cox %>% dplyr::select(c(1,12,2:11))

write.xlsx(x= list("allcause_module_cox"=allcause_module_cox,
                   "AD_module_cox"=AD_module_cox,
                   "VD_module_cox"=VD_module_cox),
           file = "E://UKB//Olink//0627//WGCNA//module_cox_model2.xlsx")


############-----------------{蛋白质网络模块hub蛋白}----------------------
olink_group <- left_join(module[,1:2],module_ref,by=c("moduleColor"="Var1")) %>% 
  left_join(.,full[,],by=c("olink"="protein"))

###spearman rank
###----------------------{1.网络模块连通性}-----------------------
library(rms)
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}

connectivity <- data.frame()
cor_net <- data.frame()
for(i in c(unique(olink_group$network))){
  cor_data <- data.frame()
  id <- olink_group %>% dplyr::filter(network == i) %>% pull(olink)
  data <- total_input[,id]
  cor_result <- rcorr(as.matrix(data),type = "spearman")
  tmp_cor <- flattenCorrMatrix(cor_result$r,cor_result$P)
  tmp_cor$Net <- i
  cor_net <- rbind(cor_net,tmp_cor)
  
  cor_matrix <- cor_result$r
  # 取相关性系数矩阵的绝对值
  abs_cor_matrix <- abs(cor_matrix)
  # 计算绝对值相关性系数矩阵每一行的和，即每个变量与其他变量的相关性系数绝对值的和
  sum_cor <- apply(abs_cor_matrix, 1, sum)
  # 计算和
  cor_data <- data.frame(Variable = id, SumCor = sum_cor,
                         Net = i)
  
  connectivity <- rbind(connectivity,cor_data)
}

connectivity <- connectivity %>%
  arrange(Net, desc(SumCor))

cor_net ###用以绘制相关性网络

####--------------------------{2.Module membership}--------------------------------
olink_module <- left_join(total_input[,c("eid",module$olink)],network,by=c("eid"="sample"))
row.names(olink_module) <- olink_module[,1]
olink_module <- olink_module[,-1]

cor_result <- rcorr(as.matrix(olink_module),type = "spearman")

MM <- flattenCorrMatrix(cor_result$r,cor_result$P)

MM <- MM %>% dplyr::filter( (row %in% olink_group[which(olink_group$network == "network1"),"olink"] & column %in% olink_group[which(olink_group$network == "network1"),"ME"]) |
                              (row %in% olink_group[which(olink_group$network == "network2"),"olink"] & column %in% olink_group[which(olink_group$network == "network2"),"ME"]) | 
                              (row %in% olink_group[which(olink_group$network == "network3"),"olink"] & column %in% olink_group[which(olink_group$network == "network3"),"ME"]) |
                              (row %in% olink_group[which(olink_group$network == "network4"),"olink"] & column %in% olink_group[which(olink_group$network == "network4"),"ME"]) |
                              (row %in% olink_group[which(olink_group$network == "network5"),"olink"] & column %in% olink_group[which(olink_group$network == "network5"),"ME"]) |
                              (row %in% olink_group[which(olink_group$network == "network6"),"olink"] & column %in% olink_group[which(olink_group$network == "network6"),"ME"]) |
                              (row %in% olink_group[which(olink_group$network == "network7"),"olink"] & column %in% olink_group[which(olink_group$network == "network7"),"ME"]) |
                              (row %in% olink_group[which(olink_group$network == "network8"),"olink"] & column %in% olink_group[which(olink_group$network == "network8"),"ME"]) |
                              (row %in% olink_group[which(olink_group$network == "network9"),"olink"] & column %in% olink_group[which(olink_group$network == "network9"),"ME"]) |
                              (row %in% olink_group[which(olink_group$network == "network10"),"olink"] & column %in% olink_group[which(olink_group$network == "network10"),"ME"]) |
                              (row %in% olink_group[which(olink_group$network == "network11"),"olink"] & column %in% olink_group[which(olink_group$network == "network11"),"ME"])|
                              (row %in% olink_group[which(olink_group$network == "network12"),"olink"] & column %in% olink_group[which(olink_group$network == "network12"),"ME"]))

MM <- MM %>% left_join(.,olink_group[,c(1,5)],by=c("row"="olink")) 

MM <- MM %>%
  mutate( abscor= abs(cor)) %>% 
  arrange(network, desc(abscor))

write.xlsx(x = list("olink module" = olink_group, 
                    "Connectivity" = connectivity,
                    "Membership"= MM,
                    "Correlation net" = cor_net),
           file = "WGCNA_network.xlsx")
write.csv(olink_group,"olink_group.csv",row.names = F)

top5_connect <- connectivity %>%
  group_by(Net) %>%
  top_n(5, SumCor)

top5_MM <- MM %>%
  group_by(network) %>%
  top_n(5, abscor)


###connectivity
top1_connect <- connectivity %>%
  group_by(Net) %>%
  top_n(1, SumCor)

top1_MM <- MM %>%
  group_by(network) %>%
  top_n(1, abscor)

hub_olink <- left_join(top1_connect,olink_group[,c(1,6,7,8)],
                       by = c("Variable"="olink"))
write.csv(hub_olink,"WGCNA/hub_olink.csv",row.names = F)

save.image(file="WGCNA/module_hub_cox.Rdata")


hub5_olink <- left_join(top5_connect,olink_group[,c(1,6,7,8)],
                       by = c("Variable"="olink"))

write.csv(hub5_olink,"WGCNA/hub5_olink.csv",row.names = F)

save.image(file="WGCNA/module_hub_cox.Rdata")





