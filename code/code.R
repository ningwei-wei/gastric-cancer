library(sigminer) ### load package
library(dplyr)
library(stringr)
library("magrittr")
library(tidyverse)
################# read sample infomation

# CNV_info<-list.files("F:/胃癌预后/segments/") %>% str_split_fixed("\\.",3)
# sample_info<-read.table("F:/胃癌预后/TCGA-STAD.survival.tsv",sep = "\t",header = T)
# over_sample<-intersect(CNV_info[,1],sample_info[,3])
# 
# sample_name<-over_sample %>% paste0(".segments.txt") 
# b<-data.frame()
# for(i in 1:length(sample_name)){
#   a<-read.table(paste0("F:/胃癌预后/segments/",sample_name[i]),sep = "\t",header = T)
#   b<-rbind(a,b)
# }
# b%<>%.[-5]
# colnames(b)[2:5]<-c("chromosome", "start", "end", "segVal")
# write.table(b,"F:/胃癌预后/process_data/segTabs.tsv",sep = "\t",quote = F,row.names = F)
####################### sigminer process
segTabs<-read.table("F:/胃癌预后/process_data/segTabs.tsv",sep = "\t",header = T)
options(sigminer.sex = "male")


cn <- read_copynumber(segTabs,
                      seg_cols = c("chromosome", "start", "end", "segVal"),
                      genome_build = "hg19", complement = FALSE, verbose = TRUE
)


options(sigminer.copynumber.max = 20)
#cn_tally_W <- sig_tally(cn, method = "W")

cn_tally_W <- sig_tally(cn, method = "W",feature_setting = sigminer::CN.features[1:7,]) ### BP10MB

sample_components<-cn_tally_W[["nmf_matrix"]]  #### sample components

############################## 
library(survival) 
library(survminer) 


cox_data<-cbind(data.frame(sample_components),sample_info[match(over_sample,sample_info[,3]),c(2,4)])

############################# single-cox
covariates <-colnames(cox_data)[1:80]
univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(OS.time,OS)~', x)))
univ_models <- lapply( univ_formulas, function(x){coxph(x, data = cox_data)})

univ_results <- lapply(univ_models,
                       function(x){ 
                         x <- summary(x)
                         #获取p值
                         p.value<-signif(x$wald["pvalue"], digits=2)
                         #获取HR
                         HR <-signif(x$coef[2], digits=2);
                         #获取95%置信区间
                         HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
                         HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
                         HR <- paste0(HR, " (", 
                                      HR.confint.lower, "-", HR.confint.upper, ")")
                         res<-c(p.value,HR)
                         names(res)<-c("p.value","HR (95% CI for HR)")
                         return(res)
                       })
res <- t(as.data.frame(univ_results, check.names = FALSE))
res<-as.data.frame(res)
######################## multipe cox
cox_name<-res %>% filter(p.value<0.05) %>% rownames(.)
a<-paste('Surv(OS.time,OS)~', cox_name[1])
for(i in 2:length(cox_name)){
  a<-paste(a,"+", covariates[i])
}
cox_formul<-as.formula(a)

cox_model<-coxph(cox_formul,data = cox_data)
ggforest(cox_model,
         data=cox_data,
         main = "Hazard ratio",        # 标题
         cpositions = c(0.01, 0.1, 0.3), #左边三列的相对距离
         fontsize = 0.8, # 字体大小
         refLabel = "reference", #显示因子的参考水平
         noDigits = 3)  #小数点位数
 
################### KM-curve
km_data<-cox_data[,c(7,81,82)]

km_data["status"]<-ifelse(km_data[,1]>0,"height","low")

fit <- survfit(Surv(OS.time,OS)~ status, data = km_data)
ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, # 添加风险表
           risk.table.col = "strata", # 根据分层更改风险表颜色
           linetype = "strata", # 根据分层更改线型
           surv.median.line = "hv", # 同时显示垂直和水平参考线
           ggtheme = theme_bw(), # 更改ggplot2的主题
           palette = c("#E7B800", "#2E9FDF"))#定义颜色

#################### check absolute copynumber and relative copynumber
library(readxl)
relative_data<-read_xlsx("F:/胃癌预后/valiadata/patient.xlsx",col_names = T)
relative_sample<-relative_data%>%filter(...2=="TCGA")%>%.[,1]
relative_sample<-as.data.frame(relative_sample)

absolute_sample<-list.files("F:/胃癌预后/segments/") %>% str_split_fixed("\\.",3) %>% .[,1]
absolute_sample<-as.data.frame(absolute_sample)

over_sample<-intersect(relative_sample[,1],absolute_sample[,1])

sample_name<-over_sample %>% paste0(".segments.txt") 
absolute_data<-data.frame(sample="1",fragment_number=2)
for(i in 1:length(sample_name)){
  a<-read.table(paste0("F:/胃癌预后/segments/",sample_name[i]),sep = "\t",header = T)
  absolute_data[i,1]<-over_sample[i]
  absolute_data[i,2]<-nrow(a)
}
absolute_data%<>%mutate(type="absolute")

relative_data<-read.table("F:/胃癌预后/valiadata/GSE77775_cbs.seg.revised_160715.txt",sep = "\t",header = T)
relative_data <- relative_data[relative_data[,1]%in%over_sample,]
write.table(relative_data,file = "F:/胃癌预后/valiadata/realtive_data.tsv",sep = "\t",quote = F)

relative_data<-as.data.frame(table(relative_data[,1]))
relative_data%<>%mutate(type="relative")
##################### plot fragment density
library(ggsci)
library(ggplot2)
library(hrbrthemes)
library(ggprism)
colnames(relative_data)<-colnames(absolute_data)
density_data<-rbind(relative_data,absolute_data)

ggplot(data=density_data, aes(x=fragment_number, color=type,fill=type)) +
  geom_density(adjust=1.5, alpha=0.5,lwd=1.5,linetype = 1) +
  theme_ipsum()+
  ggtitle("density")+
  guides(color = guide_legend(title = "type"))+theme_prism()

var.test(fragment_number ~ type, density_data, alternative = "two.sided")

######################

segTabs <- read.table("F:/胃癌预后/valiadata/realtive_data.tsv",sep = "\t")
segTabs%<>%.[,-6]
colnames(segTabs)<- c("sample","chromosome", "start", "end", "segVal")

options(sigminer.sex = "male")

cn <- read_copynumber(segTabs,
                      seg_cols = c("chromosome", "start", "end", "segVal"),
                      genome_build = "hg19", complement = FALSE, verbose = TRUE
)

options(sigminer.copynumber.max = 20)

cn_tally_W <- sig_tally(cn, method = "W",feature_setting = sigminer::CN.features[1:7,]) ### BP10MB

realtive_components<-cn_tally_W[["nmf_matrix"]]  #### sample components

################# 
absoluta_data <- sample_components%>% as.data.frame(.)%>% mutate(sample=rownames(sample_components),type="absolute")%>%.[,7:9]
realtive_data <- realtive_components%>% as.data.frame(.)%>% mutate(sample=rownames(realtive_components),type="relative")%>%.[,7:9]

over_sample<-intersect(absoluta_data[,2],realtive_data[,2])

absoluta_data<-absoluta_data[absoluta_data[,2]%in%over_sample,]
realtive_data<-realtive_data[realtive_data[,2]%in%over_sample,]

######################
plot_data<-rbind(absoluta_data,realtive_data)

library(ggsci)
library(ggplot2)
library(hrbrthemes)
library(ggprism)
library(ggsignif)

colnames(plot_data)[1]<-"BP10MB"

ggplot(data=plot_data, aes(x=BP10MB, color=type,fill=type)) +
  geom_density(adjust=1.5, alpha=0.5,lwd=1.5,linetype = 1) +
  theme_ipsum()+
  ggtitle("density")+
  guides(color = guide_legend(title = "type"))+theme_prism()

var.test(BP10MB ~ type, plot_data, alternative = "two.sided")

ggplot(plot_data, aes(x = type, y = BP10MB)) +
  geom_boxplot() +
  labs(x = "Group", y = "Value", title = "Boxplot") +
  theme_minimal()+geom_signif(comparisons = c("absolute", "relative"),
                              map_signif_level = TRUE,
                              textsize = 4, vjust = -0.8)

