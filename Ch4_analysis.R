library(tidyverse)
library(lme4)
library(performance)
library(DHARMa)
library(ggplot2)
library(glmmTMB)
library(stringr)
library(data.table)
library(lmerTest)
library(ggpubr)
library(ggeffects)

## TE landscape plot ###
TE_landscape<- read.csv("TE_landscape.csv", check.names =FALSE)
TE_landscape<- reshape2::melt(TE_landscape, id.vars=
                                "Div")
TE_landscape$value <- (TE_landscape$value/1053687097)*100
for (e in 1:length(TE_landscape$variable)) {
  TE_landscape$Class[e]<- str_split(TE_landscape$variable, "/")[[e]][1]
}

TE_ls_plot<- ggplot(TE_landscape, aes(x=Div, y=value, fill=Class))+geom_bar(position = "stack",stat = "identity", colour="black")+
  theme_bw()+labs(y="% of genome", x="Kimura substitution level (CpG adjusted)", fill="TE order", title = "Whole genome")+
  scale_x_continuous(limits = c(0,60), breaks = seq(0,60,5))+
  scale_y_continuous(limits = c(0,0.33), breaks = seq(0,0.35,0.1))

TE_landscape_macro<- read.csv("TE_landscape_macro.csv", check.names =FALSE)
TE_landscape_macro<- reshape2::melt(TE_landscape_macro, id.vars=
                                      "Div")
TE_landscape_macro$value <- (TE_landscape_macro$value/1053687097)*100
for (e in 1:length(TE_landscape_macro$variable)) {
  TE_landscape_macro$Class[e]<- str_split(TE_landscape_macro$variable, "/")[[e]][1]
}
TE_ls_plot_macro<- ggplot(TE_landscape_macro, aes(x=Div, y=value, fill=Class))+geom_bar(position = "stack",stat = "identity", colour="black", size = 0.05)+
  theme_bw()+labs(y="% of genome", x="Kimura substitution level (CpG adjusted)", fill="TE order", title = "Macro-chromosomes")+
  scale_x_continuous(limits = c(0,60), breaks = seq(0,60,5))+
  scale_y_continuous(limits = c(0,0.15), breaks = seq(0,0.35,0.1))

TE_landscape_W<- read.csv("TE_landscape_W.csv", check.names =FALSE)
TE_landscape_W<- reshape2::melt(TE_landscape_W, id.vars=
                                  "Div")
TE_landscape_W$value <- (TE_landscape_W$value/1053687097)*100
for (e in 1:length(TE_landscape_W$variable)) {
  TE_landscape_W$Class[e]<- str_split(TE_landscape_W$variable, "/")[[e]][1]
}
TE_ls_plot_W<- ggplot(TE_landscape_W, aes(x=Div, y=value, fill=Class))+geom_bar(position = "stack",stat = "identity", colour="black", size = 0.05)+
  theme_bw()+labs(y="% of genome", x="Kimura substitution level (CpG adjusted)", fill="TE order", title = "W chromosome")+
  scale_x_continuous(limits = c(0,60), breaks = seq(0,60,5))+
  scale_y_continuous(limits = c(0,0.15), breaks = seq(0,0.35,0.1))

TE_landscape_micro<- read.csv("TE_landscape_micro.csv", check.names =FALSE)
TE_landscape_micro<- reshape2::melt(TE_landscape_micro, id.vars=
                                      "Div")
TE_landscape_micro$value <- (TE_landscape_micro$value/1053687097)*100
for (e in 1:length(TE_landscape_micro$variable)) {
  TE_landscape_micro$Class[e]<- str_split(TE_landscape_micro$variable, "/")[[e]][1]
}
TE_ls_plot_micro<- ggplot(TE_landscape_micro, aes(x=Div, y=value, fill=Class))+geom_bar(position = "stack",stat = "identity", colour="black", size = 0.05)+
  theme_bw()+labs(y="% of genome", x="Kimura substitution level (CpG adjusted)", fill="TE order", title = "Micro-chromosomes")+
  scale_x_continuous(limits = c(0,60), breaks = seq(0,60,5))+
  scale_y_continuous(limits = c(0,0.15), breaks = seq(0,0.35,0.1))

TE_landscape_fig<- ggarrange(TE_ls_plot, TE_ls_plot_macro, TE_ls_plot_micro, TE_ls_plot_W, common.legend = T, nrow = 1, legend = "bottom", label.x = "Kimura substitution level (CpG adjusted)")
ggsave("Figure1_TE_landscape.png", TE_landscape_fig, width = 8, height = 3, dpi = 'retina', scale = 2)


final_mapped <- read.csv("TEs_ZF.csv")
ZF_metadata <- read.csv("TEs_ZF_metadata.csv")
final_mapped$Sample <- as.factor(final_mapped$Sample)
ZF_metadata$offspring <- as.factor(ZF_metadata$offspring)


final_mapped <- left_join(final_mapped, ZF_metadata, by=c("Sample"="offspring"))

### Check sex chromosomes ####
final_mapped_Z <- subset(final_mapped, final_mapped$V1=="NC_044241.2")
final_mapped_W <- subset(final_mapped, final_mapped$V1=="NC_045028.1")
dim(final_mapped_Z)

final_mapped_Z_ind <- final_mapped_Z %>% group_by(Sample) %>% summarise(counts = n(),sex = first(sex),gen=first(gen), length=first(length))
final_mapped_Z_ind$length[which(final_mapped_Z_ind$sex == "M")]<- final_mapped$length[1]*2


final_mapped_auto <- subset(final_mapped, final_mapped$V1!="NC_044241.2")
final_mapped_auto_ind <- final_mapped_auto %>% group_by(Sample,V1) %>% summarise(counts = n(),sex = first(sex),gen=first(gen), length=first(length))

final_mapped_auto_ind$dens <- final_mapped_auto_ind$counts/final_mapped_auto_ind$length*1000000
final_mapped_Z_ind$dens <- final_mapped_Z_ind$counts/final_mapped_Z_ind$length*1000000
final_mapped_auto_ind$type <- "Autosomal"
final_mapped_Z_ind$type <- "Z"
final_mapped_chr_typ <- rbind(final_mapped_auto_ind,final_mapped_Z_ind)
ggplot(final_mapped_chr_typ, aes(x=type, y=dens))+geom_boxplot(aes(color=sex))+
  stat_compare_means(method = "t.test", data = final_mapped_chr_typ, aes(group = sex), label = "p.signif")+geom_pwc(
    method = "t_test", label = "p.signif",
    bracket.nudge.y = 0.1
  )+xlab("Chromosome type")+ylab("Insertions per Mb")+ggtitle("Z chromosome NTEI density")+theme_bw()

final_mapped_auto_ind_micro <- subset(final_mapped_auto_ind, final_mapped_auto_ind$length < 40000000)
final_mapped_auto_ind_macro <- subset(final_mapped_auto_ind, final_mapped_auto_ind$length > 40000000)
final_mapped_chr_typ_micro <- rbind(final_mapped_auto_ind_micro, final_mapped_Z_ind)
final_mapped_chr_typ_macro <- rbind(final_mapped_auto_ind_macro, final_mapped_Z_ind)
final_mapped_chr_typ_micro$chr_size <- "Micro"
final_mapped_chr_typ_macro$chr_size <- "Macro"
final_mapped_chr_typ <- rbind(final_mapped_chr_typ_micro, final_mapped_chr_typ_macro)
Z_chr_plot<- ggplot(final_mapped_chr_typ, aes(x=type, y=dens))+geom_boxplot(aes(color=sex))+
  stat_compare_means(method = "t.test", data = final_mapped_chr_typ, aes(group = sex), label = "p.signif")+geom_pwc(
    method = "t_test", label = "p.signif",
    bracket.nudge.y = 0.1
  )+xlab("Chromosome type")+ylab("Insertions per Mbp")+ggtitle("Z chromosome NTEI density")+theme_bw()+
  facet_wrap(~chr_size)


#### NTEIs compared between families and generations ####
mapped_summary <- final_mapped %>% group_by(Sample) %>% summarise(LTR=sum(Class == "LTR"),LINE=sum(Class == "LINE"), Unknown=sum(Class == "Unknown"))

mapped_summary<- reshape2::melt(mapped_summary, id.vars=
                                  "Sample")
mapped_summary <- left_join(mapped_summary, ZF_metadata, by=c("Sample"="offspring"))

mapped_summary_1 <- subset(mapped_summary, mapped_summary$gen==1)
mapped_summary_2 <- subset(mapped_summary, mapped_summary$gen==2)

mapped_summary_1$gen_label <- "F1"
mapped_summary_2$gen_label <- "F2"
mapped_summary_1$fam_label <- NA
mapped_summary_2$fam_label <- NA
mapped_summary_1$fam_label[which(mapped_summary_1$fam==1)]<- "Fam 1"
mapped_summary_1$fam_label[which(mapped_summary_1$fam==2)]<- "Fam 2"
mapped_summary_1$fam_label[which(mapped_summary_1$fam==3)]<- "Fam 3"
mapped_summary_1$fam_label[which(mapped_summary_1$fam==4)]<- "Fam 4"

mapped_summary_2$fam_label[which(mapped_summary_1$fam==1)]<- "Fam 1"
mapped_summary_2$fam_label[which(mapped_summary_1$fam==2)]<- "Fam 2"
mapped_summary_2$fam_label[which(mapped_summary_1$fam==3)]<- "Fam 3"
mapped_summary_2$fam_label[which(mapped_summary_1$fam==4)]<- "Fam 4"

#### Plot results #####
TE_map_1_plot <- ggplot(mapped_summary_1, aes(x=Sample, y=value, fill=variable))+geom_bar(position = "stack",stat = "identity")+
  facet_grid(gen_label~fam_label, scales = "free_x")+theme_bw()+scale_y_continuous(limits = c(0,60), breaks = seq(0,60,20))+
  labs(y = "TE insertion count", x = "Sample", fill="TE Class")

TE_map_2_plot <- ggplot(mapped_summary_2, aes(x=Sample, y=value, fill=variable))+geom_bar(position = "stack",stat = "identity")+
  facet_grid(gen_label~fam_label, scales = "free_x")+theme_bw()+scale_y_continuous(limits = c(0,60), breaks = seq(0,60,20))+
  labs(y = "TE insertion count", x = "Sample", fill="TE Class")

TE_map_plot <- ggarrange(TE_map_1_plot, TE_map_2_plot, nrow = 2, ncol = 1, common.legend = TRUE, legend = "right", label.y = "Generation")
target_map_chr <- final_mapped %>% group_by(Sample,V1) %>% summarise(insertions=n(), length = first(length), sex = first(sex), gen=first(gen), fam = first(fam),
                                                                     sire = first(sire), dam=first(dam))

ggsave("Figure2_TE_summary.png", TE_map_plot, dpi = 'retina', scale = 2, width = 4, height = 3)
ggsave("FigureS1_Z_chr_plot.png", Z_chr_plot, dpi = 'retina', scale = 2, width = 4, height = 3)


#### Running models for NTEI densities ####
target_map_chr$ins_per_mb <- (target_map_chr$insertions/target_map_chr$length)*1000000
target_map_chr$V1<- as.factor(target_map_chr$V1)
target_map_chr$sex<- as.factor(target_map_chr$sex)
target_map_chr$gen<- as.factor(target_map_chr$gen)
target_map_chr$fam<- as.factor(target_map_chr$fam)
target_map_chr$dam<- as.factor(target_map_chr$dam)
target_map_chr$sire<- as.factor(target_map_chr$sire)


target_map_chr$log_size <- log(target_map_chr$length)

ins_model_map<-glmmTMB(ins_per_mb ~ log_size+fam+gen+sex*log_size+(1|sire)+(1|dam),
                       data = target_map_chr, REML = T,
                   family=Gamma(link = "log"),control = glmmTMBControl(rank_check = "adjust"))

target_map_chr_ind <- target_map_chr %>% group_by(Sample) %>% summarise(t.insertions=sum(insertions), ins_per_mb=t.insertions/sum(unique(final_mapped$length)), gen=first(gen), sex=first(sex), dam=first(dam), sire=first(sire), fam=first(fam))
ins_model_map_ind<-glmmTMB(t.insertions ~gen+sex+fam+(1|dam)+(1|sire), data = target_map_chr_ind, REML = T,
                       family="poisson",control = glmmTMBControl(rank_check = "adjust"))


summary(ins_model_map_ind)
diagnose(ins_model_map_ind)
plot(simulateResiduals(ins_model_map_ind))
check_singularity(ins_model_map_ind)
check_collinearity(ins_model_map)

#### Get model predictions ####
p.log_size_sex<- ggpredict(ins_model_map, terms = c("log_size", "sex"), bias_correction = T, ci_level = 0.95, type = "fixed", interval = "confidence")
p.sex<- ggpredict(ins_model_map, terms = c("sex"), bias_correction = T, ci_level = 0.95, type = "fixed", interval = "confidence")
p.gen<- ggpredict(ins_model_map, terms = c("gen"), bias_correction = T, ci_level = 0.95, type = "fixed", interval = "confidence")
p.gen_ind<- ggpredict(ins_model_map_ind, terms = c("gen"), bias_correction = T, ci_level = 0.95, type = "fixed", interval = "confidence")

ggplot(data=p.sex, aes(x=x,y=predicted))+geom_errorbar(aes(ymin=conf.low, ymax=conf.high), width=0.4)+geom_point(aes(x=x,y=predicted))+
  geom_violin(data = target_map_chr, aes(x=sex, y=ins_per_mb), alpha=0.2)+geom_jitter(data = target_map_chr, aes(x=sex, y=ins_per_mb), alpha=0.2)+
  labs(y="Predicted novel TE insertions per Mbp", x="Sex")+theme_bw()

gen_chr_plot<- ggplot(data=p.gen, aes(x=x,y=log(predicted)))+geom_errorbar(aes(ymin=log(conf.low), ymax=log(conf.high)),width=0.4)+geom_point(aes(x=x,y=log(predicted)))+
  geom_violin(data = target_map_chr, aes(x=gen, y=log(ins_per_mb)), alpha=0.2)+geom_jitter(data = target_map_chr, aes(x=gen, y=log(ins_per_mb)), alpha=0.2)+
  labs(y="Log of predicted novel TE insertions per Mbp", x="F")+theme_bw()

gen_ind_plot<- ggplot(data=p.gen_ind, aes(x=x,y=predicted))+geom_errorbar(aes(ymin=conf.low, ymax=conf.high),width=0.4)+geom_point(aes(x=x,y=predicted))+
  geom_violin(data = target_map_chr_ind, aes(x=gen, y=t.insertions), alpha=0.2)+geom_jitter(data = target_map_chr_ind, aes(x=gen, y=t.insertions), alpha=0.2)+
  labs(y="Predicted novel TE insertions", x="F")+theme_bw()

chr_size_plot<- ggplot()+
  geom_line(data = p.log_size_sex, aes(x=x, y = predicted, color=group)) + 
  geom_ribbon(data = p.log_size_sex, aes(x=x, ymin=conf.low, ymax=conf.high, fill=group, color=group), alpha = 0.2)+
  geom_jitter(data = target_map_chr, aes(x=log_size, y=ins_per_mb, color=sex))+
  labs(y="Predicted novel TE insertions per Mbp", x="Log chromosome size (bp)", color="Sex",fill="Sex") + theme_bw()

target_map_chr_F <- subset(target_map_chr, target_map_chr$sex=="F")
target_map_chr_M <- subset(target_map_chr, target_map_chr$sex=="M")
target_map_chr_ind_1 <- subset(target_map_chr_ind, target_map_chr_ind$gen==1)
target_map_chr_ind_2 <- subset(target_map_chr_ind, target_map_chr_ind$gen==2)
mean(target_map_chr_ind_1$t.insertions)/mean(target_map_chr_ind_2$t.insertions)
target_map_chr_ind_M <- subset(target_map_chr_ind, target_map_chr_ind$sex=="M")
target_map_chr_ind_F <- subset(target_map_chr_ind, target_map_chr_ind$sex=="F")
mean(target_map_chr_ind_F$ins_per_mb)/mean(target_map_chr_ind_M$ins_per_mb)
mean_all <- mean(c(mean(target_map_chr_ind_F$ins_per_mb),mean(target_map_chr_ind_M$ins_per_mb)))

sd(target_map_chr_ind$t.insertions)
 
plot_chr<- ggarrange(chr_size_plot,gen_chr_plot, labels = c("A","B"), ncol = 1)
#ggsave("TE_size_chr.png", chr_size_plot, dpi = 'retina', scale = 2)
ggsave("Figure4_Chr_size&gen.png", plot_chr, width = 4, height = 6, dpi = 'retina', scale = 2)
ggsave("Figure3_gen_ind_figure.png", gen_ind_plot, width = 4, height = 4, dpi = 'retina', scale = 2)

## Relative length plot to compare with recombination
chr_sizes<- read.csv("Chr_sizes.csv")
final_mapped <- inner_join(final_mapped, chr_sizes, by=c("V1"="chr"))
final_mapped$rel_pos <- final_mapped$V2/final_mapped$length_GB
final_mapped_macro<-subset(final_mapped, final_mapped$length_GB>40000000)
final_mapped_micro<-subset(final_mapped, final_mapped$length_GB<40000000)
for (p in 1:length(final_mapped_macro$V2)) {
  final_mapped_macro$acc[p]<- sum(final_mapped_macro$rel_pos < final_mapped_macro$rel_pos[p])/length(final_mapped_macro$V2)  
}

micro_chrs <- subset(chr_sizes, chr_sizes$length_GB<40000000)
macro_chrs <- subset(chr_sizes, chr_sizes$length_GB>40000000)


for (p in 1:length(final_mapped_micro$V2)) {
  final_mapped_micro$acc[p]<- sum(final_mapped_micro$rel_pos < final_mapped_micro$rel_pos[p])/length(final_mapped_micro$V2)  
}

macro_r_plot<- ggplot(data = final_mapped_macro, aes(x=rel_pos, y=acc))+
  geom_line()+geom_abline(intercept = 0, slope = 1, color="red")+
  labs(y="% Accumulated NTEIs", x="Relative position on macro-chromosome")+theme_bw()
micro_r_plot<- ggplot(data = final_mapped_micro, aes(x=rel_pos, y=acc))+
  geom_line()+geom_abline(intercept = 0, slope = 1, color="red")+
  labs(y="% Accumulated NTEIs", x="Relative position on micro-chromosome")+theme_bw()
acc_plot<- ggarrange(macro_r_plot,micro_r_plot,ncol = 2)
ggsave("Figure5_TE_acc_plot.png", acc_plot, width = 4, height = 2, dpi = 'retina', scale = 2)
ks.test(unique(final_mapped_macro$acc), "punif")
#?ks.test()

### Transposition rates per class ###
library(gggenomes)
mappable_bed <- read_bed("mappable.bed.gz")
mappable_bed$length <- mappable_bed$end- mappable_bed$start
mappable_male <- filter(mappable_bed, mappable_bed$seq_id !="NC_045028.1")
mappable_male <- sum(mappable_male$length)*2
all_size <- mappable_bed %>% group_by(seq_id) %>% summarise(length = sum(length)) 
W_size <- subset(all_size, all_size$seq_id== "NC_045028.1")
Z_size <- subset(all_size, all_size$seq_id== "NC_044241.2")
mappable_female <- mappable_male - (Z_size$length-W_size$length) 
mapped_summary_1$bp <- NA
mapped_summary_1$bp[which(mapped_summary_1$sex == "F")]<- mappable_female
mapped_summary_1$bp[which(mapped_summary_1$sex == "M")]<- mappable_male
mapped_summary_1$Rate <- mapped_summary_1$value/mapped_summary_1$bp
mapped_rates_class_1 <- mapped_summary_1 %>% group_by(variable) %>% summarise(t.value=sum(value), rate=mean(Rate))
mapped_rates_1 <- mapped_summary_1 %>% group_by(Sample) %>% summarise(t.value=sum(value), t.rate=sum(Rate))
mapped_summary_2$bp <- NA
mapped_summary_2$bp[which(mapped_summary_2$sex == "F")]<- mappable_female
mapped_summary_2$bp[which(mapped_summary_2$sex == "M")]<- mappable_male
mapped_summary_2$Rate <- mapped_summary_2$value/mapped_summary_2$bp
mapped_rates_class_2 <- mapped_summary_2 %>% group_by(variable) %>% summarise(t.value=sum(value),rate=mean(Rate))
mapped_rates_2 <- mapped_summary_2 %>% group_by(Sample) %>% summarise(t.value=sum(value), t.rate=sum(Rate))
mean(mapped_rates_2$t.rate)
sum(mapped_rates_1$t.value)

mean(c(mapped_rates_1$t.rate,mapped_rates_2$t.rate))

mapped_rates_1$copy_rate <- mapped_rates_1$t.value/560140
mapped_rates_2$copy_rate <- mapped_rates_2$t.value/560140
mean(mapped_rates_1$copy_rate)/mean(mapped_rates_2$copy_rate)
mapped_rates_class_1$copy_rate <- mapped_rates_class_1$t.value/560140
mapped_rates_class_2$copy_rate <- mapped_rates_class_2$t.value/560140


t.test(mapped_rates_2$copy_rate)
t.test(mapped_rates_2$t.rate)


