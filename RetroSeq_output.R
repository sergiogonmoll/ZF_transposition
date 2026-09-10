library(tidyverse)
setwd("C:/Users/P309374/Documents/PhD/Macquarie/ZF_scripts/")
arguments<- commandArgs(TRUE)
arguments

filename_d<- paste0(arguments[1],"RetroSeq/SRR304140", arguments[1], ".vcf")
filename_s<- paste0(arguments[2],"RetroSeq/SRR304140", arguments[2], ".vcf")
filename_f<- paste0(arguments[3],"RetroSeq/SRR304140", arguments[3], ".vcf")

filename_d
filename_s
filename_f

focal_table <- read.table(filename_f)
focal_table <- subset(focal_table, !str_detect(focal_table$V8, "TEs"))
focal_table$GQ <- unlist(str_split(focal_table$V10, ":"))[seq(2,length(unlist(str_split(focal_table$V10, ":"))),by = 6)]
focal_table$FL <- unlist(str_split(focal_table$V10, ":"))[seq(3,length(unlist(str_split(focal_table$V10, ":"))),by = 6)]
focal_table <- subset(focal_table, focal_table$FL == 8 & focal_table$GQ >= 20)

dam_table <- read.table(filename_d)
dam_table <- subset(dam_table, !str_detect(dam_table$V8, "TEs"))
dam_table$GQ <- unlist(str_split(dam_table$V10, ":"))[seq(2,length(unlist(str_split(dam_table$V10, ":"))),by = 6)]
dam_table$FL <- unlist(str_split(dam_table$V10, ":"))[seq(3,length(unlist(str_split(dam_table$V10, ":"))),by = 6)]
dam_table <- subset(dam_table, dam_table$FL == 8 & dam_table$GQ >= 20)

sire_table <- read.table(filename_s)
sire_table <- subset(sire_table, !str_detect(sire_table$V8, "TEs"))
sire_table$GQ <- unlist(str_split(sire_table$V10, ":"))[seq(2,length(unlist(str_split(sire_table$V10, ":"))),by = 6)]
sire_table$FL <- unlist(str_split(sire_table$V10, ":"))[seq(3,length(unlist(str_split(sire_table$V10, ":"))),by = 6)]
sire_table <- subset(sire_table, sire_table$FL == 8 & sire_table$GQ >= 20)

chromosomes<- unique(focal_table_het$V1)
n_shared_dam<- c()
n_shared_sire<- c()
for (chr in chromosomes) {
chr_set <- subset(focal_table, focal_table$V1==chr)  
dam_chr_set <-subset(dam_table, dam_table$V1==chr)
sire_chr_set <-subset(sire_table, sire_table$V1==chr)
  for (p in dam_chr_set$V2) {
    shared_dam<- subset(chr_set, (chr_set$V2 >= (p - 300)) & 
                          (chr_set$V2 <= (p + 300)))
    n_shared_dam <- rbind(n_shared_dam, shared_dam)
  }
  for (q in sire_chr_set$V2) {
    shared_sire<- subset(chr_set, (chr_set$V2 >= (q - 300)) & 
                           (chr_set$V2 <= (q + 300)))
    n_shared_sire <- rbind(n_shared_sire, shared_sire)
  }
}

inherited_TEs <- rbind(n_shared_dam,n_shared_sire)
private_reliable <- anti_join(focal_table, inherited_TEs)



private_reliableSample <- arguments[3]
write.csv(private_reliable, paste0(arguments[3],"private_TEs_detail.csv"))
