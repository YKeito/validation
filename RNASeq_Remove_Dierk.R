"~/Nakano_RNAseq/network_analysis/script/validation/RNASeq_Remove_Dierk.R"
#load library----
library(dplyr)
library(ggplot2)
library(stringr)
library(tidyr)
library(gridExtra)
library(Hmisc)
#input data----
cufflinks <- readRDS(file = "~/bigdata/yasue/RNASeq/CY_flg22/cufflinks/cufflinks.rds")
#check area-----
#DMSO
cufflinks %>% filter(DMSO_48h_rep1 >=10000, DMSO_48h_rep2 <=10000)  #ATCG00320
#CY15
cufflinks %>% filter(CY15_48h_rep1 >=20000, CY15_48h_rep2 <= 1000)  #AT1G09467
cufflinks %>% filter(CY15_48h_rep1 >=20000, CY15_48h_rep3 <= 1000)  #AT1G09467
#CY16
cufflinks %>% filter(CY16_48h_rep1 >= 40000)                        #ATCG00010
cufflinks %>% filter(CY16_48h_rep1 <= 1000, CY16_48h_rep2 >= 10000) #AT4G08755
cufflinks %>% filter(CY16_48h_rep2 <= 1000, CY16_48h_rep3 >= 60000) #ATCG00010
cufflinks %>% filter(CY16_48h_rep2 <= 1000, CY16_48h_rep3 >= 60000) #ATCG00010
cufflinks %>% filter(CY16_48h_rep1 <= 1000, CY16_48h_rep2 >= 10000) #AT4G08755
cufflinks %>% filter(CY16_48h_rep1 >= 7000, CY16_48h_rep2 <= 1000) #AT2G01021
cufflinks %>% filter(CY16_48h_rep2 <= 1000, CY16_48h_rep3 >= 9000) #AT2G01021
cufflinks %>% filter(CY16_48h_rep2 >= 12000, CY16_48h_rep3 <= 1000) #AT4G08755
#CY20
cufflinks %>% filter(CY20_3h_rep1 >= 30000)                        #AT1G07600
cufflinks %>% filter(CY20_12h_rep1 >= 30000)                       #AT1G07600
cufflinks %>% filter(CY20_48h_rep1 >= 15000, CY20_48h_rep2 >= 10000) #ATCG00010
cufflinks %>% filter(CY20_48h_rep1 >= 15000, CY20_48h_rep3 >= 20000) #ATCG00010
cufflinks %>% filter(CY20_48h_rep2 >= 10000, CY20_48h_rep3 >= 25000) #ATCG00010
#remove strange point----
T.AGI <- c("ATCG00320", "AT1G09467", "ATCG00010", "AT4G08755", "AT2G01021", "AT1G07600")
cufflinks <- cufflinks[-match(T.AGI, cufflinks$gene_id), ]
saveRDS(object = cufflinks, file = "~/bigdata/yasue/RNASeq/CY_flg22/cufflinks/modified_cufflinks.rds")
cuffdiff <- readRDS(file = "~/Nakano_RNAseq/network_analysis/.RData/RDS/allCY_cuffdiff.rds")
cuffdiff <- cuffdiff[-match(T.AGI, cuffdiff$AGI), ]
saveRDS(object = cuffdiff, "~/Nakano_RNAseq/network_analysis/.RData/RDS/allCY_cuffdiff_modified.rds")

T.data <- cufflinks %>% select(-matches("CY3|CY22|CY55"), -gene_id)
rownames(T.data) <- cufflinks$gene_id
temp <- T.data %>% as.matrix() %>% rcorr()
write.table(x = temp$r, file = "~/bigdata/yasue/RNASeq/MorohashiCheck_validation_Dierk.txt", sep = "\t", quote = F, row.names = T)