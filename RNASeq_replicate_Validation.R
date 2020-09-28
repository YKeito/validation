"~/Nakano_RNAseq/network_analysis/script/validation/RNASeq_replicate_Validation.R"
DirectoryInfo <- read.table("~/bigdata/yasue/RNASeq/NakanoRNASeq_Info.txt", sep = "\t", header = T, stringsAsFactors = F)
pwd <-"~/bigdata/yasue/RNASeq/CopyOf035_nakano_19021/cufflinks_results/"
cufflinks <- c()
i <- 1
for(i in i:nrow(DirectoryInfo)){
  title <- paste0(pwd, DirectoryInfo$directory[i], "_outputs/genes.fpkm_tracking")
  T.FPKM <- read.table(title, sep = "\t", fill = T, header = T, stringsAsFactors = F, quote = "")
  T.FPKM <- T.FPKM %>% select(gene_id, FPKM) %>% distinct(gene_id, .keep_all = T)
  colnames(T.FPKM)[2] <- DirectoryInfo$sample[i]
  if(i == 1){
    cufflinks <- T.FPKM
  }else{
    cufflinks <- cufflinks %>% left_join(y = T.FPKM, by = "gene_id")
  }
  print(i)
  i <- i+1
}

pwd <- "~/bigdata/yasue/RNASeq/CY_flg22/cufflinks/"
T.sample <- list.files(pwd)[grep("_cufflinks", list.files(pwd))]
i <- 1
for(i in i:length(T.sample)){
  title <- paste0(pwd, T.sample[i], "/genes.fpkm_tracking")
  T.FPKM <- read.table(title, sep = "\t", fill = T, header = T, stringsAsFactors = F, quote = "")
  T.FPKM <- T.FPKM %>% select(gene_id, FPKM) %>% distinct(gene_id, .keep_all = T)
  temp <- str_split(T.sample[i], pattern = "_", simplify = T)
  colnames(T.FPKM)[2] <- paste0(temp[, 1], "_48h_", temp[, 2])
  cufflinks <- cufflinks %>% left_join(y = T.FPKM, by = "gene_id")
  print(i)
  i <- i+1
}
library(Hmisc)
T.data <- cufflinks %>% select(-matches("CY3|CY22|CY55"), -gene_id)
rownames(T.data) <- cufflinks$gene_id
temp <- T.data %>% as.matrix() %>% rcorr()
saveRDS(object = cufflinks, file = "~/bigdata/yasue/RNASeq/CY_flg22/cufflinks/cufflinks.rds")

#load library----
library(dplyr)
library(ggplot2)
library(stringr)
library(tidyr)
library(gridExtra)
#input data----
cufflinks <- readRDS(file = "~/bigdata/yasue/RNASeq/CY_flg22/cufflinks/cufflinks.rds")
#1~24h_duplicate----
T.sample <- colnames(cufflinks)[2:ncol(cufflinks)] %>% str_split(pattern = "_rep", simplify = T)
T.sample <- T.sample[, 1] %>% unique()
temp <- str_split(T.sample, pattern = "_", simplify = T)
T.Time <- temp[, 2] %>% unique()
T.CY <- temp[, 1] %>% unique()
T.CY <- T.CY %>% setdiff(c("CY22", "CY3", "CY55"))
T.rep <- c("rep1", "rep2", "rep3")
i <- 1
for(i in i:length(T.CY)){
  j <- 1
  for(j in j:length(T.Time)){
    if(j < length(T.Time)){
      rep1.data <- cufflinks %>% select("gene_id", paste0(T.CY[i], "_", T.Time[j], "_rep1"))
      rep1.data <- rep1.data %>% gather(key = Sample, value = FPKM_rep1, -gene_id)
      rep2.data <- cufflinks %>% select("gene_id", paste0(T.CY[i], "_", T.Time[j], "_rep2"))
      rep2.data <- rep2.data %>% gather(key = Sample, value = FPKM_rep2, -gene_id)
      df <- rep1.data %>% left_join(y = rep2.data, by = "gene_id")
      g <- ggplot(data = df, aes(x = FPKM_rep1, y = FPKM_rep2))
      g <- g + geom_point()
      g <- g + geom_smooth(method = "lm")
      a <- df %>% lm(formula = FPKM_rep1~FPKM_rep2)
      a <- summary(a)
      g <- g + ggtitle(paste0(T.CY[i], "_", T.Time[j], ":", "R^2=", formatC(a$adj.r.squared, digits = 3)))
      g <- g + theme(axis.text = element_text(size=10))
      g <- g + theme(axis.title = element_text(size=10))
      assign(paste0("g", j), g)
    }else{
      k <- 1
      for(k in k:2){
        m <- k+1
        for(m in m:3){
          rep1.data <- cufflinks %>% select("gene_id", paste0(T.CY[i], "_48h_",  T.rep[k]))
          rep1.data <- rep1.data %>% gather(key = Sample, value = FPKM_rep1, -gene_id)
          rep2.data <- cufflinks %>% select("gene_id", paste0(T.CY[i], "_48h_", T.rep[m]))
          rep2.data <- rep2.data %>% gather(key = Sample, value = FPKM_rep2, -gene_id)
          df <- rep1.data %>% left_join(y = rep2.data, by = "gene_id")
          g <- ggplot(data = df, aes(x = FPKM_rep1, y = FPKM_rep2))
          g <- g + geom_point()
          g <- g + geom_smooth(method = "lm")
          a <- df %>% lm(formula = FPKM_rep1~FPKM_rep2)
          a <- summary(a)
          g <- g + ggtitle(paste0(T.CY[i], "_48h:", "R^2=", formatC(a$adj.r.squared, digits = 3)))
          g <- g + xlab(paste0("FPKM_", T.rep[k])) + ylab(paste0("FPKM_", T.rep[m]))
          g <- g + theme(axis.text = element_text(size=10))
          g <- g + theme(axis.title = element_text(size=10))
          assign(paste0("g", k, m), g)
          m <- m+1
        }
        k <- k+1
      }
    }
    j <- j+1
  }
  g <- grid.arrange(g1, g2, g3, g4, g12, g13, g23, nrow = 3)
  title <- paste0("~/bigdata/yasue/RNASeq/validation/Dierk/Check2_", T.CY[i], ".png")
  ggsave(filename = title, plot = g)
  rm(list = c("g1", "g2", "g3", "g4", "g12", "g13", "g23"))
  i <- i+1
}
