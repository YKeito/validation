"~/Nakano_RNAseq/network_analysis/script/validation/RNASeq_Tophat2_MappedInfo.R"
library(stringr)
library(dplyr)
pwd <- "/home/yasue/bigdata/yasue/RNASeq/CopyOf035_nakano_19021/"
sample <- read.table(file = "/home/yasue/bigdata/yasue/RNASeq/CopyOf035_nakano_19021/RNA-seq.txt", sep = "\t", skip = 1, row.names = 1, stringsAsFactors = F)
title <- list.files(path = pwd)
T.condition <- str_split(sample$V3, pattern = "-", simplify = T)
T.condition <- paste0(T.condition[, 2], "_", str_split(T.condition[, 3], pattern = "_", simplify = T)[, 1])
Sequence.summary <- c()
i <- 1
for(i in i:nrow(sample)){
  T.data <- read.table(file = paste0(pwd, title[grep(sample$V2[i], title)], "/tophat_out/align_summary.txt"),
                       sep = ":", fill = T, stringsAsFactors = F)
  temp <- gsub(T.data$V2[3], pattern = " ", replacement = "")
  temp <- str_split(temp, pattern = "[(]", simplify = T)
  Sequence.summary <- rbind(Sequence.summary,
                            data.frame(condition = T.condition[i],
                                       Seq_file_name = sample$V2[i],
                                       Input = gsub(T.data$V2[2], pattern = " ", replacement = ""),
                                       Mapped = temp[, 1],
                                       rate = str_split(temp[, 2], pattern = "of", simplify = T)[, 1],
                                       stringsAsFactors = F)
  )
}

write.table(x = Sequence.summary, 
            file = "~/bigdata/yasue/RNASeq/CopyOf035_nakano_19021/tophat2_summary.txt", 
            sep = "\t", quote = F, row.names = F)