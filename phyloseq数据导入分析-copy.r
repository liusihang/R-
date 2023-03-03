#phyloseq数据导入分析
library(phyloseq)
library(stringr)

path <- getwd()
file_names <- list.files(path, pattern = "_counts_matched.csv$")
file_prefixes <- str_extract(file_names, ".*(?=_counts_matched.csv)")
#otu_list <- list.files(path, pattern = "_counts_matched.csv$", full.names = TRUE)
#sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1) #获取ssr名

# 遍历文件列表
for (i in seq_len(length(file_prefixes)))  {
  # 进行操作，例如读取文件内容
  tempotu <- paste0(file_prefixes[i],'_counts.csv',collapse='')
  temptaxa <- paste0(file_prefixes[i],'_taxa.csv',collapse='')
  tempsample <- paste0(file_prefixes[i],'_counts_matched.csv',collapse='')
  otucsv <- read.table(tempotu, header = TRUE, row.names = 1, sep = ",")
  taxacsv <- read.table(temptaxa, header = TRUE, row.names = 1, sep = ",")
  samplecsv <- read.table(tempsample, header = TRUE, row.names = 1, sep = ",")
  rownames(otucsv) <- samplecsv[,1]
  rownames(samplecsv) <- samplecsv[,1]
  otucsv <- data.frame(t(otucsv))
  for (x in seq_len(length(rownames(samplecsv))))  {
    # 进行操作，例如读取文件内容
    samplecsv[x,3] <- gsub(paste0('_',samplecsv[x,1]), "", samplecsv[x,3])
    print(x)
  }# 进行其他操作
  #
  #phyloseq read
  OTU = otu_table(otucsv, taxa_are_rows = TRUE)
  TAX = tax_table(as.matrix(taxacsv))
  SAM = sample_data(samplecsv)
  temp_phyloseq <- phyloseq(OTU, TAX, SAM)
  assign(file_prefixes[i], temp_phyloseq)
  # 进行其他操作
}

