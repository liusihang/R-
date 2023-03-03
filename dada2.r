#dada2分析
#设置工作目录获取文件夹名
#放在最上级文件夹运行
path <- getwd()
setwd(path)
subdirs <- list.dirs(".", recursive = FALSE)
library(dada2)
for (i in subdirs) {
    tempdir <- i
    newdir <- gsub("\\.", "", tempdir)
    workdir <- paste(path, '/',newdir, sep = "")
    setwd(workdir)
    '#filtr data
    fnFs <- sort(list.files(".", pattern="_1.fastq.gz$", full.names = TRUE))
    fnRs <- sort(list.files(".", pattern="_2.fastq.gz$", full.names = TRUE))
    sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1) #获取ssr名
    plotQualityProfile(fnFs[1:2])
    plotQualityProfile(fnRs[1:2])
    filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
    filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
    names(filtFs) <- sample.names
    names(filtRs) <- sample.names
    out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,160),
                         maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                         compress=TRUE, multithread=TRUE)'
    #wai bu guo lv liucheng
    filtFs <- sort(list.files(".", pattern="_1.fastq.gz_paired.fq.gz", full.names = TRUE))
    filtRs <- sort(list.files(".", pattern="_2.fastq.gz_paired.fq.gz", full.names = TRUE))
    sample.names <- sapply(strsplit(basename(filtFs), "_"), `[`, 1) #获取ssr名
    errF <- learnErrors(filtFs, multithread=TRUE)
    errR <- learnErrors(filtRs, multithread=TRUE)
    derepFs <- derepFastq(filtFs, verbose=TRUE)
    derepRs <- derepFastq(filtRs, verbose=TRUE)
    #dadaFs <- dada(derepFs, err=errF, multithread=TRUE,pool="pseudo")
    dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
    #dadaRs <- dada(derepRs, err=errR, multithread=TRUE,pool="pseudo")
    dadaRs <- dada(derepRs, err=errR, multithread=TRUE)
    dadaFs[[1]]
    mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
    head(mergers[[1]])
    seqtab <- makeSequenceTable(mergers)
    dim(seqtab)
    table(nchar(getSequences(seqtab)))
    seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
    sum(seqtab.nochim)/sum(seqtab)
    taxa <- assignTaxonomy(seqtab.nochim, "/media/em/student/LSH/others/refdatabase/silva_nr99_v138_train_set.fa.gz", multithread=TRUE)
    taxa <- addSpecies(taxa, '/media/em/student/LSH/others/refdatabase/silva_species_assignment_v138.fa.gz')
    namedir <- gsub('/','',newdir)
    taxa.otu <- taxa
    rownames(taxa.otu) <- paste0('OTU',1:nrow(taxa.otu))
    #写taxa
    taxacsv <- paste(namedir,'_1taxa.csv',seq='')
    taxaonlyotucsv <- paste(namedir,'_1taxaOTU.csv',seq='')
    setwd('/media/em/student/LSH/others/nofilterdata')
    write.csv(taxa,file = taxacsv,append = FALSE, quote = FALSE , sep = " ",eol = "\n", na = "NA", dec = ".", row.names = TRUE,col.names = TRUE, qmethod = c("escape", "double"),fileEncoding = "")
    write.csv(taxa.otu,file = taxaonlyotucsv,append = FALSE, quote = FALSE , sep = " ",eol = "\n", na = "NA", dec = ".", row.names = TRUE,col.names = TRUE, qmethod = c("escape", "double"),fileEncoding = "")
    seqtable.taxa.plus <- cbind('#seq'=rownames(taxa), t(seqtab.nochim), taxa)

    # ASV表格导出
    asvcounts <- paste(namedir,'_1counts.csv',seq='')
    write.csv(seqtab.nochim, file =asvcounts,append = FALSE, quote = FALSE , sep = " ",eol = "\n", na = "NA", dec = ".", row.names = TRUE,col.names = TRUE, qmethod = c("escape", "double"),fileEncoding = "")
    #write.table(seqtab.nochim, file =asvcounts, sep="\t", quote=F, row.names = T)
    # 带注释文件的ASV表格导出
    countstaxonasv <- paste(namedir,'_1counts.taxon.species.txt',seq='')
    write.csv(seqtable.taxa.plus , file = countstaxonasv,append = FALSE, quote = FALSE , sep = " ",eol = "\n", na = "NA", dec = ".", row.names = TRUE,col.names = TRUE, qmethod = c("escape", "double"),fileEncoding = "")
    #write.table(seqtable.taxa.plus , file = countstaxonasv, sep="\t", quote=F, row.names = F)
    # track文件保存
    #trackcsv <- paste(namedir,'_track.txt',seq='')
    #write.table(track, file = trackcsv, sep="\t", quote=F, row.names = F)
    setwd(path)
}

