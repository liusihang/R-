#开口向下
ggtree(tree, branch.length="none")  +
  #geom = "text",size = 2,linesize = 0.1,nudge_x=0.5 调整字体大小
  geom_tiplab(angle=90,offset = -70,align = TRUE,linetype = 0,geom = "text",size = 2.3,family='mono',fontface=2,linesize = 0.7,nudge_x = -1.2) + 
  #画竖线label，以节点确定
  geom_cladelabel(node=72, label="VNG",color="red", 
                  angle = 0, hjust='center',offset=-12.6,offset.text=-0.5, align=TRUE,barsize=1, fontsize=5) + 
  geom_cladelabel(node=59, label="TG", angle = 0, hjust='center',
                  color="blue",offset=-12.6,offset.text=-0.5, align=TRUE) +
  xlim(NA, 10)+  layout_dendrogram()
#开口正常
ggtree(treex2, branch.length="none")  +
  #geom = "text",size = 2,linesize = 0.1,nudge_x=0.5 调整字体大小
  geom_tiplab(angle=0,offset = -70,align = TRUE,linetype = 2,geom = "text",size = 2,family='mono',fontface=1,nudge_x = 1) + 
  #画竖线label，以节点确定
  geom_cladelabel(node=114, label="VNG",color="red", 
                  angle = 90, hjust='center',offset=3.2,offset.text=0.7, align=TRUE,barsize=1, fontsize=5) + 
  geom_cladelabel(node=116, label="TG", angle = 90, hjust='center',
                  color="blue",offset=3.2,offset.text=0.7, align=TRUE) +
  geom_cladelabel(node=120, label="VNG",color="purple", 
                  angle = 90, hjust='center',offset=3.2,offset.text=0.7, align=TRUE,barsize=1, fontsize=5) + 
  geom_cladelabel(node=136, label="VNG",color="green", 
                  angle = 90, hjust='center',offset=3.2,offset.text=0.7, align=TRUE,barsize=1, fontsize=5) + 
  xlim(NA, 10)+ scale_color_continuous(low='darkgreen', high='red') +
  theme(legend.position="right")+
  geom_hilight(node=114, fill="gold") + 
  geom_hilight(node=116, fill="purple")


tempotu <- paste0('VNG_Anaerobic_Tank_Planktonic_counts.csv',collapse='')
temptaxa <- paste0('VNG_Anaerobic_Tank_Planktonic_taxa.csv',collapse='')
tempsample <- paste0('VNG_Anaerobic_Tank_Planktonic_counts_matched.csv',collapse='')
otucsv <- read.table(tempotu, header = TRUE, row.names = 1, sep = ",")
taxacsv <- read.table(temptaxa, header = TRUE, row.names = 1, sep = ",")
samplecsv <- read.table(tempsample, header = TRUE, row.names = 1, sep = ",")
rownames(otucsv) <- samplecsv[,1]
rownames(samplecsv) <- samplecsv[,1]
otucsv <- data.frame(t(otucsv))
#phyloseq read
OTU = otu_table(otucsv, taxa_are_rows = TRUE)
TAX = tax_table(as.matrix(taxacsv))
SAM = sample_data(samplecsv)
VNG_Anaerobic_Tank_Planktonic <- phyloseq(OTU, TAX, SAM)