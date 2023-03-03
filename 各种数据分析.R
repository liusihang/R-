#根据类别提取样本
AD_VNGplant <- subset_samples(AD_merged, Plantname=="VNG_Plant")
AD_VNGplant <- prune_taxa(taxa_sums(VNGplant) > 0, VNGplant)#对样本进行清理，去除多余OTU


# 基于bray-curtis距离进行计算 PERMANOVA分析
#计算PLANT对群落影响 PERMANOVA分析PERMANOVA分析
#x1 = prune_taxa(taxa_sums(all_merged) > 5, all_merged) 去除小于5次的otu
x2otu <- data.frame(t(otu_table(x2)))#转置下 对应样品id
View(x2otu)
x1sample <- data.frame(sample_data(x2))
plant.div <- adonis2(x2otu ~ Plantname, data = x2sample, permutations = 999, method="bray")

#进行两个群落间Procrustes分析
distance_AS_merged <- distance(ASprune, method = "bray")
distance_AD_merged<- distance(ADprune, method = "bray")
procrustes <- procrustes(distance_AD_merged, distance_AD_merged, scale = TRUE)
perm_results <- protest(distance_AD_merged, distance_AD_merged, permutations = 999)#带p值
procrustes_sos <- procrustes_results$SS
procrustes_rmse <- procrustes_results$MSE
p_value <- perm_results$signif[1]

#图
ADotu <- data.frame(t(otu_table(AD_merged)))
ADotu_dist <- vegdist(ADotu, method="bray", binary=F)
dune_pcoa <- cmdscale(ADotu_dist, k=3, eig=T)
dune_pcoa_points <- as.data.frame(dune_pcoa$points)
sum_eig <- sum(dune_pcoa$eig)
eig_percent <- round(dune_pcoa$eig/sum_eig*100,1)
colnames(dune_pcoa_points) <- paste0("PCoA", 1:3)
ADsample <- data.frame(sample_data(AD_merged))
dune_pcoa_result <- cbind(dune_pcoa_points, ADsample)
head(dune_pcoa_result)

ggplot(dune_pcoa_result, aes(x=PCoA1, y=PCoA2, color=Bioname)) +
      labs(x=paste("PCoA 1 (", eig_percent[1], "%)", sep=""),
                     y=paste("PCoA 2 (", eig_percent[2], "%)", sep="")) +
      geom_point(size=4) + stat_ellipse(level=0.6) +
      theme_classic()

#phyloseq自带pcoa图绘制
pcoax0 <- ordinate(x0, method = "PCoA")
plot_ordination(x0, pcoax0, color = "Plantname",shape  = 'Bioname')


#phyloseq导入流程
library("phyloseq"); packageVersion("phyloseq")
library("ggplot2"); packageVersion("ggplot2")
VNGotu <- read.table('VNG_Anaerobic_Tank_Planktonic _otu.csv', header = TRUE, row.names = 1, sep = ",")
VNGtaxa <- read.table('VNG_Anaerobic_Tank_Planktonic _taxa.csv', header = TRUE, row.names = 1, sep = ",")
#OTU可能需要转置，注意保持data。fram格式不变
VNGotu <- data.frame(t(VNGotu))
OTU = otu_table(otumat, taxa_are_rows = TRUE)
TAX = tax_table(as.matrix(VNGtaxa))

alpha_diversity <- estimate_richness(physeq, measures = c("Chao1", "Shannon")) #计算alpha_diversity


x2 = prune_taxa(names(sort(taxa_sums(physeq), TRUE))[1:200], physeq) #过滤taxa

#计算beta-diversity并画图
pcoax2 <- ordinate(x2, method = "PCoA")
plot_ordination(x2, pcoax2, color = "location") #第一个参数确定后面数据上的注释，第二个参数才起到决定性作用



# 计算 Bray-Curtis距离矩阵
distance_matrix <- distance(physeq, method = "bray")
# 使用 complete-linkage 方法聚类
tree <- hclust(as.dist(distance_matrix), method = "complete")
# 绘制聚类树
plot(tree, main = "Complete-linkage clustering")

#计算组间差异
library(vegan)
#type名指的是几因素几水平时对应的因素名，水平指的就是做的不同的样本，比如地点123，因素就是除了得到的环境数据以外其它的变量，比如处理工业地点等
adonis2res <- adonis2(distance_matrix ~ 'type名', data = sample_data(physeq))
summary(adonis2res)

#Procrustes分析可以用来比较两个微生物群落的相似性
procrustes_bc_env <- procrustes('群落1distance_matrix','群落2distance_matrix')




# 计算样本间距离矩阵
dist_matrix <- vegdist(t(otu_table(temp_phyloseq)), method = "bray")

# 使用complete-linkage clustering构建系统发育树
tree <- hclust(dist_matrix, method = "complete") #计算样本间tree
reverse_phyloseq <- phyloseq(t(otu_table(temp_phyloseq),as.phylo(tree)))
# 将系统发育树添加到phyloseq对象中
phy_tree(my_phyloseq) <- as.phylo(tree)


#hill数计算
AD_TVplantotu <- (otu_table(AD_TVplant))
AD_TGplantotuVector <- AD_TGplantotuReverse[,1]
dist_AD_VNGplant <- vegdist((otu_table(AD_VNGplant)), method = "bray")
AD_VNGplant_tree <- hclust(dist_AD_VNGplant, method = "complete")#计算otu tree
AD_VNGplant_tree <- as.phylo(AD_VNGplant_tree)
hill_div(AD_IGplantotu,1)#提供hill number及effective number of species
hill_div(AD_VNGplantotu,1,AD_VNGplant_tree)#加上了tree的信息，系统发育希尔数考虑到了这样一个事实，
#即密切相关的物种可能具有相似的生态特征，因此可能比更远的物种对群落多样性的贡献较小。
#小提琴图绘制
#testggpubr为dataframe
AD_TVplanhill <- data.frame(hill_div(AD_TVplantotu,1))# q=0 1 2 可以考虑加上进化树选项 计算后转为dataframe
colnames(AD_TVplanhill) <- 'hillnumber' #hillnumber123
AD_TVplanhill$plant <- 'TVplant' #加上地址
#分别计算四个wwtp的hillnumber
AD_combinedhill <- rbind(AD_TVplanhill, AD_TGplanhill, AD_IGplanhill, AD_VNGplanhill)
ggplot(AD_combinedhill, aes(plant, hillnumber, fill =plant)) + geom_violin(trim = FALSE,scale = 'area') + geom_jitter(shape = 16, size = 2, position = position_jitter(0.2)) + stat_compare_means() + guides(fill = FALSE) + theme_classic() + scale_y_continuous(limits = c(0, max(1200))) +theme_bw()

#使用mprr进行分析 即多响应变量置换检验。该函数可以比较两个或多个群落之间的相似性，并确定它们之间的差异是否显著。
testmrpp <- mrpp(t(otu_table(AD_merged)),sample_data(AD_merged)$Plantname,permutations = 999,distance = "bray")