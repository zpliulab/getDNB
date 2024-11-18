library(pheatmap)
library(ggsci)
library(tidyverse)
library(cowplot)
library(enrichplot)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
setwd('/home/wt/hcc_add_kegg')
################################################################################

gene_HCC <- as.data.frame(c('PLCG1', 'FOS', 'BRCA1', 'IMPDH2', 'ATF2', 'SMAD3', 
                            'E2F4', 'TCF3', 'GRB2', 'ATF5', 'MET', 'E2F5', 'PTGS2', 
                            'MYC', 'ACTL6B', 'ELK1', 'GAB1', 'FAF1', 'EGR1', 
                            'SMARCA4', 'E2F2', 'TPMT', 'SHC3', 'NFYA', 'PIK3CA', 
                            'TP53', 'ARID1B', 'GMPS', 'MZT1', 'ARID2', 'SMARCA2', 
                            'CTNNB1', 'FOSB'))
colnames(gene_HCC) <- "gene"


ego_ALL <- enrichGO(gene          = t(gene_HCC),
                    OrgDb         = org.Hs.eg.db,
                    keyType       = 'SYMBOL',
                    ont           = "ALL",  #设置为ALL时BP, CC, MF都计算
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.1)
ego_all <- data.frame(ego_ALL)
write.csv(ego_ALL,'enrichGO_all.csv')


ego_CC <- enrichGO(gene          = t(gene_HCC),
                   OrgDb         = org.Hs.eg.db,
                   keyType       = 'SYMBOL',
                   ont           = "CC",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.1)
ego_cc <- data.frame(ego_CC)
ego_cc <- ego_cc[1:10,]
ego_cc$ONTOLOGY <- "CC"
write.csv(ego_cc,'enrichGO_cc.csv') 

ego_MF <- enrichGO(gene          = t(gene_HCC),
                   OrgDb         = org.Hs.eg.db,
                   keyType       = 'SYMBOL',
                   ont           = "MF",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05)
ego_mf <- data.frame(ego_MF)
ego_mf <- ego_mf[1:10,]
ego_mf$ONTOLOGY <- "MF"
write.csv(ego_mf,'enrichGO_mf.csv') 

ego_BP <- enrichGO(gene          = t(gene_HCC),
                   OrgDb         = org.Hs.eg.db,
                   keyType       = 'SYMBOL',
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.01) 
ego_bp <- data.frame(ego_BP)
ego_bp <- ego_bp[1:10,]
ego_bp$ONTOLOGY <- "BP"
write.csv(ego_bp,'enrichGO_bp.csv') 

enrichGO_all <- rbind(ego_bp,ego_cc,ego_mf)
rownames(enrichGO_all) <- c(1:30)
GO_term_order=factor(as.integer(rownames(enrichGO_all)),labels=enrichGO_all$Description)
for(i in 1:nrow(enrichGO_all)){
  description_splite=strsplit(enrichGO_all$Description[i],split = " ")
  description_collapse=paste(description_splite[[1]][1:5],collapse = " ") #这里的5就是指5个单词的意思，可以自己更改
  enrichGO_all$Description[i]=description_collapse
  enrichGO_all$Description=gsub(pattern = "NA","",enrichGO_all$Description)
}
COLS <- c("#66C3A5", "#8DA1CB", "#FD8D62")
P1 <- ggplot(data=enrichGO_all, aes(x=GO_term_order,y=-log10(pvalue), z = Count,fill=ONTOLOGY)) +  
  geom_bar(stat="identity", width=0.8)  +   
  scale_fill_manual(values = COLS) + 
  theme_bw()  +  
  xlab("Description") + 
  ylab("-log10(P-value)") + 
  labs(title = "The Most Enriched GO Terms")+ 
  theme(axis.text.x=element_text(face = "bold", color="gray50",angle = 70,vjust = 1, hjust = 1,size=14 ))+
  theme(axis.text.y=element_text(face = "bold", color="gray50",vjust = 1, hjust = 1,size=14 ))+
  theme(axis.title.x =element_text(size=25), axis.title.y=element_text(size=25))+
  theme(legend.text=element_text(size=25))+ 
  geom_text(aes(label = Count), vjust = -0.5, color="gray30")

ggsave("GO.pdf", plot =P1, width = 15, height = 12)

d <- t(as.data.frame(gene_HCC[,1]))
gene_HCC <- bitr(d, fromType = 'SYMBOL', toType = 'ENTREZID',OrgDb = 'org.Hs.eg.db' )
ekegg <- enrichKEGG(gene = gene_HCC[,2],organism = 'hsa', pvalueCutoff = 0.5)

# 生成渐变颜色向量
hh <- as.data.frame(ego_bp)
hh <- hh[1:10,]#
rownames(hh) <- 1:nrow(hh)
hh$order=factor(rev(as.integer(rownames(hh))),labels = rev(hh$Description))
p1 <- ggplot(hh, aes(y = order, x = Count, fill = p.adjust)) +
  geom_bar(stat = "identity", width = 0.8) + ####柱子宽度
  #scale_fill_gradient(low = "#9FBA95", high = "#F4EEAC") + 
  scale_fill_gradient(low = "#ff9999", high = "#CCE5FF") + #颜色自己可以换
  labs(title = "GO Enrichment",
       x = "Counts",
       y = " Descriptions") +
  theme(title = element_text(size = 20, color = "black"),
        axis.title.x = element_text(size = 20, color = "black"),  # 设置x轴标题大小和颜色
        axis.title.y = element_text(size = 20, color = "black"),  # 设置y轴标题大小和颜色
        axis.text.x = element_text(size = 15, color = "black"),   # 设置x轴刻度文本大小和颜色
        axis.text.y = element_text(size = 15, color = "black"),   # 设置y轴刻度文本大小和颜色
        legend.text = element_text(size = 15, color = "black"),   # 设置图例文本大小和颜色
        legend.title = element_text(size = 15, color = "black"),  # 设置图例标题大小和颜色
        panel.background = element_blank(),      # 去掉背景颜色
        plot.background = element_blank())       # 去掉整个绘图区域的背景颜色
p1
plotc = p1
ggsave("enrichGO.pdf", plot = plotc, width = 14, height = 5)




