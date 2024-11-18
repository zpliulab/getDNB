library(dplyr)
library(DT)
library(SummarizedExperiment)
library(TCGAbiolinks)
BiocManager::install("TCGAbiolinks")
library(stringi)
library(clusterProfiler, lib.loc = '/home/dell/R/x86_64-conda-linux-gnu-library/4.2')
library(org.Hs.eg.db)
library(igraph, lib.loc = '/home/dell/R/x86_64-conda-linux-gnu-library/4.2')
setwd('/home/wt/HCC_TCGA')

################################################################################

TCGAbiolinks:::getProjectSummary("TCGA-LIHC")
cancer_type="TCGA-LIHC"
clinical <- GDCquery_clinic(project= cancer_type,type = "clinical")   #1098
save(clinical,file="LIHC_clinical.Rdata")
write.csv(clinical, file="TCGAbiolinks-LIHC
          -clinical.csv")

################################################################################

data_type <- "Gene Expression Quantification"
data_category <- "Transcriptome Profiling"
workflow_type <- "STAR - Counts"
query_TranscriptomeCounts <- GDCquery(project = cancer_type,
                                      data.category = data_category,
                                      data.type =  data_type, 
                                      workflow.type = workflow_type )
GDCdownload(query_TranscriptomeCounts, method = "api")

expdat <- GDCprepare(query = query_TranscriptomeCounts)
count_matrix=assay(expdat)
write.csv(count_matrix,file = "TCGAbiolinks_LIHC_counts.csv")       #60660  425

################################################################################


TCGAbiolinks_LIHC_counts <- read_csv("TCGAbiolinks_LIHC_counts.csv")
TCGAbiolinks_LIHC_counts$...1=stri_sub(TCGAbiolinks_LIHC_counts$...1,1,15)

Ensembl_ID <- TCGAbiolinks_LIHC_counts$...1
gene_symbol <- bitr(Ensembl_ID, fromType="ENSEMBL", toType=c("SYMBOL", "ENTREZID"), OrgDb="org.Hs.eg.db")
write.csv(gene_symbol,file = "gene_symbol.csv")     # 36571 3

count <- right_join(gene_symbol,TCGAbiolinks_LIHC_counts,by=c("ENSEMBL"="...1"))
count <- na.omit(count)
count <-count %>% distinct(SYMBOL, .keep_all = TRUE)
row_names <- count[,2]
data_part <- count[,-c(1:4)]
LIHC_counts <- as.data.frame(data_part)
rownames(LIHC_counts) <- t(row_names)
colnames(LIHC_counts)=stri_sub(colnames(LIHC_counts),1,15)
LIHC_counts <- log(LIHC_counts + 1)

################################################################################

clinical_4 <- TCGAbiolinks_LIHC_clinical %>% filter(grepl("Stage IV", ajcc_pathologic_stage, ignore.case = TRUE))
clinical_3 <- TCGAbiolinks_LIHC_clinical %>% filter(grepl("Stage III", ajcc_pathologic_stage, ignore.case = TRUE))
new_matrix <- TCGAbiolinks_LIHC_clinical %>% filter(!grepl("Stage III", ajcc_pathologic_stage, ignore.case = TRUE))
new_matrix <- new_matrix %>% filter(!grepl("Stage IV", ajcc_pathologic_stage, ignore.case = TRUE))
clinical_2 <- new_matrix %>% filter(grepl("Stage II", ajcc_pathologic_stage, ignore.case = TRUE))
new_matrix <- new_matrix %>% filter(!grepl("Stage II", ajcc_pathologic_stage, ignore.case = TRUE))
clinical_1 <- new_matrix %>% filter(grepl("Stage I", ajcc_pathologic_stage, ignore.case = TRUE))
new_matrix <- new_matrix %>% filter(!grepl("Stage I", ajcc_pathologic_stage, ignore.case = TRUE))
rm(new_matrix)

################################################################################

tumor_counts <- LIHC_counts[, grepl("-01",names(LIHC_counts))]
normal_counts <- LIHC_counts[, grepl("-11",names(LIHC_counts))]
tumor_counts <- cbind(colnames(tumor_counts), t(tumor_counts))
tumor_counts <- as.data.frame(tumor_counts)
tumor_counts$V1 <- stri_sub(tumor_counts$V1,1,12)
normal_counts <- cbind(colnames(normal_counts), t(normal_counts))
normal_counts <- as.data.frame(normal_counts)
normal_counts$V1 <- stri_sub(normal_counts$V1,1,12)
count0 <- normal_counts
count1 <- right_join(tumor_counts,clinical_1,by = c('V1' = 'submitter_id'))
count2 <- right_join(tumor_counts,clinical_2,by = c('V1' = 'submitter_id'))
count3 <- right_join(tumor_counts,clinical_3,by = c('V1' = 'submitter_id'))
count4 <- right_join(tumor_counts,clinical_4,by = c('V1' = 'submitter_id'))

count1 <- count1[,1:36423]
count2 <- count2[,1:36423]
count3 <- count3[,1:36423]
count4 <- count4[,1:36423]
count1 <- na.omit(count1)
count2 <- na.omit(count2)
count3 <- na.omit(count3)
count4 <- na.omit(count4)

a <- count0[,1]
count0 <- count0[,-1]
rownames(count0) <- t(a)
a <- count1[,1]
count1 <- count1[,-1]
rownames(count1) <- t(a)
a <- count2[,1]
count2 <- count2[,-1]
rownames(count2) <- t(a)
a <- count3[,1]
count3 <- count3[,-1]
rownames(count3) <- t(a)
a <- count4[,1]
count4 <- count4[,-1]
rownames(count4) <- t(a)
count0 <- t(count0)
count1 <- t(count1)
count2 <- t(count2)
count3 <- t(count3)
count4 <- t(count4)
write.csv(count0,'count0.csv')
write.csv(count1,'count1.csv')
write.csv(count2,'count2.csv')
write.csv(count3,'count3.csv')
write.csv(count4,'count4.csv')
count0 <- right_join(count0, gene_all,by = c("...1"="gene"))
count1 <- right_join(count1, gene_all,by = c("...1"="gene"))
count2 <- right_join(count2, gene_all,by = c("...1"="gene"))
count3 <- right_join(count3, gene_all,by = c("...1"="gene"))
count4 <- right_join(count4, gene_all,by = c("...1"="gene"))
new_colnames <- paste0(colnames(count0),"-00")
colnames(count0) <- new_colnames
new_colnames <- paste0(colnames(count1),"-01")
colnames(count1) <- new_colnames
new_colnames <- paste0(colnames(count2),"-02")
colnames(count2) <- new_colnames
new_colnames <- paste0(colnames(count3),"-03")
colnames(count3) <- new_colnames
new_colnames <- paste0(colnames(count4),"-04")
colnames(count4) <- new_colnames
count0 <- as.data.frame(count0)
count1 <- as.data.frame(count1)
count2 <- as.data.frame(count2)
count3 <- as.data.frame(count3)
count4 <- as.data.frame(count4)
count <- cbind(count0,count1,count2,count3,count4)
count <- as.data.frame(t(count))
a <- count[,1]
count <- count[,-1]
rownames(count) <- t(a)
gene_HCC <- as.data.frame(c('PLCG1', 'FOS', 'BRCA1', 'IMPDH2', 'ATF2', 'SMAD3', 
                            'E2F4', 'TCF3', 'GRB2', 'ATF5', 'MET', 'E2F5', 'PTGS2', 
                            'MYC', 'ACTL6B', 'ELK1', 'GAB1', 'FAF1', 'EGR1', 
                            'SMARCA4', 'E2F2', 'TPMT', 'SHC3', 'NFYA', 'PIK3CA', 
                            'TP53', 'ARID1B', 'GMPS', 'MZT1', 'ARID2', 'SMARCA2', 
                            'CTNNB1', 'FOSB'))
colnames(gene_HCC) <- "name"
HCC_splsda <- as.data.frame(c( "CLEC4G", "FCN2", "OIT3", "CLEC1B", "CAP2","COLEC10",
                               "FCN3", "CRHBP", "CFP", "ANGPTL6", "STAB2" ,  "LYVE1",
                               "CLEC4M", "ECM1", "MARCO", "ASPM", "DNASE1L3","DIRAS3", 
                               "RSPO3", "GMPS", "NAT2", "CYP1A2", "TUBE1","DNAJC6",
                               "TARBP1", "CLTRN", "TIMD4", "NVL", "ITLN1", "PABPC1",
                               "PPOX","MFSD2A","HMMR"))
colnames(HCC_splsda) <- "name"
HCC_DUBStepR <- as.data.frame(c("KIF4A","CDKN3","DNASE1L3","CDC20","RACGAP1","CRHBP","FCN3","RRM2","ANLN","MCM3","MAD2L1",  
                                "UBE2C","CCNA2","KIF20A","CENPF","PBK","BUB1B","NUF2","TOP2A","ASPM","CCNB1","SGO2",
                                "STMN1","KDM5D","MKI67","PTTG1","DDX3Y","PARPBP","TRIP13","TTK","DLGAP5","CLEC1B","GINS1"))
colnames(HCC_DUBStepR) <- "name"
HCC_SVMRFE <- as.data.frame(c('ABCC9', 'AFP', 'APCDD1', 'ARMCX3', 'CMAHP', 'CYP1B1', 'DLK1', 'DOK5',
                              'FAM50B', 'FOXP2', 'FRMD6', 'HAMP', 'HOOK1', 'IFI44L', 'IL1R2', 'IRS2',
                              'LAMB1', 'LOC105379173', 'MFAP3L', 'NAP1L3', 'NRCAM', 'PDCD6P1', 'PKIB',
                              'PLA2G2A', 'PMAIP1', 'PPARGC1A', 'PPP1R3B', 'PRMT6', 'PTPRD', 'PTX3',
                              'RGS4', 'S100A9', 'SAMD5'))
colnames(HCC_SVMRFE) <- "name"
HCC_RFRFE <- as.data.frame(c('SEMA3F','SEC62','RBM7','SEH1L','RNMT','SARS2',
                             'RNASEH2A','SERINC1','RBP2','RND3','SCAMP3',
                             'RGS4','SEM1','RAN','SEC61G','SEMA6D','SCYL1',
                             'RASA1','SANBR','RFT1','SAP30','RNF180','SAC3D1',
                             'SCARA3','RNF187','SAMD12','RIC8A','SERHL2','RELN',
                             'SDHD','RBMXL1','SEC22B'))
colnames(HCC_RFRFE) <- "name"
HCC_SERUAT <- as.data.frame(c('DKK1','PEG10','AFP','LINC01419','EPCAM','GPC3','REG3A','NQO1',
                              'S100P','MMP12','KCNN2','IGF2BP3','TOP2A','MYH4','ASPM',
                              'ALDH3A1','COX7B2','ANLN','LINC01093','CCNB1','ODAM','NEK2',
                              'MEP1A','TTK','CLEC1B','CNDP1','OIT3','FDCSP','NUF2','LIN28B','LGR5','POSTN','CLEC4G'))
colnames(HCC_SERUAT) <- "name"
HCC_V7  <- as.data.frame(c('KRT23','CTDSP2',	'STC1',	'GTPBP6',	'HSPA4L',	'DENND11',
                           'CFHR5',	'SENP2',	'LINC01806',	'JTB',	'SLC39A11',	
                           'LRRC1',	'EIF3L',	'PRR14L',	'PTPRD',	'SEC61A2',	
                           'NR2C2',	'SLC22A15',	'NEURL1B',	'TRIM31',	'CREG1',
                           'FAM168B',	'SEL1L3',	'LOXL4',	'MAGEH1',	'OAS2',
                           'COL4A3',	'RPL3',	'NR2F1',	'DLEU1'))
colnames(HCC_V7) <- 'name'
colnames(HCC_V8) <- 'name'




count <- log(count+1)
count <- cbind(rownames(count),count)
colnames(count)[1] <- 'gene'
count_1 <- right_join(count, HCC_V8, by = c('gene'= 'name'))
count_1 <- na.omit(count_1)
count_1 <- t(count_1)
write.csv(count_1,'TCGA_count_NIEE.csv')

count0 <- count0[,1:50]
count1 <- count1[,1:171]
count2 <- count2[,1:85]
count3 <- count3[,1:85]
count4 <- count4[,1:5]
################################################################################

stage0_CMI <- PCA_CMI(count0,lamda = 0.05,order0 = 1,G = adj_matrix)
write.csv(stage0_CMI[["G"]],"stage0_CMI.csv")
stage1_CMI <- PCA_CMI(count1,lamda = 0.05,order0 = 1,G = adj_matrix)
write.csv(stage1_CMI[["G"]],"stage1_CMI.csv")
stage2_CMI <- PCA_CMI(count2,lamda = 0.05,order0 = 1,G = adj_matrix)
write.csv(stage2_CMI[["G"]],"stage2_CMI.csv")
stage3_CMI <- PCA_CMI(count3,lamda = 0.05,order0 = 1,G = adj_matrix)
write.csv(stage3_CMI[["G"]],"stage3_CMI.csv")
stage4_CMI <- PCA_CMI(count4,lamda = 0.005,order0 = 1,G = adj_matrix)
write.csv(stage4_CMI[["G"]],"stage4_CMI.csv")


########################  读取CMI  #############################################

a <- stage0_CMI[,1]
stage0_CMI <- stage0_CMI[,-1]
rownames(stage0_CMI) <- t(a)
stage1_CMI <- stage1_CMI[,-1]
rownames(stage1_CMI) <- t(a)
stage2_CMI <- stage2_CMI[,-1]
rownames(stage2_CMI) <- t(a)
stage3_CMI <- stage3_CMI[,-1]
rownames(stage3_CMI) <- t(a)
stage4_CMI <- stage4_CMI[,-1]
rownames(stage4_CMI) <- t(a)

################################################################################

num_nodes <- nrow(stage0_CMI)
edge_list <- matrix(nrow = 0, ncol = 2)
for (i in 1:num_nodes) {
  for (j in (i):num_nodes) {
    if (stage0_CMI[i, j] == 1) {
      edge_list <- rbind(edge_list, c(rownames(stage0_CMI)[i], colnames(stage0_CMI)[j]))  # 减1以匹配R中的索引从1开始
    }
  }
}
g <- graph_from_data_frame(edge_list, directed = TRUE)
g_undirected <- as.undirected(g, mode = "collapse", edge.attr.comb = "mean")
g_undirected <- simplify(g_undirected, remove.multiple = FALSE, remove.loops = TRUE)
edge_list <- as_data_frame(g_undirected, what = "edges")
write.csv(edge_list,"edge_list_stage0.csv")

edge_list <- matrix(nrow = 0, ncol = 2)
for (i in 1:num_nodes) {
  for (j in (i):num_nodes) {
    if (stage1_CMI[i, j] == 1) {
      edge_list <- rbind(edge_list, c(rownames(stage1_CMI)[i], colnames(stage1_CMI)[j]))  # 减1以匹配R中的索引从1开始
    }
  }
}
g <- graph_from_data_frame(edge_list, directed = TRUE)
g_undirected <- as.undirected(g, mode = "collapse", edge.attr.comb = "mean")
g_undirected <- simplify(g_undirected, remove.multiple = FALSE, remove.loops = TRUE)
edge_list <- as_data_frame(g_undirected, what = "edges")
write.csv(edge_list,"edge_list_stage1.csv")

edge_list <- matrix(nrow = 0, ncol = 2)
for (i in 1:num_nodes) {
  for (j in (i):num_nodes) {
    if (stage2_CMI[i, j] == 1) {
      edge_list <- rbind(edge_list, c(rownames(stage2_CMI)[i], colnames(stage2_CMI)[j]))  # 减1以匹配R中的索引从1开始
    }
  }
}
g <- graph_from_data_frame(edge_list, directed = TRUE)
g_undirected <- as.undirected(g, mode = "collapse", edge.attr.comb = "mean")
g_undirected <- simplify(g_undirected, remove.multiple = FALSE, remove.loops = TRUE)
edge_list <- as_data_frame(g_undirected, what = "edges")
write.csv(edge_list,"edge_list_stage2.csv")

edge_list <- matrix(nrow = 0, ncol = 2)
for (i in 1:num_nodes) {
  for (j in (i):num_nodes) {
    if (stage3_CMI[i, j] == 1) {
      edge_list <- rbind(edge_list, c(rownames(stage3_CMI)[i], colnames(stage3_CMI)[j]))  # 减1以匹配R中的索引从1开始
    }
  }
}
g <- graph_from_data_frame(edge_list, directed = TRUE)
g_undirected <- as.undirected(g, mode = "collapse", edge.attr.comb = "mean")
g_undirected <- simplify(g_undirected, remove.multiple = FALSE, remove.loops = TRUE)
edge_list <- as_data_frame(g_undirected, what = "edges")
write.csv(edge_list,"edge_list_stage3.csv")

edge_list <- matrix(nrow = 0, ncol = 2)
for (i in 1:num_nodes) {
  for (j in (i):num_nodes) {
    if (stage4_CMI[i, j] == 1) {
      edge_list <- rbind(edge_list, c(rownames(stage4_CMI)[i], colnames(stage4_CMI)[j]))  # 减1以匹配R中的索引从1开始
    }
  }
}
g <- graph_from_data_frame(edge_list, directed = TRUE)
g_undirected <- as.undirected(g, mode = "collapse", edge.attr.comb = "mean")
g_undirected <- simplify(g_undirected, remove.multiple = FALSE, remove.loops = TRUE)
edge_list <- as_data_frame(g_undirected, what = "edges")
write.csv(edge_list,"edge_list_stage4.csv")

############################### 计算pcc ########################################

count0 <-t(count0)
matPCC_0 <- cor(count0, method='pearson')
diag(matPCC_0)
se_net_0 <- matrix(NA, ncol = 3, nrow = nrow(matPCC_0) * ncol(matPCC_0))
colnames(se_net_0) <- c("gene1", "gene2", "weight")

k <- 1
for (i in rownames(matPCC_0)) {
  for (j in colnames(matPCC_0)) {
    row_name <- i
    col_name <- j
    value <- matPCC_0[i, j]
    se_net_0[k, ] <- c(row_name, col_name, value)
    k <- k + 1
  }
}
se_net_0 <- as.data.frame(se_net_0)


count1 <-t(count1)
matPCC_1 <- cor(count1, method='pearson')
diag(matPCC_1)
se_net_1 <- matrix(NA, ncol = 3, nrow = nrow(matPCC_1) * ncol(matPCC_1))
colnames(se_net_1) <- c("gene1", "gene2", "weight")
k <- 1
for (i in rownames(matPCC_1)) {
  for (j in colnames(matPCC_1)) {
    row_name <- i
    col_name <- j
    value <- matPCC_1[i, j]
    se_net_1[k, ] <- c(row_name, col_name, value)
    k <- k + 1
  }
}
se_net_1 <- as.data.frame(se_net_1)


count2 <-t(count2)
matPCC_2 <- cor(count2, method='pearson')
diag(matPCC_2)
se_net_2 <- matrix(NA, ncol = 3, nrow = nrow(matPCC_2) * ncol(matPCC_2))
colnames(se_net_2) <- c("gene1", "gene2", "weight")
k <- 1
for (i in rownames(matPCC_2)) {
  for (j in colnames(matPCC_2)) {
    row_name <- i
    col_name <- j
    value <- matPCC_2[i, j]
    se_net_2[k, ] <- c(row_name, col_name, value)
    k <- k + 1
  }
}
se_net_2 <- as.data.frame(se_net_2)

count3 <-t(count3)
matPCC_3 <- cor(count3, method='pearson')
diag(matPCC_3)
se_net_3 <- matrix(NA, ncol = 3, nrow = nrow(matPCC_3) * ncol(matPCC_3))
colnames(se_net_3) <- c("gene1", "gene2", "weight")
k <- 1
for (i in rownames(matPCC_3)) {
  for (j in colnames(matPCC_3)) {
    row_name <- i
    col_name <- j
    value <- matPCC_3[i, j]
    se_net_3[k, ] <- c(row_name, col_name, value)
    k <- k + 1
  }
}
se_net_3 <- as.data.frame(se_net_3)

count4 <-t(count4)
matPCC_4 <- cor(count4, method='pearson')
diag(matPCC_4)
se_net_4 <- matrix(NA, ncol = 3, nrow = nrow(matPCC_4) * ncol(matPCC_4))
colnames(se_net_4) <- c("gene1", "gene2", "weight")
k <- 1
for (i in rownames(matPCC_4)) {
  for (j in colnames(matPCC_4)) {
    row_name <- i
    col_name <- j
    value <- matPCC_4[i, j]
    se_net_4[k, ] <- c(row_name, col_name, value)
    k <- k + 1
  }
}
se_net_4 <- as.data.frame(se_net_4)

################################################################################

edge_list_stage1_removed <- as.data.frame(edge_list_stage0)
edge_list_stage1_removed_1 <- cbind(edge_list_stage1_removed,paste(edge_list_stage1_removed[, 2],edge_list_stage1_removed[, 3], sep = ""))
colnames(edge_list_stage1_removed_1)[4] <- "combined"

edge_list_stage2_removed <- as.data.frame(edge_list_stage1)
edge_list_stage2_removed_1 <- cbind(edge_list_stage2_removed,paste(edge_list_stage2_removed[, 2],edge_list_stage2_removed[, 3], sep = ""))
colnames(edge_list_stage2_removed_1)[4] <- "combined"

edge_list_stage3_removed <- as.data.frame(edge_list_stage2)
edge_list_stage3_removed_1 <- cbind(edge_list_stage3_removed,paste(edge_list_stage3_removed[, 2],edge_list_stage3_removed[, 3], sep = ""))
colnames(edge_list_stage3_removed_1)[4] <- "combined"

edge_list_stage4_removed <- as.data.frame(edge_list_stage3)
edge_list_stage4_removed_1 <- cbind(edge_list_stage4_removed,paste(edge_list_stage4_removed[, 2],edge_list_stage4_removed[, 3], sep = ""))
colnames(edge_list_stage4_removed_1)[4] <- "combined"

edge_list_stage5_removed <- as.data.frame(edge_list_stage4)
edge_list_stage5_removed_1 <- cbind(edge_list_stage5_removed,paste(edge_list_stage5_removed[, 2],edge_list_stage5_removed[, 3], sep = ""))
colnames(edge_list_stage5_removed_1)[4] <- "combined"

################################################################################
se_net_S0_SE <- cbind(se_net_0, paste(se_net_0[, 1], se_net_0[, 2], sep = ""))
colnames(se_net_S0_SE)[4] <- "combined1"
se_net_S0_SE_1 <-right_join(se_net_S0_SE,edge_list_stage1_removed_1,by=c("combined1"="combined"))
se_net_S0_SE_1 <-se_net_S0_SE_1[,c(1:3)]
colnames(se_net_S0_SE_1)<- c("gene1","gene2","weight")
se_net_S0_SE_1 <- na.omit(se_net_S0_SE_1)

write.csv(se_net_S0_SE_1,"stage0.csv")


se_net_S1_SE <- cbind(se_net_1, paste(se_net_1[, 1], se_net_1[, 2], sep = ""))
colnames(se_net_S1_SE)[4] <- "combined1"
se_net_S1_SE_1 <-right_join(se_net_S1_SE,edge_list_stage2_removed_1,by=c("combined1"="combined"))
se_net_S1_SE_1 <-se_net_S1_SE_1[,c(1:3)]
colnames(se_net_S1_SE_1)<- c("gene1","gene2","weight")
se_net_S1_SE_1 <- na.omit(se_net_S1_SE_1)
write.csv(se_net_S1_SE_1,"stage1.csv")

se_net_S2_SE <- cbind(se_net_2, paste(se_net_2[, 1], se_net_2[, 2], sep = ""))
colnames(se_net_S2_SE)[4] <- "combined1"
se_net_S2_SE_1 <-right_join(se_net_S2_SE,edge_list_stage3_removed_1,by=c("combined1"="combined"))
se_net_S2_SE_1 <-se_net_S2_SE_1[,c(1:3)]
colnames(se_net_S2_SE_1)<- c("gene1","gene2","weight")
se_net_S2_SE_1 <- na.omit(se_net_S2_SE_1)
write.csv(se_net_S2_SE_1,"stage2.csv")

se_net_S3_SE <- cbind(se_net_3, paste(se_net_3[, 1], se_net_3[, 2], sep = ""))
colnames(se_net_S3_SE)[4] <- "combined1"
se_net_S3_SE_1 <-right_join(se_net_S3_SE,edge_list_stage4_removed_1,by=c("combined1"="combined"))
se_net_S3_SE_1 <-se_net_S3_SE_1[,c(1:3)]
colnames(se_net_S3_SE_1)<- c("gene1","gene2","weight")
se_net_S3_SE_1 <- na.omit(se_net_S3_SE_1)
write.csv(se_net_S3_SE_1,"stage3.csv")

se_net_S4_SE <- cbind(se_net_4, paste(se_net_4[, 1], se_net_4[, 2], sep = ""))
colnames(se_net_S4_SE)[4] <- "combined1"
se_net_S4_SE_1 <-right_join(se_net_S4_SE,edge_list_stage5_removed_1,by=c("combined1"="combined"))
se_net_S4_SE_1 <-se_net_S4_SE_1[,c(1:3)]
colnames(se_net_S4_SE_1)<- c("gene1","gene2","weight")
se_net_S4_SE_1 <- na.omit(se_net_S4_SE_1)
write.csv(se_net_S4_SE_1,"stage4.csv")

################################################################################

