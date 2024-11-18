
BiocManager::install("IOBR")
library(IOBR)
library(dplyr)
library(igraph)
setwd('/home/wt/hcc_validation')

################################################################################
a=exprSet[,1]
write.table(a, file = "a.txt", quote = F, sep="\t")

expression <- right_join(exprSet, gene_name, by = c('name'= 'initial_alias'))
expression <- na.omit(expression)
expression <- expression[,2:89]
gene_list <- as.data.frame(c('PLCG1', 'FOS', 'BRCA1', 'IMPDH2', 'ATF2', 'SMAD3', 
                             'E2F4', 'TCF3', 'GRB2', 'ATF5', 'MET', 'E2F5', 'PTGS2', 
                             'MYC', 'ACTL6B', 'ELK1', 'GAB1', 'FAF1', 'EGR1', 
                             'SMARCA4', 'E2F2', 'TPMT', 'SHC3', 'NFYA', 'PIK3CA', 
                             'TP53', 'ARID1B', 'GMPS', 'MZT1', 'ARID2', 'SMARCA2', 
                             'CTNNB1', 'FOSB'))
write.csv(gene_list,'gene list.csv')
colnames(gene_list) <- "name"
used_expression <-  right_join(expression, gene_list, by = c("converted_alias" =  "name"))
used_expression <- aggregate(.~converted_alias, mean, data = used_expression)
b <- used_expression[,1]
used_expression <- used_expression[,-1]
rownames(used_expression) <- t(b)
used_expression <- t(used_expression)
write.csv(used_expression,"used_expression.csv")

expression <-  expression%>% distinct(converted_alias, .keep_all = F)
expression = expression[!duplicated(expression$converted_alias),]#谁排第一个就取谁

c_count <- cbind(expression[,88],expression[,1:13])
colnames(c_count)[1] <- "name"
ci_count <- cbind(expression[,88],expression[,14:25])
colnames(ci_count)[1] <- "name"
LDN_count <- cbind(expression[,88],expression[,26:36])
colnames(LDN_count)[1] <- "name"
HDN_count <- cbind(expression[,88],expression[,37:47])
colnames(HDN_count)[1] <- "name"
eHCC_count <- cbind(expression[,88],expression[,48:61])
colnames(eHCC_count)[1] <- "name"
aHCC_count <- cbind(expression[,88],expression[,62:73])
colnames(aHCC_count)[1] <- "name"
vaHCC_count <- cbind(expression[,88],expression[,74:87])
colnames(vaHCC_count)[1] <- "name"



HC_count <- right_join(c_count, gene_all,by = c("name"="gene"))
ci_count <- right_join(ci_count, gene_all,by = c("name"="gene"))
LDN_count <- right_join(LDN_count, gene_all,by = c("name"="gene"))
HDN_count <- right_join(HDN_count, gene_all,by = c("name"="gene"))
eHCC_count <- right_join(eHCC_count, gene_all, by = c("name"="gene"))
aHCC_count <- right_join(aHCC_count, gene_all, by = c("name"="gene"))
vaHCC_count <- right_join(vaHCC_count, gene_all, by = c("name"="gene"))

a <- HC_count[,1]
HC_count <- HC_count[,-c(1,15)]
rownames(HC_count) <- t(a)
ci_count <- ci_count[,-c(1,14)]
rownames(ci_count) <- t(a)
LDN_count <- LDN_count[,-c(1,13)]
rownames(LDN_count) <- t(a)
HDN_count <- HDN_count[,-c(1,13)]
rownames(HDN_count) <- t(a)
eHCC_count <- eHCC_count[,-c(1,16)]
rownames(eHCC_count) <- t(a)
aHCC_count <- aHCC_count[,-c(1,14)]
rownames(aHCC_count) <- t(a)
vaHCC_count <- vaHCC_count[,-c(1,16)]
rownames(vaHCC_count) <- t(a)


HC_CMI <- PCA_CMI(HC_count,lamda = 0.05,order0 = 1,G = adj_matrix)
write.csv(HC_CMI[["G"]],"HC_CMI.csv")
ci_CMI <- PCA_CMI(ci_count,lamda = 0.05,order0 = 1,G = adj_matrix)
write.csv(ci_CMI[["G"]],"ci_CMI.csv")
LDN_CMI <- PCA_CMI(LDN_count,lamda = 0.05,order0 = 1,G = adj_matrix)
write.csv(LDN_CMI[["G"]],"LDN_CMI.csv")
HDN_CMI <- PCA_CMI(HDN_count,lamda = 0.05,order0 = 1,G = adj_matrix)
write.csv(HDN_CMI[["G"]],"HDN_CMI.csv")
eHCC_CMI <- PCA_CMI(eHCC_count, lamda = 0.05, order0 = 1, G = adj_matrix)
write.csv(eHCC_CMI[["G"]],"eHCC_CMI.csv")
aHCC_CMI <- PCA_CMI(aHCC_count, lamda = 0.05, order0 = 1, G = adj_matrix)
write.csv(aHCC_CMI[["G"]],"aHCC_CMI.csv")
vaHCC_CMI <- PCA_CMI(vaHCC_count, lamda = 0.05, order0 = 1, G = adj_matrix)
write.csv(vaHCC_CMI[["G"]],"vaHCC_CMI.csv")

a <- HC_CMI[,1]
HC_CMI <- HC_CMI[,-1]
rownames(HC_CMI) <- t(a)
ci_CMI <- ci_CMI[,-1]
rownames(ci_CMI) <- t(a)
LDN_CMI <- LDN_CMI[,-1]
rownames(LDN_CMI) <- t(a)
HDN_CMI <- HDN_CMI[,-1]
rownames(HDN_CMI) <- t(a)
eHCC_CMI <- eHCC_CMI[,-1]
rownames(eHCC_CMI) <- t(a)
aHCC_CMI <- aHCC_CMI[,-1]
rownames(aHCC_CMI) <- t(a)
vaHCC_CMI <- vaHCC_CMI[,-1]
rownames(vaHCC_CMI) <- t(a)
################################################################################
num_nodes <- nrow(HC_CMI)
edge_list <- matrix(nrow = 0, ncol = 2)
for (i in 1:num_nodes) {
  for (j in (i):num_nodes) {
    if (HC_CMI[i, j] == 1) {
      edge_list <- rbind(edge_list, c(rownames(HC_CMI)[i], colnames(HC_CMI)[j]))  # 减1以匹配R中的索引从1开始
    }
  }
}
g <- graph_from_data_frame(edge_list, directed = TRUE)
g_undirected <- as.undirected(g, mode = "collapse", edge.attr.comb = "mean")
edge_list <- as_data_frame(g_undirected, what = "edges")
write.csv(edge_list,"edge_list_stage1.csv")

edge_list <- matrix(nrow = 0, ncol = 2)
for (i in 1:num_nodes) {
  for (j in (i):num_nodes) {
    if (ci_CMI[i, j] == 1) {
      edge_list <- rbind(edge_list, c(rownames(ci_CMI)[i], colnames(ci_CMI)[j]))  # 减1以匹配R中的索引从1开始
    }
  }
}
g <- graph_from_data_frame(edge_list, directed = TRUE)
g_undirected <- as.undirected(g, mode = "collapse", edge.attr.comb = "mean")
edge_list <- as_data_frame(g_undirected, what = "edges")
write.csv(edge_list,"edge_list_stage2.csv")

edge_list <- matrix(nrow = 0, ncol = 2)
for (i in 1:num_nodes) {
  for (j in (i):num_nodes) {
    if (LDN_CMI[i, j] == 1) {
      edge_list <- rbind(edge_list, c(rownames(LDN_CMI)[i], colnames(LDN_CMI)[j]))  # 减1以匹配R中的索引从1开始
    }
  }
}
g <- graph_from_data_frame(edge_list, directed = TRUE)
g_undirected <- as.undirected(g, mode = "collapse", edge.attr.comb = "mean")
edge_list <- as_data_frame(g_undirected, what = "edges")
write.csv(edge_list,"edge_list_stage3.csv")

edge_list <- matrix(nrow = 0, ncol = 2)
for (i in 1:num_nodes) {
  for (j in (i):num_nodes) {
    if (HDN_CMI[i, j] == 1) {
      edge_list <- rbind(edge_list, c(rownames(HDN_CMI)[i], colnames(HDN_CMI)[j]))  # 减1以匹配R中的索引从1开始
    }
  }
}
g <- graph_from_data_frame(edge_list, directed = TRUE)
g_undirected <- as.undirected(g, mode = "collapse", edge.attr.comb = "mean")
edge_list <- as_data_frame(g_undirected, what = "edges")
write.csv(edge_list,"edge_list_stage4.csv")



edge_list <- matrix(nrow = 0, ncol = 2)
for (i in 1:num_nodes) {
  for (j in (i):num_nodes) {
    if (eHCC_CMI[i, j] == 1) {
      edge_list <- rbind(edge_list, c(rownames(eHCC_CMI)[i], colnames(eHCC_CMI)[j]))  # 减1以匹配R中的索引从1开始
    }
  }
}
g <- graph_from_data_frame(edge_list, directed = TRUE)
g_undirected <- as.undirected(g, mode = "collapse", edge.attr.comb = "mean")
edge_list <- as_data_frame(g_undirected, what = "edges")
write.csv(edge_list,"edge_list_stage5.csv")

edge_list <- matrix(nrow = 0, ncol = 2)
for (i in 1:num_nodes) {
  for (j in (i):num_nodes) {
    if (aHCC_CMI[i, j] == 1) {
      edge_list <- rbind(edge_list, c(rownames(aHCC_CMI)[i], colnames(aHCC_CMI)[j]))  # 减1以匹配R中的索引从1开始
    }
  }
}
g <- graph_from_data_frame(edge_list, directed = TRUE)
g_undirected <- as.undirected(g, mode = "collapse", edge.attr.comb = "mean")
edge_list <- as_data_frame(g_undirected, what = "edges")
write.csv(edge_list,"edge_list_stage6.csv")

edge_list <- matrix(nrow = 0, ncol = 2)
for (i in 1:num_nodes) {
  for (j in (i):num_nodes) {
    if (vaHCC_CMI[i, j] == 1) {
      edge_list <- rbind(edge_list, c(rownames(vaHCC_CMI)[i], colnames(vaHCC_CMI)[j]))  # 减1以匹配R中的索引从1开始
    }
  }
}
g <- graph_from_data_frame(edge_list, directed = TRUE)
g_undirected <- as.undirected(g, mode = "collapse", edge.attr.comb = "mean")
edge_list <- as_data_frame(g_undirected, what = "edges")
write.csv(edge_list,"edge_list_stage7.csv")

############################### 计算pcc ########################################

HC_count <-t(HC_count)
write.csv(HC_count,'HC_count.csv')
matPCC_HC <- cor(HC_count, method='pearson')
diag(matPCC_HC)

se_net_HC <- matrix(NA, ncol = 3, nrow = nrow(matPCC_HC) * ncol(matPCC_HC))
colnames(se_net_HC) <- c("gene1", "gene2", "weight")

k <- 1
for (i in rownames(matPCC_HC)) {
  for (j in colnames(matPCC_HC)) {
    row_name <- i
    col_name <- j
    value <- matPCC_HC[i, j]
    se_net_HC[k, ] <- c(row_name, col_name, value)
    k <- k + 1
  }
}
se_net_HC <- as.data.frame(se_net_HC)


ci_count <-t(ci_count)
write.csv(ci_count,'ci_count.csv')
matPCC_ci <- cor(ci_count, method='pearson')
diag(matPCC_ci)
se_net_ci <- matrix(NA, ncol = 3, nrow = nrow(matPCC_ci) * ncol(matPCC_ci))
colnames(se_net_ci) <- c("gene1", "gene2", "weight")
k <- 1
for (i in rownames(matPCC_ci)) {
  for (j in colnames(matPCC_ci)) {
    row_name <- i
    col_name <- j
    value <- matPCC_ci[i, j]
    se_net_ci[k, ] <- c(row_name, col_name, value)
    k <- k + 1
  }
}
se_net_ci <- as.data.frame(se_net_ci)


LDN_count <-t(LDN_count)
write.csv(LDN_count,'LDN_count.csv')
matPCC_LDN <- cor(LDN_count, method='pearson')
diag(matPCC_LDN)
se_net_LDN <- matrix(NA, ncol = 3, nrow = nrow(matPCC_LDN) * ncol(matPCC_LDN))
colnames(se_net_LDN) <- c("gene1", "gene2", "weight")
k <- 1
for (i in rownames(matPCC_LDN)) {
  for (j in colnames(matPCC_LDN)) {
    row_name <- i
    col_name <- j
    value <- matPCC_LDN[i, j]
    se_net_LDN[k, ] <- c(row_name, col_name, value)
    k <- k + 1
  }
}
se_net_LDN <- as.data.frame(se_net_LDN)


HDN_count <-t(HDN_count)
write.csv(HDN_count,'HDN_count.csv')
matPCC_HDN <- cor(HDN_count, method='pearson')
diag(matPCC_HDN)
se_net_HDN <- matrix(NA, ncol = 3, nrow = nrow(matPCC_HDN) * ncol(matPCC_HDN))
colnames(se_net_HDN) <- c("gene1", "gene2", "weight")
k <- 1
for (i in rownames(matPCC_HDN)) {
  for (j in colnames(matPCC_HDN)) {
    row_name <- i
    col_name <- j
    value <- matPCC_HDN[i, j]
    se_net_HDN[k, ] <- c(row_name, col_name, value)
    k <- k + 1
  }
}
se_net_HDN <- as.data.frame(se_net_HDN)


eHCC_count <-t(eHCC_count)
write.csv(eHCC_count,'eHCC_count.csv')
matPCC_eHCC <- cor(eHCC_count, method='pearson')
diag(matPCC_eHCC)
se_net_eHCC <- matrix(NA, ncol = 3, nrow = nrow(matPCC_eHCC) * ncol(matPCC_eHCC))
colnames(se_net_eHCC) <- c("gene1", "gene2", "weight")
k <- 1
for (i in rownames(matPCC_eHCC)) {
  for (j in colnames(matPCC_eHCC)) {
    row_name <- i
    col_name <- j
    value <- matPCC_eHCC[i, j]
    se_net_eHCC[k, ] <- c(row_name, col_name, value)
    k <- k + 1
  }
}
se_net_eHCC <- as.data.frame(se_net_eHCC)


aHCC_count <-t(aHCC_count)
write.csv(aHCC_count,'aHCC_count.csv')
matPCC_aHCC <- cor(aHCC_count, method='pearson')
diag(matPCC_aHCC)
se_net_aHCC <- matrix(NA, ncol = 3, nrow = nrow(matPCC_aHCC) * ncol(matPCC_aHCC))
colnames(se_net_aHCC) <- c("gene1", "gene2", "weight")
k <- 1
for (i in rownames(matPCC_aHCC)) {
  for (j in colnames(matPCC_aHCC)) {
    row_name <- i
    col_name <- j
    value <- matPCC_aHCC[i, j]
    se_net_aHCC[k, ] <- c(row_name, col_name, value)
    k <- k + 1
  }
}
se_net_aHCC <- as.data.frame(se_net_aHCC)


vaHCC_count <-t(vaHCC_count)
write.csv(vaHCC_count,'vaHCC_count.csv')
matPCC_vaHCC <- cor(vaHCC_count, method='pearson')
diag(matPCC_vaHCC)
se_net_vaHCC <- matrix(NA, ncol = 3, nrow = nrow(matPCC_vaHCC) * ncol(matPCC_vaHCC))
colnames(se_net_vaHCC) <- c("gene1", "gene2", "weight")
k <- 1
for (i in rownames(matPCC_vaHCC)) {
  for (j in colnames(matPCC_vaHCC)) {
    row_name <- i
    col_name <- j
    value <- matPCC_vaHCC[i, j]
    se_net_vaHCC[k, ] <- c(row_name, col_name, value)
    k <- k + 1
  }
}
se_net_vaHCC <- as.data.frame(se_net_vaHCC)

################################################################################

edge_list_stage1_removed <- as.data.frame(edge_list_stage1_removed)
edge_list_stage1_removed_1 <- cbind(edge_list_stage1_removed,paste(edge_list_stage1_removed[, 1],edge_list_stage1_removed[, 2], sep = ""))
colnames(edge_list_stage1_removed_1)[3] <- "combined"

edge_list_stage2_removed <- as.data.frame(edge_list_stage2_removed)
edge_list_stage2_removed_1 <- cbind(edge_list_stage2_removed,paste(edge_list_stage2_removed[, 1],edge_list_stage2_removed[, 2], sep = ""))
colnames(edge_list_stage2_removed_1)[3] <- "combined"

edge_list_stage3_removed <- as.data.frame(edge_list_stage3_removed)
edge_list_stage3_removed_1 <- cbind(edge_list_stage3_removed,paste(edge_list_stage3_removed[, 1],edge_list_stage3_removed[, 2], sep = ""))
colnames(edge_list_stage3_removed_1)[3] <- "combined"

edge_list_stage4_removed <- as.data.frame(edge_list_stage4_removed)
edge_list_stage4_removed_1 <- cbind(edge_list_stage4_removed,paste(edge_list_stage4_removed[, 1],edge_list_stage4_removed[, 2], sep = ""))
colnames(edge_list_stage4_removed_1)[3] <- "combined"

edge_list_stage5_removed <- as.data.frame(edge_list_stage5_removed)
edge_list_stage5_removed_1 <- cbind(edge_list_stage5_removed,paste(edge_list_stage5_removed[, 1],edge_list_stage5_removed[, 2], sep = ""))
colnames(edge_list_stage5_removed_1)[3] <- "combined"

edge_list_stage6_removed <- as.data.frame(edge_list_stage6_removed)
edge_list_stage6_removed_1 <- cbind(edge_list_stage6_removed,paste(edge_list_stage6_removed[, 1],edge_list_stage6_removed[, 2], sep = ""))
colnames(edge_list_stage6_removed_1)[3] <- "combined"

edge_list_stage7_removed <- as.data.frame(edge_list_stage7_removed)
edge_list_stage7_removed_1 <- cbind(edge_list_stage7_removed,paste(edge_list_stage7_removed[, 1],edge_list_stage7_removed[, 2], sep = ""))
colnames(edge_list_stage7_removed_1)[3] <- "combined"

################################################################################

edge_list_stage1_removed <- as.data.frame(edge_list_stage1[,2:3])
edge_list_stage1_removed_1 <- cbind(edge_list_stage1_removed,paste(edge_list_stage1_removed[, 1],edge_list_stage1_removed[, 2], sep = ""))
colnames(edge_list_stage1_removed_1)[3] <- "combined"

edge_list_stage2_removed <- as.data.frame(edge_list_stage2[,2:3])
edge_list_stage2_removed_1 <- cbind(edge_list_stage2_removed,paste(edge_list_stage2_removed[, 1],edge_list_stage2_removed[, 2], sep = ""))
colnames(edge_list_stage2_removed_1)[3] <- "combined"

edge_list_stage3_removed <- as.data.frame(edge_list_stage3[,2:3])
edge_list_stage3_removed_1 <- cbind(edge_list_stage3_removed,paste(edge_list_stage3_removed[, 1],edge_list_stage3_removed[, 2], sep = ""))
colnames(edge_list_stage3_removed_1)[3] <- "combined"

edge_list_stage4_removed <- as.data.frame(edge_list_stage4[,2:3])
edge_list_stage4_removed_1 <- cbind(edge_list_stage4_removed,paste(edge_list_stage4_removed[, 1],edge_list_stage4_removed[, 2], sep = ""))
colnames(edge_list_stage4_removed_1)[3] <- "combined"

edge_list_stage5_removed <- as.data.frame(edge_list_stage5[,2:3])
edge_list_stage5_removed_1 <- cbind(edge_list_stage5_removed,paste(edge_list_stage5_removed[, 1],edge_list_stage5_removed[, 2], sep = ""))
colnames(edge_list_stage5_removed_1)[3] <- "combined"

edge_list_stage6_removed <- as.data.frame(edge_list_stage6[,2:3])
edge_list_stage6_removed_1 <- cbind(edge_list_stage6_removed,paste(edge_list_stage6_removed[, 1],edge_list_stage6_removed[, 2], sep = ""))
colnames(edge_list_stage6_removed_1)[3] <- "combined"

edge_list_stage7_removed <- as.data.frame(edge_list_stage7[,2:3])
edge_list_stage7_removed_1 <- cbind(edge_list_stage7_removed,paste(edge_list_stage7_removed[, 1],edge_list_stage7_removed[, 2], sep = ""))
colnames(edge_list_stage7_removed_1)[3] <- "combined"
################################################################################

se_net_S1_SE <- cbind(se_net_HC, paste(se_net_HC[, 1], se_net_HC[, 2], sep = ""))
colnames(se_net_S1_SE)[4] <- "combined1"
se_net_S1_SE_1 <-right_join(se_net_S1_SE,edge_list_stage1_removed_1,by=c("combined1"="combined"))
se_net_S1_SE_1 <-se_net_S1_SE_1[,c(1:3)]
colnames(se_net_S1_SE_1)<- c("gene1","gene2","weight")
se_net_S1_SE_1 <- na.omit(se_net_S1_SE_1)
write.csv(se_net_S1_SE_1,"stage1.csv")

se_net_S2_SE <- cbind(se_net_ci, paste(se_net_ci[, 1], se_net_ci[, 2], sep = ""))
colnames(se_net_S2_SE)[4] <- "combined1"
se_net_S2_SE_1 <-right_join(se_net_S2_SE,edge_list_stage2_removed_1,by=c("combined1"="combined"))
se_net_S2_SE_1 <-se_net_S2_SE_1[,c(1:3)]
colnames(se_net_S2_SE_1)<- c("gene1","gene2","weight")
se_net_S2_SE_1 <- na.omit(se_net_S2_SE_1)
write.csv(se_net_S2_SE_1,"stage2.csv")

se_net_S3_SE <- cbind(se_net_LDN, paste(se_net_LDN[, 1], se_net_LDN[, 2], sep = ""))
colnames(se_net_S3_SE)[4] <- "combined1"
se_net_S3_SE_1 <-right_join(se_net_S3_SE,edge_list_stage3_removed_1,by=c("combined1"="combined"))
se_net_S3_SE_1 <-se_net_S3_SE_1[,c(1:3)]
colnames(se_net_S3_SE_1)<- c("gene1","gene2","weight")
se_net_S3_SE_1 <- na.omit(se_net_S3_SE_1)
write.csv(se_net_S3_SE_1,"stage3.csv")

se_net_S4_SE <- cbind(se_net_HDN, paste(se_net_HDN[, 1], se_net_HDN[, 2], sep = ""))
colnames(se_net_S4_SE)[4] <- "combined1"
se_net_S4_SE_1 <-right_join(se_net_S4_SE,edge_list_stage4_removed_1,by=c("combined1"="combined"))
se_net_S4_SE_1 <-se_net_S4_SE_1[,c(1:3)]
colnames(se_net_S4_SE_1)<- c("gene1","gene2","weight")
se_net_S4_SE_1 <- na.omit(se_net_S4_SE_1)
write.csv(se_net_S4_SE_1,"stage4.csv")

se_net_S6_SE <- cbind(se_net_eHCC, paste(se_net_eHCC[, 1], se_net_eHCC[, 2], sep = ""))
colnames(se_net_S6_SE)[4] <- "combined1"
se_net_S6_SE_1 <-right_join(se_net_S6_SE,edge_list_stage5_removed_1,by=c("combined1"="combined"))
se_net_S6_SE_1 <-se_net_S6_SE_1[,c(1:3)]
colnames(se_net_S6_SE_1)<- c("gene1","gene2","weight")
se_net_S6_SE_1 <- na.omit(se_net_S6_SE_1)
write.csv(se_net_S6_SE_1,"stage5.csv")

se_net_S7_SE <- cbind(se_net_aHCC, paste(se_net_aHCC[, 1], se_net_aHCC[, 2], sep = ""))
colnames(se_net_S7_SE)[4] <- "combined1"
se_net_S7_SE_1 <-right_join(se_net_S7_SE,edge_list_stage6_removed_1,by=c("combined1"="combined"))
se_net_S7_SE_1 <-se_net_S7_SE_1[,c(1:3)]
colnames(se_net_S7_SE_1)<- c("gene1","gene2","weight")
se_net_S7_SE_1 <- na.omit(se_net_S7_SE_1)
write.csv(se_net_S7_SE_1,"stage6.csv")

se_net_S8_SE <- cbind(se_net_vaHCC, paste(se_net_vaHCC[, 1], se_net_vaHCC[, 2], sep = ""))
colnames(se_net_S8_SE)[4] <- "combined1"
se_net_S8_SE_1 <-right_join(se_net_S8_SE,edge_list_stage7_removed_1,by=c("combined1"="combined"))
se_net_S8_SE_1 <-se_net_S8_SE_1[,c(1:3)]
colnames(se_net_S8_SE_1)<- c("gene1","gene2","weight")
se_net_S8_SE_1 <- na.omit(se_net_S8_SE_1)
write.csv(se_net_S8_SE_1,"stage7.csv")

stage_all <- rbind(stage1, stage2, stage3, stage4, stage5, stage6, stage7)
stage_all <- stage_all[2:3]
g <- graph_from_data_frame(stage_all, directed = TRUE)
g_undirected <- as.undirected(g, mode = "collapse", edge.attr.comb = "mean")
stage_all <- as_data_frame(g_undirected, what = "edges")
write.csv(stage_all, 'stage_all.csv')

################################################################################