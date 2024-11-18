setwd('/home/wt/hcc_add_kegg')
library(dplyr)
library(stringi)
library(DT)
library(tidyverse)

##########################  获取network_used和gene all #########################

Regnetwork <- read.table("human_network.txt",header = T)
colnames(anova_genes)[3] <- "pvalue"
colnames(anova_genes)[1] <- "X2"
total_gene <- subset(anova_genes,pvalue < (0.00001))
total_gene <- rbind(total_gene[,1],KEGG[,2])
total_gene <- unique(total_gene)
#
#net_used <- right_join(Regnetwork,count_1,by=c("V1"="gene"))
#net_used <- right_join(net_used,count_1,by=c("V3"="gene"))
#net_used <- na.omit(net_used)
#write.csv(net_used[,1:4],"network_used.csv")


net_used <- right_join(Regnetwork,total_gene,by=c("V1"="X2"))
net_used <- right_join(net_used,total_gene,by=c("V3"="X2"))
net_used <- na.omit(net_used)
write.csv(net_used,"network_used.csv")
gene1 <- as.data.frame(unique(net_used[,c(1,2)]))
gene2 <- as.data.frame(unique(net_used[,c(3,4)]))
colnames(gene1) <- "gene"
colnames(gene2) <- "gene"
gene_all <-unique(rbind(gene1,gene2))
colnames(gene_all)[2] <- "gene ID"
write.csv(gene_all,"gene_used.csv")
gene_prior <- right_join(gene_all,KEGG,by = c('gene'='X2'))
########################  获取adj_matrix  ######################################
network_used <- net_used[,c(1,3)]
colnames(network_used)[1] <- "gene1"
colnames(network_used)[2] <- "gene2"

filtered_rows <- matrix(NA, ncol = 2)
colnames(filtered_rows)[1] <- "gene1"
colnames(filtered_rows)[2] <- "gene2"
for(i in 1:7534)
{
  matching_row <- which(network_used[i,2] == gene_all[,1])
  
  if (length(matching_row) > 0)
  {
    filtered_rows <- rbind(filtered_rows, network_used[i, ])
  }
}
network_used <-as.data.frame(filtered_rows[-1,])
network_used <- as.data.frame(network_used)

all_nodes <- unique(c(network_used$gene1, network_used$gene2))

# 创建一个空白的邻接矩阵
adj_matrix <- matrix(0, nrow = length(all_nodes), ncol = length(all_nodes), dimnames = list(all_nodes, all_nodes))

# 将边的信息填充到邻接矩阵中
for (i in 1:nrow(network_used)) {
  source_node <- network_used$gene1[i]
  target_node <- network_used$gene2[i]
  
  # 无向图，所以对称地设置连接关系
  adj_matrix[source_node, target_node] <- 1
  adj_matrix[target_node, source_node] <- 1
}
adj_matrix <- as.data.frame(adj_matrix)

############################ 选择gene all的count ###############################

HC_count <- right_join(HC, gene_all,by = c("...1"="gene"))
ci_count <- right_join(ci, gene_all,by = c("...1"="gene"))
LDN_count <- right_join(LDN, gene_all,by = c("...1"="gene"))
HDN_count <- right_join(HDN, gene_all,by = c("...1"="gene"))
veHCC_count <- right_join(veHCC, gene_all, by = c("...1"="gene"))
eHCC_count <- right_join(eHCC, gene_all, by = c("...1"="gene"))
aHCC_count <- right_join(aHCC, gene_all, by = c("...1"="gene"))
vaHCC_count <- right_join(vaHCC, gene_all, by = c("...1"="gene"))

a <- HC_count[,1]
HC_count <- HC_count[,-c(1,12)]
rownames(HC_count) <- t(a)
ci_count <- ci_count[,-c(1,12)]
rownames(ci_count) <- t(a)
LDN_count <- LDN_count[,-c(1,12)]
rownames(LDN_count) <- t(a)
HDN_count <- HDN_count[,-c(1,9)]
rownames(HDN_count) <- t(a)
veHCC_count <- veHCC_count[,-c(1,10)]
rownames(veHCC_count) <- t(a)
eHCC_count <- eHCC_count[,-c(1,12)]
rownames(eHCC_count) <- t(a)
aHCC_count <- aHCC_count[,-c(1,9)]
rownames(aHCC_count) <- t(a)
vaHCC_count <- vaHCC_count[,-c(1,12)]
rownames(vaHCC_count) <- t(a)

count <- cbind(HC_count,ci_count,LDN_count,HDN_count,veHCC_count,eHCC_count,aHCC_count,vaHCC_count)
rownames(count) <- t(a)
write.csv(count,"count.csv")
#########################  PCA-CMI #############################################

HC_CMI <- PCA_CMI(HC_count,lamda = 0.05,order0 = 1,G = adj_matrix)
write.csv(HC_CMI[["G"]],"HC_CMI.csv")
ci_CMI <- PCA_CMI(ci_count,lamda = 0.05,order0 = 1,G = adj_matrix)
write.csv(ci_CMI[["G"]],"ci_CMI.csv")
LDN_CMI <- PCA_CMI(LDN_count,lamda = 0.05,order0 = 1,G = adj_matrix)
write.csv(LDN_CMI[["G"]],"LDN_CMI.csv")
HDN_CMI <- PCA_CMI(HDN_count,lamda = 0.05,order0 = 1,G = adj_matrix)
write.csv(HDN_CMI[["G"]],"HDN_CMI.csv")
veHCC_CMI <- PCA_CMI(veHCC_count, lamda = 0.05, order0 = 1, G = adj_matrix)
write.csv(veHCC_CMI[["G"]],"veHCC_CMI.csv")
eHCC_CMI <- PCA_CMI(eHCC_count, lamda = 0.05, order0 = 1, G = adj_matrix)
write.csv(eHCC_CMI[["G"]],"eHCC_CMI.csv")
aHCC_CMI <- PCA_CMI(aHCC_count, lamda = 0.05, order0 = 1, G = adj_matrix)
write.csv(aHCC_CMI[["G"]],"aHCC_CMI.csv")
vaHCC_CMI <- PCA_CMI(vaHCC_count, lamda = 0.05, order0 = 1, G = adj_matrix)
write.csv(vaHCC_CMI[["G"]],"vaHCC_CMI.csv")

########################  读取CMI  #############################################
a <- HC_CMI[,1]
HC_CMI <- HC_CMI[,-1]
rownames(HC_CMI) <- t(a)
ci_CMI <- ci_CMI[,-1]
rownames(ci_CMI) <- t(a)
LDN_CMI <- LDN_CMI[,-1]
rownames(LDN_CMI) <- t(a)
HDN_CMI <- HDN_CMI[,-1]
rownames(HDN_CMI) <- t(a)
veHCC_CMI <- veHCC_CMI[,-1]
rownames(veHCC_CMI) <- t(a)
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
write.csv(edge_list,"edge_list_stage1.csv")

edge_list <- matrix(nrow = 0, ncol = 2)
for (i in 1:num_nodes) {
  for (j in (i):num_nodes) {
    if (ci_CMI[i, j] == 1) {
      edge_list <- rbind(edge_list, c(rownames(ci_CMI)[i], colnames(ci_CMI)[j]))  # 减1以匹配R中的索引从1开始
    }
  }
}
write.csv(edge_list,"edge_list_stage2.csv")

edge_list <- matrix(nrow = 0, ncol = 2)
for (i in 1:num_nodes) {
  for (j in (i):num_nodes) {
    if (LDN_CMI[i, j] == 1) {
      edge_list <- rbind(edge_list, c(rownames(LDN_CMI)[i], colnames(LDN_CMI)[j]))  # 减1以匹配R中的索引从1开始
    }
  }
}
write.csv(edge_list,"edge_list_stage3.csv")

edge_list <- matrix(nrow = 0, ncol = 2)
for (i in 1:num_nodes) {
  for (j in (i):num_nodes) {
    if (HDN_CMI[i, j] == 1) {
      edge_list <- rbind(edge_list, c(rownames(HDN_CMI)[i], colnames(HDN_CMI)[j]))  # 减1以匹配R中的索引从1开始
    }
  }
}
write.csv(edge_list,"edge_list_stage4.csv")

edge_list <- matrix(nrow = 0, ncol = 2)
for (i in 1:num_nodes) {
  for (j in (i):num_nodes) {
    if (veHCC_CMI[i, j] == 1) {
      edge_list <- rbind(edge_list, c(rownames(veHCC_CMI)[i], colnames(veHCC_CMI)[j]))  # 减1以匹配R中的索引从1开始
    }
  }
}
write.csv(edge_list,"edge_list_stage5.csv")

edge_list <- matrix(nrow = 0, ncol = 2)
for (i in 1:num_nodes) {
  for (j in (i):num_nodes) {
    if (eHCC_CMI[i, j] == 1) {
      edge_list <- rbind(edge_list, c(rownames(eHCC_CMI)[i], colnames(eHCC_CMI)[j]))  # 减1以匹配R中的索引从1开始
    }
  }
}
write.csv(edge_list,"edge_list_stage6.csv")

edge_list <- matrix(nrow = 0, ncol = 2)
for (i in 1:num_nodes) {
  for (j in (i):num_nodes) {
    if (aHCC_CMI[i, j] == 1) {
      edge_list <- rbind(edge_list, c(rownames(aHCC_CMI)[i], colnames(aHCC_CMI)[j]))  # 减1以匹配R中的索引从1开始
    }
  }
}
write.csv(edge_list,"edge_list_stage7.csv")

edge_list <- matrix(nrow = 0, ncol = 2)
for (i in 1:num_nodes) {
  for (j in (i):num_nodes) {
    if (vaHCC_CMI[i, j] == 1) {
      edge_list <- rbind(edge_list, c(rownames(vaHCC_CMI)[i], colnames(vaHCC_CMI)[j]))  # 减1以匹配R中的索引从1开始
    }
  }
}
write.csv(edge_list,"edge_list_stage8.csv")

############################### 计算pcc ########################################

HC_count <-t(HC_count)
HC_count  <- log(HC_count +1)
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

veHCC_count <-t(veHCC_count)
matPCC_veHCC <- cor(veHCC_count, method='pearson')
diag(matPCC_veHCC)
se_net_veHCC <- matrix(NA, ncol = 3, nrow = nrow(matPCC_veHCC) * ncol(matPCC_veHCC))
colnames(se_net_veHCC) <- c("gene1", "gene2", "weight")
k <- 1
for (i in rownames(matPCC_veHCC)) {
  for (j in colnames(matPCC_veHCC)) {
    row_name <- i
    col_name <- j
    value <- matPCC_veHCC[i, j]
    se_net_veHCC[k, ] <- c(row_name, col_name, value)
    k <- k + 1
  }
}
se_net_veHCC <- as.data.frame(se_net_veHCC)


eHCC_count <-t(eHCC_count)
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

edge_list_stage8_removed <- as.data.frame(edge_list_stage8_removed)
edge_list_stage8_removed_1 <- cbind(edge_list_stage8_removed,paste(edge_list_stage8_removed[, 1],edge_list_stage8_removed[, 2], sep = ""))
colnames(edge_list_stage8_removed_1)[3] <- "combined"

################################################################################
setwd('/home/wt/hcc_add_kegg/add_pcc')
se_net_S1_SE <- cbind(se_net_HC, paste(se_net_HC[, 1], se_net_HC[, 2], sep = ""))
colnames(se_net_S1_SE)[4] <- "combined1"
se_net_S1_SE_1 <-right_join(se_net_S1_SE,edge_list_stage1_removed_1,by=c("combined1"="combined"))
se_net_S1_SE_1 <-se_net_S1_SE_1[,c(1:3)]
colnames(se_net_S1_SE_1)<- c("gene1","gene2","weight")
se_net_S1_SE_1 <- na.omit(se_net_S1_SE_1)
se_net_S1_SE_1 <- rbind(se_net_S1_SE_1,filtered_data)
g <- graph_from_data_frame(se_net_S1_SE_1, directed = TRUE)
g_undirected <- as.undirected(g, mode = "collapse", edge.attr.comb = "mean")
g_undirected <- simplify(g_undirected, remove.multiple = FALSE, remove.loops = TRUE)
se_net_S1_SE_1 <- as_data_frame(g_undirected, what = "edges")
write.csv(se_net_S1_SE_1,"stage1.csv")

se_net_S2_SE <- cbind(se_net_ci, paste(se_net_ci[, 1], se_net_ci[, 2], sep = ""))
colnames(se_net_S2_SE)[4] <- "combined1"
se_net_S2_SE_1 <-right_join(se_net_S2_SE,edge_list_stage2_removed_1,by=c("combined1"="combined"))
se_net_S2_SE_1 <-se_net_S2_SE_1[,c(1:3)]
colnames(se_net_S2_SE_1)<- c("gene1","gene2","weight")
se_net_S2_SE_1 <- na.omit(se_net_S2_SE_1)
se_net_S2_SE_1 <- rbind(se_net_S2_SE_1,filtered_data)
g <- graph_from_data_frame(se_net_S2_SE_1, directed = TRUE)
g_undirected <- as.undirected(g, mode = "collapse", edge.attr.comb = "mean")
se_net_S2_SE_1 <- as_data_frame(g_undirected, what = "edges")
write.csv(se_net_S2_SE_1,"stage2.csv")

se_net_S3_SE <- cbind(se_net_LDN, paste(se_net_LDN[, 1], se_net_LDN[, 2], sep = ""))
colnames(se_net_S3_SE)[4] <- "combined1"
se_net_S3_SE_1 <-right_join(se_net_S3_SE,edge_list_stage3_removed_1,by=c("combined1"="combined"))
se_net_S3_SE_1 <-se_net_S3_SE_1[,c(1:3)]
colnames(se_net_S3_SE_1)<- c("gene1","gene2","weight")
se_net_S3_SE_1 <- na.omit(se_net_S3_SE_1)
se_net_S3_SE_1 <- rbind(se_net_S3_SE_1,filtered_data)
g <- graph_from_data_frame(se_net_S3_SE_1, directed = TRUE)
g_undirected <- as.undirected(g, mode = "collapse", edge.attr.comb = "mean")
se_net_S3_SE_1 <- as_data_frame(g_undirected, what = "edges")
write.csv(se_net_S3_SE_1,"stage3.csv")

se_net_S4_SE <- cbind(se_net_HDN, paste(se_net_HDN[, 1], se_net_HDN[, 2], sep = ""))
colnames(se_net_S4_SE)[4] <- "combined1"
se_net_S4_SE_1 <-right_join(se_net_S4_SE,edge_list_stage4_removed_1,by=c("combined1"="combined"))
se_net_S4_SE_1 <-se_net_S4_SE_1[,c(1:3)]
colnames(se_net_S4_SE_1)<- c("gene1","gene2","weight")
se_net_S4_SE_1 <- na.omit(se_net_S4_SE_1)
se_net_S4_SE_1 <- rbind(se_net_S4_SE_1,filtered_data)
g <- graph_from_data_frame(se_net_S4_SE_1, directed = TRUE)
g_undirected <- as.undirected(g, mode = "collapse", edge.attr.comb = "mean")
se_net_S4_SE_1 <- as_data_frame(g_undirected, what = "edges")
write.csv(se_net_S4_SE_1,"stage4.csv")

se_net_S5_SE <- cbind(se_net_veHCC, paste(se_net_veHCC[, 1], se_net_veHCC[, 2], sep = ""))
colnames(se_net_S5_SE)[4] <- "combined1"
se_net_S5_SE_1 <-right_join(se_net_S5_SE,edge_list_stage5_removed_1,by=c("combined1"="combined"))
se_net_S5_SE_1 <-se_net_S5_SE_1[,c(1:3)]
colnames(se_net_S5_SE_1)<- c("gene1","gene2","weight")
se_net_S5_SE_1 <- na.omit(se_net_S5_SE_1)
se_net_S5_SE_1 <- rbind(se_net_S5_SE_1,filtered_data)
g <- graph_from_data_frame(se_net_S5_SE_1, directed = TRUE)
g_undirected <- as.undirected(g, mode = "collapse", edge.attr.comb = "mean")
se_net_S5_SE_1 <- as_data_frame(g_undirected, what = "edges")
write.csv(se_net_S5_SE_1,"stage5.csv")

se_net_S6_SE <- cbind(se_net_eHCC, paste(se_net_eHCC[, 1], se_net_eHCC[, 2], sep = ""))
colnames(se_net_S6_SE)[4] <- "combined1"
se_net_S6_SE_1 <-right_join(se_net_S6_SE,edge_list_stage6_removed_1,by=c("combined1"="combined"))
se_net_S6_SE_1 <-se_net_S6_SE_1[,c(1:3)]
colnames(se_net_S6_SE_1)<- c("gene1","gene2","weight")
se_net_S6_SE_1 <- na.omit(se_net_S6_SE_1)
se_net_S6_SE_1 <- rbind(se_net_S6_SE_1,filtered_data)
g <- graph_from_data_frame(se_net_S6_SE_1, directed = TRUE)
g_undirected <- as.undirected(g, mode = "collapse", edge.attr.comb = "mean")
se_net_S6_SE_1 <- as_data_frame(g_undirected, what = "edges")
write.csv(se_net_S6_SE_1,"stage6.csv")

se_net_S7_SE <- cbind(se_net_aHCC, paste(se_net_aHCC[, 1], se_net_aHCC[, 2], sep = ""))
colnames(se_net_S7_SE)[4] <- "combined1"
se_net_S7_SE_1 <-right_join(se_net_S7_SE,edge_list_stage7_removed_1,by=c("combined1"="combined"))
se_net_S7_SE_1 <-se_net_S7_SE_1[,c(1:3)]
colnames(se_net_S7_SE_1)<- c("gene1","gene2","weight")
se_net_S7_SE_1 <- na.omit(se_net_S7_SE_1)
se_net_S7_SE_1 <- rbind(se_net_S7_SE_1,filtered_data)
g <- graph_from_data_frame(se_net_S7_SE_1, directed = TRUE)
g_undirected <- as.undirected(g, mode = "collapse", edge.attr.comb = "mean")
se_net_S7_SE_1 <- as_data_frame(g_undirected, what = "edges")
write.csv(se_net_S7_SE_1,"stage7.csv")

se_net_S8_SE <- cbind(se_net_vaHCC, paste(se_net_vaHCC[, 1], se_net_vaHCC[, 2], sep = ""))
colnames(se_net_S8_SE)[4] <- "combined1"
se_net_S8_SE_1 <-right_join(se_net_S8_SE,edge_list_stage8_removed_1,by=c("combined1"="combined"))
se_net_S8_SE_1 <-se_net_S8_SE_1[,c(1:3)]
colnames(se_net_S8_SE_1)<- c("gene1","gene2","weight")
se_net_S8_SE_1 <- na.omit(se_net_S8_SE_1)
se_net_S8_SE_1 <- rbind(se_net_S8_SE_1,filtered_data)
g <- graph_from_data_frame(se_net_S8_SE_1, directed = TRUE)
g_undirected <- as.undirected(g, mode = "collapse", edge.attr.comb = "mean")
se_net_S8_SE_1 <- as_data_frame(g_undirected, what = "edges")
write.csv(se_net_S8_SE_1,"stage8.csv")