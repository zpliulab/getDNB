# PCA_CMI
#
# Input------------------------
# dada：num格式，维度 gene X sample
# lamda：互信息阈值
# order0：互信息阶数, 如果空则表示算法一直运行至最高阶！
# G（可选）：不输入则从完全图出发
# 此外，该算法会将数据格式化为离散整数再计算互信息！
#
# Output------------------------
# result$G: 去除假阳性的0-1网络
# result$Gval: 带有边强度的网络
# result$order: 最终停止的CMI阶数
#
# This code is organized by wangchuanyuan1@126.com.

PCA_CMI <- function (data,lamda = 1,order0 = 999,G = NULL)  
{ 
  library(infotheo)  
  n_gene = dim(data)[1]
  if(is.null(G))
  {
    nargin <- 3
    G <- matrix(1, ncol = n_gene, nrow = n_gene)
  }
 # G <- G + t(G)
  G[G>1] = 1
  Gval <- G
  order <- 0;
  t <- 0
  # start PCA_CMI
  while(t==0)
  {
    order <- order+1
    if (order0 < 999)
    {
      if (order>order0)
      {
        result1 <- c()
        result1$G <- G    
        result1$Gval <- Gval
        result1$order <- order-1
        
        return(result1)
      }
    }
    result <- edgereduce(G,Gval,order,data,t,lamda)
    G <- result$G
    Gval <- result$Gval
    t <- result$t
    print(t)
    if(t==0)
    {
      print('No edge is reduce! Algorithm  finished!')
      break
    }else
    {t = 0}
  }
  result1 <- c()
  result1$G <- G    
  result1$Gval <- Gval
  result1$order <- order-1
  return(result1)
}  


edgereduce = function (G,Gval,order,data,t,lamda)
{ 
  # edgereduce is pca_cmi
  # t：剩余未处理的边
  # lamda: threshold
  dis_data <- discretize(data,disc="equalwidth",(NROW(t(data))^1/2))   # equalwidth等宽 equalfreq等频
  if(order==0)   # order为0则表示只计算MI，而不是CMI
  {
    for(i in 1:(dim(G)[1]-1))
    {
      for (j in (i+1):dim(G)[1])
      {
        if (G[i,j]!=0)
        {
          mi <- mutinformation(t(dis_data[i,]),t(dis_data[j,]),method="emp")
          Gval[i,j] <- mi 
          Gval[j,i] <- mi
          if (mi<lamda) # set threshold 
          {
            G[i,j] <- 0
            G[j,i] <- 0
          }
        }
      }
    }
    t <- t+1
  }else       # 计算CMI
  {
    for (i in 1:(dim(G)[1]-1))                        # 对于每个调控者
    {
      for (j in (i+1):dim(G)[1])                      # 对于每个受体
      {
        if (G[i,j]!=0)                                # 如果存在edge
        {
          adj <- c()                                  # 两个节点公共邻居的位置
          for (k in 1:dim(G)[1])
          {
            if (G[i,k]!=0 && G[j,k]!=0)               # i、j、k三点互相都有edge
            {
              adj <- c(adj,k);
            }
          }
          if (length(adj) >= order)                   # 如果邻居数目超过order（即考个order个基因对i和j的影响）
          {
            combntnslist <- t(combn(adj,order));      #   combntns 返回所有的排列组合（组合内数目为order） Use NCHOOSEK instead. 
            combntnsrow <- dim(combntnslist)[1];      #   combntnsrow为组合总数
            cmiv <- 0
            v1 <- dis_data[i,]                        #   调控者的基因表达量
            v2 <- dis_data[j,]                        #   受体的基因表达量
            for (k in (1:combntnsrow))                #   对于每个条件组合
            {
              vcs <- dis_data[combntnslist[k,],];     #       对于每个第三方基因
              a <- condinformation(as.numeric(v1),as.numeric(v2),as.numeric(vcs),method="emp")               #       计算i、j和第三方基因的CMI
              cmiv <- max(cmiv,a)                     #       找到所有组合中的最大CMI
            }
            Gval[i,j] <- cmiv                         #   根据第三方基因重新计算edge强度
            Gval[j,i] <- cmiv                         #   根据第三方基因重新计算edge强度
            if (cmiv<lamda)                           #   根据阈值删除边
            {
              G[i,j] <- 0                             #       删除边
              G[j,i] <- 0
            }              
            t <- t+1                                    # 运行次数+1
          }
        }               
      }
    }             
  }
  result <- c()
  result$G <- G
  result$Gval <- Gval
  result$t <- t
  return(result)
}

tril = function (G, flag = -1) 
  # 这是提取矩阵上三角的函数
{
  if (flag == -1)
  {G <- t(G)}
  m <- dim(G)
  for(i in 1: m[1])
  {
    for(j in i: m[2])
    {G[i,j] <- 0}}
  G <- t(G)
  return(G)
}



