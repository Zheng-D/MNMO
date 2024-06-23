


#' MNMO scores
#'
#' @param net_gene Genes in the second layer network
#' @param network Layer-two Network
#' @param cancer_sample Gene expression of cancer samples, row names are genes, column names are samples
#' @param G_mut Gene mutation matrix, row names are samples, column names are genes
#' @param mirna_value Differentially expressed miRNAs and their scores
#' @param mir_gene_network miRNA and gene interaction network
#' @param gene_expscore_2 is the output file of construct_net function, which is the part1 of the second-layer gene control ability score
#' @param gene_expscore_3 is the output file of construct_net function, which is the part2 of the second-layer gene control ability score
#' @import psych
#' @return Ranking of potential driver genes
#' @export
#' @importFrom stats cor.test


MNMO <- function(net_gene,network,cancer_sample,G_mut,mirna_value,mir_gene_network,gene_expscore_2,gene_expscore_3){

  net_edges <- network
  mrna_exp <- cancer_sample
  M1 <- matrix(0,nrow = length(net_gene),ncol = length(net_gene))
  rownames(M1) <- net_gene
  colnames(M1) <- net_gene
  for (i in 1:nrow(M1)) {
    neib_gene <- which(net_edges[,1]==i)
    neib_list <- net_edges[neib_gene,2]
    M1[i,neib_list] <- 1
  }
  gene_exp1 <- matrix(data=0,nrow=nrow(mrna_exp),ncol=2)
  gene_exp1[,1] <- rownames(mrna_exp)
  M2 <- M1
  for(i in 1:nrow(M2))
  {
    neib_list <- which(M2[i,]==1)
    gene1_row <- which(gene_exp1[,1]==net_gene[i])
    if(length(gene1_row)==0){
      M2[i,] <- 0
      M2[,i] <- 0
      next()
    }
    for(j in 1:length(neib_list))
    {
      gene2_index <- neib_list[j]
      gene2_row <- which(gene_exp1[,1]==net_gene[gene2_index])
      if(length(gene2_row)==0){
        M2[i,gene2_index] <- 0
        M2[gene2_index,i] <- 0
        next()
      }
      x <- as.numeric(mrna_exp[gene1_row,])
      y <- as.numeric(mrna_exp[gene2_row,])
      T<- cor.test(x,y)
      if(T$p.value<0.05)
      {
        M2[i,gene2_index] <- abs(as.numeric(T$estimate))
        M2[gene2_index,i] <- abs(as.numeric(T$estimate))
      }
      else{
        M2[i,gene2_index] <- 0
        M2[gene2_index,i] <- 0
      }
    }
  }
  Mz <- M2
  for (i in 1:nrow(Mz)) {
    pcc_value <- 2-sqrt(3)
    del0 <- which(Mz[i,]< pcc_value)
    Mz[i,del0] <- 0
  }

  #ecc
  Mz_pcc <- matrix(0,nrow = length(net_gene),ncol = length(net_gene))
  rownames(Mz_pcc) <- net_gene
  colnames(Mz_pcc) <- net_gene
  for (i in 1:nrow(Mz)) {
    gene1_nei <- which(Mz[i,]!=0)
    if(length(gene1_nei)==0)next
    for (j in 1:length(gene1_nei)) {
      gene2_index <- gene1_nei[j]
      gene2_neib <- which(Mz[gene2_index,]!=0)
      inter_gene  <- intersect(gene1_nei,gene2_neib)
      if(length(inter_gene)==0)next
      l_value <- 0
      for (k in 1:length(inter_gene)) {
        gene3_index <- inter_gene[k]
        x <- Mz[i,gene3_index]
        y <- Mz[gene2_index,gene3_index]
        pcc_value <- 1- log2((x^2+y^2)/(2*x*y))
        l_value <- l_value + pcc_value
      }
      Mz_pcc[i,gene2_index] <- l_value/length(gene1_nei)
    }
  }

  #
  M_matrix_weight <- Mz_pcc
  gene_score <- matrix(0,nrow = length(net_gene),ncol = 5)
  gene_score[,1] <- net_gene
  for (i in 1:nrow(Mz_pcc)) {
    gene_score[i,2] <-sum(as.numeric(Mz_pcc[i,]))
  }
  for(i in 1:nrow(Mz_pcc))
  {
    sum1<- as.numeric(gene_score[i,2])
    if(sum1==0)next()
    neib_list1 <- which(Mz_pcc[i,]!=0)
    for(j in 1:length(neib_list1))
    {
      gene_index1 <- neib_list1[j]
      sum2 <- as.numeric(gene_score[gene_index1,2])
      M_matrix_weight[i,gene_index1] <- abs(Mz_pcc[i,gene_index1])/sqrt(sum1*sum2)
    }
  }
  M_matrix1 <- matrix(data=0,nrow=nrow(Mz_pcc),ncol=1)
  for(i in 1:length(net_gene))
  {
    index_find <- which(colnames(G_mut)==net_gene[i])
    gene_score[i,3] <- sum(G_mut[,index_find])/nrow(G_mut)
  }
  for(i in 1:nrow(Mz))
  {
    M_matrix1[i,1] <- as.numeric(gene_score[i,3])
  }

  Exp_initial <- M_matrix1
  rownames(Exp_initial) <- gene_score[,1]
  beta <- 0.4
  count2=1
  Exp_tem <- (1-beta)*M_matrix_weight%*%Exp_initial + beta*Exp_initial
  Exp_final <- (1-beta)*M_matrix_weight%*%Exp_tem + beta*Exp_initial
  Exp_text <- abs(Exp_final - Exp_tem)
  result_index <- which((Exp_text)>10^-6)
  while(length(result_index)!=0){
    Exp_tem <- Exp_final
    Exp_final <- (1-beta)*M_matrix_weight%*%Exp_tem + beta*Exp_initial
    Exp_text <- abs(Exp_final - Exp_tem)
    result_index <- which((Exp_text)>10^-6)
    count2 <- count2 + 1
    #print(Exp_final[1])
  }
  for (i in 1:nrow(gene_score)) {

    gene_score[i,4] <- Exp_final[i,1]
  }

  ##net
  for (i in 1:nrow(M1)) {
    gene_score[i,5] <- sum(as.numeric(M1[i,]))
  }
  result_gene_list <- matrix(0,nrow = nrow(gene_score),9)
  result_gene_list[,1] <- gene_score[,1]
  result_gene_list[,2] <- gene_score[,4]
  result_gene_list[,8] <- gene_score[,5]
  colnames(mir_gene_network) <- mir_gene_network[1,]
  mir_gene_network <- mir_gene_network[-1,]
  result_gene_list[,3] <- 0
  for (i in 1:nrow(mirna_value)) {
    g_list <- which(mir_gene_network[,2]==mirna_value[i,1])
    if(length(g_list)==0)next
    gn_list <- mir_gene_network[g_list,1]
    for (j in 1:length(gn_list)) {
      r_index <- which(result_gene_list[,1]==gn_list[j])
      if(length(r_index)==0)next
      value1 <- as.numeric(result_gene_list[r_index,3])
      result_gene_list[r_index,3] <- value1 + as.numeric(mirna_value[i,2])
    }
  }

  for (i in 1:nrow(result_gene_list)) {
    gene2_index <- which(gene_expscore_2[,1]==result_gene_list[i,1])
    gene3_index <- which(gene_expscore_3[,1]==result_gene_list[i,1])
    if(length(gene2_index)==0){
      result_gene_list[i,4] <- 0
    }else{
      result_gene_list[i,4] <- gene_expscore_2[gene2_index,2]
    }
    if(length(gene3_index)==0){
      result_gene_list[i,5] <- 0
    }else{
      result_gene_list[i,5] <- gene_expscore_3[gene3_index,2]
    }
    result_gene_list[i,6] <- max(as.numeric(result_gene_list[i,4]),as.numeric(result_gene_list[i,5]))
    result_gene_list[i,7] <- as.numeric(result_gene_list[i,3]) + as.numeric(result_gene_list[i,6])
  }
  result_gene_list[,1] <- 0
  result_gene_list <- apply(result_gene_list,2, as.numeric)
  m_col <-  max(result_gene_list[,2])
  mi_col <-  min(result_gene_list[,2])
  m_exp <- max(result_gene_list[,7])
  mi_exp <- min(result_gene_list[,7])
  m_net <- max(result_gene_list[,8])
  mi_net <- min(result_gene_list[,8])

  for (i in 1:nrow(result_gene_list)) {
    va1 <-  (result_gene_list[i,2] - mi_col)/m_col
    va2 <-  (result_gene_list[i,7] - mi_exp)/m_exp
    va3 <- (result_gene_list[i,8] - mi_net)/m_net
    if(va1*va2*va3!=0){
      vl <- c(va1,va2,va3)
      result_gene_list[i,9] <- harmonic.mean(vl)
    }
  }
  result_gene_list[,1] <- gene_score[,1]
  result_gene_list <- result_gene_list[order(result_gene_list[,9],decreasing = TRUE),]
  return(result_gene_list)
  #write.csv(result_gene_list,"result_gene_list.csv")

}
