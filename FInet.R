# FI-net main code
FInet <- function(input.matrix,N){
  
  X <- subset(input.matrix,select = -c(mutationScore,gene))
  X <- scale(X)
  Y <- input.matrix$mutationScore
  gene.feature.score <- as.data.frame(cbind(Y,X))
  colnames(gene.feature.score)[1] <- "mutationScore"
  rownames(gene.feature.score) <- input.matrix$gene

  h2o.init()
  
  h2o.input <- as.h2o(gene.feature.score)
  
  bptimes <- 20
  matrixPredY <- data.frame(matrix(data = NA,nrow = length(Y),ncol = bptimes))
  
  for (bp in 1:bptimes){
    model <- h2o.deeplearning(x = 2:ncol(gene.feature.score), # column numbers for predictors
                              y = 1,  # column number for label
                              training_frame = h2o.input,
                              standardize = FALSE,
                              activation = "Rectifier",
                              hidden = c(100), ## one hidden layers
                              epochs = 10,
                              rate = 0.005,
                              seed = 123)
    
    predict.result <- h2o.predict(model,h2o.input)
    predict.y <- unlist(as.data.frame(predict.result))
    matrixPredY[,bp] <- predict.y
  }
  
  meanPredY <- apply(matrixPredY,1,mean)
  gene.feature.score$estimateFI <- meanPredY

  knum <- ceiling(nrow(X)/N)
  hc <- hclust(dist(X,method = "euclidean"), method="ward.D")
  clusters <- cutree(hc, k=knum)
  table(clusters)
  
  M.cluster <- gene.feature.score
  M.cluster[,1:(ncol(gene.feature.score)-1)] <- input.matrix[,2:ncol(input.matrix)]
  M.cluster$cluster <- clusters
  M.total <- as.data.frame(matrix(data = NA,nrow = 0,ncol = ncol(M.cluster)))
  
  for (i in 1:length(unique(clusters))) {
    M.cluster.i <- M.cluster[which(M.cluster$cluster == i),]
    estimate.FIS <- M.cluster$estimateFI[which(M.cluster$cluster == i)]
    obs.FIS <- M.cluster$mutationScore[which(M.cluster$cluster == i)]
    flag.trun <- which(estimate.FIS < quantile(estimate.FIS,0.05))
    estimate.FIS.quantile <- estimate.FIS[-flag.trun]
    if(min(estimate.FIS.quantile) <= 0){
      min.FIS <- min(estimate.FIS.quantile)
      estimate.FIS.quantile <- estimate.FIS.quantile - min.FIS + 0.01
      obs.FIS <- obs.FIS - min.FIS + 0.01
    }
    
    FIS.distribution <- fitdist(estimate.FIS.quantile,"gamma",method = "mle")
    for (j in 1:nrow(M.cluster.i)) {
      if(obs.FIS[j] <= 0){
        M.cluster.i$pValue[j] <- 1
      }else{
        ks.result <- ks.test(obs.FIS[j],"pgamma",FIS.distribution$estimate[1],FIS.distribution$estimate[2],alternative = "less")
        M.cluster.i$pValue[j] <- ks.result$p.value
      }
    }
    
    M.cluster.i$qValue <- stats::p.adjust(M.cluster.i$pValue,method = "fdr",nrow(M.cluster.i))
    M.cluster.i <- M.cluster.i[order(M.cluster.i$qValue,decreasing = F),]
    M.total <- rbind(M.total,M.cluster.i)
  }

  M.total <- M.total[order(M.total$qValue,decreasing = F),]
  M.total <- cbind(rownames(M.total),M.total)
  colnames(M.total)[1] <- "gene"
  M.total <- data.frame(gene = rownames(M.total),qValue = M.total$qValue)
  return(M.total)
} 

