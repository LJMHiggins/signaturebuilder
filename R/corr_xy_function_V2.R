## Corr two datasets function from protools github ###

corr.two.data.sets.v2_NOTREADY <- function(
    x_dataset,
    y_dataset,
    common_column="cell.line",
    y_factors,
    x_factors,
    y_factor_name="y",
    x_factor_name="x",
    n_cores = 4){

  library(data.table)
  library(doParallel)


  y_factors <- intersect(y_factors,colnames(y_dataset))
  x_factors <- intersect(x_factors,colnames(x_dataset))

  df.results <- data.frame(matrix(nrow = 0,ncol=7))

  n_cores <- n_cores
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)
  clusterEvalQ(cl, {
    library(drc)
    library(dplyr)
    library(tidyr)
    library(foreach)
  })

  foreach::foreach(per = y_factors, .combine = "rbind")%dopar%{
    y <- y_dataset[,c(common_column,per)]
    colnames(y)[2] <- "per"
    xx <- merge.data.frame(y,x_dataset,by=common_column)

    dfpp <- foreach(prot = x_factors, .combine = "rbind")%do%{

      tryCatch({
        dfx <- na.omit(xx[,c("per",prot)])
        pval <- 1
        rval <- 0
        slope <- 0
        if (!is.null(dfx)){
          if (nrow(dfx)>8){
            tt <- cor.test(dfx[,1],dfx[,2],method = "pearson")
            pval <- tt$p.value
            if (pval<0.05){
              rval <-  tt$estimate
              ll <- lm.fit(cbind(1,dfx[,1]),dfx[,2])
              slope <-  ll$coefficients[2]
              return(c(per,prot, pval, rval, slope, nrow(dfx)))
            }
          }
        }
      }, error = function(e){}
      )
    }
    if (exists('dfpp')){
      if (!is.null(dfpp)){
        dfpp <- as.data.frame(dfpp)
        dfpp <- dfpp[dfpp[,3] !=1,]
        dfpp$fdr <- p.adjust(as.numeric(dfpp[,3]),method = "fdr")
        df.results <- rbind.data.frame(df.results,dfpp)
      }
    }
    return(df.results)
  }

  stopCluster(cl)
  registerDoSEQ()


  colnames(df.results) <- c(y_factor_name,x_factor_name,"pval","rval", "beta", "n_cells","fdr")
  df.results <- subset(df.results,df.results$rval!=0)
  df.results <- df.results[order(df.results$fdr),]

  return(df.results)
}
