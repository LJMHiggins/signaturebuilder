## Modified corr two datasets function from protools github ###

#' Correlate two datasets
#'
#' Function taken from protools2 and slightly altered, main change is the calculation of
#'  R2 value. Protein df scaling can be performed outside of funtion.
#'
#' @param x_dataset Intended to be independent x variable of linear function. Here,
#'  the proteomics datasets.
#' @param y_dataset Intended to be the dependent y variable of linear function. Here,
#'  the perturbation datasets.
#' @param common_column Cell line ids, for merging.
#' @param y_factors Perturbations.
#' @param x_factors Proteins
#' @param y_factor_name
#' @param x_factor_name
#' @param protein_scale Logical. whether to scale protein before fitting.
#' @param threshold_test String. Determines initial test to determine the dependence
#'  between the two variables.
#'
#' @return
#' @export
#'
#' @examples
corr.two.data.sets.v2 <- function(
    x_dataset,
    y_dataset,
    common_column = "cell.line",
    y_factors,
    x_factors,
    y_factor_name = "y",
    x_factor_name = "x",
    threshold_test = "Corr"){

  library(data.table)
  library(doParallel)

  y_factors <- intersect(y_factors, colnames(y_dataset))
  x_factors <- intersect(x_factors, colnames(x_dataset))

  ncores <- parallel::detectCores()-1
  cl <- makeCluster(ncores)
  df.results <- data.frame(matrix(nrow = 0, ncol = 8))

  for(per in y_factors){
    start_time <- Sys.time()
    registerDoParallel(cl)
    y <- y_dataset[, c(common_column, per)]
    colnames(y)[2] <- "per"
    xx <- merge.data.frame(y, x_dataset, by = common_column)
    dfpp <- foreach(prot = x_factors, .combine = 'rbind')%dopar%{
      tryCatch({
        dfx <- na.omit(xx[, c("per", prot)])
        #pval <- 1
        #rval <- 0
        slope <- 0
        r_squared <- 0
        if (!is.null(dfx)){
          if (nrow(dfx)>8){
            if (threshold_test == "Corr"){
              tt <- cor.test(dfx[, 1], dfx[, 2], method = "pearson")
              test_pval <- tt$p.value
              test_stat <- tt$estimate
            } else if (threshold_test == "MI"){
              mi <- fastmit::mi.test(x = dfx[, 1], y = dfx[, 2])
              test_pval <- mi$p.value
              test_stat <- mi$estimate
            }
            if (pval<0.05){
              rval <- tt$estimate
              ll <- lm(dfx[, 1] ~ dfx[, 2])
              slope <- ll$coefficients[2]
              r_squared <- summary(ll)$r.squared

              return(c(per,prot, test_pval, test_stat, slope, r_squared, nrow(dfx)))
            }
          }
        }
      }, error = function(e){}
      )
    }
    stopImplicitCluster()
    end_time <- Sys.time()
    print(paste(per, end_time-start_time))

    if (exists('dfpp')){
      if (!is.null(dfpp)){
        dfpp <- as.data.frame(dfpp)
        dfpp <- dfpp[dfpp[, 3] !=1, ]
        dfpp$fdr <- p.adjust(as.numeric(dfpp[, 3]), method = "fdr")
        df.results <- rbind.data.frame(df.results, dfpp)
      }
    }
  }
  if (threshold_test == "Corr"){
    colnames(df.results) <- c(y_factor_name, x_factor_name, "pval", "rval", "beta", "R_squared", "n_cells", "fdr")
    df.results <- subset(df.results, df.results$rval!=0)
  } else if (threshold_test == "MI"){
    colnames(df.results) <- c(y_factor_name, x_factor_name, "pval", "mi", "beta", "R_squared", "n_cells", "fdr")
  }
  df.results <- df.results[order(df.results$fdr), ]

  return(df.results)
}
