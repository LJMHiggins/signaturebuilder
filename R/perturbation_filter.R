#' Compute perturbagen summary statistics
#'
#' @param pert.datasets
#' @param prot.datasets
#'
#' @return
#' @export
#'
#' @examples
perturbation_summary_stats <- function(pert.datasets, prot.datasets){

  library(dplyr)
  library(doParallel)
  library(multimode) #Test modality
  library(moments)

  i <- 1
  for (per.df in pert.datasets) {
    pert.df.name <- names(pert.datasets)[i]
    i <- i + 1
    j <- 1

    perturbations <- colnames(per.df)[-1]

    for (prot.df in prot.datasets){
      prot.df.name <- names(prot.datasets)[j]
      j <- j + 1
      shared.cells <- intersect(per.df[,"cell.line"], prot.df[,"cell.line"])
      matched.pert.df <- per.df[per.df$cell.line %in% shared.cells, ]
      matched.yvar.stats <- foreach::foreach(perturb = perturbations, .combine = "rbind") %do% {
        # Not all "good.genes" are in pert cols
        if (perturb %in% colnames(matched.pert.df)){
          yvar <- matched.pert.df[,c("cell.line",perturb)]
          # Compute summary stats
          if (nrow(yvar %>% na.omit()) > 3) {

            sum.stats <- summary(yvar[,2], na.rm = TRUE)
            n.number <- nrow(yvar %>% na.omit())
            mean.val<- sum.stats[[4]]
            st.dev <- sd(yvar[,2], na.rm = TRUE)
            min.val <- sum.stats[[1]]
            Q1 <- sum.stats[[2]]
            med.val <- sum.stats[[3]]
            Q3 <- sum.stats[[5]]
            max.val <- sum.stats[[6]]
            int.q.range <- Q3 - Q1
            kurt <- moments::kurtosis(yvar[,2] %>% na.omit())
            skew <- moments::skewness(yvar[,2] %>% na.omit())

            tryCatch(
              {
                m.test <- multimode::modetest(yvar[,2] %>% na.omit())
                m.test.pval <- m.test$p.value
              }, error=function(e) m.test.pval <- NA
            )

            out <- data.frame("Perturbagen" = perturb,
                              "Sample No. with data" = n.number,
                              "Mean" = mean.val,
                              "Standard dev" = st.dev,
                              "Mininum val" = min.val,
                              "Q1" = Q1,
                              "Median" = med.val,
                              "Q3" = Q3,
                              "Maximum.val"= max.val,
                              "IQR" = int.q.range,
                              "Kurtosis" = kurt,
                              "Skew" = skew,
                              "ACR multimodality test pval" = m.test.pval)
            return(out)
          }
        }
      }
      f.name <- paste0(pert.df.name, "_", prot.df.name, "_summary_stats.csv")
      results_dir <- "computed_perturbagen_stats"
      dir.create(results_dir)
      write.csv(matched.yvar.stats,
                file = paste0(results_dir, "/", f.name))
      print(paste("Completed stats for:", f.name))
    }
  }
}

#' Filter perturbations on summary statistics
#'
#' @param summary_stats_df
#' @param gather
#' @param min_val_cut
#' @param max_val_cut
#' @param sd_cut
#' @param iqr_cut
#' @param skew_cut
#' @param kurt_cut
#'
#' @return
#' @export
#'
#' @examples
perturbation_filter <- function(summary_stats_df,
                                gather = TRUE,
                                min_val_cut,
                                max_val_cut,
                                sd_cut,
                                iqr_cut,
                                skew_cut,
                                kurt_cut){
  if (gather == TRUE){
    # Get summary stats for all datasets
    pert.sts.fnames <- list.files(path = "computed_perturbagen_stats/")
    long.df <- data.frame()
    for (f in pert.sts.fnames){
      if (grepl("summary_stats", f)) {
        df <- read.csv(paste0("Computed_perturbagen_stats/",f))
        d.set <- gsub("_summary_stats.csv","", f)
        df$dataset <- rep(d.set,nrow(df))
        long.df <- rbind(long.df, df)
      }
    }
  }
  # Subset perturbation statistics dataframe
  long.df %>%
    group_by(dataset) %>%
    filter(Mininum.val < 0.5 & Maximum.val >= 0.5) %>%
    filter((Standard.dev > 0.15) | (IQR > 0.05 &
                                      (
                                        (
                                          Skew >= 7 | Skew <= -7) & Kurtosis >= 2
                                      ))) -> df.filtered
  # Generate list of perturbations which passed QC thresholds
  datasets <- unique(df.filtered$dataset)
  list_filtered_pertubations <- list()
  for (dset in datasets){
    df.d <- df.filtered %>%
      filter(dataset==d.set)
    perturbations <- unlist(df.d$Perturbagen)
    list_filtered_pertubations[[d.set]] <- perturbations
  }
  return(list_filtered_pertubations)
}
