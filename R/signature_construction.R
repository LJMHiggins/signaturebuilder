#' @title Markers from coefficients
#'
#' @description
#' Create signatures from coefficients (linear model derived).
#'
#' @param data_files_dir Directory containing regression analysis.
#' @param corr_data_files Character vector containing file names of correlation
#'  data output.
#' @param positive_coef String indicating whether positive coefficients are
#'  indicative of sensitivity or resistance. Options are "resistance" or "sensitivity"
#' @param list_of_filtered_perts List of vectors containing perturbagens which
#'  passed QC thresholds.
#' @param scale Logical, whether to scale beta coefficients before signature
#'  construction.
#'
#' @return
#' @export
#'
#' @examples
make_coef_signatures <- function(data_files_dir,
                                 corr_data_files = NULL,
                                 positive_coef,
                                 list_of_filtered_perts = NULL,
                                 scale =  FALSE,
                                 compound_label = FALSE){

  .make.and.save.marker.df <- function(df, path, nname, q1, q2,
                                       positive_coef){

    perturbagens <- unique(df$perturbagen)
    nr <- length(perturbagens)
    proteins.sensitivity.markers <- character(nr)
    proteins.resistance.markers <- character(nr)
    genes.sensitivity.markers <- character(nr)
    genes.resistance.markers <- character(nr)
    n.sensitivity.markers <- numeric(nr)
    n.resistance.markers <- numeric(nr)
    pers <- character(nr)
    i <- 1

    for (per in perturbagens){
      dfx <- df[df$perturbagen==per,]

      if (positive_coef == "sensitive"){
        prots.sensitivity <- dfx[dfx$beta>q2, "protein"]
        prots.resistance <- dfx[dfx$beta<(q1), "protein"]
      } else if (positive_coef == "resistance"){
        prots.sensitivity <- dfx[dfx$beta<(q1), "protein"]
        prots.resistance <- dfx[dfx$beta>q2, "protein"]
      }
      if (length(prots.resistance)>2 & length(prots.resistance)>2){
        #Build new row
        proteins.resistance.markers[i] <- paste(prots.resistance, collapse = ";")
        proteins.sensitivity.markers[i] <- paste(prots.sensitivity, collapse = ";")

        pers[i] <- per
        n.resistance.markers[i] <- length(prots.resistance)
        n.sensitivity.markers[i] <- length(prots.sensitivity)
        i <- i+1
      }
    }
    df.result <- data.frame(perturbagen = pers,
                            n.resistance.markers,
                            n.sensitivity.markers,
                            proteins.resistance.markers,
                            proteins.sensitivity.markers)

    df.result <- subset(df.result, df.result$perturbagen != "")
    write.csv(df.result, file = paste0(path, nname))
    print(paste0("Finished with ", nname))
  }

  lapply(paste0("Markers_round", 1:3), dir.create)

  if (is.null(corr_data_files)){
    corr_data_files <- list.files(data_files_dir)
    corr_data_files <- corr_data_files[grepl(".csv", corr_data_files)]
    if (any(grepl("corr_values.csv", corr_data_files))){
      f.names <- gsub("corr_values.csv", "markers.csv", corr_data_files)
    } else {
      f.names <- gsub(".csv", "_markers.csv", corr_data_files)
    }
  }

  for (j in 1:length(corr_data_files)){

    f <- corr_data_files[j]
    f <- paste0(data_files_dir, "/", f)
    nname <- f.names[j]
    df_corr <- read.csv(f)
    print(paste0("Reading: ", f))
    if (compound_label == TRUE){
      df_corr <- df_corr %>%
        dplyr::rename(perturbagen = compound)
    }
    if (!is.null(list_of_filtered_perts)){
      good_perts <- list_of_filtered_perts[[j]]
      print(paste0("Got gene list for", names(list_of_filtered_perts[j])))
      df_corr <- df_corr %>%
        filter(perturbagen %in% good_perts)
    }
    if (scale == TRUE){
      df_corr$beta <- scale(df_corr$beta)
    }
    ##### Filter one, save in /Markers_round_1_crispr_2/
    f.prefix <- "Markers_round1/"
    q2 <- quantile(df_corr$beta, na.rm = T)[4]
    q1 <- quantile(df_corr$beta, na.rm = T)[2]
    .make.and.save.marker.df(df = df_corr, path = f.prefix,
                             nname = nname, q1=q1, q2=q2,
                             positive_coef = positive_coef)
    ##### Filter two using top 15% (either side - 30% total)
    f.prefix <- "Markers_round2/"
    q2 <- quantile(df_corr$beta, na.rm = T, probs = c(0.15, 0.85))[2]
    q1 <- quantile(df_corr$beta, na.rm = T, probs = c(0.15, 0.85))[1]
    .make.and.save.marker.df(df = df_corr, path = f.prefix,
                             nname = nname, q1=q1, q2=q2,
                             positive_coef = positive_coef)
    #####Filter 3, 15% total
    f.prefix <- "Markers_round3/"
    q2 <- quantile(df_corr$beta, na.rm = T, probs = c(0.075, 0.925))[2]
    q1 <- quantile(df_corr$beta, na.rm = T, probs = c(0.075, 0.925))[1]
    .make.and.save.marker.df(df = df_corr, path = f.prefix,
                             nname = nname, q1=q1, q2=q2,
                             positive_coef = positive_coef)

    print(paste0("Finished: ", nname))
  }
}


