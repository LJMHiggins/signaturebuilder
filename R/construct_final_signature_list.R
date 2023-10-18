#' Title
#'
#' @param markers_dir Directory leading to marker sheet subfolders (for different
#'  stringency levels).
#' @param all_stringencies Logical. Whether 3 marker stringency levels has been implemented.
#' @param marker_type String. Indicate whether markers are for crispr, rnai or drug response.
#' @param signature_id_suffix String. Id for signature, such as metric used in regression analysis.
#'
#' @return List containing identified markers of response (resistance/sensitivity) for
#'  each perturbation ~ proteomics dataset combo.
#' @export
#'
#' @examples
signature_list_from_markers <- function(markers_dir,
                                        all_stringencies = TRUE,
                                        marker_type = NULL,
                                        signature_id_suffix){

  marker.folders <- list.files(markers_dir)[grepl("Markers_round", list.files(markers_dir))]

  if (all_stringencies == TRUE){
    i <- 1
    for (r in c("Markers_round1", "Markers_round2", "Markers_round3")){
      marker.folders <- marker.folders[grepl(r, marker.folders)]
      list_of_signatures <- list()
      mdir <- paste0(markers_dir, "/", r)
      marker.df.files <- list.files(mdir)

      for (m.df in marker.df.files){
        f <- paste0(mdir, "/", m.df)

        xx <- read.csv(f)
        xx$perturbagen <- gsub("X.", "", xx$perturbagen,fixed = T)
        # Create signature slot name
        sig_name <- m.df
        sig_name <- gsub("markers.csv", paste0("signatures_", signature_id_suffix), sig_name)
        sig_name <- gsub("depmap_", "", sig_name)
        sig_name <- gsub("pharmaco_", "", sig_name)
        sig_name <- paste0(marker_type, "_", sig_name)
        # Append list
        list_of_signatures[[sig_name]] <-  xx[,2:ncol(xx)]
      }
      sig_file_name <- paste0("list_of_signatures_", signature_id_suffix, "_level", i)
      save(list_of_signatures, file = paste0(sig_file_name, ".rda"))
      i <- i + 1
    }
  }
}
