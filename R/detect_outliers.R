
#' Detect Outliers in Gene Expression Data using PCA and MDS
#'
#' This function detects outlier samples from gene expression data using
#' Principal Component Analysis (PCA) and Multi-Dimensional Scaling (MDS).
#' Samples deviating more than `n_sd` standard deviations from the mean
#' in any of the two dimensions are flagged as outliers. The detection
#' process is iterative until no additional outliers are found.
#'
#' @param data_matrix A numeric matrix or data frame of expression values
#'   (genes in rows, samples in columns). Row names should be gene IDs
#'   and column names sample IDs.
#' @param n_sd Numeric. Threshold in standard deviations for outlier detection
#'   (default = 5).
#'
#' @return A list containing:
#'   \item{PC_plot}{A ggplot object for PCA results.}
#'   \item{MDS_plot}{A ggplot object for MDS results.}
#'   \item{filter_table}{A data frame with sample IDs and their outlier status ("Outlier" or "Passed").}

#'
#'
#' @import ggplot2
#' @import stats
#' @import patchwork
#' @export
#'
#' @examples
#' \dontrun{
#' result <- detect_outliers_pca_mds(data_matrix, n_sd = 4)
#' }
#'

detect_outliers_pca_mds <- function(data_matrix, n_sd = 5) {


  filter_table <- data.frame(RNAsample = colnames(data_matrix),
                             filter_outliers = NA)

  # PCA function
  calculatePC <- function(matrix) {
    pca_result <- prcomp(matrix, scale. = TRUE)
    pc_values <- as.data.frame(pca_result$x[, 1:2])
    colnames(pc_values) <- c("PC1", "PC2")
    pc_stats <- c(mean1 = mean(pc_values$PC1),
                  sd1 = sd(pc_values$PC1),
                  mean2 = mean(pc_values$PC2),
                  sd2 = sd(pc_values$PC2))
    list(values = pc_values, stats = pc_stats)
  }

  # MDS function
  calculateMDS <- function(matrix) {
    distance_matrix <- dist(matrix)
    mds_result <- cmdscale(distance_matrix, k = 2)
    mds_values <- as.data.frame(mds_result)
    colnames(mds_values) <- c("MDS1", "MDS2")
    mds_stats <- c(mean1 = mean(mds_values$MDS1),
                   sd1 = sd(mds_values$MDS1),
                   mean2 = mean(mds_values$MDS2),
                   sd2 = sd(mds_values$MDS2))
    list(values = mds_values, stats = mds_stats)
  }

  # Outlier detection
  i <- 1
  matrix <- t(data_matrix)
  matrix <- matrix[, colSums(matrix) > 0]
  matrix <- matrix[, colSums(matrix >= 10) >= 5]

  pc_out <- calculatePC(matrix)
  mds_out <- calculateMDS(matrix)
  pc_ini <- pc_out
  mds_ini <- mds_out

  pc_outliers <- (abs(pc_out$values[, 1] - pc_out$stats["mean1"]) > n_sd * pc_out$stats["sd1"]) |
    (abs(pc_out$values[, 2] - pc_out$stats["mean2"]) > n_sd * pc_out$stats["sd2"])
  mds_outliers <- (abs(mds_out$values[, 1] - mds_out$stats["mean1"]) > n_sd * mds_out$stats["sd1"]) |
    (abs(mds_out$values[, 2] - mds_out$stats["mean2"]) > n_sd * mds_out$stats["sd2"])

  outliers <- unique(c(rownames(matrix)[pc_outliers], rownames(matrix)[mds_outliers]))
  discardedSamples <- outliers

  while (length(outliers) > 0) {
    i <- i + 1
    matrix <- matrix[!rownames(matrix) %in% outliers, ]
    matrix <- matrix[, colSums(matrix) > 0]
    matrix <- matrix[, colSums(matrix >= 10) >= 5]

    pc_out <- calculatePC(matrix)
    mds_out <- calculateMDS(matrix)

    pc_outliers <- (abs(pc_out$values[, 1] - pc_out$stats["mean1"]) > n_sd * pc_out$stats["sd1"]) |
      (abs(pc_out$values[, 2] - pc_out$stats["mean2"]) > n_sd * pc_out$stats["sd2"])
    mds_outliers <- (abs(mds_out$values[, 1] - mds_out$stats["mean1"]) > n_sd * mds_out$stats["sd1"]) |
      (abs(mds_out$values[, 2] - mds_out$stats["mean2"]) > n_sd * mds_out$stats["sd2"])

    outliers <- unique(c(rownames(matrix)[pc_outliers], rownames(matrix)[mds_outliers]))
    discardedSamples <- c(discardedSamples, outliers)
  }

  # Mark outliers
  pc_ini$values[["outliers"]] <- rownames(pc_ini$values) %in% discardedSamples
  mds_ini$values[["outliers"]] <- rownames(mds_ini$values) %in% discardedSamples

  pc_plot <- ggplot(pc_ini$values, aes(x = PC1, y = PC2, color = outliers)) +
    geom_point() + theme_bw() +
    scale_color_manual(values = c("skyblue4", "tomato3")) +
    ggtitle(paste0("PC analysis outliers (detected as ", n_sd, "*SD)"))

  mds_plot <- ggplot(mds_ini$values, aes(x = MDS1, y = MDS2, color = outliers)) +
    geom_point() + theme_bw() +
    scale_color_manual(values = c("skyblue4", "tomato3")) +
    ggtitle(paste0("MDS analysis outliers (detected as ", n_sd, "*SD)"))
  print(pc_plot+mds_plot)
  # Update filter table
  filter_table$filter_outliers <- ifelse(
    filter_table$RNAsample %in% discardedSamples,
    "Outlier", "Passed"
  )


  return(list(PC_plot = pc_plot,
              MDS_plot = mds_plot,
              filter_table = filter_table))
}
