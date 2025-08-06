#' MDS and MAD Plots for Variability Analysis
#'
#' This function generates plots to analyze variability in gene expression data normalized as CPM.
#' It performs a multidimensional scaling (MDS) analysis to visualize the similarity between samples
#' according to their sample type and batch. It also calculates the median absolute deviation (MAD)
#' of gene expression by sample type and by batch.
#'
#' @param cpm Matrix. CPM expression matrix with genes in rows and samples in columns.
#' @param metadata Data frame. Must contain the columns `SampleType` and `Batch`, and row names
#'        must correspond to the samples (identifiers matching the columns of the cpm matrix).
#' @param filter Numeric. Threshold to filter extreme MAD values in the density plots.
#'
#' @return A list containing:
#' \describe{
#'   \item{mds_result}{Matrix with the first two dimensions of the MDS analysis}
#'   \item{plot_mds_group}{MDS plot by sample type}
#'   \item{plot_mds_batch}{MDS plot by batch}
#'   \item{plot_mad_group}{MAD distribution by sample type}
#'   \item{plot_mad_batch}{MAD distribution by batch}
#' }
#' @importFrom matrixStats rowMads
#' @importFrom ggplot2 ggplot aes geom_point labs stat_ellipse theme_minimal geom_density
#' @importFrom dplyr filter
#' @importFrom ggsci scale_color_igv
#'
#' @export
plot_MAD_MDS <- function(cpm, metadata, filter) {

  ## Filter the 3000 most expressed genes
  gene_mads <- rowMads(cpm)
  top_3000_genes <- names(gene_mads)[order(gene_mads, decreasing = TRUE)][1:3000]
  cpm <- cpm[top_3000_genes, ]

  ## Distance matrix
  distance_matrix <- dist(t(cpm), method = "euclidean")
  mds_result <- cmdscale(distance_matrix, k = 2)

  ## MDS by SampleType
  mds_data_group <- data.frame(
    Dim1 = mds_result[, 1],
    Dim2 = mds_result[, 2],
    group = metadata$SampleType
  )

  p1 <- ggplot(mds_data_group , aes(x = Dim1, y = Dim2, color = group)) +
    geom_point(size = 4, alpha = 0.7) +
    labs(title = "MDS por grupos", x = "Dim 1", y = "Dim 2") +
    stat_ellipse() +
    theme_minimal() +
    scale_color_igv()

  ## MDS by Batch
  mds_data_batch <- data.frame(
    Dim1 = mds_result[, 1],
    Dim2 = mds_result[, 2],
    group = metadata$Batch
  )

  p2 <- ggplot(mds_data_batch , aes(x = Dim1, y = Dim2, color = group)) +
    geom_point(size = 4, alpha = 0.7) +
    labs(title = "MDS por batch", x = "Dim 1", y = "Dim 2") +
    stat_ellipse() +
    theme_minimal() +
    scale_color_igv()

  ## MAD by Sampletype
  samples <- unique(metadata$SampleType)

  mad_list <- list()

  for (group in samples) {
    group_pos <- which(metadata$SampleType == group)
    mad_values <- apply(cpm[, group_pos, drop = FALSE], 1, function(x) mad(x, center = median(x)))
    mad_list[[group]] <- mad_values
  }

  mad_data_group <- data.frame(
    mad_values = unlist(mad_list),
    group = factor(rep(names(mad_list), each = nrow(cpm)))
  )

  p3 <- ggplot(mad_data_group %>% filter(mad_values < filter), aes(x = mad_values, fill = group)) +
    geom_density(alpha = 0.5) +
    labs(title = "Densidad MAD por grupo", x = "MAD", y = "Densidad") +
    theme_minimal()

  ## MAD by Batch
  batches <- unique(metadata$Batch)

  for (b in batches) {
    pos <- grep(paste0("(_", b, ")$"), colnames(cpm))
    assign(paste0("mad_group_", b),
           apply(cpm[, pos, drop = FALSE], 1, function(x) mad(x, center = median(x))))
  }

  mad_data_batch <- data.frame(
    mad_values = unlist(lapply(batches, function(b) get(paste0("mad_group_", b)))),
    group = factor(rep(batches, each = nrow(cpm)))
  )

  p4 <- ggplot(mad_data_batch %>% filter(mad_values < filter), aes(x = mad_values, fill = group)) +
    geom_density(alpha = 0.5) +
    labs(title = "Densidad MAD por batch", x = "MAD", y = "Densidad") +
    theme_minimal()

  ## Results
  return(list(
    mds_result = mds_result,
    plot_mds_group = p1,
    plot_mds_batch = p2,
    plot_mad_group = p3,
    plot_mad_batch = p4
  ))
}
