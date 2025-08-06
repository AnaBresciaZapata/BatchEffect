#' Gráficos MDS y MAD para análisis de variabilidad
#'
#' Esta función genera gráficos para analizar la variabilidad de datos de expresión normalizados en CPM.
#' Calcula un análisis de escalamiento multidimensional (MDS) para visualizar la similitud entre muestras
#' según el tipo de muestra y batch. También calcula la dispersión mediana absoluta (MAD) por grupo y batch.
#'
#' @param cpm Matrix. Matriz de expresión en CPM con genes en filas y muestras en columnas.
#' @param metadata Data frame. Debe tener las columnas `SampleType` y `Batch`, y las filas deben corresponder
#'        a las muestras (identificadores iguales a las columnas de la matriz cpm).
#' @param filter Numeric. Umbral numérico para filtrar valores extremos de MAD en los gráficos de densidad.
#'
#' @return Una lista con:
#' \describe{
#'   \item{mds_result}{Matriz con las dos primeras dimensiones del análisis MDS}
#'   \item{plot_mds_group}{Gráfico MDS por tipo de muestra}
#'   \item{plot_mds_batch}{Gráfico MDS por batch}
#'   \item{plot_mad_group}{Distribución MAD por tipo de muestra}
#'   \item{plot_mad_batch}{Distribución MAD por batch}
#' }
#'
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
