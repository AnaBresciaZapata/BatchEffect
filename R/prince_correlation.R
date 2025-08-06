#' Correlation heatmap between metadata variables and principal components
#'
#' This function computes the association between experimental variables (metadata) and
#' the principal components of expression data using the `prince` method. It plots a heatmap
#' that shows the strength of correlation between each variable and the first 9 PCs.
#'
#' @param data A numeric matrix with genes in rows and samples in columns.
#' @param metadata Data frame. Must contain the columns `SampleType` and `Batch`, and row names
#'        must correspond to the samples (identifiers matching the columns of the data).
#'
#' @return A list with two elements:
#' \describe{
#'   \item{res}{The output of the `prince` function containing correlation results.}
#'   \item{plot}{The heatmap plot showing associations between metadata variables and PCs.}
#' }
#'
#' @importFrom grDevices heat.colors
#' @export
prince_correlation <- function(data, metadata, pcs) {
  res <- prince(g = data, o = metadata, top = pcs, imputeknn = FALSE, center = TRUE, permute = FALSE)

  plot <- prince.plot(res, label = colnames(res$o), smallest = -20, note = FALSE,
                      notecol = "black", notecex = 1,
                      breaks = seq(-20, 0, length.out = 100), col = heat.colors(99),
                      margins = c(5, 7), key = TRUE, cexRow = 1, cexCol = 1,
                      xlab = "Principal Components (Variation)", colsep = NULL,
                      rowsep = NULL, sepcolor = "black", sepwidth = c(0.05, 0.05),
                      Rsquared = FALSE, breaksRsquared = seq(0, 1, length.out = 100))

  print(plot)

  return(list(res = res, plot = plot))
}
