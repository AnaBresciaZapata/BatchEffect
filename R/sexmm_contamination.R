#' Detect sex mismatch and contamination in RNA-seq samples
#'
#' This function evaluates potential sex mismatches and contamination in RNA-seq datasets
#' based on expression of the \code{XIST} gene and Y chromosome genes. It produces
#' a scatter plot of normalized expression values, classifies samples, and updates a
#' filter table for downstream QC.
#'
#' @param counts A numeric matrix or data frame with genes in rows and samples in columns.
#'        Row names should be ENSEMBL gene IDs. The column names should be sample IDs.
#' @param metadata A data frame with at least two columns:
#'        \itemize{
#'          \item{\code{ID}:} Sample identifiers matching \code{colnames(counts)}
#'          \item{\code{Sex}:} Biological sex ("Male" or "Female")
#'        }
#' @param chrY_genes A character vector of ENSEMBL gene IDs located on chromosome Y.
#' @param contamination_slope Either `"estimated"` or a numeric value in degrees (e.g. 45).
#' @param contamination_fraction Numeric value (recommended 0.3) controlling the angular range
#'        of the contamination area.
#' @param plot_title Character. Title to display in the plot.
#'
#' @return A list with:
#' \describe{
#'   \item{plot}{A \code{ggplot2} object with the contamination/mismatch plot.}
#'   \item{status_table}{The status table with columns for contamination
#'                   and mismatch status.}
#'   \item{genes_cpm}{The CPM-normalized expression matrix.}
#' }
#'
#' @import edgeR
#' @import dplyr
#' @import ggplot2
#' @export
#'
#' @examples
#' \dontrun{
#'   result <- sexmm_contamination(
#'     counts = counts_matrix,
#'     metadata = metadata_df,
#'     chrY_genes = chrY_vector,
#'     contamination_slope = "estimated",
#'     contamination_fraction = 0.3,
#'     plot_title = "Dataset1"
#'   )
#'   print(result$plot)
#'   head(result$status_table)
#' }


Sexmm_contamination <- function(counts, metadata, chrY_genes_vector, contamination_slope, contamination_fraction, plot_title){
  if(length(grep("[.]", rownames(counts)))>0) {rownames(counts) <- sapply(strsplit(rownames(counts),"[.]"), `[`, 1) }

  #calculate cpm

  counts_cpm <- cpm(round(counts))
  rownames(counts_cpm) <- rownames(counts)
  XIST_exp <- counts_cpm["ENSG00000229807",]

  chrY_genes_exp <- colSums(counts_cpm[rownames(counts_cpm) %in% chrY_genes_vector,])
  print(max(chrY_genes_exp))
  metadata$XIST_exp <- XIST_exp[match(rownames(metadata),names(XIST_exp))]
  metadata$Y_exp <- chrY_genes_exp[match(rownames(metadata),names(chrY_genes_exp))]

  #Remove NAs
  metadata <- metadata[!is.na(metadata$XIST_exp) | !is.na(metadata$Y_exp),]

  metadata$XIST_norm <- metadata$XIST_exp/max(metadata$XIST_exp)
  metadata$Y_norm <- metadata$Y_exp/max(metadata$Y_exp)

  #preparing for the plot
  contamination_area = contamination_fraction*90

  max_norm <- max(metadata$XIST_norm, metadata$Y_norm,na.rm = TRUE)
  metadata$Sex1 <- ifelse(metadata$Sex=='Male',1,2)


  if(contamination_slope == "estimated"){
    x_estimated_median <- median(metadata[metadata$Sex1 == 1, "XIST_norm"])
    y_estimated_median <- median(metadata[metadata$Sex1 == 2, "Y_norm"])

    hip <- sqrt(x_estimated_median^2 + y_estimated_median^2)

    middle_slope <- asin(y_estimated_median/hip)
    our_contamination_slope <- atan(middle_slope) * (180 / pi)

    lower_slope <- tan((our_contamination_slope - contamination_area / 2) / 180*pi)
    upper_slope <- tan((our_contamination_slope + contamination_area / 2) / 180*pi)

  }
  else if(!is.na(as.numeric(contamination_slope))){
    contamination_slope <- as.numeric(contamination_slope)
    lower_slope <- tan((contamination_slope - contamination_area / 2) / 180*pi)
    upper_slope <- tan((contamination_slope + contamination_area / 2) / 180*pi)
    middle_slope <- tan(contamination_slope / 180*pi)
  }

  metadata$expressionSexNaive <- case_when(
    metadata$Y_norm > metadata$XIST_norm * middle_slope ~ 1,
    metadata$Y_norm < metadata$XIST_norm * middle_slope ~ 2
  )
  x_expression_median <- median(metadata[metadata$Sex1 == 1 & metadata$expressionSexNaive == 1, "XIST_norm"])
  y_expression_median <- median(metadata[metadata$Sex1 == 2 & metadata$expressionSexNaive == 2, "Y_norm"])

  x_expression_min <- min(metadata[metadata$Sex1 == 2 & metadata$expressionSexNaive == 2, "XIST_norm"])
  y_expression_min <- min(metadata[metadata$Sex1 == 1 & metadata$expressionSexNaive == 1, "Y_norm"])

  y_max <- max(metadata$Y_norm)
  x_max <- max(metadata$XIST_norm)

  x_anchor <- min(x_expression_median, x_expression_min - 0.01 * x_max)
  y_anchor <- min(y_expression_median, y_expression_min - 0.01 * y_max)

  metadata$XIST_corrected <- metadata$XIST_norm - x_anchor

  metadata$contaminated <- case_when(
    (metadata$Y_norm > ((metadata$XIST_corrected) * lower_slope + y_anchor)
     & metadata$Y_norm < ((metadata$XIST_corrected) * upper_slope + y_anchor)) ~ "yes",
    TRUE ~ "no"
  )

  metadata$expressionSex <- case_when(
    (metadata$Y_norm < ((metadata$XIST_corrected) * middle_slope + y_anchor)) ~ 2,
    (metadata$Y_norm > ((metadata$XIST_corrected) * middle_slope + y_anchor)) ~ 1
  )

  metadata$mismatch <- case_when(
    metadata$Sex1 == 0 ~ "unknown",
    metadata$expressionSex == metadata$Sex1 ~ "no",
    metadata$expressionSex != metadata$Sex1 ~ "yes"
  )

  metadata$status <- case_when(
    metadata$contaminated == "yes" & metadata$mismatch == "yes" ~ "Contaminated and sex mismatch",
    metadata$contaminated == "yes" ~ "Likely contaminated",
    metadata$mismatch == "yes" ~ "Sex mismatch",
    TRUE ~ "Passed"
  )
  status_table <- metadata[,c("ID","contaminated", "mismatch","status")]
  y_genes_zoom <- max(metadata %>% ungroup() %>%
                        filter(contaminated == "yes" | mismatch == "yes") %>%
                        summarise(max_y = max(Y_norm), max_x = max(XIST_norm)))

  exclusion_zone <- tibble(x = c(x_expression_median, max_norm)) %>%
    mutate(lower_bound = (x - x_expression_median) * lower_slope + y_expression_median,
           upper_bound = (x - x_expression_median) * upper_slope + y_expression_median,
           middle_line = (x - x_expression_median) * middle_slope + y_expression_median)

  metadata$Sex1 <- as.factor(metadata$Sex1)

  sexmm_plot <- ggplot(data = exclusion_zone, aes(x = x, ymin = lower_bound, ymax = upper_bound)) +
    geom_ribbon(alpha = 0.2) +
    geom_line(aes(x = x, y = middle_line), linetype = 2, colour = "blue") +
    geom_point(data = metadata, inherit.aes = F, aes(col = status, shape = Sex1, x = XIST_norm, y = Y_norm)) +
    scale_colour_manual(
      values = alpha(c("Passed" = "black",
                       "Likely contaminated" = "orange",
                       "Sex mismatch" = "red",
                       "Contaminated and\nsex mismatch" = "darkmagenta"),
                     0.5),
      name = "Passed checks") +
    coord_cartesian(ylim = c(0, max_norm), xlim = c(0, max_norm)) +
    theme_bw(base_size = 12) +
    xlab("XIST normalized expression")+
    ylab("Chr Y genes normalized expression")+
    ggtitle(plot_title)+
    scale_shape_manual(values = c(17,16),labels=c("Male", "Female"))
  print(sexmm_plot)

  return(list(plot=sexmm_plot,status_table=status_table, genes_cpm=counts_cpm))
}

