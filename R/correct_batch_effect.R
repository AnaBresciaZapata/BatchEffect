#' Batch effect correction in gene expression data using technical replicates
#'
#' This function performs batch effect correction on count data by using clinical replicates across different batches.
#' The correction is based on the relative dispersion between replicates, adjusting the counts accordingly by gene.
#'
#' @param ref_batch_n Integer. Batch number used as the reference batch for normalization.
#' @param counts Matrix. Raw count data rounded to integers with samples in columns and genes in rows.
#' @param metadata Data frame. With two columns: \code{SampleType} (sample classification: control, disease, replicate, etc)
#'        and \code{Batch} (batch number). Rownames must correspond to the samples and match with the colnames from \code{counts}.
#' @param replicates Character vector. Labels in \code{SampleType} to indicate replicate samples.
#'
#' @return A list with:
#' \describe{
#'   \item{counts_corrected}{Matrix count corrected}
#'   \item{adjustment_report}{Correction factors applied to each gene per batch}
#'   \item{cpm}{Original matrix in CPM}
#'   \item{cpm_corrected}{Corrected matrix in CPM}
#' }
#'
#' @importFrom edgeR cpm
#'
#' @export
correct_batch_effect <- function(ref_batch_n, counts, metadata, replicates) {

  # cpm
  cpm <- cpm(counts)

  # List of all the batches
  batches_list <- unique(metadata$Batch)

  # Take out the reference batch
  batches_list <- batches_list[batches_list != ref_batch_n]

  # Reference matrix
  indv_ref <- colnames(cpm)[colnames(cpm) %in% rownames(metadata)[metadata$Batch == ref_batch_n]]
  batch_ref <- cpm[,indv_ref]

  # Final matrix
  batch_corrected <- cpm

  # Report
  adjustment_report <- list()

  # Batch by batch
  for (batch in batches_list) {

    # New matrix
    indv <- colnames(cpm)[colnames(cpm) %in% rownames(metadata)[metadata$Batch == batch]]
    new_batch <- cpm[,indv]

    # Replicates samples
    replicates_ref <- rownames(metadata)[metadata$SampleType %in% replicates & metadata$Batch == ref_batch_n]
    replicates_new <- rownames(metadata)[metadata$SampleType %in% replicates & metadata$Batch == batch]

    overlap_sample_list <- list("ref" = c(replicates_ref),
                                "new" = c(replicates_new))

    # Correction of the new batch
    new_batch_corrected <- new_batch
    adjustment_per_gene <- numeric(length = nrow(batch_ref))
    names(adjustment_per_gene) <- rownames(batch_ref)

    # By gene
    for (gen in rownames(batch_ref)) {

      # Get replicates samples for each matrix
      samples_ref <- overlap_sample_list$ref
      samples_new <- overlap_sample_list$new

      # Extract replicate values for the gene
      bridge_ref <- batch_ref[gen, samples_ref]
      bridge_new <- new_batch[gen, samples_new]

      # Relative dispersion
      diffs <- as.numeric((bridge_new - bridge_ref) / bridge_ref)

      # Median of the differences
      adjustment <- median(diffs, na.rm = TRUE)

      # Correct new batch for the gene
      new_batch_corrected[gen, ] <- new_batch[gen, ] / (1 + adjustment)
      adjustment_per_gene[gen] <- adjustment
    }

    # Save the corrected batch
    batch_corrected[, indv] <- new_batch_corrected

    # Save the adjustment value
    adjustment_report[[paste0("Batch_", batch)]] <- adjustment_per_gene
  }

  # From cpms to counts
  total_counts_vector <- colSums(counts)
  counts_corrected <- sweep(batch_corrected, 2, total_counts_vector, FUN = function(cpm, total) {
    (cpm * total) / 1e6
  })

  # Report
  adjustment_df <- do.call(cbind, adjustment_report)

  return(list("counts_corrected" = counts_corrected,
              "adjustment_report" = adjustment_df,
              "cpm" = cpm,
              "cpm_corrected" = batch_corrected))
}


