#' introduce sequencing error in reads
#' @param miseq_reads dataframe with 4 columns: R1, R2, I1, qname character vectors
#' @param seq_error sequencing error specification
#' @return dataframe with the same columns as input df but with errors
generate_seq_error <- function(miseq_reads, seq_error) {
    data.frame(
        I1=miseq_reads$I1,
        R1=uniform_seq_error(miseq_reads$R1, seq_error),
        R2=uniform_seq_error(miseq_reads$R2, seq_error),
        qname=miseq_reads$qname
    )
}

#' introduce substitution error
#' @param sequences character vector of error-free seq
#' @param seq_error percentage of nucleotide to change
#' @return character vector of sequences with errors
uniform_seq_error <- function(sequences, seq_error) {
    stopifnot(seq_error >= 0 & seq_error <= 100)
    sequences
}
