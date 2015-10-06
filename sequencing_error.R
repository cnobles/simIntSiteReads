library(stringr)

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
    n_errors <- get_number_errors(sequences, seq_error)
    mapply(change_nucleotides, sequences, n_errors, USE.NAMES=FALSE)
}

#' how many errors need to be introduced
#' @param sequences character vector of error-free seq
#' @param seq_error percentage of nucleotide to change
#' @return integer vector of number of erros to introduce for each seq
get_number_errors <- function(sequences, seq_error) {
    n_sequences <- length(sequences)
    seq_len <- sapply(sequences, nchar, USE.NAMES=FALSE)
    n_errors <- (seq_len*seq_error)/100.00 # we can have 1.4 bases to modify
    n_errors_int <- floor(n_errors) # that takes care of int part
    n_errors_reminder <- n_errors -  n_errors_int
    add_one <- runif(n_sequences) < n_errors_reminder # that takes care of reminder
    n_errors_int + add_one # average is close to expected error rate
}

#' substitute nucleotides
#' @param sequence characters of DNA representation
#' @param n_error number of nucleotide to change
change_nucleotides <- function(sequence, n_error) {
    positions <- sample(1:nchar(sequence), n_error)
    sequence_errors <- sequence
    null <- sapply(positions, function(pos) {
        sequence_errors <<- change_nucleotide_at(sequence_errors, pos)
    })
    sequence_errors
}

UNAMBIGUOUS_DNA_ALPHABET <- c('A', 'G', 'T', 'C')

#' change one nucleotide in sequence
#' @param sequence characters of DNA representation
#' @param position at what position to introduce error
change_nucleotide_at <- function(sequence, position) {
    stopifnot(position <= nchar(sequence))
    prefix <- str_sub(sequence, start=1, end=position-1)
    nucleotide_to_change <- str_sub(sequence, start=position, end=position)
    suffix <- str_sub(sequence, start=position+1)

    neighbours <- setdiff(UNAMBIGUOUS_DNA_ALPHABET, nucleotide_to_change)
    erroneous_nucleotide <- sample(neighbours, 1)

    str_c(prefix, erroneous_nucleotide, suffix)
}
