library(stringdist)

source("sequencing_error.R")

context("Uniform random error")

seqs <- c("AAAAAAAAAA", "TTTTTTTTTT", "GATCGATCGA")
seqs_len <- sapply(seqs, nchar)

test_that("seq_error is valid", {
    expect_error(uniform_seq_error(seqs, seq_error=120))
    expect_error(uniform_seq_error(seqs, seq_error=-120))
})

test_that("seq_error of 0 % is valid and don't change sequences", {
    expect_equal(uniform_seq_error(seqs, seq_error=0), seqs)
})

test_that("seq_error of 100 % is valid and completely change sequences", {
    seq_error <- uniform_seq_error(seqs, seq_error=100)
    n_errors <- stringdist(seqs, seq_error, method="hamming")
    expect_equal(sum(n_errors), sum(seqs_len))
})

test_that("seq of 10 bases and 10% error lead to 1 error in each read", {
    seq_error <- uniform_seq_error(seqs, seq_error=10)
    n_errors <- stringdist(seqs, seq_error, method="hamming")
    expect_equal(n_errors, rep(1, length(seqs)))
})

large_seqs <- c(seqs, seqs, seqs, "AAAAAAAAAA") # 10 sequences total

test_that("seq of 10 bases and 15% error: number of errors differ ", {
    seq_error <- uniform_seq_error(seqs, seq_error=15)
    # so should have around 15 errors
    n_errors <- stringdist(seqs, seq_error, method="hamming")
    expect_true(all(n_errors >= 1), 'each read have at least one error')
    expect_true(any(n_errors >= 2), 'some reads have 2 errors')
    n_errors_total <- sum(n_errors)
    # 15 +- 3
    expect_true(n_errors_total > 12)
    expect_true(n_errors_total < 18)
})
