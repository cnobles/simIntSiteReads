library(stringdist)
library(stringr)

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
large_seqs <- rep(large_seqs, 10) # finally 100 sequences

test_that("100 seqs of 10 bases and 15% error: number of errors differ between seqs", {
    seq_error <- uniform_seq_error(large_seqs, seq_error=15)
    n_errors <- stringdist(large_seqs, seq_error, method="hamming")
    expect_true(all(n_errors >= 1), 'each read have at least one error')
    expect_true(any(n_errors >= 2), 'some reads have 2 errors')
    expect_true( ! all(n_errors >= 3), 'none reads have 3 errors')
    n_errors_total <- sum(n_errors)
    # so should have around 100 sequences * 10 basepairs * 15 % error = 150 errors
    # 150 +- 15
    expect_true(n_errors_total > 135)
    expect_true(n_errors_total < 165)
})
