source("simIntSiteReads_func.R")

context("Cut host DNA so R2 and R1 host DNA is the same")

host_seq <- "ACGTACGT"
oligo <- list(Primer="AAA", LTRBit="CCCC")

# molecule: AAA CCCC ACGTACGT

test_that("R2 is long  and host is small : no need to cut", {
    expect_equal(host_seq, cut_host_dna(host_seq, oligo, 100))
})

test_that("R2 could not be shorter than technical DNA", {
    expect_error(cut_host_dna(host_seq, oligo, 1))
})

test_that("last 3 nucleotides should be cut", {
    expect_equal("ACGTA", cut_host_dna(host_seq, oligo, 12))
})

test_that("last nucleotide should be cut", {
    expect_equal("ACGTACG", cut_host_dna(host_seq, oligo, 14))
})

test_that("last nucleotide should not be cut", {
    expect_equal(host_seq, cut_host_dna(host_seq, oligo, 15))
})

