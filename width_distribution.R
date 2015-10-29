#' generate reads from uniform distribution(possibly with duplicates)
#' 
#' @num_reads number of reads
#' @min_value minimum molecule length
#' @max_value maximum molecule length
#' @return integer vector with width values
uniform_width_distribution <- function(num_reads, min_value, max_value) {
    sample(seq(min_value, max_value), num_reads, replace=TRUE)
}

#' generate reads with Maxwell-Boltzmann distribution(emulate real 
#' gaussian-like but assymetrical molecule length distrib)
#' see: https://en.wikipedia.org/wiki/Maxwell%E2%80%93Boltzmann_distribution
#' for mean and PDF formulas (note no need for normalization)
#' 
#' @num_reads number of reads
#' @mean_value mean molecule length
#' @return integer vector with width values
maxwell_boltzmann_width_distribution <- function(num_reads, mean_value) {
    a <- sqrt(pi/2.0)*mean_value/2.0
    mb_pdf <- function(x) x**2 * exp(-x**2/(2*a**2))
    mb_dist <- AbscontDistribution(d=mb_pdf, low1=0)    
    rdist <- r(mb_dist)
    as.integer(rdist(num_reads))
}
