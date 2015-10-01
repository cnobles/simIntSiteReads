libs <- c("testthat",
          "GenomicRanges",
          "ShortRead",
          "BSgenome",
          sprintf("BSgenome.Hsapiens.UCSC.%s", "hg18"))
null <- suppressMessages(sapply(libs, library, character.only=TRUE))

source("simIntSiteReads_func.R")

## test on + strand
context('test fetching sequence downstream')
pos <- 59165408
for( chr in c("chr1") ) {
    for( pos in c(59165408, 49165408, 39165408) ) {
        context(paste('Position', chr, pos))
        context('Sequence on + strand')
        width <- c(30, 40, 50)
        expect_equal(
            get_sequence_downstream(Hsapiens, chr, pos, "+", width)$seq,
            sapply(width, function(w) as.character(Hsapiens[[chr]][1:w-1+pos])) )
        
        
        context('Sequence on - strand')
        width <- c(30, 40, 50)
        expect_equal(
            get_sequence_downstream(Hsapiens, chr, pos, "-", width)$seq,
            sapply(width, function(w) as.character(reverseComplement(Hsapiens[[chr]][1:w+pos-w]) ) ) )
        
    }
}


