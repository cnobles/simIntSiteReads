df <- data.frame(chr=site$chr,
                 start=site$position,
                 strand=site$strand,
                 width=width)

merge(data.frame(chr=site$chr,
                 start=site$position,
                 strand=site$strand),
      data.frame(width=width))



#' get random loci from reference genome
#' @param sp name of genome, Hsapiens, etc
#' @param n  number of loci
#' @return data.frame of chr, position, strand
#' @example
#' 
get_random_loci <- function(sp=Hsapiens, n=20) {
    chrLen <- seqlengths(sp)
    chrLen <- chrLen[grepl("^chr\\d+$", names(chrLen))]
    
    n.chr <- round(as.numeric(chrLen)/sum(as.numeric(chrLen))*n)
    n.chr[1] <- n.chr[1]+n-sum(n.chr)
    n.chr <- setNames(n.chr, names(chrLen))
    
    loci <- sapply(setNames(seq_along(n.chr), names(chrLen)),
                   function(i) {
                       return(as.integer(runif(n.chr[i],
                                               min = 1000,
                                               max = chrLen[i]-3000)))
                   })
    
    loci.df <- plyr::ldply(loci, function(L) data.frame(position=L), .id="chr")
    return(loci.df)
}



scorez <- paste(rep("z", max(nchar(df[[pair]]))), collapse="")
score <- substring(scorez, 1, nchar(df[[pair]])) 

score <- sapply(nchar( df[[pair]] ), function(i) paste(rep("z",i), collapse=""))
                                 
