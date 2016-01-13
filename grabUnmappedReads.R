libs <- c("stringr",
          "dplyr",
          "ggplot2",
          "GenomicRanges",
          "ShortRead",
          "BiocParallel",
          "BSgenome")
null <- suppressMessages(sapply(libs, library, character.only=TRUE))


options(stringsAsFactors=FALSE)
options(dplyr.width = Inf)
#' increase output width to console width
wideScreen <- function(howWide=as.numeric(strsplit(system('stty size', intern=T), ' ')[[1]])[2]) {
   options(width=as.integer(howWide))
}
wideScreen()

site.res <- read.csv("./truth.site.call.txt", header=TRUE, sep="\t")

site.res.nf <- dplyr::filter(site.res,
                             is.na(subjectHits.uniq) &
                             is.na(subjectHits.multi))

keyFile <- list.files(".", pattern="keys.RData", recursive=TRUE)

lostSiteName <- with(site.res.nf,
                     paste(chr,
                           ifelse(strand=="+", "p", "m"),
                           position, sep="")
                     )

lapply(keyFile, function(f) {
    f=keyFile[60] # 60, 61, 62
    keys <- get(load(f))
    keys$sitename <- str_match(keys$names, "chr.*[pm]\\d+")[,1]
    keys.sitename <- str_match(keys$names, "chr.*[pm]\\d+")
    if( !any(keys.sitename %in% lostSiteName) ) return(NULL)
    return(TRUE)
    
    key.lostSite <- unique(keys.sitename[keys.sitename %in% lostSiteName])
    
    load(sub("keys.RData", "hits.R1.RData", f))
    load(sub("keys.RData", "hits.R2.RData", f))
    
    hits.R1$siteName <- str_match(names(hits.R1), "chr.*[pm]\\d+")[,1]
    hits.R2$siteName <- str_match(names(hits.R2), "chr.*[pm]\\d+")[,1]
    
    hits.R1 <- subset(hits.R1, siteName %in% key.lostSite)
    strand(hits.R1) <- c("+","-")[3-as.integer((factor(as.factor(strand(hits.R1)), levels=c("+","-"))))]
    hits.R2 <- subset(hits.R2, siteName %in% key.lostSite)
    hits <- c(hits.R1, hits.R2)
    hits$readID <- str_match(names(hits), "\\d+$")[,1]
    
    hits <- split(hits, hits$readID)
    
    fa.R1 <- readFasta(sub("keys.RData", "R1-1.fa", f))
    fa.R2 <- readFasta(sub("keys.RData", "R2-1.fa", f))
    
    lostSite.keys <- subset(keys, sitename %in% lostSiteName)
    lostSite.keys$fa1 <- as.character(sread(fa.R1)[lostSite.keys$R1])
    lostSite.keys$fa2 <- as.character(sread(fa.R2)[lostSite.keys$R2])
    
    read.1 <- ShortRead(DNAStringSet( lostSite.keys$fa1 ),
                        BStringSet( lostSite.keys$names ) )
    read.2 <- ShortRead(DNAStringSet( lostSite.keys$fa2 ),
                        BStringSet( lostSite.keys$names ) )
    
    writeFasta(read.1, file="lostSite.R1.fa")
    writeFasta(read.2, file="lostSite.R2.fa")
    
} )


R1.gr <- subset(hits[[1]], from=="R1")
names(R1.gr) <- NULL
R2.gr <- subset(hits[[1]], from=="R2")
names(R2.gr) <- NULL

R1R2 <- merge(as.data.frame(R1.gr),
              as.data.frame(R2.gr),
              by=c("siteName", "readID", "seqnames", "strand"),
              suffixes = c(".1",".2"))

R1R2$tBaseInsert.1 <- NULL
R1R2$from.1 <- NULL
R1R2$tBaseInsert.2 <- NULL
R1R2$from.2 <- NULL

R1R2$position <- with(R1R2, ifelse(strand=="+", start.2, end.2))
R1R2$breakpoint <- with(R1R2, ifelse(strand=="+", end.1, start.1))


subset(R1R2,
       seqnames=="chr8" &
       abs(start.2 - 42377425) < 5000 )

tmp <- subset(R1R2,
       seqnames=="chr8" &
       abs(start.1 - 42377425) < 5000 )





options("showHeadLines"=7)
options("showTailLines"=2)

 ## Revert to default values
options("showHeadLines"=NULL)
options("showTailLines"=NULL)

options("showTailLines"=Inf)
options("showHeadLines"=Inf)

