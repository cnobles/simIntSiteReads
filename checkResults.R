#### check results after running intSiteCaller on simulated data ####
#' Note that reads simulated have their truth in the qnames for example
#' @M03249:1:000-SIMchr1p52699700:1:1:33:4 1:N:0:0 means
#' chr1+52699700, with length=33 down stream
#'
#' Strategy:
#' 1. read completeMetadata.RData
#' 2. load truth from Data/demultiplexedReps/GTSP0308-1_R1[2].fastq.gz
#' 3. load results from $alias/sites....
#' 4. compare alignments with different tolerance
#' 4. summarize

#' return collection of global arguments
get_args <- function() {
    suppressMessages(library(argparse))
    
    codeDir <- dirname(sub("--file=", "", grep("--file=", commandArgs(trailingOnly=FALSE), value=T)))
    if( length(codeDir)!=1 ) codeDir <- list.files(
                  path="~",
                  pattern="simIntSiteReads$",
                  recursive=TRUE, include.dirs=TRUE, full.names=TRUE)
    stopifnot(file.exists(file.path(codeDir, "simIntSiteReads.R"),
                          file.path(codeDir, "checkResults.R"),
                          file.path(codeDir, "processingParams.tsv")))
    
    parser <- ArgumentParser(formatter_class='argparse.RawTextHelpFormatter')
    parser$add_argument("-c", "--codeDir", type="character", nargs=1,
                        default=codeDir,
                        help="Directory of code")
    parser$add_argument("-p", "--workDir", type="character", nargs=1,
                        default=".",
                        help="Directory of code")
    parser$add_argument("-e", "--err", type="integer", nargs=1,
                        default=3,
                        help="tolerance for alignment")
    args <- parser$parse_args(commandArgs(trailingOnly=TRUE))
}
args <- get_args()
print(args)

libs <- c("stringr",
          "plyr",
          "dplyr",
          "GenomicRanges",
          "ShortRead")
null <- suppressMessages(sapply(libs, require, character.only=TRUE))

options(stringsAsFactors=FALSE)
options(dplyr.width = Inf)

get_metadata <- function() {
    df1 <- read.table(file.path(args$workDir, "sampleInfo.tsv"), header=TRUE)
    df2 <- read.table(file.path(args$workDir, "processingParams.tsv"), header=TRUE)
    
    df <- merge(df1, df2)
    
    ## demultiplexed R1 R2 fastq files
    df$R1fastq <- file.path(args$workDir, "Data/demultiplexedReps",
                            sprintf("%s_%s.fastq.gz", df$alias, "R1"))
    df$R2fastq <- file.path(args$workDir, "Data/demultiplexedReps",
                            sprintf("%s_%s.fastq.gz", df$alias, "R2"))
    
    ## pipeline results
    df$allSites <- file.path(args$workDir, df$alias, "allSites.RData")
    df$sites <-    file.path(args$workDir, df$alias, "sites.final.RData")
    df$rawSites <- file.path(args$workDir, df$alias, "rawSites.RData")
    df$multihit <- file.path(args$workDir, df$alias, "multihitData.RData")
    df$primerID <- file.path(args$workDir, df$alias, "primerIDData.RData")
    
    
    return(df)
}
metadata <- get_metadata()    

get_machine_file <- function(dir="Data") {
    fastqFiles <- list.files(path=dir, pattern="fastq.gz", full.names=TRUE)
    I1File <-  grep("I1", fastqFiles, value=TRUE)
    R1File <-  grep("R1", fastqFiles, value=TRUE)
    R2File <-  grep("R2", fastqFiles, value=TRUE)
}



#' Load truth from qnames in fastqfile
#' In the future will also get number of base errors
#' @param R1fastqfile
#' @return data.frame of qname, qid, chr, strand, position, width
#' @example
#' R1fastqfile="./Data/demultiplexedReps/GTSP0308-1_R1.fastq.gz"
#' load_truth_from_fastq(R1fastqfile)
load_truth_from_fastq <- function(meta=metadata) {
    stopifnot(c("R1fastq","alias") %in% colnames(meta))
    
    truth <- plyr::ldply(seq_along(nrow(meta)), function(i)
        {
            qname <- as.character(ShortRead::id(readFastq(meta$R1fastq[i])))
            info <- stringr::str_match(qname,
                                       "M03249:1:000-SIM(.*):1:1:(\\d+):(\\d+)")
            pos <- info[,2]
            pos2 <- stringr::str_match(pos, "(.*)([pm])(\\d+)")
            
            info <- cbind(info, pos2)
            
            truth <- data.frame(alias=meta[i]$alias,
                                qname =    as.character( info[,1] ),
                                qid =      as.integer(   info[,4] ),
                                chr =      as.character( info[,6] ),
                                strand =   as.character( info[,7] ),
                                position = as.integer(   info[,8] ),
                                width =    as.integer(   info[,3] ),
                                stringsAsFactors=FALSE)
            truth$strand <- sub('p', '+', truth$strand)
            truth$strand <- sub('m', '-', truth$strand)
            
            truth$breakpoint <- ifelse(truth$strand=="+",
                                       truth$position+truth$width,
                                       truth$position-truth$width)
            
            truth$qname <- sub("^M.*-", paste0(meta$alias[i],'%'), truth$qname)
            return(truth)
        } )
    
    return(truth)
}
##truth <- load_truth_from_fastq(metadata)


#' Load results from RData files, see intSiteUploader for help
#' @param meta dataframe of metadata
#' @return dataframe of siteID, multihitID, chr, strand, position, breakpoint, width, count qname, qid
#' @example
#' load_uniqueSites_from_RData()
#' load_multiSites_from_RData()
#' 
.load_uniqueSites_from_RData <- function(meta=metadata) {
    stopifnot(c("sites", "allSites") %in% colnames(meta))
    stopifnot(nrow(meta)==1)
    
    sites.final <- get(load(meta$sites))
    allSites <-    get(load(meta$allSites))
    
    sites <- plyr::ldply(1:length(sites.final), function(i)
        {
            i.gr <- allSites[unlist(sites.final[i]$revmap)]
            
            alias <- sites.final[i]$sampleName
            chr <- seqnames(sites.final[i])
            strand <- strand(sites.final[i])
            position <- start(flank(sites.final[i], -1, start=T))
            pcrBreakpoints <- start(flank(i.gr, -1, start=F))
            qname <- names(i.gr)
            qid <- as.integer(sub(".*:(\\d+)$", "\\1", qname))
            
            freq <- plyr::count(pcrBreakpoints)
            
            site <- data.frame(alias=as.character(alias),
                               siteID=i,
                               multihitID=0,
                               chr=as.character(chr),
                               strand=as.character(strand),
                               position=as.integer(position),
                               breakpoint=as.integer(pcrBreakpoints),
                               width=as.integer(abs(pcrBreakpoints-position)+1),
                               count=as.integer(1),
                               qname=as.character(qname),
                               qid=as.integer(qid),
                               stringsAsFactors=FALSE)
            
            return(site)
        } )       
    
    return(sites)
}
load_uniqueSites_from_RData <- function(meta=metadata) {
    sites <- plyr::ldply(1:nrow(meta), function(i)
        .load_uniqueSites_from_RData(meta[i,]) )
    return(sites)
}
#'
#'
#' 
.load_multiSites_from_RData <- function(meta=metadata) {
    stopifnot(c("alias", "multihit") %in% colnames(meta))
    stopifnot(nrow(meta)==1)
    multi <- get(load(meta$multihit))
    
    multiAln <- multi$unclusteredMultihits
    
    fromCluster <- findOverlaps(multiAln,
                                multi$clusteredMultihitPositions,
                                ignore.strand=FALSE,
                                maxgap=5,
                                select="first")
    
    multiAln$multihitID <- fromCluster
    
    position <- start(flank(multiAln, -1, start=TRUE))
    breakpoint <- start(flank(multiAln, -1, start=FALSE))
    qname <- names(multiAln)
    qid <- as.integer(sub(".*:(\\d+)$", "\\1", qname))
    
    msite <- data.frame(
        alias=as.character(meta$alias),
        siteID=0,
        multihitID=as.integer(multiAln$multihitID),
        chr=as.character(seqnames(multiAln)),
        strand=as.character(strand(multiAln)),
        position=as.integer(position),
        breakpoint=as.integer(breakpoint),
        width=as.integer(abs(breakpoint-position)+1),
        count=as.integer(1),
        qname=as.character(qname),
        qid=as.integer(qid),
        stringsAsFactors=FALSE)
    
    return(msite)
}
load_multiSites_from_RData <- function(meta=metadata) {
    sites <- plyr::ldply(1:nrow(meta), function(i)
        .load_multiSites_from_RData(meta[i,]) )
    return(sites)
}

truth <- load_truth_from_fastq(metadata)
res <- plyr::rbind.fill(load_uniqueSites_from_RData(),
                        load_multiSites_from_RData())

compr <- merge(truth, res, by=c("alias", "qname"), all.x=TRUE)

is_good_hit <- function(x, err=0) {
    is.good <- with(x, 
                    chr.x==chr.y &
                        strand.x==strand.y &
                            abs(position.x-position.y)<=err &
                                abs(breakpoint.x-breakpoint.y)<=err )
    is.good[is.na(is.good)] <- FALSE
    return(is.good)
}

compr$hit <- is_good_hit(compr, args$err)

compr$hit <- is_good_hit(compr, 3)
compr$hit <- is_good_hit(compr, 2)
compr$hit <- is_good_hit(compr, 1)
compr$hit <- is_good_hit(compr, 0)

compr <- (compr %>%
              group_by(qid.x) %>%
                  mutate(nAln=n(),
                         anyHit=any(hit)) )

sum(compr$hit)

sum(readRecovered$readAligned)


