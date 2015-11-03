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
                        default=5,
                        help="tolerance for alignment")
    args <- parser$parse_args(commandArgs(trailingOnly=TRUE))
}
args <- get_args()
print(args)

libs <- c("stringr",
          "plyr",
          "dplyr",
          "RMySQL",
          "GenomicRanges",
          "ShortRead")
null <- suppressMessages(sapply(libs, require, character.only=TRUE))

options(stringsAsFactors=FALSE)
options(dplyr.width = Inf)
#' increase output width to console width
wideScreen <- function(howWide=as.numeric(strsplit(system('stty size', intern=T), ' ')[[1]])[2]) {
   options(width=as.integer(howWide))
}
wideScreen()

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
##metadata <- get_metadata()    

#' get the original machine undemultiplexed I1, R1, R2 fastq.gz files
get_machine_file <- function(dir="Data") {
    fastqFiles <- list.files(path=dir, pattern="fastq.gz", full.names=TRUE)
    I1File <-  grep("I1", fastqFiles, value=TRUE)
    R1File <-  grep("R1", fastqFiles, value=TRUE)
    R2File <-  grep("R2", fastqFiles, value=TRUE)
    stopifnot(length(I1File)==1)
    stopifnot(length(R1File)==1)
    stopifnot(length(R2File)==1)
    return(data.frame(uI1=I1File, uR1=R1File, uR2=R2File) )
}
##get_machine_file(dir=file.path(args$workDir,"Data"))


#' Load truth from qnames in fastqfile
#' In the future will also get number of base errors
#' @param R1fastqfile
#' @return data.frame of qname, qid, chr, strand, position, width
#' @example
#' R1fastqfile="./Data/demultiplexedReps/GTSP0308-1_R1.fastq.gz"
#' load_truth_from_fastq(R1fastqfile)
load_truth_from_fastq <- function(meta=metadata) {
    stopifnot(c("R1fastq","alias") %in% colnames(meta))
    
    ##truth <- plyr::ldply(seq_along(nrow(meta)), function(i)
    truth <- plyr::ldply(seq_along(nrow(meta))[1], function(i)
        {
            ##qname <- as.character(ShortRead::id(readFastq(meta$R1fastq[i])))
            message("Reading ", meta$uI1[i])
            qname <- as.character(ShortRead::id(readFastq(meta$uI1[i])))
            info <- stringr::str_match(qname,
                                       "M03249:1:000-SIM(.*):1:1:(\\d+):(\\d+)")
            pos <- info[,2]
            pos2 <- stringr::str_match(pos, "(.*)([pm])(\\d+)")
            
            info <- cbind(info, pos2)
            
            ##truth <- data.frame(alias=meta$alias[i],
            truth <- data.frame(
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
            
            ##truth$qname <- sub("^M.*-", paste0(meta$alias[i],'%'), truth$qname)
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
    
    message("Loading\t", meta$sites, "\t", meta$allSites) 
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
    message("Loading\t", meta$multihit)
    multi <- get(load(meta$multihit))
    
    multiAln <- multi$unclusteredMultihits
    
    fromCluster <- findOverlaps(multiAln,
                                multi$clusteredMultihitPositions,
                                ignore.strand=FALSE,
                                maxgap=args$err,
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
    ##sites <- plyr::ldply(1:nrow(meta), function(i)
    ##    .load_multiSites_from_RData(meta[i,]) )
    ##return(sites)
    sites <- plyr::ldply(1:nrow(meta), function(i)
        {
            i.df <- try(.load_multiSites_from_RData(meta[i,]))
            if( class(i.df) == "try-error" ) i.df <- data.frame()
            return(i.df)
        } )
    return(sites)
}

metadata <- cbind(get_metadata(),
                  get_machine_file(dir=file.path(args$workDir,"Data")))

truth <- load_truth_from_fastq(metadata)

truth.site <- truth %>%
              dplyr::select(chr, strand, position) %>%
              dplyr::distinct()
truth.site.gr <-  makeGRangesFromDataFrame(truth.site,
                                           start.field="position",
                                           end="position",
                                           strand.field="strand")
    
res.uniq <- load_uniqueSites_from_RData()
res.uniq.site <- res.uniq %>%
                 dplyr::select(chr, strand, position, siteID, multihitID) %>%
                 dplyr::distinct()
res.uniq.site.gr <- makeGRangesFromDataFrame(res.uniq.site,
                                             start.field="position",
                                             end="position",
                                             strand.field="strand",
                                             keep.extra.columns=TRUE)

uniq.site.ovl <- findOverlaps(truth.site.gr, res.uniq.site.gr, maxgap=args$err)

res.multi <- load_multiSites_from_RData()
res.multi.site <- res.multi %>%
                  dplyr::select(chr, strand, position, siteID, multihitID) %>%
                  dplyr::distinct()
res.multi.site.gr <- makeGRangesFromDataFrame(res.multi.site,
                                              start.field="position",
                                              end="position",
                                              strand.field="strand",
                                              keep.extra.columns=TRUE)
multi.site.ovl <- findOverlaps(truth.site.gr, res.multi.site.gr, maxgap=args$err)

all.site.ovl <- merge(as.data.frame(uniq.site.ovl), as.data.frame(multi.site.ovl),
                      by="queryHits",
                      suffixes=c(".uniq", ".multi"),
                      all.x=TRUE, all.y=TRUE)

call.site.stat <-c(
    "all.sim"= nrow(truth.site),
    "found.any"=length(unique(all.site.ovl$queryHits)),
    "found.uniq"=length(unique(subset(all.site.ovl, !is.na(subjectHits.uniq))$queryHits)),
    "found.uniq.only"=length(unique(subset(all.site.ovl, !is.na(subjectHits.uniq) &  is.na(subjectHits.multi))$queryHits)),
    "found.multi"=length(unique(subset(all.site.ovl, !is.na(subjectHits.multi))$queryHits)),
    "found.multi.only"=length(unique(subset(all.site.ovl, !is.na(subjectHits.multi) & is.na(subjectHits.uniq)))$queryHits),
    "found.both"=length(unique(subset(all.site.ovl, !is.na(subjectHits.multi) & !is.na(subjectHits.uniq))$queryHits))
    )


print(as.data.frame(call.site.stat))
message("\nCaller site stat written callstat.txt" )
write.table(as.data.frame(call.site.stat), file="callstat.txt",
            quote=FALSE, sep="\t", col.name=FALSE)

truth.site.call <- merge(cbind(truth.site, queryHits=1:nrow(truth.site)),
                         all.site.ovl,
                         by="queryHits", all.x=TRUE)
message("\nTruth site called written truth.site.call.txt")
message("awk '$NF~/NA/ && $(NF-1)~/NA/' truth.site.call.txt ## not recovered sites")
write.table(truth.site.call, file="truth.site.call.txt",
            quote=FALSE, sep="\t", col.name=TRUE, row.name=FALSE)

q()

get_repeakMasker <- function(freeze="hg18") {
    
    connect2UCSC <- function(freeze="hg18") {
        junk <- sapply(dbListConnections(MySQL()), dbDisconnect)
        dbConn <- dbConnect(MySQL(),
                            host="genome-mysql.cse.ucsc.edu",
                            dbname=freeze, 
                            username="genome")
        return(dbConn)    
    }
    dbConn <- connect2UCSC()
    
    get_repeakMasker_tab <- function(conn=dbConn) {
        stopifnot(dbGetQuery(conn, "select 1")==1)
        freeze.tab <- dbGetQuery(conn, "SHOW TABLES")
        rmsk <- grep("chr[0-9XY]*_rmsk$", freeze.tab[,1], value=TRUE)
        return(rmsk)
    }
    
    rmsk.tab <- get_repeakMasker_tab(dbConn)
    
    df <- plyr::ldply(rmsk.tab, function(tab)
        {
            sql <- sprintf("SELECT genoName, genoStart, genoEnd, genoLeft, strand FROM %s", tab)
            message(sql)
            suppressWarnings(dbGetQuery(dbConn, sql))
        } )
    
    colnames(df) <- c("chr", "start", "end", "left", "strand")
    return(df)
}
rmsk <- get_repeakMasker()
rmsk.gr <- makeGRangesFromDataFrame(rmsk,
                                    start.field="start",
                                    end="end",
                                    strand.field="strand",
                                    keep.extra.columns=TRUE)


truth.site.call <- read.table("truth.site.call.txt", header=TRUE)
truth.site.call.gr <- makeGRangesFromDataFrame(truth.site.call,
                                               start.field="position",
                                               end="position",
                                               strand.field="strand",
                                               keep.extra.columns=TRUE)

findOverlaps(truth.site.call.gr, rmsk.gr, maxgap=args$err)

findOverlaps(subset(truth.site.call.gr, is.na(subjectHits.uniq) & is.na(subjectHits.multi)),
             subset(rmsk.gr, width>100),
             maxgap=args$err,
             ignore.strand=TRUE)





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


