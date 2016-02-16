####################################################################
#' get UCSC repeatMasker
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
                      sql <- sprintf("SELECT * FROM %s", tab)
                      message(sql)
                      suppressWarnings(dbGetQuery(dbConn, sql))
                  } )
    
    gr <- makeGRangesFromDataFrame(df,
                                   seqnames.field="genoName",
                                   start.field="genoStart",
                                   end.field="genoEnd",
                                   strand.field="strand",
                                   keep.extra.columns=TRUE)
    return(gr)
}



#' get UCSC mappability data
#' @return GRanges object with lowerLimit, dataRange, upperLimit columns
#' @example uniqueness35 <- get_wgEncodeDukeUniqueness35bp()
get_wgEncodeDukeUniqueness35bp <- function(freeze="hg18",
                                           table="wgEncodeDukeUniqueness35bp",
                                           cols=c("chrom",
                                               "chromStart",
                                               "chromEnd",
                                               "validCount",
                                               "sumData")) {
    
    connect2UCSC <- function(freeze="hg18") {
        junk <- sapply(dbListConnections(MySQL()), dbDisconnect)
        dbConn <- dbConnect(MySQL(),
                            host="genome-mysql.cse.ucsc.edu",
                            dbname=freeze, 
                            username="genome")
        return(dbConn)    
    }
    dbConn <- connect2UCSC()
    
    
    sql <- sprintf("SELECT %s FROM %s",
                   paste(cols, collapse=", "),
                   table)
    message(sql)
    df <- suppressWarnings(dbGetQuery(dbConn, sql))
    df$mappability <- df$sumData/df$validCount
    
    gr <- makeGRangesFromDataFrame(df,
                                   seqnames.field="chrom",
                                   start.field="chromStart",
                                   end.field="chromEnd",
                                   ignore.strand=TRUE,
                                   keep.extra.columns=TRUE)
    
    junk <- sapply(dbListConnections(MySQL()), dbDisconnect)
    return(gr)
}
##uniqueness35 <- get_wgEncodeDukeUniqueness35bp()



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
    sites.l <- bplapply(1:nrow(meta), function(i)
                    {
                        i.df <- try(.load_uniqueSites_from_RData(meta[i,]))
                        if( class(i.df) == "try-error" ) i.df <- data.frame()
                        return(i.df)
                    }
                        ,BPPARAM=MulticoreParam(args$nproc))
    sites <- dplyr::rbind_all(sites.l)
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
    sites.l <- bplapply(1:nrow(meta), function(i)
                    {
                        i.df <- try(.load_multiSites_from_RData(meta[i,]))
                        if( class(i.df) == "try-error" ) i.df <- data.frame()
                        return(i.df)
                    }
                        ,BPPARAM=MulticoreParam(args$nproc))
    sites <- dplyr::rbind_all(sites.l)
}   


