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
                          file.path(codeDir, "checkRunResults.R"),
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
    parser$add_argument("-n", "--nproc", type="integer", nargs=1,
                        default=5,
                        help="tolerance for alignment")
    
    args <- parser$parse_args(commandArgs(trailingOnly=TRUE))
    
    args$workDir <- normalizePath(args$workDir, mustWork=TRUE)
    return(args)
}
args <- get_args()
print(t(as.data.frame(args)), quote=FALSE)

libs <- c("stringr",
          "plyr",
          "dplyr",
          "RMySQL",
          "GenomicRanges",
          "ShortRead",
          "BiocParallel",
          "ggplot2")
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
    I1File <-  grep("_I1_", fastqFiles, value=TRUE)
    R1File <-  grep("_R1_", fastqFiles, value=TRUE)
    R2File <-  grep("_R2_", fastqFiles, value=TRUE)
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


load_run_stats <- function(workDir=".") {
    
    stats.file <- list.files(workDir, pattern="^stats.RData$",
                             recursive=TRUE, full.names=TRUE)
    
    tmp.statlist <- lapply(setNames(stats.file, stats.file), function(x) {
        a <- load(x)
        get(a)
    })
    stats <- plyr:::rbind.fill(tmp.statlist)
    stats$sample <- as.character(stats$sample)
    rownames(stats) <- NULL
    
    return( c("demultiplexed"=sum(stats$barcoded),
              "LTRed"=sum(stats$LTRed),
              "linkered"=sum(stats$linkered),
              "LTRed.linkered"=sum(stats$ltredlinkered),
              "lengthTrimed"=sum(stats$lenTrim),
              "vectorTrimed"=sum(stats$vTrimed)) )
    
}


#' check total number of decoded reads
#' @param metadata, metadata
#' @return named vector of numbers of reads in R1 and R2
check_demultiplexed_reads <- function(metadata) {
    sumR1 <- summary(readFastq(metadata$R1fastq))
    sumR2 <- summary(readFastq(metadata$R2fastq))
    return(c('nR1'=as.integer(sumR1["Length"]),
             'nR2'=as.integer(sumR2["Length"])))
}


#' check run time by looking at logs
#' time start is logged in logs/errorCorrectOutput.txt
#' time end is taken as the latest from logs/callSitesOutput*.txt
#' @return time elasped in seconds
#' 
check_runtime_by_logs <- function() {
    format <- "%a %b %d %H:%M:%S %Y"
    ##x1 <- "Mon Nov  9 10:16:39 2015"
    ##x2 <- "Mon Nov  9 10:16:58 2015"
    
    cmd <- "head -15 logs/errorCorrectOutput.txt | grep Started"
    cout <- system(cmd, intern=TRUE)
    x1 <- sub("Start.*at ", "", cout)
    t1 <- strptime(x1, format=format)
    
    cmd <- "head -15 logs/callSitesOutput*.txt | grep reported | grep Results"
    cout <- system(cmd, intern=TRUE)
    x2 <- sub("Results reported on ", "", cout)
    t2 <- head(sort(strptime(x2, format=format), decreasing=TRUE), 1)
    
    tsecs <- as.integer(difftime(t2,t1, units="secs"))
    return(tsecs)
}
##check_runtime_by_logs()


#' check blat parameter from alignmentg log
#' @return time elasped in seconds
#' 
check_blat_param_by_log <- function() {
    cmd <- "grep blat logs/alignOutput1.txt | grep psl"
    cout <- system(cmd, intern=TRUE)
    
    tileSize <- as.integer(stringr::str_match(cout, "tileSize=(\\d+)")[2])
    stepSize <- as.integer(stringr::str_match(cout, "stepSize=(\\d+)")[2])
    minIdentity <- as.integer(stringr::str_match(cout, "minIdentity=(\\d+)")[2])
    maxIntron <- as.integer(stringr::str_match(cout, "maxIntron=(\\d+)")[2])
    minScore <- as.integer(stringr::str_match(cout, "minScore=(\\d+)")[2])
    
    param <- c("tileSize"=tileSize,
               "stepSize"=stepSize,
               "minIdentity"=minIdentity,
               "maxIntron"=maxIntron,
               "minScore"=minScore)
    param.name <- names(param)
    param.default <- c(11, 11, 30, 100000, 90)
    
    param <- ifelse(is.na(param), param.default, param)
    names(param) <- param.name
    
    return(param)
}
##check_blat_param_by_log()


#### load data, truth and results ####
#### site level comparison ####

metadata <- cbind(get_metadata(),
                  get_machine_file(dir=file.path(args$workDir,"Data")))

truth <- load_truth_from_fastq(metadata)
stopifnot(!any(duplicated(truth$qid)))



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

save.image(file = "debug.checkResults.RData")


## find false positives
uniq.resIntruth.ovl <- findOverlaps(res.uniq.site.gr,
                                    truth.site.gr,
                                    maxgap=args$err)
res.uniq.site.fp.gr <- res.uniq.site.gr[-queryHits(uniq.resIntruth.ovl)]


res.multi.site.grl <- split(res.multi.site.gr, res.multi.site.gr$multihitID)
multi.resIntruth.ovl <- findOverlaps(res.multi.site.grl,
                                     truth.site.gr,
                                     maxgap=args$err)
res.multi.site.fp.grl <- res.multi.site.grl[-queryHits(multi.resIntruth.ovl)]


call.site.stat <-c(
    "all.sim"= nrow(truth.site),
    "found.any"=length(unique(all.site.ovl$queryHits)),
    "found.uniq"=length(unique(subset(all.site.ovl, !is.na(subjectHits.uniq))$queryHits)),
    "found.uniq.only"=length(unique(subset(all.site.ovl, !is.na(subjectHits.uniq) &  is.na(subjectHits.multi))$queryHits)),
    "found.multi"=length(unique(subset(all.site.ovl, !is.na(subjectHits.multi))$queryHits)),
    "found.multi.only"=length(unique(subset(all.site.ovl, !is.na(subjectHits.multi) & is.na(subjectHits.uniq)))$queryHits),
    "found.both"=length(unique(subset(all.site.ovl, !is.na(subjectHits.multi) & !is.na(subjectHits.uniq))$queryHits)),
    "falseP.uniq"=length(res.uniq.site.fp.gr),
    "falseP.multi"=length(res.multi.site.fp.grl)
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

save.image(file = "debug.checkResults.RData")


#### reads level comparison ####

res.all <- dplyr::rbind_list(res.uniq, res.multi)
stopifnot(all(res.all$qid %in% truth$qid))

compr <- merge(truth, res.all, by=c("qid"), all.x=TRUE)

is_good_hit <- function(x, err=3) {
    is.good <- with(x, 
                    chr.x==chr.y &
                        strand.x==strand.y &
                            abs(position.x-position.y)<=err &
                                abs(breakpoint.x-breakpoint.y)<=err )
    is.good[is.na(is.good)] <- FALSE
    return(is.good)
}

compr$hit <- is_good_hit(compr, args$err)

compr <- (dplyr::group_by(compr, qid) %>%
          dplyr::mutate(nAln=sum(!is.na(chr.y)),
                        anyHit=any(hit)) )

## reads aligned, not necessarily correct
reads.aligned <- (dplyr::filter(compr, nAln>0) %>%
                  dplyr::select(qid) %>%
                  dplyr::distinct() %>%
                  nrow())

## reads aligned correctly
reads.aligned.right <- (dplyr::filter(compr, hit) %>%
                        dplyr::select(qid) %>%
                        dplyr::distinct() %>%
                        nrow())


## reads aligned correctly and uniquely
reads.aligned.right.uniq <- (dplyr::filter(compr, hit & multihitID==0) %>%
                             dplyr::select(qid) %>%
                             dplyr::distinct() %>%
                             nrow())

## reads aligned correctly within one of multi hits
reads.aligned.right.multi <- (dplyr::filter(compr, hit & multihitID!=0) %>%
                              dplyr::select(qid) %>%
                              dplyr::distinct() %>%
                              nrow())

## site abundance recovery
site.abund <- (dplyr::select(compr, chr.x, strand.x, position.x, anyHit) %>%
               dplyr::distinct() %>% 
               dplyr::group_by(chr.x, strand.x, position.x) %>%
               dplyr::summarize(n=length(unique(qid)),
                                n.rec=sum(anyHit) ))


stopifnot(reads.aligned.right==reads.aligned.right.uniq+reads.aligned.right.multi)

runStat <- load_run_stats()

call.site.stat <- c(call.site.stat,
                    "reads.all"=nrow(truth),
                    runStat,
                    "reads.aligned"=reads.aligned,
                    "reads.aligned.right"=reads.aligned.right,
                    "reads.aligned.right.uniq"=reads.aligned.right.uniq,
                    "reads.aligned.right.multi"=reads.aligned.right.multi,
                    "runtime"=check_runtime_by_logs(),
                    check_blat_param_by_log())


print(as.data.frame(call.site.stat))
message("\nCaller site stat written callstat.txt" )
write.table(as.data.frame(call.site.stat), file="callstat.txt",
            quote=FALSE, sep="\t", col.name=FALSE)
save.image(file = "debug.checkResults.RData")


#### length recovered ####
x <- unlist(c(dplyr::ungroup(compr) %>%
              dplyr::filter(hit) %>%
              dplyr::select(width.y)))

widthcount.df <- merge(as.data.frame(table(width=truth$width)),
                       as.data.frame(table(width=x)),
                       by="width",
                       all.x=TRUE,
                       all.y=TRUE)

widthcount.df$width <- as.integer(as.character((widthcount.df$width)))
widthcount.df[is.na(widthcount.df)] <- 0


theme_text <- theme(text = element_text(size=14),
                    axis.text.x = element_text(size=14, face="bold"),
                    axis.title.x = element_text(size=14, face="bold"),
                    axis.text.y = element_text(size=14, face="bold"),
                    axis.title.y = element_text(size=14, face="bold"))

p <- (ggplot(widthcount.df, aes(x=width, y=Freq.x)) +
      geom_line(color="red") +
      geom_line(aes(x=width, y=Freq.y), color="blue") +
      annotate(geom="text", x=600, y=800, hjust=0,
               label=paste(sprintf("%s\t=\t%s", names(call.site.stat), call.site.stat), collapse="\n"), fontface="bold")+
      xlab("Width")+
      ylab("Read counts")+
      theme_bw() +
      theme_text )

ggsave(filename="ReadsRecoveredByLength.pdf",
       plot=p,
       width=10, height=8, units="in")
ggsave(filename="ReadsRecoveredByLength.png",
       plot=p,
       width=10, height=8, units="in")
save(widthcount.df, file="widthcount.df.RData")

#### site abundance recovered ####
cols <- c("recovered"="blue", "simulated"="red")
##cols <- factor(cols)
p2 <- (ggplot(arrange(ungroup(site.abund), n.rec), aes(x=1:length(n.rec))) +
       geom_line(aes(y=n, color="simulated")) +
       geom_line(aes(y=n.rec, color="recovered")) +
       xlab("Site index (ordered to guide the eye)")+
       ylab("Breakpoints recovered")+
       ##scale_color_manual(name="Lines", values=cols)+
       scale_color_manual(name="Lines", limits=c("simulated", "recovered"), values=cols)+
       theme_bw() +
       theme_text +
       theme(legend.justification=c(1,0), legend.position=c(1,0),
       legend.text=element_text(size=14, face="bold")) )
##p2
ggsave(filename="ReadsRecoveredBySite.pdf",
       plot=p2,
       width=10, height=8, units="in")
ggsave(filename="ReadsRecoveredBySite.png",
       plot=p2,
       width=10, height=8, units="in")
save(site.abund, file="site.abun.RData")


#########################################################
save.image(file = "debug.checkResults.RData")
##q()

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
rmsk.gr <- get_repeakMasker()


truth.site.call <- read.table("truth.site.call.txt", header=TRUE)
truth.site.call.gr <- makeGRangesFromDataFrame(truth.site.call,
                                               start.field="position",
                                               end="position",
                                               strand.field="strand",
                                               keep.extra.columns=TRUE)

rmOvl <- findOverlaps(query=truth.site.call.gr,
                      subject=rmsk.gr,
                      maxgap=args$err,
                      ignore.strand=TRUE)

truth.site.call.gr$repClass <- NA
truth.site.call.gr$repFamily <- NA

truth.site.call.gr$repClass[queryHits(rmOvl)] <- rmsk.gr$repClass[subjectHits(rmOvl)]
truth.site.call.gr$repFamily[queryHits(rmOvl)] <- rmsk.gr$repFamily[subjectHits(rmOvl)]


truth.site.call <- dplyr::select(as.data.frame(truth.site.call.gr),
                                 chr=seqnames,
                                 position=start,
                                 strand,
                                 queryHits,
                                 subjectHits.uniq,
                                 subjectHits.multi,
                                 repClass,
                                 repFamily)


message("\nAnnotating the sites with repeatMasker track")
message("File truth.site.call.txt updated")
write.table(truth.site.call, file="truth.site.call.txt",
            quote=FALSE, sep="\t", col.name=TRUE, row.name=FALSE)

save.image(file = "debug.checkResults.RData")
q()

