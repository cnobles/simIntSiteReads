#' get sonic lengths from the database
#' @param implied by ~/.my.cnf and args$group
#' @return integer vector of sonic lengths observed
#' @example
#' tmp <- get_sonicLength()
#' 
get_sonicLength <- function() {
    junk <- sapply(dbListConnections(MySQL()), dbDisconnect)
    dbConn <- dbConnect(MySQL(), group=args$group)
    sql <- "
    SELECT ABS(position - breakpoint) AS MLEN
    FROM sites JOIN pcrbreakpoints
    ON sites.siteID = pcrbreakpoints.siteID;"
    res <- dbGetQuery(dbConn,sql)
    
    return( res$MLEN )
}



#' get all sites and sonic abundance from the database
#' @param implied by ~/.my.cnf and args$group
#' @return data.frame(chr, pos, strand, estAbun)
#' @example
#' tmp <- get_sites()
#'
get_sites <- function() {
    junk <- sapply(dbListConnections(MySQL()), dbDisconnect)
    dbConn <- dbConnect(MySQL(), group=args$group)
    sql <- "
    SELECT chr, strand, position, count(*) as estAbun
    FROM sites JOIN pcrbreakpoints
    ON sites.siteID = pcrbreakpoints.siteID 
    GROUP BY sites.siteID;"
    allsites <- dbGetQuery(dbConn,sql)
    allsites <- subset(allsites, grepl("^chr\\d+$", chr))
    
    gr <- makeGRangesFromDataFrame(allsites,
                                   start.field="position",
                                   end.field="position",
                                   keep.extra.columns=TRUE)
    
    ## cobine close by sites
    gr.reduced <- reduce(gr, min.gapwidth=100, with.revmap=TRUE)
    
    ## get estAbun
    gr.estAbun <- gr$estAbun
    gr.reduced$estAbun <- sapply(gr.reduced$revmap, function(revmaplist) {
        max(gr.estAbun[revmaplist])
    } )
    gr.reduced$revmap <- NULL
    
    sites <- as.data.frame(gr.reduced)
    sites <- with(sites, data.frame(chr=as.character(seqnames),
                                    position=start,
                                    strand=as.character(strand),
                                    estAbun=estAbun) )
    sites <- sites[order(-sites$estAbun), ]
    rownames(sites) <- NULL
    
    return(sites)
}


#' get sequence from a position down stream
#' @param sp  species, Hsapiens, etc
#' @param chr chromosome
#' @param pos start position
#' @param strand strand, "+" or "-"
#' @param width number of bases from start
#' @return a dataframe with seq, chr, start, strand, width
#' @note sequences are strand specific:
#'       for + strand, it fetches sp[start, width-1]
#'       for - strand, it fetches RC(sp[start-with+1, start])
#' @require GenomicRanges, BSgenome, BSgenome.Hsapiens.UCSC.hg18, ...etc
#' @example
#' library(BSgenome)
#' library("BSgenome.Hsapiens.UCSC.hg18")
#' get_sequence_downstream(Hsapiens, "chr3", 59165408, "+", 30)
#' get_sequence_downstream(Hsapiens, "chr3", 59165408, "-", 30)
#' get_sequence_downstream(Hsapiens, "chr3", 59165408, "-", c(30,40))
get_sequence_downstream <- function(sp, chr, pos, strand, width) {
    #' debug code
    #' chr="chr3"
    #' pos=59165408
    #' strand="+"
    #' width=30
    #' sp=Hsapiens
    #'
    require(BSgenome)
    stopifnot( class(sp) == "BSgenome" )
    options(stringsAsFactors=FALSE)
    
    ## expand width
    df <- merge(data.frame(chr=chr,
                           pos=pos,
                           strand=strand),
                data.frame(width=width))
    df <- dplyr::arrange(df, chr, pos, strand, width)
    
    ## downstream along strand
    gr <- flank(GRanges(df$chr, IRanges(df$pos, df$pos), df$strand),
                width=df$width,
                start=FALSE,
                ignore.strand=FALSE)
    
    ## shift 1 upstream as flank start from the next base
    gr <- shift(gr, shift=ifelse(strand(gr)=="+", -1, +1))
    
    seq <- getSeq(sp, gr, as.character=TRUE)
    
    return( cbind(seq, df ))
}



#' Given integration sequence and oligo information make R1 R2 I1
#' @param oligo oligo information
#' @param intseq integration information, chr, start, strand, width, seq
#'               the sequences are human genome only 
#' @return dataframe of R1, R2, I1, qname character vectors
#' @note the read names are to be
#'       R1: @M03249:(#chr):(#strand):(#start):(#width):(#index) 1:N:0:0
#'       R2: @M03249:(#chr):(#strand):(#start):(#width):(#index) 2:N:0:0
#'       I1: @M03249:(#chr):(#strand):(#start):(#width):(#index) 1:N:0:0
#'       qname only contains the common part
#' @example make_miseq_reads(oligo, intseq)
make_miseq_reads <- function(oligo, intseq, R1L=175, R2L=130) {
    options(stringsAsFactors=FALSE)
    
    patch_randomGATC <- function(reads, n) {
        unname( sapply(reads, function(seq) {
                           if( nchar(seq)  >= n ) return(seq)
                           return( paste0(seq, paste0(rep("T", n-nchar(seq)), collapse="") ))
                       } ) ) }
    
    molecule_in_miseq <- paste0(oligo$P7,
                                reverseComplement(DNAStringSet(oligo$BC)),
                                oligo$Spacer,
                                oligo$SP2,
                                oligo$Primer,
                                oligo$LTRBit,
                                intseq$seq,
                                reverseComplement(DNAStringSet(oligo$Linker)),
                                reverseComplement(DNAStringSet(oligo$SP1)),
                                reverseComplement(DNAStringSet(oligo$P5)))
    
    R2seq <- substr(molecule_in_miseq, oligo$R2Start, oligo$R2Start+R2L-1)
    R2seq <- patch_randomGATC(R2seq, R2L)
    
    molecule_in_miseq.rc <- reverseComplement(DNAStringSet(molecule_in_miseq))
    R1seq <- substr(molecule_in_miseq.rc, oligo$R1Start, oligo$R1Start+R1L-1)
    R1seq <- patch_randomGATC(R1seq, R1L)
    
    I1seq <- rep(oligo$BC, length(molecule_in_miseq))
    
    ##@M03249:67:000000000-AFCKK:1:1101:14105:1552 2:N:0:0"
    qname <- sprintf("M03249:1:000-SIM%s:1:1:%s:%s",
                     paste0(intseq$chr,
                            ifelse(intseq$strand=="+", "p", "m"),
                            intseq$position),
                     intseq$width,
                     seq_along(I1seq) )
    
    return(data.frame(I1=as.character(I1seq),
                      R1=as.character(R1seq),
                      R2=as.character(R2seq),
                      qname=as.character(qname)))
}



#' make Primary Analysis Directory according to intSiteCaller readme
#' @param df data frame of I1, R1, R2 reads and the common part of qname
#' @param path string of path, default to intSiteSimulation
#' @return nothing returned, a directory is created and ready for analysis
#' @example makeInputFolder(df, "intSiteSimulation")
#'
makeInputFolder <- function(df=df, path="intSiteSimulation") {
    unlink(path, recursive=TRUE, force=TRUE)
    message("Directory ", path, " deleted")
    
    stopifnot(c("I1", "R1", "R2") %in% colnames(df))
    
    ## write fastq files in Data
    dir.create(file.path(path, "Data"),
               recursive=TRUE,
               showWarnings=FALSE)
    
    qnameComments <- c("I1"="1:N:0:0",
                       "R1"="1:N:0:0",
                       "R2"="2:N:0:0") 
    
    pairs <- c("I1", "R1", "R2")
    fastqFile <- setNames(sprintf("Undetermined_S0_L001_%s_001.fastq.gz",
                                  pairs), pairs)
    
    for( pair in pairs ) {
        message("\nWriting ", pair)
        
        read <- ShortReadQ(
            DNAStringSet( df[[pair]] ),
            FastqQuality( sapply(nchar( df[[pair]] ),
                                 function(i) paste(rep("z",i), collapse="")) ),
            BStringSet( paste(df$qname, qnameComments[pair])) )
        
        writeFastq(read,
                   file=file.path(path, "Data", fastqFile[pair]),
                   mode="w", compress=TRUE)
        message(file.path(path, "Data", fastqFile[pair]))
    }
    
    ## copy sampleInfo.tsv
    file.copy(file.path(args$codeDir, "sampleInfo.tsv"),
              file.path(path, "sampleInfo.tsv"),
              overwrite=TRUE)
    
    ## copy processingParams.tsv
    file.copy(file.path(args$codeDir, "processingParams.tsv"),
              file.path(path, "processingParams.tsv"),
              overwrite=TRUE)
    
    ## copy vector fasta
    file.copy(file.path(args$codeDir, "vector_sim.fa"),
              file.path(path, "vector_sim.fa"), overwrite=TRUE )
    
    message("Directory ", path, " created")
}

