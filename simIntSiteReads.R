#### simulate R1, R2, I1 fastq files for intSiteCaller, and save the truth file ####
#' Approach:
#' 1. Get 1000 sites, either random or from intsites_miseq database;
#' 2. Simulate template lengths ~ N(300, 100) within [30, 1000], 100 templates for each site;
#' 3. Get sequences for the templates, attach primer and linkers according to strands;
#' 4. Get R1 and R2 according to strand, simulate I1 based on barcode;
#' 5' Read names assigned to
#'    R1: @M03249:site(#id):template(#id):(#chr):(#strand):(#pos) 1:N:0:0
#'    R2: @M03249:site(#id):template(#id):(#chr):(#strand):(#pos) 2:N:0:0
#'    I1: @M03249:site(#id):template(#id):(#chr):(#strand):(#pos) 1:N:0:0
#' 6' Write to fastq.gz files
#'

#' set default and/or get arguments from command line input
#' @param commandline arguments or none for interactive session
#' @return list of default values
#'         freeze reference genome
#'         sites number to integration sites to simulate
#'         sonicLength number of templates to simulate for a site
#'         group group in the ~/.my.cnf file

options(stringsAsFactors=FALSE)

get_args <- function() {
    suppressMessages(library(argparse))
    codeDir <- dirname(sub("--file=", "", grep("--file=", commandArgs(trailingOnly=FALSE), value=T)))
    if( length(codeDir) == 0 ) codeDir <- "."
    parser <- ArgumentParser(formatter_class='argparse.RawTextHelpFormatter')
    parser$add_argument("-f", "--freeze", type="character", nargs=1,
                        default="hg18",
                        help="reference, hg18, hg19, etc")
    parser$add_argument("-o", "--outFolder", type="character", nargs=1,
                        default="intSiteSimulation",
                        help="output folder")
    parser$add_argument("-s", "--sites", type="integer", nargs=1,
                        default=1000,
                        help="number of integration sites")
    parser$add_argument("-l", "--sonicLength", type="integer", nargs=1,
                        default=100,
                        help="number of sonic lengths for each site")
    parser$add_argument("-1", "--R1L", type="integer", nargs=1,
                        default=175,
                        help="R1 read length")
    parser$add_argument("-2", "--R2L", type="integer", nargs=1,
                        default=130,
                        help="R2 read length")
    parser$add_argument("-g", "--group", type="character", nargs=1,
                        default="intsites_miseq.read",
                        help="number of sonic lengths for each site")
    parser$add_argument("-c", "--codeDir", type="character", nargs=1,
                        default=codeDir,
                        help="Directory of code")
    args <- parser$parse_args(commandArgs(trailingOnly=TRUE))
}
args <- get_args()
##print(args)

libs <- c("RMySQL",
          "ggplot2",
          "GenomicRanges",
          "ShortRead",
          "BSgenome",
          sprintf("BSgenome.Hsapiens.UCSC.%s", args$freeze))
null <- suppressMessages(sapply(libs, library, character.only=TRUE))

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

get_info_from_database <- function() {
    allSites <- get_sites()
    sonicLength <- get_sonicLength()
    return( list(site=allSites, sonicLength=sonicLength) )
}




#' get sequence from a position down stream
#' @param x reference, Hsapiens, etc
#' @param chr chromosome
#' @param start start position
#' @param strand strand, "+" or "-"
#' @param width number of bases from start
#' @return a dataframe with seq, chr, start, strand, width
#' @note sequences are strand specific:
#'       for + strand, it fetches x[start, width-1]
#'       for - strand, it fetches RC(x[start-with+1, start])
#' @require GenomicRanges, BSgenome, BSgenome.Hsapiens.UCSC.hg18, ...etc
#' @example
#' library(BSgenome)
#' library("BSgenome.Hsapiens.UCSC.hg18")
#' get_sequence_downstream(Hsapiens, "chr3", 59165408, "+", 30)
#' get_sequence_downstream(Hsapiens, "chr3", 59165408, "-", 30)
#' get_sequence_downstream(Hsapiens, "chr3", 59165408, "-", c(30,40))
get_sequence_downstream <- function(x, chr, start, strand, width) {
    require(BSgenome)
    stopifnot( class(x) == "BSgenome" )
    options(stringsAsFactors=FALSE)
    
    df <- merge(data.frame(chr=site$chr,
                           start=site$position,
                           strand=site$strand),
                data.frame(width=width))
    df <- dplyr::arrange(df, chr, start, strand, width)
    
    chr <- df$chr
    start <- df$start
    strand <- df$strand
    width <- df$width
    
    seq <- with(df, as.character( getSeq(
        x,
        shift(flank(GRanges(chr, IRanges(start, start), strand),
                    width=width, start=FALSE, ignore.strand=FALSE),
              shift= ifelse(strand=="+", -1, +1) ) )
                                 ) )
    return( cbind(seq, df ))
}



#' @note this is from one line of the sample information file
#' alias,linkerSequence,bcSeq,gender,primer,ltrBit,largeLTRFrag,vectorSeq
#' GTSP0308-1,GAACGAGCACTAGTAAGCCCNNNNNNNNNNNNCTCCGCTTAAGGGACT,GTATTCGACTTG,m,GAAAATC,TCTAGCA,TGCTAGAGATTTTCCACACTGACTAAAAGGGTCT,vector_WasLenti.fa
sampleInfo <- data.frame(alias="GTSP0308-1",
                         ##linkerSequence="GAACGAGCACTAGTAAGCCCNNNNNNNNNNNNCTCCGCTTAAGGGACT",
                         linkerSequence="GAACGAGCACTAGTAAGCCCGGGGGGTTTTTTCTCCGCTTAAGGGACT",
                         bcSeq="GTATTCGACTTG",
                         gender="m",
                         primer="GAAAATC",
                         ltrBit="TCTAGCA",
                         largeLTRFrag="TGCTAGAGATTTTCCACACTGACTAAAAGGGTCT",
                         vectorSeq="vector_WasLenti.fa")

#' @note this is specific to the integration protocol
oligo <- data.frame(P5="AATGATACGGCGACCACCGA",
                    P7="CAAGCAGAAGACGGCATACGAGAT",
                    Spacer="AGTCAGTCAGCC",
                    SP1="ATCTACACCAGGACTGACGCTATGGTAATTGT",
                    SP2="AGACCCTTTTAGTCAGTGTG",
                    IP="CACACTGACTAAAAGGGTCTGGCTGACTGACT",
                    Linker=sampleInfo$linkerSequence,
                    BC=sampleInfo$bcSeq,
                    Primer=sampleInfo$primer,
                    LTRBit=sampleInfo$ltrBit)
oligo$R2Start <- with(oligo, 1+nchar(paste0(P7, BC, Spacer, SP2))) #! 1 based
oligo$R1Start <- with(oligo, 1+nchar(paste0(P5, SP1))) #! 1 based


## validated known sites from a prep
##clone1	Clonal 293T cells with known integration at chr1+52699700
##clone2	Clonal 293T cells with known integration at chr17+77440127
##clone3	Clonal 293T cells with known integration at chr19+1330529
##clone4	Clonal 293T cells with known integration at chr1-153461600
##clone7	Clonal 293T cells with known integration at chr1+148889088
site <- data.frame(chr=c("chr1", "chr17", "chr19", "chr1", "chr1"),
                   position=c(52699700, 77440127, 1330529, 153461600, 148889088),
                   strand=c("+", "+", "+", "-", "+") )


## get sequence of integration, downstream from the point of integration
##sitesInfo <- get_info_from_database()
##site <- tail(head(sitesInfo$site, 3), 1)
##site <- sitesInfo$site[3,]
##width <- sample(sitesInfo$sonicLength, 100, replace=TRUE)
width <- sample(200:1000, 100, replace=FALSE)
width <- c(30:1000)

intseq <- get_sequence_downstream(Hsapiens,
                                  site$chr,
                                  site$position,
                                  site$strand,
                                  width)


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
make_miseq_reads <- function(oligo, intseq) {
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
    
    R2seq <- substr(molecule_in_miseq, oligo$R2Start, oligo$R2Start+args$R2L-1)
    R2seq <- patch_randomGATC(R2seq, args$R2L)
    
    molecule_in_miseq.rc <- reverseComplement(DNAStringSet(molecule_in_miseq))
    R1seq <- substr(molecule_in_miseq.rc, oligo$R1Start, oligo$R1Start+args$R1L-1)
    R1seq <- patch_randomGATC(R1seq, args$R1L)
    
    I1seq <- rep(oligo$BC, length(molecule_in_miseq))
    
    ##@M03249:67:000000000-AFCKK:1:1101:14105:1552 2:N:0:0"
    qname <- sprintf("M03249:1:000-SIM%s:1:1:%s:%s",
                     paste0(intseq$chr,
                            ifelse(intseq$strand=="+", "p", "m"),
                            intseq$start),
                     intseq$width,
                     seq_along(I1seq) )
    
    return(data.frame(I1=as.character(I1seq),
                      R1=as.character(R1seq),
                      R2=as.character(R2seq),
                      qname=as.character(qname)))
}
df <- make_miseq_reads(oligo, intseq)



#' make Primary Analysis Directory according to intSiteCaller readme
#' @param df data frame of I1, R1, R2 reads and the common part of qname 
#' @param path string of path
#' @return nothing returned, a directory is created and ready for analysis
#' @example makeInputFolder(df, "intSiteSimulation")
makeInputFolder <- function(df=df, path="intSiteSimulation") {
    unlink(path, recursive=TRUE, force=TRUE)
    
    ## make Data directory
    mapply( function(pair, comment) {
        read <- ShortReadQ(DNAStringSet( df[[pair]] ),
                           FastqQuality( sapply(nchar( df[[pair]] ), function(i) paste(rep("z",i), collapse="")) ),
                           BStringSet( paste(df$qname, comment) ) )
        
        dir.create(file.path(path, "Data"),
                   recursive=TRUE,
                   showWarnings=FALSE)
        
        file <- file.path(path, "Data",
                          sprintf("Undetermined_S0_L001_%s_001.fastq.gz", pair) )
        
        suppressWarnings(file.remove(file))
        writeFastq(read, file=file, mode="w", compress=TRUE)
    },
           pair=c("I1", "R1", "R2"),
           comment=c("1:N:0:0", "1:N:0:0", "2:N:0:0") )
    
    ## sampleInfo.tsv
    sampleInfo <- data.frame(alias="GTSP0308-1",
                             linkerSequence="GAACGAGCACTAGTAAGCCCNNNNNNNNNNNNCTCCGCTTAAGGGACT",
                             bcSeq="GTATTCGACTTG",
                             gender="m",
                             primer="GAAAATC",
                             ltrBit="TCTAGCA",
                             largeLTRFrag="TGCTAGAGATTTTCCACACTGACTAAAAGGGTCT",
                             vectorSeq="vector_sim.fa")
    write.table(sampleInfo, file.path(path, "sampleInfo.tsv"),
              quote=FALSE, sep="\t", row.names=FALSE )
    
    ## processingParams.tsv
    processingParams <- data.frame(qualityThreshold="?",
                                   badQualityBases="5",
                                   qualitySlidingWindow="10",
                                   mingDNA="30",
                                   minPctIdent="95",
                                   maxAlignStart="5",
                                   maxFragLength="2500",
                                   refGenome="hg18")
    write.table(processingParams, file.path(path, "processingParams.tsv"),
                quote=FALSE, sep="\t", row.names=FALSE )
    
    ## vector fasta
    file.copy(file.path(args$codeDir, "vector_sim.fa"),
              file.path(path, "vector_sim.fa"), overwrite=TRUE )
    
    message("Directory ", path, " created")
}
unlink(args$outFolder, recursive=TRUE, force=TRUE)
makeInputFolder(df, args$outFolder)

print(args)

