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

source(file.path(args$codeDir, "simIntSiteReads_func.R"))


get_info_from_database <- function() {
    allSites <- get_sites()
    sonicLength <- get_sonicLength()
    return( list(site=allSites, sonicLength=sonicLength) )
}
sitesInfo <- get_info_from_database()



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
                         ##vectorSeq="vector_WasLenti.fa",
                         vectorSeq="vector_sim.fa")
sampleInfo <- read.table("sampleInfo.tsv", header=TRUE)

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

I1R1R2qNamedf <- make_miseq_reads(oligo, intseq)

makeInputFolder(I1R1R2qNamedf, sampleInfo, args$outFolder)

print(args)

