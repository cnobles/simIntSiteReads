libs <- c("dplyr")
null <- suppressMessages(sapply(libs, library, character.only=TRUE))

site.res <- read.csv("./truth.site.call.txt", header=TRUE, sep="\t")

site.res.nf <- dplyr::filter(site.res,
                             is.na(subjectHits.uniq) &
                             is.na(subjectHits.multi))
message(paste(
    basename(getwd()),
    nrow(site.res.nf), # number of sites not found by caller
    collapse='\t'
    )
)
    
