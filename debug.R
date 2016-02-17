compr <- dplyr::group_by(compr, qid)
truth.read.anno <- dplyr::group_by(truth.read.anno, qid)

#### check number of reads aligned(correct or not)
df <- dplyr::summarize(compr, aligned=any(nAln>0))

table(df$aligned, useNA="ifany")

## FALSE   TRUE
##  8469 491531

dplyr::filter(truth.read.anno, aligned) %>% nrow()

## [1] 491531

## PASS


#### check number reads aligned correctly

df <- dplyr::summarize(compr,
                       aligned=any(nAln>0),
                       isAligned=any(nAln>0) & any(hit))

table(df$aligned, useNA="ifany")
## FALSE   TRUE
##  8469 491531

table(df$isAligned, useNA="ifany")
## FALSE   TRUE
## 10073 489927

## PASS

#### check distribution by if they are in repeat region

sum(table(truth.read.anno$isRep, useNA="ifany"))
##[1] 500000



#### check reads aligned incorrectly
df <- dplyr::summarize(compr,
                       aligned=any(nAln>0),
                       isAligned=any(nAln>0) & any(hit),
                       alignedNotRight = aligned & (! isAligned) )

table(df$alignedNotRight, useNA="ifany")
## FALSE   TRUE
##498396   1604


