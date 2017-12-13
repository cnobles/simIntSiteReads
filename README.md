# simIntSiteReads

To create folder with simulated reads run:
```{bash}
Rscript ./simIntSiteReads.R
```

## Command-line Arguments

We can introduce uniform random error in reads by option `-e`:
```{bash}
Rscript ./simIntSiteReads.R -e 1
```

Molecule length can be uniform or Maxwell-Boltzmann(MB).
MB emulate real skewed read distribution:
```{bash}
Rscript ./simIntSiteReads.R -w uniform
Rscript ./simIntSiteReads.R -w maxwell_boltzmann 
```

## Sample Info dependencies
Below are columns that must contain sequences in the sampleInfo to generate simulated reads:
* linkerSequence
* bcSeq
* primer
* ltrBit

## Template Structure
```
           R1 -------------->...                         I1 ->
P5 - SP1 - linkerSequence - hg - ltrBit - primer - SP2 - bcSeq - P7
                          ...<--------------- R2         
```

# Testing

In R console run:
```{r}
library(testthat)
test_dir(".")
```

Some of the test components are not deterministic,
so once in a while test can fail and we need to rerun it.

# Dependency
Library `stringdist` is required for testing.
