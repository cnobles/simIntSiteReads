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
