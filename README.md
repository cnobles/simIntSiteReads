# simIntSiteReads

To create folder with simulated reads run:
```{bash}
Rscript ./simIntSiteReads.R
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
