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

# Dependency
Library `stringdist` is required for testing.
