# Install Instructions

First, you will need to install the `Batman` package (and all its dependencies) to reproduce the results. The packages are compiled from source. On Windows, you must install [RTools](https://cran.r-project.org/bin/windows/Rtools/rtools44/rtools.html) first, while on OSX you will need at least XCode and a Fortran compiler ([see here](https://cran.r-project.org/bin/macosx/tools/)).

## Dependencies

The dependencies required are:

```{r}
deps <- c("Rcpp", "RcppParallel", "rstan", "rstantools", "StanHeaders", "BH",
          "RcppArmadillo")

install.packages(deps)
```

## Installing Batman

First, the `Batman` package needs to be installed:

```{r}
# install.packages("devtools")
devtools::install_github("theodds/Batman")
```

## Installing BART4RS

You will need to install the `BART4RS` package and its dependencies. The easiest way to do this is to open the `.Rproj` file in the `BART4RS/` directory in RStudio and then use the Build pane. Alternatively, you can run

```{r}
devtools::install_local("BART4RS")
```

# Vignettes

The folder `Code/Markdown` contains two R markdown files that illustrate the use of the proportional hazards (PH) relative survival model and the use of the non-proportional hazards (NPH) relative survival model. After installation, these files can be opened and rendered in RStudio; they illustrate basic usage of the package on the `LeukSurv` dataset from the package `spBayesSurv`.

# Reproducing Results

The file `Code/main.R` runs both our illustrative analysis on the leukemia dataset and our simulation studies. These can also be run individually by running the individual files in `Code/src/`; the files are numbered according to the order in which they should be run. Resulting figures are created in the `Code/figures` directory, and expensive computations are cached in the `Code/cache/` directory.

These files should be run with the working directory in __R__ set to the `Code/` directory.
