# flowPloidy

## Caution!

The development version of flowPloidy is working well, and corrects a
number of issues present in the stable version. However, the documentation
is now (2017-03-03) out-dated! Work to update the help files and vignettes
is ongoing.

## Installation

### Stable Version

`flowPloidy` is available in [Bioconductor](https://bioconductor.org).

To install it, you need to install the `bioconductor` R package (more
details on the [Bioconductor site ](http://bioconductor.org/install/)):

```{r}
## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite()
```

Once that's installed, you can install `flowPloidy` using the Bioconductor
tools:

```{r}
biocLite("flowPloidy")
biocLite("flowPloidyData")   # (optional) data for the examples
```

This should pull in all the package dependencies for `flowPloidy`, after
which you can load the package with the normal function
`library("flowPloidy")`.

### Development Version

Bioconductor focuses on stable software, and consequently only releases
updates once every 6 months. Development on `flowPloidy` is currently
(March 2017) proceeding at a faster rate. You may install the development
version directly from the GitHub repository if you wish to use the latest
updates:

```{r}
## Install Bioconductor tools first:
source("https://bioconductor.org/biocLite.R")
biocLite()

## Install flowCore from Bioconductor:
biocLite("flowCore")

## Install devtools so you can directly access GitHub
install.packages(devtools)
library(devtools)

## Install flowPloidy:
install_github("plantarum/flowPloidy", dependencies = TRUE, 
    build_vignettes = TRUE)
```

## Getting Started

```{r}
library("flowPloidy")
```

The `flowPloidy` workflow is documented in the vignette, which you can view
from R:

```{r}
fpVig <- vignette("flowPloidy-overview")
fpVig ## open vignette in a browser
edit(name = fpVig) ## open vignette source code in a text editor
```

## Getting Help
For general help using the package, you can post questions on
the [Bioconductor Support Site](https://support.bioconductor.org/). Use the
tag `flowploidy` to ensure your question is brought to my attention.

The development repository for `flowPloidy` is
on [Github](https://github.com/plantarum/flowPloidy), and you can file bugs
there using the **issues** tab. You are also welcome to contribute features
or bug-fixes via pull requests!
