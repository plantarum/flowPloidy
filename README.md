# flowPloidy

## Prerequisites

Before you can install `flowPloidy`, you need to install `flowCore`
from [Bioconductor](https://bioconductor.org):

```{r}
source("https://bioconductor.org/biocLite.R")
biocLite()
biocLite("flowCore")
```

and `devtools` from CRAN:

```{r}
install.packages("devtools", dependencies = TRUE)
install.package("codetools")
```

(`codetools` seems to be automatically installed sometimes, but not others.
So to be safe, install it explicitly).

In order to build the vignettes, you may need to have `pandoc` installed as
well. On Windows, the easiest way to do this is likely just
installing [RStudio](https://www.rstudio.com).

## Installation

```{r}
library(devtools)
install_github("plantarum/flowPloidy", dependencies = TRUE, 
                  build_vignettes = TRUE)
```

If you don't have the necessary tools to build the vignettes, you can omit them:

```{r}
install_github("plantarum/flowPloidy", dependencies = TRUE)
```

## Getting Started

```{r}
library(flowPloidy)
```

If you built the vignette, you can view it with:

```{r}
fpVig <- vignette("flowPloidy-overview")
fpVig ## open vignette in a browser
edit(name = fpVig) ## open vignette source code in a text editor
```
