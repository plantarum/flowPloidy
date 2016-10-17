# flowPloidy

## Installation

`flowPloidy` is now available in the devel branch of  [Bioconductor](https://bioconductor.org). 

To install it, you need to install the `bioconductor` R package (more details on the [Bioconductor site ](http://bioconductor.org/install/)):

```{r}
## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite()
```

Once that's installed, you need to select the "Devel Branch":

```{r}
library(BiocInstaller)
useDevel()
biocValid()              # checks for out of date packages
biocLite()               # (optional) updates out of date packages
```

(for more details on using the
[Bioconductor Devel Version](http://bioconductor.org/developers/how-to/useDevel/) see the linked page)

With Bioconductor set to the devel version, you can now install `flowPloidy` directly. In order to use the examples in the vignette and the help files, you'll also need to install `flowPloidyData`:

```{r}
biocLite("flowPloidy")
biocLite("flowPloidyData")   # (optional) data for the examples
```

This should pull in all the package dependencies for `flowPloidy`, after
which you can load the package with the normal function
`library("flowPloidy")`.

## Getting Started

```{r}
library("flowPloidy")
```

The `flowPloidy` workflow is documented in the vignette, which you can view from R:

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
