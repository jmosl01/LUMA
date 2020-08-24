<!-- badges: start -->
  [![Travis-CI Build Status](https://travis-ci.com/jmosl01/LUMA.svg?branch=master)](https://travis-ci.com/jmosl01/LUMA)
  [![Codecov test coverage](https://codecov.io/gh/jmosl01/LUMA/branch/master/graph/badge.svg)](https://codecov.io/gh/jmosl01/LUMA?branch=master)
  <!-- badges: end -->


# LUMA

  The `LUMA` (LCMS-Based Untargeted Metabolomics Assistant) R package is an automated data processing tool, bridging the outputs of `XCMS` and `CAMERA` to a comprehensive data matrix in the R environment for discovery-based studies. The intended users of this R package are LC-MS-based untargeted metabolomics practitioners with minimal to advanced experience in the R environment who process raw LC-MS data files in preparation for downstream statistical analysis. The `LUMA` package performs feature-reduction by minimizing single features which do not pass quality control (QC) checks and can negatively impact downstream analyses. This package contains self-contained functions (herein referred to as modules) that allow for rapid, automated workflows to perform these QC steps with minimal user input. Furthermore, to expedite manual data curation for potentially conflicting isotope and ion adduct annotations, data visualization is consolidated to a single graphic per metabolite group. This graphic contains all EIC plots and psSpectra from `CAMERA` and new correlation matrices and dendrograms for all features attributed to a single metabolite. Final processed metabolite data, containing normalized intensities and user-defined meta-data, can be exported to worksheets which are directly formatted to a number of analytical tools including `MetaboAnalyst`, while retaining traceability in the R environment.  


## Installation

You can install development version of`LUMA` from GitHub with the following:

```r
# requires devtools to install
install.packages('devtools')
library(devtools)

# install from repository
install_github('USEPA/LUMA')
library(LUMA)
```

To install from GitHub with package vignettes:
```r
library(devtools)
install_github('USEPA/LUMA', build_vignettes=TRUE)
library(LUMA)
```

To install from GitHub on Linux or MacOSX:
```r
library(devtools)
install_github('USEPA/LUMA', type = "source")
library(LUMA)
```


## Example
An overview of the `LUMA` package is provided in the users guide that is included with the package.  The documentation includes a number of examples for use of the various functions and self-contained modules to enable workflow creation.  Vignettes are also available for typical LUMA workflows.



## Package Contributions
We encourage users to submit issues and enhancement requests so we may
continue to improve our package.


## Citation

Please cite the LUMA package when using for data processing:

```{r}
citation("LUMA")

```


## Repositories

The source code for this repository is maintained at https://github.com/jmosl01/LUMA which is also mirrored at https://github.com/usepa/luma


## EPA Disclaimer

*This software/application was developed by the U.S. Environmental Protection Agency (USEPA).  No warranty expressed or implied is made regarding the accuracy or utility of the system, nor shall the act of distribution constitute any such warranty.  The USEPA has relinquished control of the information and no longer has responsibility to protect the integrity, confidentiality or availability of the information.  Any reference to specific commercial products, processes, or services by service mark, trademark, manufacturer, or otherwise, does not constitute or imply their endorsement, recommendation or favoring by the USEPA.  The USEPA seal and logo shall not be used in any manner to imply endorsement of any commercial product or activity by the USEPA or the United States Government.*

____


### License

<a rel="license" href="http://creativecommons.org/licenses/by-nc-sa/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by-nc-sa/4.0/88x31.png" /></a><br /><span xmlns:dct="http://purl.org/dc/terms/" property="dct:title">LUMA</span> by <span xmlns:cc="http://creativecommons.org/ns#" property="cc:attributionName">USEPA</span> is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by-nc-sa/4.0/">Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License</a>.

