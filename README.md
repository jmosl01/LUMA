---
title: The Lcms-based Untargeted Metabolomics Assistant (LUMA) R package.
author: 
  
- affiliation: U.S. Environmental Protection Agency
  email: evich.marina@epa.gov
  name: Marina Evich
- affiliation: U.S. Environmental Protection Agency
  email: melendez.wilson@epa.gov
  name: Wilson Melendez
- affiliation: U.S. Environmental Protection Agency
  email: mosley.jonathan@epa.gov
  name: Jonathan Mosley
output: 
  html_document: default
---

  The LCMS-Based Untargeted Metabolomics Assistant (LUMA) package is an automated data processing tool, bridging the outputs of XCMS and CAMERA to a comprehensive data matrix in the R environment for discovery-based studies. The intended users of this R package are LC-MS-based untargeted metabolomics practitioners with minimal to advanced experience in the R environment who process raw LC-MS data files in preparation for downstream statistical analysis. The LUMA package performs feature-reduction by minimizing single features which do not pass quality control (QC) checks and can negatively impact downstream analyses. LUMA contains self-contained functions (herein referred to as modules) that allow for rapid, automated workflows to perform these QC steps with minimal user input. Furthermore, to expedite manual data curation for potentially conflicting isotope and ion adduct annotations, data visualization is consolidated to a single graphic per metabolite group. This graphic contains all EIC plots and psSpectra from CAMERA and new correlation matrices and dendrograms for all features attributed to a single metabolite. Final processed metabolite data, containing normalized intensities and user-defined meta-data, can be exported to worksheets which are directly formatted to a number of analytical tools including MetaboAnalyst, while retaining traceability in the R environment.  


## Disclaimer

*This software/application was developed by the U.S. Environmental Protection Agency (USEPA).  No warranty expressed or implied is made regarding the accuracy or utility of the system, nor shall the act of distribution constitute any such warranty.  The USEPA has relinquished control of the information and no longer has responsibility to protect the integrity, confidentiality or availability of the information.  Any reference to specific commercial products, processes, or services by service mark, trademark, manufacturer, or otherwise, does not constiutite or imply their endorsement, recommendation or favoring by the USEPA.  The USEPA seal and logo shall not be used in any manner to imply endorsement of any commercial product or activity by the USEPA or the United States Government.*


____


<a rel="license" href="http://creativecommons.org/licenses/by-nc-sa/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by-nc-sa/4.0/88x31.png" /></a><br /><span xmlns:dct="http://purl.org/dc/terms/" property="dct:title">LUMA</span> by <span xmlns:cc="http://creativecommons.org/ns#" property="cc:attributionName">USEPA</span> is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by-nc-sa/4.0/">Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License</a>.



[![Travis-CI Build Status](https://travis-ci.com/USEPA/LUMA.svg?branch=master)](https://travis-ci.com/USEPA/LUMA)

### Website

LUMA is public domain and can be freely downloaded from

https://github.com/USEPA/LUMA

Repository is registered in the Reusable Component Services (RCS) system:

https://sor.epa.gov/sor_extranet/registry2/reusereg/searchandretrieve/details/general/24957
