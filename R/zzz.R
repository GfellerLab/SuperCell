.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Thank you for installing SuperCell! SuperCell is developed by the group of David Gfeller at University of Lausanne.

SuperCell can be used freely by academic groups for non-commercial purposes (see license). The product is provided free of charge, and, therefore, on an 'as is' basis, without warranty of any kind.

FOR-PROFIT USERS

If you plan to use SuperCell or any data provided with the script in any for-profit application, you are required to obtain a separate license. To do so, please contact eauffarth@licr.org at the Ludwig Institute for Cancer Research Ltd.

If required, FOR-PROFIT USERS are also expected to have proper licenses for the tools used in SuperCell, including the R packages igraph, RANN, WeightedCluster, corpora, weights, Hmisc, Matrix, ply, irlba, grDevices, patchwork, ggplot2 and velocyto.R

For scientific questions, please contact Mariia Bilous (mariia.bilous@unil.ch) or David Gfeller (David.Gfeller@unil.ch).")
}