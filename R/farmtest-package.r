#' FarmTest: Factor Adjusted Robust Multiple Testing
#'
#' This R package conducts multiple hypothesis testing of mean effects. It implements a robust procedure to estimate distribution parameters and accounts for strong dependence among coordinates via an approximate factor model. This method is particularly suitable for high-dimensional data when there are thousands of variables but only a small number of observations available. Moreover, the method is tailored to cases when the underlying distribution deviates from Gaussianity, which is commonly assumed in the literature.
#' For detailed information on how to use and install see \url{https://kbose28.github.io/FarmTest/}. See the papers on the 'FarmTest' method, Fan et al.(2017) at \url{https://arxiv.org/abs/1711.05386} and Zhou at al.(2017) at \url{https://arxiv.org/abs/1711.05381}, for detailed description of methods and further references.
#'
#'
#'@references Fan, J., Ke, Y., Sun, Q. and Zhou, W-X. (2017). "FARM-Test: Factor-adjusted robust multiple testing with false discovery control," \url{https://arxiv.org/abs/1711.05386}.
#'@references Zhou, W-X., Bose, K., Fan, J. and Liu, H. (2017). "A New Perspective on Robust M-Estimation: Finite Sample Theory and Applications to Dependence-Adjusted Multiple Testing," Annals of Statistics, to appear, \url{https://arxiv.org/abs/1711.05381}.
#' @name FarmTest
#' @docType package
NULL
