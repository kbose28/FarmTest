#' @useDynLib FarmTest
#' @importFrom  Rcpp sourceCpp
#' @importFrom  graphics mtext plot points axis par barplot
#' @import methods
#' @import utils
#' @import grDevices
#' @import stats

NULL
###################################################################################
## This is the main function that condusts the statistical test given the data
###################################################################################
#' Main function performing factor-adjusted robust test for means
#'
#' This function is used to conduct robust statistical test for means of multivariate data, after adjusting for known or unknown latent factors.
#' It uses the Huber's loss function (Huber (1964)) to robustly estimate data parameters.
#' @param X a n x p data matrix with each row being a sample.
#' You wish to test a hypothesis for the mean of each column of \code{X}.
#' @param H0 an \emph{optional} p x 1 vector of the true value of the means (or difference in means if you are performing a two sample test). The default is the zero.
#' @param fx an \emph{optional} factor matrix with each column being a factor for \code{X}. Same number of rows as \code{X}.
#' @param Kx a \emph{optional} number of factors to be estimated for \code{X}. Otherwise estimated internally.
#' @param Y an \emph{optional} data matrix that must have the same number of columns as \code{X}. You wish test the equality of means of each columns of \code{X} and \code{Y}.
#' @param fy an \emph{optional} factor matrix with each column being a factor for \code{Y}.  Same number of rows as \code{Y}. Only used for a two sample test.
#' @param Ky a \emph{optional} number of factors to be estimated for \code{Y}. Otherwise estimated internally.
#' @param alternative	an \emph{optional} character string specifying the alternative hypothesis, must be one of "two.sided" (default), "greater" or "less". You can specify just the initial letter.
#' @param alpha an \emph{optional} level for controlling the false discovery rate (in decimals). Default is 0.05. Must be in \eqn{(0,1)}.
#' @param verbose a TRUE/FALSE indicating whether to print summary of the run to console. Default is TRUE.
#' @param \dots Arguments passed to the \code{\link{farm.FDR}} function.
#' @return A list with the following items
#' \item{means}{the vector of estimated means}
#' \item{stderr}{the p x 1 vector of estimated standard errors}
#' \item{pvalue}{the p x 1 vector of unadjusted p values}
#' \item{rejected}{the indices of rejected hypotheses, along with their corresponding p values, and adjusted p values, ordered from most significant to least significant}
#' \item{alldata}{all the indices of the tested hypotheses, along with their corresponding p values, adjusted p values, and a column with 1 if declared siginificant and 0 if not}
#' \item{loadings}{estimated factor loadings}
#' \item{nfactors}{if needed, the number of estimated factors}
#' @details
#' \code{alternative = "greater"} is the alternative that \code{X} has a larger mean than \code{Y}.
#' @details
#' If some of the underlying factors are known but it is suspected that there are more confounding factors that are unobserved: Suppose we have data \eqn{X = \mu + Bf + Cg + u}, where \eqn{f} is observed and \eqn{g} is unobserved. In the first step, the user passes the data \eqn{\{X,f\}} into the main function. From the output, let us construct the residuals: \eqn{Xres = X - Bf}. Now pass \eqn{Xres} into the main function, without any factors. The output in this step is the final answer to the testing problem.
#' @details
#' Number of rows and columns of the data matrix must be at least 4 in order to be able to calculate latent factors.
#' @details For details about multiple comparison correction, see \code{\link{farm.FDR}}.
#' @examples
#' set.seed(100)
#' p = 20
#' n = 10
#' epsilon = matrix(rnorm( p*n, 0,1), nrow = n)
#' B = matrix(rnorm(p,0,1), nrow=p)
#' fx = matrix(rnorm(n, 0,1), nrow = n)
#' mu = rep(0, p)
#' mu[1:5] = 2
#' X = rep(1,n)%*%t(mu)+fx%*%t(B)+ epsilon
#' output1 = farm.test(X)
#' output = farm.test(X, alpha = 0.01,alternative = "greater")
#'
#' @references Huber, P.J. (1964). "Robust Estimation of a Location Parameter." The Annals of Mathematical Statistics, 35, 73–101.
#' @export
farm.test <- function (X, H0=NULL, fx=NULL,Kx = NULL, Y =NULL , fy=NULL, Ky  =NULL,  alternative = c("two.sided", "less", "greater"),  alpha=NULL ,verbose=TRUE, ...){
  p = NCOL(X)
  H0 <- if(is.null(H0)) rep(0,p ) else H0
  if(length(H0)!=p) stop('number of hypotheses should be the same as dimension of the data')
  alpha <- if(is.null(alpha)) 0.05 else alpha
  if(alpha>=1 || alpha <=0) stop('alpha should be between 0 and 1')

  if (!is.null(fx)){
    if(NROW(fx)!=NROW(X)) stop('number of rows in factor matrix should be the same as data matrix')
    output = farm.test.known (X, H0, fx,Y = Y , fy= fy,  alternative  = c("two.sided", "less", "greater"), alpha=alpha,...)
    if(verbose){output.call = match.call()
    cat("Call:\n")
    print(output.call)
    if(is.null(Y)){
      cat("\n One Sample Robust Test with Known Factors\n")
      cat(paste("\np = ", NCOL(X),", n = ", NROW(X), ", nfactors = ", NCOL(fx), "\n", sep = ""))
    }else
      {if(NCOL(X)!=NCOL(Y)) stop('number of rows in both data matrices must be the same')
        if(is.null(fy)) stop('must provide factors for either both or neither data matrices')
        if(NROW(fy)!=NROW(Y)) stop('number of rows in factor matrix should be the same as data matrix')
        cat("\n Two Sample Robust Test with Known Factors\n")
          cat(paste("\np = ", NCOL(X),", n.X = ", NROW(X),", n.Y = ", NROW(Y), ", X.nfactors = ", NCOL(fx), ", Y.nfactors = ", NCOL(fy), "\n", sep = ""))
      }
    cat(paste("FDR to be controlled at: ", alpha, "\n", sep = ""))
    cat(paste("alternative hypothesis: ",  match.arg(alternative), "\n",sep = ""))
    cat("hypotheses rejected:\n")
    if(is.character(output$rejected)){ cat(" no hypotheses rejected\n")} else{ cat(paste(" ", NROW(output$rejected),"\n", sep = ""))}
    }
  }
  else{
    output=farm.test.unknown (X, H0, Kx = Kx, Y=Y, Ky = Ky, alternative= c("two.sided", "less", "greater"), alpha=alpha, ...)
    if(verbose){output.call = match.call()
    cat("Call:\n")
    print(output.call)
    if(is.null(Y)){
      cat("\n One Sample Robust Test with Unknown Factors\n")
      cat(paste("\np = ", NCOL(X),", n = ", NROW(X), ", nfactors = ", output$nfactors, "\n", sep = ""))
    }else
      {if(NCOL(X)!=NCOL(Y)) stop('number of rows in both data matrices must be the same')
        cat("\n Two Sample Robust Test with Unknown Factors\n")
       cat(paste("\np = ", NCOL(X),", n.X = ", NROW(X),", n.Y = ", NROW(Y), ", X.nfactors = ", output$nfactors$X.nfactors,", Y.nfactors = ", output$nfactors$Y.nfactors,  "\n", sep = ""))}
    cat(paste("FDR to be controlled at: ",  alpha, "\n", sep = ""))
    cat(paste("alternative hypothesis: ",  match.arg(alternative),"\n", sep = ""))
    cat("hypotheses rejected:\n")
    if(is.character(output$rejected)){ cat(" no hypotheses rejected\n")} else{ cat(paste(" ", NROW(output$rejected),"\n", sep = ""))}
    }
  }
  return(output)
  }

###################################################################################
## main function (KNOWN FACTORS)
###################################################################################
farm.test.known <- function (X, H0, fx,Y  , fy,  alternative = c("two.sided", "less", "greater") , alpha,  ...){
  alternative <- match.arg(alternative)
  X = t(X)
  nx = NCOL(X)
  p = NROW(X)
  Zx <- cbind(matrix(1, nx, 1) , fx)
  Kx = NCOL(Zx)
  coefx = mu_robust_F(0.5,matrix(X, p, nx), matrix(Zx, nx, Kx))
  muhatx = coefx[1,]
  bhatx = coefx[-1,]
  thetax = mu_robust(0.5,matrix(X^2, p, nx))
  varhatx=NULL
  if(is.null( dim(bhatx))){for (j in 1:p){rvar = (thetax[j] - muhatx[j]^2)*(thetax[j] >muhatx[j]^2)+ (thetax[j])*(thetax[j]<=muhatx[j]^2)
      varf = bhatx[j] %*% stats::cov(fx) %*% bhatx[j]
      varhatx[j] = (rvar - varf)*((rvar - varf)>0) + rvar*((rvar- varf)<=0)
    }
  }else {for (j in 1:p){rvar = (thetax[j] - muhatx[j]^2)*(thetax[j] >muhatx[j]^2)+ (thetax[j])*(thetax[j]<=muhatx[j]^2)
      varf = bhatx[,j] %*% stats::cov(fx) %*% bhatx[,j]
      varhatx[j] = (rvar - varf)*((rvar - varf)>0) + rvar*((rvar- varf)<=0)
    }
  }
  sehatx = sqrt(varhatx/nx)

    if (!is.null(Y)){
        Y = t(Y)
        ny = NCOL(Y)
        Zy <- cbind(matrix(1, ny, 1) , fy)
        Ky = NCOL(Zy)
        coefy= mu_robust_F(0.5,matrix(Y, p, ny), matrix(Zy, ny, Ky))
        muhaty = coefy[1,]
        bhaty = coefy[-1,]
        thetay = mu_robust(0.5,matrix(Y^2, p, ny))
        varhaty=NULL
        if(is.null( dim(bhatx))){for (j in 1:p){rvar = (thetay[j] - muhaty[j]^2)*(thetay[j] >muhaty[j]^2)+ (thetay[j])*(thetay[j]<=muhaty[j]^2)
            varf = bhaty[j] %*% stats::cov(fy) %*% bhaty[j]
            varhaty[j] = (rvar - varf)*((rvar - varf)>0) + rvar*((rvar- varf)<=0)
            }
        }
        else{for (j in 1:p){rvar = (thetay[j] - muhaty[j]^2)*(thetay[j] >muhaty[j]^2)+ (thetay[j])*(thetay[j]<=muhaty[j]^2)
            varf = bhaty[,j] %*% stats::cov(fy) %*% bhaty[,j]
            varhaty[j] = (rvar - varf)*((rvar - varf)>0) + rvar*((rvar- varf)<=0)
        }
         }
        sehaty = sqrt(varhaty/ny)
       }

    if (is.null(Y)){
    stat=(muhatx-H0)/sehatx
    means <- muhatx
    stderr <- sehatx
    loadings = bhatx
    } else{stat=(muhatx-muhaty-H0)/sqrt(sehatx^2 + sehaty^2)
    means <- list(X.mean = muhatx, Y.mean = muhaty)
    stderr <- list(X.stderr= sehatx, Y.stderr= sehaty)
    loadings = list(X.loadings = bhatx, Y.loadings =bhaty)
    }

    if (alternative == "less"){
    pvalue = stats::pnorm(stat)
    }
    else if (alternative == "greater"){
    pvalue = stats::pnorm(stat, lower.tail = FALSE)
    }
    else {
    pvalue = 2*stats::pnorm(-abs(stat))
    }

  #Storey's procedure of pFDR
  rejected.alldata = farm.FDR(pvalue, alpha, ...)
  alldata = rejected.alldata$alldata
  rejected = rejected.alldata$rejected

  list(means = means ,  stderr=  stderr,loadings = loadings,
       pvalue = pvalue, rejected  =rejected, alldata   = alldata)
  }
#
# ###################################################################################
# ## main function (UNKNOWN FACTORS)
# ###################################################################################

farm.test.unknown <- function (X, H0,Kx, Y, Ky,  alternative = c("two.sided", "less", "greater"), alpha,... ){
  X = t(X)
  nx = NCOL(X)
  p = NROW(X)
  if(min(nx,p)<=4) stop('n and p must be at least 4')

  alternative <- match.arg(alternative)

  #for X
  muhatx = mu_robust(0.5, matrix(X, p, nx))#the first term is redundant, using CV
  covx = Cov_Huber(0.6,  X, muhatx)
  eigs = Eigen_Decomp( covx)
  values = eigs[,p+1]
  vectors = eigs[,1:p]
  #estimate nfactors
  values = pmax(values,0)
  ratio=c()
  Kx <- if (is.null(Kx)) {
    for(i in 1:(floor(min(nx,p)/2))){
      ratio=append(ratio, values[i+1]/values[i])}
    ratio = ratio[is.finite(ratio)]
    Kx = which.min(ratio)} else {Kx}
  if(Kx>=min(nx,p)/2) warning('Number of factors supplied is >= min(n,p)/2. May cause numerical inconsistencies')
  if(Kx>max(nx,p)) stop('Number of factors cannot be larger than n or p')

  Bx = matrix(NA, p, Kx)
  for (k in 1:Kx){
    Bx[,k] = sqrt(values[k])*vectors[,k]
  }
  Bx2 = apply(Bx,1, function(y) sum(y^2))
  thetax = mu_robust(0.5, matrix(X^2, p, nx))#the first term is redundant, using CV
  varhatx_0 = ( thetax - muhatx^2)* ( thetax > muhatx^2) +(thetax)* ( thetax <=muhatx^2)
  varhatx = (varhatx_0 - Bx2)* (varhatx_0 > Bx2) +(varhatx_0)* ( varhatx_0 <=Bx2)
  sehatx = sqrt(varhatx/nx)
  fx = mu_robust_F(0.5, matrix(rowMeans(X),1, p), matrix(Bx, p, Kx))


  #for Y
  if (!is.null(Y)){
  Y = t(Y)
  ny = NCOL(Y)
  if(min(ny,p)<=4) stop('n and p must be at least 4')
  muhaty = mu_robust(0.5, matrix(Y, p, ny))#the first term is redundant, using CV
  covy = Cov_Huber(0.6,  Y, muhaty)
  eigs = Eigen_Decomp( covy)
  values = eigs[,p+1]
  vectors = eigs[,1:p]
  #estimate nfactors
  values = pmax(values, 0)
  ratio=c()
  Ky <- if (is.null(Ky)) {
    for(i in 1:floor(min(p,ny)/2))
      ratio=append(ratio, values[i+1]/values[i])
    ratio = ratio[is.finite(ratio)]
    Ky = which.min(ratio)} else {Ky}
  if(Ky>=min(ny,p)/2) warning('Number of factors supplied is >= min(n,p)/2. May cause numerical inconsistencies')
  if(Ky>max(ny,p)) stop('Number of factors cannot be larger than n or p')

  By = matrix(NA, p, Ky)
  for (k in 1:Ky){
    By[,k] = sqrt(values[k])*vectors[,k]
  }
  By2 = apply(By,1, function(y) sum(y^2))
  thetay = mu_robust(0.5, matrix(Y^2, p, ny))#the first term is redundant, using CV
  varhaty_0 = ( thetay - muhaty^2)* ( thetay > muhaty^2) +(thetay)* ( thetay <=muhaty^2)
  varhaty = (varhaty_0 - By2)* (varhaty_0 > By2) +(varhaty_0)* ( varhaty_0 <=By2)
  sehaty = sqrt(varhaty/ny)
  fy = mu_robust_F(0.5, matrix(rowMeans(Y),1, p), matrix(By, p, Ky))
  }
  #test statistics
  if (is.null(Y)){
    stat=(muhatx-Bx%*%fx-H0)/sehatx
    means <- muhatx - Bx%*%fx
    stderr <- sehatx
    loadings <- Bx
    nfactors = Kx
  } else{
    stat=(muhatx - Bx%*%fx-muhaty+By%*%fy -H0)/sqrt(sehatx^2 +sehaty^2)
    means <- list(X.mean = muhatx-Bx%*%fx, Y.mean = muhaty-By%*%fy)
    stderr <- list(X.stderr= sehatx, Y.stderr= sehaty)
    loadings = list(X.loadings = Bx, Y.loadings = By)
    nfactors = list(X.nfactors= Kx, Y.nfactors =Ky)
  }
  if (alternative == "less"){
    pvalue = stats::pnorm((stat))
  }else if (alternative == "greater"){
    pvalue = stats::pnorm((stat) ,lower.tail = FALSE)
  } else {
    pvalue = 2*stats::pnorm(-abs(stat))
  }
  #Storey's procedure of pFDR
  rejected.alldata = farm.FDR(pvalue, alpha, ...)
  alldata = rejected.alldata$alldata
  rejected = rejected.alldata$rejected

  list(means = means ,  stderr=  stderr,loadings = loadings , nfactors= nfactors,
       pvalue = pvalue, rejected  =rejected, alldata = alldata)
}


# ################# #################
#' Diagnostic plots and quantities arising from estimating the number of underlying factors
#'
#' Given the data, this function  draws a scree plot and a plot of the eigenvalue ratios.
#' The eignevalue ratio test is used to estimate the number of factors. See Ahn and Horenstein(2013).
#' @param X an n x p data matrix with each row being a sample.
#' @param K.scree an \emph{optional} integer specifying the number of eigenvalues to be plotted in the scree plot. Default is min(n,p).
#' @param K.factors an \emph{optional} integer specifying the number of eigenvalues to be used for the eigenvalue ratio test. Default is min(n,p)/2.
#' @param robust a TRUE/FALSE indicating whether to use a robust covariance estimator if TRUE, or the sample covariance estimator. Default is FALSE.
#' @details The maximum eigenvalue ratio is marked differently on the plot.  The index of this maximum ratio gives the number of estimated factors.
#' @details User has to hit <Return> to see the second plot.
#' @details All the data used in the plots are output as a list.
#' @return
#' Two plots: First plot is the scree plot of the data. Second plot illustrates the eigenvalue ratio test.
#' @return A list with the data used for the plots:
#' \itemize{
#'  \item{\code{eigenvalues} }{Eigenvalues of the covariance matrix}
#'  \item{\code{proportions} }{Proportion of variance explained by the principal components}
#'  \item{\code{eigenvalue.ratios} }{Ratios calculated in the eigenvalue ratio test}
#'  \item{\code{nfactors} }{Number of factors found using the eigenvalue ratio test}
#' }
#' @examples
#' set.seed(100)
#' p = 100
#' n = 20
#' epsilon = matrix(rnorm( p*n, 0,1), nrow = n)
#' B = matrix(rnorm(p*3,0,1), nrow=p)
#' fx = matrix(rnorm(3*n, 0,1), nrow = n)
#' X = fx%*%t(B)+ epsilon
#' ouput = farm.scree(X)
#'
#' @references Ahn, S. C. and Horenstein, A. R.  (2013). "Eigenvalue Ratio Test for the Number of Factors," Econometrica, 81 (3), 1203–1227.
#' @export
farm.scree<- function(X, K.scree = NULL , K.factors = NULL , robust = FALSE){
  X = t(X)
  n = NCOL(X)
  p = NROW(X)
  if(min(n,p) <=3) stop('n and p must be at least 3')
  K.scree <- if (is.null(K.scree)) min(n,p) else K.scree
  K.factors <- if (is.null(K.factors)) (min(n,p)/2) else K.factors
  if(K.factors>min(n,p)/2) warning('Number of factors supplied is > min(n,p)/2. May cause numerical inconsistencies')
  if (robust){
    muhatx = mu_robust(0.5, matrix(X, p, n))#the first term is redundant, using CV
    covx = Cov_Huber(0.6,  X, muhatx)
    decomp = Eigen_Decomp(covx)
   eig = decomp[,p+1]
   eig = pmax(eig, 0)
  }
  else {
  pca_fit=stats::prcomp((X), center = TRUE, scale = TRUE)
  eig = (pca_fit$sdev)^2
  }
  props = eig / sum(eig)
  graphics::par(mfrow=c(1,1), mex=0.5,oma=c(0,0,4,0),mar=c(6,7,5,6))

  #plot first n eigenvalues
  grid =seq(1,K.scree)
  graphics::barplot(props[1:K.scree], main="Scree plot of the data",
          xlab=paste("Top", K.scree ,"principle components", sep=" "), ylab="Proportion of variance explained", lwd = 2, cex.lab=1, cex.axis=1, cex.main=1)
  graphics::par(new=T)
  graphics::plot(grid, eig[1:K.scree], type="b", pch=19, axes = FALSE, xlab="", ylab="")
  graphics::axis(1, at=pretty(range(grid)))
  graphics::axis(side = 4, at = pretty(range(eig[1:K.scree])))
  graphics::mtext("Eigenvalues",side=4,col="black",line=3)
  oldpar = graphics::par(ask = TRUE)
  on.exit(graphics::par(oldpar))
  #eigenvalue ratio test
  ratio=c()
  for(i in 1:K.factors)
    ratio=append(ratio, eig[i]/eig[i+1])
  ratio = ratio[is.finite(ratio)]
    graphics::plot(ratio,type="b",pch=19, ylim=c(min(ratio),max(ratio)),main = paste("Eigenvalue ratio plot:\n", which.max(ratio), "factor(s) found"),cex.main=1)
  graphics::points(x = which.max(ratio), max(ratio), col = "red",bg = "red", pch = 23,cex = 1)
  list(eigenvalues = eig, proportions = props ,  eigenvalue.ratios=  ratio, nfactors = which.max(ratio))
  }

# ################# rejections using storeys method#################
#' Control FDR given a list of pvalues
#'
#' Given a list of p-values, this function conducts multiple testing and outputs the indices of the rejected hypothesis. Uses an adaptive Benjamini-Hochberg (BH) procedure where the proportion of true nulls \eqn{\pi_0} is estimated.
#' This estimation is done based on the \href{https://www.rdocumentation.org/packages/qvalue/versions/2.4.2/topics/pi0est}{pi0est} function in the \href{http://bioconductor.org/packages/release/bioc/html/qvalue.html}{qvalue} package. See Storey(2015).
#' @param pvalue a vector of p-values obtained from multiple testing
#' @param alpha an \emph{optional} significance level for testing (in decimals). Default is 0.05. Must be in \eqn{(0,1)}.
#' @param type an \emph{optional} character string specifying the type of test. The default is the modified BH procedure (type = "mBH"). The usual BH procedure is also available (type = "BH"). See Benjamini and Hochberg (1995).
#' @param lambda an \emph{optional} threshold for estimating the proportion of true null hypotheses \eqn{\pi_0}. Must be in \eqn{[0,1)}.
#' @param pi0.method \emph{optional}, either "smoother" or "bootstrap"; the method for automatically choosing tuning parameter in the estimation of \eqn{\pi_0}, the proportion of true null hypotheses.
#' @param smooth.df an \emph{optional} number of degrees-of-freedom to use when estimating \eqn{\pi_0} with a smoother.
#' @param smooth.log.pi0 an \emph{optional} TRUE/FALSE. If TRUE and pi0.method = "smoother", \eqn{\pi_0} will be estimated by applying a smoother to a scatterplot of \eqn{\log(\pi_0)} estimates against the tuning parameter lambda. Default is FALSE.
#' @return
#' \item{rejected}{the indices of rejected hypotheses, along with their corresponding p values, and adjusted p values, ordered from most significant to least significant}
#' \item{alldata}{all the indices of the tested hypotheses, along with their corresponding p values, adjusted p values, and a column with 1 if declared siginificant and 0 if not}
#' @details The "mBH" procedure is simply the regular Benjamini-Hochberg pocedure, but in the rejection threshold the denominator \eqn{p} is replaced by  \eqn{\pi_0 * p}. This is a less conservative approach. See Storey (2002).
#' @examples
#' set.seed(100)
#' Y = matrix(rnorm(1000, 0, 1),10)
#' pval = apply(Y, 1, function(x) t.test(x)$p.value)
#' farm.FDR(pval, 0.05)
#' farm.FDR(pval, 0.01, type = "BH")
#'
#' @references Benjamini, Y. and Hochberg, Y. (1995). "Controlling the False Discovery Rate: A Practical and PowerfulApproach to Multiple Testing." Journal of the Royal Statistical Society B, 51, 289–300.
#' @references Storey, J.D. (2015). "qvalue: Q-value estimation for false discovery rate control. R package version 2.8.0, \url{https://github.com/jdstorey/qvalue}.
#' @references Storey, J.D. (2002). " Direct Approach to False Discovery Rates." Journal of the Royal Statistical Society B, 64(3), 479–498.

#' @export
farm.FDR<- function(pvalue, alpha= NULL , type = c("mBH", "BH"),lambda = seq(0.05,0.95,0.05), pi0.method = c("smoother", "bootstrap"),
smooth.df = 3, smooth.log.pi0 = FALSE){
  if((sum(pvalue<0)!=0)|| (sum(pvalue>1)!=0)) stop("pvalues must be between 0 and 1")
  type = match.arg(type)
  alpha <- if (is.null(alpha)) 0.05 else alpha
  if(alpha <=0|| alpha>=1) stop("alpha must be between 0 and 1")
  if (type == "BH"){
    p = length(pvalue)
    index = order(pvalue)
    pvalue_adj = stats::p.adjust(pvalue, method = "BH")
    alldata  = cbind(1:p, pvalue, pvalue_adj)
  }
  else {
    hatp0 = mypi0est(pvalue, lambda)[1]
    p = length(pvalue)
    pvalue_sort = sort(pvalue)
    index = order(pvalue)
    pvalue_adj = pmin(as.numeric(hatp0)* p*pvalue_sort/c(1:p), 1)
    alldata  = cbind(index, pvalue_sort, pvalue_adj)
    alldata = alldata[order(alldata[,1]),]
  }
  if (min(pvalue_adj-alpha)>0){
    rejected = "no hypotheses rejected"
    alldata = cbind( alldata, rep(0,p))
    colnames(alldata) = c("index", "pvalue", "pvalue adjusted", "significant?")
  } else{
    selected = index[1:max(which(pvalue_adj<=alpha))]
    rejected = cbind(selected, pvalue[selected], pvalue_adj[1:max(which(pvalue_adj<=alpha))])
    alldata = cbind( alldata, rep(0,p))
    alldata[selected,4] = 1
    colnames(rejected)= c("index", "pvalue", "pvalue adjusted")
    colnames(alldata) = c("index", "pvalue", "pvalue adjusted", "significant?")

  }

 list(rejected = rejected , alldata = alldata)
}

mypi0est <- function(p, lambda = seq(0.05,0.95,0.05), pi0.method = c("smoother", "bootstrap"),
                                       smooth.df = 3, smooth.log.pi0 = FALSE, ...) {
                      # Check input arguments
                      rm_na <- !is.na(p)
                      p <- p[rm_na]
                      pi0.method = match.arg(pi0.method)
                      m <- length(p)
                      lambda <- sort(lambda) # guard against user input

                      ll <- length(lambda)
                      if (min(p) < 0 || max(p) > 1) {
                        stop("ERROR: p-values not in valid range [0, 1].")
                      } else if (ll > 1 && ll < 4) {
                        stop(cat("ERROR:", paste("length(lambda)=", ll, ".", sep=""),
                                 "If length of lambda greater than 1,",
                                 "you need at least 4 values."))
                      } else if (min(lambda) < 0 || max(lambda) >= 1) {
                        stop("ERROR: Lambda must be within [0, 1).")
                      }
                      # Determines pi0
                      if (ll == 1) {
                        pi0 <- mean(p >= lambda)/(1 - lambda)
                        pi0.lambda <- pi0
                        pi0 <- min(pi0, 1)
                        pi0Smooth <- NULL
                      } else {
                        pi0 <- sapply(lambda, function(l) mean(p >= l) / (1 - l))
                        pi0.lambda <- pi0
                        # Smoother method approximation
                        if (pi0.method == "smoother") {
                          if (smooth.log.pi0) {
                            pi0 <- log(pi0)
                            spi0 <- stats::smooth.spline(lambda, pi0, df = smooth.df)
                            pi0Smooth <- exp(stats::predict(spi0, x = lambda)$y)
                            pi0 <- min(pi0Smooth[ll], 1)
                          } else {
                            spi0 <- stats::smooth.spline(lambda, pi0, df = smooth.df)
                            pi0Smooth <- stats::predict(spi0, x = lambda)$y
                            pi0 <- min(pi0Smooth[ll], 1)
                          }
                        } else if (pi0.method == "bootstrap") {
                          # Bootstrap method closed form solution by David Robinson
                          minpi0 <- stats::quantile(pi0, prob = 0.1)
                          W <- sapply(lambda, function(l) sum(p >= l))
                          mse <- (W / (m ^ 2 * (1 - lambda) ^ 2)) * (1 - W / m) + (pi0 - minpi0) ^ 2
                          pi0 <- min(pi0[mse == min(mse)], 1)
                          pi0Smooth <- NULL
                        } else {
                          stop('ERROR: pi0.method must be one of "smoother" or "bootstrap".')
                        }
                      }
                      if (pi0 <= 0) {
                        stop("ERROR: The estimated pi0 <= 0. Check that you have valid p-values or use a different range of lambda. Alternatively, to not estimate pi0, use type = \"BH\".")
                      }
                      return(list(pi0 = pi0, pi0.lambda = pi0.lambda,
                                  lambda = lambda, pi0.smooth = pi0Smooth))
      }
