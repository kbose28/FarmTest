#' @useDynLib FarmTest
#' @importFrom  Rcpp sourceCpp
#' @importFrom  graphics mtext plot points axis par barplot
#' @import methods
#' @import utils
#' @import grDevices
#' @import stats
NULL
###################################################################################
## This is the main function that conducts the statistical test given the data
###################################################################################
#' Main function performing factor-adjusted robust test for means
#'
#' This function is used to conduct robust statistical test for means of multivariate data, after adjusting for known or unknown latent factors using the methods in Fan et al.(2017) and Zhou et al.(2017).
#' It uses the Huber's loss function (Huber (1964)) to robustly estimate data parameters.
#' @param X a n x p data matrix with each row being a sample.
#' You wish to test a hypothesis for the mean of each column of \code{X}.
#' @param H0 an \emph{optional} p x 1 vector of the true value of the means (or difference in means if you are performing a two sample test). The default is the zero.
#' @param fx an \emph{optional} factor matrix with each column being a factor for \code{X}. Same number of rows as \code{X}.
#' @param Kx a \emph{optional} number of factors to be estimated for \code{X}. Otherwise estimated internally. Kx>=0
#' @param Y an \emph{optional} data matrix that must have the same number of columns as \code{X}. You wish test the equality of means of each columns of \code{X} and \code{Y}.
#' @param fy an \emph{optional} factor matrix with each column being a factor for \code{Y}.  Same number of rows as \code{Y}. Only used for a two sample test.
#' @param Ky a \emph{optional} number of factors to be estimated for \code{Y}. Otherwise estimated internally.
#' @param alternative	an \emph{optional} character string specifying the alternate hypothesis, must be one of "two.sided" (default), "greater" or "lesser". You can specify just the initial letter.
#' @param alpha an \emph{optional} level for controlling the false discovery rate (in decimals). Default is 0.05. Must be in \eqn{(0,1)}.
#' @param robust a boolean, specifying whether or not to use robust estimators for mean and variance. Default is TRUE.
#' @param cv a boolean, specifying whether or  not to run cross-validation for the tuning parameter. Default is TRUE. Only used if \code{robust} is TRUE.
#' @param tau \code{>0}, multiplier for the tuning parameter for Huber loss function. Default is 2. Only used if \code{robust} is TRUE and \code{cv} is FALSE. See details.
#' @param verbose a boolean specifying whether to print runtime updates to the console. Default is TRUE.
#' @param \dots Arguments passed to the \code{\link{farm.FDR}} function.
#' @return An object with S3 class \code{farm.test} containing:
#' \item{means}{estimated means}
#' \item{stderr}{estimated standard errors}
#' \item{pvalue}{unadjusted p values}
#' \item{rejected}{the indices of rejected hypotheses, along with their corresponding p values, and adjusted p values, ordered from most significant to least significant}
#' \item{alldata}{all the indices of the tested hypotheses, along with their corresponding p values, adjusted p values, and a column with 1 if declared siginificant and 0 if not}
#' \item{loadings}{estimated factor loadings}
#' \item{nfactors}{the number of (estimated) factors}
#' \item{significant}{the number of means that are found significant}
#' \item{\dots}{further arguments passed to methods. For complete list use the function \code{\link{names}} on the output object}
#' @details
#' \code{alternative = "greater"} is the alternative that \code{X} has a larger mean than \code{Y}.
#' @details If some of the underlying factors are known but it is suspected that there are more confounding factors that are unobserved: Suppose we have data \eqn{X = \mu + Bf + Cg + u}, where \eqn{f} is observed and \eqn{g} is unobserved. In the first step, the user passes the data \eqn{\{X,f\}} into the main function. From the output, let us construct the residuals: \eqn{Xres = X - Bf}. Now pass \eqn{Xres} into the main function, without any factors. The output in this step is the final answer to the testing problem.
#' @details For two-sample test, the output values \code{means}, \code{stderr}, \code{n}, \code{nfactors},\code{loadings} are all lists containing two items, each pertaining to \code{X} and \code{Y}, indicated by a prefix \code{X.} and \code{Y.} respectively.
#' @details Number of rows and columns of the data matrix must be at least 4 in order to be able to calculate latent factors.
#' @details For details about multiple comparison correction, see \code{\link{farm.FDR}}.
#' @details The tuning parameter \code{= tau *  sigma * optimal rate } where \code{optimal rate } is the optimal rate for the tuning parameter. For details, see Fan et al.(2017). \code{sigma} is the standard deviation of the data.
#' @seealso \code{\link{farm.FDR}}, \code{\link{print.farm.test}}
#' @examples
#' set.seed(100)
#' p = 100
#' n = 50
#' epsilon = matrix(rnorm( p*n, 0,1), nrow = n)
#' B = matrix(runif(p*3,-2,2), nrow=p)
#' fx = matrix(rnorm(3*n, 0,1), nrow = n)
#' mu = rep(0, p)
#' mu[1:5] = 2
#' X = rep(1,n)%*%t(mu)+fx%*%t(B)+ epsilon
#' output = farm.test(X, cv=FALSE)#robust, no cross-validation
#' output
#'
#' #other robustification options
#' output = farm.test(X, robust = FALSE, verbose=FALSE) #non-robust
#' output = farm.test(X, tau = 3, cv=FALSE, verbose=FALSE) #robust, no cross-validation, specified tau
#' #output = farm.test(X) #robust, cross-validation, longer running
#'
#' #two sample test
#' n2 = 25
#' epsilon = matrix(rnorm( p*n2, 0,1), nrow = n2)
#' B = matrix(rnorm(p*3,0,1), nrow=p)
#' fy = matrix(rnorm(3*n2, 0,1), nrow = n2)
#' Y = fy%*%t(B)+ epsilon
#' output = farm.test(X=X,Y=Y, robust=FALSE)
#' names(output$means)
#'
#' @references Huber, P.J. (1964). "Robust Estimation of a Location Parameter." The Annals of Mathematical Statistics, 35, 73–101.
#' @references Fan, J., Ke, Y., Sun, Q. and Zhou, W-X. (2017). "FARM-Test: Factor-Adjusted Robust Multiple Testing with False Discovery Control", \url{https://arxiv.org/abs/1711.05386}.
#' @references Zhou, W-X., Bose, K., Fan, J. and Liu, H. (2017). "A New Perspective on Robust M-Estimation: Finite Sample Theory and Applications to Dependence-Adjusted Multiple Testing," Annals of Statistics, to appear, \url{https://arxiv.org/abs/1711.05381}.
#' @export
farm.test <- function (X, H0=NULL, fx=NULL,Kx = NULL, Y =NULL , fy=NULL, Ky  =NULL,  alternative = c("two.sided", "lesser", "greater"),  alpha=NULL ,robust=TRUE,cv = TRUE, tau=2,verbose = FALSE,...){
  p = NCOL(X)
  alternative <- match.arg(alternative)
  H0 <- if(is.null(H0)) rep(0,p ) else H0
  alpha <- if(is.null(alpha)) 0.05 else alpha
  if(tau<=0) stop('tau should be a positive number')
  #error checking
  if(length(H0)!=p) stop('number of hypotheses should be the same as dimension of the data')
  if(alpha>=1 || alpha <=0) stop('alpha should be between 0 and 1')

  if (!is.null(fx)){
    if(NROW(fx)!=NROW(X)) stop('number of rows in factor matrix should be the same as data matrix')
    if(!is.null(Y)){
      if(NCOL(X)!=NCOL(Y)) stop('number of rows in both data matrices must be the same')
      if(is.null(fy)) stop('must provide factors for either both or neither data matrices')
      if(NROW(fy)!=NROW(Y)) stop('number of rows in factor matrix should be the same as data matrix')
      }
    }else{
    if(!is.null(Y)){
      if(!is.null(fy)) stop('must provide factors for either both or neither data matrices')
      if(NCOL(X)!=NCOL(Y)) stop('number of rows in both data matrices must be the same')
      }
    }
  #call main function
  if(!is.null(fx)){output = farm.testknown (X, H0, fx,Y = Y , fy= fy,  alternative  =alternative, alpha=alpha,robust, cv, tau,verbose=verbose,...)}
  else{output=farm.testunknown (X, H0, Kx = Kx, Y=Y, Ky = Ky, alternative=alternative, alpha=alpha,robust, cv, tau,verbose=verbose,...)}
  value = (output)
  attr(value, "class") <- "farm.test"
  value
}
###################################################################################
## Print method
###################################################################################
#' Summarize and print the results of the multiple testing
#'
#' Print method for \code{farm.test} objects
#' @param x A \code{farm.test} object.
#' @param \dots Further arguments passed to or from other methods.
#' @return A list with the following items:
#' \item{means}{estimated means}
#' \item{stderr}{estimated standard errors}
#' \item{pvalue}{unadjusted p values}
#' \item{rejected}{the indices of rejected hypotheses, along with their corresponding p values, and adjusted p values, ordered from most significant to least significant}
#' \item{alldata}{all the indices of the tested hypotheses, along with their corresponding p values, adjusted p values, and a column with 1 if declared siginificant and 0 if not}
#' \item{loadings}{estimated factor loadings}
#' \item{nfactors}{number of (estimated) factors}
#' \item{n}{number of observations}
#' \item{p}{number of dimensions}
#' \item{alpha}{level at which FDR was controlled}
#' \item{H0}{null hypothesis}
#' \item{alternative}{alternate hypothesis}
#' \item{robust}{whether robust parameters were used}
#' \item{type}{whether the test is one or two-sided}
#' \item{significant}{the number of means that are found significant}
#' @seealso \code{\link{farm.test}}
#' @examples
#' set.seed(100)
#' p = 50
#' n = 100
#' X = matrix(rnorm( p*n, 0,1), nrow = n)
#' output = farm.test(X)
#' output
#' names(output)
#' @export
print.farm.test<-function(x,...){
  if (x$type=="known"){
    if(length(x$n)==1){
     cat(paste("\n One Sample",if(x$robust) "Robust", "Test with Known Factors\n"))
     cat(paste("\np = ", x$p,", n = ", x$n, ", nfactors = ", x$nfactors, "\n", sep = ""))
   }else{
    cat(paste("\n Two Sample",if(x$robust) "Robust", "Test with Known Factors\n"))
    cat(paste("\np = ", x$p,", nX = ", x$n$X.n,", nY = ", x$n$Y.n, ", X.nfactors = ",x$nfactors$X.nfactors, ", Y.nfactors = ", x$nfactors$Y.nfactors, "\n", sep = ""))
  }
  cat(paste("FDR to be controlled at: ", x$alpha, "\n", sep = ""))
  cat(paste("alternative hypothesis: ",  x$alternative, "\n",sep = ""))
  cat("hypotheses rejected:\n")
  if(x$significant==0){ cat(" no hypotheses rejected\n")} else{ cat(paste(" ", x$significant,"\n", sep = ""))}
  }else{
  if(length(x$n)==1){
    cat(paste("\n One Sample",if(x$robust) "Robust", "Test with Unknown Factors\n"))
    cat(paste("\np = ", x$p,", n = ", x$n, ", nfactors = ", x$nfactors, "\n", sep = ""))
  }else{cat(paste("\n Two Sample",if(x$robust) "Robust", "Test with Unknown Factors\n"))
  cat(paste("\np = ", x$p,", nX = ",x$n$X.n,", nY = ", x$n$Y.n, ", X.nfactors = ", x$nfactors$X.nfactors,", Y.nfactors = ", x$nfactors$Y.nfactors,  "\n", sep = ""))
  }
  cat(paste("FDR to be controlled at: ", x$alpha, "\n", sep = ""))
  cat(paste("alternative hypothesis: ",  x$alternative, "\n",sep = ""))
  cat("hypotheses rejected:\n")
  if(x$significant==0){ cat(" no hypotheses rejected\n")} else{ cat(paste(" ", x$significant,"\n", sep = ""))}
  }
}

###################################################################################
## Plot method for farm.scree
###################################################################################
#' Diagnostic plots from factor-finding
#'
#' Plot method for \code{farm.scree} objects. Plots the eigenvalue ratio plot and the scree plot.
#' @param x A "\code{farm.scree}" object.
#' @param scree.plot \emph{optional} indicating whether to show the scree plot. Default TRUE
#' @param ratio.plot \emph{optional} indicating whether to show the scree plot. Default TRUE.
#' @param col Controls the color of the maximim eigenvalue dot. Defaut "red".
#' @param \dots graphical parameters to \code{\link{plot}}.
#' @return Two plots: First plot is the scree plot of the data. Second plot illustrates the eigenvalue ratio test.
#' @details By default, two plots are output with default options. To customize plots, plot one at a time and customize.
#' @seealso \code{\link{farm.scree}} and \code{\link{print.farm.scree}}
#' @examples
#' set.seed(100)
#' p = 100
#' n = 20
#' epsilon = matrix(rnorm( p*n, 0,1), nrow = n)
#' B = matrix(rnorm(p*3,0,1), nrow=p)
#' fx = matrix(rnorm(3*n, 0,1), nrow = n)
#' X = fx%*%t(B)+ epsilon
#' output = farm.scree(X)
#' plot(output)
#' plot(output, scree.plot=FALSE, col="blue", main="Customized plot")
#'
#' @export
plot.farm.scree<-function(x, scree.plot=TRUE, ratio.plot=TRUE, col="red", ...){


  col <- if(is.null(col))"red" else col
  new.args = list(...)
  if(scree.plot==TRUE &ratio.plot==TRUE){
    graphics::par(mfrow=c(1,1), mex=0.5,oma=c(0,0,4,0),mar=c(6,7,5,6))
    #plot first n eigenvalues
    grid =seq(1,x$K.scree)
    graphics::barplot(x$proportions[1:x$K.scree], main="Scree plot of the data",
                      xlab=paste("Top", x$K.scree ,"principle components", sep=" "), ylab="Proportion of variance explained", lwd = 2, cex.lab=1, cex.axis=1, cex.main=1)
    graphics::par(new=T)
    graphics::plot(grid, x$eigenvalues[1:x$K.scree], type="b", pch=19, axes = FALSE, xlab="", ylab="")
    graphics::axis(1, at=pretty(range(grid)))
    graphics::axis(side = 4, at = pretty(range(x$eigenvalues[1:x$K.scree])))
    graphics::mtext("Eigenvalues",side=4,col="black",line=3)
    oldpar = graphics::par(ask = TRUE)
    on.exit(graphics::par(oldpar))
    graphics::plot(x$eigenvalue.ratios,type="b",pch=19,ylab="ratio", ylim=c(min(x$eigenvalue.ratios),max(x$eigenvalue.ratios)),main = paste("Eigenvalue ratio plot:\n", x$nfactors, "factor(s) found"),cex.main=1)
    graphics::points(x = x$nfactors, max(x$eigenvalue.ratios), col = col,bg = col, pch = 23,cex = 1)
  }else if(scree.plot == TRUE & ratio.plot==FALSE){
    plot.args = list(x=seq(1,x$K.scree), y=x$proportions[1:x$K.scree], main="Scree plot of the data", xlab=paste("Top", length(x$eigenvalues) ,"principle components", sep=" "),ylab="Proportion of variance explained", lwd = 2, cex.lab=1, cex.axis=1, cex.main=1)
    if (length(new.args)) plot.args[names(new.args)] = new.args
    do.call("plot", plot.args)
    }else if(scree.plot == FALSE & ratio.plot==TRUE){
    plot.args = list(x=seq(1,length(x$eigenvalue.ratios)),y= x$eigenvalue.ratios, ylim=c(min(x$eigenvalue.ratios),max(x$eigenvalue.ratios)), xlab="Index",type="b", main= paste("Eigenvalue ratio plot:\n", x$nfactors, "factor(s) found"),pch=19,ylab="ratio", cex.main=1)
    if (length(new.args)) plot.args[names(new.args)] = new.args
    do.call("plot", plot.args)
    graphics::points(x = x$nfactors, max(x$eigenvalue.ratios), col = col,bg = col, pch = 23,cex = 1)
 }
}

###################################################################################
## Print method for farm.scree
###################################################################################
#' Summarize and print the results of the eignevalue ratio test
#'
#' Print method for farm.scree objects.
#' @param x A "farm.scree" object.
#' @param \dots Further arguments passed to or from other methods.
#' @return Summarizes the results of the factor-finding step.
#' @seealso \code{\link{farm.scree}} and \code{\link{plot.farm.scree}}
#' @examples
#' set.seed(100)
#' p = 100
#' n = 20
#' epsilon = matrix(rnorm( p*n, 0,1), nrow = n)
#' B = matrix(rnorm(p*3,0,1), nrow=p)
#' fx = matrix(rnorm(3*n, 0,1), nrow = n)
#' X = fx%*%t(B)+ epsilon
#' output = farm.scree(X)
#' output
#' @export
print.farm.scree<-function(x,...){
  cat("Summary of eigenvalue ratio test\n")
  cat(paste("  Number of factors found: ",  x$nfactors, "\n", sep = ""))
  cat(paste("  Proportion of variation explained by the top ", x$nfactors, " principal components: ", paste(round(100*sum(x$proportions[1:x$nfactors]), 2), "%", sep=""),"\n",sep = ""))
}

###################################################################################
## main function (KNOWN FACTORS)
###################################################################################
farm.testknown <- function (X, H0, fx,Y  , fy,  alternative = alternative, alpha,robust, cv, tau, verbose,...){
  X = t(X)
  nx = NCOL(X)
  p = NROW(X)
  Zx <- cbind(matrix(1, nx, 1) , fx)
  Kx = NCOL(Zx)
  if(robust==TRUE){
    if(cv==TRUE){
    coefx = mu_robust_F(matrix(X, p, nx), matrix(Zx, nx, Kx))
    muhatx = coefx[1,]
    bhatx = coefx[-1,]
    thetax = mu_robust(matrix(X^2, p, nx))
    }else{olsres = apply(X, 1, function(x) lm(x~Zx-1)$residuals)
    CT = tau*apply(olsres,2,sd)*sqrt(nx/log(p*nx))
      coefx = mu_robust_F_noCV(matrix(X, p, nx), matrix(Zx, nx, Kx), matrix(CT, p,1) )
      muhatx = coefx[1,]
      bhatx = coefx[-1,]
      CT = tau*apply(X^2,1,sd)*sqrt(nx/log(p*nx))
      thetax = mu_robust_noCV(matrix(X^2, p, nx), matrix(CT, p,1))
    }
  }else{
    coefx= apply((X), 1, function(x) lm(x~Zx-1)$coefficients)
    muhatx = coefx[1,]
    bhatx = coefx[-1,]
    thetax = rowMeans(X^2)
  }


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


    if(robust==TRUE){
      if(cv==TRUE){
        coefy = mu_robust_F(matrix(Y, p, ny), matrix(Zy, ny, Ky))
        muhaty = coefy[1,]
        bhaty = coefy[-1,]
        thetay = mu_robust(matrix(Y^2, p, ny))
      }else{olsres = apply(Y,1, function(x) lm(x~Zy-1)$residuals)
        CT = tau*apply(olsres,2,sd)*sqrt(ny/log(p*ny))
      coefy = mu_robust_F_noCV(matrix(Y, p, ny), matrix(Zy, ny, Ky), matrix(CT, p,1) )
      muhaty = coefy[1,]
      bhaty = coefy[-1,]
      CT = tau*apply(Y^2,1,sd)*sqrt(ny/log(p*ny))
      thetay = mu_robust_noCV(matrix(Y^2, p, ny), matrix(CT, p,1))
      }
    }else{
      coefy= apply((Y), 1, function(x) lm(x~Zy-1)$coefficients)
      muhaty = coefy[1,]
      bhaty = coefy[-1,]
      thetay = rowMeans(Y^2)
    }


    varhaty=NULL
    if(is.null( dim(bhaty))){for (j in 1:p){rvar = (thetay[j] - muhaty[j]^2)*(thetay[j] >muhaty[j]^2)+ (thetay[j])*(thetay[j]<=muhaty[j]^2)
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
    nfactors  = Kx-1
    n = nx
  } else{stat=(muhatx-muhaty-H0)/sqrt(sehatx^2 + sehaty^2)
  means <- list(X.mean = muhatx, Y.mean = muhaty)
  stderr <- list(X.stderr= sehatx, Y.stderr= sehaty)
  loadings = list(X.loadings = bhatx, Y.loadings =bhaty)
  nfactors = list(X.nfactors= Kx-1, Y.nfactors =Ky-1)
  n = list(X.n = nx,Y.n = ny)
  }

  if (alternative == "lesser"){
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
  significant= rejected.alldata$significant

  val<-list(means = means ,  stderr=  stderr,loadings = loadings , nfactors= nfactors,
                      pvalue = pvalue, rejected  =rejected, alldata = alldata, alternative = alternative,
                      H0 = H0, robust = robust, n = n, p = p, alpha = alpha, type = "known", significant = significant)
  return(val)

}
#
# ###################################################################################
# ## main function (UNKNOWN FACTORS)
# ###################################################################################

farm.testunknown <- function (X, H0,Kx, Y, Ky,  alternative = alternative, alpha,robust, cv, tau,verbose, ... ){
  X = t(X)
  nx = NCOL(X)
  p = NROW(X)
  if(min(nx,p)<=4) stop('n and p must be at least 4')

  #for X
  if(!is.null(Kx)){
    if(Kx < 0)stop('number of factors cannot be < 0')
    if(Kx==0){
      if(robust==TRUE){
        if(cv==TRUE){
          muhatx = mu_robust( matrix(X, p, nx))
          thetax = mu_robust( matrix(X^2, p, nx))
          }else{
            CT = tau*apply(X,1,sd)*sqrt(nx/log(p*nx))
            muhatx = mu_robust_noCV( matrix(X, p, nx), matrix(CT, p,1))
            CT = tau*apply(X^2,1,sd)*sqrt(nx/log(p*nx))
            thetax = mu_robust_noCV( matrix(X^2, p, nx), matrix(CT,p,1))
            }
        }else{
          muhatx = rowMeans(X)
          thetax = rowMeans(X^2)
          }
      varhatx = ( thetax - muhatx^2)* ( thetax > muhatx^2) +(thetax)* ( thetax <=muhatx^2)
      sehatx = sqrt(varhatx/nx)
      }else{
        if(Kx>min(nx,p)/2) warning('Number of factors supplied is >= min(n,p)/2. May cause numerical inconsistencies')
        if(Kx>max(nx,p)) stop('Number of factors cannot be larger than n or p')
        if(robust==TRUE){
          if(cv==TRUE){
            muhatx = mu_robust( matrix(X, p, nx))
            thetax = mu_robust( matrix(X^2, p, nx))
            if(verbose){cat("calculating covariance matrix for X...\n")}
            covx = Cov_Huber( X, matrix(muhatx, p, 1))
          }else{
            CT = tau*apply(X,1,sd)*sqrt(nx/log(p*nx))
            muhatx = mu_robust_noCV( matrix(X, p, nx), matrix(CT, p,1))
            CT = tau*apply(X^2,1,sd)*sqrt(nx/log(p*nx))
            thetax = mu_robust_noCV( matrix(X^2, p, nx), matrix(CT,p,1))
            if(verbose){cat("calculating covariance matrix for X...\n")}
            CT = Cov_Huber_tune(X, tau)
            covx =  Cov_Huber_noCV(matrix((X),p,nx), matrix(muhatx, p, 1), matrix(CT,p,p))
          }
        }else{
          muhatx = rowMeans(X)
          covx = cov(t(X))
          thetax = rowMeans(X^2)
        }
        eigs = Eigen_Decomp( covx)
        values = eigs[,p+1]
        vectors = eigs[,1:p]
        #estimate nfactors
        values = pmax(values,0)
        Bx = matrix(NA, p, Kx)
        for (k in 1:Kx){Bx[,k] = sqrt(values[k])*vectors[,k]}
        Bx2 = apply(Bx,1, function(y) sum(y^2))
        varhatx_0 = ( thetax - muhatx^2)* ( thetax > muhatx^2) +(thetax)* ( thetax <=muhatx^2)
        varhatx = (varhatx_0 - Bx2)* (varhatx_0 > Bx2) +(varhatx_0)* ( varhatx_0 <=Bx2)
        sehatx = sqrt(varhatx/nx)
        if(robust == TRUE){
          if(cv==TRUE){
            fx = mu_robust_F( matrix(rowMeans(X),1, p), matrix(Bx, p, Kx))
          }else{
            olsres = lm( matrix(rowMeans(X),p, 1) ~  matrix(Bx, p, Kx)-1)$residuals
            CT = tau*sd(olsres)*sqrt(p/log(nx))
            fx = mu_robust_F_noCV( matrix(rowMeans(X),1, p), matrix(Bx, p, Kx), matrix(CT, 1,1))
          }
        }else{
          fx = coef(lm(rowMeans(X)~Bx-1))
          }
        }
    }else{
      if(robust==TRUE){
        if(cv==TRUE){
          muhatx = mu_robust( matrix(X, p, nx))
          thetax = mu_robust( matrix(X^2, p, nx))
          if(verbose){cat("calculating covariance matrix for X...\n")}
          covx = Cov_Huber( X, matrix(muhatx, p, 1))
        }else{
          CT = tau*apply(X,1,sd)*sqrt(nx/log(p*nx))
          muhatx = mu_robust_noCV( matrix(X, p, nx), matrix(CT, p,1))
          CT = tau*apply(X^2,1,sd)*sqrt(nx/log(p*nx))
          thetax = mu_robust_noCV( matrix(X^2, p, nx), matrix(CT,p,1))
          if(verbose){cat("calculating covariance matrix for X...\n")}
          CT = Cov_Huber_tune(X, tau)
          covx =  Cov_Huber_noCV(matrix((X),p,nx), matrix(muhatx, p, 1), matrix(CT,p,p))
        }
      }else{
        muhatx = rowMeans(X)
        covx = cov(t(X))
        thetax = rowMeans(X^2)
      }
      eigs = Eigen_Decomp( covx)
      values = eigs[,p+1]
      vectors = eigs[,1:p]
      #estimate nfactors
      values = pmax(values,0)
      ratio=c()
      for(i in 1:(floor(min(nx,p)/2))){ratio=append(ratio, values[i+1]/values[i])}
      ratio = ratio[is.finite(ratio)]
      Kx = which.min(ratio)
      Bx = matrix(NA, p, Kx)
      for (k in 1:Kx){Bx[,k] = sqrt(values[k])*vectors[,k]}
      Bx2 = apply(Bx,1, function(y) sum(y^2))
      varhatx_0 = ( thetax - muhatx^2)* ( thetax > muhatx^2) +(thetax)* ( thetax <=muhatx^2)
      varhatx = (varhatx_0 - Bx2)* (varhatx_0 > Bx2) +(varhatx_0)* ( varhatx_0 <=Bx2)
      sehatx = sqrt(varhatx/nx)
      if(robust == TRUE){
        if(cv==TRUE){
          fx = mu_robust_F( matrix(rowMeans(X),1, p), matrix(Bx, p, Kx))
        }else{
          olsres = lm( matrix(rowMeans(X),p, 1) ~  matrix(Bx, p, Kx)-1)$residuals
          CT = tau*sd(olsres)*sqrt(p/log(nx))
          fx = mu_robust_F_noCV( matrix(rowMeans(X),1, p), matrix(Bx, p, Kx), matrix(CT, 1,1))
        }
      }else{
        fx = coef(lm(rowMeans(X)~Bx-1))
      }
    }

  #for Y
  if (!is.null(Y)){
    Y = t(Y)
    ny = NCOL(Y)
    if(min(ny,p)<=4) stop('n and p must be at least 4')
      if(!is.null(Ky)){
        if(Ky < 0)stop('number of factors cannot be < 0')
        if(Ky==0){
          if(robust==TRUE){
            if(cv==TRUE){
              muhaty = mu_robust( matrix(Y, p, ny))
              thetay = mu_robust( matrix(Y^2, p, ny))
              }else{
                CT = tau*apply(Y,1,sd)*sqrt(ny/log(p*ny))
                muhaty = mu_robust_noCV( matrix(Y, p, ny), matrix(CT, p,1))
                CT = tau*apply(Y^2,1,sd)*sqrt(ny/log(p*ny))
                thetay = mu_robust_noCV( matrix(Y^2, p, ny), matrix(CT, p,1))
              }
            }else{
              muhaty = rowMeans(Y)
              thetay = rowMeans(Y^2)
              }
          varhaty = ( thetay - muhaty^2)* ( thetay > muhaty^2) +(thetay)* ( thetay <=muhaty^2)
          sehaty = sqrt(varhaty/ny)
          }else{
            if(Ky>min(ny,p)/2) warning('Number of factors supplied is >= min(n,p)/2. May cause numerical inconsistencies')
            if(Ky>max(ny,p)) stop('Number of factors cannot be larger than n or p')
            By = matrix(NA, p, Ky)
            for (k in 1:Ky){By[,k] = sqrt(values[k])*vectors[,k]}
            By2 = apply(By,1, function(y) sum(y^2))
            varhaty_0 = ( thetay - muhaty^2)* ( thetay > muhaty^2) +(thetay)* ( thetay <=muhaty^2)
            varhaty = (varhaty_0 - By2)* (varhaty_0 > By2) +(varhaty_0)* ( varhaty_0 <=By2)
            sehaty = sqrt(varhaty/ny)
            if(robust == TRUE){
              if(cv==TRUE){
                fy = mu_robust_F( matrix(rowMeans(Y),1, p), matrix(By, p, Ky))
              }else
              {
                olsres = lm( matrix(rowMeans(Y),p, 1)~  matrix(By, p, Ky)-1)$residuals
                CT = tau*sd(olsres)*sqrt(p/log(ny))
                fy = mu_robust_F_noCV( matrix(rowMeans(Y),1, p), matrix(By, p, Ky), matrix(CT, 1,1))
              }
            }else{
              fy = coef(lm(rowMeans(Y)~By-1))
              }
            }
        }else{
              if(robust==TRUE){
               if(cv==TRUE){
              muhaty = mu_robust( matrix(Y, p, ny))
              thetay = mu_robust( matrix(Y^2, p, ny))
              if(verbose){cat("calculating covariance matrix for Y...\n")}
              covy = Cov_Huber( Y, matrix(muhaty, p, 1))
            }else{
              CT = tau*apply(Y,1,sd)*sqrt(ny/log(p*ny))
              muhaty = mu_robust_noCV( matrix(Y, p, ny), matrix(CT, p,1))
              CT = tau*apply(Y^2,1,sd)*sqrt(ny/log(p*ny))
              thetay = mu_robust_noCV( matrix(Y^2, p, ny), matrix(CT, p,1))
              if(verbose){cat("calculating covariance matrix for Y...\n")}
              CT = Cov_Huber_tune(Y, tau)
              covy =  Cov_Huber_noCV(matrix((Y),p,ny), matrix(muhaty, p, 1), matrix(CT,p,p))
            }
          }else{
            muhaty = rowMeans(Y)
            covy = cov(t(Y))
            thetay = rowMeans(Y^2)
          }
          eigs = Eigen_Decomp( covy)
          values = eigs[,p+1]
          vectors = eigs[,1:p]
          #estimate nfactors
          values = pmax(values, 0)
          ratio=c()
          for(i in 1:floor(min(p,ny)/2)){ratio=append(ratio, values[i+1]/values[i])}
          ratio = ratio[is.finite(ratio)]
          Ky = which.min(ratio)
          By = matrix(NA, p, Ky)
          for (k in 1:Ky){By[,k] = sqrt(values[k])*vectors[,k]}
          By2 = apply(By,1, function(y) sum(y^2))
          varhaty_0 = ( thetay - muhaty^2)* ( thetay > muhaty^2) +(thetay)* ( thetay <=muhaty^2)
          varhaty = (varhaty_0 - By2)* (varhaty_0 > By2) +(varhaty_0)* ( varhaty_0 <=By2)
          sehaty = sqrt(varhaty/ny)
          if(robust == TRUE){
            if(cv==TRUE){
              fy = mu_robust_F( matrix(rowMeans(Y),1, p), matrix(By, p, Ky))
            }else
            {
              olsres = lm( matrix(rowMeans(Y),p, 1)~  matrix(By, p, Ky)-1)$residuals
              CT = tau*sd(olsres)*sqrt(p/log(ny))
              fy = mu_robust_F_noCV( matrix(rowMeans(Y),1, p), matrix(By, p, Ky), matrix(CT, 1,1))
            }
          }else{
            fy = coef(lm(rowMeans(Y)~By-1))
          }
            }

  }


  #test statistics
  if (is.null(Y)){
    if(Kx==0){
      stat=(muhatx-H0)/sehatx
      means <- muhatx
      stderr <- sehatx
      loadings <- NULL
      nfactors=Kx
      n = nx
    }else{
    stat=(muhatx-Bx%*%fx-H0)/sehatx
    means <- muhatx - Bx%*%fx
    stderr <- sehatx
    loadings <- Bx
    nfactors=Kx
    n = nx}
  } else{
    if(Kx==0 & Ky!=0){
      stat=(muhatx - muhaty+By%*%fy -H0)/sqrt(sehatx^2 +sehaty^2)
      means <- list(X.mean = muhatx, Y.mean = muhaty-By%*%fy)
      stderr <- list(X.stderr= sehatx, Y.stderr= sehaty)
      loadings = list(X.loadings = NULL, Y.loadings = By)
      nfactors = list(X.nfactors= Kx, Y.nfactors =Ky)
      n = list(X.n=nx, Y.n=ny)
    } else if(Ky==0 & Kx!=0){
      stat=(muhatx - Bx%*%fx-muhaty -H0)/sqrt(sehatx^2 +sehaty^2)
      means <- list(X.mean = muhatx-Bx%*%fx, Y.mean = muhaty)
      stderr <- list(X.stderr= sehatx, Y.stderr= sehaty)
      loadings = list(X.loadings = Bx, Y.loadings = NULL)
      nfactors = list(X.nfactors= Kx, Y.nfactors =Ky)
      n = list(X.n=nx, Y.n=ny)
    }
    else if(Kx==0 & Ky==0){
      stat=(muhatx -muhaty-H0)/sqrt(sehatx^2 +sehaty^2)
      means <- list(X.mean = muhatx, Y.mean = muhaty)
      stderr <- list(X.stderr= sehatx, Y.stderr= sehaty)
      loadings = list(X.loadings =NULL, Y.loadings = NULL)
      nfactors = list(X.nfactors= Kx, Y.nfactors =Ky)
      n = list(X.n=nx, Y.n=ny)
    }else{
    stat=(muhatx - Bx%*%fx-muhaty+By%*%fy -H0)/sqrt(sehatx^2 +sehaty^2)
    means <- list(X.mean = muhatx-Bx%*%fx, Y.mean = muhaty-By%*%fy)
    stderr <- list(X.stderr= sehatx, Y.stderr= sehaty)
    loadings = list(X.loadings = Bx, Y.loadings = By)
    nfactors = list(X.nfactors= Kx, Y.nfactors =Ky)
    n = list(X.n=nx, Y.n=ny)
    }
    }
  if (alternative == "lesser"){
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
  significant= rejected.alldata$significant


  val<-list(means = means ,  stderr=  stderr,loadings = loadings , nfactors= nfactors,
                 pvalue = pvalue, rejected  =rejected, alldata = alldata, alternative = alternative,
                 H0 = H0, robust = robust, n = n, p = p, alpha = alpha, type = "unknown",significant = significant)
  return(val)

}



# ################# #################
#' Diagnostic plots and quantities arising from estimating the number of underlying factors
#'
#' Given the data, this function  draws a scree plot and a plot of the eigenvalue ratios.
#' The eignevalue ratio test is used to estimate the number of factors. See Ahn and Horenstein(2013).
#' @param X an n x p data matrix with each row being a sample.
#' @param K.scree an \emph{optional} integer specifying the number of eigenvalues to be plotted in the scree plot. Default is min(n,p).
#' @param K.factors an \emph{optional} integer specifying the number of eigenvalues to be used for the eigenvalue ratio test. Default is min(n,p)/2.
#' @param robust a TRUE/FALSE indicating whether to use a robust covariance estimator if TRUE, or the sample covariance estimator. Default is TRUE.
#' @param cv a boolean, specifying whether or  not to run cross-validation for the tuning parameter. Default is TRUE. Only used if \code{robust} is TRUE.
#' @param tau \code{>0}, multiplier for the tuning parameter for Huber loss function. Default is 2. Only used if \code{robust} is TRUE and \code{cv} is FALSE. See details.
#' @param show.plot a TRUE/FALSE indicating whether to show the resulting plots. Default is FALSE.
#' @details The maximum eigenvalue ratio is marked differently on the plot.  The index of this maximum ratio gives the number of estimated factors.
#' @details If \code{show.plots=TRUE}, plots are output and user has to hit <Return> to see the second plot. Alternatively, one may use the plot method for this class.
#' @details The tuning parameter \code{= tau *  sigma * optimal rate } where \code{optimal rate } is the optimal rate for the tuning parameter. For details, see Fan et al.(2017). \code{sigma} is the standard deviation of the data.
#' @return An object with S3 class \code{farm.scree} containing:
#' \itemize{
#'  \item{\code{eigenvalues} }{Eigenvalues of the covariance matrix}
#'  \item{\code{proportions} }{Proportion of variance explained by the principal components}
#'  \item{\code{eigenvalue.ratios} }{Ratios calculated in the eigenvalue ratio test}
#'  \item{\code{nfactors} }{Number of factors found using the eigenvalue ratio test}
#' }
#' @return If \code{show.plots=TRUE} function returns two plots: First plot is the scree plot of the data. Second plot illustrates the eigenvalue ratio test.
#' @seealso \code{\link{plot.farm.scree}} and \code{\link{print.farm.scree}}
#' @examples
#' set.seed(100)
#' p = 100
#' n = 20
#' epsilon = matrix(rnorm( p*n, 0,1), nrow = n)
#' B = matrix(rnorm(p*3,0,1), nrow=p)
#' fx = matrix(rnorm(3*n, 0,1), nrow = n)
#' X = fx%*%t(B)+ epsilon
#' output = farm.scree(X,show.plot = TRUE)
#' output = farm.scree(X,show.plot = FALSE, cv=FALSE, K.scree=5, K.factors =10)
#' output
#' plot(output, scree.plot=FALSE, col="blue", main="Customized plot")
#'
#' @references Ahn, S. C. and Horenstein, A. R.  (2013). "Eigenvalue Ratio Test for the Number of Factors," Econometrica, 81 (3), 1203–1227.
#' @references Fan, J., Ke, Y., Sun, Q. and Zhou, W-X. (2017). "FARM-Test: Factor-Adjusted Robust Multiple Testing with False Discovery Control", \url{https://arxiv.org/abs/1711.05386}.
#' @references Zhou, W-X., Bose, K., Fan, J. and Liu, H. (2017). "A New Perspective on Robust M-Estimation: Finite Sample Theory and Applications to Dependence-Adjusted Multiple Testing," Annals of Statistics, to appear, \url{https://arxiv.org/abs/1711.05381}.
#' @export
farm.scree<- function(X, K.scree = NULL , K.factors = NULL , robust = TRUE,cv=TRUE,tau=2, show.plot=FALSE){
  X = t(X)
  n = NCOL(X)
  p = NROW(X)
  if(tau<=0) stop('tau should be a positive number')

  if(min(n,p) <=4) stop('n and p must be at least 4')
  K.scree <- if (is.null(K.scree)) floor(min(n,p)) else K.scree
  if(K.scree <= 0)stop('number of eigenvalues to be plotted cannot be <= 0')

  K.factors <- if (is.null(K.factors)) floor(min(n,p)/2) else K.factors
  if(K.factors>min(n,p)/2) warning('Number of eigenvalue ratios is > min(n,p)/2. May cause numerical inconsistencies')
  if(K.factors < 2)stop('number of eigenvalues in the ratio test cannot be < 2')


  if(robust==TRUE){
    if(cv==TRUE){
      muhatx = mu_robust( matrix(X, p, n))
      thetax = mu_robust( matrix(X^2, p, n))
      covx = Cov_Huber( X, matrix(muhatx, p, 1))
    }
    else{CT = tau*apply(X,1,sd)*sqrt(n/log(p*n))
    muhatx = mu_robust_noCV( matrix(X, p, n), matrix(CT, p,1))
    CT = tau*apply(X^2,1,sd)*sqrt(n/log(p*n))
    thetax = mu_robust_noCV( matrix(X^2, p, n), matrix(CT,p,1))
    CT = Cov_Huber_tune(X, tau)
    covx =  Cov_Huber_noCV(matrix((X),p,n), matrix(muhatx, p, 1), matrix(CT,p,p))
    }
    decomp = Eigen_Decomp(covx)
    eig = decomp[,p+1]
    eig = pmax(eig, 0)
  }
  else {
    pca_fit=stats::prcomp((X), center = TRUE, scale = TRUE)
    eig = (pca_fit$sdev)^2
  }
  props = eig / sum(eig)
  if(show.plot){
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
    on.exit(graphics::par(oldpar))}
  #eigenvalue ratio test
  ratio=c()
  for(i in 1:K.factors)
    ratio=append(ratio, eig[i]/eig[i+1])
  ratio = ratio[is.finite(ratio)]
  if(show.plot){
    graphics::plot(ratio,type="b",pch=19, ylim=c(min(ratio),max(ratio)),main = paste("Eigenvalue ratio plot:\n", which.max(ratio), "factor(s) found"),cex.main=1)
    graphics::points(x = which.max(ratio), max(ratio), col = "red",bg = "red", pch = 23,cex = 1)}
  value   =list(eigenvalues = eig, proportions = props ,  eigenvalue.ratios=  ratio, nfactors = which.max(ratio), K.scree=K.scree, K.factors=K.factors)

  attr(value, "class") <- "farm.scree"
  value
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
#' \item{significant}{The number of hypotheses rejected}
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

  list(rejected = rejected , alldata = alldata, significant = NROW(rejected))
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


#################### huber covariance calculation ##############################################
#' Covariance estimation with Huber's loss function
#'
#' This function estimates covariance of multivariate data using the Huber's loss. The tuning parameter is chosen by cross validation.
#' @param X an n x p data matrix with each row being a sample.
#' @param cv a boolean, specifying whether or  not to run cross-validation for the tuning parameter. Default is TRUE.
#' @param tau \code{>0}, multiplier for the tuning parameter for Huber loss function. Default is 2. Only used if \code{cv} is FALSE. See details.
#' @param verbose a boolean specifying whether to print runtime updates to the console. Default is TRUE.
#' @details The tuning parameter \code{= tau *  sigma * optimal rate } where \code{optimal rate } is the optimal rate for the tuning parameter. For details, see Fan et al.(2017). \code{sigma} is the standard deviation of the data.
#' @return A list with the following items
#' \item{covhat}{the covariance matrix}
#' @examples
#' set.seed(100)
#' p = 20
#' n = 10
#' X = matrix(rnorm( p*n, 0,1), nrow = n)
#' covhat = farm.cov(X)
#'
#' @references Huber, P.J. (1964). "Robust Estimation of a Location Parameter." The Annals of Mathematical Statistics, 35, 73–101.
#' @references Fan, J., Ke, Y., Sun, Q. and Zhou, W-X. (2017). "FARM-Test: Factor-Adjusted Robust Multiple Testing with False Discovery Control", \url{https://arxiv.org/abs/1711.05386}.
#' @references Zhou, W-X., Bose, K., Fan, J. and Liu, H. (2017). "A New Perspective on Robust M-Estimation: Finite Sample Theory and Applications to Dependence-Adjusted Multiple Testing," Annals of Statistics, to appear, \url{https://arxiv.org/abs/1711.05381}.
#' @export
farm.cov <- function (X, cv=TRUE, tau=2, verbose=FALSE){
  X = t(X)
  p  = NROW(X)
  n = NCOL(X)
  if(tau<=0) stop('tau should be a positive number')
    if(cv==TRUE){
      muhatx = mu_robust( matrix(X, p, n))
      if(verbose){cat("calculating covariance matrix for X...\n")}
      covx = Cov_Huber( X, matrix(muhatx, p, 1))
    }
    else{CT = tau*apply(X,1,sd)*sqrt(n/log(p*n))
    muhatx = mu_robust_noCV( matrix(X, p, n), matrix(CT, p,1))
    CT = Cov_Huber_tune(X, tau)
    if (verbose){cat("calculating covariance matrix for X...\n")}
    covx =  Cov_Huber_noCV(matrix((X),p,n), matrix(muhatx, p, 1), matrix(CT,p,p))
    }
  return(covx)
}

#################### huber mean calculation ##############################################
#' Mean estimation with Huber's loss function
#'
#' This function estimates mean of multivariate data using the Huber's loss. The tuning parameter is chosen by cross validation.
#' @param X a n x p data matrix with each row being a sample.
#' @param cv a boolean, specifying whether or  not to run cross-validation for the tuning parameter. Default is TRUE.
#' @param tau \code{>0}, multiplier for the tuning parameter for Huber loss function. Default is 2. Only used if \code{cv} is FALSE. See details.
#' @param verbose a boolean specifying whether to print runtime updates to the console. Default is TRUE.
#' @return A list with the following items
#' \item{muhat}{the mean vector}
#' @details The tuning parameter \code{= tau *  sigma * optimal rate } where \code{optimal rate } is the optimal rate for the tuning parameter. For details, see Fan et al.(2017). \code{sigma} is the standard deviation of the data.
#' @examples
#' set.seed(100)
#' p = 20
#' n = 10
#' X = matrix(rnorm( p*n, 0,1), nrow = n)
#' muhat = farm.mean(X)
#'
#' @references Huber, P.J. (1964). "Robust Estimation of a Location Parameter." The Annals of Mathematical Statistics, 35, 73–101.
#' @references Fan, J., Ke, Y., Sun, Q. and Zhou, W-X. (2017). "FARM-Test: Factor-Adjusted Robust Multiple Testing with False Discovery Control", \url{https://arxiv.org/abs/1711.05386}.
#' @references Zhou, W-X., Bose, K., Fan, J. and Liu, H. (2017). "A New Perspective on Robust M-Estimation: Finite Sample Theory and Applications to Dependence-Adjusted Multiple Testing," Annals of Statistics, to appear, \url{https://arxiv.org/abs/1711.05381}.
#' @export
farm.mean <- function(X, cv=TRUE, tau=2, verbose=FALSE){
  X = t(X)
  p  = NROW(X)
  n = NCOL(X)
  if(tau<=0) stop('tau should be a positive number')
  if(cv==TRUE){
      muhatx = mu_robust( matrix(X, p, n))
    }
    else{CT = tau*apply(X,1,sd)*sqrt(n/log(p*n))
    muhatx = mu_robust_noCV( matrix(X, p, n), matrix(CT, p,1))
    }

  return(muhatx)
}


