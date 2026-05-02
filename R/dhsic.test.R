##' Independence test based on dHSIC
##'
##' Hypothesis test for finding statistically significant evidence of
##' dependence between several variables. Uses the d-variable Hilbert
##' Schmidt independence criterion (dHSIC) as measure of
##' dependence. Several types of hypothesis tests are included. The
##' null hypothesis (H_0) is that all variables are jointly
##' independent.
##' 
##' @title Independence test based on dHSIC
##' @param X either a list of at least two numeric matrices or a
##'   single numeric matrix. The rows of a matrix correspond to the
##'   observations of a variable. It is always required that there are
##'   an equal number of observations for all variables (i.e. all
##'   matrices have to have the same number of rows). If \code{X} is a
##'   single numeric matrix than one has to specify the second
##'   variable as \code{Y} or set \code{matrix.input} to "TRUE". See
##'   below for more details.
##' @param Y a numeric matrix if \code{X} is also a numeric matrix and
##'   omitted if \code{X} is a list.
##' @param K a list of the gram matrices corresponding to each
##'   variable. If \code{K} specified the the following inputs
##'   \code{X}, \code{Y}, \code{kernel}, \code{pairwise},
##'   \code{bandwidth} and \code{matrix.input} will be ignored.
##' @param alpha a numeric value in (0,1) specifying the confidence
##'   level of the hypothesis test.
##' @param method a character string specifying the type of hypothesis
##'   test used. The available options are: "gamma" (gamma
##'   approximation based test), "permutation" (permutation test
##'   (slow)), "bootstrap" (bootstrap test (slow)) and "eigenvalue"
##'   (eigenvalue based test).
##' @param kernel a vector of character strings specifying the kernels
##'   for each variable. There exist two pre-defined kernels:
##'   "gaussian" (Gaussian kernel with median heuristic as bandwidth)
##'   and "discrete" (discrete kernel). User defined kernels can also
##'   be used by passing the function name as a string, which will
##'   then be matched using \code{\link[base]{match.fun}}. If the
##'   length of \code{kernel} is smaller than the number of variables
##'   the kernel specified in \code{kernel[1]} will be used for all
##'   variables.
##' @param B an integer value specifying the number of Monte-Carlo
##'   iterations made in the permutation and bootstrap test. Only
##'   relevant if \code{method} is set to "permutation" or to
##'   "bootstrap".
##' @param pairwise a logical value indicating whether one should use
##'   HSIC with pairwise comparisons instead of dHSIC. Can only be
##'   true if there are more than two variables.
##' @param bandwidth a numeric value specifying the size of the
##'   bandwidth used for the Gaussian kernel. Only used if
##'   kernel="gaussian.fixed".
##' @param matrix.input a boolean. If \code{matrix.input} is "TRUE"
##'   the input \code{X} is assumed to be a matrix in which the
##'   columns correspond to the variables.
##'
##' @details The d-variable Hilbert Schmidt independence criterion is
##'   a direct extension of the standard Hilbert Schmidt independence
##'   criterion (HSIC) from two variables to an arbitrary number of
##'   variables. It is 0 if and only if the variables are jointly
##'   independent.
##'
##' 4 different statistical hypothesis tests are implemented all with
##' null hypothesis (H_0: \code{X[[1]]},...,\code{X[[d]]} are jointly
##' independent) and alternative hypothesis (H_A:
##' \code{X[[1]]},...,\code{X[[d]]} are not jointly independent):
##' 1. Permutation test for dHSIC: exact level, slow 2. Bootstrap test
##' for dHSIC: pointwise asymptotic level and pointwise consistent,
##' slow 3. Gamma approximation based test for dHSIC: only
##' approximate, fast 4. Eigenvalue based test for dHSIC: pointwise
##' asymptotic level and pointwise consistent, medium
##'
##' The null hypothesis is rejected if \code{statistic} is strictly
##' greater than \code{crit.value}.
##'
##' If \code{X} is a list with d matrices, the function tests for
##' joint independence of the corresponding d random vectors. If
##' \code{X} is a matrix and \code{matrix.input} is "TRUE" the
##' functions tests the independence between the columns of
##' \code{X}. If \code{X} is a matrix and \code{matrix.input} is
##' "FALSE" then \code{Y} needs to be a matrix, too; in this case, the
##' function tests the (pairwise) independence between the
##' corresponding two random vectors.
##'
##' For more details see the references.
##' 
##' @return A list containing the following components:
##'
##' \item{statistic}{the value of the test statistic}
##' \item{crit.value}{critical value of the hypothesis test. The null
##' hypothesis (H_0: joint independence) is rejected if
##' \code{statistic} is greater than \code{crit.value}.}
##' \item{p.value}{p-value of the hypothesis test, i.e. the
##' probability that a random version of the test statistic is greater
##' than \code{statistic} under the calculated null hypothesis (H_0:
##' joint independence) based on the data.}  \item{time}{numeric
##' vector containing computation times. \code{time[1]} is time to
##' compute Gram matrix, \code{time[2]} is time to compute dHSIC and
##' \code{time[3]} is the time to compute \code{crit.value} and
##' \code{p.value.}}  \item{bandwidth}{bandwidth used during the
##' computation. Only relevant if Gaussian kernel was used.}
##' 
##' @export
##'
##' @import stats
##'
##' @useDynLib dHSIC
##' @importFrom Rcpp sourceCpp
##' 
##' @author Niklas Pfister and Jonas Peters
##'
##' @references
##' Gretton, A., K. Fukumizu, C. H. Teo, L. Song, B. Sch{\"o}lkopf and A. J. Smola (2007).
##' A kernel statistical test of independence.
##' In Advances in Neural Information Processing Systems (pp. 585-592).
##' 
##' Pfister, N., P. B{\"u}hlmann, B. Sch{\"o}lkopf and J. Peters (2018).
##' Kernel-based Tests for Joint Independence.
##' Journal of the Royal Statistical Society, Series B.
##'
##' @seealso In order to only compute the test statistic without
##'   p-values, use the function \code{\link{dhsic}}.
##'
##' @keywords nonparametric, htest
##'
##' @examples
##'
##' ### pairwise independent but not jointly independent (pairwise HSIC vs dHSIC)
##' set.seed(0)
##' x <- matrix(rbinom(100,1,0.5),ncol=1)
##' y <- matrix(rbinom(100,1,0.5),ncol=1)
##' z <- matrix(as.numeric((x+y)==1)+rnorm(100),ncol=1)
##' X <- list(x,y,z)
##' 
##' dhsic.test(X, method="permutation",
##'            kernel=c("discrete", "discrete", "gaussian"),
##'            pairwise=TRUE, B=1000)$p.value
##' dhsic.test(X, method="permutation",
##'            kernel=c("discrete", "discrete", "gaussian"),
##'            pairwise=FALSE, B=1000)$p.value


dhsic.test <- function(X, Y, K,
                       alpha=0.05,
                       method = "permutation",
                       kernel = "gaussian",
                       B = 1000,
                       pairwise = FALSE,
                       bandwidth=1,
                       matrix.input=FALSE){

  ############
  # Pairwise #
  ############
  
  if(pairwise & missing(K)){
    # Deal with matrix input
    if(matrix.input){
      if(!is.matrix(X)){
        stop("X needs to be a matrix if matrix.input=TRUE")
      }
      else{
        X <- split(X, rep(1:ncol(X), each = nrow(X)))
      }
    }
    if(!missing(Y)||!is.list(X)||length(X)<=2){
      stop("pairwise only makes sense if number of variables is greater than 2")
    }
    d <- length(X)
    for(j in 1:d){
      if(is.matrix(X[[j]])==FALSE){
        X[[j]] <- as.matrix(X[[j]])
      }
    }   
    len <- nrow(X[[1]])
    
    # Case: len<2d
    if(len<2*d){
      warning("Sample size is smaller than twice the number of variables. Test is trivial.")
      test=list(statistic=0,
                crit.value = Inf,
                p.value = 1,
                time = c(GramMat=0,dHSIC=0,CritVal=0))
      return(test)
    }
    
    pValVec <- rep(0,d-1)
    statVec <- rep(0,d-1)
    critVec <- rep(0,d-1)
    ptm <- proc.time()
    for(j in (d:2)){
      resTmp <- dhsic.test(do.call("cbind", X[1:(j-1)]), X[[j]], alpha=alpha/(d-1), method = method, kernel = kernel, B = B, pairwise=FALSE, bandwidth=bandwidth)
      pValVec[d-j+1] <- resTmp$p.value
      statVec[d-j+1] <- resTmp$statistic
      critVec[d-j+1] <- resTmp$crit.value
    }
    ind <- which.min(pValVec)
    stat <- statVec[ind]
    critical_value <- critVec[ind]
    # bonferroni correction for p-value
    p_value = min(pValVec[ind]*(d-1),1)
    timeCritVal <- as.numeric((proc.time()-ptm)[1])
    timeGramMat <- as.numeric(resTmp$time[1])
    timeHSIC <- as.numeric(resTmp$time[2])
  }
  
  ###################
  # d-variable HSIC #
  ###################
  
  else{
    # Check whether gram matrices have been specified in K and compute them if not
    if(missing(K)){

      ###
      # Prerequisites
      ###

      # Deal with matrix input
      if(matrix.input){
        if(!is.matrix(X)){
          stop("X needs to be a matrix if matrix.input=TRUE")
        }
        else{
          X <- split(X, rep(1:ncol(X), each = nrow(X)))
        }
      }
    
      # Check if Y has been passed as argument (implies two variable case)
      if(!missing(Y)){
        X <- list(X,Y)
      }
    
      # Set d, convert to matricies and set len
      d <- length(X)
      for(j in 1:d){
        if(is.matrix(X[[j]])==FALSE){
          X[[j]] <- as.matrix(X[[j]])
        }
      }     
      len <- nrow(X[[1]])
      
      # Case: len<2d
      if(len<2*d){
        warning("Sample size is smaller than twice the number of variables. Test is trivial.")
        test=list(statistic=0,
                  crit.value = Inf,
                  p.value = 1,
                  time = c(GramMat=0,dHSIC=0,CritVal=0),
                  bandwidth=NA)
        return(test)
      }
      
      # Check if enough kernels where set given the number of variables (else only used first kernel)
      if(length(kernel)<d){
        kernel <- rep(kernel[1],d)
      }
      if(length(bandwidth)<d){
        bandwidth <- rep(bandwidth[1],d)
      }
    
      # Define median heuristic bandwidth function
      median_bandwidth <- function(x){
        bandwidth <- median_bandwidth_rcpp(x[sample(len, replace=FALSE),,drop=FALSE],len,ncol(x))
        if(bandwidth==0){
          bandwidth <- 0.001
        }
        return(bandwidth)
      }

      # Define custom kernel Gram-matrix
      custom_grammat <- function(x,fun){
        KX <- matrix(nrow=len,ncol=len)
        for (i in 1:len){
          for (j in i:len){
            KX[i,j] <- match.fun(fun)(x[i,],x[j,])
            KX[j,i] <- KX[i,j]
          }
        }
        return(KX)
      }
      
      ###
      # Compute Gram-matrix (using specified kernels)
      ###
      K <- vector("list", d)
      ptm <- proc.time()
      for(j in 1:d){
        if(kernel[j]=="gaussian"){
          bandwidth[j] <- median_bandwidth(X[[j]])
          K[[j]] <- gaussian_grammat_rcpp(X[[j]],bandwidth[j],len,ncol(X[[j]]))
        }
        else if(kernel[j]=="gaussian.fixed"){
          K[[j]] <- gaussian_grammat_rcpp(X[[j]],bandwidth[j],len,ncol(X[[j]]))
        }
        else if(kernel[j]=="discrete"){
          bandwidth[j] <- NA
          K[[j]] <- discrete_grammat_rcpp(X[[j]],len,ncol(X[[j]]))
        }
        else{
          bandwidth[j] <- NA
          K[[j]] <- custom_grammat(X[[j]],kernel[j])
        }
      }
      timeGramMat <- as.numeric((proc.time() - ptm)[1])
    }
    else{
      # check if K is a list
      if(is.list(K)){
        d <- length(K)
        len <- nrow(K[[1]])
        timeGramMat <- NA
      }
      else{
        stop("K needs to be a list of matrices")
      }
      
    }
    
    
    ###
    # Main Computation
    ###
  
    ### Compute dHSIC
    ptm <- proc.time()
    term1 <- 1
    term2 <- 1
    term3 <- 2/len
    for (j in 1:d){
      term1 <- term1*K[[j]]
      term2 <- 1/len^2*term2*sum(K[[j]])
      term3 <- 1/len*term3*colSums(K[[j]])
    }
    term1 <- sum(term1)
    term3 <- sum(term3)
    dHSIC=1/len^2*term1+term2-term3
    timeHSIC <- as.numeric((proc.time() - ptm)[1])
    
    ### Using permutation test for dHSIC
    ptm <- proc.time()
    if(method == "permutation"){
      # Define permutation function
      dhsic_perm_fun <- function(i){
        term1 <- K[[1]]
        term2 <- 1/len^2*sum(K[[1]])
        term3 <- 2/len^2*colSums(K[[1]])
        for (j in 2:d){
          perm <- sample(1:len,replace=FALSE)
          K_perm <- K[[j]][perm, perm]
          term1 <- term1*K_perm
          term2 <- 1/len^2*term2*sum(K_perm)
          term3 <- 1/len*term3*colSums(K_perm)
        }
        term1 <- sum(term1)
        term3 <- sum(term3)
        return(1/len^2*term1+term2-term3)
      }
      # Perform performing
      dHSIC_perm <- sapply(1:B,dhsic_perm_fun)
      
      # Compute critical value and p-value
      sortdHSIC <- sort(len*dHSIC_perm)
      Bind <- sum(len*dHSIC==sortdHSIC)+ceiling((1-alpha)*(B+1))
      if(Bind<=B){
        critical_value <- sortdHSIC[Bind]
      }
      else{
        critical_value <- Inf
      }
      p_value <- (sum(dHSIC_perm>=dHSIC)+1)/(B+1)
    }
  
    ### Using bootstrap test for dHSIC
    if (method == "bootstrap"){
      # Define bootstrap function
      dhsic_boot_fun <- function(i){
        term1 <- K[[1]]
        term2 <- 1/len^2*sum(K[[1]])
        term3 <- 2/len^2*colSums(K[[1]])
        for (j in 2:d){
          boot <- sample(0:(len-1),replace=TRUE)
          K_boot <- K[[j]][boot, boot]
          term1 <- term1*K_boot
          term2 <- 1/len^2*term2*sum(K_boot)
          term3 <- 1/len*term3*colSums(K_boot)
        }
        term1 <- sum(term1)
        term3 <- sum(term3)
        return(1/len^2*term1+term2-term3)
      }
      # Perform bootstrapping
      dHSIC_boot <- sapply(1:B,dhsic_boot_fun)
      
      # Compute critical value and p-value
      sortdHSIC <- sort(len*dHSIC_boot)
      Bind <- sum(len*dHSIC==sortdHSIC)+ceiling((1-alpha)*(B+1))
      if(Bind<=B){
        critical_value <- sortdHSIC[Bind]
      }
      else{
        critical_value <- Inf
      }
      p_value <- (sum(dHSIC_boot>=dHSIC)+1)/(B+1)
    }
    
    ### Gamma approximation based test for dHSIC
    else if (method == "gamma"){
      # estimators
      est.a <- vector("numeric", d)
      est.b <- vector("numeric", d)
      est.c <- vector("numeric", d)
      for (j in 1:d){
        est.a[j] <- 1/(len^2)*sum(K[[j]])
        est.b[j] <- 1/(len^2)*sum(K[[j]]^2)
        est.c[j] <- 1/len^3*sum(rowSums(K[[j]])^2)
      }
      prod.a <- prod(est.a)
      prod.b <- prod(est.b)
      prod.c <- prod(est.c)
      oneoutprod.a <- vector("numeric",d)
      oneoutprod.b <- vector("numeric",d)
      oneoutprod.c <- vector("numeric",d)
      for(j in 1:d){
        oneoutprod.a[j] <- prod.a/est.a[j]
        oneoutprod.b[j] <- prod.b/est.b[j]
        oneoutprod.c[j] <- prod.c/est.c[j]
      }
      est.d <- est.a^2
      prod.d <- prod.a^2
      oneoutprod.d <- oneoutprod.a^2
      # exp-estimator
      exp.est <- (1-sum(oneoutprod.a)+(d-1)*prod.a)/len
      # var-estimtor
      term1 <- prod.b
      term2 <- (d-1)^2*prod.d
      term3 <- 2*(d-1)*prod.c
      term4 <- 0
      term5 <- 0
      term6 <- 0
      term7 <- 0
      for (r in 1:(d-1)){
        term4 <- term4+est.b[r]*oneoutprod.d[r]
        term5 <- term5+est.b[r]*oneoutprod.c[r]
        term6 <- term6+est.c[r]*oneoutprod.d[r]
        for (s in (r+1):d){
          term7 <- term7+2*est.c[r]*est.c[s]*oneoutprod.d[r]/est.d[s]
        }
      }
      term4 <- term4+est.b[d]*oneoutprod.d[d]
      term5 <- -2*(term5+est.b[d]*oneoutprod.c[d])
      term6 <- -2*(d-1)*(term6+est.c[d]*oneoutprod.d[d])
      
      factor1 <- len-2*d
      factor2 <- len*(len-1)*(len-2)
      for (j in 1:(2*d-3)){
        factor1 <- factor1*(len-2*d-j)
        factor2 <- factor2*(len-2-j)
      }
      var.est=2*factor1/factor2*(term1+term2+term3+term4+term5+term6+term7)
      # Compute critical value and p-value
      a <- (exp.est^2)/var.est
      b <- len*var.est/exp.est
      critical_value <- qgamma(1-alpha,shape=a,scale=b)
      p_value <- pgamma(len*dHSIC,shape=a,scale=b, lower.tail=FALSE)
    }
    
    
    ### Eigenvalue based test for dHSIC
    else if (method == "eigenvalue"){
      ## Use eigenvalue approximation to calculate the critical value
      # calculate H2
      est1 <- vector("numeric",d)
      est2 <- matrix(NA,len,d)
      for (j in 1:d){
        est1[j] <- 1/(len*(len-1))*sum(K[[j]]-diag(diag(K[[j]])))
        est2[,j] <- 1/len*rowSums(K[[j]])
      }
      est1_prod <- prod(est1)
      est2_prod <- apply(est2,1,prod)
      
      a1 <- matrix(1,len,len)
      a2 <- (d-1)^2*est1_prod
      a3 <- (d-1)*est2_prod
      a5 <- matrix(0,len,len)
      a6 <- matrix(0,len,len)
      a8 <- matrix(0,len,len)
      a9 <- matrix(0,len,1)
      for (j in 1:d){
        a1 <- a1*K[[j]]
        a5 <- a5+K[[j]]*est1_prod/est1[j]
        a6 <- a6+K[[j]]*matrix(rep(est2_prod/est2[,j],len),len,len)
        a9 <- a9+est2[,j,drop=FALSE]*est1_prod/est1[j]
        i <- j+1
        while ( i<=d ){
          a8 <- a8+est2[,j,drop=FALSE]%*%t(est2[,i,drop=FALSE])*est1_prod/(est1[j]*est1[i])+est2[,i,drop=FALSE]%*%t(est2[,j,drop=FALSE])*est1_prod/(est1[j]*est1[i])
          i <- i+1
        }
      }
      a3 <- matrix(rep(a3,len),len,len)
      a4 <- t(a3)
      a6 <- -a6
      a7 <- t(a6)
      a8 <- a8
      a9 <- matrix(rep((1-d)*a9,len),len,len)
      a10 <- t(a9)
      
      H2 <- 1/(d*(2*d-1))*(a1+a2+a3+a4+a5+a6+a7+a8+a9+a10)
      
      eigenvalues <- eigen(H2,only.values = TRUE,symmetric = TRUE)$values/len
      
      # simulate the sum of chi-squared distribution
      M <- 5000
      Z <- rnorm(M*len)^2
      chi_dist <- d*(2*d-1)*colSums(matrix(Z*eigenvalues,len,M))
      critical_value <- as.numeric(quantile(chi_dist,1-alpha))
      p_value <- 1-ecdf(chi_dist)(dHSIC*len)
    }
    timeCritVal <- as.numeric((proc.time()-ptm)[1])
    stat <- dHSIC*len
  }
  
  ###
  # Collect result
  ###
  test=list(statistic=stat,
            crit.value = critical_value,
            p.value = p_value,
            time = c(GramMat=timeGramMat,dHSIC=timeHSIC,CritVal=timeCritVal),
            bandwidth=bandwidth)
  
  return(test)
}
