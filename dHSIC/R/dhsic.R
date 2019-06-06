##' d-variable Hilbert Schmidt independence criterion - dHSIC
##'
##' The d-variable Hilbert Schmidt independence criterion (dHSIC) is a
##' non-parametric measure of dependence between an arbitrary number
##' of variables. In the large sample limit the value of dHSIC is 0 if
##' thevariables are jointly independent and positive if there is
##' adependence. It is therefore able to detect any type of dependence
##' given a sufficient amount of data.
##' 
##' @title d-variable Hilbert Schmidt independence criterion - dHSIC
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
##'   variable. If \code{K} specified the other inputs will have no
##'   effect on the computations.
##' @param kernel a vector of character strings specifying the kernels
##'   for each variable. There exist two pre-defined kernels:
##'   "gaussian" (Gaussian kernel with median heuristic as bandwidth)
##'   and "discrete" (discrete kernel). User defined kernels can also
##'   be used by passing the function name as a string, which will
##'   then be matched using \code{\link[base]{match.fun}}. If the
##'   length of \code{kernel} is smaller than the number of variables
##'   the kernel specified in \code{kernel[1]} will be used for all
##'   variables.
##' @param bandwidth a numeric value specifying the size of the
##'   bandwidth used for the Gaussian kernel. Only used if
##'   kernel="gaussian.fixed".
##' @param matrix.input a boolean. If \code{matrix.input} is "TRUE"
##'   the input \code{X} is assumed to be a matrix in which the
##'   columns correspond to the variables.
##' 
##' @return A list containing the following components:
##'
##' \item{dHSIC}{the value of the empirical estimator of dHSIC}
##' \item{time}{numeric vector containing computation
##' times. \code{time[1]} is time to compute Gram matrix and
##' \code{time[2]} is time to compute dHSIC.}
##' \item{bandwidth}{bandwidth used during computations. Only relevant
##' if Gaussian kernel was used.}
##' 
##' @export
##'
##' @import stats
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
##' @seealso In order to perform hypothesis tests based on dHSIC use
##'   the function \code{\link{dhsic.test}}.
##'
##' @keywords nonparametric
##'
##' @examples
##'
##' ### Three different input methods
##' set.seed(0)
##' x <- matrix(rnorm(200),ncol=2)
##' y <- matrix(rbinom(100,30,0.1),ncol=1)
##' # compute dHSIC of x and y (x is taken as a single variable)
##' dhsic(list(x,y),kernel=c("gaussian","discrete"))$dHSIC
##' dhsic(x,y,kernel=c("gaussian","discrete"))$dHSIC
##' # compute dHSIC of x[,1], x[,2] and y
##' dhsic(cbind(x,y),kernel=c("gaussian","discrete"), matrix.input=TRUE)$dHSIC
##' 
##' ### Using a user-defined kernel (here: sigmoid kernel)
##' set.seed(0)
##' x <- matrix(rnorm(500),ncol=1)
##' y <- x^2+0.02*matrix(rnorm(500),ncol=1)
##' sigmoid <- function(x_1,x_2){
##'   return(tanh(sum(x_1*x_2)))
##' }
##' dhsic(x,y,kernel="sigmoid")$dHSIC


dhsic <- function(X, Y, K,
                  kernel="gaussian",
                  bandwidth=1,
                  matrix.input=FALSE){

  ###
  # Prerequisites
  ###

  # Check whether gram matrices have been specified in K and compute them if not
  if(missing(K)){
    # Check if Y has been passed as argument (implies two variable case)
    if(!missing(Y)){
      X <- list(X,Y)
    }
    
    # Deal with matrix input
    if(matrix.input){
      if(!is.matrix(X)){
        stop("X needs to be a matrix if matrix.input=TRUE")
      }
      else{
        X <- split(X, rep(1:ncol(X), each = nrow(X)))
      }
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
      warning("Sample size is smaller than twice the number of variables. dHSIC is trivial.")
      result=list(dHSIC = 0,
                  time = c(GramMat=0,HSIC=0))
      return(result)
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
      bandwidth <- median_bandwidth_rcpp(x[sample(1:len),,drop=FALSE],len,ncol(x))
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
  # Compute dHSIC
  ###
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

  ###
  # Collect result
  ###
  result=list(dHSIC = dHSIC,
              time = c(GramMat=timeGramMat,HSIC=timeHSIC),
              bandwidth=bandwidth)
  
  return(result)

}
