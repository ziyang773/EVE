#' @importFrom stats lm qnorm
NULL
#' Equivariant Variance Estimator (EVE)
#'
#' \code{eve} computes the Equivariant Variance Estimator (EVE) for a time
#' series under frequent mean changes. The method is based on fitting a linear
#' regression to the sequence \eqn{Y_k = T_k/2}, where \eqn{T_k} is constructed
#' from shifted squared differences of the input series. The variance estimate
#' is given by the intercept of the fitted regression.
#'
#' Users may specify a fixed parameter \eqn{K}; alternatively, if
#' \code{K = NULL}, the optimal \eqn{K} is selected automatically using a
#' data-driven procedure with maximum search bound \code{Kmax}.
#'
#' @param x Numeric vector. Input time vector.
#' @param K Integer (or vector of integers). The length of the \eqn{Y} vector
#'   to be used in estimating the variance. If \code{NULL}, \eqn{K} is selected
#'   automatically. Default: \code{NULL}.
#' @param Kmax Integer. Upper bound for the search when \code{K = NULL}. Default: 20.
#'
#' @return A list with:
#' \itemize{
#'   \item \code{sigmahat} — Named vector of estimated standard deviation(s) using EVE wits selected \eqn{K}s.
#'   \item \code{K} — The specified or automatically selected value(s) of \eqn{K}.
#' }
#'
#' @details
#' When \code{K = NULL}, \code{eve} internally calls \code{eve_sel_K()} to select
#' \eqn{K} using the standardized prediction error criterion. When \code{K} is
#' provided, regression is performed separately for each value of \eqn{K}.
#' A warning is issued if the estimated variance is negative.
#'
#' @references
#' Hao, N., Niu, Y. S., & Xiao, H. (2023).
#' \emph{Equivariant variance estimation for multiple change-point model}.
#' Electronic Journal of Statistics, 17(2), 3811–3853. arXiv:2108.09431.
#'
#' @examples
#' set.seed(1234)
#' mu <- rep(c(rep(5, 50), rep(0, 50)), 20)
#' x <- rnorm(1000) + mu # true standard deviation is 1.0
#' eve(x)              # automatic selection of K
#' eve(x, K = 5)       # user-specified K
#' eve(x, K = c(5,10,15,20)) # compare different K values
#'
#' @export
eve <- function(x, K = NULL, Kmax=20){ #If K is NULL, we will select K with Kmax = Kmax
  if (is.null(K)) return(eve_sel_K(x, Kmax = Kmax))
  n <- length(x)
  Kmax <- max(K)
  xx <- c(x,x[1:Kmax])                           # glue top Kmax obvervations to tail
  Y <- numeric(Kmax)
  for (k in 1:Kmax){
    Y[k] <- mean((x-xx[(k+1):(n+k)])^2)/2        # calculate Y_k=T_k/2
  }
  sigmahat <- numeric(length(K))
  for (i in 1:length(K)){
    cK <- K[i]
    obj  <- lm(Y[1:cK]~c(1:cK))
    variancehat <- obj$coefficients[1]
    if (variancehat < 0) {
      print("warning: K is too small or there is a linear trend in data.") # give a warning if negative
      variancehat <- 0
    }
    sigmahat[i] <- sqrt(variancehat)                   # our method
  }
  names(sigmahat) = "EVE"
  if (length(K) > 1) names(sigmahat) = K
  return(list(sigmahat = sigmahat, K = K))
}

#' EVE variance estimator with data-driven K (internal)
#'
#' \code{eve_sel_K} computes the EVE variance estimator using a data-driven
#' truncation parameter \eqn{K}. For each \eqn{k = 1,\dots,K_{\max}}, it constructs
#' the sequence \eqn{Y_k = T_k/2} from the shifted squared differences of
#' the input series, and then selects \eqn{K} via \code{selK()}. The variance is
#' estimated by regressing \eqn{Y_1,\dots,Y_K} on \eqn{1, \dots, K} and taking the
#' intercept. A warning is issued if the estimated variance becomes negative.
#'
#' This function serves as the backend for the exported interface \code{EVE()} in the package.
#'
#' @param x Numeric vector. Input time series.
#' @param Kmax Integer. Maximum candidate value of \eqn{K}. Default = 20.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{sigmahat} — Estimated standard deviation by EVE with selected \eqn{K}.
#'   \item \code{K} — The selected value of \eqn{K}.
#' }
#'
#' @keywords internal
eve_sel_K  <-  function(x, Kmax = 20){
  n <- length(x)
  xx  <-  c(x,x[1:Kmax])                             # glue top Kmax obvervations to tail
  Y  <-  numeric(Kmax)
  for (k in 1:Kmax){
    Y[k] <- mean((x-xx[(k+1):(n+k)])^2)/2            # calculate Y_k=T_k/2
  }
  K <-  selK(Y)                                      # select K
  obj <- lm(Y[1:K]~c(1:K))
  variancehat <- obj$coefficients[1]
  if (variancehat < 0) {
    print("warning: K is too small or there is a linear trend in data.") # give a warning if negative
    variancehat <- 0
  }
  sigmahat  <-  sqrt(variancehat)               # our method
  names(sigmahat) <- "EVEsel"
  return(list(sigmahat = sigmahat, K = K))
}


#' Select the tuning parameter K for EVE (internal)
#'
#' \code{selK} determines the truncation parameter \eqn{K} in the EVE variance
#' estimator. For each candidate \eqn{k = 5, \dots, K_{\max} - 1}, it fits a linear
#' regression to \eqn{Y_1,\dots,Y_k}, evaluates the standardized one-step-ahead
#' prediction error for \eqn{Y_{k+1}}, and chooses the \eqn{k} that maximizes it.
#' If all standardized scores fall below the 0.99 quantile of the standard normal
#' distribution, the default \code{K = length(Y)} is returned. This function is internally called within
#' \code{eve_sel_K()}.
#'
#' @param Y Numeric vector of \eqn{Y_k} values used to construct the EVE estimator.
#'
#' @return Integer. The selected value of \eqn{K}. If \code{length(Y) < 6},
#'   the function returns \code{length(Y)}.
#'
#' @keywords internal
selK  <-  function(Y){
  Kmax <- length(Y)
  if (Kmax<6) {print("Warning: K is too small"); return(Kmax)}
  tempfit  <-  rep(0,Kmax)
  for (k in 5:(Kmax-1)){
    tempY  <-  Y[1:k]
    tempX  <-  1:k
    obj  <-  lm(tempY~tempX)
    tempv  <-  mean(obj$residuals^2)
    tempfit[k]  <-  abs(Y[k+1]-sum(obj$coefficients*c(1,(k+1))))/sqrt(tempv)
  }
  if (max(tempfit) < qnorm(0.99)) {
    #print("Warning: default K = Kmax is used")
    return (Kmax)          # nothing significant, return Kmax
  }
  return(which.max(tempfit))
}

#' Müller–Stadtmüller variance estimator (internal)
#'
#' \code{MSve} computes the classical difference-based variance estimator of
#' Müller and Stadtmüller. For \eqn{k = 1,\dots,K_{\max}}, it constructs
#' \eqn{Y_k = \frac{1}{2n} \sum_{t=1}^{n-k} (x_t - x_{t+k})^2} and estimates the
#' variance by regressing \eqn{Y_1,\dots,Y_K} on \eqn{1,\dots,K}, taking the
#' intercept. A warning is issued if the estimated variance is negative.
#'
#' This function is provided for comparison with the EVE estimator.
#'
#' @param x Numeric vector. Input time series.
#' @param K Integer (or vector of integers). The length of the \eqn{Y} vector
#'   to be used. Default: \code{10}.
#'
#' @return A list containing \code{sigmahat} and \code{K}.
#' @importFrom stats filter na.omit
#'
#' @keywords internal
MSve <- function(x, K = 10){
  n <- length(x)
  Kmax <- max(K)
  YS = numeric(Kmax)
  for (k in 1:Kmax){
    YS[k]=sum((x[1:(n-k)]-x[(k+1):n])^2)/(2*n)
  }
  sigmahat <- numeric(length(K))
  for (i in 1:length(K)){
    cK <- K[i]
    obj  <- lm(YS[1:cK]~c(1:cK))
    variancehat <- obj$coefficients[1]
    if (variancehat < 0) {
      print("warning: K is too small or there is a linear trend in data.") # give a warning if negative
      variancehat <- 0
    }
    sigmahat[i] <- sqrt(variancehat)
  }
  names(sigmahat) = "MS"
  if (length(K) > 1) names(sigmahat) = K
  return(list(sigmahat = sigmahat, K = K))
}
