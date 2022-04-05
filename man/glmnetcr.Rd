\name{glmnetcr}
\alias{glmnetcr}
\title{Fit a Penalized Constrained Continuation Ratio Model Using Lasso or Elasticnet Regularization Via 'glmnet'}
\description{
   Fit a constrained continuation ratio model via penalized maximum likelihood.  The
  regularization path is computed for the lasso or elasticnet penalty at a grid
  of values for the regularization parameter lambda. Can deal with all
  shapes of data, including very large sparse data matrices.
}
\usage{
glmnetcr(x, y, method = "backward", weights = NULL, offset = NULL, alpha = 1,
nlambda = 100, lambda.min.ratio = NULL, lambda = NULL, standardize = TRUE,
thresh = 1e-04,  exclude = NULL, penalty.factor = NULL, maxit = 100,
dfmax = nvars + 1, pmax = min(dfmax * 2 + 20, nvars),
type.logistic = c("Newton","modified.Newton"),
trace.it = 0 )
}
\arguments{
  \item{x}{
input data, of dimension nobs x nvars; each row is an observation vector. Can be in sparse column format (class \code{"dgCMatrix"} as in package \code{Matrix})}
  \item{y}{
ordinal response variable}
  \item{method}{
select between fitting a backward (default) versus a forward continuation ratio model}
  \item{weights}{
observation weights. Can be total counts if responses are proportion matrices. Default is 1 for each observation}
  \item{offset}{
A vector of length \code{nobs} that is included in the linear predictor (a \code{nobs x nc} matrix for the \code{"multinomial"} family). Useful for the \code{"poisson"} family (e.g. log of exposure time), or for refining a model by starting at a current fit. Default is \code{NULL}. If supplied, then values must also be supplied to the \code{predict.glmnet.cr} function.}
  \item{alpha}{
The elasticnet mixing parameter, with
    \eqn{0<\alpha\le 1}. The penalty is defined
    as \deqn{(1-\alpha)/2||\beta||_2^2+\alpha||\beta||_1.} \code{alpha=1}
    is the lasso penalty.}
  \item{nlambda}{
The number of \code{lambda} values - default is 100.}
  \item{lambda.min.ratio}{
Smallest value for \code{lambda}, as a fraction of
    \code{lambda.max}, the (data derived) entry value (i.e. the smallest
  value for which all coefficients are zero). The default depends on the
  sample size \code{nobs} relative to the number of variables
  \code{nvars}. If \code{nobs > nvars}, the default is \code{0.0001},
  close to zero.  If \code{nobs < nvars}, the default is \code{0.01}.
  A very small value of
  \code{lambda.min.ratio} will lead to a saturated fit in the \code{nobs <
  nvars} case. This is undefined for
  \code{"binomial"} and \code{"multinomial"} models, and \code{glmnet}
  will exit gracefully when the percentage deviance explained is almost
  1.}
  \item{lambda}{
A user supplied \code{lambda} sequence. Typical usage
    is to have the
    program compute its own \code{lambda} sequence based on
    \code{nlambda} and \code{lambda.min.ratio}. Supplying a value of
    \code{lambda} overrides this. WARNING: use with care. Do not supply
  a single value for \code{lambda} (for predictions after CV use \code{predict()}
  instead).  Supply instead
    a decreasing sequence of \code{lambda} values. \code{glmnet} relies
  on its warms starts for speed, and its often faster to fit a whole
  path than compute a single fit.}
  \item{standardize}{
Logical flag for variable standardization, prior to
    fitting the model sequence. The coefficients are always returned on
    the original scale. Default is \code{standardize=TRUE}.
  If variables are in the same units already, you might not wish to standardize.}
  \item{thresh}{Convergence threshold for coordinate descent. Each inner
  coordinate-descent loop continues until the maximum change in the
  objective after any coefficient update is less than \code{thresh}
  times the null deviance. Defaults value is \code{1E-6}.}
  \item{exclude}{
Indices of variables to be excluded from the
    model. Default is none. Equivalent to an infinite penalty factor
    (next item).}
  \item{penalty.factor}{
Separate penalty factors can be applied to each
    coefficient. This is a number that multiplies \code{lambda} to allow
    differential shrinkage. Can be 0 for some variables, which implies
    no shrinkage, and that variable is always included in the
    model. Default is 1 for all variables (and implicitly infinity for
    variables listed in \code{exclude}).}
  \item{maxit}{
Maximum number of outer-loop iterations for
      \code{"binomial"} or \code{"multinomial"} families. Default is 100.}
   \item{dfmax}{Limit the maximum number of variables in the model. Useful for very large \code{nvars}, if a partial path is desired.}
   \item{pmax}{Limit the maximum number of variables ever to be nonzero.}
   \item{type.logistic}{If \code{"Newton"} then the exact hessian is used (default), while \code{"modified.Newton"} uses an upper-bound on the hessian, and can be faster.}
   \item{trace.it}{If \code{trace.it=1}, then a progress bar is displayed; useful for big models that take a long time to fit.}
}
\details{
  The sequence of models implied by \code{lambda} is fit by coordinate
  descent. This is a lasso or elasticnet regularization path
  for fitting the generalized linear regression
  paths, by maximizing the appropriate penalized log-likelihood. Sometimes the sequence is truncated before \code{nlambda}
  values of \code{lambda} have been used, because of instabilities in
  the inverse link functions near a saturated fit. This function invokes \code{glmnet(...,family="binomial")}
  which fits a traditional logistic regression model for the
  log-odds. The
  penalties take care of redundancies.
}
\value{
  \item{call}{the call that produced this object}
  \item{a0}{Intercept sequence of length \code{length(lambda)}}
  \item{beta}{a \code{nvars x
      length(lambda)} matrix of coefficients}
  \item{lambda}{The actual sequence of \code{lambda} values used}
  \item{dev}{The fraction of (null) deviance explained (for \code{"elnet"}, this
      is the R-square). The deviance calculations incorporate weights if
  present in the model. The deviance is defined to be -2*(loglike_sat -
  loglike), where loglike_sat is the log-likelihood for the saturated
  model (a model with a free parameter per observation). Hence dev.fraction=1-dev/nulldev.}
    \item{nulldev}{Null deviance (per observation). This is defined to
  be  -2*(loglike_sat -loglike(Null)); The NULL model refers to the
  intercept model.}
  \item{df}{The number of nonzero coefficients for each value of
    \code{lambda}.}
  \item{dfmat}{Appliable for \code{"multnet"} only so does not apply to the CR model. A matrix consisting of the
    number of nonzero coefficients per class}
  \item{dim}{dimension of coefficient matrix (ices)}
  \item{nobs}{number of observations}
  \item{npasses}{total passes over the data summed over all lambda
    values}
  \item{offset}{a logical variable indicating whether an offset was included in the model}
  \item{jerr}{error flag, for warnings and errors (largely for internal debugging).}
 \item{x}{returns the original data matrix. Useful for \code{predict} and \code{fitted} functions}
 \item{y}{returns the original ordinal response. Useful for \code{predict} and \code{fitted} functions}
 \item{method}{returns the continuation ratio method (forward or backward). Useful for \code{predict} and \code{fitted} functions}
}
\author{
Kellie J. Archer, \email{archer.43@osu.edu}}
\note{
This function depends heavily upon the \code{glmnet} function in the \pkg{glmnet} package written by Jerome Friedman, Trevor Hastie and Rob Tibshirani
}

\seealso{
See Also as \code{\link{glmnet}}, \code{\link{select.glmnetcr}}, \code{\link{fitted.glmnetcr}}
}
\examples{
data(diabetes)

x <- diabetes[, 2:dim(diabetes)[2]]

y <- diabetes$y

glmnet.fit <- glmnetcr(x, y)

fitted(glmnet.fit, s = select.glmnetcr(glmnet.fit))
}
\keyword{ models }
\keyword{ classif }
