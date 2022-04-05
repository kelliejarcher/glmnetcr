\name{cr.backward}
\alias{cr.backward}
\title{
Backward Continuation Ratio Restructure Function}
\description{
This is an internal function used by \code{glmnetcr} to restructure an ordinal dataset to represent the K-1 conditionally independent likelihoods needed for fitting a backrward continuation ratio model}
\usage{
cr.backward(x, y, weights)
}
\arguments{
  \item{x}{covariates to be included in the model}
  \item{y}{ordinal outcome}
   \item{weights}{observation weights. Can be total counts if responses are proportion matrices. Default is 1 for each observation}
}
\value{
  \item{object}{Yields an object where the first column y is the ordinal outcome that has been dichotomously coded and the expanded covariate matrix for the K-1 conditionally independent likelihoods}
}
\author{
Kellie J. Archer}
\seealso{
See Also as \code{\link{cr.forward}}, \code{\link{glmnetcr}}
}
\examples{
data(diabetes)
x <- diabetes[,2:dim(diabetes)[2]]
y <- diabetes$y
fit <- glmnetcr(x, y, method = "backward")
}
\keyword{ internal }
