\name{print.glmnetcr}
\alias{print.glmnetcr}
\title{
Print a 'glmnetcr' Object}
\description{
Print a summary of the \code{glmnetcr} path at each step along the path.
}
\usage{
\method{print}{glmnetcr}(x, digits = max(3, getOption("digits") - 3), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{fitted \code{glmnetcr} object}
  \item{digits}{significant digits in printout}
  \item{\dots}{additional print arguments}
}

\details{
  The call that produced the object \code{x} is printed, followed by a three-column
  matrix
  with columns \code{Df}, \code{\%dev} and \code{Lambda}. The \code{Df}
    column is the number of nonzero coefficients (Df is a reasonable
    name only for lasso fits). \code{\%dev} is the percent deviance
    explained (relative to the null deviance).
}
\value{The matrix above is silently returned}
\references{Friedman, J., Hastie, T. and Tibshirani, R. (2008)
  \emph{Regularization Paths for Generalized Linear Models via Coordinate
  Descent}}
\author{Jerome Friedman, Trevor Hastie and Rob Tibshirani}
\note{This function is essentially the same as \code{print.glmnet} from the \pkg{glmnet} package by but was edited to operate on a returned \code{glmnetcr} object.}

\seealso{
See Also as \code{\link{glmnetcr}}
}
\examples{
data(diabetes)
x <- diabetes[, 2:dim(diabetes)[2]]
y <- diabetes$y
glmnet.fit <- glmnetcr(x, y)
print(glmnet.fit)
}
\keyword{ methods }
