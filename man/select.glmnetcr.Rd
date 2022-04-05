\name{select.glmnetcr}
\alias{select.glmnetcr}
\title{Select Step of Optimal Fitted AIC or BIC CR Model}
\description{
This function returns the step number associated with the minimum AIC or BIC models for a \code{glmnet.cr} object and so can be used in selecting the step number in calls in functions such as \code{coef.glmnet.cr}, \code{nonzero.glmnet.cr}, and \code{fitted.glmnet.cr}}
\usage{
select.glmnetcr(fit, which = "BIC")
}
\arguments{
  \item{fit}{a \code{glmnetcr} object}
  \item{which}{determines whether the step is associated with the optimal \code{which="BIC"} model or \code{which="AIC"} model}
}

\value{
A numeric value of length one representing the step number.
}
\author{
Kellie J. Archer, \email{archer.43@osu.edu}}


\seealso{
See Also as \code{\link{coef.glmnetcr}}, \code{\link{nonzero.glmnetcr}}, \code{\link{fitted.glmnetcr}}
}
\examples{
data(diabetes)
x <- diabetes[, 2:dim(diabetes)[2]]
y <- diabetes$y
glmnet.fit <- glmnetcr(x, y)
AIC <- select.glmnetcr(glmnet.fit, which = "AIC")
AIC
nonzero.glmnetcr(glmnet.fit, s = AIC)
}
\keyword{ misc }
