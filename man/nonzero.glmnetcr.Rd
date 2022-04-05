\name{nonzero.glmnetcr}
\alias{nonzero.glmnetcr}
\title{
Extract Non-Zero Model Coefficients}
\description{
The \code{nonzero.glmnetcr} function returns only those non-zero coefficient estimates for a selected model}
\usage{
nonzero.glmnetcr(fit, s)
}
\arguments{
  \item{fit}{a \code{glmnetcr} object}
  \item{s}{the step at which the non-zero coefficient estimates are desired}
}
\value{
  \item{a0}{intercept estimate}
  \item{beta}{non-zero estimates for variables and ordinal thresholds}
}
\author{
 Kellie J. Archer, \email{archer.43@osu.edu}
}
\seealso{
 See Also as \code{\link{glmnetcr}}, \code{\link{coef.glmnetcr}}, \code{\link{select.glmnetcr}}
}
\examples{
data(diabetes)
x <- diabetes[, 2:dim(diabetes)[2]]
y <- diabetes$y
glmnet.fit <- glmnetcr(x, y)
AIC.step <- select.glmnetcr(glmnet.fit, which = "AIC")
nonzero.glmnetcr(glmnet.fit, s = AIC.step)
}
\keyword{ misc }
