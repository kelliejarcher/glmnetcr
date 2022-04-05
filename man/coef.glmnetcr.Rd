\name{coef.glmnetcr}
\alias{coef.glmnetcr}
\title{
Extract All Model Coefficients
}
\description{
The \code{coef.glmnetcr} function returns all estimated coefficients for a \code{glmnetcr} fitted model, where the model selected is indicated by step s.}
\usage{
\method{coef}{glmnetcr}(object, s, ...)
}
\arguments{
  \item{object}{a \code{glmnetcr} object}
  \item{s}{the step at which the coefficient estimates are desired}
  \item{\dots}{additional optional arguments}
}
\value{
\item{a0}{the intercept}
\item{beta}{estimated coefficients for the variables and the ordinal thresholds}
}
\author{
Kellie J. Archer, \email{archer.43@osu.edu}
}
\seealso{
See Also as \code{\link{glmnetcr}}, \code{\link{nonzero.glmnetcr}}, \code{\link{select.glmnetcr}}
}
\examples{
data(diabetes)
x <- diabetes[, 2:dim(diabetes)[2]]
y <- diabetes$y
glmnet.fit <- glmnetcr(x, y)
BIC.model <- select.glmnetcr(glmnet.fit)
estimates <- coef(glmnet.fit, s = BIC.model)
}
\keyword{ methods }
