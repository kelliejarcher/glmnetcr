\name{fitted.glmnetcr}
\alias{fitted.glmnetcr}
\title{
AIC, BIC, Predicted Class, and Fitted Probabilities of Class Membership
}
\description{
For a given step, returns the AIC, BIC, predicted class, and the fitted probabilities for the K classes.
}
\usage{
\method{fitted}{glmnetcr}(object, newx = NULL, s, ...)
}
\arguments{
  \item{object}{a \code{glmnetcr} object}
  \item{newx}{a data matrix representing the predictor variables, if missing defaults to original data used in fitting the model}
  \item{s}{the step at which the fitted probabilities and class are desired}
  \item{\dots}{additional optional arguments}
}
\value{
  \item{AIC }{AIC at step s}
  \item{BIC }{BIC at step s}
  \item{class}{a vector of length n indicating the predicted class for each observation in newx at step s}
  \item{probs}{a matrix with n rows and K columns indicating the fitted class probabilities for each observation and class in new at step x}
}
\author{
Kellie J. Archer, \email{archer.43@osu.edu}
}
\seealso{
See Also as \code{\link{predict.glmnetcr}}, \code{\link{select.glmnetcr}}
}
\examples{
data(diabetes)
x <- diabetes[, 2:dim(diabetes)[2]]
y <- diabetes$y
glmnet.fit <- glmnetcr(x, y)
fitted(glmnet.fit, s = select.glmnetcr(glmnet.fit))
}
\keyword{ methods }