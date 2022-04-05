\name{predict.glmnetcr}
\alias{predict.glmnetcr}
\title{
AIC, BIC, Predicted Class, and Fitted Probabilities for All Models}
\description{
Returns the AIC, BIC, predicted class, and the fitted probabilities for the K classes for all steps along the path.
}
\usage{
\method{predict}{glmnetcr}(object, newx = NULL, ...)
}
\arguments{
  \item{object}{a \code{glmnetcr} object}
  \item{newx}{a data matrix representing the predictor variables, if missing defaults to original data used in fitting the model}
  \item{\dots}{additional optional arguments}
}
\value{
  \item{AIC }{a numeric vector representing the AIC for each step}
  \item{BIC }{a numeric vector representing the BIC for each step}
  \item{class}{a matrix with n rows indicating the predicted class for each observation in newx for each stepp}
  \item{probs}{an array where for each step, a matrix of n rows and K columns indicates the fitted class probabilities for each observation and class}
}
\author{
Kellie J. Archer, \email{archer.43@osu.edu}
}
\seealso{
 See Also as \code{\link{glmnetcr}}, \code{\link{fitted.glmnetcr}}
}
\examples{
data(diabetes)
x <- diabetes[, 2:dim(diabetes)[2]]
y <- diabetes$y
glmnet.fit <- glmnetcr(x, y)
result <- predict(glmnet.fit)
names(result)
}
\keyword{ methods }