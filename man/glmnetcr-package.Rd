\name{glmnetcr-package}
\alias{glmnetcr-package}
\docType{package}
\title{
Penalized Constrained Continuation Ratio Models for Ordinal Response Prediction using 'glmnet'
}
\description{
This package provides a function \code{glmnetcr} for fitting penalized constrained continuation ratio models for predicting an ordinal response and associated methods for plotting, printing, extracting predicted classes and probabilities, and extracting estimated coefficients for selected models in the regularization path. }
\details{
The DESCRIPTION file:
\packageDESCRIPTION{glmnetcr}
\packageIndices{glmnetcr}
This package contains functions for fitting penalized constrained continuation ratio models and extracting estimated coefficients, predicted class, and fitted probabilities. The model and methods can be used when the response to be predicted is ordinal, and is particularly relevant when there are more covariates than observations.
}
\author{
\packageAuthor{glmnetcr} Kellie J. Archer <archer.43@osu.edu>

Maintainer: \packageMaintainer{glmnetcr} Kellie J. Archer <archer.43@osu.edu>
}
\references{
Archer K.J., Williams A.A.A. (2012) L1 penalized continuation ratio models for ordinal response prediction using high-dimensional datasets. \emph{Statistics in Medicine}, 31(14), 1464-74.
}
\keyword{ package }
\keyword{ models }
\keyword{ classif }

\seealso{
 See also \code{\link[glmnet:glmnet]{glmnet}} ~~
}
\examples{
data(diabetes)
x <- diabetes[, 2:dim(diabetes)[2]]
y <- diabetes$y
glmnet.fit <- glmnetcr(x, y)
AIC <- select.glmnetcr(glmnet.fit, which="AIC")
fitted(glmnet.fit, s=AIC)
}
