\name{plot.glmnetcr}
\alias{plot.glmnetcr}
\title{
Plots the Regularization Path Computed
}
\description{
This function takes a \code{glmnetcr} object and visualizes the regularization path. The horizontal axis can be \code{norm}, \code{lambda} or \code{step}. The vertical axis can be \code{coefficients}, \code{aic} or \code{bic}.
}
\usage{
\method{plot}{glmnetcr}(x, xvar = c("norm", "lambda", "step"),
type = c("coefficients", "aic", "bic"), omit.zero = TRUE, breaks = TRUE, mar = NULL,
eps = .Machine$double.eps, main = NULL, ...)
}
\arguments{
  \item{x}{
a \code{glmnetcr} object
}
  \item{xvar}{
horizontal axis. \code{xvar=norm} plots against the L1 norm of the coefficients (to which L1 norm penalty was applied); \code{xvar=lambda} plots against \eqn{\lambda}; and \code{xvar=step} plots against the number of steps taken. Default is "norm".}
  \item{type}{
type of the plot, or the vertical axis. Default is "coefficients".}
  \item{omit.zero}{
If \code{TRUE} and \code{type="coefficients"}, the predictors that were never in the active set are omitted.}
  \item{breaks}{
If \code{TRUE}, vertical lines are drawn at the points where the active set changes and numbered with the degrees of freedom.}
  \item{mar}{
margin relative to the current font size
}
  \item{eps}{
an effective zero
}
  \item{main}{
title of the plot
}
  \item{\dots}{
other options for the plot}
}
\author{
Kellie J. Archer \email{archer.43@osu.edu}
}
\note{
This function was adapted from the \pkg{glmpath} package by Park and Hastie for a \code{glmnetcr} object.
}

\seealso{
See also \code{\link{glmnetcr}}, \code{\link{predict.glmnetcr}}
}
\examples{
library(glmnetcr)
data(diabetes)
x<-diabetes[,2:dim(diabetes)[2]]
y<-diabetes$y
summary(y)
fit<-glmnetcr(x,y)
plot(fit)
plot(fit,xvar="step",type="bic")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Graphics }
