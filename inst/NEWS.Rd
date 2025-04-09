\name{NEWS}
\title{News for Package 'glmnetcr'}
\section{Changes in version 1.0.3}{
\itemize{
  \item Change syntax to remove the period from all relevant functions. The main modeling fitting function is now glmnetcr which yields an object of class glmnetcr.
  \item Other functions are: nonzero.glmnetcr, select.glmnetcr, as well as the coef.glmnetcr, fitted.glmnetcr, plot.glmnetcr, predict.glmnetcr, print.glmnetcr.
  \item glmnetcr function now handles weights.
  \item Housekeeping: Cleaned up the DESCRIPTION, NAMESPACE, and URLs in help files.
  }
}
\section{Changes in version 1.0.4}{
\itemize{
   \item Edited glmnetcr-package.Rd to ensure it would be in sync with the rest of the package.
}
}
\section{Changes in version 1.0.5}{
\itemize{
   \item Changed code to be compatible with R 4.0.0 class function returning a vector of length greater than 1 when the object is a matrix.
   \item Added functionality of additional parameters from glmnet package to glmnetcr call.
}
}
\section{Changes in version 1.0.6}{
\itemize{
   \item Changed glmnetcr function to include the default of exclude = NULL which had not been previously implemented.
}
}
\section{Changes in version 1.0.7}{
\itemize{
   \item Fixed help Rd files to include link to package functions outside of this package.
}
}
