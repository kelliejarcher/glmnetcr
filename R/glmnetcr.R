glmnetcr <-
function (x, y, method = "backward", weights=NULL, offset=NULL,
alpha = 1, nlambda = 100, lambda.min.ratio=NULL, lambda=NULL,
standardize = TRUE, thresh = 1e-04,
exclude = NULL, penalty.factor = NULL, maxit=100,
dfmax = nvars + 1,
pmax = min(dfmax * 2 + 20, nvars),
type.logistic = c("Newton","modified.Newton"),
trace.it = 0 )
{
    if (length(unique(y))==2) stop("Binary response: Use glmnet with family='binomial' parameter")
    penalty.factor <- unique(c(penalty.factor, exclude))
    n <- nobs <- dim(x)[1]
    p <- m <- nvars <- dim(x)[2]
    if (is.null(penalty.factor) & !is.null(exclude)) {
      penalty.factor<-rep(1,p)
      penalty.factor[exclude]<-Inf
    } else if (!is.null(penalty.factor) & !is.null(exclude)) {
      penalty.factor[exclude]<-Inf
    }
    k <- length(unique(y))
    x <- as.matrix(x)
	if (is.null(penalty.factor)) penalty.factor<-rep(1, nvars) else penalty.factor<-penalty.factor
	if (is.null(lambda.min.ratio)) lambda.min.ratio<-ifelse(nobs<nvars,0.01,0.0001)
	if (is.null(weights)) weights<-rep(1, length(y))
    if (c("backward", "forward")[charmatch(method, c("backward",
													 "forward"))] == "backward") {
        restructure <- cr.backward(x = x, y = y, weights=weights)
    }
    if (c("backward", "forward")[charmatch(method, c("backward",
													 "forward"))] == "forward") {
        restructure <- cr.forward(x = x, y = y, weights=weights)
    }
    glmnet.data <- list(x = restructure[, -c(1,2)], y = restructure[,
						"y"], weights = restructure[, "weights"])
    object <- glmnet(glmnet.data$x, glmnet.data$y, family = "binomial", weights=glmnet.data$weights, offset=offset, alpha = alpha, nlambda = nlambda, lambda.min.ratio = lambda.min.ratio, lambda=lambda,
					 standardize = standardize, thresh = thresh, dfmax = dfmax, pmax=pmax, exclude=exclude, penalty.factor = c(penalty.factor,rep(0,k-1)), maxit=maxit, type.gaussian=ifelse(nvars<500,"covariance","naive"), trace.it = trace.it, type.logistic = type.logistic)
	object$x<-x
	object$y<-y
	object$method<-method
    class(object) <- "glmnetcr"
    object
}
