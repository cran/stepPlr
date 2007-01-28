\name{predict.stepplr}
\alias{predict.stepplr}
\title{
  prediction function for step.plr
}
\description{
  This function computes the linear predictors, probability estimates,
  or the class labels for new data, using a \code{stepplr} object.
}
\usage{
  predict.stepplr(object, x = NULL, newx = NULL,
                  type = c("link", "response", "class"), ...)
}
\arguments{
  \item{object}{
    a \code{stepplr} object
  }
  \item{x}{
    the matrix of features used for fitting \code{object.} If
    \code{newx} is provided, \code{x} must be provided as well.
  }
  \item{newx}{
    a matrix of features at which the predictions are made. If
    \code{newx=NULL,} predictions for the training data are returned.
  }
  \item{type}{
    If \code{type=link,} the linear predictors are returned; if
    \code{type=response,} the probability estimates are returned; and if
    \code{type=class,} the class labels are returned. Default is
    \code{type=link.}
  }
  \item{...}{
    other options for the prediction
  }
}
\references{
  Mee Young Park and Trevor Hastie (2006) Penalized Logistic Regression
  for Detecting Gene Interactions - available at the authors' websites,
  \url{http://stat.stanford.edu/~mypark} or
  \url{http://stat.stanford.edu/~hastie/pub.htm}.
}
\author{Mee Young Park and Trevor Hastie}
\seealso{
stepplr
}
\examples{
n <- 100
p <- 5
x0 <- matrix(sample(seq(3),n*p,replace=TRUE),nrow=n)
x0 <- cbind(rnorm(n),x0)
y <- sample(c(0,1),n,replace=TRUE)
level <- vector("list",length=6)
for (i in 2:6) level[[i]] <- seq(3)
fit <- step.plr(x0,y,level=level)
x1 <- matrix(sample(seq(3),n*p,replace=TRUE),nrow=n)
x1 <- cbind(rnorm(n),x1)
pred1 <- predict(fit,x0,x1,type="link")
pred2 <- predict(fit,x0,x1,type="response")
pred3 <- predict(fit,x0,x1,type="class")
}
\keyword{models}
\keyword{regression}
