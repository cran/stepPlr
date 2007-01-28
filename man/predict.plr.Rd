\name{predict.plr}
\alias{predict.plr}
\title{
  prediction function for plr
}
\description{
  This function computes the linear predictors, probability estimates,
  or the class labels for new data, using a \code{plr} object.
}
\usage{
  predict.plr(object, newx = NULL,
              type = c("link", "response", "class"), ...)
}
\arguments{
  \item{object}{
    a \code{plr} object
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
plr
}
\examples{
n <- 100

p <- 10
x0 <- matrix(rnorm(n*p),nrow=n)
y <- sample(c(0,1),n,replace=TRUE)
fit <- plr(x0,y,lambda=1)
x1 <- matrix(rnorm(n*p),nrow=n)
pred1 <- predict(fit,x1,type="link")
pred2 <- predict(fit,x1,type="response")
pred3 <- predict(fit,x1,type="class")

p <- 3
z <- matrix(sample(seq(3),n*p,replace=TRUE),nrow=n)
x0 <- data.frame(x1=factor(z[ ,1]),x2=factor(z[ ,2]),x3=factor(z[ ,3]))
y <- sample(c(0,1),n,replace=TRUE)
fit <- plr(x0,y,lambda=1)
z <- matrix(sample(seq(3),n*p,replace=TRUE),nrow=n)
x1 <- data.frame(x1=factor(z[ ,1]),x2=factor(z[ ,2]),x3=factor(z[ ,3]))
pred1 <- predict(fit,x1,type="link")
pred2 <- predict(fit,x1,type="response")
pred3 <- predict(fit,x1,type="class")
}
\keyword{models}
\keyword{regression}
