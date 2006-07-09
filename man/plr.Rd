\name{plr}
\alias{plr}
\title{
  Logistic regression with a quadratic penalization on the coefficients
}
\description{
  This function fits a logistic regression model penalizing the size of
  the L2 norm of the coefficients.
}
\usage{
  plr(x, y, weights = rep(1,length(y)), lambda = 1e-4, cp = "bic")
}
\arguments{
  \item{x}{
    matrix of features
  }
  \item{y}{
    binary response 
  }
  \item{weights}{
    an optional vector of weights for observations
  }
  \item{lambda}{
    regularization parameter for the L2 norm of the coefficients. The
    minimizing criterion in \code{plr} is
    -log-likelihood+\eqn{\lambda*\|\beta\|^2}. Default is
    \code{lambda=1e-4.}
  }
  \item{cp}{
    complexity parameter to be used when computing the
    score. \code{score=deviance+cp*df.} If \code{cp="aic"} or
    \code{cp="bic",} these are converted to \code{cp=2} and
    \code{cp=log(sample size),} respectively. Default is
    \code{cp="bic".}
  }
}
\value{
  A \code{plr} object is returned. \code{predict, print,} and
  \code{summary} functions can be applied.  
  \item{coefficients}{
    vector of the coefficient estimates
  }
  \item{covariance}{
    a sandwich estimate of the covariance matrix for the coefficients
  }
  \item{deviance}{
    deviance of the fitted model
  }
  \item{null.deviance}{
    deviance of the null model
  }
  \item{df}{
    degrees of freedom of the fitted model
  }
  \item{score}{
    deviance + cp*df
  }
  \item{nobs}{
    number of observations
  }
  \item{cp}{
    complexity parameter used when computing the score
  }
  \item{fitted.values}{
    fitted probabilities
  }
  \item{linear.predictors}{
    linear predictors computed with the estimated coefficients
  }
  \item{level}{
    If any categorical factors are input, level - the list of level sets
    - is automatically generated and returned. See \code{step.plr} for
    details of how it is generated.
  }
}
\details{
  We proposed using logistic regression with a quadratic penalization on
  the coefficients for detecting gene interactions as described in
  "Penalized Logistic Regression for Detecting Gene Interactions (2006)"
  by Park and Hastie. However, this function \code{plr} may be used for
  a general purpose. 
    
  We thank Michael Saunders of SOL, Stanford University for providing
  the solver used for the convex optimization in this function.
}
\references{
  Mee Young Park and Trevor Hastie (2006) Penalized Logistic Regression
  for Detecting Gene Interactions - available at the authors' websites,
  \url{http://www.stanford.edu/~mypark} or
  \url{http://stat.stanford.edu/~hastie/pub.htm}.
}
\author{Mee Young Park and Trevor Hastie}
\seealso{
predict.plr, step.plr
}
\examples{
n <- 100

p <- 10
x <- matrix(rnorm(n*p),nrow=n)
y <- sample(c(0,1),n,replace=TRUE)
fit <- plr(x,y,lambda=1)

p <- 3
z <- matrix(sample(seq(3),n*p,replace=TRUE),nrow=n)
x <- data.frame(x1=factor(z[ ,1]),x2=factor(z[ ,2]),x3=factor(z[ ,3]))
y <- sample(c(0,1),n,replace=TRUE)
fit <- plr(x,y,lambda=1)
# 'level' is automatically generated. Check 'fit$level'.
}
\keyword{models}
\keyword{regression}
