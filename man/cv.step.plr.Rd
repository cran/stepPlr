\name{cv.step.plr}
\alias{cv.step.plr}
\title{
  Computes cross-validated deviance or prediction errors for step.plr
}
\description{  
  This function computes cross-validated deviance or prediction errors
  for \code{step.plr.} The parameters that can be cross-validated are
  \code{lambda} and \code{cp}. 
}
\usage{
  cv.step.plr(x, y, weights = rep(1, length(y)), nfold = 5,
              folds = NULL, lambda = c(1e-4, 1e-2, 1),
              cp = c("aic", "bic"), cv.type=c("deviance", "class"),
              trace = TRUE, ...)
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
  \item{nfold}{
    number of folds to be used in cross-validation. Default is
    \code{nfold=5.}
  }
  \item{folds}{
    the list of cross-validation folds. Its length must be nfold. If
    \code{NULL,} the folds are randomly generated.
  }
  \item{lambda}{
    vector of the candidate values for \code{lambda} in \code{step.plr} 
  }
  \item{cp}{
    vector of the candidate values for \code{cp} in \code{step.plr} 
  }
  \item{cv.type}{
    If \code{cv.type=deviance,} cross-validated deviances are returned. If
    \code{cv.type=class,} cross-validated prediction errors are returned.
  }
  \item{trace}{
    If \code{TRUE,} the steps are printed out.
  }
  \item{...}{
    other options for \code{step.plr}
  }
}
\details{
  This function computes cross-validated deviance or prediction errors
  for \code{step.plr.} The parameters that can be cross-validated are
  \code{lambda} and \code{cp}. If both are input as vectors (of length
  greater than 1), then a two-dimensional cross-validation is done. If
  either one is input as a single value, then the cross-validation is
  done only on the parameter with multiple inputs.
}
\references{
  Mee Young Park and Trevor Hastie (2006) Penalized Logistic Regression
  for Detecting Gene Interactions - available at the authors' websites,
  \url{http://www.stanford.edu/~mypark} or
  \url{http://stat.stanford.edu/~hastie/pub.htm}.
}
\author{Mee Young Park and Trevor Hastie}
\seealso{
step.plr
}
\examples{
n <- 100
p <- 5
x <- matrix(sample(seq(3),n*p,replace=TRUE),nrow=n)
y <- sample(c(0,1),n,replace=TRUE)
level <- vector("list",length=p)
for (i in 1:p) level[[i]] <- seq(3)
cvfit1 <- cv.step.plr(x,y,level=level,lambda=c(1e-4,1e-2,1),cp="bic")
cvfit2 <- cv.step.plr(x,y,level=level,lambda=1e-4,cp=c(2,3,4))
cvfit3 <- cv.step.plr(x,y,level=level,lambda=c(1e-4,1e-2,1),cp=c(2,3,4))
}
\keyword{models}
\keyword{regression}
