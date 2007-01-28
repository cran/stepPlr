\name{step.plr}
\alias{step.plr}
\title{
  Forward stepwise selection procedure for penalized logistic
  regression
}
\description{
  This function fits a series of L2 penalized logistic regression models
  selecting variables through the forward stepwise selection procedure.
}
\usage{
  step.plr(x, y, weights = rep(1,length(y)), fix.subset = NULL,
           level = NULL, lambda = 1e-4, cp = "bic",
           max.terms = 5, type = c("both", "forward"), trace = FALSE)  
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
  \item{fix.subset}{
    a vector of indices for the variables that are forced to be in the
    model
  }
  \item{level}{
    a list of length \code{ncol(x).} The j-th element corresponds to
    the j-th column of \code{x.} If the j-th column of \code{x} is
    discrete, \code{level[[j]]} is the set of levels for the
    categorical factor. If the j-th column of \code{x} is continuous,
    \code{level[[j]] = NULL.} \code{level} is automatically
    generated in the function; however, if any levels of the
    categorical factors are not observed, but still need to be included
    in the model, then the user must provide the complete sets of the
    levels through \code{level.} If a numeric column needs to be
    considered discrete, it can be done by manually providing
    \code{level} as well.
  }
  \item{lambda}{
    regularization parameter for the L2 norm of the
    coefficients. The minimizing criterion in \code{plr} is
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
  \item{max.terms}{
    the maximum number of terms to be added in the forward selection
    procedure. Default is \code{max.terms=5.}
  }
  \item{type}{
    If \code{type="both",} the forward selection is followed by a
    backward deletion. If \code{type="forward",} only a forward
    selection is done. Default is \code{"both".}
  }
  \item{trace}{
    If \code{TRUE,} the variable selection procedure prints out its
    progress.
  }
}
\value{
  A \code{stepplr} object is returned. \code{anova, predict, print,} and
  \code{summary} functions can be applied.
  \item{fit}{
    a \code{plr} object for the optimal model selected
  }
  \item{action}{
    a list that stores the selection order of the terms in the optimal
    model. 
  }
  \item{action.name}{
    a list of the names of the sequentially added terms - in the same
    order as in \code{action}
  }
  \item{deviance}{
    deviance of the fitted model
  }
  \item{df}{
    residual degrees of freedom of the fitted model
  }
  \item{score}{
    deviance + cp*df, where df is the model degrees of freedom
  }
  \item{group}{
    a vector of the counts for the dummy variables, to be used in
    \code{predict.stepplr}
  }
  \item{y}{
    response variable used
  }
  \item{weight}{
    weights used
  }
  \item{fix.subset}{
    fix.subset used
  }
  \item{level}{
    level used
  }
  \item{lambda}{
    lambda used
  }
  \item{cp}{
    complexity parameter used when computing the score
  }
  \item{type}{
    type used
  }
  \item{xnames}{
    column names of \code{x}
  }
}
\details{
  This function implements an L2 penalized logistic regression along
  with the stepwise variable selection procedure, as described in
  "Penalized Logistic Regression for Detecting Gene Interactions (2006)"
  by Park and Hastie.

  If \code{type="forward",} \code{max.terms} terms are sequentially
  added to the model, and the model that minimizes \code{score} is
  selected as the optimal fit. If \code{type="both",} a backward
  deletion is done in addition, which provides a series of models with a
  different combination of the selected terms. The optimal model
  minimizing \code{score} is chosen from the second list.
  
  We thank Michael Saunders of SOL, Stanford University for providing
  the solver used for the convex optimization in this function.
}
\references{
  Mee Young Park and Trevor Hastie (2006) Penalized Logistic Regression
  for Detecting Gene Interactions - available at the authors' websites,
  \url{http://stat.stanford.edu/~mypark} or
  \url{http://stat.stanford.edu/~hastie/pub.htm}.
}
\author{Mee Young Park and Trevor Hastie}
\seealso{
cv.step.plr, plr, predict.stepplr
}
\examples{
n <- 100

p <- 3
z <- matrix(sample(seq(3),n*p,replace=TRUE),nrow=n)
x <- data.frame(x1=factor(z[ ,1]),x2=factor(z[ ,2]),x3=factor(z[ ,3]))
y <- sample(c(0,1),n,replace=TRUE)
fit <- step.plr(x,y)
# 'level' is automatically generated. Check 'fit$level'.

p <- 5
x <- matrix(sample(seq(3),n*p,replace=TRUE),nrow=n)
x <- cbind(rnorm(n),x)
y <- sample(c(0,1),n,replace=TRUE)
level <- vector("list",length=6)
for (i in 2:6) level[[i]] <- seq(3)
fit1 <- step.plr(x,y,level=level,cp="aic")
fit2 <- step.plr(x,y,level=level,cp=4)
fit3 <- step.plr(x,y,level=level,type="forward")
fit4 <- step.plr(x,y,level=level,max.terms=10)
# This is an example in which 'level' was input manually.
# level[[1]] should be either 'NULL' or 'NA' since the first factor is continuous.
}
\keyword{models}
\keyword{regression}
