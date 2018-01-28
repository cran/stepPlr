\name{plr-internal}
\title{
  Internal plr functions
}
\alias{anova.stepplr}
\alias{cross.imat}
\alias{get.imat}
\alias{imat}
\alias{print.plr}
\alias{print.stepplr}
\alias{summary.plr}
\alias{summary.stepplr}
\alias{term.match}
\usage{
  \method{anova}{stepplr}(object, \dots)
  cross.imat(imat1, imat2)
  get.imat(x, level = NULL)
  imat(object, x)
  \method{print}{plr}(x, \dots)
  \method{print}{stepplr}(x, \dots)
  \method{summary}{plr}(object, \dots)
  \method{summary}{stepplr}(object, \dots)
  term.match(term1, term2, termlist, len)
}
\description{
  Internal plr functions
}
\references{
  Mee Young Park and Trevor Hastie (2008) Penalized Logistic Regression
  for Detecting Gene Interactions
}
\author{Mee Young Park and Trevor Hastie}
\keyword{internal}
