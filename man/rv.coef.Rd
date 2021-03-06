% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{rv.coef}
\alias{rv.coef}
\title{Computes the RV coefficient}
\usage{
rv.coef(S1, S2)
}
\arguments{
\item{S1}{Configuration matrix corresponding to dataset 1}

\item{S2}{Configuration matrix corresponding to dataset 2}
}
\value{
The RV coefficient between dataset 1 and dataset 2
}
\description{
Computes the RV coefficient between dataset 1 and dataset 2. You'll typically won't need to call this function directly,
but should use rv.cor.matrix() instead, as it will make determining partial RV coefficients, p-values and confidence
intervals easier later on.
}
\examples{
set.seed(2)
n = 100
p = 100
x1 = matrix(rnorm(n*p), n, p)
x2 = x1 + matrix(rnorm(n*p), n, p)
S1 = compute.config.matrix(x1)
S2 = compute.config.matrix(x1)
rv.coef(S1, S2)
}
