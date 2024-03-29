% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/logistlogF.R
\name{logistlogF}
\alias{logistlogF}
\title{Logistic regression with log-F(m,m) penalty}
\usage{
logistlogF(form, dat, m, control = glm.control())
}
\arguments{
\item{form}{R formula for the model}

\item{dat}{dataframe of data}

\item{m}{degree-of-freedome parameter}

\item{control}{control convergence}
}
\value{
The fitted penalized logistic regression
}
\description{
Logistic regression with log-F(m,m) penalty
}
\examples{
data(DES); DES$fmatched <- factor(DES$matched.set)
form = formula(case~fmatched+DES+matern.smoke)
fit = logF(form,DES,m=2)
coefficients(fit)
}
