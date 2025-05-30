% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/desmon-extra.R
\name{bin1samp_power}
\alias{bin1samp_power}
\alias{bin1samp_sim}
\title{Power for one-sample exact binomial designs}
\usage{
bin1samp_power(p0, pa, n, r)

bin1samp_sim(p0, pa, n, r = seq.int(n), plot = TRUE)
}
\arguments{
\item{p0, pa}{probability of success under the null and alternative
hypotheses, respectively}

\item{n}{sample size, a single value or vector}

\item{r}{a vector of critical values, typically the minimum number of
successes required to reject \code{pa}}

\item{plot}{logical; if \code{TRUE}, a the sequence of \code{r} versus
type-I and type-II errors is plotted}
}
\value{
\code{bin1samp_power} returns a vector or matrix with the following:

\item{\code{type1}}{the overall type-I error}
\item{\code{type2}}{the overall type-II error}

\code{bin1samp_sim} returns a data frame with the following columns:

\item{\code{r}}{critical values}
\item{\code{type1}}{the overall type-I errors}
\item{\code{type2}}{the overall type-II errors}
}
\description{
Determines the power and significance level for a one-sided, one-sample
exact binomial test.

\code{bin1samp_power} accepts a single value or vector for \code{n} and/or
\code{r} and will return a matrix with results for each combination. If
only one value is given for each, a vector is returned.
}
\examples{
p0 <- 0.1
pa <- 0.3
des <- desmon2:::bin1samp(p0, pa)
bin1samp_power(p0, pa, des['n'], des['r'] + 1)
## compare
des[c('size', 'type2')]

bin1samp_power(p0, pa, des['n'], des['r'] + -2:2)
bin1samp_power(p0, pa, des['n'] + 0:1, des['r'] + -2:2)

## simulate over critical values
bin1samp_sim(p0, pa, des['n'])

des <- desmon2:::bin1samp(pa, p0)
bin1samp_power(pa, p0, des['n'], des['r'])
## compare
des[c('size', 'type2')]

}
\seealso{
\code{\link{twostg_power}}; \code{\link{twostg_sim}}
}
