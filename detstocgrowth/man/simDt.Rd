% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BayesHybridMod.R
\name{simDt}
\alias{simDt}
\title{Hybrid model expressed as number of cells at time t1, after starting at t0.}
\usage{
simDt(K = 1000, r = 1, N0 = 1, NSwitch = 100, t0 = 0, t1 = 1)
}
\arguments{
\item{K}{- Carrying Capacity}

\item{r}{- Growth Rate}

\item{N0}{- Population Start Size}

\item{t0}{- Initial time point}

\item{t1}{- Time points for which to do inference; single value or vector}

\item{Nswitch}{- Switch from stochastic to deterministic Model}
}
\value{
A vector of cell count(s) for the specified input time(s) is returned.
}
\description{
Hybrid model expressed as number of cells at time t1, after starting at t0.
}

