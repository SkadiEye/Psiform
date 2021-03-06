% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/2-1-psiform.R
\name{psiform}
\alias{psiform}
\title{A tail-based test statistics (TTS) for differential expression detection.}
\usage{
psiform(
  ...,
  y,
  G = NULL,
  dim_u = 5,
  weight = "none",
  scale = FALSE,
  lambda = NULL,
  a = 3.7,
  param_init = NULL,
  n_iter = 100,
  epsilon = 10^-8,
  conv_thresh = 10^-6,
  method = 1
)
}
\arguments{
\item{...}{One or multiple matrices of expression levels from multiple platforms.}

\item{y}{Subject labels.}

\item{G}{Pathway information.}

\item{dim_u}{Number of dimension of u.}

\item{weight}{Type of weight on the pathway. Possible values include "sqrt-degree",
"none" or "degree". The default is "none".}

\item{scale}{Whether to scale input matrices.}

\item{lambda}{Hyper-parameters in SCAD, TLP, and L1 of phenotype differences.}

\item{a}{Hyper-parameter in SCAD. Default is 3.7.}

\item{param_init}{If not NULL, should be a list including c, beta, u, v and sigma_sq.}

\item{n_iter}{Number of maximum iterations.}

\item{epsilon}{A small number to avoid zero-denominators.}

\item{conv_thresh}{Threshold for convergence.}

\item{method}{Type of TLP.}
}
\value{
Returns a \code{Psiform} object.
}
\description{
A tail-based test statistics (TTS) for differential expression detection.
}
