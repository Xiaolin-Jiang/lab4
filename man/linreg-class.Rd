% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/linear_regression.R
\docType{class}
\name{linreg-class}
\alias{linreg-class}
\alias{linreg}
\title{calculate the statistics using ordinary least squares}
\description{
calculate the statistics using ordinary least squares
}
\section{Fields}{

\describe{
\item{\code{beta:}}{regression coefficient}

\item{\code{y_hat:}}{fitted value}

\item{\code{e_hat}}{: residuals (e with a hat)}

\item{\code{df}}{: degree of fredom}

\item{\code{sigma_square:}}{variance of residuals}

\item{\code{variance_beta:}}{variance of regression coefficient}

\item{\code{t_beta:}}{t value of regression coefficient}

\item{\code{Form:}}{input formula}

\item{\code{data}}{: input dataframe(the data reference of formula)}

\item{\code{dfname:}}{name of dataframe}

\item{\code{Formulaf}}{: string of formula(without bracket)}

\item{\code{X:}}{created by model.matrix()}

\item{\code{y:}}{true value in dataframe not estimated value}
}}


