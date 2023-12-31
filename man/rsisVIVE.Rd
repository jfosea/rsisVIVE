% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rsisVIVE.R
\name{rsisVIVE}
\alias{rsisVIVE}
\title{Robustified Some Valid Some Invalid Instrumental Variable Estimator}
\usage{
rsisVIVE(Y, D, Z, method = c("BIC", "PE_MIN", "PE_SE"), ncores = 1)
}
\arguments{
\item{Y}{matrix array nx1 specifying the outcome values}

\item{D}{matrix array nx1 specifying the exposure values}

\item{Z}{matrix array nxL specifying the instruments where L is the number of instruments}

\item{method}{character string specifying the type of method to obtain the step one (alpha) coefficients. If \code{BIC} is chosen, it uses the robust BIC as specified in the Sparse Least Trimmed Squares method (Alfons et. al., 2013). If \code{PE_MIN} or \code{PE_SE} is chosen, the function uses a 10-fold cross-validation minimizing the Root Trimmed Mean Squared Prediction Error and selecting the minimum or the one standard error, correspondingly. Note that the selection procedure follows the majority rule as outlined in the}

\item{ncores}{numeric indicating the number of cores used for the cross-validation}
}
\value{
description
}
\description{
This function provides the Robustified Some Valid Some Invalid Instrumental Variable Estimator (rsisVIVE) coefficients for step one (alpha) and step two (beta) of the instrumental variable estimation. The algorithm closely follows Kang et. al.'s some invalid some valid instrumental variable estimator (2016) however, the rsisVIVE can tolerate large proportions of contamination both in the exposure and/or outcome variables in addition to tolerating some invalid instrumental variables.
}
\examples{
n <- 50
L <- 10
s <- 3
contam <- 0.2
Sigma <- diag(1, nrow = L)

a0 <- 1
a1 <- rep(0.05, L)

b0 <- 1
b1 <- c(rep(1,s), rep(0, L-s)) # alpha
b2 <- 1                        # beta

error <- MASS::mvrnorm(n, c(0, 0), matrix(c(1, 0.9, 0.9, 1), 2, 2))

Z <- MASS::mvrnorm(n, mu = rep(0, L), Sigma)
D <- a0 + Z \%*\% a1 + error[, 1]
Y <- b0 + Z \%*\% b1 + D * b2 + error[, 2]

rsisVIVE(Y, D, Z, method = 'PE_SE', ncores = 1)
}
