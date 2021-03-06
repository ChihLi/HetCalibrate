% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hetCalibrate.R
\name{rho_AN}
\alias{rho_AN}
\title{# Rho function for SiNK prediction, anistropic case
## @param covtype covariance kernel type, either 'Gaussian' or 'Matern5_2'}
\usage{
rho_AN(
  xx,
  X0,
  theta_g,
  g,
  sigma = 1,
  type = "Gaussian",
  SiNK_eps = 1e-04,
  eps = sqrt(.Machine$double.eps),
  mult
)
}
\description{
# Rho function for SiNK prediction, anistropic case
## @param covtype covariance kernel type, either 'Gaussian' or 'Matern5_2'
}
