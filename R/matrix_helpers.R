# Internal matrix helpers. No dependence on hetGP's private native symbols.
distance_cpp <- function(X1, X2 = NULL) {
  if (is.null(X2)) X2 <- X1
  out <- matrix(0, nrow(X1), nrow(X2))
  for (j in seq_len(ncol(X1))) out <- out + outer(X1[, j], X2[, j], "-")^2
  out
}
add_diag <- function(A, B) { diag(A) <- diag(A) + B; A }
fast_diag <- function(A, B) rowSums(A * t(B))
trace_sym <- function(A, B) sum(A * B)
fast_tUY2 <- function(mult, Y2) {
  as.vector(rowsum(Y2, rep(seq_along(mult), mult), reorder = FALSE))
}
find_corres <- function(X0, X) {
  vapply(seq_len(nrow(X)), function(i) {
    which(rowSums(sweep(X0, 2, X[i, ], "!=") ) == 0)[1L]
  }, integer(1))
}
# Log-correlation derivative with respect to a common lengthscale.
matern_theta_derivative <- function(X1, X2, theta, nu) {
  out <- matrix(0, nrow(X1), nrow(X2))
  for (j in seq_len(ncol(X1))) {
    r <- sqrt(2 * nu) * abs(outer(X1[, j], X2[, j], "-"))/theta
    out <- out + if (nu == 1.5) r^2/(theta * (1 + r)) else
      r^2 * (1 + r)/(theta * (3 + 3*r + r^2))
  }
  out
}
partial_d_dist_dX1_i1_i2_X2 <- function(X1, X2, i1, i2) {
  out <- matrix(0, nrow(X1), nrow(X2))
  out[i1, ] <- -2 * (X1[i1, i2] - X2[, i2])
  out
}
partial_d_dist_dX_i1_i2 <- function(X1, i1, i2) {
  out <- partial_d_dist_dX1_i1_i2_X2(X1, X1, i1, i2)
  out + t(out)
}
partial_d_dist_abs_dX1_i1_i2_X2 <- function(X1, X2, i1, i2) {
  out <- matrix(0, nrow(X1), nrow(X2))
  delta <- X1[i1, i2] - X2[, i2]
  r <- sqrt(5) * abs(delta)
  out[i1, ] <- -sqrt(5) * sign(delta) * r * (1+r)/(3+3*r+r^2)
  out
}
partial_d_dist_abs_dX_i1_i2 <- function(X1, i1, i2) {
  out <- partial_d_dist_abs_dX1_i1_i2_X2(X1, X1, i1, i2)
  out + t(out)
}
matern5_2_1args <- function(X1) hetGP::cov_gen(X1, theta = 1, type = "Matern5_2")
matern5_2_2args <- function(X1, X2) hetGP::cov_gen(X1, X2, theta = 1, type = "Matern5_2")
d_matern5_2_1args_theta_k <- function(X1, theta) matern_theta_derivative(X1, X1, theta, 2.5)
d_matern5_2_1args_theta_k_iso <- function(X1, theta) matern_theta_derivative(X1, X1, theta, 2.5)
d_matern5_2_2args_theta_k_iso <- function(X1, X2, theta) matern_theta_derivative(X1, X2, theta, 2.5)
d_matern5_2_1args_kthetag <- function(X1, kt) matern_theta_derivative(X1, X1, kt, 2.5)
d_matern5_2_2args_kthetag <- function(X1, X2, kt) matern_theta_derivative(X1, X2, kt, 2.5)
matern3_2_1args <- function(X1) hetGP::cov_gen(X1, theta = 1, type = "Matern3_2")
matern3_2_2args <- function(X1, X2) hetGP::cov_gen(X1, X2, theta = 1, type = "Matern3_2")
d_matern3_2_1args_theta_k <- function(X1, theta) matern_theta_derivative(X1, X1, theta, 1.5)
d_matern3_2_1args_theta_k_iso <- function(X1, theta) matern_theta_derivative(X1, X1, theta, 1.5)
d_matern3_2_2args_theta_k_iso <- function(X1, X2, theta) matern_theta_derivative(X1, X2, theta, 1.5)
d_matern3_2_1args_kthetag <- function(X1, kt) matern_theta_derivative(X1, X1, kt, 1.5)
d_matern3_2_2args_kthetag <- function(X1, X2, kt) matern_theta_derivative(X1, X2, kt, 1.5)
