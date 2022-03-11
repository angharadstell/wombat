.make_omega_sampler <- function(
  measurement_model,
  process_model,
  X = .make_X_omega(process_model, measurement_model),
  Sigma_epsilon = .make_Sigma_epsilon(measurement_model),
  Xt_Q_epsilon_X = .make_Xt_Q_epsilon_X(
    X,
    measurement_model,
    Sigma_epsilon = Sigma_epsilon
  ),
  mu_omega = .make_mu_omega(process_model, measurement_model),
  Q_omega = .make_Q_omega(process_model, measurement_model),
  chol_Q_omega_conditional = .make_chol_Q_omega_conditional(
    process_model,
    measurement_model,
    X = X,
    Xt_Q_epsilon_X = Xt_Q_epsilon_X,
    Q_omega = Q_omega
  ),
  Z2_tilde = calculate(measurement_model, 'Z2_tilde', process_model)
) {
  omega_unpack <- .make_omega_unpack(process_model, measurement_model)

  function(current) {
    # original WOMBAT paper, Appendix A1
    chol_Q_omega_conditional_i <- chol_Q_omega_conditional(current)
    mu_omega_conditional <- as.vector(.chol_solve(
      chol_Q_omega_conditional_i,
      crossprod(X, solve(Sigma_epsilon(current), Z2_tilde))
      + Q_omega(current) %*% mu_omega(current)
    ))
    # draw new samples for omegas (alphas and betas)
    # ORIG omega <- (
    # ORIG   mu_omega_conditional +
    # ORIG   .sample_normal_precision_chol(chol_Q_omega_conditional_i)
    # ORIG )
    # limit omegas to -1 or greater
    # this will also limit betas as well as alphas!
    # So don't want to use this if solving for betas
    # Arg 1: The current value of the chain
    # Arg 2: The conditional mean
    # Arg 3: The conditional covariance matrix
    # Arg 4: Tuning parameter, see hmc-exact-truncated.R
    omega <- sample_alpha_truncated(
      current$alpha,
      mu = mu_omega_conditional,
      sigma = chol2inv(chol_Q_omega_conditional_i),
      T = pi / 2)

    # check omegas have been correctly adjusted to greater than or equal to -1
    log_debug(paste('Minimum omega value:', min(omega)))
    log_debug(paste('Ratio of omegas more than or equal to -1:',
              sum(omega >= -1) / length(omega)))

    omega_unpack(current, omega)
  }
}
