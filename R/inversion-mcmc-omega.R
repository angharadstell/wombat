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
    # page 38, Appendix A1
    # Q_epsilon is the same as the paper? X is the basis functions
    # chol_Q_omega_conditional_i is the choleski factorisation of (Q_alpha + t(basis) * Q_epsilon * basis) - aka measure of uncertainty on alpha 
    chol_Q_omega_conditional_i <- chol_Q_omega_conditional(current)
    # mu_omega_conditional is the mean of alpha
    # this isnt quite the same as the paper?
    # what actually is Z2_tilde?
    # X isn't transposed, Q_omega(current) %*% mu_omega(current) isn't in paper
    mu_omega_conditional <- as.vector(.chol_solve(
      chol_Q_omega_conditional_i,
      crossprod(X, solve(Sigma_epsilon(current), Z2_tilde))
      + Q_omega(current) %*% mu_omega(current)
    ))
    # new sample is from gaussian distrib with mean mu_omega_conditional, choleski chol_Q_omega_conditional_i
    omega <- (
      mu_omega_conditional
      + .sample_normal_precision_chol(chol_Q_omega_conditional_i)
    )
    # attempt to limit at -1
    # omega <- TruncatedNormal::rtmvnorm(n = 1,
    #                                    mu = mu_omega_conditional,
    #                                    sigma = chol2inv(chol_Q_omega_conditional_i),
    #                                    lb = rep(-1, ncol(chol_Q_omega_conditional_i)))

    # check omegas have been correctly adjusted to greater than or equal to -1
    log_debug(paste("Minimum omega value:", min(omega)))
    log_debug(paste("Ratio of omegas more than or equal to -1:", sum(omega >= -1) / length(omega)))

    omega_unpack(current, omega)
  }
}
