hmc_exact <- function(L, F, g, mu_r, Sigma, Q, initial_X, T = pi / 2) {
  # Returns samples from a d-dimensional Gaussian with m constraints given
  # by F*X+g>0
  # If cov == true
  # then M is the covariance and the mean is mu = mu_r
  # if cov== false
  # then M is the precision matrix and the log-density is -1/2 X'*M*X + r'*X
  # Input
  # F:          m x d matrix
  # g:          m x 1 vector
  # M           d x d matrix, must be symmmetric and definite positive
  # mu_r        d x 1 vector.
  # cov:        see explanation above
  # L:          number of samples desired
  # initial_X   d x 1 vector. Must satisfy the constraint.
  # Output
  # Xs:      d x L matrix, each column is a sample
  # bounce_count:  number of times the particle bounced

  # go to a whitened frame

  d <- length(initial_X)
  m <- length(g)

  if (!missing(Sigma)) {
    mu <- mu_r
    g <- g + F %*% mu
    R <- chol(Sigma)
    F <- tcrossprod(F, R)
    initial_X <- initial_X - mu
    initial_X <- backsolve(R, initial_X, transpose = TRUE)
  } else {
    mu <- mu_r
    R <- chol(Q)
    g <- g + F %*% mu
    F <- F %*% solve(R)
    initial_X <- initial_X - mu
    initial_X <- R %*% initial_X
  }

  bounce_count <- 0
  nearzero <- 10000 * .Machine$double.eps


  # Verify that initial_X is feasible
  c <- F %*% initial_X + g
  if (any(c < 0)) {
    stop('error: inconsistent initial condition')
  }

  # squared norm of the rows of F, needed for reflecting the velocity
  F2 <- rowSums(F ^ 2)
  Ft <- t(F)

  # Sampling loop
  last_X <- initial_X
  Xs <- matrix(NA, d, L)
  Xs[, 1] <- initial_X

  i <- 2
  while (i <= L) {
    stop <- 0
    j <- 0
    V0 <- rnorm(d)  # initial velocity
    X <- last_X

    tt <- 0  # records how much time the particle already moved

    while (1) {
      a <- V0
      b <- X
      A <- cbind(a, b)
      f <- F %*% A

      U <- sqrt(f[, 1] ^ 2 + f[, 2] ^ 2)
      phi <- atan2(-f[, 1], f[, 2])  # -pi < phi < +pi

      pn <- abs(g / U) <= 1  # these are the walls that may be hit

      # find the first time constraint becomes zero
      if (any(pn)) {
        inds <- which(pn)
        phn <- phi[pn]

        # time at which coordinates hit the walls
        # this expression always gives the correct result because
        # of U*cos(phi+t)+g>=0
        t1 <- -phn + acos(-g[pn] / U[pn])

        # if there was a previous reflection (j>0)
        # and there is a potential reflection at the sample plane
        # make sure that a new reflection at j is not found because of
        # numerical error
        if (j > 0) {
          if (pn[j] == 1) {
            cs <- cumsum(pn)
            indj <- cs[j]
            tt1 <- t1[indj]
            if (abs(tt1) < nearzero || abs(tt1 - 2 * pi) < nearzero) {
              t1[indj] <- Inf
            }
          }
        }

        m_ind <- which.min(t1)
        mt <- t1[m_ind]

        # find the reflection plane
        # j is an index in the full vector of dim-m,
        # not in the restriced vector determined by pn.
        j <- inds[m_ind]

      } else {  # if pn[i] =0 for all i
        mt <- T
      }

      tt <- tt + mt
      if (tt >= T) {
        mt <- mt - (tt - T)
        stop <- 1
      }

      # move the particle a time mt

      X <- a * sin(mt) + b * cos(mt)
      V <- a * cos(mt) - b * sin(mt)

      if (stop) {
        break
      }

      # compute reflected velocity

      qj <- sum(F[j, ] * V) / F2[j]
      V0 <- V - 2 * qj * Ft[, j]
      bounce_count <- bounce_count + 1
    }

    # at this point we have a sampled value X
    Xs[, i] <- X
    last_X <- X
    i <- i + 1
  }

  # transform back to the unwhitened frame
  if (!missing(Sigma)) {
    Xs <- mu + crossprod(R, Xs)
  } else {
    Xs <- mu + backsolve(R, Xs)
  }

  list(Xs = Xs, bounce_count = bounce_count)
}

sample_alpha_truncated <- function(alpha_current, mu, sigma, T = pi / 2) {
  n_alpha <- length(alpha_current)
  hmc_exact(
    2,
    F = diag(n_alpha),
    g = rep(1, n_alpha),
    mu_r = mu,
    Sigma = sigma,
    initial_X = alpha_current,
    T = T
  )$Xs[, 2]
}
