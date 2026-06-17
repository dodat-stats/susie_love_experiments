############################################################
## SuSiE-style CAVI + VEB for inverse-Wishart column model
## with summary-statistic likelihood covariance Rbar / N
############################################################

logsumexp <- function(x) {
  m <- max(x)
  m + log(sum(exp(x - m)))
}

safe_normalize_logweights <- function(logw) {
  lse <- logsumexp(logw)
  exp(logw - lse)
}

chol_solve <- function(ch, b) {
  # ch is upper-triangular Cholesky factor from chol(A), so A = t(ch) %*% ch
  backsolve(ch, forwardsolve(t(ch), b))
}

quad_form_chol <- function(ch, x) {
  y <- chol_solve(ch, x)
  sum(x * y)
}

############################################################
## One-dimensional quadrature for q(s)
##
## q(s) proportional to
## s^{-(a+1)} exp{-b/s + D s - 0.5 H s^2}, s > 0.
##
## With s = exp(xi),
## q(xi) proportional to
## exp{-a xi - b exp(-xi) + D exp(xi) - 0.5 H exp(2 xi)}.
############################################################

qs_moments_grid <- function(a, b, D, H,
                            grid_n = 121,
                            width = 12,
                            mode_lower = -40,
                            mode_upper = 40) {
  H <- max(H, 1e-14)
  b <- max(b, 1e-14)

  logf <- function(xi) {
    -a * xi - b * exp(-xi) + D * exp(xi) - 0.5 * H * exp(2 * xi)
  }

  opt <- optimize(logf, interval = c(mode_lower, mode_upper), maximum = TRUE)
  mode <- opt$maximum

  # Expand grid if boundary mass is non-negligible
  for (attempt in 1:5) {
    grid <- seq(mode - width, mode + width, length.out = grid_n)
    lf <- logf(grid)
    m <- max(lf)
    w <- exp(lf - m)

    # trapezoid weights
    wt <- w
    wt[1] <- wt[1] / 2
    wt[length(wt)] <- wt[length(wt)] / 2

    p <- wt / sum(wt)

    boundary_mass <- p[1] + p[length(p)]
    if (boundary_mass < 1e-7) break
    width <- width * 1.5
  }

  dx <- grid[2] - grid[1]
  logZ <- m + log(sum(wt) * dx)

  Es     <- sum(exp(grid) * p)
  Es2    <- sum(exp(2 * grid) * p)
  ElogS  <- sum(grid * p)
  EInvS  <- sum(exp(-grid) * p)

  list(
    Es = Es,
    Es2 = Es2,
    ElogS = ElogS,
    EInvS = EInvS,
    logZ = logZ,
    mode = mode,
    width = width
  )
}

############################################################
## Main fitting function
############################################################

iw_susie_veb <- function(x,
                         Rbar,
                         L,
                         N = 1,
                         sigma2 = rep(1, L),
                         nu0_init = 10,
                         max_iter = 50,
                         inner_iter = 3,
                         update_nu0 = TRUE,
                         update_sigma2 = FALSE,
                         nu0_lower = 2,
                         ridge = 1e-8,
                         tol = 1e-5,
                         qs_grid_n = 121,
                         init_gamma = NULL,
                         init_beta = NULL,
                         verbose = TRUE) {

  x <- as.numeric(x)
  J <- length(x)

  stopifnot(nrow(Rbar) == J, ncol(Rbar) == J)
  stopifnot(length(N) == 1, is.finite(N), N > 0)
  stopifnot(length(sigma2) == L)
  if (!is.null(init_gamma)) {
    stopifnot(length(init_gamma) == L)
    stopifnot(all(init_gamma == as.integer(init_gamma)))
    stopifnot(all(init_gamma >= 1), all(init_gamma <= J))
  }
  if (!is.null(init_beta)) {
    stopifnot(length(init_beta) == L)
    stopifnot(all(is.finite(init_beta)))
  }

  # Ridge for numerical stability
  Rbar <- (Rbar + t(Rbar)) / 2
  if (ridge > 0) {
    Rbar <- Rbar + ridge * diag(J)
  }

  rho <- diag(Rbar)

  if (any(rho <= 0)) {
    stop("All diagonal entries of Rbar must be positive.")
  }

  chR <- chol(Rbar)
  logdetR <- 2 * sum(log(diag(chR)))

  # G[, j] = Rbar[, j] / Rbar[j, j]
  # This is E[u_j(z)] at the prior mean.
  G <- sweep(Rbar, 2, rho, "/")

  logdetQ <- logdetR - log(rho)

  # Initialize nu0
  nu0 <- nu0_init

  # Helper to get prior moments for s under current nu0
  prior_s_moments <- function(nu0) {
    a <- (nu0 + 2) / 2
    b <- nu0 * rho / 2

    Es <- b / (a - 1)

    if (a <= 2) {
      Es2 <- rep(Inf, J)
    } else {
      Es2 <- b^2 / ((a - 1)^2 * (a - 2))
    }

    ElogS <- log(b) - digamma(a)
    EInvS <- a / b

    list(a = a, b = b, Es = Es, Es2 = Es2,
         ElogS = ElogS, EInvS = EInvS)
  }

  smom <- prior_s_moments(nu0)

  # Variational parameters
  pi_mat <- matrix(1 / J, nrow = L, ncol = J)

  mb <- matrix(0, nrow = L, ncol = J)
  vb <- matrix(rep(sigma2, each = J), nrow = L, byrow = TRUE)

  Es <- matrix(rep(smom$Es, times = L), nrow = L, byrow = TRUE)
  Es2 <- matrix(rep(smom$Es2, times = L), nrow = L, byrow = TRUE)
  ElogS <- matrix(rep(smom$ElogS, times = L), nrow = L, byrow = TRUE)
  EInvS <- matrix(rep(smom$EInvS, times = L), nrow = L, byrow = TRUE)

  # q_s natural parameters, needed to evaluate entropy of q(s)
  qsa <- matrix(smom$a, nrow = L, ncol = J)
  qsb <- matrix(rep(smom$b, times = L), nrow = L, byrow = TRUE)

  D_s <- matrix(0, nrow = L, ncol = J)
  H_s <- matrix(0, nrow = L, ncol = J)

  # log normalizer of q_s
  # If initialized at prior, Z = Gamma(a) / b^a
  logZs <- matrix(
    rep(lgamma(smom$a) - smom$a * log(smom$b), times = L),
    nrow = L,
    byrow = TRUE
  )

  # omega moments
  c0 <- (nu0 + 3) / 2
  Eomega <- matrix(nu0 + 3, nrow = L, ncol = J)
  Elogomega <- matrix(digamma(c0) + log(2), nrow = L, ncol = J)

  # q_omega parameters
  aw <- matrix(c0, nrow = L, ncol = J)
  bw <- matrix(0.5, nrow = L, ncol = J)

  # z-related scalar summaries
  Cz <- matrix(NA_real_, nrow = L, ncol = J)
  alpha_z <- matrix(0, nrow = L, ncol = J)
  A_z <- matrix(NA_real_, nrow = L, ncol = J)
  B_z <- matrix(NA_real_, nrow = L, ncol = J)
  Delta_z <- matrix(NA_real_, nrow = L, ncol = J)

  # Effect posterior means m_l = E[beta_l R_{gamma_l}]
  m_eff <- matrix(0, nrow = J, ncol = L)

  if (!is.null(init_gamma)) {
    init_beta <- if (is.null(init_beta)) rep(0, L) else as.numeric(init_beta)
    for (l in seq_len(L)) {
      j <- init_gamma[l]
      pi_mat[l, ] <- 1e-12 / (J - 1)
      pi_mat[l, j] <- 1 - 1e-12
      mb[l, j] <- init_beta[l]
      m_eff[, l] <- init_beta[l] * Rbar[, j]
    }
  }

  ############################################################
  ## Helper: refresh z scalar summaries for one effect
  ############################################################

  refresh_z_scalars <- function(l, y) {
    yKy <- quad_form_chol(chR, y)

    # d_j(y) = y^T K y - y_j^2 / rho_j
    d <- yKy - y^2 / rho
    d <- pmax(d, 0)

    Eb2 <- vb[l, ] + mb[l, ]^2
    T1 <- mb[l, ] * Es[l, ]
    T2 <- Eb2 * Es2[l, ]

    Cz_l <- N * T2 + rho * Eomega[l, ]
    Cz_l <- pmax(Cz_l, 1e-14)

    alpha_l <- N * T1 / Cz_l

    B_l <- y / rho + alpha_l * d

    A_l <- 1 / rho + alpha_l^2 * d + (J - 1) / Cz_l

    Delta_l <- rho * (alpha_l^2 * d + (J - 1) / Cz_l)

    list(
      yKy = yKy,
      d = d,
      Eb2 = Eb2,
      T1 = T1,
      T2 = T2,
      Cz = Cz_l,
      alpha = alpha_l,
      A = A_l,
      B = B_l,
      Delta = Delta_l
    )
  }

  ############################################################
  ## Helper: prior-q terms for each candidate j
  ############################################################

  local_prior_minus_q <- function(l, sigma2_l, nu0_current) {
    a_cur <- (nu0_current + 2) / 2
    b_cur <- nu0_current * rho / 2
    c_cur <- (nu0_current + 3) / 2

    Eb2 <- vb[l, ] + mb[l, ]^2

    ## beta
    Elogp_beta <- -0.5 * (log(2 * pi * sigma2_l) + Eb2 / sigma2_l)
    Elogq_beta <- -0.5 * (log(2 * pi * vb[l, ]) + 1)

    ## s prior under current nu0
    Elogp_s <- a_cur * log(b_cur) -
      lgamma(a_cur) -
      (a_cur + 1) * ElogS[l, ] -
      b_cur * EInvS[l, ]

    ## q(s), using stored q_s natural parameters
    Elogq_s <- -(qsa[l, ] + 1) * ElogS[l, ] -
      qsb[l, ] * EInvS[l, ] +
      D_s[l, ] * Es[l, ] -
      0.5 * H_s[l, ] * Es2[l, ] -
      logZs[l, ]

    ## omega prior under current nu0
    Elogp_omega <- c_cur * log(0.5) -
      lgamma(c_cur) +
      (c_cur - 1) * Elogomega[l, ] -
      0.5 * Eomega[l, ]

    ## q(omega)
    Elogq_omega <- aw[l, ] * log(bw[l, ]) -
      lgamma(aw[l, ]) +
      (aw[l, ] - 1) * Elogomega[l, ] -
      bw[l, ] * Eomega[l, ]

    ## z prior
    p_dim <- J - 1

    Elogp_z <- -0.5 * p_dim * log(2 * pi) -
      0.5 * logdetQ +
      0.5 * p_dim * (Elogomega[l, ] + log(rho)) -
      0.5 * Eomega[l, ] * Delta_z[l, ]

    ## q(z), where V_z = Q_j / C_z
    Elogq_z <- -0.5 * (
      p_dim * (log(2 * pi) + 1) +
        logdetQ -
        p_dim * log(Cz[l, ])
    )

    out <- (Elogp_beta - Elogq_beta) +
      (Elogp_s - Elogq_s) +
      (Elogp_omega - Elogq_omega) +
      (Elogp_z - Elogq_z)

    out
  }

  ############################################################
  ## Helper: update one SuSiE effect
  ############################################################

  update_effect <- function(l, y) {
    a_cur <- (nu0 + 2) / 2
    b_cur <- nu0 * rho / 2

    for (inner in seq_len(inner_iter)) {
      ## 1. Refresh z summaries
      zsum <- refresh_z_scalars(l, y)

      Cz[l, ] <- zsum$Cz
      alpha_z[l, ] <- zsum$alpha
      A_z[l, ] <- zsum$A
      B_z[l, ] <- zsum$B
      Delta_z[l, ] <- zsum$Delta

      ## 2. Update beta
      vb[l, ] <<- 1 / (1 / sigma2[l] + N * Es2[l, ] * A_z[l, ])
      mb[l, ] <<- vb[l, ] * N * Es[l, ] * B_z[l, ]

      ## 3. Refresh z summaries after beta update
      zsum <- refresh_z_scalars(l, y)

      Cz[l, ] <<- zsum$Cz
      alpha_z[l, ] <<- zsum$alpha
      A_z[l, ] <<- zsum$A
      B_z[l, ] <<- zsum$B
      Delta_z[l, ] <<- zsum$Delta

      ## 4. Update s by one-dimensional quadrature
      Eb2 <- vb[l, ] + mb[l, ]^2
      D_vec <- N * mb[l, ] * B_z[l, ]
      H_vec <- N * Eb2 * A_z[l, ]

      for (j in seq_len(J)) {
        qs <- qs_moments_grid(
          a = a_cur,
          b = b_cur[j],
          D = D_vec[j],
          H = H_vec[j],
          grid_n = qs_grid_n
        )

        Es[l, j] <<- qs$Es
        Es2[l, j] <<- qs$Es2
        ElogS[l, j] <<- qs$ElogS
        EInvS[l, j] <<- qs$EInvS
        logZs[l, j] <<- qs$logZ

        qsa[l, j] <<- a_cur
        qsb[l, j] <<- b_cur[j]
        D_s[l, j] <<- D_vec[j]
        H_s[l, j] <<- H_vec[j]
      }

      ## 5. Refresh z summaries after s update
      zsum <- refresh_z_scalars(l, y)

      Cz[l, ] <<- zsum$Cz
      alpha_z[l, ] <<- zsum$alpha
      A_z[l, ] <<- zsum$A
      B_z[l, ] <<- zsum$B
      Delta_z[l, ] <<- zsum$Delta

      ## 6. Update omega
      aw_new <- (nu0 + J + 2) / 2
      bw_new <- (1 + Delta_z[l, ]) / 2

      aw[l, ] <<- aw_new
      bw[l, ] <<- bw_new

      Eomega[l, ] <<- aw[l, ] / bw[l, ]
      Elogomega[l, ] <<- digamma(aw[l, ]) - log(bw[l, ])
    }

    ## Final refresh before scoring
    zsum <- refresh_z_scalars(l, y)

    Cz[l, ] <<- zsum$Cz
    alpha_z[l, ] <<- zsum$alpha
    A_z[l, ] <<- zsum$A
    B_z[l, ] <<- zsum$B
    Delta_z[l, ] <<- zsum$Delta

    Eb2 <- vb[l, ] + mb[l, ]^2
    T1 <- mb[l, ] * Es[l, ]
    T2 <- Eb2 * Es2[l, ]

    ## Likelihood part, dropping constants common over j:
    ## E log p(y | theta_l) = N * (T1 B - 0.5 T2 A) + const.
    lik_score <- N * (T1 * B_z[l, ] - 0.5 * T2 * A_z[l, ])

    prior_score <- local_prior_minus_q(l, sigma2[l], nu0)

    logw <- lik_score + prior_score

    pi_mat[l, ] <<- safe_normalize_logweights(logw)

    ## Update effect mean m_l.
    ##
    ## E[u_j(z)] = g_j + alpha_j (y - y_j g_j)
    ## where g_j = Rbar[, j] / rho_j.
    ##
    ## Sum_j w_j E[u_j] with w_j = pi_lj E[beta] E[s].
    w <- pi_mat[l, ] * mb[l, ] * Es[l, ]

    coef_g <- w * (1 - alpha_z[l, ] * y)
    coef_y <- sum(w * alpha_z[l, ])

    m_eff[, l] <<- as.numeric(G %*% coef_g) + coef_y * y

    invisible(NULL)
  }

  ############################################################
  ## ELBO
  ############################################################

  compute_elbo <- function() {
    mu_total <- rowSums(m_eff)
    resid <- x - mu_total

    resid_quad <- quad_form_chol(chR, resid)

    var_quad <- 0

    for (l in seq_len(L)) {
      Eb2 <- vb[l, ] + mb[l, ]^2

      # E[theta_l^T K theta_l]
      E_quad_l <- sum(pi_mat[l, ] * Eb2 * Es2[l, ] * A_z[l, ])

      # m_l^T K m_l
      m_quad_l <- quad_form_chol(chR, m_eff[, l])

      var_quad <- var_quad + E_quad_l - m_quad_l
    }

    var_quad <- max(var_quad, 0)

    Eloglik <- -0.5 * (
      J * log(2 * pi) +
        logdetR -
        J * log(N) +
        N * (resid_quad + var_quad)
    )

    Eprior_minus_q <- 0

    for (l in seq_len(L)) {
      lpq <- local_prior_minus_q(l, sigma2[l], nu0)

      gamma_term <- log(1 / J) - log(pmax(pi_mat[l, ], 1e-300))

      Eprior_minus_q <- Eprior_minus_q +
        sum(pi_mat[l, ] * (gamma_term + lpq))
    }

    Eloglik + Eprior_minus_q
  }

  ############################################################
  ## VEB update for nu0
  ############################################################

  M_nu0 <- function(nu0_candidate) {
    if (nu0_candidate <= nu0_lower) return(-Inf)

    a_cand <- (nu0_candidate + 2) / 2
    b_cand <- nu0_candidate * rho / 2

    c_cand <- (nu0_candidate + 3) / 2

    total <- 0

    for (l in seq_len(L)) {
      Elogp_s <- a_cand * log(b_cand) -
        lgamma(a_cand) -
        (a_cand + 1) * ElogS[l, ] -
        b_cand * EInvS[l, ]

      Elogp_omega <- c_cand * log(0.5) -
        lgamma(c_cand) +
        (c_cand - 1) * Elogomega[l, ] -
        0.5 * Eomega[l, ]

      total <- total + sum(pi_mat[l, ] * (Elogp_s + Elogp_omega))
    }

    total
  }

  update_nu0_EB <- function() {
    # nu0 = nu0_lower + exp(eta)
    obj <- function(eta) {
      nu <- nu0_lower + exp(eta)
      -M_nu0(nu)
    }

    opt <- optimize(
      obj,
      interval = c(log(1e-8), log(1e4))
    )

    nu0_lower + exp(opt$minimum)
  }

  ############################################################
  ## Main loop
  ############################################################

  elbo_trace <- numeric(max_iter)

  for (iter in seq_len(max_iter)) {
    old_elbo <- if (iter == 1) -Inf else elbo_trace[iter - 1]

    for (l in seq_len(L)) {
      current_sum <- rowSums(m_eff)
      y_l <- x - (current_sum - m_eff[, l])

      update_effect(l, y_l)
    }

    if (update_sigma2) {
      for (l in seq_len(L)) {
        Eb2_l <- vb[l, ] + mb[l, ]^2
        sigma2[l] <- sum(pi_mat[l, ] * Eb2_l)
        sigma2[l] <- max(sigma2[l], 1e-10)
      }
    }

    if (update_nu0) {
      nu0 <- update_nu0_EB()
    }

    elbo_trace[iter] <- compute_elbo()

    if (verbose) {
      cat(sprintf(
        "iter %03d | ELBO %.6f | nu0 %.6f\n",
        iter,
        elbo_trace[iter],
        nu0
      ))
    }

    if (iter > 1) {
      rel_change <- abs(elbo_trace[iter] - old_elbo) /
        (abs(old_elbo) + 1e-8)

      if (rel_change < tol) {
        elbo_trace <- elbo_trace[seq_len(iter)]
        break
      }
    }
  }

  PIP <- 1 - apply(1 - pi_mat, 2, prod)

  list(
    pi = pi_mat,
    PIP = PIP,
    nu0 = nu0,
    sigma2 = sigma2,
    N = N,
    mb = mb,
    vb = vb,
    Es = Es,
    Es2 = Es2,
    ElogS = ElogS,
    EInvS = EInvS,
    Eomega = Eomega,
    Elogomega = Elogomega,
    alpha_z = alpha_z,
    A_z = A_z,
    B_z = B_z,
    Delta_z = Delta_z,
    m_effects = m_eff,
    fitted = rowSums(m_eff),
    elbo = elbo_trace,
    Rbar_used = Rbar
  )
}


set.seed(1)

J <- 300
L <- 2
sample_size <- 10000
nu0_true <- 10000
sigma2_fit <- rep(1, L)

# 1. Choose a correlation matrix Rbar with eigenvalues bounded below.
rho_ar <- 0.35
idx <- seq_len(J)
Rbar <- rho_ar^abs(outer(idx, idx, "-"))
eig_min <- min(eigen(Rbar, symmetric = TRUE, only.values = TRUE)$values)
stopifnot(eig_min > 0.1)

# 2. Generate R ~ IW(nu0_true * Rbar, nu0_true + J + 1).
riwish <- function(scale, df) {
  W <- rWishart(1, df = df, Sigma = solve(scale))[, , 1]
  R <- solve(W)
  (R + t(R)) / 2
}

R_true <- riwish(
  scale = nu0_true * Rbar,
  df = nu0_true + J + 1
)

# 3. Choose one column of R and form x_bar.
gamma_true <- c(35, J / 2)
beta_true <- c(1.0, 0.9) * 1e-1
x_bar <- as.numeric(R_true[, gamma_true, drop = FALSE] %*% beta_true)

# 4. Generate x from N(x_bar, Rbar / N).
z <- rnorm(J)
x <- x_bar + as.numeric(t(chol(Rbar)) %*% z) / sqrt(sample_size)

# 5. Fit and check recovery of nu0 and selected coordinates.
fit <- iw_susie_veb(
  x = x,
  Rbar = Rbar,
  L = L,
  N = sample_size,
  sigma2 = sigma2_fit,
  nu0_init = 100,
  max_iter = 1000,
  inner_iter = 3,
  update_nu0 = TRUE,
  update_sigma2 = FALSE,
  tol = 1e-6,
  verbose = FALSE
)

top_by_effect <- apply(fit$pi, 1, which.max)
top_pip <- order(fit$PIP, decreasing = TRUE)[seq_len(L)]
coordinate_recovered <- all(gamma_true %in% top_pip)
nu0_relative_error <- abs(fit$nu0 - nu0_true) / nu0_true
recovery <- data.frame(
  effect = seq_len(L),
  gamma_true = gamma_true,
  beta_true = beta_true,
  top_gamma = top_by_effect,
  prob_true_gamma = fit$pi[cbind(seq_len(L), gamma_true)],
  prob_top_gamma = fit$pi[cbind(seq_len(L), top_by_effect)]
)

summarize_fit <- function(fit_obj, label) {
  top_by_effect_obj <- apply(fit_obj$pi, 1, which.max)
  top_pip_obj <- order(fit_obj$PIP, decreasing = TRUE)[seq_len(L)]
  recovery_obj <- data.frame(
    effect = seq_len(L),
    gamma_true = gamma_true,
    beta_true = beta_true,
    top_gamma = top_by_effect_obj,
    prob_true_gamma = fit_obj$pi[cbind(seq_len(L), gamma_true)],
    prob_top_gamma = fit_obj$pi[cbind(seq_len(L), top_by_effect_obj)]
  )

  cat(sprintf("\n%s recovery summary\n", label))
  cat(sprintf("nu0: %.3f | final ELBO: %.6f\n", fit_obj$nu0, tail(fit_obj$elbo, 1)))
  cat("top PIP coordinates:\n")
  print(top_pip_obj)
  cat(sprintf("all true coordinates in top %d PIPs: %s\n",
              L, all(gamma_true %in% top_pip_obj)))
  print(recovery_obj)

  invisible(list(
    top_by_effect = top_by_effect_obj,
    top_pip = top_pip_obj,
    recovery = recovery_obj
  ))
}

cat("\nSimulation recovery summary\n")
cat(sprintf("true nu0: %.3f | estimated nu0: %.3f\n", nu0_true, fit$nu0))
cat(sprintf("min eigenvalue of Rbar: %.3f\n", eig_min))
cat(sprintf("relative nu0 error: %.3f\n", nu0_relative_error))
cat("true gamma:\n")
print(gamma_true)
cat("top PIP coordinates:\n")
print(top_pip)
cat(sprintf("all true coordinates in top %d PIPs: %s\n", L, coordinate_recovered))
print(recovery)

plot(fit$elbo, type = "b", xlab = "Iteration", ylab = "ELBO")
plot(fit$PIP, type = "h", xlab = "Index j", ylab = "PIP")
abline(v = gamma_true, col = "red", lty = 2)


plot(seq_len(J), x, type = "l", xlab = "Index j", ylab = "x")

first_true_part <- as.numeric(R_true[, gamma_true[1]] * beta_true[1])
first_estimated_part <- fit$m_effects[, 1]

second_true_part <- as.numeric(R_true[, gamma_true[2]] * beta_true[2])
second_estimated_part <- fit$m_effects[, 2]

part_ylim <- range(
  first_true_part,
  first_estimated_part,
  second_true_part,
  second_estimated_part
)

plot(
  seq_len(J),
  first_true_part,
  type = "l",
  col = "black",
  lwd = 2,
  ylim = part_ylim,
  xlab = "Index j",
  ylab = "Contribution",
  main = "First Effect Contribution"
)
lines(seq_len(J), first_estimated_part, col = "blue", lwd = 2, lty = 2)
abline(v = gamma_true[1], col = "red", lty = 3)
legend(
  "topright",
  legend = c("true_part", "estimated_part"),
  col = c("black", "blue"),
  lty = c(1, 2),
  lwd = 2,
  bty = "n"
)

plot(
  seq_len(J),
  second_true_part,
  type = "l",
  col = "black",
  lwd = 2,
  ylim = part_ylim,
  xlab = "Index j",
  ylab = "Contribution",
  main = "Second Effect Contribution"
)
lines(seq_len(J), second_estimated_part, col = "blue", lwd = 2, lty = 2)
abline(v = gamma_true[2], col = "red", lty = 3)
legend(
  "topright",
  legend = c("true_part", "estimated_part"),
  col = c("black", "blue"),
  lty = c(1, 2),
  lwd = 2,
  bty = "n"
)

# Fit the same data with nu0 fixed at the truth. This checks whether the
# coordinate and contribution recovery improve when the R-column uncertainty is
# constrained to the data-generating level.
fit_fixed_nu0 <- iw_susie_veb(
  x = x,
  Rbar = Rbar,
  L = L,
  N = sample_size,
  sigma2 = sigma2_fit,
  nu0_init = nu0_true,
  max_iter = 1000,
  inner_iter = 3,
  update_nu0 = FALSE,
  update_sigma2 = FALSE,
  tol = 1e-6,
  verbose = FALSE
)

fixed_summary <- summarize_fit(fit_fixed_nu0, "Fixed true nu0")

cat("\nELBO comparison\n")
cat(sprintf("EB nu0 fit final ELBO: %.6f\n", tail(fit$elbo, 1)))
cat(sprintf("Fixed true nu0 fit final ELBO: %.6f\n", tail(fit_fixed_nu0$elbo, 1)))
cat(sprintf("ELBO fixed - EB: %.6f\n", tail(fit_fixed_nu0$elbo, 1) - tail(fit$elbo, 1)))

plot(
  fit$elbo,
  type = "l",
  col = "blue",
  lwd = 2,
  xlab = "Iteration",
  ylab = "ELBO",
  main = "ELBO Comparison"
)
lines(fit_fixed_nu0$elbo, col = "black", lwd = 2, lty = 2)
legend(
  "bottomright",
  legend = c("EB nu0", "fixed true nu0"),
  col = c("blue", "black"),
  lty = c(1, 2),
  lwd = 2,
  bty = "n"
)

first_estimated_part_fixed <- fit_fixed_nu0$m_effects[, 1]
second_estimated_part_fixed <- fit_fixed_nu0$m_effects[, 2]

fixed_part_ylim <- range(
  first_true_part,
  first_estimated_part_fixed,
  second_true_part,
  second_estimated_part_fixed
)

plot(
  seq_len(J),
  first_true_part,
  type = "l",
  col = "black",
  lwd = 2,
  ylim = fixed_part_ylim,
  xlab = "Index j",
  ylab = "Contribution",
  main = "First Effect Contribution, Fixed True nu0"
)
lines(seq_len(J), first_estimated_part_fixed, col = "blue", lwd = 2, lty = 2)
abline(v = gamma_true[1], col = "red", lty = 3)
legend(
  "topright",
  legend = c("true_part", "estimated_part_fixed"),
  col = c("black", "blue"),
  lty = c(1, 2),
  lwd = 2,
  bty = "n"
)

plot(
  seq_len(J),
  second_true_part,
  type = "l",
  col = "black",
  lwd = 2,
  ylim = fixed_part_ylim,
  xlab = "Index j",
  ylab = "Contribution",
  main = "Second Effect Contribution, Fixed True nu0"
)

lines(seq_len(J), second_estimated_part_fixed, col = "blue", lwd = 2, lty = 2)
abline(v = gamma_true[2], col = "red", lty = 3)
legend(
  "topright",
  legend = c("true_part", "estimated_part_fixed"),
  col = c("black", "blue"),
  lty = c(1, 2),
  lwd = 2,
  bty = "n"
)

# Warm-start a fixed-nu0 fit from the EB-recovered coordinates. The starting
# effect sizes are GLS coefficients using Rbar columns at the top PIP indices.
init_gamma_warm <- top_pip
X_warm <- Rbar[, init_gamma_warm, drop = FALSE]
KX_warm <- chol_solve(chol(Rbar), X_warm)
Kx_warm <- chol_solve(chol(Rbar), x)
init_beta_warm <- as.numeric(solve(t(X_warm) %*% KX_warm, t(X_warm) %*% Kx_warm))

fit_fixed_nu0_warm <- iw_susie_veb(
  x = x,
  Rbar = Rbar,
  L = L,
  N = sample_size,
  sigma2 = sigma2_fit,
  nu0_init = nu0_true,
  max_iter = 1000,
  inner_iter = 3,
  update_nu0 = FALSE,
  update_sigma2 = FALSE,
  tol = 1e-6,
  init_gamma = init_gamma_warm,
  init_beta = init_beta_warm,
  verbose = FALSE
)

fixed_warm_summary <- summarize_fit(
  fit_fixed_nu0_warm,
  "Fixed true nu0, EB-coordinate warm start"
)

cat("\nWarm-start initialization\n")
cat("init gamma:\n")
print(init_gamma_warm)
cat("init beta:\n")
print(init_beta_warm)

cat("\nELBO comparison including warm start\n")
cat(sprintf("EB nu0 fit final ELBO: %.6f\n", tail(fit$elbo, 1)))
cat(sprintf("Fixed true nu0 cold final ELBO: %.6f\n", tail(fit_fixed_nu0$elbo, 1)))
cat(sprintf("Fixed true nu0 warm final ELBO: %.6f\n", tail(fit_fixed_nu0_warm$elbo, 1)))
cat(sprintf("ELBO warm fixed - EB: %.6f\n",
            tail(fit_fixed_nu0_warm$elbo, 1) - tail(fit$elbo, 1)))

plot(
  fit$elbo,
  type = "l",
  col = "blue",
  lwd = 2,
  xlab = "Iteration",
  ylab = "ELBO",
  main = "ELBO Comparison Including Warm Start"
)
lines(fit_fixed_nu0$elbo, col = "black", lwd = 2, lty = 2)
lines(fit_fixed_nu0_warm$elbo, col = "darkgreen", lwd = 2, lty = 3)
legend(
  "bottomright",
  legend = c("EB nu0", "fixed true nu0 cold", "fixed true nu0 warm"),
  col = c("blue", "black", "darkgreen"),
  lty = c(1, 2, 3),
  lwd = 2,
  bty = "n"
)

first_match_warm <- match(gamma_true[1], init_gamma_warm)
second_match_warm <- match(gamma_true[2], init_gamma_warm)
first_estimated_part_fixed_warm <- fit_fixed_nu0_warm$m_effects[, first_match_warm]
second_estimated_part_fixed_warm <- fit_fixed_nu0_warm$m_effects[, second_match_warm]

warm_part_ylim <- range(
  first_true_part,
  first_estimated_part_fixed_warm,
  second_true_part,
  second_estimated_part_fixed_warm
)

plot(
  seq_len(J),
  first_true_part,
  type = "l",
  col = "black",
  lwd = 2,
  ylim = warm_part_ylim,
  xlab = "Index j",
  ylab = "Contribution",
  main = "First Effect Contribution, Fixed nu0 Warm Start"
)
lines(seq_len(J), first_estimated_part_fixed_warm, col = "blue", lwd = 2, lty = 2)
abline(v = gamma_true[1], col = "red", lty = 3)
legend(
  "topright",
  legend = c("true_part", "estimated_part_fixed_warm"),
  col = c("black", "blue"),
  lty = c(1, 2),
  lwd = 2,
  bty = "n"
)

plot(
  seq_len(J),
  second_true_part,
  type = "l",
  col = "black",
  lwd = 2,
  ylim = warm_part_ylim,
  xlab = "Index j",
  ylab = "Contribution",
  main = "Second Effect Contribution, Fixed nu0 Warm Start"
)
lines(seq_len(J), second_estimated_part_fixed_warm, col = "blue", lwd = 2, lty = 2)
abline(v = gamma_true[2], col = "red", lty = 3)
legend(
  "topright",
  legend = c("true_part", "estimated_part_fixed_warm"),
  col = c("black", "blue"),
  lty = c(1, 2),
  lwd = 2,
  bty = "n"
)
