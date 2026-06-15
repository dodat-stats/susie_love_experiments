library(susielove)
library(susieR)

compute_wen_stephens <- function(R_0, N_ref) {
  # Number of haplotypes (diploid population)
  d <- 2 * N_ref

  # Calculate Watterson's estimator constant (Harmonic sum)
  harmonic_sum <- sum(1 / (1:(d - 1)))
  theta        <- 1 / harmonic_sum

  # Calculate the exact mutation shrinkage parameter lambda
  lambda <- theta / (theta + 2 * N_ref)

  # Apply local uniform shrinkage: Psi = (1-lambda)^2 * R_0 + (1 - (1-lambda)^2) * I
  p <- ncol(R_0)
  R_ws <- (1 - lambda)^2 * R_0 + (1 - (1 - lambda)^2) * diag(p)

  return(R_ws)
}


compute_lw_linear_from_R0 <- function(R_0, N_ref) {
  p <- ncol(R_0)

  # 1. Sum of the squared off-diagonal elements of R_0
  # (This is the Frobenius norm of the off-diagonals)
  sum_sq_off_diag <- sum(R_0^2) - p

  # 2. Estimate the sampling variance (noise) parameter under normality
  # In correlation space, the sum of variances of r_ij is approximated by:
  beta_sq <- (1 / N_ref) * (p + sum_sq_off_diag)

  # 3. Calculate the distance to the identity target matrix
  # Since diag(R_0) is all 1s and diag(I) is all 1s, the distance is just
  # the sum of the squared off-diagonal elements.
  delta_sq <- sum_sq_off_diag

  # 4. Compute the optimal shrinkage intensity alpha
  # Guard against division by zero if R_0 is already the identity matrix
  if (delta_sq == 0) {
    alpha <- 1
  } else {
    alpha <- min(1, beta_sq / delta_sq)
  }

  # 5. Construct the shrunk matrix: alpha * I + (1 - alpha) * R_0
  R_lw_linear <- (1 - alpha) * R_0 + alpha * diag(p)

  # Print the learned alpha for transparency
  cat("Ledoit-Wolf Linear Alpha estimated from R_0:", alpha, "\n")

  return(R_lw_linear)
}

sim_asn_eur <- function(chrom, start, end, n_gwas, n_ref, J, h2, seed) {

  stdpopsim <- reticulate::import("stdpopsim")
  msprime   <- reticulate::import("msprime")

  set.seed(seed)

  # 1. Get Realistic Genetic Map from stdpopsim
  species <- stdpopsim$get_species("HomSap")
  contig  <- species$get_contig(chrom, left = as.integer(start), right = as.integer(end))

  # 2. Get OutOfAfrica_3G09 Demography (YRI, CEU, CHB)
  demo_model <- species$get_demographic_model("OutOfAfrica_3G09")
  demography <- demo_model$model

  # 3. Run Simulation
  ts <- msprime$sim_ancestry(
    samples = reticulate::dict(CHB = as.integer(n_gwas), CEU = as.integer(n_ref)),
    demography = demography,
    recombination_rate = contig$recombination_map,
    random_seed = as.integer(seed)
  )

  ts <- msprime$sim_mutations(ts, rate = 1.29e-8, random_seed = as.integer(seed))

  # 4. Extract Genotypes
  chb_id <- NULL
  ceu_id <- NULL
  for (p in reticulate::iterate(ts$populations())) {
    if (!is.null(p$metadata) && p$metadata$name == "CHB") chb_id <- p$id
    if (!is.null(p$metadata) && p$metadata$name == "CEU") ceu_id <- p$id
  }

  samples_gwas <- ts$samples(population = as.integer(chb_id))
  samples_ref  <- ts$samples(population = as.integer(ceu_id))

  G_all <- ts$genotype_matrix()

  to_diploid <- function(hap_idx) {
    idx <- as.integer(hap_idx) + 1
    g_hap <- G_all[, idx, drop = FALSE]
    g_dip <- g_hap[, seq(1, ncol(g_hap), by = 2)] + g_hap[, seq(2, ncol(g_hap), by = 2)]
    return(t(g_dip))
  }

  G_gwas <- to_diploid(samples_gwas)
  G_ref  <- to_diploid(samples_ref)

  # 5. Filter for Common Variants (MAF > 5%)
  maf_gwas <- colMeans(G_gwas) / 2
  maf_ref  <- colMeans(G_ref) / 2

  common_mask <- (maf_gwas > 0.05) & (maf_gwas < 0.95) & (maf_ref > 0.05) & (maf_ref < 0.95)
  n_common <- sum(common_mask)
  if (n_common < J) stop(paste0("Only ", n_common, " common variants found but J = ", J, ". Try a larger region."))
  G_gwas <- G_gwas[, common_mask, drop = FALSE]
  G_ref  <- G_ref[, common_mask, drop = FALSE]

  J_all <- min(ncol(G_gwas), J)
  max_start <- J_all - J
  start_idx <- if (max_start > 0) sample(0:max_start, 1) else 0
  cols_to_keep <- (start_idx + 1):(start_idx + J)
  G  <- G_gwas[, cols_to_keep, drop = FALSE]
  G0 <- G_ref[, cols_to_keep, drop = FALSE]

  # 6. Generate summary stats
  X <- sweep(G, 2, colMeans(G), "-")
  N <- nrow(X)

  G0_centered <- sweep(G0, 2, colMeans(G0), "-")

  R0 <- stats::cov(G0_centered)
  R0[is.na(R0)] <- 0

  R <- stats::cov(X)
  R[is.na(R)] <- 0

  num_causal <- 2
  first_causal_SNP <- sample(1:J, 1)
  second_causal_SNP <- which.min(R[first_causal_SNP, ]^2)
  causal_SNPs_loc <- c(first_causal_SNP, second_causal_SNP)

  beta_true <- rep(0, J)
  beta_true[causal_SNPs_loc] <- stats::runif(num_causal) + 0.01

  var_expl <- sum((X %*% beta_true)^2 / N)
  scale <- h2 / var_expl
  beta_true <- beta_true * sqrt(scale)

  y <- (X %*% beta_true) + sqrt(1 - h2) * stats::rnorm(N)
  y_std <- sqrt(sum((y - mean(y))^2) / N)
  y <- (y - mean(y)) / y_std

  v <- (t(X) %*% y) / N

  return(list(
    v = as.vector(v),
    R = R,
    R0 = R0,
    beta_true = beta_true,
    causal_SNPs_loc = causal_SNPs_loc
  ))
}


seed = 2
chrom <- "chr22"
start <- 20e6
end   <- 21e6
N     <- 25000L
N0 = 500L
J     <- 500L
h2    <- 0.01

res <- sim_asn_eur(chrom, start, end,
                   n_gwas = N,
                   n_ref = N0,
                   J = J,
                   h2 = h2,
                   seed = seed)

v               <- res$v
R               <- res$R
R0              <- res$R0
beta_true       <- res$beta_true
causal_SNPs_loc <- res$causal_SNPs_loc

## standardized
d = 1 / sqrt(diag(R))
v = d * v
R = cov2cor(R)
R0 = cov2cor(R0)

plot(v, pch = 16)
points(causal_SNPs_loc, v[causal_SNPs_loc], pch = 16, col = "red")

##########
par(mfrow = c(1, 2))

image(R0[1:100, 1:100], main = "Out-of-sample (EUR)", col = heat.colors(15))
image(R[1:100, 1:100], main = "In-sample corr (CHB)", col = heat.colors(15))
par(mfrow = c(1, 1))

#########
fit_R <- susie_suff_stat(
  XtX = N * R,
  Xty = N * v,
  n = N,
  yty = N,
  L = 5
)
susie_plot(fit_R, y = "PIP", b = beta_true, main = "In-sample cov mat")

#########
fit_R_silly <- susie_suff_stat(
  XtX = N * R0,
  Xty = N * v,
  n = N,
  yty = N,
  L = 5,
  residual_variance = 1
)
susie_plot(fit_R_silly, y = "PIP", b = beta_true, main = "Silly corr mat")



R_tilde = (1 - 1 / N0) * R0 + 1 / N0 * diag(J)
fit_R_tilde <- susie_suff_stat(
  XtX = N * R_tilde,
  Xty = N * v,
  n = N,
  yty = N,
  L = 5
)
susie_plot(fit_R_tilde, y = "PIP", b = beta_true, main = "Silly corr mat")

print(summary(fit_R_tilde)$cs)


## ---- susielove
R_out  <- compute_wen_stephens(R_0 = R0, N_ref = N0)
# R_out = diag(J)
# R_out = R_tilde
k      <- min(N0, J)
diag_R = diag(R)

love_res_WS <- susie_love_lambda(v, R_out, diag_R,
                                 lambda = 0,
                                 N, L = 5,
                                 k,
                                 nu_delta_init = 500,
                                 step_susie = 20,
                                 step_R = 30,
                                 step_nu = 10,
                                 max_iter = 30,
                                 tol = 1e-4)

b_res <- love_res_WS$b_res
R_q   <- love_res_WS$R_q

susie_plot_alpha(
  alpha = b_res$alpha,
  R = cov2cor(R_q),
  b = beta_true,
  main = paste0("SuSiE-LovE with WS shrinkage")
)


loss_val = susielove:::compute_loss(love_res_WS$b_res,
                                    love_res_WS$Z_final,
                                    love_res_WS$nu_q,
                                    love_res_WS$nu_delta,
                                    v, N, R_out)
loss_val

loss_list1 = susielove:::loss_decompose(love_res_WS$b_res,
                                    love_res_WS$Z_final,
                                    love_res_WS$nu_q,
                                    love_res_WS$nu_delta,
                                    v, N, R_out)
loss_list1$loss_total
cat("loss fit", loss_list1$loss_fit)
cat("loss b KL", loss_list1$loss_b_kl)

cat("loss cov", loss_list1$loss_cov)
cat("loss logdet term", loss_list1$loss_logdet)
cat("loss trace term", loss_list1$loss_trace)
cat("loss nu term", loss_list1$loss_nu)




###  ---- susielove with higher lambda
R_out  <- compute_wen_stephens(R_0 = R0, N_ref = N0)
R_out = (1 - 0.1) * R_out + 0.1 * diag(J)
k      <- min(N0, J)
diag_R = diag(R)

love_res2 <- susie_love_lambda(v, R_out, diag_R,
                                 lambda = 0,
                                 N, L = 5,
                                 k,
                                 nu_delta_init = 500,
                                 step_susie = 20,
                                 step_R = 30,
                                 step_nu = 10,
                                 max_iter = 30,
                                 tol = 1e-4)

b_res <- love_res2$b_res
R_q   <- love_res2$R_q

susie_plot_alpha(
  alpha = b_res$alpha,
  R = cov2cor(R_q),
  b = beta_true,
  main = paste0("SuSiE-LovE with WS shrinkage")
)

loss_list2 = susielove:::loss_decompose(love_res2$b_res,
                                       love_res2$Z_final,
                                       love_res2$nu_q,
                                       love_res2$nu_delta,
                                       v, N, R_out)
loss_list2$loss_total
cat("loss fit", loss_list2$loss_fit)
cat("loss b KL", loss_list2$loss_b_kl)

cat("loss cov", loss_list2$loss_cov)
cat("loss logdet term", loss_list2$loss_logdet)
cat("loss trace term", loss_list2$loss_trace)
cat("loss nu term", loss_list2$loss_nu)
loss_list2$loss_trace + loss_list2$loss_nu

cat("loss cov", loss_list1$loss_cov)
cat("loss cov", loss_list2$loss_cov)


cat("loss fit term Rq", loss_list1$loss_fit_Rq)
cat("loss fit term Rq", loss_list2$loss_fit_Rq)
cat("loss fit term null", loss_list1$loss_fit_null)
cat("loss fit term null", loss_list2$loss_fit_null)

