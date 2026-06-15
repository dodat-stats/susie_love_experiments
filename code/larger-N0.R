

library(susieR)
library(susielove)

seed  <- 1
T_split <- 50
chrom <- "chr22"
start <- 20e6
end   <- 21e6
N     <- 25000L
N0    <- 1000L
J     <- 1000L
h2    <- 0.005

res <- sim_2pop(chrom, start, end,
                T_split = T_split,
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

plot(v, pch = 16)
points(causal_SNPs_loc, v[causal_SNPs_loc], pch = 16, col = "red")



#####
fit <- susie_suff_stat(
  XtX = N * R,
  Xty = N * v,
  n = N,
  yty = N
)
susie_plot(fit, y = "PIP", b = beta_true, main = "In-sample cov mat")


#####
fit_R0 <- susie_suff_stat(
  XtX = N * R0,
  Xty = N * v,
  n = N,
  yty = N,
  L = 5
)
susie_plot(fit_R0, y = "PIP", b = beta_true, main = "Out-of-sample cov mat")

## susie-love-lambda
nu_delta_init = 1000
lambda <- 1 / (N0 + 1)
diag_R <- diag(R)
R_out  <- R0
k      <- min(N0, J)

love_res <- susie_love_lambda(v, R_out, diag_R,
                              lambda,
                              N, L = 5,
                              k,
                              nu_delta_init = nu_delta_init,
                              step_susie = 20,
                              step_R = 30,
                              step_nu = 10,
                              max_iter = 30,
                              tol = 1e-4)

b_res <- love_res$b_res
R_q   <- love_res$R_q

susie_plot_alpha(
  alpha = b_res$alpha,
  R = cov2cor(R_q),
  b = beta_true,
  main = paste0("SuSiE-LovE with lambda = ", lambda)
)










