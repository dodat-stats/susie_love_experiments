library(data.table)
library(susieR)

data_dir <- "/Users/td250/Dat/ukb-height"
N0 <- 100  # out-of-sample size for LD reference
set.seed(15)

# --- Load phenotype + covariates ---
pheno <- fread(file.path(data_dir, "height.csv.gz"))

# --- Load genotype dosage matrix (PLINK --recode A output) ---
geno <- fread(file.path(data_dir, "height.chr10.104285594.raw.gz"))

# Merge on individual ID
pheno[, IID := id]
merged <- merge(geno, pheno, by = "IID")

# --- Extract dosage matrix ---
dosage_cols <- names(geno)[-(1:6)]  # skip FID IID PAT MAT SEX PHENOTYPE
X <- as.matrix(merged[, ..dosage_cols])

# --- Residualise height on covariates ---
covar_cols <- c("sex", "age", "assessment_centre", "genotype_measurement_batch",
                "missingness", paste0("pc_genetic", 1:40))

covar_formula <- as.formula(paste("height ~", paste(covar_cols, collapse = " + ")))
fit <- lm(covar_formula, data = merged, na.action = na.exclude)
y <- residuals(fit)

# Also residualise each column of X on the same covariates
covar_mat <- model.matrix(~ ., data = merged[, ..covar_cols])[, -1]

# Drop rows with NA residuals
keep <- !is.na(y)
y <- y[keep]
X <- X[keep, ]
covar_mat <- covar_mat[keep, ]

# Residualise X on covariates via projection
Q <- qr.Q(qr(cbind(1, covar_mat)))
X <- X - Q %*% crossprod(Q, X)

# Drop variants with zero variance
col_var <- colSums(X^2)
X <- X[, col_var > 0]

# Mean-impute any remaining missing genotypes
for (j in seq_len(ncol(X))) {
  na_idx <- is.na(X[, j])
  if (any(na_idx)) X[na_idx, j] <- mean(X[, j], na.rm = TRUE)
}

rm(geno, merged, covar_mat, Q)
gc()

# ============================================================
# Split into in-sample (GWAS) and out-of-sample (LD reference)
# ============================================================
n_total <- length(y)
idx_oos <- sample(n_total, N0)
idx_in  <- setdiff(seq_len(n_total), idx_oos)

X_in  <- X[idx_in, ]
y_in  <- y[idx_in]
X_oos <- X[idx_oos, ]

rm(X)
gc()

n_in <- length(y_in)

# Restrict to variants 900-2000 for tractable fine-mapping
var_idx <- 900:min(2000, ncol(X_in))
X_in  <- X_in[, var_idx]
X_oos <- X_oos[, var_idx]

cat(sprintf("In-sample: %d individuals, Out-of-sample LD ref: %d individuals, %d variants\n",
            n_in, N0, ncol(X_in)))

# --- Standardize in-sample genotypes and compute sufficient statistics ---
X_in_scaled <- scale(X_in)
XtX_in <- crossprod(X_in_scaled)
Xty    <- as.numeric(crossprod(X_in_scaled, y_in))
yty    <- sum(y_in^2)
bhat   <- Xty / (n_in - 1)

# --- In-sample LD (correlation) matrix ---
cat("Computing in-sample LD matrix...\n")
R_in <- XtX_in / (n_in - 1)
gc()

rm(X_in, X_in_scaled)

# --- Out-of-sample LD matrix ---
cat("Computing out-of-sample LD matrix...\n")
X_oos_scaled <- scale(X_oos)
R_oos <- crossprod(X_oos_scaled) / (nrow(X_oos_scaled) - 1)
rm(X_oos, X_oos_scaled)
gc()

# ============================================================
# Fit 1: SuSiE with in-sample LD
# ============================================================
cat("\n=== Fitting SuSiE with in-sample LD ===\n")
fit_in <- susie_suff_stat(XtX = XtX_in, Xty = Xty, yty = yty, n = n_in,
                           L = 10, standardize = FALSE,
                           estimate_residual_variance = TRUE,
                           estimate_prior_variance = TRUE,
                           verbose = TRUE)

cs_in <- susie_get_cs(fit_in, Xcorr = R_in, coverage = 0.95)
cs_in_list <- if (!is.null(cs_in$cs)) cs_in$cs else cs_in
pip_in <- susie_get_pip(fit_in)

cat(sprintf("In-sample LD: %d credible sets\n", length(cs_in_list)))

# ============================================================
# Fit 2: SuSiE with out-of-sample LD
# ============================================================
cat("\n=== Fitting SuSiE with out-of-sample LD ===\n")
fit_oos <- susie_suff_stat(XtX = R_oos * (n_in - 1), Xty = Xty, yty = yty, n = n_in,
                            L = 10, standardize = FALSE,
                            estimate_residual_variance = TRUE,
                            estimate_prior_variance = TRUE,
                            verbose = TRUE)

cs_oos <- susie_get_cs(fit_oos, Xcorr = R_oos, coverage = 0.95)
cs_oos_list <- if (!is.null(cs_oos$cs)) cs_oos$cs else cs_oos
pip_oos <- susie_get_pip(fit_oos)

cat(sprintf("Out-of-sample LD: %d credible sets\n", length(cs_oos_list)))

# ============================================================
# Compare results
# ============================================================
pvar <- fread(file.path(data_dir, "height.chr10.104285594.pvar"))
names(pvar) <- c("CHROM", "POS", "ID", "REF", "ALT")

variant_names <- colnames(R_in)
rsids <- sub("_[^_]+$", "", variant_names)

comparison <- data.table(
  rsid    = rsids,
  bhat    = bhat,
  pip_in  = pip_in,
  pip_oos = pip_oos,
  cs_in   = NA_integer_,
  cs_oos  = NA_integer_
)

for (i in seq_along(cs_in_list))  comparison[cs_in_list[[i]],  cs_in  := i]
for (i in seq_along(cs_oos_list)) comparison[cs_oos_list[[i]], cs_oos := i]

comparison <- merge(comparison, pvar[, .(ID, CHROM, POS)],
                    by.x = "rsid", by.y = "ID", all.x = TRUE)
setorder(comparison, -pip_in)

# ============================================================
# Plot comparison
# ============================================================
par(mfrow = c(2, 1), mar = c(4, 4, 3, 1))
susie_plot(fit_in,  y = "PIP", main = "In-sample LD")
susie_plot(fit_oos, y = "PIP", main = "Out-of-sample LD")

cat(sprintf("\nPIP correlation (in vs out-of-sample LD): %.4f\n",
            cor(pip_in, pip_oos)))

cat("\n=== Credible set sizes ===\n")
cat("In-sample LD:\n")
for (i in seq_along(cs_in_list))  cat(sprintf("  CS%d: %d variants\n", i, length(cs_in_list[[i]])))
cat("Out-of-sample LD:\n")
for (i in seq_along(cs_oos_list)) cat(sprintf("  CS%d: %d variants\n", i, length(cs_oos_list[[i]])))





par(mfrow = c(2, 1), mar = c(4, 4, 3, 1))

### small shrinkage
J = dim(R_oos)[1]
lambda <- 1 / (N0 + 1)
R_tilde = (1- lambda) * R_oos + lambda * diag(R_oos)

fit_tilde <- susie_suff_stat(XtX = R_tilde * (n_in - 1), Xty = Xty, yty = yty, n = n_in,
                             L = 10, standardize = FALSE,
                             estimate_residual_variance = TRUE,
                             estimate_prior_variance = TRUE,
                             verbose = TRUE)

susie_plot(fit_tilde,  y = "PIP", main = "Regularized LD")


###
library(susielove)

J = dim(R_oos)[1]
lambda <- .1 / (N0 + 1)
diag_R <- diag(R_oos)
R_out  <- R_oos
k      <- min(N0, J)

# Standardize v to var(y)=1 scale (susie_love assumes unit residual variance)
s_y <- sqrt(yty / (n_in - 1))
bhat_love <- bhat / s_y

love_res <- susie_love_lambda(bhat_love, R_out, diag_R,
                              lambda,
                              n_in, L = 5,
                              k,
                              nu_delta_init = 1000,
                              step_susie = 10,
                              step_R = 30,
                              step_nu = 10,
                              max_iter = 30,
                              max_iter_init = 1000,
                              tol = 1e-4)

b_res <- love_res$b_res
R_q   <- love_res$R_q

susie_plot_alpha(
  alpha = b_res$alpha,
  R = cov2cor(R_q),
  main = paste0("SuSiE-LovE with lambda = ", lambda)
)



## WS shrinkage
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

diag_R <- diag(R_oos)
R_out  <- compute_wen_stephens(R_0 = R_oos, N_ref = N0)
k      <- min(N0, J)

love_res_WS <- susie_love_lambda(bhat_love, R_out, diag_R,
                                 lambda = 0,
                                 N = n_in, L = 5,
                                 k,
                                 nu_delta_init = 1000,
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
  main = paste0("SuSiE-LovE with WS shrinakge")
)






### test oos with the "correct" L
fit_oos <- susie_suff_stat(XtX = R_oos * (n_in - 1), Xty = Xty, yty = yty, n = n_in,
                           L = 2, standardize = FALSE,
                           estimate_residual_variance = TRUE,
                           estimate_prior_variance = TRUE,
                           verbose = TRUE)
susie_plot(fit_oos,  y = "PIP", main = "Oos with correct L")
