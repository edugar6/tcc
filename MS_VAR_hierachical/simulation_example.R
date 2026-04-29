##---------------
rm(list = ls())

current_working_dir = dirname(rstudioapi::getSourceEditorContext()$path)
setwd(current_working_dir)

# =============================================================================
# simulation_example.R
# Simulate data from the MS-VAR model and verify estimation via utils.R
# Asset-specific idiosyncratic states + transition matrix deviation analysis
# =============================================================================

source("utils.R")

set.seed(2024)

# -----------------------------------------------------------------------------
# SIMULATION PARAMETERS
# -----------------------------------------------------------------------------

TT <- 500
n  <- 3
K  <- 3    # bear=1, stable=2, bull=3
M  <- 2    # idiosyncratic states per asset

# True transition matrix (structural zeros: bear<->bull)
P_true <- matrix(c(
  0.85, 0.15, 0.00,
  0.05, 0.90, 0.05,
  0.00, 0.15, 0.85
), nrow = K, byrow = TRUE)

# True intercepts per global regime + idiosyncratic states
# mu[[k]]   : n-vector for global regime k
# mu[[K+m]] : n-vector for idiosyncratic state m (used row-by-row per asset)
mu_true <- list(
  c(-0.30, -0.25, -0.35),   # bear
  c( 0.00,  0.01, -0.01),   # stable
  c( 0.30,  0.25,  0.35),   # bull
  c(-0.10, -0.08, -0.12),   # idiosyncratic state 1 (asset-specific rows)
  c( 0.10,  0.08,  0.12)    # idiosyncratic state 2
)

# True AR matrices
A_bear <- matrix(c(
  0.10,  0.00,  0.00,
  0.20,  0.05,  0.00,   # asset 2 influenced by asset 1 in bear
  0.15,  0.00,  0.05    # asset 3 influenced by asset 1 in bear
), nrow = n, byrow = TRUE)

A_stable <- diag(0.05, n)

A_bull <- matrix(c(
  0.05,  0.00,  0.15,   # asset 1 influenced by asset 3 in bull
  0.00,  0.05,  0.20,   # asset 2 influenced by asset 3 in bull
  0.00,  0.00,  0.10
), nrow = n, byrow = TRUE)

# Idiosyncratic AR matrices (weak, asset-specific dynamics)
A_idio1 <- diag(0.03, n)
A_idio2 <- diag(0.03, n)

A_true <- list(A_bear, A_stable, A_bull, A_idio1, A_idio2)

# True participation probabilities rho[i, k]
rho_true <- matrix(c(
  0.85, 0.90, 0.85,
  0.80, 0.90, 0.82,
  0.88, 0.92, 0.87
), nrow = n, byrow = TRUE)

# True emission params: list[[i]][[s]], i=1..n, s=1..(K+M)
# Global states: distinct params per asset
# Idiosyncratic states: asset-specific local dynamics

sigma_glob <- list(
  bear   = c(2.0, 1.8, 2.2),
  stable = c(0.5, 0.4, 0.6),
  bull   = c(1.0, 0.9, 1.1)
)
alpha_glob <- list(
  bear   = c(-0.6, -0.5, -0.7),
  stable = c( 0.0,  0.0,  0.0),
  bull   = c( 0.4,  0.3,  0.5)
)
nu_glob <- list(
  bear   = c(5, 6, 5),
  stable = c(10, 10, 10),
  bull   = c(7, 8, 6)
)

# Idiosyncratic: each asset has own local params for state m=1,2
sigma_idio <- matrix(c(
  1.5, 1.3, 1.6,   # state K+1 per asset
  1.4, 1.2, 1.5    # state K+2 per asset
), nrow = M, byrow = TRUE)

alpha_idio <- matrix(0, nrow = M, ncol = n)
nu_idio    <- matrix(6,  nrow = M, ncol = n)

# Build params_emit[[i]][[s]]
params_true <- lapply(seq_len(n), function(i) {
  states <- vector("list", K + M)
  # Global states
  states[[1]] <- list(sigma = sigma_glob$bear[i],
                      alpha = alpha_glob$bear[i],
                      nu    = nu_glob$bear[i])
  states[[2]] <- list(sigma = sigma_glob$stable[i],
                      alpha = alpha_glob$stable[i],
                      nu    = nu_glob$stable[i])
  states[[3]] <- list(sigma = sigma_glob$bull[i],
                      alpha = alpha_glob$bull[i],
                      nu    = nu_glob$bull[i])
  # Idiosyncratic states
  for (m in seq_len(M)) {
    states[[K + m]] <- list(sigma = sigma_idio[m, i],
                            alpha = alpha_idio[m, i],
                            nu    = nu_idio[m, i])
  }
  states
})

# -----------------------------------------------------------------------------
# SIMULATE REGIME SEQUENCE z_t
# -----------------------------------------------------------------------------

z    <- integer(TT)
z[1] <- 2

for (t in 2:TT) {
  z[t] <- sample(K, 1, prob = P_true[z[t - 1], ])
}

cat("True regime frequencies:\n")
print(round(table(z) / TT, 3))

# -----------------------------------------------------------------------------
# SIMULATE RETURNS y_t
# -----------------------------------------------------------------------------

Y <- matrix(0, TT, n)

for (t in 2:TT) {
  k   <- z[t]
  yt1 <- Y[t - 1, ]
  
  for (i in seq_len(n)) {
    s_it <- rbinom(1, 1, rho_true[i, k])
    
    if (s_it == 1) {
      # Global regime: full AR row
      mu_val  <- mu_true[[k]][i]
      ar_val  <- A_true[[k]][i, ]
      p_emit  <- params_true[[i]][[k]]
    } else {
      # Idiosyncratic: asset-specific local state
      m_loc  <- sample(M, 1)
      s_idx  <- K + m_loc
      mu_val <- mu_true[[s_idx]][i]
      ar_val <- A_true[[s_idx]][i, ]
      p_emit <- params_true[[i]][[s_idx]]
    }
    
    eps_i   <- rst(1, xi = 0, omega = p_emit$sigma,
                   alpha = p_emit$alpha, nu = p_emit$nu)
    Y[t, i] <- mu_val + sum(ar_val * yt1) + eps_i
  }
}

colnames(Y) <- paste0("Asset_", seq_len(n))

cat("\nDescriptive statistics of simulated returns:\n")
print(round(apply(Y, 2, function(x)
  c(mean = mean(x), sd   = sd(x),
    skew = mean((x - mean(x))^3) / sd(x)^3,
    kurt = mean((x - mean(x))^4) / sd(x)^4 - 3)), 3))

# -----------------------------------------------------------------------------
# PLOT: RETURNS + TRUE REGIMES
# -----------------------------------------------------------------------------

par(mfrow = c(2, 1), mar = c(3, 4, 2, 1))
regime_col <- c("tomato", "grey60", "steelblue")

matplot(Y, type = "l", lty = 1,
        col = c("steelblue", "tomato", "forestgreen"),
        main = "Simulated Returns (3 assets)",
        ylab = "Return", xlab = "")
abline(h = 0, lty = 2, col = "grey80")
legend("topright", colnames(Y),
       col = c("steelblue", "tomato", "forestgreen"),
       lty = 1, bty = "n", cex = 0.8)

plot(z, type = "s", col = regime_col[z],
     main = "True Regime Sequence",
     ylab = "Regime", xlab = "Time", yaxt = "n", lwd = 1.5)
axis(2, at = 1:3, labels = c("Bear", "Stable", "Bull"))
par(mfrow = c(1, 1))

# -----------------------------------------------------------------------------
# ESTIMATE MODEL
# -----------------------------------------------------------------------------

cat("\n--- Starting EM estimation ---\n")

fit <- em_msvar(
  Y        = Y,
  K        = K,
  M        = M,
  max_iter = 100,
  tol      = 1e-5,
  P_zeros  = list(c(1, 3), c(3, 1)),
  verbose  = TRUE
)

# -----------------------------------------------------------------------------
# LOG-LIKELIHOOD CONVERGENCE
# -----------------------------------------------------------------------------

plot(fit$history, type = "l", col = "steelblue", lwd = 2,
     main = "EM Log-Likelihood Convergence",
     xlab = "Iteration", ylab = "Log-Likelihood")
abline(v = length(fit$history), lty = 2, col = "tomato")

# -----------------------------------------------------------------------------
# REGIME CLASSIFICATION ACCURACY (label-switching correction)
# -----------------------------------------------------------------------------

z_hat <- apply(fit$xi_filt, 1, which.max)

# Try all K! permutations
perms    <- combinat::permn(seq_len(K))
best_acc <- 0
best_map <- seq_len(K)

for (perm in perms) {
  acc <- mean(perm[z_hat] == z)
  if (acc > best_acc) {
    best_acc <- acc
    best_map <- perm
  }
}

z_aligned <- best_map[z_hat]
cat(sprintf("\nRegime classification accuracy: %.1f%%\n", best_acc * 100))

cat("\nConfusion matrix (rows=True, cols=Estimated):\n")
print(table(True = z, Estimated = z_aligned))

# Filtered probability plot
par(mfrow = c(3, 1), mar = c(2, 4, 2, 1))
state_names <- c("Bear", "Stable", "Bull")
for (k in seq_len(K)) {
  plot(fit$xi_filt[, best_map[k]], type = "l",
       col = regime_col[k], ylim = c(0, 1),
       main = paste("Filtered P(z_t =", state_names[k], ")"),
       ylab = "Probability", xlab = "")
  # Overlay true regime indicator
  lines(as.numeric(z == k) * 0.95, col = "grey40",
        lty = 2, lwd = 0.8)
}
par(mfrow = c(1, 1))

# -----------------------------------------------------------------------------
# TRANSITION MATRIX: COMPARISON + DEVIATION
# -----------------------------------------------------------------------------

# Re-permute estimated P to match label alignment
P_est_aligned        <- matrix(0, K, K)
for (j in seq_len(K)) {
  for (k in seq_len(K)) {
    P_est_aligned[best_map[j], best_map[k]] <- fit$P_trans[j, k]
  }
}

cat("\n--- Transition Matrix Comparison ---\n")
cat("\nTrue P:\n")
print(round(P_true, 4))
cat("\nEstimated P (label-aligned):\n")
print(round(P_est_aligned, 4))

# Element-wise absolute deviation
dev_matrix <- abs(P_true - P_est_aligned)
cat("\nAbsolute deviation |P_true - P_est|:\n")
print(round(dev_matrix, 4))

# Summary statistics of deviation
cat(sprintf("\nMax deviation  : %.4f\n", max(dev_matrix)))
cat(sprintf("Mean deviation : %.4f\n", mean(dev_matrix[P_true > 0])))
cat(sprintf("Frobenius norm : %.4f\n", norm(dev_matrix, type = "F")))

# Relative deviation (only for non-zero true entries)
rel_dev <- ifelse(P_true > 0, dev_matrix / P_true, NA)
cat("\nRelative deviation |P_true - P_est| / P_true:\n")
print(round(rel_dev, 4))
cat(sprintf("Mean relative deviation (non-zero): %.2f%%\n",
            100 * mean(rel_dev, na.rm = TRUE)))

# Visual: heatmap of deviation matrix
par(mfrow = c(1, 3), mar = c(4, 4, 3, 2))

image(t(P_true[K:1, ]), col = hcl.colors(20, "Blues", rev = TRUE),
      main = "True P", axes = FALSE)
axis(1, at = seq(0, 1, length = K), labels = state_names)
axis(2, at = seq(0, 1, length = K), labels = rev(state_names), las = 1)
for (j in seq_len(K)) for (k in seq_len(K))
  text((k-1)/(K-1), 1-(j-1)/(K-1),
       sprintf("%.3f", P_true[j, k]), cex = 1.1)

image(t(P_est_aligned[K:1, ]),
      col = hcl.colors(20, "Blues", rev = TRUE),
      main = "Estimated P", axes = FALSE)
axis(1, at = seq(0, 1, length = K), labels = state_names)
axis(2, at = seq(0, 1, length = K), labels = rev(state_names), las = 1)
for (j in seq_len(K)) for (k in seq_len(K))
  text((k-1)/(K-1), 1-(j-1)/(K-1),
       sprintf("%.3f", P_est_aligned[j, k]), cex = 1.1)

image(t(dev_matrix[K:1, ]),
      col = hcl.colors(20, "Reds", rev = TRUE),
      main = "|Deviation|", axes = FALSE)
axis(1, at = seq(0, 1, length = K), labels = state_names)
axis(2, at = seq(0, 1, length = K), labels = rev(state_names), las = 1)
for (j in seq_len(K)) for (k in seq_len(K))
  text((k-1)/(K-1), 1-(j-1)/(K-1),
       sprintf("%.3f", dev_matrix[j, k]), cex = 1.1)

par(mfrow = c(1, 1))

# -----------------------------------------------------------------------------
# AR COEFFICIENT DEVIATION per regime
# -----------------------------------------------------------------------------

cat("\n--- AR Matrix Deviation per Regime ---\n")

for (k in seq_len(K)) {
  k_est   <- which(best_map == k)
  A_dev   <- abs(A_true[[k]] - fit$A[[k_est]])
  cat(sprintf("\nRegime %s — Frobenius norm of A deviation: %.4f\n",
              state_names[k], norm(A_dev, type = "F")))
  cat("  True A:\n"); print(round(A_true[[k]], 3))
  cat("  Estimated A:\n"); print(round(fit$A[[k_est]], 3))
  cat("  |Deviation|:\n"); print(round(A_dev, 3))
}

# -----------------------------------------------------------------------------
# PARTICIPATION PROBABILITIES rho: TRUE vs ESTIMATED
# -----------------------------------------------------------------------------

cat("\n--- Participation Probabilities rho ---\n")
cat("True rho (n x K):\n")
print(round(rho_true, 3))
cat("Estimated rho (n x K, label-aligned):\n")
rho_aligned <- fit$rho[, best_map]
print(round(rho_aligned, 3))
cat("Absolute deviation:\n")
print(round(abs(rho_true - rho_aligned), 3))

# -----------------------------------------------------------------------------
# LEADERSHIP ANALYSIS
# -----------------------------------------------------------------------------

cat("\n--- Leadership Analysis ---\n")

W_list <- influence_matrix(fit$A[seq_len(K)], alpha_sig = 0.10, K = K)

for (k in seq_len(K)) {
  k_est <- which(best_map == k)
  cat(sprintf("\nInfluence matrix W^(%s):\n", state_names[k]))
  rownames(W_list[[k_est]]) <- colnames(W_list[[k_est]]) <- colnames(Y)
  print(round(W_list[[k_est]], 3))
}

h_mat <- dynamic_leadership(W_list, fit$xi_filt, K = K)
r_mat <- regime_pagerank(W_list, fit$xi_filt, K = K)
colnames(h_mat) <- colnames(r_mat) <- colnames(Y)

cat("\nMean dynamic leadership score h_{it}:\n")
print(round(colMeans(h_mat), 4))
cat("\nMean PageRank score r_t:\n")
print(round(colMeans(r_mat), 4))

# Leadership over time
matplot(h_mat, type = "l", lty = 1,
        col = c("steelblue", "tomato", "forestgreen"), lwd = 1.5,
        main = "Dynamic Leadership Score h_{it}",
        ylab = expression(h[it]), xlab = "Time")
legend("topright", colnames(Y),
       col = c("steelblue", "tomato", "forestgreen"),
       lty = 1, bty = "n", cex = 0.8)

# -----------------------------------------------------------------------------
# WALD TESTS
# -----------------------------------------------------------------------------

cat("\n--- Wald Tests — Bear regime ---\n")
wald_bear <- wald_test_leadership(Y, fit$xi_smooth, k = 1, lags_nw = 5)
cat("P-values [row i leads column l]:\n")
rownames(wald_bear$p_values) <- colnames(wald_bear$p_values) <- colnames(Y)
print(round(wald_bear$p_values, 3))

cat("\n--- Wald Tests — Bull regime ---\n")
wald_bull <- wald_test_leadership(Y, fit$xi_smooth, k = 3, lags_nw = 5)
cat("P-values [row i leads column l]:\n")
rownames(wald_bull$p_values) <- colnames(wald_bull$p_values) <- colnames(Y)
print(round(wald_bull$p_values, 3))

# -----------------------------------------------------------------------------
# BIC COMPARISON K=2 vs K=3
# -----------------------------------------------------------------------------

cat("\n--- Model Selection via BIC ---\n")

bic_k3 <- compute_bic(fit$log_lik, K = 3, n = n, p = 1, TT = TT)
cat(sprintf("BIC (K=3): %.2f\n", bic_k3))

fit_k2 <- em_msvar(Y, K = 2, M = M, max_iter = 60,
                   P_zeros = NULL, verbose = FALSE)
bic_k2 <- compute_bic(fit_k2$log_lik, K = 2, n = n, p = 1, TT = TT)
cat(sprintf("BIC (K=2): %.2f\n", bic_k2))
cat(sprintf("Selected model: K=%d\n", ifelse(bic_k3 < bic_k2, 3, 2)))





# Para cada regime, mostra true vs estimado para elementos off-diagonal
cat("=== Verificação de liderança por coeficientes A_k ===\n")

for (k in seq_len(K)) {
  k_est <- which(best_map == k)
  cat(sprintf("\n-- Regime %s --\n", state_names[k]))
  
  A_t <- A_true[[k]]
  A_e <- fit$A[[k_est]]
  
  # Extrai apenas off-diagonais: par (i -> l) significa A[l, i]
  df_coef <- data.frame(
    from     = integer(0),
    to       = integer(0),
    true_val = numeric(0),
    est_val  = numeric(0),
    abs_dev  = numeric(0),
    same_sign = logical(0)
  )
  
  for (i in seq_len(n)) {
    for (l in seq_len(n)) {
      if (i == l) next
      df_coef <- rbind(df_coef, data.frame(
        from      = i,
        to        = l,
        true_val  = round(A_t[l, i], 4),
        est_val   = round(A_e[l, i], 4),
        abs_dev   = round(abs(A_t[l, i] - A_e[l, i]), 4),
        same_sign = sign(A_t[l, i]) == sign(A_e[l, i])
      ))
    }
  }
  print(df_coef)
}

cat("\n=== Verificação do ranking de out-strength ===\n")

for (k in seq_len(K)) {
  k_est <- which(best_map == k)
  
  # Out-strength verdadeiro: soma dos |A_true[l,i]| para l != i
  os_true <- sapply(seq_len(n), function(i)
    sum(abs(A_true[[k]][-i, i])))
  
  # Out-strength estimado
  os_est <- rowSums(W_list[[k_est]])
  
  rank_true <- rank(-os_true)
  rank_est  <- rank(-os_est)
  
  cat(sprintf("\nRegime %s:\n", state_names[k]))
  print(data.frame(
    asset       = colnames(Y),
    os_true     = round(os_true, 4),
    rank_true   = rank_true,
    os_est      = round(os_est, 4),
    rank_est    = rank_est,
    rank_match  = rank_true == rank_est
  ))
}

cat("\n=== Correlação de Spearman dos rankings ===\n")

for (k in seq_len(K)) {
  k_est <- which(best_map == k)
  
  os_true <- sapply(seq_len(n), function(i)
    sum(abs(A_true[[k]][-i, i])))
  os_est  <- rowSums(W_list[[k_est]])
  
  rho_sp <- cor(os_true, os_est, method = "spearman")
  cat(sprintf("Regime %s — Spearman rho: %.4f\n", state_names[k], rho_sp))
}

cat("\n=== Acurácia do teste de Wald vs. verdade ===\n")

alpha_test <- 0.10

wald_results <- list(
  wald_test_leadership(Y, fit$xi_smooth, k = 1, lags_nw = 5),
  wald_test_leadership(Y, fit$xi_smooth, k = 2, lags_nw = 5),
  wald_test_leadership(Y, fit$xi_smooth, k = 3, lags_nw = 5)
)

for (k in seq_len(K)) {
  k_est  <- which(best_map == k)
  pv     <- wald_results[[k]]$p_values
  A_t    <- A_true[[k]]
  
  tp <- fp <- tn <- fn <- 0
  
  for (i in seq_len(n)) {
    for (l in seq_len(n)) {
      if (i == l || is.na(pv[i, l])) next
      true_nonzero <- abs(A_t[l, i]) > 1e-6
      rejected     <- pv[i, l] < alpha_test
      
      if (true_nonzero  &&  rejected)  tp <- tp + 1
      if (!true_nonzero &&  rejected)  fp <- fp + 1
      if (!true_nonzero && !rejected)  tn <- tn + 1
      if (true_nonzero  && !rejected)  fn <- fn + 1
    }
  }
  
  precision <- ifelse(tp + fp > 0, tp / (tp + fp), NA)
  recall    <- ifelse(tp + fn > 0, tp / (tp + fn), NA)
  f1        <- ifelse(!is.na(precision) & !is.na(recall) &
                        (precision + recall) > 0,
                      2 * precision * recall / (precision + recall), NA)
  
  cat(sprintf(
    "\nRegime %s | TP=%d FP=%d TN=%d FN=%d | Precision=%.2f Recall=%.2f F1=%.2f\n",
    state_names[k], tp, fp, tn, fn,
    ifelse(is.na(precision), 0, precision),
    ifelse(is.na(recall), 0, recall),
    ifelse(is.na(f1), 0, f1)
  ))
}

par(mfrow = c(2, K), mar = c(3, 3, 3, 1))

for (k in seq_len(K)) {
  k_est  <- which(best_map == k)
  
  # W verdadeiro (sem threshold)
  W_true_k <- matrix(0, n, n)
  for (i in seq_len(n))
    for (l in seq_len(n))
      if (i != l) W_true_k[i, l] <- abs(A_true[[k]][l, i])
  
  image(t(W_true_k[n:1, ]),
        col  = hcl.colors(20, "Blues", rev = TRUE),
        main = sprintf("W true — %s", state_names[k]),
        axes = FALSE)
  axis(1, at = seq(0,1,len=n), labels = colnames(Y), cex.axis = 0.9)
  axis(2, at = seq(0,1,len=n), labels = rev(colnames(Y)),
       las = 1, cex.axis = 0.9)
  for (i in seq_len(n)) for (l in seq_len(n))
    text((l-1)/(n-1), 1-(i-1)/(n-1),
         sprintf("%.2f", W_true_k[i,l]), cex = 0.9)
}

for (k in seq_len(K)) {
  k_est <- which(best_map == k)
  Wk    <- W_list[[k_est]]
  
  image(t(Wk[n:1, ]),
        col  = hcl.colors(20, "Blues", rev = TRUE),
        main = sprintf("W est — %s", state_names[k]),
        axes = FALSE)
  axis(1, at = seq(0,1,len=n), labels = colnames(Y), cex.axis = 0.9)
  axis(2, at = seq(0,1,len=n), labels = rev(colnames(Y)),
       las = 1, cex.axis = 0.9)
  for (i in seq_len(n)) for (l in seq_len(n))
    text((l-1)/(n-1), 1-(i-1)/(n-1),
         sprintf("%.2f", Wk[i,l]), cex = 0.9)
}

par(mfrow = c(1, 1))