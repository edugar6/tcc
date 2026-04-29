##---------------
rm(list = ls())

current_working_dir = dirname(rstudioapi::getSourceEditorContext()$path)
setwd(current_working_dir)

# =============================================================================
# simulation_example.R — versão corrigida
# Parâmetros com maior separação entre regimes, T=1500, warm start via K-means
# =============================================================================

source("utils.R")
library(ggplot2)
library(reshape2)
library(combinat)

set.seed(2024)

# -----------------------------------------------------------------------------
# PARÂMETROS DA SIMULAÇÃO
# -----------------------------------------------------------------------------

TT <- 1500
n  <- 3
K  <- 3    # bear=1, stable=2, bull=3
M  <- 2

# Matriz de transição verdadeira
P_true <- matrix(c(
  0.85, 0.15, 0.00,
  0.05, 0.90, 0.05,
  0.00, 0.15, 0.85
), nrow = K, byrow = TRUE)

# Interceptos — maior separação entre regimes
mu_true <- list(
  c(-0.60, -0.50, -0.70),   # bear
  c( 0.00,  0.01, -0.01),   # stable
  c( 0.60,  0.50,  0.70),   # bull
  c(-0.05, -0.04, -0.06),   # idiossincrático 1 (low vol, ~zero)
  c( 0.05,  0.04,  0.06)    # idiossincrático 2 (low vol, ~zero)
)

# Matrizes AR — coeficientes off-diagonal pronunciados
# Bear: Asset_1 lidera Asset_2 e Asset_3
A_bear <- matrix(c(
  0.15,  0.00,  0.00,
  0.40,  0.08,  0.00,
  0.35,  0.00,  0.08
), nrow = n, byrow = TRUE)

# Stable: dinâmica fraca, sem liderança clara
A_stable <- diag(0.05, n)

# Bull: Asset_3 lidera Asset_1 e Asset_2
A_bull <- matrix(c(
  0.08,  0.00,  0.35,
  0.00,  0.08,  0.40,
  0.00,  0.00,  0.15
), nrow = n, byrow = TRUE)

# Idiossincráticos: fraco, sem estrutura de liderança
A_idio1 <- diag(0.03, n)
A_idio2 <- diag(0.03, n)

A_true <- list(A_bear, A_stable, A_bull, A_idio1, A_idio2)

# Probabilidades de participação verdadeiras
rho_true <- matrix(c(
  0.85, 0.90, 0.85,
  0.80, 0.90, 0.82,
  0.88, 0.92, 0.87
), nrow = n, byrow = TRUE)

# Parâmetros de emissão skew-t — volatilidades bem separadas
sigma_glob <- list(
  bear   = c(3.0, 2.8, 3.2),
  stable = c(0.5, 0.4, 0.6),
  bull   = c(1.2, 1.0, 1.4)
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

# Idiossincráticos: vol baixa e claramente distinta dos globais
sigma_idio <- matrix(c(
  0.8, 0.7, 0.9,
  0.8, 0.7, 0.9
), nrow = M, byrow = TRUE)
alpha_idio <- matrix(0,  nrow = M, ncol = n)
nu_idio    <- matrix(10, nrow = M, ncol = n)

# Monta params_emit[[i]][[s]] — asset-specific
params_true <- lapply(seq_len(n), function(i) {
  states      <- vector("list", K + M)
  states[[1]] <- list(sigma = sigma_glob$bear[i],
                      alpha = alpha_glob$bear[i],
                      nu    = nu_glob$bear[i])
  states[[2]] <- list(sigma = sigma_glob$stable[i],
                      alpha = alpha_glob$stable[i],
                      nu    = nu_glob$stable[i])
  states[[3]] <- list(sigma = sigma_glob$bull[i],
                      alpha = alpha_glob$bull[i],
                      nu    = nu_glob$bull[i])
  for (m in seq_len(M)) {
    states[[K + m]] <- list(sigma = sigma_idio[m, i],
                            alpha = alpha_idio[m, i],
                            nu    = nu_idio[m, i])
  }
  states
})

# -----------------------------------------------------------------------------
# SIMULA SEQUÊNCIA DE REGIMES z_t
# -----------------------------------------------------------------------------

z    <- integer(TT)
z[1] <- 2

for (t in 2:TT) {
  z[t] <- sample(K, 1, prob = P_true[z[t - 1], ])
}

cat("Frequências de regime verdadeiras:\n")
print(round(table(z) / TT, 3))

# -----------------------------------------------------------------------------
# SIMULA RETORNOS Y
# -----------------------------------------------------------------------------

Y <- matrix(0, TT, n)

for (t in 2:TT) {
  k   <- z[t]
  yt1 <- Y[t - 1, ]
  
  for (i in seq_len(n)) {
    s_it <- rbinom(1, 1, rho_true[i, k])
    
    if (s_it == 1) {
      mu_val <- mu_true[[k]][i]
      ar_val <- A_true[[k]][i, ]
      p_emit <- params_true[[i]][[k]]
    } else {
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

cat("\nEstatísticas descritivas dos retornos simulados:\n")
print(round(apply(Y, 2, function(x)
  c(mean = mean(x), sd   = sd(x),
    skew = mean((x - mean(x))^3) / sd(x)^3,
    kurt = mean((x - mean(x))^4) / sd(x)^4 - 3)), 3))

# -----------------------------------------------------------------------------
# VISUALIZAÇÃO: RETORNOS + REGIMES VERDADEIROS
# -----------------------------------------------------------------------------

par(mfrow = c(2, 1), mar = c(3, 4, 2, 1))
regime_col <- c("tomato", "grey60", "steelblue")

matplot(Y, type = "l", lty = 1,
        col = c("steelblue", "tomato", "forestgreen"),
        main = "Retornos Simulados (3 ativos)",
        ylab = "Retorno", xlab = "")
abline(h = 0, lty = 2, col = "grey80")
legend("topright", colnames(Y),
       col = c("steelblue", "tomato", "forestgreen"),
       lty = 1, bty = "n", cex = 0.8)

plot(z, type = "s", col = regime_col[z],
     main = "Sequência de Regimes Verdadeira",
     ylab = "Regime", xlab = "Tempo", yaxt = "n", lwd = 1.5)
axis(2, at = 1:3, labels = c("Bear", "Stable", "Bull"))
par(mfrow = c(1, 1))

# -----------------------------------------------------------------------------
# INICIALIZAÇÃO VIA K-MEANS (warm start)
# -----------------------------------------------------------------------------

cat("\nComputando inicialização via K-means...\n")

km            <- kmeans(Y, centers = K, nstart = 50, iter.max = 100)
cluster_order <- order(km$centers[, 1])   # ordena por retorno médio crescente

init_mu <- lapply(seq_len(K + M), function(s) {
  if (s <= K) as.numeric(km$centers[cluster_order[s], ])
  else        rep(0, n)
})

init_A <- lapply(seq_len(K + M), function(s) {
  if (s > K) return(diag(0.02, n))
  idx <- which(km$cluster == cluster_order[s])
  idx <- idx[idx > 1]
  if (length(idx) < n + 2) return(diag(0.05, n))
  Xc  <- cbind(1, Y[idx - 1, , drop = FALSE])
  Yc  <- Y[idx, , drop = FALSE]
  B   <- tryCatch(
    solve(t(Xc) %*% Xc + diag(1e-6, ncol(Xc)), t(Xc) %*% Yc),
    error = function(e) matrix(0, ncol(Xc), n)
  )
  t(B[2:(n + 1), , drop = FALSE])
})

# Inicializa rho alto (próximo do verdadeiro)
init_rho <- matrix(0.85, n, K)

# Inicializa P_trans com zeros estruturais
init_P <- matrix(c(
  0.80, 0.20, 0.00,
  0.10, 0.80, 0.10,
  0.00, 0.20, 0.80
), nrow = K, byrow = TRUE)

# Inicializa params_emit com sigma calibrado por cluster
init_params <- lapply(seq_len(n), function(i) {
  states <- vector("list", K + M)
  for (s in seq_len(K)) {
    idx_s  <- which(km$cluster == cluster_order[s])
    sig_s  <- if (length(idx_s) > 1) sd(Y[idx_s, i]) else 1.0
    states[[s]] <- list(sigma = sig_s, alpha = 0.0, nu = 5)
  }
  for (m in seq_len(M)) {
    states[[K + m]] <- list(sigma = 0.8, alpha = 0.0, nu = 10)
  }
  states
})

# -----------------------------------------------------------------------------
# ESTIMAÇÃO DO MODELO
# -----------------------------------------------------------------------------

cat("\n--- Iniciando estimação EM ---\n")

fit <- em_msvar(
  Y        = Y,
  K        = K,
  M        = M,
  max_iter = 150,
  tol      = 1e-6,
  P_zeros  = list(c(1, 3), c(3, 1)),
  init     = list(
    P_trans     = init_P,
    pi0         = rep(1 / K, K),
    mu          = init_mu,
    A           = init_A,
    rho         = init_rho,
    params_emit = init_params
  ),
  verbose  = TRUE
)

# -----------------------------------------------------------------------------
# CONVERGÊNCIA
# -----------------------------------------------------------------------------

plot(fit$history, type = "l", col = "steelblue", lwd = 2,
     main = "Convergência da Log-verossimilhança (EM)",
     xlab = "Iteração", ylab = "Log-verossimilhança")
abline(v = length(fit$history), lty = 2, col = "tomato")

# -----------------------------------------------------------------------------
# ACURÁCIA DE CLASSIFICAÇÃO (correção de label switching)
# -----------------------------------------------------------------------------

state_names <- c("Bear", "Stable", "Bull")
z_hat       <- apply(fit$xi_filt, 1, which.max)

perms    <- combinat::permn(seq_len(K))
best_acc <- 0
best_map <- seq_len(K)

for (perm in perms) {
  acc <- mean(perm[z_hat] == z)
  if (acc > best_acc) { best_acc <- acc; best_map <- perm }
}

z_aligned <- best_map[z_hat]
cat(sprintf("\nAcurácia de classificação: %.1f%%\n", best_acc * 100))

cat("\nMatriz de confusão (linhas=Verdadeiro, colunas=Estimado):\n")
print(table(True = z, Estimated = z_aligned))

# Probabilidades filtradas por regime
par(mfrow = c(3, 1), mar = c(2, 4, 2, 1))
for (k in seq_len(K)) {
  plot(fit$xi_filt[, best_map[k]], type = "l",
       col = regime_col[k], ylim = c(0, 1),
       main = paste("P filtrada(z_t =", state_names[k], ")"),
       ylab = "Probabilidade", xlab = "")
  lines(as.numeric(z == k) * 0.95, col = "grey40",
        lty = 2, lwd = 0.8)
}
par(mfrow = c(1, 1))

# -----------------------------------------------------------------------------
# COMPARAÇÃO DA MATRIZ DE TRANSIÇÃO
# -----------------------------------------------------------------------------

P_est_aligned <- matrix(0, K, K)
for (j in seq_len(K))
  for (k in seq_len(K))
    P_est_aligned[best_map[j], best_map[k]] <- fit$P_trans[j, k]

cat("\n--- Comparação da Matriz de Transição ---\n")
cat("\nP verdadeira:\n");        print(round(P_true, 4))
cat("\nP estimada (alinhada):\n"); print(round(P_est_aligned, 4))

dev_matrix <- abs(P_true - P_est_aligned)
cat("\nDesvio absoluto |P_true - P_est|:\n"); print(round(dev_matrix, 4))
cat(sprintf("\nNorma de Frobenius : %.4f\n", norm(dev_matrix, type = "F")))
cat(sprintf("Desvio máximo      : %.4f\n", max(dev_matrix)))
cat(sprintf("Desvio médio       : %.4f\n", mean(dev_matrix[P_true > 0])))

rel_dev <- ifelse(P_true > 0, dev_matrix / P_true, NA)
cat(sprintf("Desvio relativo médio: %.2f%%\n",
            100 * mean(rel_dev, na.rm = TRUE)))

# Heatmaps: P_true | P_est | |Desvio|
par(mfrow = c(1, 3), mar = c(4, 4, 3, 2))
plot_P_heatmap <- function(mat, title) {
  image(t(mat[K:1, ]), col = hcl.colors(20, "Blues", rev = TRUE),
        main = title, axes = FALSE)
  axis(1, at = seq(0, 1, len = K), labels = state_names)
  axis(2, at = seq(0, 1, len = K), labels = rev(state_names), las = 1)
  for (j in seq_len(K)) for (k in seq_len(K))
    text((k-1)/(K-1), 1-(j-1)/(K-1),
         sprintf("%.3f", mat[j, k]), cex = 1.1)
}
plot_P_heatmap(P_true,        "P verdadeira")
plot_P_heatmap(P_est_aligned, "P estimada")

image(t(dev_matrix[K:1, ]), col = hcl.colors(20, "Reds", rev = TRUE),
      main = "|Desvio|", axes = FALSE)
axis(1, at = seq(0, 1, len = K), labels = state_names)
axis(2, at = seq(0, 1, len = K), labels = rev(state_names), las = 1)
for (j in seq_len(K)) for (k in seq_len(K))
  text((k-1)/(K-1), 1-(j-1)/(K-1),
       sprintf("%.3f", dev_matrix[j, k]), cex = 1.1)
par(mfrow = c(1, 1))

# -----------------------------------------------------------------------------
# DESVIO DAS MATRIZES A_k
# -----------------------------------------------------------------------------

cat("\n--- Desvio das Matrizes A_k por Regime ---\n")

for (k in seq_len(K)) {
  k_est <- which(best_map == k)
  A_dev <- abs(A_true[[k]] - fit$A[[k_est]])
  cat(sprintf("\nRegime %s — Frobenius: %.4f\n",
              state_names[k], norm(A_dev, type = "F")))
  cat("  A verdadeira:\n");   print(round(A_true[[k]], 3))
  cat("  A estimada:\n");     print(round(fit$A[[k_est]], 3))
  cat("  |Desvio|:\n");       print(round(A_dev, 3))
}

# -----------------------------------------------------------------------------
# RHO: VERDADEIRO vs. ESTIMADO
# -----------------------------------------------------------------------------

cat("\n--- Probabilidades de Participação rho ---\n")
rho_aligned <- fit$rho[, best_map]
cat("rho verdadeiro:\n");   print(round(rho_true, 3))
cat("rho estimado:\n");     print(round(rho_aligned, 3))
cat("Desvio absoluto:\n");  print(round(abs(rho_true - rho_aligned), 3))

# -----------------------------------------------------------------------------
# ANÁLISE DE LIDERANÇA
# -----------------------------------------------------------------------------

cat("\n--- Análise de Liderança ---\n")

W_list <- influence_matrix(fit$A[seq_len(K)], alpha_sig = 0.10, K = K)

for (k in seq_len(K)) {
  k_est <- which(best_map == k)
  rownames(W_list[[k_est]]) <- colnames(W_list[[k_est]]) <- colnames(Y)
  cat(sprintf("\nW^(%s):\n", state_names[k]))
  print(round(W_list[[k_est]], 3))
}

h_mat <- dynamic_leadership(W_list, fit$xi_filt, K = K)
r_mat <- regime_pagerank(W_list, fit$xi_filt, K = K)
colnames(h_mat) <- colnames(r_mat) <- colnames(Y)

cat("\nScore médio de liderança h_{it}:\n")
print(round(colMeans(h_mat), 4))
cat("\nPageRank médio r_t:\n")
print(round(colMeans(r_mat), 4))

# -----------------------------------------------------------------------------
# VERIFICAÇÃO 1 — COEFICIENTES A_k
# -----------------------------------------------------------------------------

cat("\n=== Verificação de liderança por coeficientes A_k ===\n")

for (k in seq_len(K)) {
  k_est <- which(best_map == k)
  A_t   <- A_true[[k]]
  A_e   <- fit$A[[k_est]]
  cat(sprintf("\n-- Regime %s --\n", state_names[k]))
  
  rows <- list()
  for (i in seq_len(n)) for (l in seq_len(n)) {
    if (i == l) next
    rows[[length(rows) + 1]] <- data.frame(
      from      = i, to = l,
      true_val  = round(A_t[l, i], 4),
      est_val   = round(A_e[l, i], 4),
      abs_dev   = round(abs(A_t[l, i] - A_e[l, i]), 4),
      same_sign = sign(A_t[l, i]) == sign(A_e[l, i])
    )
  }
  print(do.call(rbind, rows))
}

# -----------------------------------------------------------------------------
# VERIFICAÇÃO 2 — RANKING DE OUT-STRENGTH
# -----------------------------------------------------------------------------

cat("\n=== Ranking de out-strength ===\n")

for (k in seq_len(K)) {
  k_est   <- which(best_map == k)
  os_true <- sapply(seq_len(n), function(i) sum(abs(A_true[[k]][-i, i])))
  os_est  <- rowSums(W_list[[k_est]])
  cat(sprintf("\nRegime %s:\n", state_names[k]))
  print(data.frame(
    asset      = colnames(Y),
    os_true    = round(os_true, 4),
    rank_true  = rank(-os_true),
    os_est     = round(os_est, 4),
    rank_est   = rank(-os_est),
    rank_match = rank(-os_true) == rank(-os_est)
  ))
}

# -----------------------------------------------------------------------------
# VERIFICAÇÃO 3 — CORRELAÇÃO DE SPEARMAN
# -----------------------------------------------------------------------------

cat("\n=== Correlação de Spearman dos rankings ===\n")

for (k in seq_len(K)) {
  k_est   <- which(best_map == k)
  os_true <- sapply(seq_len(n), function(i) sum(abs(A_true[[k]][-i, i])))
  os_est  <- rowSums(W_list[[k_est]])
  rho_sp  <- suppressWarnings(cor(os_true, os_est, method = "spearman"))
  cat(sprintf("Regime %s — Spearman rho: %s\n",
              state_names[k],
              ifelse(is.na(rho_sp), "NA (desvio zero)", round(rho_sp, 4))))
}

# -----------------------------------------------------------------------------
# VERIFICAÇÃO 4 — WALD + PRECISION / RECALL / F1
# -----------------------------------------------------------------------------

cat("\n=== Acurácia do teste de Wald vs. verdade ===\n")

alpha_test   <- 0.10
wald_results <- lapply(seq_len(K), function(k)
  wald_test_leadership(Y, fit$xi_smooth, k = k, lags_nw = 5))

for (k in seq_len(K)) {
  k_est <- which(best_map == k)
  pv    <- wald_results[[k]]$p_values
  A_t   <- A_true[[k_est]]
  tp <- fp <- tn <- fn <- 0
  
  for (i in seq_len(n)) for (l in seq_len(n)) {
    if (i == l || is.na(pv[i, l])) next
    true_nz  <- abs(A_t[l, i]) > 1e-6
    rejected <- pv[i, l] < alpha_test
    if ( true_nz &&  rejected) tp <- tp + 1
    if (!true_nz &&  rejected) fp <- fp + 1
    if (!true_nz && !rejected) tn <- tn + 1
    if ( true_nz && !rejected) fn <- fn + 1
  }
  
  prec <- ifelse(tp + fp > 0, tp / (tp + fp), 0)
  rec  <- ifelse(tp + fn > 0, tp / (tp + fn), 0)
  f1   <- ifelse(prec + rec > 0, 2 * prec * rec / (prec + rec), 0)
  cat(sprintf(
    "Regime %s | TP=%d FP=%d TN=%d FN=%d | Precision=%.2f Recall=%.2f F1=%.2f\n",
    state_names[k], tp, fp, tn, fn, prec, rec, f1))
}

# -----------------------------------------------------------------------------
# VERIFICAÇÃO 5 — HEATMAP W^(k) VERDADEIRO vs. ESTIMADO
# -----------------------------------------------------------------------------

par(mfrow = c(2, K), mar = c(3, 3, 3, 1))

for (k in seq_len(K)) {
  W_true_k <- matrix(0, n, n)
  for (i in seq_len(n)) for (l in seq_len(n))
    if (i != l) W_true_k[i, l] <- abs(A_true[[k]][l, i])
  
  image(t(W_true_k[n:1, ]), col = hcl.colors(20, "Blues", rev = TRUE),
        main = sprintf("W true — %s", state_names[k]), axes = FALSE)
  axis(1, at = seq(0, 1, len = n), labels = colnames(Y), cex.axis = 0.9)
  axis(2, at = seq(0, 1, len = n), labels = rev(colnames(Y)),
       las = 1, cex.axis = 0.9)
  for (i in seq_len(n)) for (l in seq_len(n))
    text((l-1)/(n-1), 1-(i-1)/(n-1),
         sprintf("%.2f", W_true_k[i, l]), cex = 0.9)
}

for (k in seq_len(K)) {
  k_est <- which(best_map == k)
  Wk    <- W_list[[k_est]]
  image(t(Wk[n:1, ]), col = hcl.colors(20, "Blues", rev = TRUE),
        main = sprintf("W est — %s", state_names[k]), axes = FALSE)
  axis(1, at = seq(0, 1, len = n), labels = colnames(Y), cex.axis = 0.9)
  axis(2, at = seq(0, 1, len = n), labels = rev(colnames(Y)),
       las = 1, cex.axis = 0.9)
  for (i in seq_len(n)) for (l in seq_len(n))
    text((l-1)/(n-1), 1-(i-1)/(n-1),
         sprintf("%.2f", Wk[i, l]), cex = 0.9)
}
par(mfrow = c(1, 1))

# -----------------------------------------------------------------------------
# HIERARQUIA DINÂMICA AO LONGO DO TEMPO
# -----------------------------------------------------------------------------

regime_dominant <- apply(fit$xi_filt, 1, which.max)
regime_label    <- factor(best_map[regime_dominant],
                          levels = 1:K, labels = state_names)
df_regime <- data.frame(t = seq_len(TT), regime = regime_label)
regime_fills <- c(Bear = "#FDDCDC", Stable = "#E8E8E8", Bull = "#D6EAF8")

# Score h_{it}
df_h <- as.data.frame(h_mat)
df_h$t <- seq_len(TT)
df_h_long <- melt(df_h, id.vars = "t",
                  variable.name = "Asset", value.name = "Score")

p1 <- ggplot() +
  geom_tile(data = df_regime,
            aes(x = t, y = 0, fill = regime, height = Inf),
            alpha = 0.35, width = 1) +
  scale_fill_manual(values = regime_fills, name = "Regime") +
  geom_line(data = df_h_long,
            aes(x = t, y = Score, colour = Asset), linewidth = 0.7) +
  scale_colour_manual(values = c("steelblue", "tomato", "forestgreen")) +
  labs(title    = "Dynamic Leadership Score h_{it}",
       subtitle = "Fundo = regime dominante",
       x = "Tempo", y = expression(h[it])) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())

print(p1)

# Ranking dinâmico (bump chart)
rank_mat <- t(apply(h_mat, 1, function(row) rank(-row, ties.method = "min")))
colnames(rank_mat) <- colnames(Y)
df_rank <- as.data.frame(rank_mat)
df_rank$t <- seq_len(TT)
df_rank_long <- melt(df_rank, id.vars = "t",
                     variable.name = "Asset", value.name = "Rank")

p2 <- ggplot() +
  geom_tile(data = df_regime,
            aes(x = t, y = 0, fill = regime, height = Inf),
            alpha = 0.30, width = 1) +
  scale_fill_manual(values = regime_fills, name = "Regime") +
  geom_line(data = df_rank_long,
            aes(x = t, y = Rank, colour = Asset, group = Asset),
            linewidth = 0.8) +
  geom_point(data = df_rank_long[df_rank_long$t %% 100 == 0, ],
             aes(x = t, y = Rank, colour = Asset), size = 2) +
  scale_colour_manual(values = c("steelblue", "tomato", "forestgreen")) +
  scale_y_reverse(breaks = seq_len(n),
                  labels = paste("Rank", seq_len(n))) +
  labs(title    = "Ranking Dinâmico de Liderança",
       subtitle = "Rank 1 = maior liderança no período t",
       x = "Tempo", y = "Ranking") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())

print(p2)

# Proporção como líder por regime
df_rank_regime <- merge(df_rank_long, df_regime, by = "t")
prop_top <- aggregate(Rank == 1 ~ Asset + regime,
                      data = df_rank_regime, FUN = mean)
names(prop_top)[3] <- "Prop_Rank1"

p3 <- ggplot(prop_top, aes(x = regime, y = Prop_Rank1, fill = Asset)) +
  geom_col(position = "dodge", colour = "white", linewidth = 0.4) +
  scale_fill_manual(values = c("steelblue", "tomato", "forestgreen")) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(title    = "Proporção do Tempo como Líder por Regime",
       subtitle = "Fração de períodos com Rank 1 por ativo e regime",
       x = "Regime", y = "% dos períodos") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

print(p3)

# Heatmap suavizado de h_{it}
smooth_h <- function(x, w = 30)
  as.numeric(stats::filter(x, rep(1/w, w), sides = 2))

df_hs <- as.data.frame(apply(h_mat, 2, smooth_h, w = 30))
df_hs$t <- seq_len(TT)
df_hs   <- na.omit(df_hs)
df_hs_long <- melt(df_hs, id.vars = "t",
                   variable.name = "Asset", value.name = "Score_smooth")

transitions <- which(diff(as.integer(regime_label)) != 0)

p4 <- ggplot(df_hs_long, aes(x = t, y = Asset, fill = Score_smooth)) +
  geom_tile() +
  scale_fill_gradient2(low  = "white", mid = "steelblue", high = "navy",
                       midpoint = median(df_hs_long$Score_smooth),
                       name = expression(h[it])) +
  geom_vline(xintercept = transitions,
             colour = "grey40", linewidth = 0.25, alpha = 0.5) +
  labs(title    = "Heatmap do Score de Liderança (MA 30 períodos)",
       subtitle = "Linhas verticais = transições de regime",
       x = "Tempo", y = "Ativo") +
  theme_bw(base_size = 12) +
  theme(legend.position = "right")

print(p4)

# Sumário numérico de estabilidade
cat("\n=== Estabilidade da hierarquia ao longo do tempo ===\n")

inversions <- mean(apply(rank_mat, 1, function(r) any(diff(r) != 0)))
cat(sprintf("Frequência de mudanças de ranking: %.1f%%\n", 100 * inversions))

cat("\nFração do tempo em cada rank por ativo:\n")
for (a in seq_len(n)) {
  cat(sprintf("  %s: ", colnames(Y)[a]))
  for (r in seq_len(n))
    cat(sprintf("Rank%d=%.1f%%  ", r, 100 * mean(rank_mat[, a] == r)))
  cat("\n")
}

cat("\nFração como líder (Rank 1) por regime:\n")
print(round(
  tapply(df_rank_regime$Rank == 1,
         list(df_rank_regime$Asset, df_rank_regime$regime), mean) * 100, 1))

# -----------------------------------------------------------------------------
# BIC: K=2 vs K=3
# -----------------------------------------------------------------------------

cat("\n--- Seleção de Modelo via BIC ---\n")
bic_k3 <- compute_bic(fit$log_lik, K = 3, n = n, p = 1, TT = TT)
cat(sprintf("BIC (K=3): %.2f\n", bic_k3))

fit_k2 <- em_msvar(Y, K = 2, M = M, max_iter = 80,
                   P_zeros = NULL, verbose = FALSE)
bic_k2 <- compute_bic(fit_k2$log_lik, K = 2, n = n, p = 1, TT = TT)
cat(sprintf("BIC (K=2): %.2f\n", bic_k2))
cat(sprintf("Modelo selecionado: K=%d\n", ifelse(bic_k3 < bic_k2, 3, 2)))