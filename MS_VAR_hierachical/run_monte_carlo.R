##---------------
rm(list = ls())

setwd("C:\\Users\\egarc\\OneDrive\\Documentos\\Edu\\IA\\tcc\\MS_VAR_hierachical")

# =============================================================================
# monte_carlo.R
# Estudo de Monte Carlo para o MS-VAR com hierarquia dinâmica
# Compatível com Windows, Linux e macOS
# Resultados salvos a cada lote — pode ser interrompido e retomado
# =============================================================================

source("utils.R")
library(combinat)
library(parallel)
library(ggplot2)
library(reshape2)

# -----------------------------------------------------------------------------
# CONFIGURAÇÃO DO EXPERIMENTO
# -----------------------------------------------------------------------------

MC_CONFIG <- list(
  R        = 500,
  T_grid   = c(500, 750, 1000, 1500, 2000),
  sig_grid = c("fraco", "moderado", "forte"),
  rho_grid = c(0.70, 0.80, 0.90, 0.95),
  K        = 3,
  M        = 2,
  n        = 3,
  n_cores  = max(1L, detectCores() - 1L),
  out_dir  = "mc_results",
  seed     = 42
)

dir.create(MC_CONFIG$out_dir, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# PARÂMETROS FIXOS
# -----------------------------------------------------------------------------

P_TRUE <- matrix(c(
  0.85, 0.15, 0.00,
  0.05, 0.90, 0.05,
  0.00, 0.15, 0.85
), nrow = 3, byrow = TRUE)

MU_TRUE <- list(
  c(-0.60, -0.50, -0.70),
  c( 0.00,  0.01, -0.01),
  c( 0.60,  0.50,  0.70),
  c(-0.05, -0.04, -0.06),
  c( 0.05,  0.04,  0.06)
)

AR_OFFDIAG <- list(
  fraco    = list(lead = c(0.10, 0.15), follow = c(0.05, 0.08)),
  moderado = list(lead = c(0.35, 0.40), follow = c(0.08, 0.08)),
  forte    = list(lead = c(0.45, 0.55), follow = c(0.10, 0.10))
)

SIGMA_GLOB <- list(
  bear   = c(3.0, 2.8, 3.2),
  stable = c(0.5, 0.4, 0.6),
  bull   = c(1.2, 1.0, 1.4)
)
ALPHA_GLOB <- list(
  bear   = c(-0.6, -0.5, -0.7),
  stable = c( 0.0,  0.0,  0.0),
  bull   = c( 0.4,  0.3,  0.5)
)
NU_GLOB <- list(
  bear   = c(5, 6, 5),
  stable = c(10, 10, 10),
  bull   = c(7, 8, 6)
)

# -----------------------------------------------------------------------------
# FUNÇÕES AUXILIARES
# -----------------------------------------------------------------------------

make_A_true <- function(signal) {
  od <- AR_OFFDIAG[[signal]]
  n  <- 3
  
  A_bear <- matrix(c(
    0.15,         0.00,           0.00,
    od$lead[2],   od$follow[1],   0.00,
    od$lead[1],   0.00,           od$follow[1]
  ), nrow = n, byrow = TRUE)
  
  A_stable <- diag(0.05, n)
  
  A_bull <- matrix(c(
    od$follow[1],  0.00,          od$lead[1],
    0.00,          od$follow[1],  od$lead[2],
    0.00,          0.00,          0.15
  ), nrow = n, byrow = TRUE)
  
  list(A_bear, A_stable, A_bull, diag(0.03, n), diag(0.03, n))
}

make_params_true <- function(n, K, M) {
  lapply(seq_len(n), function(i) {
    states      <- vector("list", K + M)
    states[[1]] <- list(sigma = SIGMA_GLOB$bear[i],
                        alpha = ALPHA_GLOB$bear[i],
                        nu    = NU_GLOB$bear[i])
    states[[2]] <- list(sigma = SIGMA_GLOB$stable[i],
                        alpha = ALPHA_GLOB$stable[i],
                        nu    = NU_GLOB$stable[i])
    states[[3]] <- list(sigma = SIGMA_GLOB$bull[i],
                        alpha = ALPHA_GLOB$bull[i],
                        nu    = NU_GLOB$bull[i])
    for (m in seq_len(M))
      states[[K + m]] <- list(sigma = 0.8, alpha = 0.0, nu = 10)
    states
  })
}

simulate_data <- function(TT, n, K, M, A_true, rho_val, params_true) {
  rho_mat <- matrix(rho_val, n, K)
  z       <- integer(TT); z[1] <- 2
  for (t in 2:TT)
    z[t] <- sample(K, 1, prob = P_TRUE[z[t - 1], ])
  
  Y <- matrix(0, TT, n)
  for (t in 2:TT) {
    k   <- z[t]
    yt1 <- Y[t - 1, ]
    for (i in seq_len(n)) {
      if (rbinom(1, 1, rho_mat[i, k]) == 1) {
        mu_v <- MU_TRUE[[k]][i]
        ar_v <- A_true[[k]][i, ]
        p_e  <- params_true[[i]][[k]]
      } else {
        s_idx <- K + sample(M, 1)
        mu_v  <- MU_TRUE[[s_idx]][i]
        ar_v  <- A_true[[s_idx]][i, ]
        p_e   <- params_true[[i]][[s_idx]]
      }
      Y[t, i] <- mu_v + sum(ar_v * yt1) +
        rst(1, xi = 0, omega = p_e$sigma,
            alpha = p_e$alpha, nu = p_e$nu)
    }
  }
  list(Y = Y, z = z)
}

kmeans_init <- function(Y, K, M, n) {
  km            <- kmeans(Y, centers = K, nstart = 30, iter.max = 100)
  cluster_order <- order(km$centers[, 1])
  
  init_mu <- lapply(seq_len(K + M), function(s) {
    if (s <= K) as.numeric(km$centers[cluster_order[s], ])
    else        rep(0, n)
  })
  init_A <- lapply(seq_len(K + M), function(s) {
    if (s > K) return(diag(0.02, n))
    idx <- which(km$cluster == cluster_order[s])
    idx <- idx[idx > 1]
    if (length(idx) < n + 2) return(diag(0.05, n))
    Xc <- cbind(1, Y[idx - 1, , drop = FALSE])
    Yc <- Y[idx, , drop = FALSE]
    B  <- tryCatch(
      solve(t(Xc) %*% Xc + diag(1e-6, ncol(Xc)), t(Xc) %*% Yc),
      error = function(e) matrix(0, ncol(Xc), n)
    )
    t(B[2:(n + 1), , drop = FALSE])
  })
  init_params <- lapply(seq_len(n), function(i) {
    states <- vector("list", K + M)
    for (s in seq_len(K)) {
      idx_s <- which(km$cluster == cluster_order[s])
      sig_s <- if (length(idx_s) > 1) sd(Y[idx_s, i]) else 1.0
      states[[s]] <- list(sigma = sig_s, alpha = 0.0, nu = 5)
    }
    for (m in seq_len(M))
      states[[K + m]] <- list(sigma = 0.8, alpha = 0.0, nu = 10)
    states
  })
  list(
    P_trans     = matrix(c(0.80, 0.20, 0.00,
                           0.10, 0.80, 0.10,
                           0.00, 0.20, 0.80), 3, 3, byrow = TRUE),
    pi0         = rep(1 / K, K),
    mu          = init_mu,
    A           = init_A,
    rho         = matrix(0.85, n, K),
    params_emit = init_params
  )
}

align_labels <- function(z_hat, z, K) {
  perms    <- combinat::permn(seq_len(K))
  best_acc <- 0
  best_map <- seq_len(K)
  for (perm in perms) {
    acc <- mean(perm[z_hat] == z)
    if (acc > best_acc) { best_acc <- acc; best_map <- perm }
  }
  list(acc = best_acc, map = best_map)
}

# -----------------------------------------------------------------------------
# FUNÇÃO DE UMA REPLICAÇÃO
# -----------------------------------------------------------------------------

one_rep <- function(rep_id, TT, signal, rho_val, K, M, n, seed_base) {
  set.seed(seed_base + rep_id * 1000)
  
  A_true      <- make_A_true(signal)
  params_true <- make_params_true(n, K, M)
  
  sim <- simulate_data(TT, n, K, M, A_true, rho_val, params_true)
  Y   <- sim$Y
  z   <- sim$z
  
  init <- tryCatch(kmeans_init(Y, K, M, n), error = function(e) NULL)
  
  fit <- tryCatch(
    em_msvar(Y, K = K, M = M, max_iter = 150, tol = 1e-6,
             P_zeros  = list(c(1, 3), c(3, 1)),
             init     = init, verbose = FALSE),
    error = function(e) NULL
  )
  
  if (is.null(fit))
    return(list(converged = FALSE, rep_id = rep_id,
                TT = TT, signal = signal, rho_val = rho_val))
  
  n_iter <- length(fit$history)
  z_hat  <- apply(fit$xi_filt, 1, which.max)
  al     <- align_labels(z_hat, z, K)
  acc    <- al$acc
  map    <- al$map
  
  P_est_al <- matrix(0, K, K)
  for (j in seq_len(K)) for (k in seq_len(K))
    P_est_al[map[j], map[k]] <- fit$P_trans[j, k]
  frob_P <- norm(abs(P_TRUE - P_est_al), type = "F")
  
  frob_A <- rmse_A <- numeric(K)
  for (k in seq_len(K)) {
    k_est     <- which(map == k)
    dev       <- fit$A[[k_est]] - A_true[[k]]
    frob_A[k] <- norm(dev, type = "F")
    rmse_A[k] <- sqrt(mean(dev^2))
  }
  
  rho_true <- matrix(rho_val, n, K)
  rho_al   <- fit$rho[, map]
  rmse_rho <- sqrt(mean((rho_al - rho_true)^2))
  
  W_list   <- influence_matrix(fit$A[seq_len(K)], alpha_sig = 0.10, K = K)
  spearman <- prob_top1 <- numeric(K)
  for (k in seq_len(K)) {
    k_est   <- which(map == k)
    os_true <- sapply(seq_len(n), function(i) sum(abs(A_true[[k]][-i, i])))
    os_est  <- rowSums(W_list[[k_est]])
    sp <- suppressWarnings(cor(os_true, os_est, method = "spearman"))
    spearman[k]  <- ifelse(is.na(sp), NA, sp)
    prob_top1[k] <- as.numeric(which.min(rank(-os_true)) ==
                                 which.min(rank(-os_est)))
  }
  
  f1_wald <- numeric(K)
  for (k in seq_len(K)) {
    k_est <- which(map == k)
    wt    <- tryCatch(
      wald_test_leadership(Y, fit$xi_smooth, k = k, lags_nw = 5),
      error = function(e) NULL
    )
    if (is.null(wt)) { f1_wald[k] <- NA; next }
    pv  <- wt$p_values
    A_t <- A_true[[k_est]]
    tp <- fp <- tn <- fn <- 0
    for (i in seq_len(n)) for (l in seq_len(n)) {
      if (i == l || is.na(pv[i, l])) next
      tnz <- abs(A_t[l, i]) > 1e-6
      rej <- pv[i, l] < 0.10
      if ( tnz &&  rej) tp <- tp + 1
      if (!tnz &&  rej) fp <- fp + 1
      if (!tnz && !rej) tn <- tn + 1
      if ( tnz && !rej) fn <- fn + 1
    }
    prec       <- ifelse(tp + fp > 0, tp / (tp + fp), 0)
    rec        <- ifelse(tp + fn > 0, tp / (tp + fn), 0)
    f1_wald[k] <- ifelse(prec + rec > 0,
                         2 * prec * rec / (prec + rec), 0)
  }
  
  bic_k3 <- compute_bic(fit$log_lik, K = 3, n = n, p = 1, TT = TT)
  fit_k2 <- tryCatch(
    em_msvar(Y, K = 2, M = M, max_iter = 80,
             P_zeros = NULL, verbose = FALSE),
    error = function(e) NULL
  )
  bic_k2      <- if (!is.null(fit_k2))
    compute_bic(fit_k2$log_lik, K = 2, n = n, p = 1, TT = TT) else NA
  bic_correct <- if (!is.na(bic_k2)) as.numeric(bic_k3 < bic_k2) else NA
  
  acc_per_regime <- sapply(seq_len(K), function(k) {
    idx <- which(z == k)
    if (length(idx) == 0) return(NA)
    mean((map[z_hat])[idx] == k)
  })
  
  list(
    converged      = TRUE,
    rep_id         = rep_id,
    TT             = TT,
    signal         = signal,
    rho_val        = rho_val,
    acc            = acc,
    acc_per_regime = acc_per_regime,
    frob_P         = frob_P,
    bias_P         = as.vector(P_est_al - P_TRUE),
    frob_A         = frob_A,
    rmse_A         = rmse_A,
    rmse_rho       = rmse_rho,
    spearman       = spearman,
    prob_top1      = prob_top1,
    f1_wald        = f1_wald,
    bic_k3         = bic_k3,
    bic_k2         = bic_k2,
    bic_correct    = bic_correct,
    n_iter         = n_iter,
    label_switch   = as.numeric(!all(map == seq_len(K)))
  )
}

# -----------------------------------------------------------------------------
# BACKEND PARALELO — compatível com Windows, Linux e macOS
# -----------------------------------------------------------------------------

make_cluster <- function(n_cores) {
  os_type  <- .Platform$OS.type
  max_safe <- if (os_type == "windows") min(n_cores, 10L) else n_cores
  if (max_safe < n_cores)
    message(sprintf(
      "Windows: limitando workers de %d para %d (limite seguro PSOCK)",
      n_cores, max_safe))
  
  utils_path <- normalizePath("utils.R", mustWork = TRUE)
  main_wd    <- getwd()
  
  # Valida ambiente com 1 worker antes de abrir o cluster completo
  cl_test <- tryCatch({
    ct <- makeCluster(1L, type = "PSOCK")
    clusterCall(ct, function(wd) setwd(wd), main_wd)
    clusterEvalQ(ct, {
      library(sn); library(igraph)
      library(Matrix); library(combinat)
    })
    clusterExport(ct, "utils_path", envir = environment())
    clusterEvalQ(ct, source(utils_path))
    stopCluster(ct)
    TRUE
  }, error = function(e) {
    message("Falha no worker de teste: ", e$message)
    FALSE
  })
  
  if (!cl_test) stop("Nao foi possivel inicializar workers. Verifique utils.R.")
  
  cl <- makeCluster(max_safe, type = "PSOCK",
                    outfile = file.path(main_wd, "mc_worker.log"))
  
  setup_ok <- tryCatch({
    clusterCall(cl, function(wd) setwd(wd), main_wd)
    clusterEvalQ(cl, {
      library(sn); library(igraph)
      library(Matrix); library(combinat)
    })
    clusterExport(cl, "utils_path", envir = environment())
    clusterEvalQ(cl, source(utils_path))
    TRUE
  }, error = function(e) {
    message("Erro ao configurar workers: ", e$message)
    tryCatch(stopCluster(cl), error = function(e) NULL)
    FALSE
  })
  
  if (!setup_ok) stop("Falha na configuracao dos workers.")
  
  clusterExport(cl, varlist = c(
    "P_TRUE", "MU_TRUE", "AR_OFFDIAG",
    "SIGMA_GLOB", "ALPHA_GLOB", "NU_GLOB",
    "MC_CONFIG",
    "make_A_true", "make_params_true",
    "simulate_data", "kmeans_init",
    "align_labels", "one_rep"
  ), envir = globalenv())
  
  message(sprintf("Cluster iniciado com %d workers.", max_safe))
  cl
}

par_run <- function(cl, batch, TT, signal, rho_val, cfg) {
  tryCatch({
    clusterExport(cl, varlist = c(
      "P_TRUE", "MU_TRUE", "AR_OFFDIAG",
      "SIGMA_GLOB", "ALPHA_GLOB", "NU_GLOB"
    ), envir = globalenv())
    parLapplyLB(cl, batch, function(r) {
      one_rep(r, TT, signal, rho_val,
              cfg$K, cfg$M, cfg$n, cfg$seed)
    })
  }, error = function(e) {
    warning("Paralelo falhou (", e$message, ") — executando serial.")
    lapply(batch, function(r)
      one_rep(r, TT, signal, rho_val,
              cfg$K, cfg$M, cfg$n, cfg$seed))
  })
}

# -----------------------------------------------------------------------------
# GRADE DE CÉLULAS
# -----------------------------------------------------------------------------

cells <- expand.grid(
  TT      = MC_CONFIG$T_grid,
  signal  = MC_CONFIG$sig_grid,
  rho_val = MC_CONFIG$rho_grid,
  stringsAsFactors = FALSE
)

cat(sprintf("Grade: %d células x %d replicações = %d execuções totais\n",
            nrow(cells), MC_CONFIG$R, nrow(cells) * MC_CONFIG$R))
cat(sprintf("Núcleos utilizados: %d\n\n", MC_CONFIG$n_cores))

# -----------------------------------------------------------------------------
# LOOP PRINCIPAL COM CHECKPOINT
# -----------------------------------------------------------------------------

cl <- make_cluster(MC_CONFIG$n_cores)
on.exit(tryCatch(stopCluster(cl), error = function(e) NULL), add = TRUE)

for (cell_idx in seq_len(nrow(cells))) {
  
  TT      <- cells$TT[cell_idx]
  signal  <- cells$signal[cell_idx]
  rho_val <- cells$rho_val[cell_idx]
  
  cell_tag  <- sprintf("T%d_sig%s_rho%.2f", TT, signal, rho_val)
  ckpt_file <- file.path(MC_CONFIG$out_dir,
                         paste0("ckpt_", cell_tag, ".rds"))
  
  if (file.exists(ckpt_file)) {
    ckpt     <- readRDS(ckpt_file)
    done_ids <- ckpt$done_ids
    results  <- ckpt$results
    cat(sprintf("[%d/%d] %s — retomando (%d/%d)\n",
                cell_idx, nrow(cells), cell_tag,
                length(done_ids), MC_CONFIG$R))
  } else {
    done_ids <- integer(0)
    results  <- vector("list", MC_CONFIG$R)
    cat(sprintf("[%d/%d] %s — iniciando\n",
                cell_idx, nrow(cells), cell_tag))
  }
  
  pending <- setdiff(seq_len(MC_CONFIG$R), done_ids)
  if (length(pending) == 0) { cat("  Já completa.\n"); next }
  
  batches <- split(pending,
                   ceiling(seq_along(pending) / 50))
  
  for (batch in batches) {
    batch_res <- par_run(cl, batch, TT, signal, rho_val, MC_CONFIG)
    
    for (i in seq_along(batch))
      results[[batch[i]]] <- batch_res[[i]]
    
    done_ids <- c(done_ids, batch)
    saveRDS(list(done_ids = done_ids, results = results), ckpt_file)
    cat(sprintf("  %d/%d replicações concluídas\n",
                length(done_ids), MC_CONFIG$R))
  }
  
  cat(sprintf("  Célula %s completa.\n\n", cell_tag))
}

tryCatch(stopCluster(cl), error = function(e) NULL)
cat("Todas as células concluídas.\n")

# -----------------------------------------------------------------------------
# COMPILAÇÃO DOS RESULTADOS
# -----------------------------------------------------------------------------

cat("Compilando resultados...\n")

compile_results <- function(out_dir, cells, R) {
  rows <- list()
  for (cell_idx in seq_len(nrow(cells))) {
    TT      <- cells$TT[cell_idx]
    signal  <- cells$signal[cell_idx]
    rho_val <- cells$rho_val[cell_idx]
    cell_tag  <- sprintf("T%d_sig%s_rho%.2f", TT, signal, rho_val)
    ckpt_file <- file.path(out_dir, paste0("ckpt_", cell_tag, ".rds"))
    if (!file.exists(ckpt_file)) next
    res <- readRDS(ckpt_file)$results
    
    for (r in seq_len(R)) {
      rv <- res[[r]]
      if (is.null(rv) || !isTRUE(rv$converged)) {
        rows[[length(rows) + 1]] <- data.frame(
          TT = TT, signal = signal, rho_val = rho_val, rep_id = r,
          converged     = FALSE,
          acc           = NA, frob_P      = NA,
          frob_A_bear   = NA, frob_A_stable = NA, frob_A_bull = NA,
          rmse_A_bear   = NA, rmse_A_stable = NA, rmse_A_bull = NA,
          rmse_rho      = NA,
          sp_bear       = NA, sp_stable   = NA, sp_bull     = NA,
          top1_bear     = NA, top1_stable = NA, top1_bull   = NA,
          f1_bear       = NA, f1_stable   = NA, f1_bull     = NA,
          bic_correct   = NA, n_iter      = NA, label_switch = NA,
          acc_bear      = NA, acc_stable  = NA, acc_bull    = NA
        )
        next
      }
      rows[[length(rows) + 1]] <- data.frame(
        TT            = TT,
        signal        = signal,
        rho_val       = rho_val,
        rep_id        = r,
        converged     = TRUE,
        acc           = rv$acc,
        frob_P        = rv$frob_P,
        frob_A_bear   = rv$frob_A[1],
        frob_A_stable = rv$frob_A[2],
        frob_A_bull   = rv$frob_A[3],
        rmse_A_bear   = rv$rmse_A[1],
        rmse_A_stable = rv$rmse_A[2],
        rmse_A_bull   = rv$rmse_A[3],
        rmse_rho      = rv$rmse_rho,
        sp_bear       = rv$spearman[1],
        sp_stable     = rv$spearman[2],
        sp_bull       = rv$spearman[3],
        top1_bear     = rv$prob_top1[1],
        top1_stable   = rv$prob_top1[2],
        top1_bull     = rv$prob_top1[3],
        f1_bear       = rv$f1_wald[1],
        f1_stable     = rv$f1_wald[2],
        f1_bull       = rv$f1_wald[3],
        bic_correct   = rv$bic_correct,
        n_iter        = rv$n_iter,
        label_switch  = rv$label_switch,
        acc_bear      = rv$acc_per_regime[1],
        acc_stable    = rv$acc_per_regime[2],
        acc_bull      = rv$acc_per_regime[3]
      )
    }
  }
  do.call(rbind, rows)
}

df_all <- compile_results(MC_CONFIG$out_dir, cells, MC_CONFIG$R)
saveRDS(df_all, file.path(MC_CONFIG$out_dir, "mc_all_results.rds"))
write.csv(df_all, file.path(MC_CONFIG$out_dir, "mc_all_results.csv"),
          row.names = FALSE)

cat(sprintf("Total de linhas: %d\n", nrow(df_all)))
cat(sprintf("Taxa de convergência: %.1f%%\n",
            100 * mean(df_all$converged, na.rm = TRUE)))

# -----------------------------------------------------------------------------
# TABELA-RESUMO
# -----------------------------------------------------------------------------

summarise_mc <- function(df) {
  df_conv <- df[isTRUE(df$converged) | df$converged == TRUE, ]
  metrics <- c("acc", "frob_P",
               "rmse_A_bear", "rmse_A_bull",
               "rmse_rho",
               "sp_bear", "sp_bull",
               "top1_bear", "top1_bull",
               "f1_bear", "f1_bull",
               "bic_correct", "label_switch", "n_iter",
               "acc_bear", "acc_stable", "acc_bull")
  
  result <- do.call(rbind, lapply(
    split(df_conv,
          list(df_conv$TT, df_conv$signal, df_conv$rho_val)),
    function(sub) {
      if (nrow(sub) == 0) return(NULL)
      row <- data.frame(TT      = sub$TT[1],
                        signal  = sub$signal[1],
                        rho_val = sub$rho_val[1],
                        n_conv  = nrow(sub))
      for (m in metrics) {
        x <- sub[[m]]
        row[[paste0(m, "_mean")]] <- mean(x,                  na.rm = TRUE)
        row[[paste0(m, "_sd")]]   <- sd(x,                    na.rm = TRUE)
        row[[paste0(m, "_p10")]]  <- quantile(x, 0.10,        na.rm = TRUE)
        row[[paste0(m, "_p90")]]  <- quantile(x, 0.90,        na.rm = TRUE)
      }
      row
    }
  ))
  result[order(result$TT, result$signal, result$rho_val), ]
}

df_summary <- summarise_mc(df_all)
saveRDS(df_summary, file.path(MC_CONFIG$out_dir, "mc_summary.rds"))
write.csv(df_summary, file.path(MC_CONFIG$out_dir, "mc_summary.csv"),
          row.names = FALSE)

cat("\n=== TABELA-RESUMO PRINCIPAL ===\n")
cols_print <- c("TT", "signal", "rho_val", "n_conv",
                "acc_mean", "frob_P_mean",
                "rmse_A_bear_mean", "rmse_A_bull_mean",
                "sp_bear_mean", "sp_bull_mean",
                "top1_bear_mean", "top1_bull_mean",
                "f1_bear_mean", "f1_bull_mean",
                "bic_correct_mean", "rmse_rho_mean",
                "label_switch_mean")
print(round(df_summary[, cols_print], 3), row.names = FALSE)

# -----------------------------------------------------------------------------
# FIGURAS
# -----------------------------------------------------------------------------

theme_mc <- theme_bw(base_size = 11) +
  theme(legend.position  = "bottom",
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey92"))

sig_levels  <- c("fraco", "moderado", "forte")
sig_colours <- c("fraco"    = "#E07B54",
                 "moderado" = "#4A90D9",
                 "forte"    = "#27AE60")
df_summary$signal <- factor(df_summary$signal, levels = sig_levels)

# Figura 1 — Acurácia de classificação
p1 <- ggplot(df_summary,
             aes(x = TT, y = acc_mean, colour = signal,
                 ymin = acc_mean - acc_sd,
                 ymax = acc_mean + acc_sd)) +
  geom_ribbon(aes(fill = signal), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.9) + geom_point(size = 2) +
  facet_wrap(~rho_val, labeller = label_both) +
  scale_colour_manual(values = sig_colours, name = "Sinal") +
  scale_fill_manual(values   = sig_colours, name = "Sinal") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     limits = c(0.5, 1)) +
  labs(title = "Acurácia de Classificação de Regime vs T",
       x = "T", y = "Acurácia (%)") +
  theme_mc
ggsave(file.path(MC_CONFIG$out_dir, "fig1_accuracy.png"),
       p1, width = 10, height = 7, dpi = 150)

# Figura 2 — Spearman Bear e Bull
df_sp <- rbind(
  data.frame(df_summary[, c("TT","signal","rho_val")],
             spearman = df_summary$sp_bear_mean,
             sd       = df_summary$sp_bear_sd, regime = "Bear"),
  data.frame(df_summary[, c("TT","signal","rho_val")],
             spearman = df_summary$sp_bull_mean,
             sd       = df_summary$sp_bull_sd, regime = "Bull")
)
df_sp$signal <- factor(df_sp$signal, levels = sig_levels)

p2 <- ggplot(df_sp,
             aes(x = TT, y = spearman, colour = signal,
                 ymin = spearman - sd, ymax = spearman + sd)) +
  geom_ribbon(aes(fill = signal), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.9) + geom_point(size = 2) +
  geom_hline(yintercept = 0.8, linetype = "dashed", colour = "grey40") +
  facet_grid(regime ~ rho_val, labeller = label_both) +
  scale_colour_manual(values = sig_colours, name = "Sinal") +
  scale_fill_manual(values   = sig_colours, name = "Sinal") +
  labs(title    = "Correlação de Spearman do Ranking de Liderança vs T",
       subtitle = "Linha tracejada = limiar 0.8",
       x = "T", y = "Spearman rho") +
  theme_mc
ggsave(file.path(MC_CONFIG$out_dir, "fig2_spearman.png"),
       p2, width = 10, height = 7, dpi = 150)

# Figura 3 — P(líder correto)
df_top <- rbind(
  data.frame(df_summary[, c("TT","signal","rho_val")],
             top1 = df_summary$top1_bear_mean, regime = "Bear"),
  data.frame(df_summary[, c("TT","signal","rho_val")],
             top1 = df_summary$top1_bull_mean, regime = "Bull")
)
df_top$signal <- factor(df_top$signal, levels = sig_levels)

p3 <- ggplot(df_top,
             aes(x = TT, y = top1, colour = signal)) +
  geom_line(linewidth = 0.9) + geom_point(size = 2) +
  geom_hline(yintercept = 0.9, linetype = "dashed", colour = "grey40") +
  facet_grid(regime ~ rho_val, labeller = label_both) +
  scale_colour_manual(values = sig_colours, name = "Sinal") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     limits = c(0, 1)) +
  labs(title    = "P(Identificação Correta do Líder) vs T",
       subtitle = "Linha tracejada = 90%",
       x = "T", y = "P(Rank 1 correto)") +
  theme_mc
ggsave(file.path(MC_CONFIG$out_dir, "fig3_top1.png"),
       p3, width = 10, height = 7, dpi = 150)

# Figura 4 — F1 do Wald
df_f1 <- rbind(
  data.frame(df_summary[, c("TT","signal","rho_val")],
             f1 = df_summary$f1_bear_mean,
             sd = df_summary$f1_bear_sd, regime = "Bear"),
  data.frame(df_summary[, c("TT","signal","rho_val")],
             f1 = df_summary$f1_bull_mean,
             sd = df_summary$f1_bull_sd, regime = "Bull")
)
df_f1$signal <- factor(df_f1$signal, levels = sig_levels)

p4 <- ggplot(df_f1,
             aes(x = TT, y = f1, colour = signal,
                 ymin = f1 - sd, ymax = f1 + sd)) +
  geom_ribbon(aes(fill = signal), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.9) + geom_point(size = 2) +
  facet_grid(regime ~ rho_val, labeller = label_both) +
  scale_colour_manual(values = sig_colours, name = "Sinal") +
  scale_fill_manual(values   = sig_colours, name = "Sinal") +
  scale_y_continuous(limits = c(0, 1)) +
  labs(title = "F1 do Teste de Wald (Liderança) vs T",
       x = "T", y = "F1") +
  theme_mc
ggsave(file.path(MC_CONFIG$out_dir, "fig4_f1_wald.png"),
       p4, width = 10, height = 7, dpi = 150)

# Figura 5 — Taxa BIC correto
p5 <- ggplot(df_summary,
             aes(x = TT, y = bic_correct_mean, colour = signal)) +
  geom_line(linewidth = 0.9) + geom_point(size = 2) +
  geom_hline(yintercept = 0.9, linetype = "dashed", colour = "grey40") +
  facet_wrap(~rho_val, labeller = label_both) +
  scale_colour_manual(values = sig_colours, name = "Sinal") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     limits = c(0, 1)) +
  labs(title    = "Taxa de Seleção Correta K=3 via BIC vs T",
       subtitle = "Linha tracejada = 90%",
       x = "T", y = "P(BIC seleciona K=3)") +
  theme_mc
ggsave(file.path(MC_CONFIG$out_dir, "fig5_bic.png"),
       p5, width = 10, height = 7, dpi = 150)

# Figura 6 — RMSE de rho
p6 <- ggplot(df_summary,
             aes(x = TT, y = rmse_rho_mean, colour = signal)) +
  geom_line(linewidth = 0.9) + geom_point(size = 2) +
  facet_wrap(~rho_val, labeller = label_both) +
  scale_colour_manual(values = sig_colours, name = "Sinal") +
  labs(title    = "RMSE de rho (Participação) vs T",
       subtitle = "Avalia o viés de atenuação do M-step",
       x = "T", y = "RMSE(rho)") +
  theme_mc
ggsave(file.path(MC_CONFIG$out_dir, "fig6_rmse_rho.png"),
       p6, width = 10, height = 7, dpi = 150)

# Figura 7 — Heatmap Spearman Bear (rho=0.90)
df_heat <- df_summary[df_summary$rho_val == 0.90, ]
df_heat$signal <- factor(df_heat$signal, levels = sig_levels)

p7 <- ggplot(df_heat,
             aes(x = factor(TT), y = signal, fill = sp_bear_mean)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.2f", sp_bear_mean)), size = 3.5) +
  scale_fill_gradient2(low = "#E74C3C", mid = "#F7DC6F",
                       high = "#27AE60", midpoint = 0.7,
                       limits = c(0, 1), name = "Spearman") +
  labs(title = "Heatmap Spearman — Regime Bear (rho=0.90)",
       x = "T", y = "Força do Sinal") +
  theme_mc
ggsave(file.path(MC_CONFIG$out_dir, "fig7_heatmap_spearman.png"),
       p7, width = 8, height = 3.5, dpi = 150)

cat("\nFiguras salvas em:", MC_CONFIG$out_dir, "\n")
cat("\nArquivos gerados:\n")
cat("  mc_all_results.rds/.csv   — dados brutos de todas as replicações\n")
cat("  mc_summary.rds/.csv       — tabela resumo agregada\n")
cat("  fig1_accuracy.png         — acurácia de classificação\n")
cat("  fig2_spearman.png         — Spearman do ranking\n")
cat("  fig3_top1.png             — P(líder correto)\n")
cat("  fig4_f1_wald.png          — F1 do Wald\n")
cat("  fig5_bic.png              — seleção K via BIC\n")
cat("  fig6_rmse_rho.png         — RMSE de rho\n")
cat("  fig7_heatmap_spearman.png — heatmap Spearman\n")