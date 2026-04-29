# =============================================================================
# MGHMM - Mixture Gaussian Hidden Markov Model
# Dias, Vermunt & Ramos (2015) - European Journal of Operational Research
#
# Aplicação: 30 ações do Dow Jones Industrial Average (DJIA)
# TCC – Identificação de Regimes de Mercado via MGHMM
# =============================================================================
#
# ORDEM DE EXECUÇÃO:
#   1. Rode download_dados_djia.R  (apenas uma vez, gera dados_djia.RData)
#   2. Rode este script
#
# PACOTES NECESSÁRIOS:
#   install.packages(c("quantmod", "zoo", "ggplot2", "reshape2", "gridExtra",
#                      "tseries", "moments"))
# =============================================================================

library(quantmod)
library(zoo)
library(ggplot2)
library(reshape2)
library(gridExtra)




# =============================================================================
# 1. FUNÇÕES AUXILIARES
# =============================================================================


gaussian_density <- function(y, mu, sigma2) {
  sigma2 <- pmax(sigma2, 1e-8)
  dnorm(y, mean = mu, sd = sqrt(sigma2))
}


transition_probs <- function(y_prev, gamma0, gamma1, K) {
  P <- matrix(0, K, K)
  for (j in 1:K) {
    eta    <- gamma0[j, ] + gamma1[j, ] * y_prev
    eta[K] <- 0
    eta    <- eta - max(eta)
    P[j, ] <- exp(eta) / sum(exp(eta))
  }
  return(P)
}




# =============================================================================
# 2. FORWARD-BACKWARD ESTENDIDO
# =============================================================================


forward_backward_mghmm <- function(y, S, K, pi_w, delta_w,
                                   gamma0, gamma1, mu, sigma2) {
  T       <- length(y)
  alpha   <- array(0,  dim = c(T, S, K))
  beta    <- array(1,  dim = c(T, S, K))
  c_scale <- matrix(0, T, S)

  # Forward t = 1
  for (s in 1:S) {
    alpha[1, s, ] <- delta_w[s, ] * gaussian_density(y[1], mu, sigma2)
    cs <- sum(alpha[1, s, ])
    if (!is.finite(cs) || cs <= 0) {
      c_scale[1, s] <- 1e-300; alpha[1, s, ] <- rep(1/K, K)
    } else {
      c_scale[1, s] <- cs;     alpha[1, s, ] <- alpha[1, s, ] / cs
    }
  }

  # Forward t = 2,...,T
  for (t in 2:T) {
    for (s in 1:S) {
      P_t <- transition_probs(y[t-1], gamma0[s,,], gamma1[s,,], K)
      for (k in 1:K) {
        alpha[t, s, k] <- gaussian_density(y[t], mu[k], sigma2[k]) *
          sum(alpha[t-1, s, ] * P_t[, k])
      }
      cs <- sum(alpha[t, s, ])
      if (!is.finite(cs) || cs <= 0) {
        c_scale[t, s] <- 1e-300; alpha[t, s, ] <- rep(1/K, K)
      } else {
        c_scale[t, s] <- cs;     alpha[t, s, ] <- alpha[t, s, ] / cs
      }
    }
  }

  # Backward t = T
  beta[T, , ] <- 1

  # Backward t = T-1,...,1
  for (t in (T-1):1) {
    for (s in 1:S) {
      P_t1 <- transition_probs(y[t], gamma0[s,,], gamma1[s,,], K)
      for (k in 1:K) {
        beta[t, s, k] <- sum(
          P_t1[k, ] * gaussian_density(y[t+1], mu, sigma2) * beta[t+1, s, ]
        )
      }
      bs <- sum(beta[t, s, ])
      if (!is.finite(bs) || bs <= 0) {
        beta[t, s, ] <- rep(1/K, K)
      } else {
        beta[t, s, ] <- beta[t, s, ] / bs
      }
    }
  }

  # Verossimilhança marginal por classe
  iw <- numeric(S)
  for (s in 1:S) {
    log_scale_sum <- sum(log(pmax(c_scale[, s], 1e-300)))
    ab_sum        <- sum(alpha[T, s, ] * beta[T, s, ])
    val           <- pi_w[s] * exp(log_scale_sum) * ab_sum
    iw[s]         <- ifelse(is.finite(val) && val > 0, val, 0)
  }
  fy <- sum(iw)
  if (!is.finite(fy) || fy <= 0) fy <- 1e-300

  # Posteriors
  post_w <- iw / fy
  post_w[!is.finite(post_w)] <- 1/S

  post_wzt <- array(0, dim = c(T, S, K))
  for (t in 1:T) {
    for (s in 1:S) {
      ab    <- alpha[t, s, ] * beta[t, s, ]
      denom <- sum(ab)
      if (!is.finite(denom) || denom <= 0) {
        post_wzt[t, s, ] <- (iw[s] / fy) * rep(1/K, K)
      } else {
        post_wzt[t, s, ] <- (iw[s] / fy) * ab / denom
      }
    }
  }

  post_wzt1zt <- array(0, dim = c(T-1, S, K, K))
  for (t in 1:(T-1)) {
    for (s in 1:S) {
      P_t1  <- transition_probs(y[t], gamma0[s,,], gamma1[s,,], K)
      denom <- sum(alpha[t, s, ] * beta[t, s, ])
      if (!is.finite(denom) || denom <= 0) next
      for (j in 1:K) {
        for (k in 1:K) {
          val <- (iw[s] / fy) * alpha[t, s, j] * P_t1[j, k] *
            gaussian_density(y[t+1], mu[k], sigma2[k]) *
            beta[t+1, s, k] / denom
          post_wzt1zt[t, s, j, k] <- ifelse(is.finite(val), val, 0)
        }
      }
    }
  }

  log_lik <- log(fy) + sum(log(pmax(c_scale, 1e-300)))

  return(list(
    alpha = alpha, beta = beta,
    post_w = post_w, post_wzt = post_wzt,
    post_wzt1zt = post_wzt1zt,
    log_lik = log_lik
  ))
}




# =============================================================================
# 3. M-STEP
# =============================================================================


mstep_mghmm <- function(y_list, S, K, fb_list) {
  n <- length(y_list)
  T <- length(y_list[[1]])

  # pi_w
  post_w_mat <- do.call(rbind, lapply(fb_list, function(fb) fb$post_w))
  pi_w_new   <- colMeans(post_w_mat)
  pi_w_new   <- pmax(pi_w_new, 1e-10)
  pi_w_new   <- pi_w_new / sum(pi_w_new)

  # delta_w
  delta_w_new <- matrix(0, S, K)
  for (s in 1:S) {
    d <- colMeans(do.call(rbind,
                          lapply(fb_list, function(fb) fb$post_wzt[1, s, ])))
    d <- pmax(d, 1e-10)
    delta_w_new[s, ] <- d / sum(d)
  }

  # mu e sigma2
  mu_new     <- numeric(K)
  sigma2_new <- numeric(K)
  all_y      <- unlist(y_list)

  for (k in 1:K) {
    w_k <- numeric(n * T)
    idx <- 1
    for (i in 1:n) {
      w_k[idx:(idx+T-1)] <- rowSums(
        matrix(fb_list[[i]]$post_wzt[, , k], T, S))
      idx <- idx + T
    }
    w_k <- pmax(w_k, 0)
    sw  <- sum(w_k)
    if (sw < 1e-10) {
      mu_new[k] <- 0; sigma2_new[k] <- 1
    } else {
      mu_new[k]     <- sum(w_k * all_y) / sw
      sigma2_new[k] <- max(sum(w_k * (all_y - mu_new[k])^2) / sw, 1e-8)
    }
  }

  # gamma0 e gamma1 via optim
  gamma0_new <- array(0, dim = c(S, K, K))
  gamma1_new <- array(0, dim = c(S, K, K))

  if (K > 1) {
    for (s in 1:S) {
      for (j in 1:K) {
        w_agg  <- matrix(0, T-1, K)
        yp_agg <- y_list[[1]][1:(T-1)]
        for (i in 1:n) {
          w_agg <- w_agg +
            matrix(fb_list[[i]]$post_wzt1zt[, s, j, ], T-1, K)
        }
        w_j <- rowSums(w_agg)
        if (sum(w_j) < 1e-10) next
        p_obs <- w_agg / pmax(w_j, 1e-15)

        obj_fn <- function(params) {
          g0  <- c(params[1:(K-1)], 0)
          g1  <- c(params[K:(2*(K-1))], 0)
          ll  <- 0
          for (t in seq_len(T-1)) {
            if (w_j[t] < 1e-15) next
            eta <- g0 + g1 * yp_agg[t]
            eta <- eta - max(eta)
            lp  <- eta - log(sum(exp(eta)))
            ll  <- ll - w_j[t] * sum(p_obs[t, ] * lp)
          }
          if (!is.finite(ll)) ll <- 1e10
          return(ll)
        }

        init <- rep(0, 2*(K-1))
        opt  <- tryCatch(
          optim(init, obj_fn, method = "BFGS",
                control = list(maxit = 100, reltol = 1e-6)),
          error = function(e) list(par = init)
        )
        gamma0_new[s, j, 1:(K-1)] <- opt$par[1:(K-1)]
        gamma1_new[s, j, 1:(K-1)] <- opt$par[K:(2*(K-1))]
      }
    }
  }

  return(list(pi_w = pi_w_new, delta_w = delta_w_new,
              gamma0 = gamma0_new, gamma1 = gamma1_new,
              mu = mu_new, sigma2 = sigma2_new))
}




# =============================================================================
# 4. EM COMPLETO
# =============================================================================


fit_mghmm <- function(y_list, S = 2, K = 3,
                      max_iter = 200, tol = 1e-5,
                      n_init = 10, verbose = TRUE) {
  n       <- length(y_list)
  best_ll <- -Inf
  best_params <- NULL

  for (init_run in 1:n_init) {
    if (verbose) cat("Inicialização", init_run, "de", n_init, "\n")

    set.seed(init_run * 42)
    pi_w    <- rep(1/S, S)
    delta_w <- matrix(1/K, S, K)
    mu      <- sort(rnorm(K, 0, 1))
    sigma2  <- runif(K, 0.5, 3)
    gamma0  <- array(rnorm(S*K*K, 0, 0.1),  dim = c(S, K, K))
    gamma1  <- array(rnorm(S*K*K, 0, 0.01), dim = c(S, K, K))
    gamma0[,,K] <- 0; gamma1[,,K] <- 0

    ll_prev <- -Inf
    ll_curr <- -Inf

    for (iter in 1:max_iter) {

      all_fb <- lapply(y_list, function(y) {
        forward_backward_mghmm(y, S, K, pi_w, delta_w,
                               gamma0, gamma1, mu, sigma2)
      })

      ll_curr <- sum(sapply(all_fb, function(fb) fb$log_lik))

      if (!is.finite(ll_curr)) {
        if (verbose) cat("  LL não finito, abortando inicialização\n")
        ll_curr <- -Inf
        break
      }

      if (verbose && iter %% 10 == 0)
        cat("  Iter", iter, "| LL =", round(ll_curr, 4), "\n")

      if (is.finite(ll_prev) && abs(ll_curr - ll_prev) < tol) {
        if (verbose) cat("  Convergiu na iteração", iter, "\n")
        break
      }
      ll_prev <- ll_curr

      params_new <- mstep_mghmm(y_list, S, K, all_fb)
      pi_w    <- params_new$pi_w
      delta_w <- params_new$delta_w
      mu      <- params_new$mu
      sigma2  <- params_new$sigma2
      gamma0  <- params_new$gamma0
      gamma1  <- params_new$gamma1
    }

    if (is.finite(ll_curr) && ll_curr > best_ll) {
      best_ll  <- ll_curr
      N_params <- (S-1) + S*(K-1) + 2*S*K*(K-1) + 2*K
      best_params <- list(
        pi_w     = pi_w,
        delta_w  = delta_w,
        gamma0   = gamma0,
        gamma1   = gamma1,
        mu       = mu,
        sigma2   = sigma2,
        log_lik  = ll_curr,
        n_iter   = iter,
        n_params = N_params,
        BIC      = -2 * ll_curr + N_params * log(n),
        S        = S,
        K        = K,
        all_fb   = all_fb
      )
    }
  }

  if (is.null(best_params)) stop("Nenhuma inicialização convergiu.")
  return(best_params)
}




# =============================================================================
# 5. SELEÇÃO DE MODELO VIA BIC
# =============================================================================


select_model_mghmm <- function(y_list, S_range = 1:2, K_range = 2:3,
                               n_init = 5, verbose = FALSE) {
  results <- data.frame()
  for (S in S_range) {
    for (K in K_range) {
      cat("Ajustando S =", S, "| K =", K, "...\n")
      fit <- tryCatch(
        fit_mghmm(y_list, S = S, K = K,
                  n_init = n_init, verbose = verbose),
        error = function(e) { cat("  Erro:", e$message, "\n"); NULL }
      )
      if (!is.null(fit)) {
        results <- rbind(results, data.frame(
          S        = S,
          K        = K,
          LogLik   = fit$log_lik,
          n_params = fit$n_params,
          BIC      = fit$BIC
        ))
      }
    }
  }
  results[order(results$BIC), ]
}




# =============================================================================
# 6. PREVISÃO 1-PASSO-À-FRENTE
# =============================================================================


predict_mghmm <- function(y, params, class_idx = NULL) {
  S <- params$S
  K <- params$K
  T <- length(y)

  if (is.null(class_idx)) {
    fb        <- forward_backward_mghmm(y, S, K, params$pi_w, params$delta_w,
                                        params$gamma0, params$gamma1,
                                        params$mu, params$sigma2)
    class_idx <- which.max(fb$post_w)
  }

  e        <- params$delta_w[class_idx, ]
  y_hat    <- numeric(T)
  y_hat[1] <- sum(e * params$mu)

  for (t in 2:T) {
    P_t      <- transition_probs(y[t-1], params$gamma0[class_idx,,],
                                 params$gamma1[class_idx,,], K)
    e        <- e %*% P_t
    y_hat[t] <- sum(e * params$mu)
  }
  return(y_hat)
}




# =============================================================================
# 7. SINCRONIZAÇÃO ENTRE MERCADOS
# =============================================================================


compute_synchronization <- function(fb_list, series_names = NULL,
                                    regime_k = 1) {
  n        <- length(fb_list)
  T        <- dim(fb_list[[1]]$post_wzt)[1]
  S        <- dim(fb_list[[1]]$post_wzt)[2]
  K        <- dim(fb_list[[1]]$post_wzt)[3]
  regime_k <- min(regime_k, K)

  build_prob_mat <- function(rk) {
    pm <- matrix(0, T, n)
    for (i in 1:n) {
      p          <- rowSums(matrix(fb_list[[i]]$post_wzt[, , rk], T, S))
      p          <- pmin(pmax(p, 1e-6), 1 - 1e-6)
      pm[, i]    <- log(p / (1 - p))
    }
    pm
  }

  prob_mat <- build_prob_mat(regime_k)
  sd_cols  <- apply(prob_mat, 2, sd)
  ok       <- is.finite(sd_cols) & sd_cols > 1e-8

  if (sum(ok) < 2) {
    avg_occ  <- sapply(1:K, function(k)
      mean(sapply(fb_list, function(fb)
        mean(rowSums(matrix(fb$post_wzt[,,k], T, S))))))
    regime_k <- which.max(avg_occ)
    message("Regime ", regime_k,
            " usado para sincronização (maior ocupação média = ",
            round(avg_occ[regime_k], 3), ")")
    prob_mat <- build_prob_mat(regime_k)
    sd_cols  <- apply(prob_mat, 2, sd)
    ok       <- is.finite(sd_cols) & sd_cols > 1e-8
  }

  sync <- matrix(NA, n, n)
  diag(sync) <- 1
  if (!is.null(series_names)) rownames(sync) <- colnames(sync) <- series_names

  if (sum(ok) >= 2) {
    sync_sub     <- cor(prob_mat[, ok, drop = FALSE], use = "complete.obs")
    sync[ok, ok] <- sync_sub
  } else {
    warning("Nenhum regime apresenta variação suficiente para calcular sincronização.")
  }
  return(sync)
}




# =============================================================================
# 8. VISUALIZAÇÕES
# =============================================================================


plot_regimes <- function(y, fb, series_name = "Série", K = NULL) {
  T  <- length(y)
  S  <- dim(fb$post_wzt)[2]
  if (is.null(K)) K <- dim(fb$post_wzt)[3]

  cores <- c("firebrick","steelblue","darkgreen","orange","purple")[1:K]

  post_df <- as.data.frame(
    t(sapply(1:T, function(t)
      colSums(matrix(fb$post_wzt[t,,], S, K)))))
  colnames(post_df) <- paste0("Regime_", 1:K)
  post_df$t <- 1:T
  post_long <- melt(post_df, id.vars = "t",
                    variable.name = "Regime",
                    value.name    = "Probabilidade")

  modal <- apply(post_df[, 1:K], 1, which.max)

  p1 <- ggplot(data.frame(t = 1:T, retorno = y, regime = factor(modal)),
               aes(x = t, y = retorno, color = regime)) +
    geom_line(linewidth = 0.35) +
    scale_color_manual(values = cores, labels = paste("Regime", 1:K)) +
    labs(title = paste("Retornos e Regimes –", series_name),
         x = "Tempo", y = "Retorno (%)", color = NULL) +
    theme_minimal(base_size = 11)

  p2 <- ggplot(post_long, aes(x = t, y = Probabilidade, fill = Regime)) +
    geom_area(alpha = 0.75, position = "identity") +
    scale_fill_manual(values = cores) +
    labs(title = paste("Probabilidades Posteriores –", series_name),
         x = "Tempo", y = "Probabilidade", fill = NULL) +
    theme_minimal(base_size = 11)

  grid.arrange(p1, p2, ncol = 1)
}


plot_transition_probs <- function(params, class_idx = 1) {
  S      <- params$S
  K      <- params$K
  y_grid <- seq(-10, 10, length.out = 300)
  cores  <- c("firebrick","steelblue","darkgreen","orange","purple")[1:K]

  df_all <- do.call(rbind, lapply(1:K, function(j) {
    do.call(rbind, lapply(1:K, function(k) {
      probs <- sapply(y_grid, function(yp) {
        transition_probs(yp, params$gamma0[class_idx,,],
                         params$gamma1[class_idx,,], K)[j, k]
      })
      data.frame(y_prev = y_grid, prob = probs,
                 from = paste("De Regime", j),
                 to   = paste("Para Regime", k))
    }))
  }))

  ggplot(df_all, aes(x = y_prev, y = prob, color = to)) +
    geom_line(linewidth = 0.8) +
    facet_wrap(~ from, ncol = K) +
    scale_color_manual(values = cores) +
    labs(title  = paste("Probabilidades de Transição – Classe", class_idx),
         x = expression(y[t-1]), y = "Probabilidade", color = NULL) +
    theme_minimal(base_size = 11)
}


plot_bic <- function(model_sel) {
  ggplot(model_sel,
         aes(x = factor(K), y = BIC,
             group = factor(S), color = factor(S))) +
    geom_line(linewidth = 1) +
    geom_point(size = 3) +
    labs(title = "Seleção de Modelo via BIC",
         x = "Número de Regimes (K)", y = "BIC", color = "Classes (S)") +
    theme_minimal(base_size = 11)
}


plot_sync <- function(sync_matrix) {
  df <- melt(sync_matrix, na.rm = TRUE)
  ggplot(df, aes(x = Var1, y = Var2, fill = value)) +
    geom_tile(color = "white") +
    geom_text(aes(label = round(value, 2)), size = 3.5) +
    scale_fill_gradient2(low = "steelblue", mid = "white", high = "firebrick",
                         midpoint = 0, limits = c(-1, 1)) +
    labs(title = "Sincronização entre Ações",
         x = NULL, y = NULL, fill = "Correlação") +
    theme_minimal(base_size = 9) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}




# =============================================================================
# 8B. VISUALIZAÇÕES – ANÁLISE EXPLORATÓRIA
# =============================================================================


plot_descriptive_stats <- function(y_list) {
  if (!requireNamespace("moments", quietly = TRUE))
    stop("Instale o pacote 'moments': install.packages('moments')")
  library(moments)

  stats_df <- do.call(rbind, lapply(names(y_list), function(nm) {
    y  <- y_list[[nm]]
    jb <- jarque.test(y)
    data.frame(
      Acao       = nm,
      Media      = mean(y),
      Mediana    = median(y),
      DP         = sd(y),
      Assimetria = skewness(y),
      Curtose    = kurtosis(y) - 3,
      JB_stat    = jb$statistic,
      JB_pvalor  = jb$p.value,
      stringsAsFactors = FALSE
    )
  }))

  cat("\n=== Estatísticas Descritivas dos Retornos ===\n")
  print(round(stats_df[, -1], 4), row.names = FALSE)

  p_media <- ggplot(stats_df, aes(x = reorder(Acao, Media), y = Media)) +
    geom_col(fill = "steelblue") +
    coord_flip() +
    labs(title = "Retorno Médio Diário por Ação (%)", x = NULL, y = "Média (%)") +
    theme_minimal(base_size = 10)

  p_dp <- ggplot(stats_df, aes(x = reorder(Acao, DP), y = DP)) +
    geom_col(fill = "firebrick") +
    coord_flip() +
    labs(title = "Risco (Desvio-Padrão) por Ação (%)", x = NULL, y = "DP (%)") +
    theme_minimal(base_size = 10)

  grid.arrange(p_media, p_dp, ncol = 2)
  return(invisible(stats_df))
}


plot_return_series <- function(y_list, n_show = 6, window = 30) {
  nomes   <- names(y_list)[seq_len(min(n_show, length(y_list)))]
  roll_fn <- function(x, w) {
    res <- rep(NA, length(x))
    for (i in w:length(x)) res[i] <- mean(x[(i-w+1):i])
    res
  }
  roll_sd <- function(x, w) {
    res <- rep(NA, length(x))
    for (i in w:length(x)) res[i] <- sd(x[(i-w+1):i])
    res
  }

  plots_ret  <- list()
  plots_mean <- list()
  plots_sd   <- list()

  for (nm in nomes) {
    y  <- y_list[[nm]]
    df <- data.frame(t = seq_along(y), ret = y,
                     rm = roll_fn(y, window), rsd = roll_sd(y, window))
    plots_ret[[nm]] <- ggplot(df, aes(t, ret)) +
      geom_line(color = "steelblue", linewidth = 0.3) +
      labs(title = nm, x = NULL, y = "Retorno (%)") +
      theme_minimal(base_size = 9)
    plots_mean[[nm]] <- ggplot(df, aes(t, rm)) +
      geom_line(color = "darkgreen", linewidth = 0.4) +
      labs(title = nm, x = NULL, y = "Média móvel") +
      theme_minimal(base_size = 9)
    plots_sd[[nm]] <- ggplot(df, aes(t, rsd)) +
      geom_line(color = "firebrick", linewidth = 0.4) +
      labs(title = nm, x = NULL, y = "DP móvel") +
      theme_minimal(base_size = 9)
  }

  cat("\n--- Retornos (primeiras", n_show, "ações) ---\n")
  do.call(grid.arrange, c(plots_ret,  list(ncol = 3)))
  cat("\n--- Médias Móveis (janela =", window, "dias) ---\n")
  do.call(grid.arrange, c(plots_mean, list(ncol = 3)))
  cat("\n--- Desvios-Padrão Móveis (janela =", window, "dias) ---\n")
  do.call(grid.arrange, c(plots_sd,   list(ncol = 3)))
}


plot_return_correlation <- function(y_list, method = "pearson") {
  T_min   <- min(sapply(y_list, length))
  ret_mat <- do.call(cbind, lapply(y_list, tail, T_min))
  colnames(ret_mat) <- names(y_list)

  cor_mat <- cor(ret_mat, method = method, use = "complete.obs")

  cat("\n=== Correlação dos Retornos (", method, ") ===\n")
  print(round(cor_mat, 3))

  df_cor <- melt(cor_mat)
  p <- ggplot(df_cor, aes(x = Var1, y = Var2, fill = value)) +
    geom_tile(color = "white") +
    geom_text(aes(label = sprintf("%.2f", value)), size = 2.2) +
    scale_fill_gradient2(low = "steelblue", mid = "white", high = "firebrick",
                         midpoint = 0, limits = c(-1, 1)) +
    labs(title = "Correlação dos Retornos Logarítmicos (Pearson)",
         subtitle = "Retornos log são estacionários → correlação economicamente válida",
         x = NULL, y = NULL, fill = "Correlação") +
    theme_minimal(base_size = 8) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  print(p)
  return(invisible(cor_mat))
}


test_stationarity <- function(y_list) {
  if (!requireNamespace("tseries", quietly = TRUE))
    stop("Instale o pacote 'tseries': install.packages('tseries')")
  library(tseries)

  cat("\n=== Teste ADF de Estacionariedade (H0: raiz unitária) ===\n")
  cat(sprintf("%-22s  %10s  %8s  %s\n",
              "Ação", "Estatística", "p-valor", "Resultado"))
  cat(strrep("-", 58), "\n")

  adf_df <- do.call(rbind, lapply(names(y_list), function(nm) {
    res       <- adf.test(y_list[[nm]], alternative = "stationary")
    resultado <- ifelse(res$p.value < 0.05, "Estacionária", "Raiz unitária")
    cat(sprintf("%-22s  %10.4f  %8.4f  %s\n",
                nm, res$statistic, res$p.value, resultado))
    data.frame(Acao = nm, ADF = res$statistic, pvalor = res$p.value,
               Estacionario = res$p.value < 0.05,
               stringsAsFactors = FALSE)
  }))

  n_ok <- sum(adf_df$Estacionario)
  cat(sprintf("\n%d de %d séries rejeitam H0 (p < 0.05) → estacionárias.\n",
              n_ok, nrow(adf_df)))
  return(invisible(adf_df))
}


analyze_synchronization <- function(fit, series_names) {
  K           <- fit$K
  regime_bear <- which.min(fit$mu)
  regime_bull <- which.max(fit$mu)

  rotulo <- function(rk) {
    if (rk == regime_bear) "Bear"
    else if (rk == regime_bull) "Bull"
    else "Estável"
  }

  cat("\n=== Sincronização via Probabilidades Posteriores ===\n")

  for (rk in 1:K) {
    lbl <- rotulo(rk)
    cat("\n--- Regime", rk, "(", lbl, ") ---\n")
    sync <- compute_synchronization(fit$all_fb,
                                    series_names = series_names,
                                    regime_k     = rk)
    print(round(sync, 3))
    print(plot_sync(sync) +
            labs(title = paste("Sincronização – Regime", rk, "(", lbl, ")")))
  }
}




# =============================================================================
# 9. EXECUÇÃO PRINCIPAL — DADOS DO DJIA
# =============================================================================


run_djia <- function(arquivo_dados    = "dados_djia.RData",
                     S_range          = 1:3,
                     K_range          = 2:4,
                     n_init_selecao   = 5,
                     n_init_final     = 10,
                     n_show_series    = 6,
                     janela_movel     = 30) {

  # ---------------------------------------------------------------------------
  # CARREGAMENTO DOS DADOS
  # ---------------------------------------------------------------------------
  cat("=== Carregando dados do DJIA ===\n")

  if (!file.exists(arquivo_dados))
    stop(
      "Arquivo '", arquivo_dados, "' nao encontrado.\n",
      "Execute primeiro o script download_dados_djia.R para gerar os dados."
    )

  load(arquivo_dados)   # carrega: y_list, preco_list, tickers, nomes,
                        #          data_inicio, data_fim

  n_acoes <- length(y_list)
  T_obs   <- length(y_list[[1]])

  cat("Acoes carregadas :", n_acoes, "\n")
  cat("Periodo          :", data_inicio, "a", data_fim, "\n")
  cat("Observacoes (T)  :", T_obs, "dias uteis por acao\n\n")


  # ===========================================================================
  # BLOCO A — ANÁLISE EXPLORATÓRIA (pré-MGHMM)
  # ===========================================================================

  cat("\n========================================\n")
  cat("  BLOCO A - ANALISE EXPLORATORIA\n")
  cat("========================================\n")

  # A1. Estatísticas descritivas — replica Tabela 1 do artigo base
  cat("\n[A1] Estatisticas descritivas\n")
  desc <- plot_descriptive_stats(y_list)

  # A2. Séries de retornos + médias e desvios móveis — replica Figs. 2-4
  cat("\n[A2] Series temporais de retornos\n")
  plot_return_series(y_list, n_show = n_show_series, window = janela_movel)

  # A3. Testes de estacionariedade ADF
  # H0: raiz unitária. Esperamos rejeitar H0 para os retornos log.
  # Isso justifica NÃO usar preços em nível (não-estacionários → correlação
  # espúria) e SIM usar retornos log (estacionários → inferência válida).
  cat("\n[A3] Testes de estacionariedade ADF\n")
  adf_results <- test_stationarity(y_list)

  # A4. Correlação dos retornos logarítmicos
  # Retornos estacionários → correlação de Pearson economicamente interpretável.
  cat("\n[A4] Correlacao dos retornos logaritmicos\n")
  cor_ret <- plot_return_correlation(y_list, method = "pearson")


  # ===========================================================================
  # SELEÇÃO DE MODELO VIA BIC
  # ===========================================================================

  cat("\n========================================\n")
  cat("  SELECAO DE MODELO VIA BIC\n")
  cat("========================================\n")
  cat("Testando S =", paste(range(S_range), collapse = ":"),
      "| K =", paste(range(K_range), collapse = ":"), "\n\n")

  model_sel <- select_model_mghmm(y_list,
                                  S_range = S_range,
                                  K_range = K_range,
                                  n_init  = n_init_selecao,
                                  verbose = FALSE)
  cat("\n--- Tabela BIC completa ---\n")
  print(model_sel)
  print(plot_bic(model_sel))

  best_S <- model_sel$S[1]
  best_K <- model_sel$K[1]
  cat("\nMelhor modelo: S =", best_S, "| K =", best_K,
      "| BIC =", round(model_sel$BIC[1], 2), "\n\n")


  # ===========================================================================
  # ESTIMAÇÃO FINAL
  # ===========================================================================

  cat("\n========================================\n")
  cat("  ESTIMACAO FINAL DO MGHMM\n")
  cat("========================================\n")
  cat("S =", best_S, "| K =", best_K,
      "| Inicializacoes =", n_init_final, "\n\n")

  fit <- fit_mghmm(y_list, S = best_S, K = best_K,
                   n_init = n_init_final, verbose = TRUE)

  # Parâmetros estimados
  cat("\n=== Parametros estimados ===\n")
  cat("Medias por regime (mu):        ", round(fit$mu, 4), "\n")
  cat("Variancias por regime (sigma2):", round(fit$sigma2, 4), "\n")
  cat("Tamanhos das classes (pi_w):   ", round(fit$pi_w, 4), "\n")
  cat("Log-verossimilhanca:           ", round(fit$log_lik, 2), "\n")
  cat("No parametros:                 ", fit$n_params, "\n")
  cat("BIC:                           ", round(fit$BIC, 2), "\n")

  # Identificação dos regimes
  regime_bear    <- which.min(fit$mu)
  regime_bull    <- which.max(fit$mu)
  regime_estaveis <- setdiff(1:best_K, c(regime_bear, regime_bull))

  cat("\n--- Identificacao dos regimes ---\n")
  cat("Regime bear (menor retorno): Regime", regime_bear,
      "| mu =", round(fit$mu[regime_bear], 4),
      "| sigma2 =", round(fit$sigma2[regime_bear], 4), "\n")
  cat("Regime bull (maior retorno): Regime", regime_bull,
      "| mu =", round(fit$mu[regime_bull], 4),
      "| sigma2 =", round(fit$sigma2[regime_bull], 4), "\n")
  if (length(regime_estaveis) > 0)
    cat("Regime(s) estavel(is):       Regime(s)",
        paste(regime_estaveis, collapse = ", "), "\n")

  # Ocupação média por regime
  cat("\n--- Ocupacao media por regime ---\n")
  occ <- sapply(1:best_K, function(k)
    mean(sapply(fit$all_fb, function(fb) {
      T_ <- dim(fb$post_wzt)[1]
      S_ <- dim(fb$post_wzt)[2]
      mean(rowSums(matrix(fb$post_wzt[,,k], T_, S_)))
    })))
  names(occ) <- paste("Regime", 1:best_K)
  print(round(occ, 4))

  # Visualização de regimes por ação
  cat("\n--- Visualizacao de regimes por acao ---\n")
  for (i in seq_along(y_list)) {
    plot_regimes(y_list[[i]], fit$all_fb[[i]],
                 series_name = names(y_list)[i], K = best_K)
  }

  # Probabilidades de transição por classe
  for (s in 1:best_S)
    print(plot_transition_probs(fit, class_idx = s))

  # Classe modal por ação
  if (best_S > 1) {
    cat("\n--- Classe modal por acao ---\n")
    cat(sprintf("%-22s  %6s  %s\n", "Acao", "Classe", "P(w)"))
    cat(strrep("-", 55), "\n")
    for (i in seq_along(y_list)) {
      cls <- which.max(fit$all_fb[[i]]$post_w)
      prb <- round(fit$all_fb[[i]]$post_w, 3)
      cat(sprintf("%-22s  %6d  %s\n",
                  names(y_list)[i], cls, paste(prb, collapse = " / ")))
    }
  }


  # ===========================================================================
  # BLOCO B — ANÁLISE DE SINCRONIZAÇÃO (pós-MGHMM)
  # ===========================================================================

  cat("\n========================================\n")
  cat("  BLOCO B - ANALISE DE SINCRONIZACAO\n")
  cat("========================================\n")

  # B1. Sincronização via probabilidades posteriores — replica Tabela 5
  # Medida: correlação dos logits das probabilidades posteriores.
  # Robusta a outliers (cf. Dias et al., 2015, p. 9).
  cat("\n[B1] Sincronizacao por regime (equivalente a Tabela 5 do artigo base)\n")
  analyze_synchronization(fit, series_names = names(y_list))

  # B2. Comparação retornos (A4) vs. sincronização de regimes (B1)
  # Discrepâncias revelam dependências não lineares e de cauda capturadas
  # pelo MGHMM mas invisíveis na correlação linear de Pearson.
  cat("\n[B2] Comparacao: correlacao de retornos vs. sincronizacao de regimes\n")
  cat("  Discrepancias entre A4 e B1 indicam dependencias de cauda\n")
  cat("  capturadas pelo MGHMM mas nao pela correlacao linear.\n")


  # ===========================================================================
  # PREVISÃO 1-PASSO-À-FRENTE
  # ===========================================================================

  cat("\n========================================\n")
  cat("  PREVISAO 1-PASSO-A-FRENTE\n")
  cat("========================================\n")
  cat(sprintf("%-22s  %8s  %8s\n", "Acao", "RMSE", "MAE"))
  cat(strrep("-", 42), "\n")

  prev_df <- do.call(rbind, lapply(seq_along(y_list), function(i) {
    nm    <- names(y_list)[i]
    y_hat <- predict_mghmm(y_list[[i]], fit)
    rmse  <- sqrt(mean((y_list[[i]] - y_hat)^2))
    mae   <- mean(abs(y_list[[i]] - y_hat))
    cat(sprintf("%-22s  %8.4f  %8.4f\n", nm, rmse, mae))
    data.frame(Acao = nm, RMSE = rmse, MAE = mae)
  }))

  cat(sprintf("%-22s  %8.4f  %8.4f\n",
              "MEDIA GERAL",
              mean(prev_df$RMSE), mean(prev_df$MAE)))


  # ===========================================================================
  # RETORNO
  # ===========================================================================

  sync_bear <- compute_synchronization(fit$all_fb,
                                       series_names = names(y_list),
                                       regime_k     = regime_bear)

  return(invisible(list(
    fit       = fit,
    model_sel = model_sel,
    desc      = desc,
    adf       = adf_results,
    cor_ret   = cor_ret,
    sync_bear = sync_bear,
    prev      = prev_df
  )))
}


# =============================================================================
# EXECUTAR
# =============================================================================
result <- run_djia()
