# =============================================================================
# utils.R
# MS-VAR with Partially Shared Regimes — asset-specific idiosyncratic states
# Regimes: K=3 (bear=1, stable=2, bull=3)
# Each asset i has M idiosyncratic local states with own params
# Innovations: skew-t (Azzalini-Capitanio)
# =============================================================================

library(sn)
library(igraph)
library(Matrix)

# -----------------------------------------------------------------------------
# ESTRUTURA DE PARÂMETROS
#
# params_emit : lista de comprimento n, onde cada elemento i é uma lista
#               de comprimento K + M:
#               [[k]]   para k=1..K : estado global (compartilhado entre ativos,
#                                     mas cada ativo pode ter sigma/alpha/nu próprios)
#               [[K+m]] para m=1..M : estado idiossincrático do ativo i
#
# mu  : lista de comprimento K + M, cada elemento é vetor de comprimento n
#       Para estados idiossincráticos, apenas mu[[K+m]][i] é usado para ativo i
#
# A   : lista de comprimento K + M, cada elemento é matriz n x n
#       Para estados idiossincráticos, apenas A[[K+m]][i,] é usado para ativo i
#
# rho : matriz n x K  — rho[i,k] = P(s_{it}=1 | z_t=k, asset i)
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# 1. SKEW-T LOG-DENSITY
# -----------------------------------------------------------------------------

skewt_logdens <- function(x, sigma, alpha, nu) {
  dst(x, xi = 0, omega = sigma, alpha = alpha, nu = nu, log = TRUE)
}

# -----------------------------------------------------------------------------
# 2. COMPUTE VAR RESIDUALS
# -----------------------------------------------------------------------------

compute_residuals <- function(yt, yt1, mu_k, A_k) {
  as.numeric(yt - mu_k - A_k %*% yt1)
}

# -----------------------------------------------------------------------------
# 3. JOINT EMISSION LOG-DENSITY — asset-specific idiosyncratic states
#
# params_emit[[i]] : list of length K+M for asset i
#   params_emit[[i]][[k]]   : global state k params for asset i
#   params_emit[[i]][[K+m]] : idiosyncratic state m params for asset i
# -----------------------------------------------------------------------------

joint_logdens <- function(yt, yt1, k, mu, A, rho, params_emit,
                          K = 3, M = 2) {
  n        <- length(yt)
  eps_glob <- compute_residuals(yt, yt1, mu[[k]], A[[k]])
  log_dens <- 0
  
  for (i in seq_len(n)) {
    
    # --- Component 1: s_{it}=1, follows global regime k ---
    p_glob  <- params_emit[[i]][[k]]
    f_glob  <- exp(skewt_logdens(eps_glob[i],
                                 p_glob$sigma, p_glob$alpha, p_glob$nu))
    
    # --- Component 2: s_{it}=0, follows asset-i idiosyncratic state ---
    f_idio <- 0
    for (m in seq_len(M)) {
      # Asset-i idiosyncratic intercept and AR row
      mu_im  <- mu[[K + m]][i]
      a_im   <- A[[K + m]][i, ]
      eps_im <- yt[i] - mu_im - sum(a_im * yt1)
      
      p_idio <- params_emit[[i]][[K + m]]
      f_idio <- f_idio + (1 / M) * exp(
        skewt_logdens(eps_im, p_idio$sigma, p_idio$alpha, p_idio$nu)
      )
    }
    
    f_mix    <- rho[i, k] * f_glob + (1 - rho[i, k]) * f_idio
    log_dens <- log_dens + log(max(f_mix, .Machine$double.xmin))
  }
  
  log_dens
}

# -----------------------------------------------------------------------------
# 4. HAMILTON FILTER
# -----------------------------------------------------------------------------

hamilton_filter <- function(Y, P_trans, pi0, mu, A, rho,
                            params_emit, K = 3, M = 2) {
  TT       <- nrow(Y)
  xi_filt  <- matrix(0, TT, K)
  xi_pred  <- matrix(0, TT, K)
  xi_prev  <- pi0
  log_lik  <- 0
  
  for (t in seq_len(TT)) {
    yt  <- Y[t, ]
    yt1 <- if (t == 1) rep(0, ncol(Y)) else Y[t - 1, ]
    
    xi_p <- as.vector(t(P_trans) %*% xi_prev)
    xi_pred[t, ] <- xi_p
    
    f_k <- vapply(seq_len(K), function(k) {
      exp(joint_logdens(yt, yt1, k, mu, A, rho, params_emit, K, M))
    }, numeric(1))
    
    denom        <- max(sum(xi_p * f_k), .Machine$double.xmin)
    xi_filt[t, ] <- (xi_p * f_k) / denom
    log_lik      <- log_lik + log(denom)
    xi_prev      <- xi_filt[t, ]
  }
  
  list(xi_filt = xi_filt, xi_pred = xi_pred, log_lik = log_lik)
}

# -----------------------------------------------------------------------------
# 5. BAUM-WELCH BACKWARD PASS
# -----------------------------------------------------------------------------

baum_welch_smooth <- function(xi_filt, xi_pred, P_trans, K = 3) {
  TT         <- nrow(xi_filt)
  xi_smooth  <- matrix(0, TT, K)
  xi_smooth[TT, ] <- xi_filt[TT, ]
  
  for (t in (TT - 1):1) {
    ratio <- xi_smooth[t + 1, ] /
      pmax(xi_pred[t + 1, ], .Machine$double.xmin)
    for (j in seq_len(K)) {
      xi_smooth[t, j] <- xi_filt[t, j] * sum(P_trans[j, ] * ratio)
    }
    s <- sum(xi_smooth[t, ])
    if (s > 0) xi_smooth[t, ] <- xi_smooth[t, ] / s
  }
  xi_smooth
}

# -----------------------------------------------------------------------------
# 6. SMOOTHED JOINT PROBABILITIES eta[t, j, k]
# -----------------------------------------------------------------------------

smooth_joint_probs <- function(xi_filt, xi_pred, xi_smooth,
                               P_trans, K = 3) {
  TT  <- nrow(xi_filt)
  eta <- array(0, dim = c(TT, K, K))
  
  for (t in 2:TT) {
    ratio <- xi_smooth[t, ] /
      pmax(xi_pred[t, ], .Machine$double.xmin)
    for (j in seq_len(K)) {
      for (k in seq_len(K)) {
        eta[t, j, k] <- xi_filt[t - 1, j] * P_trans[j, k] * ratio[k]
      }
    }
    s <- sum(eta[t, , ])
    if (s > 0) eta[t, , ] <- eta[t, , ] / s
  }
  eta
}

# -----------------------------------------------------------------------------
# 7. M-STEP: TRANSITION MATRIX
# -----------------------------------------------------------------------------

mstep_P <- function(eta, K = 3, P_fixed_zeros = NULL) {
  num   <- apply(eta[2:dim(eta)[1], , , drop = FALSE], c(2, 3), sum)
  P_new <- num / pmax(rowSums(num), .Machine$double.xmin)
  
  if (!is.null(P_fixed_zeros)) {
    for (idx in P_fixed_zeros) P_new[idx[1], idx[2]] <- 0
    rs <- rowSums(P_new)
    for (j in seq_len(K)) {
      if (rs[j] > 0) P_new[j, ] <- P_new[j, ] / rs[j]
    }
  }
  P_new
}

# -----------------------------------------------------------------------------
# 8. M-STEP: mu_k AND A_k via WLS
# -----------------------------------------------------------------------------

mstep_mu_A <- function(Y, xi_smooth, k) {
  TT <- nrow(Y)
  n  <- ncol(Y)
  w  <- xi_smooth[2:TT, k]
  
  X  <- cbind(1, Y[1:(TT - 1), ])
  Yw <- Y[2:TT, ]
  
  XtW <- t(X) %*% diag(w)
  B   <- solve(XtW %*% X + diag(1e-8, ncol(X)), XtW %*% Yw)
  
  list(mu = as.numeric(B[1, ]),
       A  = t(B[2:(n + 1), ]))
}

# -----------------------------------------------------------------------------
# 9. M-STEP: rho_{ik} — asset-specific participation probabilities
# -----------------------------------------------------------------------------

mstep_rho <- function(Y, xi_smooth, mu, A, params_emit,
                      K = 3, M = 2) {
  TT  <- nrow(Y)
  n   <- ncol(Y)
  
  # Acumuladores separados: numerador e denominador
  num_rho <- matrix(0, n, K)
  den_rho <- matrix(0, n, K)
  
  for (t in 2:TT) {
    yt  <- Y[t, ]
    yt1 <- Y[t - 1, ]
    
    for (k in seq_len(K)) {
      eps_glob <- compute_residuals(yt, yt1, mu[[k]], A[[k]])
      w_k      <- xi_smooth[t, k]
      
      for (i in seq_len(n)) {
        # Likelihood under global state k for asset i
        p_g <- params_emit[[i]][[k]]
        f_g <- exp(skewt_logdens(eps_glob[i],
                                 p_g$sigma, p_g$alpha, p_g$nu))
        f_g <- max(f_g, .Machine$double.xmin)
        
        # Likelihood under idiosyncratic mixture for asset i
        f_id <- 0
        for (m in seq_len(M)) {
          mu_im  <- mu[[K + m]][i]
          a_im   <- A[[K + m]][i, ]
          eps_im <- yt[i] - mu_im - sum(a_im * yt1)
          p_id   <- params_emit[[i]][[K + m]]
          f_id   <- f_id + (1 / M) * exp(
            skewt_logdens(eps_im, p_id$sigma, p_id$alpha, p_id$nu)
          )
        }
        f_id <- max(f_id, .Machine$double.xmin)
        
        # E[s_{it} | data, z_t=k] usando rho da iteração anterior
        # Extraído directamente de params_emit (evita rho=0 na 1ª iter)
        rho_prev <- 0.5   # prior neutro — será sobrescrito após 1ª iter
        # Recuperar rho actual do acumulador se já há massa suficiente
        if (den_rho[i, k] > 1e-6) {
          rho_prev <- num_rho[i, k] / den_rho[i, k]
        }
        rho_prev <- max(min(rho_prev, 0.99), 0.01)
        
        # Posterior de s_{it} = 1
        num_s <- rho_prev * f_g
        den_s <- num_s + (1 - rho_prev) * f_id
        den_s <- max(den_s, .Machine$double.xmin)
        e_s   <- num_s / den_s
        
        # Acumula ponderado por xi_{t|T}^(k)
        num_rho[i, k] <- num_rho[i, k] + w_k * e_s
        den_rho[i, k] <- den_rho[i, k] + w_k
      }
    }
  }
  
  # rho_{ik} = E[s_{it}=1 | data] marginalizado sobre t
  rho <- num_rho / pmax(den_rho, .Machine$double.xmin)
  rho[is.nan(rho) | is.na(rho)] <- 0.5
  pmin(pmax(rho, 0.01), 0.99)
}

# -----------------------------------------------------------------------------
# 10. M-STEP: SKEW-T PARAMS — asset-specific, per effective state
# Returns updated params_emit (list of n assets x (K+M) states)
# -----------------------------------------------------------------------------

mstep_skewt_asset <- function(Y, xi_smooth, mu, A, rho,
                              params_emit, K = 3, M = 2) {
  TT <- nrow(Y)
  n  <- ncol(Y)
  
  # Accumulate residuals and weights per asset per state
  res_arr <- vector("list", n)
  wgt_arr <- vector("list", n)
  for (i in seq_len(n)) {
    res_arr[[i]] <- vector("list", K + M)
    wgt_arr[[i]] <- vector("list", K + M)
    for (s in seq_len(K + M)) {
      res_arr[[i]][[s]] <- numeric(0)
      wgt_arr[[i]][[s]] <- numeric(0)
    }
  }
  
  for (t in 2:TT) {
    yt  <- Y[t, ]
    yt1 <- Y[t - 1, ]
    
    for (k in seq_len(K)) {
      eps_glob <- compute_residuals(yt, yt1, mu[[k]], A[[k]])
      w_k      <- xi_smooth[t, k]
      
      for (i in seq_len(n)) {
        # Weight for global state: xi_{t|T}^(k) * rho_{ik}
        w_ik <- w_k * rho[i, k]
        res_arr[[i]][[k]] <- c(res_arr[[i]][[k]], eps_glob[i])
        wgt_arr[[i]][[k]] <- c(wgt_arr[[i]][[k]], w_ik)
      }
    }
    
    for (m in seq_len(M)) {
      s_idx <- K + m
      for (i in seq_len(n)) {
        mu_im  <- mu[[s_idx]][i]
        a_im   <- A[[s_idx]][i, ]
        eps_im <- yt[i] - mu_im - sum(a_im * yt1)
        
        # Weight for idiosyncratic state: mean xi * (1 - mean rho_i)
        w_im <- mean(xi_smooth[t, ]) * (1 - mean(rho[i, ]))
        res_arr[[i]][[s_idx]] <- c(res_arr[[i]][[s_idx]], eps_im)
        wgt_arr[[i]][[s_idx]] <- c(wgt_arr[[i]][[s_idx]], w_im)
      }
    }
  }
  
  # Optimise skew-t params per asset per state
  params_new <- vector("list", n)
  for (i in seq_len(n)) {
    params_new[[i]] <- vector("list", K + M)
    for (s in seq_len(K + M)) {
      r <- res_arr[[i]][[s]]
      w <- wgt_arr[[i]][[s]]
      sw <- sum(w)
      
      if (sw < 1e-6 || length(r) < 5) {
        params_new[[i]][[s]] <- params_emit[[i]][[s]]
        next
      }
      
      neg_ll <- function(par) {
        sigma <- exp(par[1])
        alpha <- par[2]
        nu    <- exp(par[3]) + 2
        -sum(w * dst(r, xi = 0, omega = sigma,
                     alpha = alpha, nu = nu, log = TRUE))
      }
      
      p0   <- params_emit[[i]][[s]]
      init <- c(log(p0$sigma), p0$alpha, log(max(p0$nu - 2, 1)))
      opt  <- tryCatch(
        optim(init, neg_ll, method = "BFGS",
              control = list(maxit = 200, reltol = 1e-6)),
        error = function(e) list(par = init, convergence = 1)
      )
      
      par <- opt$par
      params_new[[i]][[s]] <- list(
        sigma = exp(par[1]),
        alpha = par[2],
        nu    = exp(par[3]) + 2
      )
    }
  }
  
  params_new
}

# -----------------------------------------------------------------------------
# 11. FULL EM ALGORITHM
# -----------------------------------------------------------------------------

em_msvar <- function(Y,
                     K          = 3,
                     M          = 2,
                     max_iter   = 100,
                     tol        = 1e-5,
                     P_zeros    = list(c(1, 3), c(3, 1)),
                     init       = NULL,
                     verbose    = TRUE) {
  TT <- nrow(Y)
  n  <- ncol(Y)
  
  # --- Initialisation -------------------------------------------------------
  if (is.null(init)) {
    set.seed(42)
    
    P_trans <- matrix(1 / K, K, K)
    for (idx in P_zeros) P_trans[idx[1], idx[2]] <- 0
    P_trans <- P_trans / rowSums(P_trans)
    
    pi0 <- rep(1 / K, K)
    
    # Global regime intercepts and AR matrices
    mu_vals <- c(-0.3, 0.0, 0.3, -0.1, 0.1)
    mu <- lapply(seq_len(K + M), function(s) rep(mu_vals[s], n))
    A  <- lapply(seq_len(K + M), function(s) diag(0.05, n))
    
    rho <- matrix(0.8, n, K)
    
    # Asset-specific emission params: list[[i]][[s]]
    sigma_vals <- c(2.0, 0.5, 1.0, 1.5, 1.5)
    alpha_vals <- c(-0.5, 0.0, 0.3, 0.0, 0.0)
    params_emit <- lapply(seq_len(n), function(i) {
      lapply(seq_len(K + M), function(s) {
        list(sigma = sigma_vals[s],
             alpha = alpha_vals[s],
             nu    = 5)
      })
    })
    
  } else {
    P_trans     <- init$P_trans
    pi0         <- init$pi0
    mu          <- init$mu
    A           <- init$A
    rho         <- init$rho
    params_emit <- init$params_emit
  }
  
  log_lik_prev <- -Inf
  history      <- numeric(max_iter)
  
  # --- EM Loop --------------------------------------------------------------
  for (iter in seq_len(max_iter)) {
    
    # E-STEP
    filt   <- hamilton_filter(Y, P_trans, pi0, mu, A, rho,
                              params_emit, K, M)
    smooth <- baum_welch_smooth(filt$xi_filt, filt$xi_pred,
                                P_trans, K)
    eta    <- smooth_joint_probs(filt$xi_filt, filt$xi_pred,
                                 smooth, P_trans, K)
    
    log_lik       <- filt$log_lik
    history[iter] <- log_lik
    
    if (verbose) cat(sprintf("Iter %3d | log-lik = %.4f\n", iter, log_lik))
    
    if (abs(log_lik - log_lik_prev) < tol && iter > 5) {
      if (verbose) cat("Converged.\n")
      break
    }
    log_lik_prev <- log_lik
    
    # M-STEP
    P_trans <- mstep_P(eta, K, P_zeros)
    
    for (k in seq_len(K)) {
      res_k   <- mstep_mu_A(Y, smooth, k)
      mu[[k]] <- res_k$mu
      A[[k]]  <- res_k$A
    }
    
    rho         <- mstep_rho(Y, smooth, mu, A, params_emit, K, M)
    params_emit <- mstep_skewt_asset(Y, smooth, mu, A, rho,
                                     params_emit, K, M)
    pi0 <- smooth[1, ]
  }
  
  list(
    P_trans     = P_trans,
    pi0         = pi0,
    mu          = mu,
    A           = A,
    rho         = rho,
    params_emit = params_emit,
    xi_filt     = filt$xi_filt,
    xi_smooth   = smooth,
    log_lik     = log_lik,
    history     = history[1:iter]
  )
}

# -----------------------------------------------------------------------------
# 12. INFLUENCE MATRIX W^(k)
# -----------------------------------------------------------------------------

influence_matrix <- function(A_list, alpha_sig = 0.05,
                             vcov_list = NULL, K = 3) {
  n <- nrow(A_list[[1]])
  W <- vector("list", K)
  for (k in seq_len(K)) {
    Wk <- matrix(0, n, n)
    for (i in seq_len(n)) {
      for (l in seq_len(n)) {
        if (i == l) next
        coef_val <- A_list[[k]][l, i]
        sig <- if (!is.null(vcov_list)) {
          se <- sqrt(vcov_list[[k]][l + (i-1)*n, l + (i-1)*n])
          2 * (1 - pnorm(abs(coef_val / se))) < alpha_sig
        } else TRUE
        Wk[i, l] <- abs(coef_val) * as.numeric(sig)
      }
    }
    W[[k]] <- Wk
  }
  W
}

# -----------------------------------------------------------------------------
# 13. DYNAMIC LEADERSHIP SCORE h_{it}
# -----------------------------------------------------------------------------

dynamic_leadership <- function(W_list, xi_filt, K = 3) {
  n            <- nrow(W_list[[1]])
  out_strength <- sapply(seq_len(K), function(k) rowSums(W_list[[k]]))
  # out_strength: n x K; xi_filt: T x K
  h <- xi_filt %*% t(out_strength)   # T x n
  h
}

# -----------------------------------------------------------------------------
# 14. REGIME-WEIGHTED PAGERANK r_t
# -----------------------------------------------------------------------------

regime_pagerank <- function(W_list, xi_filt, K = 3) {
  pr_k <- sapply(seq_len(K), function(k) {
    g  <- graph_from_adjacency_matrix(
      W_list[[k]], mode = "directed", weighted = TRUE)
    page_rank(g, directed = TRUE)$vector
  })
  # pr_k: n x K
  r <- xi_filt %*% t(pr_k)   # T x n
  r
}

# -----------------------------------------------------------------------------
# 15. WALD TEST — Newey-West SE (corrigido)
# -----------------------------------------------------------------------------

wald_test_leadership <- function(Y, xi_smooth, k, lags_nw = 5) {
  TT  <- nrow(Y)
  n   <- ncol(Y)
  TT2 <- TT - 1
  
  w   <- xi_smooth[2:TT, k]
  X   <- cbind(1, Y[1:(TT - 1), ])
  Yw  <- Y[2:TT, ]
  
  XtW <- t(X) %*% diag(w)
  B   <- solve(XtW %*% X + diag(1e-8, ncol(X)), XtW %*% Yw)
  
  p_matrix <- matrix(NA, n, n)
  
  for (l in seq_len(n)) {
    resid_l <- as.numeric(Yw[, l] - X %*% B[, l])
    score   <- X * (w * resid_l)   # TT2 x (n+1)
    
    S <- t(score) %*% score / TT2
    for (lag in seq_len(lags_nw)) {
      idx1      <- 1:(TT2 - lag)
      idx2      <- (1 + lag):TT2
      gamma_lag <- t(score[idx1, , drop = FALSE]) %*%
        score[idx2, , drop = FALSE] / TT2
      w_bart    <- 1 - lag / (lags_nw + 1)
      S         <- S + w_bart * (gamma_lag + t(gamma_lag))
    }
    
    XtWX_inv <- solve(XtW %*% X + diag(1e-8, ncol(X)))
    vcov_l   <- XtWX_inv %*% S %*% XtWX_inv * TT2
    
    for (i in seq_len(n)) {
      if (i == l) next
      coef_idx  <- i + 1
      coef_val  <- B[coef_idx, l]
      se_val    <- sqrt(max(vcov_l[coef_idx, coef_idx], 1e-12))
      wald_stat <- (coef_val / se_val)^2
      p_matrix[i, l] <- 1 - pchisq(wald_stat, df = 1)
    }
  }
  
  list(p_values = p_matrix, B = B)
}

# -----------------------------------------------------------------------------
# 16. BIC
# -----------------------------------------------------------------------------

compute_bic <- function(log_lik, K, n, p = 1, TT) {
  N_params <- K * (K - 1) +
    K * n^2 * p +
    K * n * (n + 1) / 2 +
    K * n +
    n * K +
    2 * n * K * (1 + 1)  # skew-t params agora por ativo
  -2 * log_lik + log(TT) * N_params
}