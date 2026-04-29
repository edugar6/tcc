# =============================================================================
# HIERARQUIA DINÂMICA AO LONGO DO TEMPO
# =============================================================================

library(ggplot2)
library(reshape2)

# -----------------------------------------------------------------------------
# 1. RANK DINÂMICO: converte h_{it} em ranking por período
# rank=1 significa maior liderança no período t
# -----------------------------------------------------------------------------

rank_mat <- t(apply(h_mat, 1, function(row) rank(-row, ties.method = "min")))
colnames(rank_mat) <- colnames(Y)

# Converte para data.frame longo
df_rank <- as.data.frame(rank_mat)
df_rank$t <- seq_len(nrow(df_rank))
df_rank_long <- melt(df_rank, id.vars = "t",
                     variable.name = "Asset",
                     value.name    = "Rank")

# -----------------------------------------------------------------------------
# 2. SCORE DINÂMICO h_{it} em formato longo
# -----------------------------------------------------------------------------

df_h <- as.data.frame(h_mat)
df_h$t <- seq_len(nrow(df_h))
df_h_long <- melt(df_h, id.vars = "t",
                  variable.name = "Asset",
                  value.name    = "Score")

# -----------------------------------------------------------------------------
# 3. REGIME DOMINANTE por período (para faixa de fundo)
# -----------------------------------------------------------------------------

regime_dominant <- apply(fit$xi_filt, 1, which.max)
# Mapeia para label alinhado
regime_label    <- factor(best_map[regime_dominant],
                          levels = 1:K,
                          labels = state_names)

df_regime <- data.frame(
  t     = seq_len(TT),
  regime = regime_label
)

regime_fills <- c(Bear = "#FDDCDC", Stable = "#E8E8E8", Bull = "#D6EAF8")

# -----------------------------------------------------------------------------
# 4. FIGURA 1: Score de liderança h_{it} com fundo de regime
# -----------------------------------------------------------------------------

p1 <- ggplot() +
  # Faixa de fundo por regime dominante
  geom_tile(data = df_regime,
            aes(x = t, y = 0, fill = regime, height = Inf),
            alpha = 0.35, width = 1) +
  scale_fill_manual(values = regime_fills, name = "Regime") +
  
  # Linhas de score
  geom_line(data = df_h_long,
            aes(x = t, y = Score, colour = Asset),
            linewidth = 0.7) +
  scale_colour_manual(values = c("steelblue", "tomato", "forestgreen")) +
  
  labs(title    = "Dynamic Leadership Score h_{it} over Time",
       subtitle = "Background shading = dominant regime",
       x = "Time", y = expression(h[it])) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom",
        panel.grid.minor = element_blank())

print(p1)

# -----------------------------------------------------------------------------
# 5. FIGURA 2: Ranking dinâmico (bump chart)
# -----------------------------------------------------------------------------

p2 <- ggplot() +
  geom_tile(data = df_regime,
            aes(x = t, y = 0, fill = regime, height = Inf),
            alpha = 0.30, width = 1) +
  scale_fill_manual(values = regime_fills, name = "Regime") +
  
  geom_line(data = df_rank_long,
            aes(x = t, y = Rank, colour = Asset, group = Asset),
            linewidth = 0.8) +
  geom_point(data = df_rank_long[df_rank_long$t %% 50 == 0, ],
             aes(x = t, y = Rank, colour = Asset),
             size = 2) +
  
  scale_colour_manual(values = c("steelblue", "tomato", "forestgreen")) +
  scale_y_reverse(breaks = seq_len(n),
                  labels = paste("Rank", seq_len(n))) +
  
  labs(title    = "Dynamic Leadership Ranking over Time",
       subtitle = "Rank 1 = highest leadership in period t",
       x = "Time", y = "Leadership Rank") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom",
        panel.grid.minor = element_blank())

print(p2)

# -----------------------------------------------------------------------------
# 6. FIGURA 3: Proporção do tempo em cada rank por regime
# Responde: "qual ativo lidera MAIS em cada regime?"
# -----------------------------------------------------------------------------

df_rank_regime <- merge(df_rank_long,
                        df_regime, by = "t")

# Proporção de períodos em rank=1 por ativo e regime
prop_top <- aggregate(
  Rank == 1 ~ Asset + regime,
  data = df_rank_regime,
  FUN  = mean
)
names(prop_top)[3] <- "Prop_Rank1"

p3 <- ggplot(prop_top,
             aes(x = regime, y = Prop_Rank1,
                 fill = Asset)) +
  geom_col(position = "dodge", colour = "white", linewidth = 0.4) +
  scale_fill_manual(values = c("steelblue", "tomato", "forestgreen")) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(title    = "Proportion of Time as Top Leader by Regime",
       subtitle = "Fraction of periods where each asset holds Rank 1",
       x = "Regime", y = "% of periods as leader") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

print(p3)

# -----------------------------------------------------------------------------
# 7. FIGURA 4: Heatmap de h_{it} — tempo x ativo (com regimes no eixo x)
# -----------------------------------------------------------------------------

# Suaviza h_{it} com média móvel de 20 períodos para visualização
smooth_h <- function(x, w = 20) stats::filter(x, rep(1/w, w), sides = 2)

df_h_smooth <- as.data.frame(apply(h_mat, 2, smooth_h, w = 20))
df_h_smooth$t <- seq_len(nrow(df_h_smooth))
df_h_smooth   <- na.omit(df_h_smooth)

df_hs_long <- melt(df_h_smooth, id.vars = "t",
                   variable.name = "Asset",
                   value.name    = "Score_smooth")

p4 <- ggplot(df_hs_long,
             aes(x = t, y = Asset, fill = Score_smooth)) +
  geom_tile() +
  scale_fill_gradient2(low  = "white",
                       mid  = "steelblue",
                       high = "navy",
                       midpoint = median(df_hs_long$Score_smooth,
                                         na.rm = TRUE),
                       name = "h_{it}\n(smoothed)") +
  geom_vline(data = data.frame(
    t      = which(diff(as.integer(regime_label)) != 0),
    regime = regime_label[which(diff(as.integer(regime_label)) != 0)]
  ),
  aes(xintercept = t),
  colour = "grey40", linewidth = 0.3, alpha = 0.5) +
  labs(title    = "Leadership Score Heatmap (20-period MA)",
       subtitle = "Vertical lines = regime transitions",
       x = "Time", y = "Asset") +
  theme_bw(base_size = 12) +
  theme(legend.position = "right",
        axis.text.y     = element_text(size = 11))

print(p4)

# -----------------------------------------------------------------------------
# 8. SUMÁRIO NUMÉRICO: estabilidade da hierarquia
# -----------------------------------------------------------------------------

cat("\n=== Estabilidade da hierarquia ao longo do tempo ===\n")

# Frequência de inversões de ranking entre t e t+1
inversions <- sum(apply(rank_mat, 1, function(r) {
  any(diff(r) != 0)
})) / TT

cat(sprintf("Frequência de mudanças de ranking: %.1f%%\n",
            100 * inversions))

# Fracção do tempo em que cada ativo ocupa cada rank
cat("\nFracção do tempo em cada rank:\n")
for (a in seq_len(n)) {
  cat(sprintf("  %s: ", colnames(Y)[a]))
  for (r in seq_len(n)) {
    cat(sprintf("Rank%d=%.1f%%  ", r,
                100 * mean(rank_mat[, a] == r)))
  }
  cat("\n")
}

# Por regime
cat("\nFracção do tempo como líder (Rank 1) por regime:\n")
print(round(
  tapply(df_rank_regime$Rank == 1,
         list(df_rank_regime$Asset, df_rank_regime$regime),
         mean) * 100, 1
))