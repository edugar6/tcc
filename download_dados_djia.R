# =============================================================================
# DOWNLOAD E PREPARAÇÃO DOS DADOS – 30 AÇÕES DO DJIA
# =============================================================================
# Execute este script UMA VEZ para baixar e salvar os dados localmente.
# Depois, o script principal (mghmm_djia.R) carrega o arquivo salvo em vez
# de baixar da internet a cada execução.
# =============================================================================

# install.packages(c("quantmod", "zoo"))
library(quantmod)
library(zoo)

# =============================================================================
# 1. TICKERS E NOMES DAS 30 AÇÕES DO DJIA
# =============================================================================
# Composição do DJIA vigente (atualizada para 2024).

tickers <- c(
  "AAPL",  # Apple
  "AMGN",  # Amgen
  "AXP",   # American Express
  "BA",    # Boeing
  "CAT",   # Caterpillar
  "CRM",   # Salesforce
  "CSCO",  # Cisco
  "CVX",   # Chevron
  "DIS",   # Disney
  "DOW",   # Dow Inc.
  "GS",    # Goldman Sachs
  "HD",    # Home Depot
  "HON",   # Honeywell
  "IBM",   # IBM
  "INTC",  # Intel
  "JNJ",   # Johnson & Johnson
  "JPM",   # JPMorgan Chase
  "KO",    # Coca-Cola
  "MCD",   # McDonald's
  "MMM",   # 3M
  "MRK",   # Merck
  "MSFT",  # Microsoft
  "NKE",   # Nike
  "PG",    # Procter & Gamble
  "TRV",   # Travelers
  "UNH",   # UnitedHealth
  "V",     # Visa
  "VZ",    # Verizon
  "WBA",   # Walgreens Boots Alliance
  "WMT"    # Walmart
)

nomes <- c(
  "Apple", "Amgen", "AmericanExpress", "Boeing", "Caterpillar",
  "Salesforce", "Cisco", "Chevron", "Disney", "Dow",
  "GoldmanSachs", "HomeDepot", "Honeywell", "IBM", "Intel",
  "JohnsonJohnson", "JPMorgan", "CocaCola", "McDonalds", "3M",
  "Merck", "Microsoft", "Nike", "ProcterGamble", "Travelers",
  "UnitedHealth", "Visa", "Verizon", "Walgreens", "Walmart"
)

# =============================================================================
# 2. PARÂMETROS DO PERÍODO
# =============================================================================
# Período idêntico ao do artigo base: ~15 anos de dados diários.
# DOW Inc. (ticker DOW) só existe como empresa independente desde abril/2019;
# por isso o período começa em 2019-06-01 para garantir que todas as 30
# ações tenham dados disponíveis.

data_inicio <- "2019-06-03"   # primeiro pregão completo do DOW Inc.
data_fim    <- "2024-12-31"

# =============================================================================
# 3. DOWNLOAD
# =============================================================================

cat("=== Baixando dados do Yahoo Finance ===\n")
cat("Período:", data_inicio, "a", data_fim, "\n\n")

y_list     <- list()   # retornos log diários (%)
preco_list <- list()   # preços ajustados (para referência)
falhas     <- character(0)

for (i in seq_along(tickers)) {
  tk <- tickers[i]
  nm <- nomes[i]
  cat(sprintf("  [%2d/30] %-6s (%s) ... ", i, tk, nm))

  resultado <- tryCatch({
    raw <- getSymbols(tk, src = "yahoo",
                      from = data_inicio, to = data_fim,
                      auto.assign = FALSE)

    # Interpola NAs internos (feriados de mercados diferentes)
    raw <- na.approx(raw, na.rm = FALSE)
    px  <- as.numeric(Ad(raw))          # preço ajustado por dividendos/splits
    px  <- as.numeric(na.approx(px, na.rm = FALSE))
    px  <- na.omit(px)

    # Retorno log diário em percentual (igual ao artigo base)
    ret <- 100 * diff(log(px))
    ret <- as.numeric(na.omit(ret))

    stopifnot(all(is.finite(ret)), length(ret) > 100)

    y_list[[nm]]     <<- ret
    preco_list[[nm]] <<- px
    cat("OK | T =", length(ret), "\n")
    TRUE
  }, error = function(e) {
    cat("ERRO:", conditionMessage(e), "\n")
    falhas <<- c(falhas, tk)
    FALSE
  })
}

# =============================================================================
# 4. ALINHAMENTO PELO COMPRIMENTO MÍNIMO
# =============================================================================
# Garante que todas as séries tenham o mesmo número de observações,
# condição necessária para o MGHMM.

if (length(y_list) > 0) {
  T_cada  <- sapply(y_list, length)
  T_min   <- min(T_cada)
  T_max   <- max(T_cada)

  cat("\n=== Comprimentos antes do alinhamento ===\n")
  print(T_cada)
  cat("Mínimo:", T_min, "| Máximo:", T_max, "\n")

  # Usa os T_min dias mais recentes de cada série
  y_list     <- lapply(y_list,     tail, T_min)
  preco_list <- lapply(preco_list, tail, T_min + 1)  # +1 pois px tem 1 obs a mais

  cat("Séries alinhadas em T =", T_min, "observações.\n")
}

# =============================================================================
# 5. SUMÁRIO RÁPIDO
# =============================================================================

cat("\n=== Sumário dos retornos ===\n")
cat(sprintf("%-20s  %7s  %7s  %7s\n", "Ação", "Média", "DP", "Min/Max"))
cat(strrep("-", 50), "\n")
for (nm in names(y_list)) {
  y <- y_list[[nm]]
  cat(sprintf("%-20s  %7.4f  %7.4f  %6.2f / %5.2f\n",
              nm, mean(y), sd(y), min(y), max(y)))
}

if (length(falhas) > 0) {
  cat("\n⚠ Tickers com falha no download:", paste(falhas, collapse = ", "), "\n")
} else {
  cat("\n✓ Todos os 30 tickers baixados com sucesso.\n")
}

# =============================================================================
# 6. SALVAR EM ARQUIVO .RData
# =============================================================================
# O arquivo gerado pode ser carregado no script principal com:
#   load("dados_djia.RData")

save(y_list, preco_list, tickers, nomes, data_inicio, data_fim,
     file = "dados_djia.RData")

cat("\n✓ Dados salvos em: dados_djia.RData\n")
cat("  Objetos disponíveis após load():\n")
cat("    y_list     – lista com os retornos log (%) de cada ação\n")
cat("    preco_list – lista com os preços ajustados de cada ação\n")
cat("    tickers    – vetor com os tickers do DJIA\n")
cat("    nomes      – vetor com os nomes completos das empresas\n")
cat("    data_inicio, data_fim – período da amostra\n")

# =============================================================================
# 7. COMO USAR NO SCRIPT PRINCIPAL
# =============================================================================
# Substitua o bloco de download dentro de run_example() por:
#
#   load("dados_djia.RData")
#
# e passe y_list diretamente para as funções de análise, por exemplo:
#
#   desc     <- plot_descriptive_stats(y_list)
#   adf      <- test_stationarity(y_list)
#   cor_ret  <- plot_return_correlation(y_list)
#   fit      <- fit_mghmm(y_list, S = 2, K = 3, n_init = 10)
