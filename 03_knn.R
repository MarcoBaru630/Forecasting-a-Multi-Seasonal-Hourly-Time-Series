# =============================================================================
# 03_knn.R
# Modello di Machine Learning: KNN regression.
#
# Impostazione del report (sezione 2.5):
#   - target = differenza stagionale a 24 ore su scala log, cioe' la variazione
#     rispetto alla stessa ora del giorno precedente. Cosi' la stagionalita'
#     giornaliera, che e' la componente dominante, e' gia' tolta dal target e il
#     modello non deve impararla dalle distanze.
#   - feature: lag 168 (stagionalita' settimanale) + regressori ciclici sin/cos
#     su ora, giorno della settimana e giorno dell'anno.
#   - la stagionalita' annuale NON usa un lag a 8760 ore: costerebbe un anno
#     intero di dati e aggiungerebbe rumore. Viene catturata dai ciclici sul
#     giorno dell'anno.
#   - standardizzazione stimata solo sul training e riapplicata a test e
#     previsione con gli stessi centri e scale.
#
# Output: output/pred_knn.csv, output/plots/knn_*.png
# =============================================================================

source("00_common.R")

suppressPackageStartupMessages({
  library(FNN)
})

K_GRID   <- 1:20
LAG_DAY  <- 24
LAG_WEEK <- 168

# Valuta anche il test set in modo ricorsivo (previsione a 720 passi con i lag
# ricostruiti dalle previsioni stesse). E' molto piu' lento ma e' l'unico
# confronto onesto con la previsione finale; a FALSE resta solo la valutazione
# a un passo.
RUN_RECURSIVE_HOLDOUT <- TRUE

# =============================================================================
# 1. Timeline completa e feature
# =============================================================================
# Storico e finestra da prevedere in un unico data frame ordinato: i lag a 24 e
# 168 ore delle prime ore da prevedere pescano dagli ultimi valori osservati,
# quindi le due parti devono stare nella stessa serie.

all_df <- rbind(
  train[,  c("time", "date", "hour", "value")],
  future[, c("time", "date", "hour", "value")]
)
all_df <- all_df[order(all_df$time), ]
rownames(all_df) <- NULL

all_df$value_log <- ifelse(is.na(all_df$value), NA_real_, log(all_df$value + 1))

all_df$hour_n  <- lubridate::hour(all_df$time)
all_df$weekday <- lubridate::wday(all_df$time, week_start = 1)
all_df$doy     <- lubridate::yday(all_df$time)

all_df$sin_h <- sin(2 * pi * all_df$hour_n  / 24)
all_df$cos_h <- cos(2 * pi * all_df$hour_n  / 24)
all_df$sin_w <- sin(2 * pi * all_df$weekday / 7)
all_df$cos_w <- cos(2 * pi * all_df$weekday / 7)
all_df$sin_y <- sin(2 * pi * all_df$doy     / 365.25)
all_df$cos_y <- cos(2 * pi * all_df$doy     / 365.25)

all_df$lag24  <- lag_vec(all_df$value_log, LAG_DAY)
all_df$lag168 <- lag_vec(all_df$value_log, LAG_WEEK)
all_df$target <- all_df$value_log - all_df$lag24

# Ordine delle colonne definito una volta sola e riusato ovunque. Nella versione
# originale la matrice di training veniva costruita con select(-...) mentre le
# righe di previsione erano assemblate a mano in una funzione separata: bastava
# aggiungere una colonna da una parte per disallineare silenziosamente le due
# matrici e mandare in vacca le distanze.
FEATURE_COLS <- c("hour_n", "weekday", "doy",
                  "sin_h", "cos_h", "sin_w", "cos_w", "sin_y", "cos_y",
                  "lag24", "lag168")

make_features <- function(d) as.matrix(d[, FEATURE_COLS, drop = FALSE])

# =============================================================================
# 2. Split training / test
# =============================================================================

n_train   <- nrow(train)
test_rows <- holdout_idx
fit_rows  <- setdiff(seq_len(n_train), test_rows)

train_knn <- all_df[fit_rows, ]
train_knn <- train_knn[complete.cases(train_knn[, c(FEATURE_COLS, "target")]), ]

test_knn <- all_df[test_rows, ]
test_knn <- test_knn[complete.cases(test_knn[, c(FEATURE_COLS, "target")]), ]

cat(sprintf("Righe di training utilizzabili: %d\n", nrow(train_knn)))
cat(sprintf("Righe di test utilizzabili:     %d\n", nrow(test_knn)))

X_train <- scale(make_features(train_knn))
center  <- attr(X_train, "scaled:center")
scale_v <- attr(X_train, "scaled:scale")
y_train <- train_knn$target

apply_scaling <- function(X) scale(X, center = center, scale = scale_v)
X_test <- apply_scaling(make_features(test_knn))

# Distribuzione del target: mostra perche' la differenza stagionale e' una buona
# variabile da modellare (centrata su zero, senza trend residuo).
save_plot("knn_01_target_distribuzione.png", {
  par(mfrow = c(1, 2))
  hist(y_train, breaks = 80, col = COL_ACTUAL, border = "white",
       main = "Target: diff. stagionale 24h (scala log)", xlab = "Target")
  acf(y_train, lag.max = 200, main = "ACF del target")
  par(mfrow = c(1, 1))
}, width = 1300, height = 450)

# =============================================================================
# 3. Selezione di k
# =============================================================================
# Valutazione a un passo: per ogni ora del test si usa il lag24 vero, che nella
# previsione reale non sarebbe noto. E' la procedura del report (Tabella 6) e
# viene mantenuta per riprodurre quei numeri, ma va letta come un limite
# superiore ottimistico rispetto alla previsione ricorsiva.

results <- data.frame(k = K_GRID, MAE = NA_real_, MSE = NA_real_, R2 = NA_real_)

for (i in seq_along(K_GRID)) {
  k <- K_GRID[i]
  pred_diff <- FNN::knn.reg(train = X_train, test = X_test,
                            y = y_train, k = k)$pred
  pred <- exp(pred_diff + test_knn$lag24) - 1
  err  <- test_knn$value - pred

  results$MAE[i] <- mean(abs(err))
  results$MSE[i] <- mean(err^2)
  results$R2[i]  <- 1 - sum(err^2) / sum((test_knn$value - mean(test_knn$value))^2)
}

cat("\nSelezione di k (valutazione a un passo):\n")
print(results, row.names = FALSE)

best_k <- results$k[which.min(results$MAE)]
cat(sprintf("\nk selezionato (MAE minimo): %d\n", best_k))

# Nota: lo script originale selezionava k = 7 sulla griglia ma poi lanciava la
# previsione finale con k = 15. Qui best_k viene usato ovunque.

# Curva di selezione: e' il grafico che rende leggibile la Tabella 6 del report.
# La forma a U e' il trade-off bias-varianza: k piccolo sovradatta, k grande
# media su vicini troppo diversi.
save_plot("knn_02_selezione_k.png", {
  par(mfrow = c(1, 3))
  plot(results$k, results$MAE, type = "b", pch = 16, col = COL_ACTUAL, lwd = 2,
       main = "MAE al variare di k", xlab = "k", ylab = "MAE")
  points(best_k, min(results$MAE), col = COL_PRED, pch = 16, cex = 2)
  plot(results$k, results$MSE, type = "b", pch = 16, col = COL_ACTUAL, lwd = 2,
       main = "MSE al variare di k", xlab = "k", ylab = "MSE")
  plot(results$k, results$R2, type = "b", pch = 16, col = COL_ACTUAL, lwd = 2,
       main = "R2 al variare di k", xlab = "k", ylab = "R2")
  par(mfrow = c(1, 1))
}, width = 1400, height = 450)

# Fitted vs original a un passo, con il k scelto.
pred_1step <- postprocess_counts(
  exp(FNN::knn.reg(train = X_train, test = X_test,
                   y = y_train, k = best_k)$pred + test_knn$lag24) - 1
)

save_plot("knn_03_test_1passo.png", {
  plot_fit_vs_actual(test_knn$value, pred_1step,
                     main = sprintf("KNN (k=%d) - previsione a un passo sul test set", best_k),
                     xlab = "Ore", ylab = "Conteggio persone")
}, width = 1400, height = 500)

save_plot("knn_04_test_1passo_zoom.png", {
  idx <- 1:(24 * 7)
  plot_fit_vs_actual(test_knn$value[idx], pred_1step[idx],
                     main = "KNN - prima settimana del test set (un passo)",
                     xlab = "Ore", ylab = "Conteggio persone")
}, width = 1200, height = 450)

save_plot("knn_05_scatter.png", {
  plot(test_knn$value, pred_1step, pch = 16,
       col = adjustcolor(COL_ACTUAL, alpha.f = 0.25),
       main = "KNN - previsto vs osservato (test set)",
       xlab = "Osservato", ylab = "Previsto")
  abline(0, 1, col = COL_PRED, lwd = 2)
}, width = 700, height = 700)

# =============================================================================
# 4. Previsione ricorsiva
# =============================================================================
# Le ore da prevedere sono consecutive: dalla 25esima in poi il lag a 24 ore cade
# dentro la finestra da prevedere e non e' osservato. Va quindi ricostruito
# dalle previsioni gia' fatte, un passo alla volta.

forecast_recursive <- function(series_log, idx_to_predict, k) {
  s <- series_log
  for (i in idx_to_predict) {
    lag24  <- s[i - LAG_DAY]
    lag168 <- s[i - LAG_WEEK]
    if (is.na(lag24) || is.na(lag168)) next

    x_row <- data.frame(
      hour_n  = all_df$hour_n[i],
      weekday = all_df$weekday[i],
      doy     = all_df$doy[i],
      sin_h   = all_df$sin_h[i], cos_h = all_df$cos_h[i],
      sin_w   = all_df$sin_w[i], cos_w = all_df$cos_w[i],
      sin_y   = all_df$sin_y[i], cos_y = all_df$cos_y[i],
      lag24   = lag24,
      lag168  = lag168
    )

    x_scaled  <- apply_scaling(make_features(x_row))
    pred_diff <- FNN::knn.reg(train = X_train, test = x_scaled,
                              y = y_train, k = k)$pred
    s[i] <- pred_diff + lag24
  }
  s
}

if (RUN_RECURSIVE_HOLDOUT) {
  cat("\nValutazione ricorsiva sul test set (piu' lenta)...\n")
  s_cv <- all_df$value_log
  s_cv[test_rows] <- NA_real_
  s_cv <- forecast_recursive(s_cv, test_rows, best_k)

  pred_cv <- postprocess_counts(exp(s_cv[test_rows]) - 1)
  metrics_holdout <- eval_forecast(train$value[test_rows], pred_cv,
                                   train$value[fit_rows])
  cat("\nMetriche out-of-sample ricorsive (ultimo mese):\n")
  print(metrics_holdout)

  save_plot("knn_06_test_ricorsivo.png", {
    plot_fit_vs_actual(train$value[test_rows], pred_cv,
                       main = sprintf("KNN (k=%d) - previsione ricorsiva sul test set", best_k),
                       xlab = "Ore", ylab = "Conteggio persone")
  }, width = 1400, height = 500)

  # Degrado dell'errore con l'orizzonte: e' la differenza fra la valutazione a
  # un passo e quella ricorsiva, e spiega perche' la Tabella 6 e' ottimistica.
  save_plot("knn_07_errore_per_orizzonte.png", {
    err_abs <- abs(train$value[test_rows] - pred_cv)
    giorno  <- rep(1:ceiling(length(err_abs) / 24), each = 24)[seq_along(err_abs)]
    mae_g   <- tapply(err_abs, giorno, mean)
    plot(as.integer(names(mae_g)), mae_g, type = "b", pch = 16,
         col = COL_ACTUAL, lwd = 2,
         main = "KNN - MAE per giorno di orizzonte (previsione ricorsiva)",
         xlab = "Giorni avanti", ylab = "MAE")
  }, width = 1100, height = 450)
}

cat("\nPrevisione ricorsiva sulla finestra richiesta...\n")
future_rows <- (n_train + 1):nrow(all_df)
s_fut <- forecast_recursive(all_df$value_log, future_rows, best_k)

pred_knn <- postprocess_counts(exp(s_fut[future_rows]) - 1)
stopifnot(length(pred_knn) == H)

save_plot("knn_08_previsione.png", {
  plot_forecast_tail(train$value, pred_knn,
                     main = "KNN - ultime 2 settimane osservate e previsione")
}, width = 1400, height = 500)

save_plot("knn_09_previsione_completa.png", {
  plot(pred_knn, type = "l", col = COL_PRED, lwd = 1.2,
       main = sprintf("KNN - previsione sulle %d ore richieste", H),
       xlab = "Ore dall'inizio della finestra", ylab = "Conteggio persone")
}, width = 1400, height = 450)

# --- Salvataggio -------------------------------------------------------------

out <- data.frame(time = future$time, KNN = pred_knn)
write.csv(out, file.path(OUTPUT_DIR, "pred_knn.csv"), row.names = FALSE)

cat(sprintf("\nSalvate %d previsioni in %s\n",
            nrow(out), file.path(OUTPUT_DIR, "pred_knn.csv")))
print(summary(pred_knn))
cat(sprintf("\nGrafici salvati in %s\n", PLOTS_DIR))
