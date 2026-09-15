

source("00_common.R")

suppressPackageStartupMessages({
  library(forecast)
})

# --- Iperparametri -----------------------------------------------------------

K_YEAR   <- 10           
K_WEEK   <- 5            
P_YEAR   <- 24 * 365.25  
P_WEEK   <- 168          
ORDER    <- c(2, 0, 0)
SEASONAL <- list(order = c(0, 1, 1), period = 24)

FIT_MODEL2 <- FALSE

# =============================================================================
# Scelta parametri
# =============================================================================


lz   <- train$value_log
lz24 <- diff(lz, lag = 24)

save_plot("arima_00_acf_diff24.png", {
  acf(lz24, lag.max = 210,
      main = "ACF di diff(log(1+z), s=24) - scelta di q e Q")
}, width = 1200, height = 500)

save_plot("arima_01_pacf_diff24.png", {
  pacf(lz24, lag.max = 210,
       main = "PACF di diff(log(1+z), s=24) - scelta di p e P")
}, width = 1200, height = 500)


save_plot("arima_02_pacf_diff24_zoom.png", {
  pacf(lz24, lag.max = 30,
       main = "PACF di diff(log(1+z), s=24) - primi 30 lag")
}, width = 900, height = 500)



build_xreg <- function(d, t_index) {
  frt_y <- outer(t_index, seq_len(K_YEAR)) * 2 * pi / P_YEAR
  co_y  <- cos(frt_y); colnames(co_y) <- paste0("cos_year", seq_len(K_YEAR))
  si_y  <- sin(frt_y); colnames(si_y) <- paste0("sin_year", seq_len(K_YEAR))

  frt_w <- outer(t_index, seq_len(K_WEEK)) * 2 * pi / P_WEEK
  co_w  <- cos(frt_w); colnames(co_w) <- paste0("cos_week", seq_len(K_WEEK))
  si_w  <- sin(frt_w); colnames(si_w) <- paste0("sin_week", seq_len(K_WEEK))

  dow   <- factor(lubridate::wday(d$date, week_start = 1), levels = 1:7)
  X_dow <- model.matrix(~ dow - 1)[, -1, drop = FALSE]
  colnames(X_dow) <- paste0("dow", 2:7)

  cbind(co_y, si_y, co_w, si_w, X_dow,
        holiday = d$holiday, vacation = d$vacation)
}

n_train  <- nrow(train)
X_train  <- build_xreg(train,  seq_len(n_train))
X_future <- build_xreg(future, n_train + seq_len(H))
stopifnot(identical(colnames(X_train), colnames(X_future)))


rk <- qr(diff(X_train, lag = 24))$rank
if (rk < ncol(X_train)) {
  stop(sprintf("xreg singolare dopo diff(s=24): rango %d su %d colonne",
               rk, ncol(X_train)))
}


y <- train$value + 1
stopifnot(all(y > 0))

# =============================================================================
#  Valutazione 
# =============================================================================

library(forecast)
fit_idx <- setdiff(seq_len(n_train), holdout_idx)

cat("\n--- Stima sul training ridotto (test = ultimo mese) ---\n")
mod_cv <- Arima(
  y[fit_idx],
  order            = ORDER,
  seasonal         = SEASONAL,
  include.constant = FALSE,   
  method           = "CSS",   
  xreg             = X_train[fit_idx, ],
  lambda           = 0
)

fc_cv     <- forecast(mod_cv, h = length(holdout_idx),
                      xreg = X_train[holdout_idx, ])
pred_cv   <- postprocess_counts(as.numeric(fc_cv$mean) - 1)
actual_cv <- train$value[holdout_idx]

metrics_holdout <- eval_forecast(actual_cv, pred_cv, train$value[fit_idx])
cat("\nMetriche out-of-sample (ultimo mese):\n")
print(metrics_holdout)


save_plot("arima_03_test_fitted_vs_original.png", {
  plot_fit_vs_actual(actual_cv, pred_cv,
                     main = "ARIMA - osservato vs previsto sul test set (ultimo mese)",
                     xlab = "Ore", ylab = "Conteggio persone")
}, width = 1400, height = 500)

save_plot("arima_04_test_zoom_settimana.png", {
  idx <- 1:(24 * 7)
  plot_fit_vs_actual(actual_cv[idx], pred_cv[idx],
                     main = "ARIMA - prima settimana del test set",
                     xlab = "Ore", ylab = "Conteggio persone")
}, width = 1200, height = 450)

save_plot("arima_05_test_errori.png", {
  err <- actual_cv - pred_cv
  plot(err, type = "h", col = adjustcolor(COL_PRED, alpha.f = 0.7),
       main = "ARIMA - errori di previsione sul test set",
       xlab = "Ore", ylab = "Osservato - Previsto")
  abline(h = 0, col = "grey30")
}, width = 1400, height = 400)

# =============================================================================
# Stima finale su tutto lo storico
# =============================================================================

cat("\n--- Stima finale su tutto lo storico ---\n")
mod <- Arima(
  y,
  order            = ORDER,
  seasonal         = SEASONAL,
  include.constant = FALSE,
  method           = "CSS",
  xreg             = X_train,
  lambda           = 0
)

print(summary(mod))


cat("\nMetriche in-sample (accuracy sul training):\n")
print(accuracy(mod))

# --- Diagnostica dei residui -------------------------------------------------


res <- residuals(mod)

save_plot("arima_06_residui_acf.png", {
  acf(res, lag.max = 500, main = "ARIMA - ACF dei residui")
}, width = 1400, height = 500)

save_plot("arima_07_residui_pacf.png", {
  pacf(res, lag.max = 500, main = "ARIMA - PACF dei residui")
}, width = 1400, height = 500)

save_plot("arima_08_residui_acf_zoom.png", {
  acf(res, lag.max = 200, main = "ARIMA - ACF dei residui (primi 200 lag)")
}, width = 1200, height = 500)

save_plot("arima_09_residui_serie.png", {
  plot(as.numeric(res), type = "l", col = adjustcolor(COL_ACTUAL, alpha.f = 0.6),
       main = "ARIMA - residui nel tempo", xlab = "Ore", ylab = "Residuo")
  abline(h = 0, col = COL_PRED)
}, width = 1400, height = 400)

save_plot("arima_10_residui_istogramma.png", {
  hist(as.numeric(res), breaks = 80, col = COL_ACTUAL, border = "white",
       main = "ARIMA - distribuzione dei residui", xlab = "Residuo")
}, width = 900, height = 500)

save_plot("arima_11_residui_qq.png", {
  qqnorm(as.numeric(res), main = "ARIMA - QQ plot dei residui",
         col = adjustcolor(COL_ACTUAL, alpha.f = 0.4), pch = 16)
  qqline(as.numeric(res), col = COL_PRED, lwd = 2)
}, width = 800, height = 600)

# Test di Ljung-Box: versione numerica di quello che l'ACF mostra a occhio.
cat("\nLjung-Box sui residui (lag 48 e 336):\n")
print(Box.test(res, lag = 48,  type = "Ljung-Box"))
print(Box.test(res, lag = 336, type = "Ljung-Box"))

# Fitted vs original in-sample sulle ultime due settimane dello storico.
save_plot("arima_12_insample_fitted.png", {
  idx <- (n_train - 24 * 14 + 1):n_train
  plot_fit_vs_actual(train$value[idx], as.numeric(fitted(mod))[idx] - 1,
                     main = "ARIMA - fitted vs original (ultime 2 settimane di training)",
                     xlab = "Ore", ylab = "Conteggio persone")
}, width = 1200, height = 450)

# =============================================================================
#  Previsione
# =============================================================================

fc <- forecast(mod, h = H, xreg = X_future)
pred_arima <- postprocess_counts(as.numeric(fc$mean) - 1)
stopifnot(length(pred_arima) == H)

save_plot("arima_13_previsione.png", {
  plot_forecast_tail(train$value, pred_arima,
                     main = "ARIMA - ultime 2 settimane osservate e previsione")
}, width = 1400, height = 500)

save_plot("arima_14_previsione_completa.png", {
  plot(pred_arima, type = "l", col = COL_PRED, lwd = 1.2,
       main = sprintf("ARIMA - previsione sulle %d ore richieste", H),
       xlab = "Ore dall'inizio della finestra", ylab = "Conteggio persone")
}, width = 1400, height = 450)


save_plot("arima_15_previsione_intervalli.png", {
  plot(fc, include = 24 * 14,
       main = "ARIMA - previsione con intervalli all'80% e 95%",
       xlab = "Ore", ylab = "value + 1")
}, width = 1400, height = 500)

# --- Salvataggio -------------------------------------------------------------

out <- data.frame(time = future$time, Arima = pred_arima)
write.csv(out, file.path(OUTPUT_DIR, "pred_arima.csv"), row.names = FALSE)

cat(sprintf("\nSalvate %d previsioni in %s\n",
            nrow(out), file.path(OUTPUT_DIR, "pred_arima.csv")))
print(summary(pred_arima))

# =============================================================================
# Modello 2 
# =============================================================================

if (FIT_MODEL2) {
  cat("\n--- Modello 2 (confronto) ---\n")

  lz168 <- diff(lz, lag = 168)

  save_plot("arima_m2_00_acf_diff168.png", {
    acf(lz168, lag.max = 1000, main = "ACF di diff(log(1+z), s=168)")
  }, width = 1400, height = 500)

  save_plot("arima_m2_01_pacf_diff168.png", {
    pacf(lz168, lag.max = 1000, main = "PACF di diff(log(1+z), s=168)")
  }, width = 1400, height = 500)

  K_YEAR2 <- 15
  build_xreg2 <- function(d, t_index) {
    frt <- outer(t_index, seq_len(K_YEAR2)) * 2 * pi / P_YEAR
    co  <- cos(frt); colnames(co) <- paste0("cos_year", seq_len(K_YEAR2))
    si  <- sin(frt); colnames(si) <- paste0("sin_year", seq_len(K_YEAR2))
    cbind(co, si, holiday = d$holiday, vacation = d$vacation)
  }

  X2_train  <- build_xreg2(train,  seq_len(n_train))
  X2_future <- build_xreg2(future, n_train + seq_len(H))

  mod2 <- Arima(
    y[fit_idx],
    order            = c(24, 0, 0),
    seasonal         = list(order = c(0, 1, 1), period = 168),
    include.constant = FALSE,
    method           = "CSS",
    xreg             = X2_train[fit_idx, ],
    lambda           = 0
  )
  print(summary(mod2))

  res2 <- residuals(mod2)
  save_plot("arima_m2_02_residui_acf.png", {
    acf(res2, lag.max = 500, main = "Modello 2 - ACF dei residui")
  }, width = 1400, height = 500)

  save_plot("arima_m2_03_residui_pacf.png", {
    pacf(res2, lag.max = 500, main = "Modello 2 - PACF dei residui")
  }, width = 1400, height = 500)

  fc2_cv   <- forecast(mod2, h = length(holdout_idx),
                       xreg = X2_train[holdout_idx, ])
  pred2_cv <- postprocess_counts(as.numeric(fc2_cv$mean) - 1)

  save_plot("arima_m2_04_test_fitted_vs_original.png", {
    plot_fit_vs_actual(actual_cv, pred2_cv,
                       main = "Modello 2 - osservato vs previsto sul test set",
                       xlab = "Ore", ylab = "Conteggio persone")
  }, width = 1400, height = 500)

  cat("\nConfronto out-of-sample sul test set:\n")
  print(rbind(
    "Modello 1 (2,0,0)(0,1,1)[24]"   = metrics_holdout,
    "Modello 2 (24,0,0)(0,1,1)[168]" = eval_forecast(actual_cv, pred2_cv,
                                                     train$value[fit_idx])
  ))
}

cat(sprintf("\nGrafici salvati in %s\n", PLOTS_DIR))

