
source("00_common.R")

library(KFAS)


N_HARMONICS <- 16 

PLOT_HOUR <- 18

fit_one_hour <- function(h) {

  hist_h <- train[train$hour == h, c("date", "value", "holiday", "vacation")]
  fut_h  <- future[future$hour == h, c("date", "value", "holiday", "vacation")]

  wdf <- rbind(hist_h, fut_h)
  wdf <- wdf[order(wdf$date), ]
  rownames(wdf) <- NULL

  test_dates <- unique(train$date[holdout_idx])
  wdf$y_na <- wdf$value
  wdf$y_na[wdf$date %in% test_dates] <- NA_real_
  wdf$y_na[is.na(wdf$value)]         <- NA_real_

  vy <- var(wdf$y_na, na.rm = TRUE)

  mod <- SSModel(
    y_na ~ holiday + vacation +
      SSMtrend(2, list(NA, NA)) +
      SSMseasonal(7,   NA, "dummy") +
      SSMseasonal(365, NA, "trig", harmonics = 1:N_HARMONICS),
    H    = NA,
    data = wdf
  )


  mod$a1["level", 1] <- mean(wdf$y_na, na.rm = TRUE)
  diag(mod$P1)       <- vy
  diag(mod$P1inf)    <- 0

  updt <- function(pars, model) {
    model$Q[1, 1, 1] <- exp(pars[1])                 # varianza livello
    model$Q[2, 2, 1] <- exp(pars[2])                 # varianza slope
    model$Q[3, 3, 1] <- exp(pars[3])                 # varianza dummy s=7
    diag(model$Q[-(1:3), -(1:3), 1]) <- exp(pars[4]) # varianza trig s=365
    model$H[1, 1, 1] <- exp(pars[5])                 # varianza osservazione
    model
  }

  inits <- log(c(vy / 20, vy / 10000, vy / 10000, vy / 100, vy / 10))

  fit <- fitSSM(mod, inits, updt, method = "BFGS")

  if (fit$optim.out$convergence != 0) {
    warning(sprintf("Ora %d: convergenza non raggiunta (code %d)",
                    h, fit$optim.out$convergence))
  }

  kfs <- KFS(fit$model, smoothing = c("state", "signal"))
  wdf$pred <- as.numeric(kfs$muhat)

  list(hour = h, data = wdf, kfs = kfs,
       convergence = fit$optim.out$convergence,
       variances = exp(fit$optim.out$par))
}

# =============================================================================
# Stima delle 24 fasce orarie
# =============================================================================

results <- vector("list", 24)

for (h in 1:24) {
  t0 <- Sys.time()
  results[[h]] <- tryCatch(fit_one_hour(h), error = function(e) {
    warning(sprintf("Ora %d: stima fallita (%s)", h, conditionMessage(e)))
    NULL
  })
  cat(sprintf("Ora %2d/24 completata in %5.1f s%s\n", h,
              as.numeric(difftime(Sys.time(), t0, units = "secs")),
              if (is.null(results[[h]])) "  [FALLITA]" else ""))
}

failed <- which(vapply(results, is.null, logical(1)))
if (length(failed) > 0) {
  stop("Stima fallita per le ore: ", paste(failed, collapse = ", "))
}

cat("\nCodici di convergenza (0 = ok):\n")
print(setNames(vapply(results, function(r) r$convergence, numeric(1)), 1:24))


var_tab <- do.call(rbind, lapply(results, function(r) r$variances))
colnames(var_tab) <- c("level", "slope", "seas7", "seas365", "epsilon")
rownames(var_tab) <- paste0("h", 1:24)
cat("\nVarianze stimate per fascia oraria:\n")
print(round(var_tab, 4))

# =============================================================================
# Scomposizione nelle componenti 
# =============================================================================

r   <- results[[PLOT_HOUR]]
wdf <- r$data
ah  <- r$kfs$alphahat

save_plot(sprintf("ucm_01_livello_h%d.png", PLOT_HOUR), {
  plot(wdf$date, wdf$value, type = "l",
       col = adjustcolor(COL_ACTUAL, alpha.f = 0.55),
       main = sprintf("UCM ora %d - serie osservata e livello smoothed", PLOT_HOUR),
       xlab = "Data", ylab = "Conteggio persone")
  lines(wdf$date, ah[, "level"], col = COL_PRED, lwd = 2)
  legend("topright", legend = c("Osservato", "Livello"),
         col = c(COL_ACTUAL, COL_PRED), lwd = c(1, 2), bty = "n")
}, width = 1400, height = 500)

save_plot(sprintf("ucm_02_slope_h%d.png", PLOT_HOUR), {
  plot(wdf$date, ah[, "slope"], type = "l", col = COL_PRED, lwd = 1.5,
       main = sprintf("UCM ora %d - componente slope", PLOT_HOUR),
       xlab = "Data", ylab = "Slope")
  abline(h = 0, col = "grey40", lty = 2)
}, width = 1400, height = 400)


seas7_name <- grep("^sea_dummy1$", colnames(ah), value = TRUE)
if (length(seas7_name) == 1) {
  save_plot(sprintf("ucm_03_stagionalita_settimanale_h%d.png", PLOT_HOUR), {
    idx <- which(wdf$date >= as.Date("2019-02-01") &
                   wdf$date <= as.Date("2019-04-28"))
    plot(wdf$date[idx], ah[idx, seas7_name], type = "l",
         col = COL_PRED, lwd = 1.6,
         main = sprintf("UCM ora %d - stagionalita' settimanale (dummy)", PLOT_HOUR),
         xlab = "Data", ylab = "Effetto")
    abline(h = 0, col = "grey40", lty = 2)
  }, width = 1200, height = 450)
}

trig_names <- grep("^sea_trig[0-9]+$", colnames(ah), value = TRUE)
if (length(trig_names) > 0) {
  seas365 <- rowSums(ah[, trig_names, drop = FALSE])

  save_plot(sprintf("ucm_04_stagionalita_annuale_h%d.png", PLOT_HOUR), {
    plot(wdf$date, seas365, type = "l", col = COL_PRED, lwd = 1.5,
         main = sprintf("UCM ora %d - stagionalita' annuale (%d armoniche)",
                        PLOT_HOUR, length(trig_names) / 2),
         xlab = "Data", ylab = "Effetto")
    abline(h = 0, col = "grey40", lty = 2)
  }, width = 1400, height = 450)

  
  has_cal <- all(c("holiday", "vacation") %in% colnames(ah))
  if (has_cal) {
    cal_effect <- ah[, "holiday"] * wdf$holiday + ah[, "vacation"] * wdf$vacation
    save_plot(sprintf("ucm_05_annuale_e_calendario_h%d.png", PLOT_HOUR), {
      plot(wdf$date, seas365 + cal_effect, type = "l", col = COL_PRED, lwd = 1.4,
           main = sprintf("UCM ora %d - stagionalita' annuale + effetti calendario",
                          PLOT_HOUR),
           xlab = "Data", ylab = "Effetto")
      abline(h = 0, col = "grey40", lty = 2)
    }, width = 1400, height = 450)

    cat(sprintf("\nEffetti di calendario stimati (ora %d): holiday %.2f, vacation %.2f\n",
                PLOT_HOUR, tail(ah[, "holiday"], 1), tail(ah[, "vacation"], 1)))
  }
}


save_plot(sprintf("ucm_06_segnale_h%d.png", PLOT_HOUR), {
  idx <- which(wdf$date >= as.Date("2019-01-01"))
  plot(wdf$date[idx], wdf$value[idx], type = "l",
       col = adjustcolor(COL_ACTUAL, alpha.f = 0.6),
       main = sprintf("UCM ora %d - osservato vs segnale stimato (2019)", PLOT_HOUR),
       xlab = "Data", ylab = "Conteggio persone")
  lines(wdf$date[idx], wdf$pred[idx], col = COL_PRED, lwd = 1.6)
  legend("topright", legend = c("Osservato", "Segnale"),
         col = c(COL_ACTUAL, COL_PRED), lwd = c(1, 1.6), bty = "n")
}, width = 1400, height = 500)


save_plot(sprintf("ucm_07_residui_h%d.png", PLOT_HOUR), {
  resid_h <- wdf$value - wdf$pred
  par(mfrow = c(1, 2))
  plot(wdf$date, resid_h, type = "h",
       col = adjustcolor(COL_PRED, alpha.f = 0.6),
       main = sprintf("UCM ora %d - residui", PLOT_HOUR),
       xlab = "Data", ylab = "Residuo")
  abline(h = 0, col = "grey30")
  acf(resid_h, na.action = na.pass, lag.max = 60,
      main = "ACF dei residui")
  par(mfrow = c(1, 1))
}, width = 1400, height = 450)

# =============================================================================
# Ricomposizione della serie oraria
# =============================================================================

all_pred <- do.call(rbind, lapply(results, function(r) {
  data.frame(date = r$data$date, hour = r$hour, pred = r$data$pred)
}))

# --- Valutazione sul test set ------------------------------------------------

test_actual <- train[holdout_idx, c("date", "hour", "value")]
test_eval   <- merge(test_actual, all_pred, by = c("date", "hour"))
test_eval   <- test_eval[order(test_eval$date, test_eval$hour), ]

pred_test <- postprocess_counts(test_eval$pred)
fit_rows  <- setdiff(seq_len(nrow(train)), holdout_idx)

metrics_holdout <- eval_forecast(test_eval$value, pred_test, train$value[fit_rows])
cat("\nMetriche out-of-sample (ultimo mese):\n")
print(metrics_holdout)

save_plot("ucm_08_test_fitted_vs_original.png", {
  plot_fit_vs_actual(test_eval$value, pred_test,
                     main = "UCM - osservato vs previsto sul test set (ultimo mese)",
                     xlab = "Ore", ylab = "Conteggio persone")
}, width = 1400, height = 500)

save_plot("ucm_09_test_zoom_settimana.png", {
  idx <- 1:(24 * 7)
  plot_fit_vs_actual(test_eval$value[idx], pred_test[idx],
                     main = "UCM - prima settimana del test set",
                     xlab = "Ore", ylab = "Conteggio persone")
}, width = 1200, height = 450)

save_plot("ucm_10_test_errori.png", {
  err <- test_eval$value - pred_test
  plot(err, type = "h", col = adjustcolor(COL_PRED, alpha.f = 0.7),
       main = "UCM - errori di previsione sul test set",
       xlab = "Ore", ylab = "Osservato - Previsto")
  abline(h = 0, col = "grey30")
}, width = 1400, height = 400)

# =============================================================================
#  Previsione
# =============================================================================

fut_key <- future[, c("time", "date", "hour")]
fut_key <- merge(fut_key, all_pred, by = c("date", "hour"), all.x = TRUE)
fut_key <- fut_key[order(fut_key$time), ]

stopifnot(nrow(fut_key) == H)
stopifnot(!any(is.na(fut_key$pred)))

pred_ucm <- postprocess_counts(fut_key$pred)

save_plot("ucm_11_previsione.png", {
  plot_forecast_tail(train$value, pred_ucm,
                     main = "UCM - ultime 2 settimane osservate e previsione")
}, width = 1400, height = 500)

save_plot("ucm_12_previsione_completa.png", {
  plot(pred_ucm, type = "l", col = COL_PRED, lwd = 1.2,
       main = sprintf("UCM - previsione sulle %d ore richieste", H),
       xlab = "Ore dall'inizio della finestra", ylab = "Conteggio persone")
}, width = 1400, height = 450)

# --- Salvataggio -------------------------------------------------------------

out <- data.frame(time = fut_key$time, UCM = pred_ucm)
write.csv(out, file.path(OUTPUT_DIR, "pred_ucm.csv"), row.names = FALSE)

cat(sprintf("\nSalvate %d previsioni in %s\n",
            nrow(out), file.path(OUTPUT_DIR, "pred_ucm.csv")))
print(summary(pred_ucm))
cat(sprintf("\nGrafici salvati in %s\n", PLOTS_DIR))
