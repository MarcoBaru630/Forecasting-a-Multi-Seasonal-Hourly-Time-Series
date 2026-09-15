
source("00_common.R")

suppressPackageStartupMessages({
  library(tseries)
})

z  <- train$value
lz <- train$value_log

# =============================================================================
#  Andamento della serie
# =============================================================================

save_plot("eda_01_serie_intera.png", {
  plot(train$time, z, type = "l", col = COL_ACTUAL,
       main = "Serie oraria completa",
       xlab = "Tempo", ylab = "Conteggio persone")
}, width = 1400, height = 500)

save_plot("eda_02_zoom_3_settimane.png", {
  idx <- 1:(24 * 21)
  plot(train$time[idx], z[idx], type = "l", col = COL_ACTUAL,
       main = "Prime tre settimane - ciclo giornaliero e settimanale",
       xlab = "Tempo", ylab = "Conteggio persone")
}, width = 1400, height = 450)

# =============================================================================
#  Eteroschedasticita': mean-sd plot
# =============================================================================

nday  <- rep(1:ceiling(length(z) / 24), each = 24)[seq_along(z)]
means <- tapply(z, nday, mean)
sds   <- tapply(z, nday, sd)

save_plot("eda_03_mean_sd.png", {
  plot(means, sds, pch = 16, col = adjustcolor(COL_ACTUAL, alpha.f = 0.4),
       main = "Mean-sd plot giornaliero (serie originale)",
       xlab = "Media giornaliera", ylab = "Deviazione standard giornaliera")
  abline(lm(sds ~ means), col = COL_PRED, lwd = 2)
}, width = 800, height = 600)

means_l <- tapply(lz, nday, mean)
sds_l   <- tapply(lz, nday, sd)

save_plot("eda_04_mean_sd_log.png", {
  plot(means_l, sds_l, pch = 16, col = adjustcolor(COL_ACTUAL, alpha.f = 0.4),
       main = "Mean-sd plot giornaliero (dopo log(1+z))",
       xlab = "Media giornaliera", ylab = "Deviazione standard giornaliera")
  abline(lm(sds_l ~ means_l), col = COL_PRED, lwd = 2)
}, width = 800, height = 600)

cat(sprintf("\nPendenza mean-sd: originale %.4f -> log %.4f\n",
            coef(lm(sds ~ means))[2], coef(lm(sds_l ~ means_l))[2]))

# =============================================================================
# Serie trasformata
# =============================================================================

save_plot("eda_05_serie_log.png", {
  plot(train$time, lz, type = "l", col = COL_ACTUAL,
       main = "Serie dopo la trasformazione log(1 + z)",
       xlab = "Tempo", ylab = "log(1 + value)")
}, width = 1400, height = 500)

# =============================================================================
# Test di stazionarieta'
# =============================================================================


cat("\n--- Test di stazionarieta' sulla serie originale ---\n")
print(suppressWarnings(kpss.test(z, null = "Level")))
print(suppressWarnings(kpss.test(z, null = "Trend")))
print(suppressWarnings(adf.test(z)))

cat("\n--- Test di stazionarieta' sulla serie log ---\n")
print(suppressWarnings(kpss.test(lz, null = "Level")))
print(suppressWarnings(kpss.test(lz, null = "Trend")))

# =============================================================================
# stagionalita'
# =============================================================================


lz24  <- diff(lz, lag = 24)
lz168 <- diff(lz, lag = 168)


save_plot("eda_06_acf_s24.png", {
  acf(lz, lag.max = 72, main = "ACF log(1+z): picchi ogni 24 lag (giornaliera)")
})

save_plot("eda_07_pacf_s24.png", {
  pacf(lz, lag.max = 72, main = "PACF log(1+z)")
})


save_plot("eda_08_acf_s168.png", {
  acf(lz24, lag.max = 504,
      main = "ACF dopo diff(s=24): picchi ogni 168 lag (settimanale)")
}, width = 1400, height = 500)


save_plot("eda_09_acf_s8760.png", {
  acf(lz168, lag.max = 10000,
      main = "ACF dopo diff(s=168): struttura residua annuale")
}, width = 1400, height = 500)

save_plot("eda_10_profilo_orario.png", {
  prof_h <- tapply(z, train$hour, mean)
  plot(as.integer(names(prof_h)), prof_h, type = "b", pch = 16,
       col = COL_ACTUAL, lwd = 2,
       main = "Profilo medio orario", xlab = "Ora (1-24)", ylab = "Media")
}, width = 900, height = 500)

save_plot("eda_11_profilo_settimanale.png", {
  wd <- lubridate::wday(train$date, week_start = 1)
  prof_w <- tapply(z, wd, mean)
  barplot(prof_w, names.arg = c("Lun", "Mar", "Mer", "Gio", "Ven", "Sab", "Dom"),
          col = COL_ACTUAL, border = NA,
          main = "Profilo medio per giorno della settimana", ylab = "Media")
}, width = 900, height = 500)

save_plot("eda_12_profilo_mensile.png", {
  mth <- as.integer(format(train$date, "%m"))
  prof_m <- tapply(z, mth, mean)
  barplot(prof_m, names.arg = month.abb, col = COL_ACTUAL, border = NA,
          main = "Profilo medio mensile", ylab = "Media")
}, width = 900, height = 500)

cat(sprintf("\nGrafici salvati in %s\n", PLOTS_DIR))

