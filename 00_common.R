


library(lubridate)
library(timeDate)


DATA_PATH <- "C:/Users/famil/Downloads/TimeSeries.csv"

OUTPUT_DIR <- "output"
PLOTS_DIR  <- file.path(OUTPUT_DIR, "plots")
for (d in c(OUTPUT_DIR, PLOTS_DIR)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}


HOLDOUT_DAYS  <- 30
HOLDOUT_HOURS <- HOLDOUT_DAYS * 24



save_plot <- function(filename, expr, width = 1100, height = 550) {
  path <- file.path(PLOTS_DIR, filename)
  png(path, width = width, height = height)
  on.exit(dev.off(), add = TRUE)
  force(expr)
  invisible(path)
}

COL_ACTUAL <- "steelblue4"
COL_PRED   <- "firebrick2"


plot_fit_vs_actual <- function(actual, predicted, main,
                               xlab = "Index", ylab = "Value") {
  rng <- range(c(actual, predicted), na.rm = TRUE)
  plot(actual, type = "l", col = COL_ACTUAL, lwd = 1.4,
       ylim = rng, main = main, xlab = xlab, ylab = ylab)
  lines(predicted, col = COL_PRED, lwd = 1.4)
  legend("topright", legend = c("Original", "Fitted"),
         col = c(COL_ACTUAL, COL_PRED), lwd = 1.4, bty = "n")
}


plot_forecast_tail <- function(observed, predicted, main,
                               n_tail = 24 * 14) {
  obs <- tail(observed, n_tail)
  y   <- c(obs, predicted)
  plot(y, type = "l", col = COL_ACTUAL, lwd = 1.2,
       main = main, xlab = "Ore", ylab = "Value")
  lines(x = seq_along(y)[-seq_along(obs)], y = predicted,
        col = COL_PRED, lwd = 1.2)
  abline(v = length(obs), col = "grey40", lty = 2)
  legend("topleft", legend = c("Osservato", "Previsto"),
         col = c(COL_ACTUAL, COL_PRED), lwd = 1.2, bty = "n")
}



stopifnot(file.exists(DATA_PATH))
raw <- read.csv(DATA_PATH, header = TRUE, stringsAsFactors = FALSE)
stopifnot(all(c("time", "value") %in% names(raw)))

raw$time <- as.POSIXct(raw$time, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")

# --- Timestamp mancanti ------------------------------------------------------
# 2016-03-27, 2017-03-26, 2018-03-25, 2019-03-31)
na_time <- which(is.na(raw$time))
if (length(na_time) > 0) {
  stopifnot(min(na_time) > 1)
  for (i in na_time) raw$time[i] <- raw$time[i - 1] + 3600
}
stopifnot(!any(is.na(raw$time)))

raw <- raw[order(raw$time), ]
rownames(raw) <- NULL

raw$date <- as.Date(raw$time)
raw$hour <- as.integer(format(raw$time, "%H")) + 1L

n_dup <- sum(table(raw$date, raw$hour) > 1)
if (n_dup > 0) warning("Trovate ", n_dup, " coppie (date, hour) duplicate.")


is_future <- is.na(raw$value)
stopifnot(any(is_future), any(!is_future))
stopifnot(all(diff(which(is_future)) == 1))
stopifnot(min(which(is_future)) > max(which(!is_future)))

train  <- raw[!is_future, ]
future <- raw[ is_future, ]
rownames(train)  <- NULL
rownames(future) <- NULL

H <- nrow(future) 

cat(sprintf("Osservazioni storiche: %d (%s -> %s)\n",
            nrow(train), min(train$date), max(train$date)))
cat(sprintf("Ore da prevedere:      %d (%s -> %s)\n",
            H, min(future$date), max(future$date)))

# --- Festivita' nazionali australiane ----------------------------------------


YEARS <- 2015:2020

fixed_hdays <- data.frame(
  name  = c("New Year", "Australia Day", "ANZAC Day", "Christmas", "Boxing Day"),
  month = c(1, 1, 4, 12, 12),
  day   = c(1, 26, 25, 25, 26)
)

hdays_fixed <- unlist(lapply(YEARS, function(y) {
  as.Date(paste(y, fixed_hdays$month, fixed_hdays$day, sep = "-"))
}))
hdays_fixed   <- as.Date(hdays_fixed, origin = "1970-01-01")
easter_monday <- as.Date(timeDate::EasterMonday(YEARS))
hdays_aus     <- sort(unique(c(hdays_fixed, easter_monday)))

# --- Vacanze scolastiche -----------------------------------------------------


vacations <- list(
  c("2015-04-03", "2015-04-19"), c("2015-06-27", "2015-07-12"),
  c("2015-09-19", "2015-10-05"), c("2015-12-17", "2016-01-24"),
  c("2016-04-09", "2016-04-26"), c("2016-07-02", "2016-07-18"),
  c("2016-09-24", "2016-10-09"), c("2016-12-21", "2017-01-26"),
  c("2017-04-08", "2017-04-25"), c("2017-07-01", "2017-07-17"),
  c("2017-09-23", "2017-10-08"), c("2017-12-16", "2018-01-29"),
  c("2018-04-14", "2018-04-29"), c("2018-07-07", "2018-07-22"),
  c("2018-09-29", "2018-10-14"), c("2018-12-22", "2019-01-28"),
  c("2019-04-13", "2019-04-28"), c("2019-07-06", "2019-07-21"),
  c("2019-09-28", "2019-10-13"), c("2019-12-21", "2020-01-27")
)
vacations <- lapply(vacations, function(r) as.Date(r))

is_in_vacation <- function(x, ranges) {
  as.numeric(Reduce(`|`, lapply(ranges, function(r) x >= r[1] & x <= r[2])))
}

add_calendar <- function(d) {
  d$holiday  <- as.numeric(d$date %in% hdays_aus)
  d$vacation <- is_in_vacation(d$date, vacations)
  d
}

train  <- add_calendar(train)
future <- add_calendar(future)

cat(sprintf("Giorni festivi nello storico: %d | ore in vacanza: %d (%.1f%%)\n",
            sum(tapply(train$holiday, train$date, max)),
            sum(train$vacation), 100 * mean(train$vacation)))

# --- Test set ----------------------------------------------------------------

holdout_idx <- (nrow(train) - HOLDOUT_HOURS + 1):nrow(train)
stopifnot(min(holdout_idx) > 0)

cat(sprintf("Test set: %s -> %s (%d ore)\n",
            min(train$date[holdout_idx]), max(train$date[holdout_idx]),
            length(holdout_idx)))

# --- Serie trasformata -------------------------------------------------------

train$value_log <- log(train$value + 1)

# --- Metriche ----------------------------------------------------------------

scale_mase <- function(y_train, m = 24) {
  mean(abs(diff(y_train, lag = m)), na.rm = TRUE)
}

eval_forecast <- function(actual, predicted, y_train, m = 24) {
  err <- actual - predicted
  data.frame(
    ME   = mean(err, na.rm = TRUE),
    RMSE = sqrt(mean(err^2, na.rm = TRUE)),
    MAE  = mean(abs(err), na.rm = TRUE),
    MASE = mean(abs(err), na.rm = TRUE) / scale_mase(y_train, m)
  )
}


postprocess_counts <- function(x) {
  x <- round(x)
  x[is.na(x)] <- 0
  pmax(x, 0)
}


lag_vec <- function(x, n) c(rep(NA_real_, n), head(x, -n))

