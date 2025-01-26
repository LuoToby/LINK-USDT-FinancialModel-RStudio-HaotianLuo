################################################################################
#   CHAINLINK_TIME_SERIES_ARIMA_GARCH.R
#   PROFESSIONAL SCRIPT - MINIMAL EXTRA TEXT
#
#   STEPS:
#   1) LOAD PACKAGES
#   2) DATA FETCH & AGGREGATION
#   3) BASIC EDA & OUTLIERS
#   4) STATIONARITY & SEASONALITY
#   5) SEASONAL ARIMA (AUTO & GRID)
#   6) CHECK ARCH -> GARCH
#   7) FORECAST (ARIMA + GARCH)
#   8) HOLD-OUT EVALUATION
#   9) ROLLING EVALUATION
#   10) COMPARISON
#
#   ~400+ LINES OF CODE
################################################################################

# 1) LOAD PACKAGES
################################################################################

if(!require("quantmod"))  install.packages("quantmod")
if(!require("forecast"))  install.packages("forecast")
if(!require("tseries"))   install.packages("tseries")
if(!require("rugarch"))   install.packages("rugarch")
if(!require("KFAS"))      install.packages("KFAS")
if(!require("ggplot2"))   install.packages("ggplot2")
if(!require("tsoutliers")) install.packages("tsoutliers")
install.packages("FinTS")        # install if needed

library(FinTS)
library(quantmod)
library(forecast)
library(tseries)
library(rugarch)
library(KFAS)
library(ggplot2)
library(tsoutliers)

################################################################################
# 2) DATA FETCH & AGGREGATION
################################################################################

start_date <- as.Date("2020-01-01")
end_date   <- as.Date("2025-01-01")
getSymbols("LINK-USD", src="yahoo", from=start_date, to=end_date)
link_daily <- Ad(`LINK-USD`)
link_monthly_xts <- to.monthly(link_daily, indexAt="lastof", drop.time=TRUE)
link_monthly <- Cl(link_monthly_xts)
start_yr <- as.numeric(format(start(link_monthly_xts), "%Y"))
start_mo <- as.numeric(format(start(link_monthly_xts), "%m"))
link_ts <- ts(link_monthly, start=c(start_yr, start_mo), frequency=12)
N_all <- length(link_ts)

################################################################################
# 3) EXPLORATORY & OUTLIERS
################################################################################

plot(link_ts, main="LINK Monthly", col="blue")
summary(link_ts)
link_log <- log(link_ts)
link_diff <- diff(link_log)
box_stats <- boxplot.stats(na.omit(link_diff))
outs <- box_stats$out
num_outs <- length(outs)
num_outs     ##result = 0 which means the data are correct with original one
box_u <- box_stats$stats[5]
box_l <- box_stats$stats[1]
link_diff_cap <- link_diff
link_diff_cap[link_diff_cap>box_u] <- box_u
link_diff_cap[link_diff_cap<box_l] <- box_l

################################################################################
# 4) STATIONARITY & SEASONALITY
################################################################################

d_adf <- ndiffs(link_log,  test="adf")
d_kpss <- ndiffs(link_log, test="kpss")

# Convert link_log to a time series object with frequency = 12 (12 months in a year)
link_log_ts <- ts(link_log, frequency = 12)

# Use the updated object in nsdiffs
D_seas <- nsdiffs(link_log_ts, test = "ocsb")

d_adf ##The value 1 indicates that one non-seasonal difference
d_kpss#The value 0 indicates that no non-seasonal differencing is needed
D_seas#The value 0 indicates that no seasonal differencing
adf_1 <- adf.test(na.omit(link_diff))
adf_1##result p-value is 0.1162 fail to reject H0
kpss_1 <- kpss.test(na.omit(link_diff))
kpss_1##result p-value is 0.1 fail to reject H0
# Convert `link_log` to a time series object with correct frequency
ts_data <- ts(link_log, frequency = 12, start = c(2020, 1))  # Adjust start year if needed

# Add colors and labels to make the plot easier to read
seasonplot(ts_data,
           year.labels = TRUE,        # Add year labels for better understanding
           year.labels.left = TRUE,   # Place labels on the left side as well
           col = rainbow(5),          # Use rainbow colors (adjust for your number of years)
           main = "Seasonal Plot for Your Data with Highlights",
           pch = 19                   # Use solid points to show exact values
)
###########the result shows the seasonality but not every year

##Perform decomposition to separate the trend, seasonality, and residuals for clearer insights.
decomposed <- decompose(ts(link_log, frequency = 12))
plot(decomposed)

acf(ts(link_log, frequency = 12))


################################################################################
# 5) SEASONAL ARIMA SEARCH (AUTO & GRID)
################################################################################

auto_sarima <- auto.arima(link_log, 
                          seasonal=TRUE,
                          stepwise=FALSE,
                          approximation=FALSE,
                          trace=TRUE)

summary(auto_sarima)##the result is ARIMA(1,0,0)

checkresiduals(auto_sarima)
### autoarima is not the best, so i try to use manually to catch season
best_aic    <- Inf
best_model  <- NULL
max_pq      <- 3
max_PQ      <- 2
d_use       <- max(d_adf, d_kpss)
D_use       <- D_seas

for(p in 0:max_pq){
  for(q in 0:max_pq){
    for(P in 0:max_PQ){
      for(Q in 0:max_PQ){
        fit_tmp <- tryCatch({
          Arima(link_log, order=c(p,d_use,q), 
                seasonal=c(P,D_use,Q),
                method="ML")
        }, error=function(e) NULL)
        if(!is.null(fit_tmp)){
          aic_tmp <- AIC(fit_tmp)
          if(aic_tmp < best_aic){
            best_aic   <- aic_tmp
            best_model <- fit_tmp
          }
        }
      }
    }
  }
}

best_aic
best_model
if(!is.null(best_model)){
  if(AIC(best_model) < AIC(auto_sarima)){
    final_sarima <- best_model
  } else {
    final_sarima <- auto_sarima
  }
} else {
  final_sarima <- auto_sarima
}
summary(final_sarima)
checkresiduals(final_sarima)

######getting better in p q part

# Force both regular (d=1) and seasonal (D=1) differencing,
# raise the maximum allowed orders for p, q, P, Q,
# and perform a non‐stepwise, non‐approximate search:


fit_auto_new <- auto.arima(
  link_log,
  d = 1,       # Force 1 regular difference
  D = 0,       # Force 0 seasonal differences
  max.p = 10,   # Allow AR up to order 5
  max.q = 10,   # Allow MA up to order 5
  max.P = 10,   # Allow seasonal AR up to order 2
  max.Q = 10,   # Allow seasonal MA up to order 2
  seasonal = TRUE,  # Still model seasonality, but no differencing
  method = "ML",    # Could also try "CSS"
  stepwise = FALSE, # Do a thorough search
  approximation = FALSE,
  trace = TRUE,     # Show each tried model’s AICc
  ic = "aicc"       # Use AICc to decide
)

summary(fit_auto_new)
checkresiduals(fit_auto_new)
### but we find last is still the best

################################################################################
# 6) RESIDUAL CHECK (ARCH) & BASIC UCM
################################################################################

resid_sarima <- residuals(final_sarima)
arch_tst <- ArchTest(resid_sarima, lags=12)
arch_tst

ucm_fit <- StructTS(link_log, type="BSM")
ucm_fit


