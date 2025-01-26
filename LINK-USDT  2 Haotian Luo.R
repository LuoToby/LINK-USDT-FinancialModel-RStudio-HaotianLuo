################################################################################
# 0) Setup & Get Data
################################################################################
# install.packages(c("quantmod","rugarch","FinTS","forecast")) # If not installed

library(quantmod)
library(rugarch)
library(FinTS)
library(forecast)
library(xts)


# 0A) Fetch daily LINK data from Yahoo
start_date <- as.Date("2020-01-01")
end_date   <- as.Date("2025-01-01")
getSymbols("LINK-USD", src="yahoo", from=start_date, to=end_date)

# 0B) Create log-returns
link_daily <- Ad(`LINK-USD`)
link_log   <- log(link_daily)
link_ret   <- diff(link_log)
link_ret   <- na.omit(link_ret)

# 0C) [Optional] Check ARMA orders with auto.arima on returns
arma_search <- auto.arima(link_ret, stationary=TRUE, seasonal=FALSE, trace=TRUE)
arma_search 
# Suppose we decide ARMA(1,1) is best.

################################################################################
# 1) Define GARCH Model Search
################################################################################
garch_models <- c("sGARCH", "eGARCH", "gjrGARCH")  # GARCH families
p_candidates <- 1:2  # GARCH p in {1,2}
q_candidates <- 1:2  # GARCH q in {1,2}

# We fix ARMA(1,1) in the mean. If you want to also search ARMA, add more loops.

################################################################################
# 2) Data Structures to Store Results
################################################################################
################################################################################
# 1) Define GARCH Model Search Grid
################################################################################
garch_models <- c("sGARCH", "eGARCH", "gjrGARCH")  # GARCH families
distributions <- c("std", "norm", "sstd")          # Student-t, Normal, Skewed Student-t
p_candidates <- 1:3                                # GARCH p in {1,2,3}
q_candidates <- 1:3                                # GARCH q in {1,2,3}

# We'll use ARMA(1,1) as the mean model in all attempts
arma_order <- c(1, 1)

################################################################################
# 2) Data Structures to Store Results
################################################################################
results_df <- data.frame(
  ModelType   = character(),  # GARCH model type: sGARCH, eGARCH, gjrGARCH
  Distribution = character(), # Distribution: std, norm, sstd
  p_garch     = integer(),    # GARCH p-order
  q_garch     = integer(),    # GARCH q-order
  Converged   = logical(),    # Whether the model converged
  AIC         = numeric(),    # Akaike Information Criterion
  BIC         = numeric(),    # Bayesian Information Criterion
  LogLik      = numeric(),    # Log-Likelihood
  ARCH_LM_p   = numeric(),    # p-value of ARCH-LM test
  LjungSq_p   = numeric(),    # p-value of Ljung-Box test on squared residuals
  stringsAsFactors = FALSE
)

# List to store model fits for later inspection
fits_list <- list()

################################################################################
# 3) Loop Through Model Variants and Fit Each
################################################################################
n_obs <- length(link_ret)  # Total observations in data (for AIC/BIC scaling)

for (garch_type in garch_models) {
  for (dist in distributions) {  # Loop through distributions
    for (p_garch in p_candidates) {  # Loop through GARCH p-orders
      for (q_garch in q_candidates) {  # Loop through GARCH q-orders
        
        # Create a unique model label
        model_label <- paste0(garch_type, "_", dist, "_p", p_garch, "_q", q_garch)
        
        # Define the GARCH specification
        spec_tmp <- ugarchspec(
          variance.model = list(
            model = garch_type,        # GARCH model type (e.g., sGARCH)
            garchOrder = c(p_garch, q_garch)  # (p, q) GARCH orders
          ),
          mean.model = list(
            armaOrder = arma_order,    # Fixed ARMA(1,1) mean model
            include.mean = TRUE        # Include the mean in the model
          ),
          distribution.model = dist    # Dynamic distribution (e.g., std, norm, sstd)
        )
        
        # Try fitting the model and handle errors gracefully
        fit_tmp <- tryCatch(
          ugarchfit(spec = spec_tmp, data = link_ret),
          error = function(e) NULL
        )
        
        if (!is.null(fit_tmp)) {
          # Extract Information Criteria and Diagnostics
          ic_vals <- infocriteria(fit_tmp)  # Per-observation AIC/BIC
          aic_raw <- ic_vals["Akaike"] * n_obs  # Convert to raw AIC
          bic_raw <- ic_vals["Bayes"] * n_obs   # Convert to raw BIC
          loglik <- fit_tmp@fit$LLH            # Log-Likelihood
          
          # Residual Diagnostics
          z_resid <- residuals(fit_tmp, standardize = TRUE)
          arch_test <- ArchTest(z_resid, lags = 12)
          lb_test <- Box.test(z_resid^2, lag = 12, type = "Ljung-Box")
          
          # Append Results to Data Frame
          results_df <- rbind(
            results_df,
            data.frame(
              ModelType = garch_type,
              Distribution = dist,
              p_garch = p_garch,
              q_garch = q_garch,
              Converged = TRUE,
              AIC = aic_raw,
              BIC = bic_raw,
              LogLik = loglik,
              ARCH_LM_p = arch_test$p.value,
              LjungSq_p = lb_test$p.value,
              stringsAsFactors = FALSE
            )
          )
          
          # Save the fit for later retrieval
          fits_list[[model_label]] <- fit_tmp
          
          cat(sprintf("[OK] %s => AIC=%.3f, BIC=%.3f, LogLik=%.3f, ARCH.p=%.3f, LjungSq.p=%.3f\n",
                      model_label, aic_raw, bic_raw, loglik, arch_test$p.value, lb_test$p.value))
        } else {
          # If model fails to converge, append NA values
          results_df <- rbind(
            results_df,
            data.frame(
              ModelType = garch_type,
              Distribution = dist,
              p_garch = p_garch,
              q_garch = q_garch,
              Converged = FALSE,
              AIC = NA,
              BIC = NA,
              LogLik = NA,
              ARCH_LM_p = NA,
              LjungSq_p = NA,
              stringsAsFactors = FALSE
            )
          )
          cat(sprintf("[FAIL] %s did not converge.\n", model_label))
        }
      }
    }
  }
}

# Filter converged models
conv_df <- subset(results_df, Converged == TRUE)


######Inspect Results After running the loop, sort the results 
######by AIC or BIC to find the best model:
# Sort by AIC (or BIC)
conv_df_sorted <- conv_df[order(conv_df$AIC), ]
print(conv_df_sorted)


#####Analyze the Best Model 
best_row <- conv_df_sorted[1, ]  # Best model by AIC
best_label <- paste0(best_row$ModelType, "_", best_row$Distribution, 
                     "_p", best_row$p_garch, "_q", best_row$q_garch)
best_fit <- fits_list[[best_label]]

show(best_fit)

#####Residual Diagnostics
z_resid_best <- residuals(best_fit, standardize = TRUE)
par(mfrow = c(2, 1))
plot(z_resid_best, type = "l", col = "blue", main = "Standardized Residuals")
plot(z_resid_best^2, type = "l", col = "red", main = "Squared Residuals")
par(mfrow = c(1, 1))

adf.test(link_ret)
kpss.test(link_ret)

################################################################################
# Final Forecast with Best Model
################################################################################

# Forecast horizon (e.g., 10 days ahead)
forecast_horizon <- 10

# Ensure best_fit contains the optimal model
if (!is.null(best_fit)) {
  # Generate the forecast
  garch_forecast <- ugarchforecast(best_fit, n.ahead = forecast_horizon)
  
  # Display the forecast
  cat("\n===== GARCH Forecast =====\n")
  show(garch_forecast)
  
  # Extract forecasted mean returns
  forecasted_mean <- garch_forecast@forecast$seriesFor
  cat("\nForecasted Mean Returns:\n")
  print(forecasted_mean)
  
  # Extract forecasted volatility (conditional standard deviations)
  forecasted_volatility <- garch_forecast@forecast$sigmaFor
  cat("\nForecasted Volatility (Standard Deviations):\n")
  print(forecasted_volatility)
  
  # Visualize the forecast
  par(mfrow = c(2, 1))
  plot(forecasted_mean, type = "l", col = "blue", main = "Forecasted Mean Returns", xlab = "Days Ahead", ylab = "Returns")
  plot(forecasted_volatility, type = "l", col = "red", main = "Forecasted Volatility", xlab = "Days Ahead", ylab = "Volatility")
  par(mfrow = c(1, 1))
} else {
  cat("Best fit model is NULL. Ensure the model has been identified before forecasting.")
}



###########comparing with real data


library(xts)


# 确保必要的数据
price_data <- link_daily  # 历史价格数据
forecasted_volatility <- garch_forecast@forecast$sigmaFor  # GARCH 模型预测的波动率

# 确保基准时间从固定预测日期开始，例如 2025-01-01
base_date <- as.Date("2025-01-01")  # 设置预测起点
forecast_dates <- seq.Date(from = base_date, by = "days", length.out = length(forecasted_volatility))

# 计算历史波动率（用对数收益计算）
log_returns <- diff(log(price_data))  # 对数收益
log_returns <- na.omit(log_returns)  # 移除 NA 值

# 滚动计算历史波动率
rolling_window <- 20  # 定义滚动窗口大小，例如 20 天
historical_volatility <- rollapply(
  data = log_returns,
  width = rolling_window,
  FUN = sd,
  align = "right",
  fill = NA
)

# 确保历史波动率与日期对齐
historical_volatility <- na.omit(historical_volatility)  # 去掉滚动计算中的 NA 值
historical_volatility_xts <- xts(historical_volatility, order.by = index(log_returns)[rolling_window:length(log_returns)])

# 创建预测波动率的 xts 对象
forecasted_volatility_xts <- xts(as.numeric(forecasted_volatility), order.by = forecast_dates)

# 可视化比较历史波动率和预测波动率
plot(index(historical_volatility_xts), historical_volatility_xts, type = "l", col = "blue", lwd = 2,
     xlab = "Date", ylab = "Volatility", main = "Historical vs Forecasted Volatility")
lines(index(forecasted_volatility_xts), forecasted_volatility_xts, col = "red", lwd = 2)
legend("topright", legend = c("Historical Volatility", "Forecasted Volatility"), col = c("blue", "red"), lwd = 2)


######monte carlo simulation
# 假设您已经拟合了一个 GARCH 模型
spec <- ugarchspec(mean.model = list(armaOrder = c(1, 1)),
                   variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
                   distribution.model = "std")

fit <- ugarchfit(spec = spec, data = link_ret)

# 1. 使用拟合结果进行蒙特卡洛模拟
set.seed(123)
n.sim <- 1000  # 模拟路径数
simulations <- ugarchsim(fit, n.sim = length(link_ret), m.sim = n.sim)

# 2. 提取模拟数据
simulated_returns <- fitted(simulations)  # 模拟的收益率矩阵
simulated_volatility <- sigma(simulations)  # 模拟的波动率矩阵

# 3. 计算模拟波动率的统计特性
# 计算模拟数据每条路径的平均波动率
simulated_mean_vol <- apply(simulated_volatility, 2, mean)
simulated_sd_vol <- apply(simulated_volatility, 2, sd)

# 实际数据的波动率
realized_volatility <- sigma(fit)  # 模型拟合的历史波动率
realized_mean_vol <- mean(realized_volatility)
realized_sd_vol <- sd(realized_volatility)

# 4. 可视化比较
# 绘制模拟的波动率分布与实际波动率
hist(simulated_mean_vol, breaks = 50, col = "blue", main = "Simulated vs Realized Volatility",
     xlab = "Volatility", xlim = range(c(simulated_mean_vol, realized_mean_vol)))
abline(v = realized_mean_vol, col = "red", lwd = 2, lty = 2)  # 实际波动率均值
legend("topright", legend = c("Simulated", "Realized Mean"), col = c("blue", "red"), lty = c(1, 2), lwd = 2)



