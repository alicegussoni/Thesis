### PORTFOLIO OPTIMIZATION ###

### 4. MODELS IMPLEMENTATION ###

# Equal Weight Portfolio #

equal_weight <- function(R_window){
  
  p <- ncol(R_window)
  
  rep(1/p, p)
  
}

# Plug-in: classical Markowitz, long-only, fully invested #

plugin_weights <- function(R_window, lambda = 5) {
  #plug-in estimator
  p         <- ncol(R_window)
  mu_hat    <- colMeans(R_window)
  #Small ridge adjustment to improve numerical stability
  Sigma_hat <- cov(R_window) + diag(1e-6, p)
  
  #mean-variance optimization
  Dmat <- 2 * lambda * Sigma_hat
  dvec <- mu_hat
  
  #constraints
  Amat <- cbind(rep(1, p), diag(p))
  bvec <- c(1, rep(0, p))
  meq  <- 1
  
  sol <- tryCatch(
    solve.QP(Dmat = Dmat, dvec = dvec, Amat = Amat, bvec = bvec, meq = meq),
    error = function(e) NULL
  )
  if (is.null(sol)) return(rep(1 / p, p))
  w <- pmax(sol$solution, 0)
  w / sum(w)
}

# Bayesian NIW: Normal-Inverse-Wishart conjugate prior #

bayes_weights <- function(R_window, lambda = 5,
                          kappa0    = 25,
                          nu0_extra = 20) {
  X  <- as.matrix(R_window)
  n  <- nrow(X)
  p  <- ncol(X)
  
  xbar     <- colMeans(X)
  S_sample <- cov(X)
  
  # Weakly informative Normal-Inverse-Wishart prior
  mu0 <- rep(0, p)
  nu0 <- p + nu0_extra
  
  # Empirical Bayes prior scale matrix
  prior_var <- mean(diag(S_sample))
  S0 <- prior_var * diag(p) * (nu0 - p - 1)
  
  #  Posterior update 
  kappa_n <- kappa0 + n
  nu_n    <- nu0 + n
  mu_n    <- (kappa0 * mu0 + n * xbar) / kappa_n
  diff_mu <- matrix(xbar - mu0, ncol = 1)
  S_n     <- S0 +
    (n - 1) * S_sample +
    (kappa0 * n / kappa_n) * (diff_mu %*% t(diff_mu))
  
  # Posterior predictive moments
  scale_factor <- (kappa_n + 1) / (kappa_n * (nu_n - p - 1))
  # Small ridge adjustment for numerical stability
  Sigma_pred   <- scale_factor * S_n + diag(1e-6, p)
  
  # Mean-variance optimization using posterior predictive moments
  Dmat <- 2 * lambda * Sigma_pred
  dvec <- mu_n
  Amat <- cbind(rep(1, p), diag(p))
  bvec <- c(1, rep(0, p))
  meq  <- 1
  
  sol <- tryCatch(
    solve.QP(Dmat = Dmat, dvec = dvec, Amat = Amat, bvec = bvec, meq = meq),
    error = function(e) NULL
  )
  
  if (is.null(sol)) return(rep(1 / p, p))
  w <- pmax(sol$solution, 0)
  w / sum(w)
}


# Machine Learning algorithm: Random Forest #

#prova 7/07: funziona
library(ranger)
library(future.apply)

plan(multisession,
     workers = parallel::detectCores()-1)

# Random Forest dataset #
create_rf_dataset <- function(asset_return,
                              market_return,
                              market_vol,
                              n_lags = 6){
  
  n <- length(asset_return)
  
  df <- data.frame(
    Asset_Lag1 = asset_return[n_lags:(n-1)]
  )
  
  # Asset lags
  if(n_lags > 1){
    
    for(l in 2:n_lags){
      
      df[[paste0("Asset_Lag",l)]] <-
        asset_return[(n_lags-l+1):(n-l)]
    }
  }
  
  # Market lags
  for(l in 1:n_lags){
    df[[paste0("Market_Lag",l)]] <-
      market_return[(n_lags-l+1):(n-l)]
  }
  
  # Market volatility
  df$Market_Vol <-
    market_vol[(n_lags+1):n]
  
  # Target
  df$Target <-
    asset_return[(n_lags+1):n]
  
  return(df)
}

#check
test <- create_rf_dataset(
  R_all[1:60,1],
  rowMeans(R_all[1:60,]),
  rep(sd(rowMeans(R_all[1:60,])),60),
  n_lags = 6
)

dim(test)
head(test)


# Random Forest weights + feature importance #
rf_weights <- function(R_window,
                       lambda = 5,
                       n_lags = 6,
                       num_trees = 300,
                       shrinkage = 0.2){
  
  
  R_window <- as.matrix(R_window)
  
  p <- ncol(R_window)
  
  # Market variables
  market_return <- rowMeans(R_window)

  market_vol <- zoo::rollapply(
    market_return,
    width = 12,
    FUN = sd,
    fill = NA,
    align = "right"
  )
  
  # Replace initial NA values
  market_vol[is.na(market_vol)] <- 
    mean(market_vol, na.rm = TRUE)
  
# Expected returns + feature importance #
  mu_rf <- numeric(p)
  
  importance_list <- list()
  
  for(asset in 1:p){
    
    asset_return <- R_window[,asset]
    
    rf_data <- create_rf_dataset(
      asset_return,
      market_return,
      market_vol,
      n_lags
    )
    
    
    # Train 
    rf_model <- ranger(
      Target ~ .,
      data = rf_data,
      num.trees = num_trees,
      mtry = max(
        1,
        floor(sqrt(ncol(rf_data)-1))
      ),
      min.node.size = 10,
      importance = "impurity",
      seed = 123
    )
    importance_list[[asset]] <- ranger::importance(rf_model)
    
    # Forecast #
    new_x <- data.frame(
      matrix(
        nrow = 1,
        ncol = 2*n_lags + 1
      )
    )
    
    
    colnames(new_x) <- c(
      paste0("Asset_Lag",1:n_lags),
      paste0("Market_Lag",1:n_lags),
      "Market_Vol"
    )
    
    for(l in 1:n_lags){
      new_x[[paste0("Asset_Lag",l)]] <-
        asset_return[length(asset_return)-l+1]
    }
    
    for(l in 1:n_lags){
      new_x[[paste0("Market_Lag",l)]] <-
        market_return[length(market_return)-l+1]
    }
    
    new_x$Market_Vol <- tail(market_vol,1)
    
    mu_rf[asset] <-
      predict(
        rf_model,
        data = new_x
      )$predictions
  }
  
  # Shrinkage expected returns #
 mu_hist <- colMeans(R_window)
  
  mu_hat <-
    shrinkage * mu_rf +
    (1-shrinkage)*mu_hist
  

  # Covariance matrix #
  Sigma_hat <- cov(R_window)
  
  Sigma_hat <- Sigma_hat +
    diag(1e-6,p)
  
  # Markowitz optimization #
  Dmat <- 2 * lambda * Sigma_hat
  
  dvec <- mu_hat
  
  Amat <- cbind(
    rep(1,p),
    diag(p)
  )
  
  bvec <- c(
    1,
    rep(0,p)
  )
  
  sol <- tryCatch(
    solve.QP(
      Dmat = Dmat,
      dvec = dvec,
      Amat = Amat,
      bvec = bvec,
      meq = 1
    ),
    error = function(e) NULL
  )
  
  if(is.null(sol)){
    weights <- rep(1/p,p)
  } else {
    w <- pmax(sol$solution,0)
    if(sum(w)==0){
      weights <- rep(1/p,p)
    } 
    else {
      weights <- w/sum(w)
    }
  }
  
  # Attach feature importance #
  attr(weights,"importance") <- importance_list

  return(weights)
}


### 5. SENSITIVITY ANALYSIS ###
run_backtest <- function(window_size,
                         R_all,
                         dates_all,
                         lambda = 5){
  
  T <- nrow(R_all)
  p <- ncol(R_all)

  
  ret_plugin <- rep(NA, T)
  ret_bayes  <- rep(NA, T)
  ret_rf     <- rep(NA, T)
  ret_equal <- rep(NA, T)
  
  W_plugin <- matrix(NA, T, p)
  W_bayes  <- matrix(NA, T, p)
  W_rf     <- matrix(NA, T, p)
  
  
  RF_importance <- list()
  
  pb <- txtProgressBar(
    min = window_size + 1,
    max = T,
    style = 3
  )
  
  # Rolling window #
  for(t in (window_size+1):T){
    R_window <- R_all[(t-window_size):(t-1), ]
        r_next <- R_all[t, ]
        
    ## Portfolio weights
    w_plugin <- plugin_weights(
      R_window,
      lambda = lambda
    )
    
    w_bayes <- bayes_weights(
      R_window,
      lambda = lambda
    )
    
    w_equal <- equal_weight(R_window)
    
    w_rf <- rf_weights(
      R_window,
      lambda = lambda,
      n_lags = 6,
      num_trees = 300,
      shrinkage = 0.30
    )
    RF_importance[[t]] <- attr(w_rf,"importance")
    
    # Save weights
    W_plugin[t, ] <- w_plugin
    W_bayes[t, ]  <- w_bayes
    W_rf[t, ]     <- w_rf
    
    # Portfolio returns
    ret_plugin[t] <- log(1 + sum(w_plugin * (exp(r_next)-1)))
    ret_bayes[t]  <- log(1 + sum(w_bayes  * (exp(r_next)-1)))
    ret_rf[t]     <- log(1 + sum(w_rf     * (exp(r_next)-1)))
    ret_equal[t] <- log(1 + sum(w_equal * (exp(r_next)-1)))
    
    setTxtProgressBar(pb,t)
  }

  close(pb)
  
  keep <- (window_size+1):T
  
  ret_plugin <- ret_plugin[keep]
  ret_bayes  <- ret_bayes[keep]
  ret_rf     <- ret_rf[keep]
  ret_equal  <- ret_equal[keep]
  
  W_plugin <- W_plugin[keep,]
  W_bayes  <- W_bayes[keep,]
  W_rf     <- W_rf[keep,]
  
  backtest_dates <- dates_all[keep]
  
  # Performance metrics #
  m_plugin <- performance_metrics(ret_plugin)
  m_bayes  <- performance_metrics(ret_bayes)
  m_rf     <- performance_metrics(ret_rf)
  m_equal <- performance_metrics(ret_equal)

  # a.Turnover
  turnover_plugin <- compute_turnover(W_plugin)
  turnover_bayes  <- compute_turnover(W_bayes)
  turnover_rf     <- compute_turnover(W_rf)

  ## Results table
  results <- data.frame(
    Window = window_size,
    Method = c(
      "Equal Weight",
      "Plugin",
      "Bayesian",
      "Random Forest"
    ),
    
    Annual_Return = c(
      m_equal$Annualized_Return,
      m_plugin$Annualized_Return,
      m_bayes$Annualized_Return,
      m_rf$Annualized_Return
    ),
    
    Annual_Volatility = c(
      m_equal$Annualized_Volatility,
      m_plugin$Annualized_Volatility,
      m_bayes$Annualized_Volatility,
      m_rf$Annualized_Volatility
    ),
    
    Sharpe = c(
      m_equal$Sharpe_Ratio,
      m_plugin$Sharpe_Ratio,
      m_bayes$Sharpe_Ratio,
      m_rf$Sharpe_Ratio
    ),
    
    Max_Drawdown = c(
      m_equal$Max_Drawdown,
      m_plugin$Max_Drawdown,
      m_bayes$Max_Drawdown,
      m_rf$Max_Drawdown
    ),
    
    Turnover = c(
      0,
      turnover_plugin,
      turnover_bayes,
      turnover_rf
    )
    
  )
  
  W_equal <- matrix(
    rep(1/p, length(keep)*p),
    ncol = p,
    byrow = TRUE
  )
  return(
    list(
      summary = results,
      returns = data.frame(
        Date = backtest_dates,
        Equal_Weight = ret_equal,
        Plugin = ret_plugin,
        Bayesian = ret_bayes,
        RandomForest = ret_rf
      ),
      
      weights = list(
        Equal_Weight = W_equal,
        Plugin = W_plugin,
        Bayesian = W_bayes,
        RandomForest = W_rf
      ),
      RF_importance = RF_importance,
      dates = backtest_dates
    )
  )
}


performance_metrics <- function(r){
  ann_return <- mean(r) * 12
  ann_vol <- sd(r) * sqrt(12)
  sharpe <- ifelse(
    ann_vol > 0,
    ann_return / ann_vol,
    NA
  )
  
  wealth <- c(1, exp(cumsum(r)))
  
  running_max <- cummax(wealth)
  
  drawdown <- (wealth - running_max) / running_max
  
  max_dd <- min(drawdown)
  
  data.frame(
    Annualized_Return = ann_return,
    Annualized_Volatility = ann_vol,
    Sharpe_Ratio = sharpe,
    Max_Drawdown = max_dd
  )
}


compute_turnover <- function(W){
  if(nrow(W) < 2){
    return(NA)
  }
  mean(
    rowSums(
      abs(diff(W))
    )
  )
}


# Run sensitivity analysis #
windows <- c(60, 96, 120)

backtests <- lapply(
  windows,
  function(w)
    run_backtest(
      window_size = w,
      R_all = R_all,
      dates_all = dates_all,
      lambda = 5
    )
)

sensitivity_results <- do.call(
  rbind,
  lapply(backtests, function(x) x$summary)
)

print(sensitivity_results)

write.csv(
  sensitivity_results,
  "output/sensitivity_analysis.csv",
  row.names = FALSE
)

final_backtest <- backtests[[3]]

# Sensitivity analysis for shrinkage parameter #
run_backtest_shrikage <- function(window_size,
                                  R_all,
                                  dates_all,
                                  lambda = 5,
                                  shrinkage = 0.2){
  
  T <- nrow(R_all)
  p <- ncol(R_all)
  
  ret_plugin <- rep(NA, T)
  ret_bayes  <- rep(NA, T)
  ret_rf     <- rep(NA, T)
  ret_equal <- rep(NA, T)
  
  W_plugin <- matrix(NA, T, p)
  W_bayes  <- matrix(NA, T, p)
  W_rf     <- matrix(NA, T, p)
  
  
  RF_importance <- list()
  
  pb <- txtProgressBar(
    min = window_size + 1,
    max = T,
    style = 3
  )
  
  # Rolling window
  for(t in (window_size+1):T){
    R_window <- R_all[(t-window_size):(t-1), ]
    r_next <- R_all[t, ]
    
    ## Portfolio weights
    w_plugin <- plugin_weights(
      R_window,
      lambda = lambda
    )
    
    w_bayes <- bayes_weights(
      R_window,
      lambda = lambda
    )
    
    w_equal <- equal_weight(R_window)
    
    w_rf <- rf_weights(
      R_window,
      lambda = lambda,
      n_lags = 6,
      num_trees = 300,
      shrinkage = shrinkage
    )
    RF_importance[[t]] <- attr(w_rf,"importance")
    
    # weights
    W_plugin[t, ] <- w_plugin
    W_bayes[t, ]  <- w_bayes
    W_rf[t, ]     <- w_rf
    
    ## Portfolio returns
    ret_plugin[t] <- log(1 + sum(w_plugin * (exp(r_next)-1)))
    ret_bayes[t]  <- log(1 + sum(w_bayes  * (exp(r_next)-1)))
    ret_rf[t]     <- log(1 + sum(w_rf     * (exp(r_next)-1)))
    ret_equal[t] <- log(1 + sum(w_equal * (exp(r_next)-1)))
    
    setTxtProgressBar(pb,t)
  }
  
  close(pb)
  
  keep <- (window_size+1):T
  
  ret_plugin <- ret_plugin[keep]
  ret_bayes  <- ret_bayes[keep]
  ret_rf     <- ret_rf[keep]
  ret_equal  <- ret_equal[keep]
  
  W_plugin <- W_plugin[keep,]
  W_bayes  <- W_bayes[keep,]
  W_rf     <- W_rf[keep,]
  
  backtest_dates <- dates_all[keep]
  
  # Performance metrics #
  m_plugin <- performance_metrics(ret_plugin)
  m_bayes  <- performance_metrics(ret_bayes)
  m_rf     <- performance_metrics(ret_rf)
  m_equal <- performance_metrics(ret_equal)
  
  # Turnover #
  turnover_plugin <- compute_turnover(W_plugin)
  turnover_bayes  <- compute_turnover(W_bayes)
  turnover_rf     <- compute_turnover(W_rf)

  results <- data.frame(
    Window = window_size,
    Method = c(
      "Equal Weight",
      "Plugin",
      "Bayesian",
      "Random Forest"
    ),
    Annual_Return = c(
      m_equal$Annualized_Return,
      m_plugin$Annualized_Return,
      m_bayes$Annualized_Return,
      m_rf$Annualized_Return
    ),
    Annual_Volatility = c(
      m_equal$Annualized_Volatility,
      m_plugin$Annualized_Volatility,
      m_bayes$Annualized_Volatility,
      m_rf$Annualized_Volatility
    ),
    Sharpe = c(
      m_equal$Sharpe_Ratio,
      m_plugin$Sharpe_Ratio,
      m_bayes$Sharpe_Ratio,
      m_rf$Sharpe_Ratio
    ),
    Max_Drawdown = c(
      m_equal$Max_Drawdown,
      m_plugin$Max_Drawdown,
      m_bayes$Max_Drawdown,
      m_rf$Max_Drawdown
    ),
    Turnover = c(
      0,
      turnover_plugin,
      turnover_bayes,
      turnover_rf
    )
  )
 
  W_equal <- matrix(
    rep(1/p, length(keep)*p),
    ncol = p,
    byrow = TRUE
  )
  
  return(
    list(
      summary = results,
      returns = data.frame(
        Date = backtest_dates,
        Equal_Weight = ret_equal,
        Plugin = ret_plugin,
        Bayesian = ret_bayes,
        RandomForest = ret_rf
      ),
      
      weights = list(
        Equal_Weight = W_equal,
        Plugin = W_plugin,
        Bayesian = W_bayes,
        RandomForest = W_rf
      ),
      RF_importance = RF_importance,
      dates = backtest_dates
    )
  )
}

# Performance metrics #
performance_metrics <- function(r){
  ann_return <- mean(r) * 12
  ann_vol <- sd(r) * sqrt(12)
  sharpe <- ifelse(
    ann_vol > 0,
    ann_return / ann_vol,
    NA
  )
  
  wealth <- c(1, exp(cumsum(r)))
  
  running_max <- cummax(wealth)
  
  drawdown <- (wealth - running_max) / running_max
  
  max_dd <- min(drawdown)
  
  data.frame(
    Annualized_Return = ann_return,
    Annualized_Volatility = ann_vol,
    Sharpe_Ratio = sharpe,
    Max_Drawdown = max_dd
  )
}

# Portfolio turnover #
compute_turnover <- function(W){
  if(nrow(W) < 2){
    return(NA)
  }
  mean(
    rowSums(
      abs(diff(W))
    )
  )
}

# Run sensitivity analysis #
windows <- c(60, 96, 120)

backtests <- lapply(
  windows,
  function(w)
    run_backtest(
      window_size = w,
      R_all = R_all,
      dates_all = dates_all,
      lambda = 5
    )
)

sensitivity_results <- do.call(
  rbind,
  lapply(backtests, function(x) x$summary)
)

print(sensitivity_results)

write.csv(
  sensitivity_results,
  "output/sensitivity_analysis.csv",
  row.names = FALSE
)

final_backtest <- backtests[[3]]

# Shrinkage sensitivity analysis #
alpha_grid <- c(0, 0.2, 0.4, 0.6, 0.8, 1)

alpha_backtests <- lapply(
  alpha_grid,
  function(a){
    
    cat("\n============================\n")
    cat("Shrinkage alpha =", a, "\n")
    cat("============================\n")
    
    run_backtest_shrikage(
      window_size = 120,
      R_all = R_all,
      dates_all = dates_all,
      lambda = 5,
      shrinkage = a
    )
  }
)

alpha_results <- do.call(
  rbind,
  lapply(seq_along(alpha_grid), function(i){
    tmp <- alpha_backtests[[i]]$summary
    tmp <- tmp[tmp$Method=="Random Forest",]
    tmp$Alpha <- alpha_grid[i]
    tmp
  })
)

alpha_results <- alpha_results[,
                               c("Alpha",
                                 "Annual_Return",
                                 "Annual_Volatility",
                                 "Sharpe",
                                 "Max_Drawdown",
                                 "Turnover")]

print(alpha_results)

write.csv(
  alpha_results,
  "output/shrinkage_sensitivity.csv",
  row.names = FALSE
)


### 6. BACKTESTING (120-MONTH WINDOW) ###
returns_final <- final_backtest$returns

weights_plugin <- final_backtest$weights$Plugin
weights_bayes  <- final_backtest$weights$Bayesian
weights_rf     <- final_backtest$weights$RandomForest

backtest_dates <- final_backtest$dates

plugin_ret <- returns_final$Plugin
bayes_ret  <- returns_final$Bayesian
rf_ret     <- returns_final$RandomForest
equal_ret <- returns_final$Equal_Weight


# Final performance metrics #
m_plugin <- performance_metrics(plugin_ret)
m_bayes  <- performance_metrics(bayes_ret)
m_rf     <- performance_metrics(rf_ret)
m_equal <- performance_metrics(equal_ret)

results_final <- data.frame(
  Method = c(
    "Equal Weight",
    "Plugin",
    "Bayesian",
    "Random Forest"
  ),
  
  Annual_Return = c(
    m_equal$Annualized_Return,
    m_plugin$Annualized_Return,
    m_bayes$Annualized_Return,
    m_rf$Annualized_Return
  ),
  
  Annual_Volatility = c(
    m_equal$Annualized_Volatility,
    m_plugin$Annualized_Volatility,
    m_bayes$Annualized_Volatility,
    m_rf$Annualized_Volatility
  ),
  
  Sharpe = c(
    m_equal$Sharpe_Ratio,
    m_plugin$Sharpe_Ratio,
    m_bayes$Sharpe_Ratio,
    m_rf$Sharpe_Ratio
  ),
  
  Max_Drawdown = c(
    m_equal$Max_Drawdown,
    m_plugin$Max_Drawdown,
    m_bayes$Max_Drawdown,
    m_rf$Max_Drawdown
  ),
  
  Turnover = c(
    0,
    compute_turnover(weights_plugin),
    compute_turnover(weights_bayes),
    compute_turnover(weights_rf)
  )
)

print(results_final)

write.csv(
  results_final,
  "output/final_results.csv",
  row.names = FALSE
)



