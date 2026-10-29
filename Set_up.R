### PORTFOLIO OPTIMIZATION ###

### 0. SETUP ###

# m = numero di asset (50);
# T = numero totale di osservazioni nel dataset;
# t = indice temporale (t=1,…,T);
# n = numero di osservazioni rolling window
# i = indice dell'asset.

required_packages <- c("tidyverse", "moments", "corrplot", "quadprog", "ggplot2", "scales", "zoo")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
}

library(tidyverse)
library(moments)
library(corrplot)
library(quadprog)
library(rvest)
library(dplyr)
library(ggplot2)
library(scales)
library(zoo)
library(knitr)
library(kableExtra)


dir.create("output",  showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

COL_PLUGIN <- "#d73027"
COL_BAYES  <- "#2166ac"
COL_RF     <- "#762a83"
COL_EQUAL  <- "#1b9e77"
BG_COLOR   <- "#FFFFFF"

theme_portfolio <- theme_minimal(base_size = 13) +
  theme(
    legend.position   = "top",
    plot.title        = element_text(face = "bold"),
    plot.subtitle     = element_text(colour = "grey40"),
    plot.background   = element_rect(fill = BG_COLOR, color = NA),
    panel.background  = element_rect(fill = BG_COLOR, color = NA),
    legend.background = element_rect(fill = BG_COLOR, color = NA)
  )

### 1. LOAD DATA ###

stocks      <- read.csv("~/Downloads/Assets.csv", stringsAsFactors = FALSE)
stocks$Date <- as.Date(stocks$Date)
stocks      <- stocks %>% arrange(Date)

cat("\n=== Data ===\n")
cat("Rows:   ", nrow(stocks), "\n")
cat("Assets: ", ncol(stocks) - 1, "\n")
cat("Period: ", as.character(min(stocks$Date)),
    "to", as.character(max(stocks$Date)), "\n")

# NA check #
na_check <- colSums(is.na(stocks %>% select(-Date)))
na_check <- na_check[na_check > 0]
if (length(na_check) == 0) {
  cat("NA check: no missing values found.\n\n")
} else {
  cat("NAs found — applying forward fill:\n")
  print(na_check)
  cat("\n")
  stocks <- stocks %>%
    mutate(across(where(is.numeric), ~ zoo::na.locf(., na.rm = FALSE)))
}

#stocks table #
selected_tickers <- colnames(stocks)[-1]
url <- "https://en.wikipedia.org/wiki/List_of_S%26P_500_companies"

sp500 <- url %>%
  read_html() %>%
  html_table() %>%
  .[[1]]

head(sp500)

assets_table <- sp500 %>%
  filter(Symbol %in% selected_tickers) %>%
  select(Symbol, Security, `GICS Sector`)

colnames(assets_table) <- c(
  "Ticker",
  "Company",
  "Sector"
)


# Check #

nrow(assets_table)

sector_table <- assets_table %>%
  group_by(Sector) %>%
  summarise(
    Number_of_assets = n()
  ) %>%
  arrange(desc(Number_of_assets))

sector_table <- sector_table %>%
  mutate(
    Percentage = round(
      Number_of_assets / sum(Number_of_assets) * 100,
      1
    )
  )

kable(
  assets_table,
  format = "latex",
  booktabs = TRUE,
  caption = "Selected assets included in the empirical analysis"
)

kable(
  sector_table,
  format = "latex",
  booktabs = TRUE,
  caption = "Sector composition of the selected asset universe"
)



### 2. LOG RETURNS ###

returns <- stocks %>%
  mutate(across(where(is.numeric), ~ log(. / lag(.)))) %>%
  drop_na()

R_all       <- returns %>% select(-Date) %>% as.matrix()
dates_all   <- returns$Date
N_assets    <- ncol(R_all)
asset_names <- colnames(R_all)
N_obs       <- nrow(R_all)

cat("Monthly obs: ", N_obs, "\n")
cat("Assets:      ", N_assets, "\n\n")


# Cumulative Market Return Chart: S&P 500 Performance from 2010 to 2025 #
market_return <- rowMeans(R_all)
market_wealth <- exp(cumsum(market_return))

market_df <- data.frame(
  Date = dates_all,
  Wealth = market_wealth
)

p_market <- ggplot(
  market_df,
  aes(Date, Wealth)
) +
  
  geom_line(
    colour = COL_PLUGIN,
    linewidth = 0.8
  ) +
  
  labs(
    title = "Cumulative Market Performance",
    x = NULL,
    y = "Growth of $1 Invested"
  ) +
  
  theme_portfolio


# Histogram of all monthly returns #
all_returns <- as.vector(R_all)

mu <- mean(all_returns)
sigma <- sd(all_returns)

p_hist <- ggplot(
  data.frame(Return = all_returns),
  aes(x = Return)
) +
  
  geom_histogram(
    aes(y = after_stat(density)),
    bins = 40,
    fill = COL_PLUGIN,
    colour = "white",
    linewidth = 0.3,
    alpha = 0.85
  ) +
  
  geom_density(
    colour = COL_BAYES,
    linewidth = 0.8
  ) +
  
  stat_function(
    fun = dnorm,
    args = list(mean = mu, sd = sigma),
    colour = COL_RF,
    linewidth = 0.8,
    linetype = "dashed"
  ) +
  
  geom_vline(
    xintercept = 0,
    colour = "grey50",
    linewidth = 0.5,
    linetype = "dotted"
  ) +
  
  labs(
    title = "Distribution of Monthly Log Returns",
    x = "Monthly Log Return",
    y = "Density"
  ) +
  
  theme_portfolio +
  
  theme(
    legend.position = "none"
  )

