### PORTFOLIO OPTIMIZATION ###

### 3 SUMMARY STATISTICS ###

# Excess kurtosis = kurtosis(x) - 3

summary_stats <- data.frame(
  Asset          = asset_names,
  Mean           = unname(colMeans(R_all)),
  SD             = unname(apply(R_all, 2, sd)),
  Min            = unname(apply(R_all, 2, min)),
  Max            = unname(apply(R_all, 2, max)),
  Skewness       = unname(apply(R_all, 2, skewness)),
  Excess_Kurtosis = unname(apply(R_all, 2, function(x) kurtosis(x) - 3)),
  row.names      = NULL
)
print(summary_stats[, sapply(summary_stats, is.numeric)] %>% round(5))
write.csv(summary_stats, "output/summary_stats.csv", row.names = FALSE)

# Mean distributions #

mean_df <- data.frame(
  Mean = summary_stats$Mean
)

p_mean <- ggplot(
  mean_df,
  aes(x = Mean)
) +
  
  geom_histogram(
    bins = 12,
    fill = COL_PLUGIN,
    colour = "white",
    alpha = 0.85
  ) +
  
  geom_vline(
    xintercept = 0,
    colour = "grey40",
    linetype = "dashed"
  ) +
  
  labs(
    x = "Average Monthly Return",
    y = "Number of Assets"
  ) +
  
  theme_portfolio


# standard deviation #

sd_rank <- summary_stats %>%
  arrange(SD) %>%
  mutate(Asset = factor(Asset, levels = Asset))

p_sd_rank <- ggplot(
  sd_rank,
  aes(x = Asset, y = SD)
) +
  
  geom_point(
    colour = COL_BAYES,
    size = 2
  ) +
  
  coord_flip() +
  
  labs(
    title = "Monthly Volatility by Asset",
    x = NULL,
    y = "Standard Deviation"
  ) +
  
  theme_portfolio

ggsave(
  "figures/sd_returns_distribution_ggplot.png",
  p_sd_rank,
  width = 8,
  height = 5,
  dpi = 150
)


# Skewness table #
skewness_table <- data.frame(
  Asset    = asset_names,
  Skewness = unname(apply(R_all, 2, skewness)),
  row.names = NULL
)

# Shapiro-Wilk normality test #
# H0: normal distribution -->  reject H0 if p-value < 0.05 #

shapiro_table <- data.frame(
  Asset           = asset_names,
  Shapiro_p_value = unname(apply(R_all, 2, function(x) shapiro.test(x)$p.value)),
  row.names       = NULL
)
cat("Normal assets (Shapiro p > 0.05):",
    sum(shapiro_table$Shapiro_p_value > 0.05), "/", N_assets, "\n")

shapiro_plot <- shapiro_table %>%
  arrange(Shapiro_p_value) %>%
  mutate(Asset = factor(Asset, levels = Asset))

p_shapiro_rank <- ggplot(
  shapiro_plot,
  aes(x = Asset, y = Shapiro_p_value)
) +
  
  geom_point(
    colour = COL_RF,
    size = 2
  ) +
  
  geom_hline(
    yintercept = 0.05,
    colour = "red3",
    linewidth = 0.8,
    linetype = "dashed"
  ) +
  
  coord_flip() +
  
  labs(
    title = "Shapiro-Wilk Test by Asset",
    x = NULL,
    y = "p-value"
  ) +
  
  theme_portfolio


shapiro_summary <- data.frame(
  Result = c("Normality not rejected", "Normality rejected"),
  Number = c(
    sum(shapiro_table$Shapiro_p_value > 0.05),
    sum(shapiro_table$Shapiro_p_value <= 0.05)
  )
)

shapiro_summary

p_shapiro_summary <- ggplot(
  shapiro_summary,
  aes(x = Result, y = Number)
) +
  
  geom_col(
    fill = COL_RF,
    width = 0.6
  ) +
  
  geom_text(
    aes(label = Number),
    vjust = -0.5,
    size = 5
  ) +
  
  labs(
    title = "Shapiro-Wilk Normality Test Results",
    x = NULL,
    y = "Number of assets"
  ) +
  
  theme_portfolio


# Correlation and covariance (full sample) #
correlation_matrix <- cor(R_all)
covariance_matrix  <- cov(R_all)
write.csv(as.data.frame(correlation_matrix), "output/correlation_matrix.csv", row.names = TRUE)
write.csv(as.data.frame(covariance_matrix),  "output/covariance_matrix.csv",  row.names = TRUE)

png("figures/correlation_matrix.png", width = 1000, height = 800)
par(bg = BG_COLOR)
corrplot(correlation_matrix, method = "color", type = "upper", tl.cex = 0.7,
         mar   = c(0, 0, 2, 0),
         col   = colorRampPalette(c(COL_PLUGIN, "white", COL_BAYES))(200))
dev.off()

png("figures/covariance_matrix.png", width = 1000, height = 800)
par(bg = BG_COLOR)
corrplot(covariance_matrix, 
         method  = "color", 
         type    = "upper", 
         tl.cex  = 0.7,
         is.corr = FALSE,
         na.label = " ", 
         mar     = c(0, 0, 2, 0),
         col     = colorRampPalette(c("white", COL_BAYES))(200)) 
dev.off()


summary_table <- data.frame(
  Statistic = c(
    "Number of assets",
    "Mean pairwise correlation",
    "Assets rejecting normality (5%)",
    "Mean skewness",
    "Median skewness",
    "Mean excess kurtosis",
    "Median excess kurtosis",
    "Condition number"
  ),
  Value = c(
    ncol(R_all),
    round(mean(cor(R_all)[upper.tri(cor(R_all))]), 3),
    paste0(sum(shapiro_table$Shapiro_p_value < 0.05),
           " (",
           round(100*mean(shapiro_table$Shapiro_p_value < 0.05),1),
           "%)"),
    round(mean(summary_stats$Skewness),3),
    round(median(summary_stats$Skewness),3),
    round(mean(summary_stats$Excess_Kurtosis),3),
    round(median(summary_stats$Excess_Kurtosis),3),
    round(kappa(cov(R_all)),1)
  )
)

summary_table

dir.create("tables", showWarnings = FALSE)

latex_table <-
  kable(summary_table,
        format = "latex",
        booktabs = TRUE,
        caption = "Descriptive statistics of the dataset") %>%
  kable_styling(
    latex_options = c("hold_position")
  )


corr_values <- correlation_matrix[upper.tri(correlation_matrix)]

data.frame(
  Statistic = c("Mean","Median","Minimum","Maximum"),
  Value = c(
    mean(corr_values),
    median(corr_values),
    min(corr_values),
    max(corr_values)
  )
)
