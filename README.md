
<!-- README.md is generated from README.Rmd. Please edit that file -->

# cforecast

**cforecast** is an R package for conducting conditional forecasts and
scenario analysis using vector autoregressive (VAR) models. It
implements the Kalman filtering methodology proposed by Clarida and
Coyle (1984) and Banbura, Giannone, and Lenza (2015), allowing users to
simulate forecast paths under imposed constraints on future values of
selected variables.

## Installation

You can install the development version from GitHub using:

``` r
# install.packages("devtools")
devtools::install_github("timginker/cforecast")
```

## Example

``` r
library(cforecast)
library(vars)
library(ggplot2)

# Load data and fit VAR model
data(Canada)
fit <- VAR(Canada, p = 2, type = "const")

# Define a scenario: increase unemployment (4th variable) to 15 for 3 periods
cond_path <- matrix(rep(15, 3), ncol = 1)
cond_var <- 4

# Generate conditional forecast
conditional_forecast <- cforecast(fit = fit, cond_path = cond_path, cond_var = cond_var)

# Extract forecast and standard deviation
forecast <- conditional_forecast$forecast[, 1]
std_dev <- sqrt(conditional_forecast$mse[1, 1, ])
horizon <- length(forecast)
time <- seq_len(horizon)

df <- data.frame(
  Time = time,
  Forecast = forecast,
  Lower = forecast - 1.645 * std_dev,
  Upper = forecast + 1.645 * std_dev,
  Series = "Forecast"
)

# Plot
ggplot(df, aes(x = Time, y = Forecast)) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = "grey70", alpha = 0.4) +
  geom_line(color = "black", linewidth = 0.6) +
  labs(
    x = "Horizon",
    y = "Forecasted Value",
    title = "Conditional Forecast with 90% Confidence Interval"
  ) +
  theme_minimal(base_size = 12, base_family = "serif") +
  theme(
    legend.position = "none",
    axis.title = element_text(face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()
  )
```

<img src="man/figures/README-unnamed-chunk-2-1.png" width="100%" />

## References

- Clarida, R., & Coyle, D. (1984). *Conditional Projection by Means of
  Kalman Filtering*.
- Banbura, M., Giannone, D., & Lenza, M. (2015). *Conditional forecasts
  and scenario analysis with vector autoregressions for large
  cross-sections*. *International Journal of Forecasting*, 31(3),
  739–756.

# Disclaimer

The views expressed here are solely of the author and do not necessarily
represent the views of the Bank of Israel.

Please note that `cforecast` is still under development and may contain
bugs or other issues that have not yet been resolved. While we have made
every effort to ensure that the package is functional and reliable, we
cannot guarantee its performance in all situations.

We strongly advise that you regularly check for updates and install any
new versions that become available, as these may contain important bug
fixes and other improvements. By using this package, you acknowledge and
accept that it is provided on an “as is” basis, and that we make no
warranties or representations regarding its suitability for your
specific needs or purposes.
