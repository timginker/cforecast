#' U.S. Macroeconomic Dataset from FRED
#'
#' Source: Federal Reserve Economic Data (FRED)
#' \url{https://fred.stlouisfed.org/}
#'
#' Data vintage: 2026-01-27.
#'
#' All variables are obtained from the Federal Reserve Economic Data (FRED)
#' database and converted to quarterly frequency. Observations are aligned
#' to quarter-end dates to ensure consistency across series and compatibility
#' with quarterly VAR estimation and forecasting.
#'
#' Interest-rate and spread variables are averaged within each quarter,
#' while price indices use end-of-period quarterly values.
#'
#' @details
#' \strong{Variable definitions:}
#'
#' \describe{
#'   \item{GDPC1}{Real Gross Domestic Product. Inflation-adjusted U.S. GDP,
#'   expressed as the quarter-to-quarter percent change.}
#'
#'   \item{PCEPILFE}{Core PCE Price Index. Personal Consumption Expenditures
#'   price index excluding food and energy (“core PCE”), expressed as the
#'   quarter-to-quarter percent change using end-of-period quarterly values.}
#'
#'   \item{FEDFUNDS}{Federal Funds Effective Rate. The effective federal
#'   funds rate, the primary short-term monetary policy rate in the United
#'   States. Constructed as the quarterly average and expressed in
#'   annualized percentage points (levels).}
#'
#'   \item{BAA10YM}{Corporate Credit Spread. The spread between Moody’s
#'   seasoned Baa corporate bond yield and the 10-year U.S. Treasury
#'   constant maturity yield. Averaged within each quarter and expressed
#'   in percentage points.}
#'
#'   \item{DCOILWTICO}{Crude Oil Price (WTI). Spot price of West Texas
#'   Intermediate (WTI) crude oil at Cushing, Oklahoma. Quarterly average growth rate (quarter-to-quarter
#'   percent change).}
#' }
#'
#' \strong{Date convention:}
#' All quarterly observations are aligned to the last calendar day of the
#' quarter to ensure temporal consistency across variables and compatibility
#' with quarterly VAR and forecasting analysis.
#'
#' @keywords datasets
"fred_macro"
