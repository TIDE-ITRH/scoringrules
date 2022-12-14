
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Evaluating Probabilistic Forecasts for Maritime Engineering Operations

<!-- badges: start -->

<!-- badges: end -->

This is the GitHub repo for the translational paper ‘’Evaluating
Probabilistic Forecasts for Maritime Engineering Operations’’ submitted
to *Data-centric Engineering*. The main purpose of the paper is to make
existing statistical methodology available to an engineering audience
with the hope that it may be more widely adopted. The methodology is
then supported by a case study of forecasting surface winds at a
location in Australia’s North–West Shelf. This repository is broken into
two sections: the first, and main section, shows code and methods for
calculating three chosen scoring rules given a probabilistic forecast;
the second section documents the code for the case study in the paper.
Some data is provided alongside. Note, we heavily use
[tidyverse](https://www.tidyverse.org) packages and syntax; in
particular, [pipes](https://style.tidyverse.org/pipes.html).

## Scoring Rules

To fully assess a probabilistic forecast we advocate three scoring
rules: the squared error (SE); the Dawid-Sebastiani score (DSS); and the
continuous ranked probability score (CRPS) or its multivariate analogue
the energy score (ES). Roughly, the SE assesses the forecast mean, the
DSS assesses the forecast tails, and the CRPS/ES assesses the main body
of the forecast. Definitions and discussions are in the paper. First we
introduce some data and the forecasts of a probabilistic forecasting
model, then we show how the scores may be calculated.

### Data and Forecasting Output

We use the predictions from the Model Output Statistics (MOS) model made
at 00:00:00AWST 06-07-2018 (this is the same time as shown in Figure 3
of the paper). The MOS model is a Bayesian forecasting model estimated
via Markov chain Monte Carlo (MCMC). As such, parametric forecast
distributions are not available at each prediction horizon; rather, we
have a collection of samples. The tibble (the particular data structure
we use to store the forecasts) of the forecasts is shown below; it
documents the time issued (*t*), the time forecast (*t+h*), the
prediction horizon (*h*), the sample from the forecasting model (*j*),
the true measured eastings and northings at the time forecast, and the
forecast eastings and northings.

``` r
library(tidyverse)

mos_predictions <- readRDS("data/mos_predictions.RDS")
mos_predictions
#> # A tibble: 56,000 × 8
#>    time_issued         time_forecast       horizon sample meas_east meas_north
#>    <dttm>              <dttm>                <dbl>  <int>     <dbl>      <dbl>
#>  1 2018-07-06 00:00:00 2018-07-06 00:00:00       0      1 -6.58e- 2      -3.77
#>  2 2018-07-06 00:00:00 2018-07-06 01:00:00       1      1  2.25e- 1      -3.22
#>  3 2018-07-06 00:00:00 2018-07-06 02:00:00       2      1  5.76e- 1      -3.63
#>  4 2018-07-06 00:00:00 2018-07-06 03:00:00       3      1  1.36e+ 0      -1.62
#>  5 2018-07-06 00:00:00 2018-07-06 04:00:00       4      1  4.24e- 1      -1.58
#>  6 2018-07-06 00:00:00 2018-07-06 05:00:00       5      1  8.14e-17       1.33
#>  7 2018-07-06 00:00:00 2018-07-06 06:00:00       6      1  1.09e- 1      -3.11
#>  8 2018-07-06 00:00:00 2018-07-06 07:00:00       7      1 -3.78e- 1      -3.08
#>  9 2018-07-06 00:00:00 2018-07-06 08:00:00       8      1 -1.24e+ 0      -3.80
#> 10 2018-07-06 00:00:00 2018-07-06 09:00:00       9      1 -3.58e- 1      -2.54
#> # … with 55,990 more rows, and 2 more variables: fore_east <dbl>,
#> #   fore_north <dbl>
```

### Calculating Scoring Rules

**Squared error** calculates the distance between the mean of the
forecast distribution and the observed values. The forecast mean is
approximated by the sample mean and the SE is calculated at each
prediction horzion as per Equation 3.1.

``` r
mos_mean <- mos_predictions %>%
  group_by(time_issued, time_forecast, meas_east, meas_north) %>%
  summarise(
    mu_east = mean(fore_east),
    mu_north = mean(fore_north)
  ) %>%
  ungroup()

# The forecast means
mos_mean
#> # A tibble: 56 × 6
#>    time_issued         time_forecast       meas_east meas_north mu_east mu_north
#>    <dttm>              <dttm>                  <dbl>      <dbl>   <dbl>    <dbl>
#>  1 2018-07-06 00:00:00 2018-07-06 00:00:00 -6.58e- 2      -3.77   1.40    -2.55 
#>  2 2018-07-06 00:00:00 2018-07-06 01:00:00  2.25e- 1      -3.22   1.71    -2.22 
#>  3 2018-07-06 00:00:00 2018-07-06 02:00:00  5.76e- 1      -3.63   2.09    -1.95 
#>  4 2018-07-06 00:00:00 2018-07-06 03:00:00  1.36e+ 0      -1.62   2.53    -1.18 
#>  5 2018-07-06 00:00:00 2018-07-06 04:00:00  4.24e- 1      -1.58   2.56    -0.363
#>  6 2018-07-06 00:00:00 2018-07-06 05:00:00  8.14e-17       1.33   2.33     0.152
#>  7 2018-07-06 00:00:00 2018-07-06 06:00:00  1.09e- 1      -3.11   1.79    -0.137
#>  8 2018-07-06 00:00:00 2018-07-06 07:00:00 -3.78e- 1      -3.08   1.36    -0.400
#>  9 2018-07-06 00:00:00 2018-07-06 08:00:00 -1.24e+ 0      -3.80   0.920   -0.476
#> 10 2018-07-06 00:00:00 2018-07-06 09:00:00 -3.58e- 1      -2.54   0.913   -0.270
#> # … with 46 more rows

# The SEs for each prediction horizon
mos_mean %>%
  mutate(
    SE = (meas_east - mu_east)^2 + (meas_north - mu_north)^2
  )
#> # A tibble: 56 × 7
#>    time_issued         time_forecast       meas_east meas_north mu_east mu_north
#>    <dttm>              <dttm>                  <dbl>      <dbl>   <dbl>    <dbl>
#>  1 2018-07-06 00:00:00 2018-07-06 00:00:00 -6.58e- 2      -3.77   1.40    -2.55 
#>  2 2018-07-06 00:00:00 2018-07-06 01:00:00  2.25e- 1      -3.22   1.71    -2.22 
#>  3 2018-07-06 00:00:00 2018-07-06 02:00:00  5.76e- 1      -3.63   2.09    -1.95 
#>  4 2018-07-06 00:00:00 2018-07-06 03:00:00  1.36e+ 0      -1.62   2.53    -1.18 
#>  5 2018-07-06 00:00:00 2018-07-06 04:00:00  4.24e- 1      -1.58   2.56    -0.363
#>  6 2018-07-06 00:00:00 2018-07-06 05:00:00  8.14e-17       1.33   2.33     0.152
#>  7 2018-07-06 00:00:00 2018-07-06 06:00:00  1.09e- 1      -3.11   1.79    -0.137
#>  8 2018-07-06 00:00:00 2018-07-06 07:00:00 -3.78e- 1      -3.08   1.36    -0.400
#>  9 2018-07-06 00:00:00 2018-07-06 08:00:00 -1.24e+ 0      -3.80   0.920   -0.476
#> 10 2018-07-06 00:00:00 2018-07-06 09:00:00 -3.58e- 1      -2.54   0.913   -0.270
#> # … with 46 more rows, and 1 more variable: SE <dbl>
```

**The Dawid-Sebastiani Score** requires the forecast means and
variances. It requires slightly more calculation (due to estimating a
sample variance matrix). First, we define the DSS loss function. Then we
use <tt>lapply</tt> to loop over each prediction horizon and calculate
DSS.

``` r

# DSS loss function
calc_DSS <- function(obs, mu, Sigma){
  ld <- log(det(Sigma))
  DSS <- ld + t(mu - obs) %*% solve(Sigma) %*% (mu - obs)
  return(as.numeric(DSS))
}

# Define the prediction horizons to loop over (this is the prediction domain)
horizons <- mos_predictions$horizon %>% unique()

lapply(horizons, function(i){
  
  # Subset the predictions at h
  mos_h <- mos_predictions %>%
  filter(horizon == i)

  # Define the observed value
  obs_h <- as.matrix(c(mos_h$meas_east[1], mos_h$meas_north[1]))

  # Collate the forecasts into a matrix
  forecasts_h <- matrix(
    c(mos_h$fore_east, mos_h$fore_north), 
    ncol = 2
  )

  # Calculate the empirical mean and variance
  fmean <- colMeans(forecasts_h) %>% matrix()
  fvar <- var(forecasts_h)

  # Calculate DSS
  DSS_calc <- calc_DSS(obs_h, fmean, fvar)

  #  Collate back into a tibble and return
  mos_h[1,] %>%
    dplyr::select(time_issued, time_forecast, horizon) %>%
    mutate(DSS = DSS_calc)
}) %>%
  bind_rows()
#> # A tibble: 56 × 4
#>    time_issued         time_forecast       horizon   DSS
#>    <dttm>              <dttm>                <dbl> <dbl>
#>  1 2018-07-06 00:00:00 2018-07-06 00:00:00       0  3.66
#>  2 2018-07-06 00:00:00 2018-07-06 01:00:00       1  3.62
#>  3 2018-07-06 00:00:00 2018-07-06 02:00:00       2  4.12
#>  4 2018-07-06 00:00:00 2018-07-06 03:00:00       3  3.01
#>  5 2018-07-06 00:00:00 2018-07-06 04:00:00       4  4.55
#>  6 2018-07-06 00:00:00 2018-07-06 05:00:00       5  4.52
#>  7 2018-07-06 00:00:00 2018-07-06 06:00:00       6  6.08
#>  8 2018-07-06 00:00:00 2018-07-06 07:00:00       7  5.65
#>  9 2018-07-06 00:00:00 2018-07-06 08:00:00       8  7.17
#> 10 2018-07-06 00:00:00 2018-07-06 09:00:00       9  4.55
#> # … with 46 more rows
```

**The continuous ranked probability score** and the **energy score**
require the cumulative distribution function (CDF) of the forecast. When
forecasts are described via samples, Equations 3.4 and 3.6 can be used
to instead calculate CRPS and ES, respectively. In <tt>R</tt>, the
<tt>scoringRules</tt> package may be used to efficiently calculate this
via either <tt>crps\_sample</tt> or <tt>es\_sample</tt>. Note, similar
packages exist in [python](https://pypi.org/project/properscoring/) and
[MatLab](https://www.mathworks.com/matlabcentral/fileexchange/77203-acps-package).

``` r
library(scoringRules)
#> Warning: package 'scoringRules' was built under R version 3.6.2

lapply(horizons, function(i){
  
  # Subset the predictions at h
  mos_h <- mos_predictions %>%
  filter(horizon == i)

  # Define the observed value
  obs_h <- c(mos_h$meas_east[1], mos_h$meas_north[1])

  # Collate the forecasts into a matrix
  forecasts_h <- matrix(
    c(mos_h$fore_east, mos_h$fore_north), 
    ncol = 2
  )

  # Calculate ES
  ES_calc <- es_sample(obs_h, t(forecasts_h))

  #  Collate back into a tibble and return
  mos_h[1,] %>%
    dplyr::select(time_issued, time_forecast, horizon) %>%
    mutate(ES = ES_calc)
}) %>%
  bind_rows()
#> # A tibble: 56 × 4
#>    time_issued         time_forecast       horizon    ES
#>    <dttm>              <dttm>                <dbl> <dbl>
#>  1 2018-07-06 00:00:00 2018-07-06 00:00:00       0 1.28 
#>  2 2018-07-06 00:00:00 2018-07-06 01:00:00       1 1.22 
#>  3 2018-07-06 00:00:00 2018-07-06 02:00:00       2 1.48 
#>  4 2018-07-06 00:00:00 2018-07-06 03:00:00       3 0.952
#>  5 2018-07-06 00:00:00 2018-07-06 04:00:00       4 1.64 
#>  6 2018-07-06 00:00:00 2018-07-06 05:00:00       5 1.72 
#>  7 2018-07-06 00:00:00 2018-07-06 06:00:00       6 2.36 
#>  8 2018-07-06 00:00:00 2018-07-06 07:00:00       7 2.14 
#>  9 2018-07-06 00:00:00 2018-07-06 08:00:00       8 2.80 
#> 10 2018-07-06 00:00:00 2018-07-06 09:00:00       9 1.69 
#> # … with 46 more rows
```

## Case Study of Surface Wind Prediction

We
