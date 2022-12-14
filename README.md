
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
Some data is provided alongside, although these data have been
anonymised and are not original data shown in the paper.

## Scoring Rules

To fully assess a probabilistic forecast we advocate three scoring
rules: the squared error (SE), the Dawid-Sebastiani score (DSS), and the
continuous ranked probability score (CRPS). Roughly, the SE assesses the
forecast mean, the DSS assesses the forecast tails, and the CRPS
assesses the main body of the forecast. Definitions and discussions are
in the paper. First we introduce some data and the forecasts of a
probabilistic forecasting model, then we show how the scores may be
calculated.

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
simple calculated by the samples as \(\mu = 1/n \sum_{j} y_j\)

## Case Study of Surface Wind Prediction

<tt>bayeslinear</tt> is currently in development; install the latest
version from [GitHub](https://github.com/) with:
