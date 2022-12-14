
setwd("~/Documents/PostDoc/20220218_wind_time_series/github/scoringrules")

library(tidyverse)
library(lubridate)
library(scoringRules)
library(reshape2)

mos_predict <- readRDS("../../inter/mos_predict.RDS")

load(file = "data/wind_meas_val.rda")
load(file = "data/wind_fore_val.rda")

wind_meas_val <- wind_meas_val %>%
  group_by(time) %>%
  summarize(
    ew_meas = mean(ew_meas),
    nw_meas = mean(nw_meas)
  )

times_issued <- wind_fore_val$time_issued %>% unique()
horizons <- wind_fore_val$horizon %>% unique() %>% sort()

length(times_issued)
str(mos_predict)

time_select <- ymd_hms("2018-07-06 00:00:00", tz = "Australia/Perth")

idx <- which(time_select == times_issued)

time <- wind_fore_val %>%
	filter(
		time_issued == time_select,
		time_predict < time_select + days(5)
	) %>% pull(time_predict)

dim(mos_predict$y_predict[idx, , 1:56, ])

east_tmp <- reshape2::melt(
	mos_predict$y_predict[idx, 1, 1:56, ]
) %>% 
	as_tibble() %>%
	rename(horizon_tmp = Var1, sample = Var2, fore_east = value)

north_tmp <- reshape2::melt(
	mos_predict$y_predict[idx, 2, 1:56, ]
) %>% 
	as_tibble() %>%
	rename(horizon_tmp = Var1, sample = Var2, fore_north = value)

mos_mcmc_predict <- left_join(east_tmp, north_tmp) %>%
	mutate(
		horizon = horizons[horizon_tmp],
		time_issued = time_select,
		time_forecast = time_select + hours(horizon)
	) %>%
	left_join(wind_meas_val, by = c("time_forecast" = "time")) %>%
	dplyr::select(
		time_issued, time_forecast, horizon, sample, 
		meas_east = ew_meas, meas_north = nw_meas, fore_east, fore_north
	)

saveRDS(mos_mcmc_predict, "data/mos_predictions.RDS")

# SQUARED ERROR

mos_mean <- mos_mcmc_predict %>%
	group_by(time_issued, time_forecast, meas_east, meas_north) %>%
	summarise(
		mu_east = mean(fore_east),
		mu_north = mean(fore_north)
	) %>%
	ungroup()

mos_mean %>%
	mutate(
		SE = (meas_east - mu_east)^2 + (meas_north - mu_north)^2
	)

# DSS

horizons <- mos_mcmc_predict$horizon %>% unique()

lapply(horizons, function(i){
	tmp <- mos_mcmc_predict %>%
	filter(horizon == i)

	forecasts_tmp <- matrix(c(tmp$fore_east, tmp$fore_north), ncol = 2)

	fmean <- colMeans(forecasts_tmp) %>% matrix()
	fvar <- var(forecasts_tmp)

	obs_tmp <- as.matrix(c(tmp$meas_east[1], tmp$meas_north[1]))

	DSS_calc <- as.numeric(log(det(fvar)) + t(fmean - obs_tmp) %*% solve(fvar) %*% (fmean - obs_tmp))

	tmp[1,] %>%
		dplyr::select(time_issued, time_forecast, horizon) %>%
		mutate(DSS = DSS_calc)
}) %>%
	bind_rows()

# ?es_sample

# lapply(horizons, function(i){
	tmp <- mos_mcmc_predict %>%
	filter(horizon == 1)

	forecasts_tmp <- matrix(c(tmp$fore_east, tmp$fore_north), ncol = 2)

	obs_tmp <- c(tmp$meas_east[1], tmp$meas_north[1])

	es_sample(obs_tmp, t(forecasts_tmp))

	tmp[1,] %>%
		dplyr::select(time_issued, time_forecast, horizon) %>%
		mutate(DSS = DSS_calc)
# }) %>%
# 	bind_rows()







mos_mcmc_predict %>%
	filter(sample < 100) %>%
	ggplot() + 
		geom_line(
			aes(x = time_forecast, y = fore_east, group = sample),
			alpha = 0.1
		) + 
		geom_line(
			aes(x = time_forecast, y = meas_east, group = sample), 
			colour = "blue"
		) +
		theme_bw()


wind_meas_val



mos_fore <- t(apply(mos_predict$y_predict[idx, , 1:56, ], c(1,2), quantile, 0.1))
colnames(mos_fore) <- c("Easting", "Northing")
mos_qmin_tbl <- as_tibble(mos_fore) %>%
	mutate(time_predict = time) %>%
	pivot_longer(-time_predict) %>%
	rename(qmin = value)







numerical_forecasts <- readRDS("data/wind_forecats.RDS")

mos_forecasts <- readRDS("data/")

mos_forecasts <- readRDS("../../inter/mos_predict.RDS")


y_predict <- mos_forecasts$y_predict

n_mcmc <- mos_forecasts$n_mcmc

K <- mos_forecasts$K

horizons <- mos_forecasts$horizons
H <- length(horizons)

times_issued <- numerical_forecasts$time_issued %>% unique()


wind_meas_all <- readRDS("~/Documents/PostDoc/201911_stochvol_mcmc/real_data/metocean_meas.RDS") %>%
  dplyr::select(time, ew_meas, nw_meas)

wind_meas_val <- wind_meas_all %>%
  filter(time > ymd_hms("2018-07-01 00:00:00", tz = "Australia/Perth"))

head(times_issued)

length(times_issued)

str(y_predict)




times_issued[100]

str(mos_forecasts)





load_all()

metra_rms <- calc_error_horizon(
	wind_meas_val, c("ew_meas", "nw_meas"),
	wind_fore_val, c("ew_fore", "nw_fore")
) %>%
	mutate(
		Model = "NWP", metric = "Root SE Skill"
	)

metra_dss <- wind_meas_val %>%
  left_join(wind_fore_val, by = c("time" = "time_predict")) %>%
  mutate(
  	eresid = ew_fore - ew_meas,
  	nresid = nw_fore - nw_meas
  )

dss_cov <- var(as.matrix(dplyr::select(metra_dss, eresid, nresid)), na.rm = TRUE)

metra_dss <- metra_dss %>%
	rowwise() %>%
	mutate(dss = as.numeric(-dmvnorm(matrix(c(eresid, nresid)) %>% t(), sigma = as.matrix(dss_cov), log = TRUE))) %>%
	ungroup() %>%
	dplyr::select(times_issued = time_issued, horizon, dss) %>%
	group_by(horizon) %>%
  summarise(score = mean(dss, na.rm = TRUE)) %>%
	mutate(
		metric = "DSS Skill",
		Model = "NWP"
	)

metra_es <- wind_meas_val %>%
  left_join(wind_fore_val, by = c("time" = "time_predict")) %>%
  rowwise() %>%
  mutate(
  	ES = es_sample(c(ew_meas, nw_meas), matrix(c(ew_fore, nw_fore)))
  ) %>%
  group_by(horizon) %>%
  summarise(score = mean(ES, na.rm = TRUE)) %>%
	mutate(
		metric = "ES Skill",
		Model = "NWP"
	)

metra_error <- bind_rows(metra_rms, metra_dss, metra_es)

mos_diagnostics <- readRDS("../inter/mos_diagnostics.RDS") %>%
	mutate(
		metric = ifelse(metric == "RMSE", "Root SE Skill", metric),
		metric = ifelse(metric == "ES", "ES Skill", metric),
		metric = ifelse(metric == "DSS", "DSS Skill", metric)
	) %>%
	rename(Model = model)

var_diagnostics <- readRDS("../inter/var_diagnostics.RDS") %>%
	bind_rows(
		tibble(
			horizon = c(0,0), 
			score = c(0,0), 
			metric = c("RMSE", "ES"), 
			model = c("VAR", "VAR")
		)
	) %>%
	mutate(
		metric = ifelse(metric == "RMSE", "Root SE Skill", metric),
		metric = ifelse(metric == "ES", "ES Skill", metric),
		metric = ifelse(metric == "DSS", "DSS Skill", metric)
	) %>%
	rename(Model = model)

blank_tbl <- tibble(
	y = c(0, 0, 3, 3.5, 3, 5), 
	x = rep(0, 6), 
	metric = c(
		"Root SE Skill", "ES Skill", "DSS Skill", 
		"Root SE Skill", "ES Skill", "DSS Skill"
	)
)

diag_plot <- bind_rows(
	metra_error,
	mos_diagnostics,
	var_diagnostics
) %>%
	ggplot() +
		geom_line(aes(x = horizon, y = score, linetype = Model)) +
		geom_blank(data = blank_tbl, aes(x = x, y = y)) +
		scale_y_continuous("Model Skill [m/s]", expand = c(0, 0)) +
		scale_x_continuous("Prediction Horizon [hrs]", expand = c(0, 0), limits = c(0, 48)) +
		facet_wrap(~metric, ncol = 1, scales = "free") +
		theme_bw() +
		theme(
			legend.position = "bottom", 
			axis.title = element_text(size = 8),
			axis.text = element_text(size = 8),
			legend.title = element_text(size = 8),
			legend.text = element_text(size = 8),
			strip.text = element_text(size = 8)
		)

ggsave(
	"../latex/overleaf/images/model_horizons.png", diag_plot,
	width = 15, height = 12, units = "cm", dpi = 300
)	
