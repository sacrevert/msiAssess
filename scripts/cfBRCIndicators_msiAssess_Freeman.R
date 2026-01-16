## Test msiAssess Freeman (ovm = 4) versus BRCIndicators
#source("scripts/00_newSimsFramework_30102025.R") # indicatorTo2023_msiAssess_ovm4_251203.rds run
source("scripts/00_newSimsFramework_15012026.R") # 1512026 run
#library(devtools)
library(BRCindicators)
#install_github(repo = 'biologicalrecordscentre/BRCindicators')
offDat <- read.csv(file = "data/fionaDat/d4a_spp_ts_engto23_Smoothed_unsmVP.csv")
head(offDat)
offDatSE <- read.csv(file = "data/fionaDat/d4_spp_timeseries_engto23.csv")
head(offDatSE)
head(offDat[offDat$species=="acanthis cabaret",])

# following the DataLabs process
data <- offDat
logSmoothing <- TRUE
if(logSmoothing){
  data$index <- data$index2 # log scale
} else {
  data$index <- data$index1
}
# because no fish at all for 2020 (shouldn't matter for msiAssess, but should check this at some point)
data <- rbind(data[,c("species", "year", "group", "index")], 
              data.frame(species = "sprattus sprattus", 
                         year = 2020,
                         group = "fish",
                         index = NA))
inddat <- data[,c("species","year", "index")]

taxa <- unique(data$group)
#ypk <- c(10, 3) # two smoothing options explored
# I can just look at one for the moment
ypk <- 3L
p <- expand.grid(ypk = ypk, gr = setdiff(taxa, c("moths", "fish")))
n.yr    <- length(unique(inddat$year))
n.knots <- round(n.yr / ypk) # n.knots = 18
spTrends = FALSE
nit = 5000

# ## Note that this function does return Rhats for some parameters given that incl.model = T (default)
# n.chains hardcoded to 3 in bma()
# mod <- BRCindicators::bma(data = inddat, 
#            num.knots = n.knots,
#            parallel = TRUE,
#            plot = FALSE,
#            model = 'smooth', # default
#            seFromData = FALSE, # default
#            Y1perfect = TRUE, # default
#            errorY1 = FALSE, # default
#            n.iter = nit,
#            n.thin = 5, # default
#            incl.model = TRUE, # default
#            save.sppars = spTrends) # FALSE
# saveRDS(mod, file=paste0("data/fionaDat/all_",
#                          ypk,"ypk_L_",
#                          nit/1000,"kit_",
#                          format(Sys.Date(),"%y%m%d"),".rds"))
# mod <- readRDS("data/fionaDat/all_3ypk_L_5kit_251128.rds")
# plot(mod$Year, mod$Index.Mprime) # matches DataLabs plot for ypk = 3
modRhats <- attributes(mod)$model$Rhat # rhats of (hard-coded) parameters (note that general advice is that rhat > 1.1 = not converged/bad)
modq50 <- attributes(mod)$model$q50 # median estimates

## indices for msiAssess
inddat$year <- as.integer(as.character(inddat$year))
msiDat <- reshape(inddat, direction = 'wide', idvar = 'species', timevar = 'year')
useful::topleft(msiDat)
row.names(msiDat) <- msiDat[,1]
msiDat <- msiDat[,-1]
colnames(msiDat) <- as.numeric(sub("index\\.", "", colnames(msiDat)))
msiDat <- msiDat[, order(as.numeric(colnames(msiDat))), drop = FALSE]

# --------------------------------------------------------------------------------------
# test as_sim_from_empirical(), but don't use output as this function is called within
# run_full_analysis() at the beginning
#msiDat_ <- as_sim_from_empirical(as.matrix(msiDat))
# run with 30102025 version first
modMsiAssess <- run_full_analysis(data_source = "empirical", # simulated data or empirical?
                         empirical_index_mat = as.matrix(msiDat),
                         fit_models = c("freeman"),
                         ## Model-specific settings (match BRCIndicators defaults)
                         jags_freeman = list(obs_var_model=4, n_iter=5000, n_thin=5,
                                             n_chains=3, n_burnin=2500), # floor(n_iter/2)) # latter is BRCInd default
                         ## Cross-model settings (match Freeman settings above)
                         brc_opts = list(num_knots = 18, seFromData = FALSE,
                                         Y1perfect = TRUE, m.scale = "loge"),
                         ## Growth-rate presentation scale (model returns both anyway, this is for plot)
                         growth_scale = "log", quiet = FALSE) # suppress JAGS output or not
# saveRDS(modMsiAssess, file = paste0("outputs/indicatorTo2023_msiAssess_ovm4_",
#                          format(Sys.Date(),"%y%m%d"),".rds"))
# modMsiAssess <- readRDS(file = "outputs/indicatorTo2023_msiAssess_ovm4_251203.rds")
modMsiAssess$checks # worth noting that theta has an Rhat of 2.33 here
print(modMsiAssess$plots$MSI) # M_prime
print(modMsiAssess$plots$CumLog) # m
evaluate_inclusion_process(modMsiAssess$sim)
#fyList <- evaluate_inclusion_process(modMsiAssess$sim)$entry_stats$FY
#foo<-environment(modMsiAssess[["fits"]][["freeman"]][["model"]][["data"]])[["data"]][["y1perfect_mask"]]
#apply(foo,1,sum) # numbers of sp with FY in a given year
#plot_species_trends(modMsiAssess, n_sample = 12, scale = "log")

# --------------------------------------------------------------------------------------
# 15-01-2026 v. run with updated theta parameterisation
modMsiAssess_150126 <- run_full_analysis(data_source = "empirical", # simulated data or empirical?
                                  empirical_index_mat = as.matrix(msiDat),
                                  fit_models = c("freeman"),
                                  ## Model-specific settings (match BRCIndicators defaults)
                                  jags_freeman = list(obs_var_model=4, n_iter=5000, n_thin=5,
                                                      n_chains=3, n_burnin=2500), # floor(n_iter/2)) # latter is BRCInd default
                                  ## Cross-model settings (match Freeman settings above)
                                  brc_opts = list(num_knots = 18, seFromData = FALSE,
                                                  Y1perfect = TRUE, m.scale = "loge"),
                                  ## Growth-rate presentation scale (model returns both anyway, this is for plot)
                                  growth_scale = "log", quiet = FALSE) # suppress JAGS output or not
saveRDS(modMsiAssess_150126, file = paste0("outputs/indicatorTo2023_msiAssess_ovm4_",
                         format(Sys.Date(),"%y%m%d"),".rds"))
# modMsiAssess_150126 <- readRDS(file = "outputs/indicatorTo2023_msiAssess_ovm4_260115.rds")
modMsiAssess_150126$checks # theta is 1.15 now
print(modMsiAssess_150126$plots$MSI) # M_prime
print(modMsiAssess_150126$plots$CumLog) # m
evaluate_inclusion_process(modMsiAssess_150126$sim)
#fyList <- evaluate_inclusion_process(modMsiAssess_150126$sim)$entry_stats$FY
#foo2<-environment(modMsiAssess_150126[["fits"]][["freeman"]][["model"]][["data"]])[["data"]][["y1perfect_mask"]]
#apply(foo2,1,sum) # numbers of sp with FY in a given year
#plot_species_trends(modMsiAssess_150126, n_sample = 12, scale = "log")

### Simulated examples
## Run BRCindicators using sim from msiAssess (example 1 in README)
sim_args1 <- list(n_species = 40, n_years = 30, seed = 232680,
                 ## Growth rate DGP
                 dgp_mode = "spline", # "spline","rw","ar1","changepoint","seasonal","mixture"
                 ## DGP options
                 # spline options
                 df_mu = 6, # default
                 # RW / AR1
                 sigma_eta = 0.05, phi = 0.7,
                 # changepoint
                 n_cp = 2, jump_sd = 0.15, piecewise_linear = FALSE,
                 # seasonal
                 A = 0.15, period = 8, phi0 = runif(1, 0, 2*pi),
                 # mixture of random walks
                 K_guilds = 5, # number of groups
                 # species time-series entry mode (guild_staggered for mixture DGP)
                 entry_mode = c("none"),
                 ## Simulate species options
                 # species/state variation
                 sigma_sp = 0.05, sigma_delta = 0.05, sd_alpha0 = 0.4,
                 innov_dist = "normal", df_u = 3, # "normal" or "student_t" (df_u for t)
                 sp_trend = "none", sigma_gamma = 0.05, log_sd_se = 0.35, # "none" or "random_slope"
                 use_delta = FALSE, # Use cross-species annual shocks?
                 # observation model
                 sigma_obs_mean = 0.2,
                 prop_missing = 0, # set to 0 if using approx_a0_for_target_pi() for MNAR
                 ## Inclusion MNAR (additional to prop_missing)
                 # setting use_delta = T means cross-spp shocks enter inclusion selection predictor too
                 inclusion_bias = list(enabled = F, # include inclusion bias?
                                       a0 = -0.4, # lower = P(inc) down at avg growth
                                       a1 = 5.0, # higher = P(inc) sensitive to recent growth
                                       # rho1 - rho0 up = increase stickiness
                                       rho1 = 1.2,
                                       rho0 = 0.1,
                                       p_init = 0.2)) # initial inc prob
out <- run_full_analysis(data_source = "simulate", # simulated data or empirical?
                         sim_args = sim_args1, # as above
                         fit_models = NULL)# sim/empirical outs only
# assume out$sim$y  is S × T  matrix of log indices
Y <- out$sim$y
years <- out$sim$years
species <- rownames(Y) %||% paste0("sp", seq_len(nrow(Y)))
colnames(Y) <- paste0("year.", years)
df_wide <- data.frame(species = species, Y)
df_long <- reshape(df_wide,
                   direction = "long",
                   varying = list(names(df_wide)[-1]),
                   v.names = "index",
                   timevar = "year",
                   times = years)
row.names(df_long) <- NULL
head(df_long)
# run data through BRCindicators
mod <- BRCindicators::bma(data = df_long, 
                          num.knots = 12,
                          parallel = TRUE,
                          plot = FALSE,
                          n.iter = 500)
plot(mod$Year, mod$Index.Mprime)

# second README example
sim_args2 <- list(n_species = 100, n_years = 30, seed = 232680,
                 ## Growth rate DGP
                 dgp_mode = "mixture",
                 ## DGP options for mixture
                 # RW / AR1
                 sigma_eta = 0.05, phi = 0.7,
                 # mixture of random walks
                 K = 5, # number of groups
                 # time-series entry mode
                 entry_mode = c("none"),
                 ## Simulate species options
                 # species/state variation
                 sigma_sp = 0.05, sigma_delta = 0.05, sd_alpha0 = 0.4,
                 innov_dist = "normal", df_u = 3, # "normal" or "student_t" (df_u for t)
                 sp_trend = "random_slope", sigma_gamma = 0.05, log_sd_se = 0.35, # random slopes this time
                 use_delta = FALSE, # cross-species annual shocks?
                 # observation model
                 sigma_obs_mean = 0.05, prop_missing = 0.2,
                 # Inclusion MNAR (additional to prop_missing)
                 inclusion_bias = list(enabled = TRUE, # include inclusion bias?
                                       a0 = -0.5, 
                                       a1 = 3.0,
                                       rho1 = 1.5,
                                       rho0 = 0.1,
                                       p_init = 1))
out2 <- run_full_analysis(data_source = "simulate", # simulated data or empirical?
                         sim_args = sim_args2, # as above
                         fit_models = NULL)
Y <- out2$sim$y
years <- out2$sim$years
species <- rownames(Y) %||% paste0("sp", seq_len(nrow(Y)))
colnames(Y) <- paste0("year.", years)
df_wide <- data.frame(species = species, Y)
df_long2 <- reshape(df_wide,
                   direction = "long",
                   varying = list(names(df_wide)[-1]),
                   v.names = "index",
                   timevar = "year",
                   times = years)
row.names(df_long2) <- NULL
head(df_long2)
# run data through BRCindicators
mod <- BRCindicators::bma(data = df_long2, 
                          num.knots = 12,
                          parallel = TRUE,
                          plot = FALSE,
                          seFromData = FALSE,
                          Y1perfect = TRUE,
                          n.iter = 1000)
plot(mod$Year, mod$Index.Mprime)
plot(mod$Year, mod$Index.M)

out2cf <- run_full_analysis(data_source = "simulate",
                          sim_args = sim_args2, # as above
                          fit_models = "freeman",
                          # same settings as BRCindicator defaults
                          jags_freeman = list(obs_var_model = 4,
                                              n_iter = 500, n_burnin = 100, n_chains = 3, n_thin = 5, seed = NULL),
                          brc_opts = list(num_knots = 12, seFromData = FALSE,
                                          Y1perfect = TRUE, m.scale = "loge"))
print(out2cf$plots$MSI) # M_prime equivalents
print(out2cf$plots$CumLog)