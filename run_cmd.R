
install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))

library(cmdstanr)
check_cmdstan_toolchain(fix = TRUE, quiet = TRUE)
install_cmdstan()


load("matrices_Joaquin_ag22.rdata")

#load("matrices_Joaquin_abril20b1.RData")
#load("data_Joaquin_may20.RData")
#library("rstan")
#Sys.setenv(LOCAL_CPPFLAGS = '-march=native')
#rstan_options(auto_write = TRUE)
#options(mc.cores = parallel::detectCores())


stan_dat <- list(
  N = dim(y)[1],
  T = dim(y)[2],
  E = n.e,
  Y = y,
  visits = visits,
  K = dim(X)[2]+ dim(XE)[2],
  N_J = max(j),
  L_J = dim(TT_pres)[2],
  N_1 = n.e,
  J = j,
  J_1 = est,
  Kx = dim(X)[2],
  X = X,
  TT = TT_pres,
  C = C, 
  ones = numeric(max(j)) + 1,
  Kp = 1,
  L_Jp = dim(TT_det)[2], 
  TTp = TT_det, 
  est = est,
  K_e = dim(XE)[2],
  XE = XE,
  Dmat = D,
  years = years,
  NY = max(years))

pars <- c("Omega", "tau", "betas", "rho", "z",
          "taup", "ps", "zp", "rhop", "sigmae",
          "rhosq", "etasq", "sigmaee", "pa")


mod <- cmdstan_model('model_joaquin_apr22.stan')

fit <- mod$sample(
  data = stan_dat, 
  #seed = 123, 
  chains = 4, 
  parallel_chains = 4,
  output_dir = "C:\\Users\\jm361n\\uruguay",
  iter_warmup = 2000,
  iter_sampling = 10000,
  thin = 10,
  refresh = 200 # print update every 500 iters
)

# fit_vb = mod$variational(
#   data = stan_dat, 
#   grad_samples = 100,
#   seed = 123
# )

fit_summary = fit$summary(pars)
write.csv(fit_summary, file = "fit_summary1.csv")
save(fit_summary, file = "fit_summary.RData")

# Rhat Neff ---------------------------------------------------------------

op <- par(mfrow = c(1,2))
hist(fit_summary$rhat, main = "R-hat")
hist(fit_summary$ess_bulk + fit_summary$ess_tail, main = "n-eff" )
par(op)

summary(fit_summary$rhat)
summary(fit_summary$ess_bulk + fit_summary$ess_tail)
#--------------------------------------------------------------------------

#load("data_Joaquin_may20.RData")
load("matrices_Joaquin_ag22.RData")
library(sfsmisc)
library(ggplot2)
library(rstan)
library(ggpubr)



# Probability of detection ------------------------------------------------

pdet <- fit_summary[grepl("ps", fit_summary$variable),]

pdet <- fit_summary[grepl("ps", rownames(fit_summary)),]
prob_det = data.frame(spp = as.character(spp),
                      tag = as.character(tag), 
                      probs = plogis(pdet[,1]),#convierto el fit en plogis porque est? en logit
                      mean = pdet,
                      q005 = pdet[,5], 
                      q095 = pdet[,7],
                      Neff = pdet[,9],
                      Rhat = pdet[,10])
write.csv(prob_det, "detection_probability_per_spp.csv", row.names = F)







fit$summary("rho")
fit$summary("betas")
