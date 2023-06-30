
# install.packages("cmdstanr", 
# repos = c("https://mc-stan.org/r-packages/", getOption("repos")))

library(cmdstanr)
check_cmdstan_toolchain(fix = TRUE, quiet = TRUE)
install_cmdstan()


load("matrices_Joaquin_sep22.rdata")

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
  #output_dir = "C:\\Users\\jm361n\\uruguay",
  iter_warmup = 2000,
  iter_sampling = 10000,
  thin = 10,
  refresh = 200 
)

fit_summary = fit$summary(pars)
write.csv(fit_summary, file = "fit_summary.csv")
save(fit, file = "fit.RData")

# Rhat and Neff ---------------------------------------------------------------

op <- par(mfrow = c(1,2))
hist(fit_summary$rhat, main = "R-hat")
hist(fit_summary$ess_bulk + fit_summary$ess_tail, main = "n-eff" )
par(op)

summary(fit_summary$rhat)
summary(fit_summary$ess_bulk) 
