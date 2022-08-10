
# install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))

library(cmdstanr)
check_cmdstan_toolchain(fix = TRUE, quiet = TRUE)
#install_cmdstan()


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
  seed = 123, 
  chains = 4, 
  parallel_chains = 4,
  iter_warmup = 500,
  iter_sampling = 500,
  thin = 1,
  refresh = 100 # print update every 500 iters
)
