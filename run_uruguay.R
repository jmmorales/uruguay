
load("matrices_joaquin_ag22.rdata")

#load("matrices_Joaquin_abril20b1.RData")
#load("data_Joaquin_may20.RData")
library("rstan")
#Sys.setenv(LOCAL_CPPFLAGS = '-march=native')
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())


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


fit <- stan(file = 'model_joaquin_apr22.stan',
            data = stan_dat,
            pars = pars,
            iter = 10000,
            warmup = 5000,
            thin = 3,
            chains = 4)

save(fit, file="fit_uruguay.RData")