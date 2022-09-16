# Posterior predictive check ----------------------------------------------

library(mvtnorm)
library(cmdstanr)

load("fit.RData")

post_bs = fit$draws('betas')

nn = nrow(post_bs) * ncol(post_bs)

pnames = c("Height", "brush", "tree", "land_use", "Intercept", "Landscape")
np = length(pnames)
# N = nrow(post_bs)
# id = sample(1:N, nn, replace = T)#iteraciones que he sampleado

resp = vector("list", np)
for(p in 1:np)
{
  init = ((p-1)*n.s)+1
  fin = init+(n.s-1)
  
  tmp = as.vector(post_bs[,,init:fin])
  dim(tmp) = c(dim(post_bs)[1] * dim(post_bs)[2], n.s)
  resp[[p]]=  tmp#guardo las betas para las simulaciones
}

sigma = as.vector(fit$draws('sigmae')) #[id]#sd de la normal del intercepto (centrado en establecimieto)
sigma2 = as.vector(fit$draws('sigmaee')) #[id]#sd del efecto a?o

tmp = as.vector(plogis(fit$draws('ps')))
dim(tmp) = c(dim(post_bs)[1] * dim(post_bs)[2], n.s)
det_prob = tmp

#par?metros de autocorrelaci?n especial

rhosq_ = as.vector(fit$draws('rhosq')) #[id]
etasq_ = as.vector(fit$draws('etasq')) #[id]
delta_ = as.vector(fit$draws('delta')) #[id]


X_ = X[1:n.su, ]
est_ = est[1:n.su]
visits_ = visits[1:n.su]
year_ = years[1:n.su]


# Simulate responses ------------------------------------------------------

#-------------------------------------------------------------------------------

# Betas plot ordered by occupancy (most common, rarest) -------------------

Occ_obs = rep(NA, n.s)
id = which(grepl("Visita", colnames(Y_data)) == TRUE)
for(s in 1:n.s)
{
  tmp = Y_data[which(Y_data$Especie == spp[s]), id]
  tmp2 = apply(tmp, 1, function(x) sum(na.omit(x)))
  Occ_obs[s] = length(which(tmp2 >0))/n.su
}
Occ_sorted = sort(Occ_obs, decreasing = T)
length(Occ_sorted) == n.s

#create a vector identifying spp

id_sorted = c()
names_sorted = c()
names_sorted2 = c()
for(s in 1:n.s)
{
  tmp = which(Occ_obs == Occ_sorted[s])
  id_sorted = c(id_sorted, tmp)
  names_sorted = c(names_sorted, tag[tmp])
  names_sorted2 = c(names_sorted2, spp[tmp])
}
id_sorted = unique(id_sorted)
names_sorted = unique(names_sorted)
names_sorted2 = unique(names_sorted2)