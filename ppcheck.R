# Posterior predictive check ----------------------------------------------

library(mvtnorm)
library(cmdstanr)

load("data_Joaquin_sep22.RData")
load("matrices_Joaquin_sep22.rdata")
Alt = scale(log(X_data$Altura_p))#altura del pasto
Paj = rep(1, length(X_data$Cob_paj)) #  scale(X_data$Cob_paj)#cobertura de pajonal
Paj[which(X_data$Cob_paj==0)] = 0
Arb = scale(sqrt(X_data$Cob_arb) )#cobertura arb?rea
Uso_suelo = rep(NA, nrow(X_data))
Uso_suelo[which(X_data$Uso_suelo == "Campo_natural")] = 0
Uso_suelo[which(X_data$Uso_suelo == "Mejorado")] = 1
summary(factor(Uso_suelo))

#La forma de la funci?n va a ser b0+b1*altura+b2*pajonal+b3*tree+b4*uso_suelo
betas_ = c("Height", "brush", "tree", "land_use")
X = cbind(rep(Alt, n.s), rep(Paj, n.s), rep(Arb, n.s), rep(Uso_suelo, n.s))


load("fitc.RData")

post_bs = fit$draws('betas')

nn = nrow(post_bs) * ncol(post_bs)

pnames = c("grass_height", "tussock", "tree_cover", "land_use", "Intercept", "Landscape")

np = length(pnames)
n.s = 69
# N = nrow(post_bs)
# id = sample(1:N, nn, replace = T)#iteraciones que he sampleado

resp = vector("list", np)
for(p in 1:np)
{
  init = ((p-1)*n.s)+1
  fin = init+(n.s-1)
  
  tmp = as.vector(post_bs[,,init:fin])
  dim(tmp) = c(dim(post_bs)[1] * dim(post_bs)[2], n.s)
  resp[[p]]=  tmp #guardo las betas para las simulaciones
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



#Voy a guardar la ocupancia de cada especie su tasa de ocupaci?n (naive), es decir 
#su presencia y detecci?n


OCC = array(NA, c(length(sigma), n.s))

for(i in 1:length(sigma))# length(id))
{
  print(i/length(sigma))
  #######################################
  ##Simulo par?metros que ocurren a escala de comunidad#####
  
  #Covariaci?n espacial
  K = array(NA, c(n.e, n.e))
  for(e1 in 1:n.e)
  {
    for(e2 in 1:n.e)
    {
      K[e1,e2] =  etasq_[i]*exp(-rhosq_[i]*D[e1,e2]^2)
    }
  }
  gamma = rmvnorm(1, mean = rep(0, n.e), K)
  
  #Intercepto por cada establecimiento
  beta0 = array(NA, c(n.s,n.e))
  for(ss in 1:n.s)
  {
    for(e in 1:n.e)
    {
      tmp = resp[[5]][i,ss]+resp[[6]][i,ss]*XE[e,2]
      beta0[ss,e] = rnorm(1, mean = tmp, sd = sigma[i])
    }
  }
  #Intercepto por a?o
  beta00 = rep(NA, ny)
  for(yy in 1:ny)
  {
    beta00[yy] = rnorm(1, mean = 0, sd = sigma2[i])
  }
  
  
  
  
  #########################################
  #Simulo para cada especie ocupaci?n & detecci?n
  occ = rep(NA, n.s)
  for(s in 1:n.s)
  {
    tmp_Y = array(NA, c(n.su, max(visits)))
    for(n in 1:n.su)
    {
      p_pres = plogis(beta0[s,est_[n]]+ beta00[year_[n]]+ X_[n,1]*resp[[1]][i,s]+ X_[n,2]*resp[[2]][i,s]+
                        X_[n,3]*resp[[3]][i,s]+X_[n,4]*resp[[4]][i,s]+ gamma[est_[n]]) 
      z = rbinom(1, 1, p_pres)#presente o no
      tmp_Y[n,1:visits_[n]] = rbinom(visits_[n], 1, z*det_prob[i,s])#detectado en la visita
    }
    tmp = apply(tmp_Y, 1, function(x) sum(na.omit(x)))
    occ[s] = length(tmp[tmp>0])/n.su
  }#end loop spp
  ############################################
  OCC[i,]  = occ
}#end loop simulations

save(OCC, file = "Occ_PPCc.RData")

#Calcular las tasas de ocupaci?n de las especies
Occ_obs = rep(NA, n.s)
id = which(grepl("Visita", colnames(Y_data)) == TRUE)
for(s in 1:n.s)
{
  tmp = Y_data[which(Y_data$Especie == spp[s]), id]
  tmp2 = apply(tmp, 1, function(x) sum(na.omit(x)))
  Occ_obs[s] = length(which(tmp2 >0))/n.su
}

#Hago la base de datos y el plot

L = apply(OCC, 2, function(x) quantile(x, probs  = 0.025))
U = apply(OCC, 2, function(x) quantile(x, probs  = 0.975))
idx = sort(Occ_obs, index.retur = TRUE)$ix
df = data.frame(spp = as.character(spp[idx]), 
                tags = as.character(tag[idx]), 
                obs = Occ_obs[idx], 
                L = L[idx], 
                U = U[idx])

library(ggplot2)

ggplot(df, aes(x = factor(spp, levels = spp), y = obs)) +
  geom_hline(yintercept = 0, color = "grey", size = 1)+
  geom_linerange(aes(ymin =L , ymax = U), color = "darkgrey")+
  #geom_point(aes(shape=vuln, color=f), size = 4)+
  geom_point(size = 2, stroke = 1, color = "grey", alpha = 0.5) +
  theme_minimal()+
  ylab("Observed occupancy" )+
  xlab("") +
  coord_flip()+
  theme(legend.position = "none") 

# ggplot(df, aes(x = tags, y = obs))+
#   theme_classic()+
#   geom_errorbar(aes(ymin=L, ymax=U), width=.05, size = 1, col = "darkgrey")+
#   geom_point(size = 4, col = "blue")+
#   ylab("Occupancy (naive)")+
#   xlab("")+
#   theme(axis.text.x = element_text(colour = "black", size = 6, angle = 75), 
#         axis.text.y = element_text(colour = "black", size = 12))+
#   theme(axis.line.x = element_line(color="black"),
#         axis.line.y = element_line(color="black"))+
#   theme(axis.title=element_text(size=16))+  
#   theme(axis.text.x = element_text(margin = margin(t = 20)))
# 
# ggsave("PPC_naive_occupancy.png", dpi = 300, scale = 1.5)


mean_pred = apply(OCC, 2, mean)
df2 = data.frame(x = Occ_obs, y = mean_pred)

ggplot(df2, aes(x = x, y =y)) +
  theme_minimal()+
  ylab("Predicted occupancy")+
  xlab("Observed occupancy")+
  # theme(axis.text.x = element_text(colour = "black", size = 14), 
  #       axis.text.y = element_text(colour = "black", size = 14))+
  # theme(axis.line.x = element_line(color="black"),
  #       axis.line.y = element_line(color="black"))+
  theme(axis.title=element_text(size=14))+
  geom_line(aes(x=x, y=x), color = "darkgrey", alpha = 0.5, size = 1, linetype = "solid")+
  geom_point(size = 3, color = "grey", alpha = 0.5) 

ggsave("PPC_naive_occupancy_scatter.png", dpi = 300, scale = 1.5)
