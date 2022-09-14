
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

# fit <- mod$sample(
#   data = stan_dat, 
#   #seed = 123, 
#   chains = 4, 
#   parallel_chains = 4,
#   output_dir = "C:\\Users\\jm361n\\uruguay",
#   iter_warmup = 2000,
#   iter_sampling = 10000,
#   thin = 10,
#   refresh = 200 # print update every 500 iters
# )

# fit_vb = mod$variational(
#   data = stan_dat, 
#   grad_samples = 100,
#   seed = 123
# )

# fit_summary = fit$summary(pars)
# write.csv(fit_summary, file = "fit_summary.csv")
# save(fit, file = "fit.RData")

# Rhat Neff ---------------------------------------------------------------

op <- par(mfrow = c(1,2))
hist(fit_summary$rhat, main = "R-hat")
hist(fit_summary$ess_bulk + fit_summary$ess_tail, main = "n-eff" )
par(op)

summary(fit_summary$rhat)
summary(fit_summary$ess_bulk + fit_summary$ess_tail)
#--------------------------------------------------------------------------

load("data_Joaquin_ag22.RData")
load("matrices_Joaquin_ag22.rdata")

library(sfsmisc)
library(ggplot2)
#library(rstan)
library(ggpubr)



# Probability of detection ------------------------------------------------

pdet <- fit_summary[grepl("ps", fit_summary$variable),]

prob_det = data.frame(spp = as.character(spp),
                      tag = as.character(tag), 
                      probs = plogis(pdet$mean),
                      mean = pdet$mean,
                      q005 = pdet$q5, 
                      q095 = pdet$q95,
                      ess_bulk = pdet$ess_bulk,
                      ess_tail = pdet$ess_tail,
                      Rhat = pdet$rhat)

write.csv(prob_det, "detection_probability_per_spp.csv", row.names = F)

df = data.frame(x = as.character(prob_det$tag), y = prob_det$probs,
                L = plogis(prob_det$q005), U = plogis(prob_det$q095))

plot_det1 = ggplot(df, aes(x = x, y = y)) +
  geom_linerange(aes(ymin = L, ymax = U), color = "black")+
  geom_point(size = 4) +
  theme_classic()+
  ylab("Prob.det.")+
  xlab("")+
  theme(axis.text.x = element_text(colour = "black", size = 6, angle = 75), 
        axis.text.y = element_text(colour = "black", size = 12))+
  theme(axis.line.x = element_line(color="black"),
        axis.line.y = element_line(color="black"))+
  theme(axis.title=element_text(size=16))+  
  #geom_hline(yintercept = 0, color = "grey", size = 1.2, linetype = "dashed")+
  theme(axis.text.x = element_text(margin = margin(t = 20)))
plot_det1
ggsave("det_prob.png", dpi = 300, scale = 1.5)


# Trait effects on probability of detection -------------------------------
tr_det = colnames(TT_det)
zps <- fit_summary[grepl("zp", fit_summary$variable),]
#calculo las f
tmp = fit$draws('zp')
f = rep(NA, dim(tmp)[3])
for(i in 1:dim(tmp)[3])
{
  tmp1 = density(tmp[,,i], from = min(tmp[,,i])*1.1, to = max(tmp[,,i]*1.1))
  if(min(tmp[,,i])*1.1 > 0)
  {tmp1 = density(tmp[,,i], from = 0, to = max(tmp[,,i]*1.1))}
  dd = data.frame(x = tmp1$x, y = tmp1$y)
  if(mean(tmp[,,i]) > 0)
  {
    f[i] = integrate.xy(dd$x, dd$y, a = 0, b = max(dd$x))
  }else{
    f[i] = integrate.xy(dd$x, dd$y, a = min(dd$x), b = 0)
  }
}


df <- data.frame(x =  tr_det,
                 fz = zps$mean,
                 L = zps$q5,
                 U = zps$q95, 
                 f = f)
df2 = df[-1,]

plot_det2 = ggplot(df2, aes(x = x, y = fz, color=f)) +
  geom_linerange(aes(ymin =L , ymax = U), color = "black")+
  geom_point(size = 4) +
  scale_colour_gradient(low = "lightgrey", high = "black",limits = c(0.5,1.1))+
  theme_classic()+
  ylab("Effects on prob. det")+
  xlab("")+
  theme(axis.text.x = element_text(colour = "black", size = 14, angle = 45), 
        axis.text.y = element_text(colour = "black", size = 12))+
  theme(axis.line.x = element_line(color="black"),
        axis.line.y = element_line(color="black"))+
  theme(axis.title=element_text(size=14))+  
  geom_hline(yintercept = 0, color = "grey", size = 1.2, linetype = "dashed")+
  theme(axis.text.x = element_text(margin = margin(t = 10)))
plot_det2
ggsave("Traits_prob_det.png", dpi = 300, scale = 1.5)

zps$variable = c("imtercept", "antidep", "gregarious", "size", "vocal")
write.csv(zps, "traits_prob_det.csv", row.names = F)

# Rho para la probabilidad de deteccion -----------------------------------
rhop <- fit_summary[grepl("rhop", fit_summary$variable),]

rhop = as.vector(rhop)

plot(density(fit$draws('rhop')), main = "rho detection probability", xlab = "")
write.csv(rhop, "rho_pd.csv")

# Covariate effects (betas) -----------------------------------------------
bs <- fit_summary[grepl("betas", fit_summary$variable),]

post_bs = fit$draws('betas')
##Los betas vienen ordenados as? b1sp1...b1spj,,, bpsp1...bpspj
np = ncol(X)+ ncol(XE)
pnames = c("Height", "brush", "tree", "land_use", "Intercept", "Landscape")
np = length(pnames)
ylabs = c("Pasture height (mean)", "Brush cover (prop)", "Tree cover (prop)", 
          "Antrop. pasture", "Intercept", "Highland")
betas_plot = list()
betas_df = list()
for(p in 1:np)
{
  init = ((p-1)*n.s)+1
  fin = init+(n.s-1)
  tmp_b = data.frame(bs[init:fin,])
  tmp_b$spp = as.character(tag)
  ####calculo las f
  tmp = post_bs[,,init:fin]
  f = rep(NA, dim(tmp)[3])
  for(i in 1:dim(tmp)[3])
  {
    tmp1 = density(tmp[,,i], from = min(tmp[,,i])*1.1, to = max(tmp[,,i]*1.1))
    if(min(tmp[,,i])*1.1 > 0)
    {tmp1 = density(tmp[,,i], from = 0, to = max(tmp[,,i]*1.1))}
    dd = data.frame(x = tmp1$x, y = tmp1$y)
    if(mean(tmp[,,i]) > 0)
    {
      f[i] = integrate.xy(dd$x, dd$y, a = 0, b = max(dd$x))
    }else{
      if(max(dd$x) < 0)
      {
        f[i] = integrate.xy(dd$x, dd$y, a = min(dd$x), b = max(dd$x)) 
      }else{
        f[i] = integrate.xy(dd$x, dd$y, a = min(dd$x), b = 0)
      }
      
    }
  } 
  
  df = data.frame(x = as.character(tag), 
                  y = tmp_b$mean, 
                  L = tmp_b$q5, 
                  U = tmp_b$q95, 
                  f= f)
  x = ggplot(df, aes(x = x, y = y, color=f)) +
    geom_linerange(aes(ymin =L , ymax = U), color = "black")+
    geom_point(size = 4) +
    scale_colour_gradient(low = "lightgrey", high = "black",limits = c(0.5,1.1))+
    theme_classic()+
    ylab(as.character(ylabs[[p]]))+
    xlab("")+
    theme(axis.text.x = element_text(colour = "black", size = 8, angle =90 ), 
          axis.text.y = element_text(colour = "black", size = 12))+
    theme(axis.line.x = element_line(color="black"),
          axis.line.y = element_line(color="black"))+
    theme(axis.title=element_text(size=14))+  
    geom_hline(yintercept = 0, color = "grey", size = 1.2, linetype = "dashed")+
    theme(axis.text.x = element_text(margin = margin(t = 0)))
  betas_plot[[p]] = x
  betas_df[[p]]= data.frame(param = rep(as.character(pnames[p]), nrow(tmp_b)),
                            spp = as.character(spp), 
                            tag = as.character(tag),
                            mean = tmp_b$mean, 
                            q5 = tmp_b$q5, 
                            q95 = tmp_b$q95,
                            f = f,
                            bulk = tmp_b$ess_bulk, 
                            tail = tmp_b$ess_tail,
                            Rhat = tmp_b$rhat)
}

betas_plot[[1]]
ggsave("Altura_pasto.png", dpi = 300, scale = 2)
betas_plot[[2]]
ggsave("Cobertura_pajonal.png", dpi = 300, scale = 2)
betas_plot[[3]]
ggsave("Cobertura_arborea.png", dpi = 300, scale = 2)
betas_plot[[4]]
ggsave("Mejorado.png", dpi = 300, scale = 2)
betas_plot[[6]]
ggsave("Sierra.png", dpi = 300, scale = 2)


ggarrange(betas_plot[[1]], betas_plot[[2]], betas_plot[[3]], betas_plot[[4]],
          labels = c("A", "B", "C", "D"),  
          ncol = 2, nrow = 2, common.legend = T, 
          font.label = list(size = 12, color = "black", face = "bold"),
          vjust = 1)
ggsave("Covariates_local.png", dpi = 300, scale = 3)
##Guardo las betas
OUT = NULL
for(p in 1:np)
{
  OUT = rbind(OUT, betas_df[[p]])
}
View(OUT)
OUT = data.frame(OUT)
write.csv(OUT, "betas.csv", row.names = F)



# Trait effects en los betas (Z) ----------------------------------------------
zs <- fit_summary[grepl("z", fit_summary$variable),]

tmp = which(grepl("zp",  zs$variable) == TRUE)
zs = zs[-tmp,]#quito las zp de las probabilidades de detecci?n
#View(zs)
post_zs = fit$draws('z')
#En el vector z est? la informaci?n organizada como p1tr1, p1tr2...p1trt, pptr1...pptrt

tr2 =colnames(TT_pres)
nt = length(tr2)
ylabs
zs_plot = list()
zs_df = list()
pnames
for(p in 1:length(pnames))
{
  init = ((p-1)*nt)+1
  fin = init+ (nt -1)
  tmp_z = data.frame(zs[init:fin,])
  tmp_z$traits = as.character(tr2)
  ###Calculo las f
  tmp = post_zs[,,init:fin]
  f = rep(NA, dim(tmp)[3])
  for(i in 1:dim(tmp)[3])
  {
    tmp1 = density(tmp[,,i], from = min(tmp[,,i])*1.1, to = max(tmp[,,i]*1.1))
    if(min(tmp[,,i])*1.1 > 0)
    {tmp1 = density(tmp[,,i], from = 0, to = max(tmp[,,i]*1.1))}
    dd = data.frame(x = tmp1$x, y = tmp1$y)
    if(mean(tmp[,,i]) > 0)
    {
      f[i] = integrate.xy(dd$x, dd$y, a = 0, b = max(dd$x))
    }else{
      if(max(dd$x) < 0)
      {
        f[i] = integrate.xy(dd$x, dd$y, a = min(dd$x), b = max(dd$x)) 
      }else{
        f[i] = integrate.xy(dd$x, dd$y, a = min(dd$x), b = 0)
      }
      
    }
  }
  df = data.frame(x = as.character(tr2),
                  y = tmp_z$mean, 
                  L = tmp_z$q5, 
                  U = tmp_z$q95, 
                  f = f)
  df2 = df[-1,]
  x = ggplot(df2, aes(x = x, y = y, color=f)) +
    geom_linerange(aes(ymin =L , ymax = U), color = "black")+
    geom_point(size = 4) +
    scale_colour_gradient(low = "lightgrey", high = "black",limits = c(0.5,1.1))+
    theme_classic()+
    ylab(as.character(ylabs[[p]]))+
    xlab("")+
    theme(axis.text.x = element_text(colour = "black", size = 14, angle = 45), 
          axis.text.y = element_text(colour = "black", size = 12))+
    theme(axis.line.x = element_line(color="black"),
          axis.line.y = element_line(color="black"))+
    theme(axis.title=element_text(size=14))+  
    geom_hline(yintercept = 0, color = "grey", size = 1.2, linetype = "dashed")+
    theme(axis.text.x = element_text(margin = margin(t = 15)))
  zs_plot[[p]] = x
  colnames(df2) = c("Trait", "Mean", "q5", "q95", "f")
  zs_df[[p]] = df2
}
p = 1
zs_plot[[p]]



zs_plot[[1]]
ggsave("Altura_pasto_traits.png", dpi = 300, scale = 2)
zs_plot[[2]]
ggsave("Cobertura_pajonal_traits.png", dpi = 300, scale = 2)
zs_plot[[3]]
ggsave("Cobertura_arborea_traits.png", dpi = 300, scale = 2)
zs_plot[[4]]
ggsave("Mejorado_traits.png", dpi = 300, scale = 2)

ggarrange(zs_plot[[1]], zs_plot[[2]], zs_plot[[3]], zs_plot[[4]],
          labels = c("A", "B", "C", "D"),  
          ncol = 2, nrow = 2, common.legend = T, 
          font.label = list(size = 12, color = "black", face = "bold"),
          vjust = 1)
ggsave("Covariates_local_traits.png", dpi = 300, scale = 3)


# Scatterplots y boxplots de los traits que han tenido un efecto ----------

betas = read.csv("betas.csv")
uruguay_traits$Greg = as.character(uruguay_traits$Greg)
uruguay_traits$Greg[which(uruguay_traits$Greg == "SI")] = "YES"
df = data.frame(b_height= betas$mean[which(betas$param == "Height")],
                b_brush = betas$mean[which(betas$param == "brush")],
                b_use = betas$mean[which(betas$param == "land_use")],
                b_lands = betas$mean[which(betas$param == "Landscape")],
                ground = scale(uruguay_traits$Ground_stratum), 
                size= log(uruguay_traits$body_size), 
                greg = uruguay_traits$Greg)

plots = vector("list", 4)

plots[[1]] = ggplot(df, aes(x = size, y = b_height)) +
  geom_point(size = 4) +
  theme_classic()+
  ylab("Pasture height response")+
  xlab("Body size (log)")+
  theme(axis.text.x = element_text(colour = "black", size = 14), 
        axis.text.y = element_text(colour = "black", size = 14))+
  theme(axis.line.x = element_line(color="black"),
        axis.line.y = element_line(color="black"))+
  theme(axis.title=element_text(size=16))+
  geom_smooth(method=lm,linetype="dashed",
              color="darkgray", fill="lightgray")

plots[[2]] = ggplot(df, aes(x = ground, y = b_height)) +
  geom_point(size = 4) +
  theme_classic()+
  ylab("Pasture height response")+
  xlab("Ground use (scaled)")+
  theme(axis.text.x = element_text(colour = "black", size = 14), 
        axis.text.y = element_text(colour = "black", size = 14))+
  theme(axis.line.x = element_line(color="black"),
        axis.line.y = element_line(color="black"))+
  theme(axis.title=element_text(size=16))+
  geom_smooth(method=lm,linetype="dashed",
              color="darkgray", fill="lightgray")


plots[[3]] = ggplot(df, aes(x = size, y = b_brush)) +
  geom_point(size = 4) +
  theme_classic()+
  ylab("Brush cover response")+
  xlab("Body size (log)")+
  theme(axis.text.x = element_text(colour = "black", size = 14), 
        axis.text.y = element_text(colour = "black", size = 14))+
  theme(axis.line.x = element_line(color="black"),
        axis.line.y = element_line(color="black"))+
  theme(axis.title=element_text(size=16))+
  geom_smooth(method=lm,linetype="dashed",
              color="darkgray", fill="lightgray")

plots[[4]] = ggplot(df, aes(x=greg, y=b_use)) + 
  geom_boxplot()+
  theme_classic()+
  ylab("Antrop. response")+
  xlab("Gregarious")+
  theme(axis.text.x = element_text(colour = "black", size = 14), 
        axis.text.y = element_text(colour = "black", size = 14))+
  theme(axis.line.x = element_line(color="black"),
        axis.line.y = element_line(color="black"))+
  theme(axis.title=element_text(size=16))

plots[[1]]
ggsave("Altura_pasto_size_scatter.png", dpi = 300, scale = 2)
plots[[2]]
ggsave("Altura_pasto_ground_scatter.png", dpi = 300, scale = 2)
plots[[3]]
ggsave("Pajonal_size_scatter.png", dpi = 300, scale = 2)
plots[[4]]
ggsave("Mejorado_greg_boxplot.png", dpi = 300, scale = 2)

ggarrange(plots[[1]], plots[[2]], plots[[3]], plots[[4]],
          labels = c("A", "B", "C", "D"),  
          ncol = 2, nrow = 2, common.legend = T, 
          font.label = list(size = 12, color = "black", face = "bold"),
          vjust = 1)
ggsave("Traits_betas.png", dpi = 300, scale = 3)

#Efectos de si es monte o sierra
zs_plot[[6]]
ggsave("sierra_traits.png", dpi = 300, scale = 2)


ggplot(df, aes(x = size, y = b_lands)) +
  geom_point(size = 4) +
  theme_classic()+
  ylab("Highlands response")+
  xlab("Body size (log)")+
  theme(axis.text.x = element_text(colour = "black", size = 14), 
        axis.text.y = element_text(colour = "black", size = 14))+
  theme(axis.line.x = element_line(color="black"),
        axis.line.y = element_line(color="black"))+
  theme(axis.title=element_text(size=16))+
  geom_smooth(method=lm,linetype="dashed",
              color="darkgray", fill="lightgray")

ggsave("Highland_size_scatter.png", dpi = 300, scale = 2)


# Guardo la base de datos 


OUT = NULL
for(p in 1:np)
{
  tmp = zs_df[[p]]
  out  = cbind(as.character(pnames[p], nrow(tmp)), tmp)
  OUT = rbind(OUT,out )
}
#View(OUT)
OUT = data.frame(OUT)
write.csv(OUT, "trait_effects_on_betas.csv", row.names= F)

# rho  presencia ----------------------------------------------------------

rho <- fit_summary[which(fit_summary$variable == "rho"),]

plot(density(fit$draws('rho')), main = "rho presence", xlab = "")

write.csv(rho, "rho_presence.csv")


# Autocorrelacion espacial ------------------------------------------------

op = par(mfrow=c(2,2))
plot(density(fit$draws('etasq')), main = "eta")
plot(density(fit$draws('rhosq')), main = "rho_spatial")
plot(density(fit$draws('delta')), main = "delta")
par(op)


#Hago un plot con la autocorrelaci?n
x = seq(from = 0, to = max(D), by = 0.01)

curve(median(fit$draws('etasq'))*exp(-median(fit$draws('rhosq'))*x^2), from = 0, to = 1,
      xlab = "distance 100 km", ylab = "covariance", ylim = c(0,0.3), col = "black",
      lwd = 2)


etas = as.vector(fit$draws('etasq'))
rhosqs = as.vector(fit$draws('rhosq'))

nn =length(etas)

tmp_y = array(NA, c(length(x), nn))

for(i in 1:nn)
{
  tmp_y[,i] = etas[i]*exp(-rhosqs[i]*x^2)
}

mean_cor = apply(tmp_y, 1, mean)
L = apply(tmp_y, 1, function(x) quantile(x, probs = 0.025))
U = apply(tmp_y, 1, function(x) quantile(x, probs = 0.975))

df = data.frame(x = x, y = mean_cor, L = L, U)

ggplot(df, aes(x=x, y=y, L = L, U = U)) + 
  geom_line(size = 1, colour = "black")+
  theme_classic()+
  geom_ribbon(aes(ymin=L, ymax=U), linetype=2, alpha=0.1)+
  ylab("Covariance")+
  xlab("Distance (100 km)")+
  theme(axis.text.x = element_text(colour = "black", size = 12), 
        axis.text.y = element_text(colour = "black", size = 12))+
  theme(axis.line.x = element_line(color="black"),
        axis.line.y = element_line(color="black"))+
  theme(axis.title=element_text(size=14))
ggsave("covariance_distance.png", dpi = 300)


write.csv(fit_summary, "fit_summary.csv")
fit_summary = data.frame(fit_summary)
tmp = which(fit_summary$Rhat > 1.1)
rownames(fit_summary)[tmp]
save.image("analisis_Joaquin.RData")


# Posterior predictive check ----------------------------------------------
library(mvtnorm)


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

load("matrices_Joaquin_ag22.rdata")
#Voy a guardar la ocupancia de cada especie su tasa de ocupaci?n (naive), es decir 
#su presencia y detecci?n


OCC = array(NA, c(nn, n.s))

for(i in 1:nn)# length(id))
{
  #print(i/length(id))
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

save(OCC, file = "Occ_PPC.RData")

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
df = data.frame(spp = as.character(spp), tags = as.character(tag), obs = Occ_obs, L = L, U = U)


ggplot(df, aes(x = tags, y = obs))+
  theme_classic()+
  geom_errorbar(aes(ymin=L, ymax=U), width=.05, size = 1, col = "darkgrey")+
  geom_point(size = 4, col = "blue")+
  ylab("Occupancy (naive)")+
  xlab("")+
  theme(axis.text.x = element_text(colour = "black", size = 6, angle = 75), 
        axis.text.y = element_text(colour = "black", size = 12))+
  theme(axis.line.x = element_line(color="black"),
        axis.line.y = element_line(color="black"))+
  theme(axis.title=element_text(size=16))+  
  theme(axis.text.x = element_text(margin = margin(t = 20)))

ggsave("PPC_naive_occupancy.png", dpi = 300, scale = 1.5)


mean_pred = apply(OCC, 2, mean)
df2 = data.frame(x = Occ_obs, y = mean_pred)

ggplot(df2, aes(x = x, y =y)) +
  geom_point(size = 4) +
  theme_classic()+
  ylab("Occ. predicted (naive)")+
  xlab("Occ. obs (naive)")+
  theme(axis.text.x = element_text(colour = "black", size = 14), 
        axis.text.y = element_text(colour = "black", size = 14))+
  theme(axis.line.x = element_line(color="black"),
        axis.line.y = element_line(color="black"))+
  theme(axis.title=element_text(size=16))+
  geom_line(aes(x=x, y=x), color = "darkred", size = 1.5, linetype = "dashed") 

ggsave("PPC_naive_occupancy_scatter.png", dpi = 300, scale = 1.5)

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

# Vulnerable spp ----------------------------------------------------------

vuln = c("Limnoctites_rectirostris", "Spartonoica_maluroides", "Xanthopsar_flavus", "Xolmis_dominicanus")
VUL = rep(0, length(names_sorted))
VUL[which(names_sorted2 %in% vuln)] = 1

#bs <- fit_summary[grepl("betas", rownames(fit_summary)),]
#post_bs = draws$betas
##Los betas vienen ordenados as? b1sp1...b1spj,,, bpsp1...bpspj
np = ncol(X)+ ncol(XE)
pnames = c("Height", "brush", "tree", "land_use", "Intercept", "Landscape")
np = length(pnames)
ylabs = c("Pasture height (mean)", "Brush cover (prop)", "Tree cover (prop)", 
          "Antrop. pasture", "Intercept", "Highland")
plot_list = vector("list", length(ylabs))


for(p in 1:length(ylabs))
{
  init = ((p-1)*n.s)+1
  fin = init+(n.s-1)
  tmp_b = data.frame(bs[init:fin,])
  tmp_b$spp = as.character(tag)
  ####calculo las f
  tmp = post_bs[,,init:fin]
  f = rep(NA, dim(tmp)[3])
  for(i in 1:dim(tmp)[3])
  {
    tmp1 = density(tmp[,,i], from = min(tmp[,,i])*1.1, to = max(tmp[,,i]*1.1))
    if(min(tmp[,,i])*1.1 > 0)
    {tmp1 = density(tmp[,,i], from = 0, to = max(tmp[,,i]*1.1))}
    dd = data.frame(x = tmp1$x, y = tmp1$y)
    if(mean(tmp[,,i]) > 0)
    {
      f[i] = integrate.xy(dd$x, dd$y, a = 0, b = max(dd$x))
    }else{
      if(max(dd$x) < 0)
      {
        f[i] = integrate.xy(dd$x, dd$y, a = min(dd$x), b = max(dd$x)) 
      }else{
        f[i] = integrate.xy(dd$x, dd$y, a = min(dd$x), b = 0)
      }
      
    }
  }
  
  #reordeno en funci?n de su abundancia
  tmp_mean = tmp_b$mean[id_sorted]
  tmp_L = tmp_b$q5[id_sorted]
  tmp_U = tmp_b$q95[id_sorted]
  tmp_f = f[id_sorted]
  names2 = tolower(names_sorted)
  df = data.frame(x= as.character(names2), 
                  y = tmp_mean, 
                  L = tmp_L,
                  U = tmp_U, 
                  f= tmp_f, 
                  vuln = as.factor(VUL))
  x =  ggplot(df, aes(x = factor(x,levels = names2), y = y, group = vuln, color=f)) +
    geom_linerange(aes(ymin =L , ymax = U), color = "darkgrey")+
    #geom_point(aes(shape=vuln, color=f), size = 4)+
    geom_point(size = 4)+
    theme(axis.text.x = element_text(colour = "black", size = 6, angle =90 ), 
          axis.text.y = element_text(colour = "black", size = 12))+
    scale_colour_gradient(low = "lightgrey", high = "black",limits = c(0.5,1.1))+
    theme_classic()+
    ylab(as.character(ylabs[[p]]))+
    xlab("")+
    theme(axis.text.x = element_text(colour = "black", size = 12), 
          axis.text.y = element_text(colour = "black", size = 9))+
    theme(axis.line.x = element_line(color="black"),
          axis.line.y = element_line(color="black"))+
    theme(axis.title=element_text(size=14))+  
    geom_hline(yintercept = 0, color = "grey", size = 1.2, linetype = "dashed")+
    theme(axis.text.x = element_text(margin = margin(t = 0)))+
    coord_flip()+
    theme(legend.position = "none") 
  plot_list[[p]] =  x+theme(aspect.ratio = 1) 
}

plot_list[[4]]

plot_list[[1]] 
ggsave("Altura_pasto2.png", dpi = 300, scale = 2)
plot_list[[2]] 
ggsave("Cobertura_pajonal2.png", dpi = 300, scale = 2)
plot_list[[3]] 
ggsave("Cobertura_arborea2.png", dpi = 300, scale = 2)
plot_list[[4]] 
ggsave("Mejorado2.png", dpi = 300, scale = 2)

ggarrange(plot_list[[1]], plot_list[[2]], plot_list[[3]],plot_list[[4]],
          labels = c("A", "B", "C", "D"),  
          ncol = 2, nrow = 2, common.legend = F, 
          font.label = list(size = 12, color = "black", face = "bold"),
          vjust = 1)
ggsave("Cov_2.png", dpi = 300, scale = 2)










