# Load data and libraries -------------------------------------------------
rm(list = ls())
library(raster)
load("data_Joaquin_may20.RData")

#Compruebo que los datos est?n en los mismos ?rdenes
length(spp) == ncol(phy)
length(spp) == nrow(phy)
nrow(uruguay_traits) == length(spp)
length(unique(Y_data$TransectaID)) == length(unique(X_data$TransectaID))


# Arrange matrices --------------------------------------------------------

#Lo primero que tengo que hacer es colocar las matrices de datos en el formato que deber?an estar para el modelo 
#que hemos programado
n.s = length(spp)#N?mero de especies
n.su = nrow(X_data)#N?mero de sampling units (transectas)
n.e = length(unique(X_data$Establecimiento))#N?mero de establecimientos
nV = 4 #N?mero de visitas como m?ximo 


##Vectores que identifiquen: las especies (j), las unidades de muestreo (su) y los establecimientos (est)
##vectores que van a servir para identificar las especies de p?jaros,  las transectas y las estancias
j = rep(1:n.s, each= n.su)#repite cada especie un n?mero de veces igual  al n?mero de potreros
su = rep(1:n.su, n.s)#repite las sampling units un n?mero igual al n?mero de especies
j[1:n.su]
su[1:n.su]

tmp = unique(X_data$Establecimiento)#hago una lista de los establecimientos
tag_est = tmp
tmp_e = rep(NA, n.su)
for(i in 1:n.su)
{
  tmp_e[i] = which(tmp == X_data$Establecimiento[i]) 
}
tmp_e
View(cbind(as.character(X_data$Establecimiento), tmp_e))
tag_establecimiento = cbind(as.character(X_data$Establecimiento), tmp_e)
est = rep(tmp_e, n.s)
est[1:n.su]



# Create Y matrix ---------------------------------------------------------




#La matriz y (observaciones) tiene que estar organizada de la siguiente forma [n.s*n.su, nv]. Cada fila
#corresponde a la combinaci?n de una especie con una transecta y en ella hay hasta 4 visitas 1 si se ha observado y 0 si no
#las filas van organizadas de la siguiente manera sp1site1, ..., sp1siten.su, ...., spn.s, siten.su Es decir
#la primera especie las observaciones para todas las transectas,
#la segunda especie las observaciones para todas las transectas...
#y as? sucesivamente

ids = unique(X_data$TransectaID)

ERR = NULL
y  = NULL
tmp2 = grep("Visita", colnames(Y_data))
for(s in 1:n.s)
{
  tmp_y = array(NA, c(length(ids),4))
  for(j in 1:length(ids))
  {
    tmp = Y_data[which(Y_data$Especie == spp[s] &
                         Y_data$TransectaID== ids[j]),]
    if(nrow(tmp) == 0)
    {
      tmp_y[j,] = rep(0, nV)
    }else{
      if(nrow(tmp) >2)
      { err = c(tmp$A.o, tmp$Proyecto, as.character(tmp$Pot), as.character(tmp$Especie))
      ERR = rbind(ERR, err)
      }else{
        tmp_y[j,] = as.numeric(as.vector(tmp[,tmp2]))
      }
    }
  }
  y = rbind(y, tmp_y)
}
ERR

##Compruebo que la matriz se ha contruido bien
ERR = NULL
for(i in 1:length(spp))
{
  init = ((i-1)*n.su)+1
  fin = init + (n.su-1)
  tmp_y = y[init:fin,]
  for(j in 1:length(ids))
  {
    tmp = Y_data[which(Y_data$Especie == as.character(spp[i]) &
                         Y_data$TransectaID == as.character(ids[j])),]
    if(nrow(tmp) == 0)
    {nObs = 0
    }else{
      nObs = sum(na.omit(as.numeric(as.vector(tmp[,tmp2]))))
    }
    nM = sum(na.omit(tmp_y[j,]))
    if(nM != nObs)
    {err = c(i,j)
    ERR =rbind(ERR, err)}
  }
}
ERR#Est? todo OK

tmp2 = grep("Visita", colnames(Y_data))

visits = rep(NA, length(su))#Vector donde voy a ir guardando las visitas

for(i in 1:length(ids))
{
  tmp = Y_data[which(Y_data$TransectaID==  ids[i]),]
  tmp = as.vector(tmp[1,tmp2])
  nv = nV - length(which(is.na(tmp)))
  visits[which(su == i)] = nv
}
which(is.na(visits))

##Tengo que convertir los NA de y en 0, porque stan no permite que haya NA (pero eso no 
#pasa nada porque est? contemplado en el vector visits)
tmp  = which(is.na(y), arr.ind=T)
for(i in 1:nrow(tmp))
{
  y[tmp[i,1], tmp[i,2]] = 0
}

sum(colSums(y))
tmp2 = grep("Visita", colnames(Y_data))
tmp3 = as.vector(Y_data[,tmp2])
tmp3 = apply(tmp3, 2, function(x) sum(na.omit(x)))
sum(tmp3)
rm(tmp)
rm(tmp2)
rm(tmp3)


# X_matrices --------------------------------------------------------------

#Va a haber 2: una que tengan las variables ambientales de las transectas y otra que tenga las variables 
#a escala de establecimiento

#Escala de transecta
Alt = scale(X_data$Altura_p)#altura del pasto
Paj = scale(X_data$Cob_paj)#cobertura de pajonal
Arb = scale(X_data$Cob_arb)#cobertura arb?rea
Uso_suelo = rep(NA, nrow(X_data))
Uso_suelo[which(X_data$Uso_suelo == "Campo_natural")] = 0
Uso_suelo[which(X_data$Uso_suelo == "Mejorado")] = 1
summary(factor(Uso_suelo))

#La forma de la funci?n va a ser b0+b1*altura+b2*pajonal+b3*tree+b4*uso_suelo
betas_ = c("Height", "brush", "tree", "land_use")
X = cbind(rep(Alt, n.s), rep(Paj, n.s), rep(Arb, n.s), rep(Uso_suelo, n.s))#repito un n?mero de veces = n.s para que tenga la misma
#estructura que la matrix Y (sp1site1...sp1siteN,..., spSsite1....spSsiteN)
X = as.matrix(X)
summary(X)
hist(X[,1])
hist(X[,2])
hist(X[,3])
table(X[,4])
nrow(X)/n.su == n.s

#Escala de establecimiento 
tmp = unique(X_data$Establecimiento)#hago una lista de los establecimientos

tmp = unique(X_data$Establecimiento)
Paisaje = rep(NA, length(tmp))
for(i in 1:length(tmp))
{
  Paisaje[i]= as.character(X_data$Unidad_pais[which(X_data$Establecimiento == tmp[i])][1])
}
Paisaje[which(Paisaje == "Valle")] = 0
Paisaje[which(Paisaje == "Sierra")]  = 1
Paisaje = as.numeric(as.character(Paisaje))
XE = cbind(rep(1, n.e), Paisaje)


# Matriz de los traits y la filogenia ----------------------------------------------------

head(uruguay_traits)
tr = c("Size", "Ground", "Insect", "Antidep", "Greg", "Vocal")
Dep = rep(NA, n.s)
Dep[which(uruguay_traits$Antidep == "Camuflaje")] = 0
Dep[which(uruguay_traits$Antidep == "Deteccion")] = 1
Greg = rep(NA, n.s)
Greg[which(uruguay_traits$Greg == "SI")] = 1
Greg[which(uruguay_traits$Greg == "NO")] = 0
Vocal = rep(NA, n.s)
Vocal[which(uruguay_traits$Voice == "BAJA")] = 0
Vocal[which(uruguay_traits$Voice == "ALTA")] = 1
uruguay_traits$body_size = as.numeric(as.character(uruguay_traits$body_size))

TT_pres = as.matrix(cbind(scale(log(uruguay_traits$body_size)), scale(uruguay_traits$Ground),
                          scale(log(uruguay_traits$Insect+1)),Greg))
TT_pres = cbind(rep(1, nrow(TT_pres)), TT_pres)
colnames(TT_pres)= c("Intercept", "Size", "Ground",  "Insect", "Greg")
View(TT_pres)

TT_det = as.matrix(cbind(scale(log(uruguay_traits$body_size)), Dep, Greg, Vocal))
TT_det = cbind(rep(1, nrow(TT_det)), TT_det)
colnames(TT_det) = c("Intercept", "Size", "Antidep", "Greg", "Vocal")


hist(TT_pres[,2])
hist(TT_pres[,3])
hist(TT_pres[,4])



C = as.matrix(phy)#matriz de filogenia


# Matriz de distancias entre transectas -----------------------------------

#Voy a hacer el promedio de las transectas de cada uno de los establecimientos para que
#la parte de autocorrelaci?n espacial sea a una escala m?s grande y m?s manejable (46*46, vs 1079x1079 la matriz de distancias)
coord = array(NA, c(n.e, 2))

for(i in 1:length(tag_est))
{
  tmp = X_data[which(X_data$Establecimiento == tag_est[i]),]
  coord[i,] = c(mean(tmp$X), mean(tmp$Y))
}


D  = array(NA, c(length(tag_est), length(tag_est)))

for(i in 1:length(tag_est))
{
  coord_i = coord[i,]
  for(j in 1:length(tag_est))
  {
    coord_j = coord[j,]
    D[i,j] = pointDistance(coord_i, coord_j, lonlat = F)/1000#En km
  }
}
hist(D)
D = D/100
hist(D)#centenas de km

j = rep(1:n.s, each= n.su)#repite cada especie un n?mero de veces igual  al n?mero de potreros


# A?os --------------------------------------------------------------------

tmp_y = sort(unique(Y_data$A.o))
ny = length(tmp_y)
X_data$Year2 = rep(NA, nrow(X_data))
for(i in 1:length(tmp_y))
{
  tmp = which(X_data$A.o == tmp_y[i])
  X_data$Year2[tmp] = i
}
View(cbind(X_data$A.o, X_data$Year2))

years = rep(X_data$Year2, n.s)

# save input data ---------------------------------------------------------------

save(y, X, XE, TT_det, TT_pres, C, n.su, n.s, su, est, n.e, tag_establecimiento, visits, j, D, 
     ny, years, file  = "matrices_Joaquin_abril20b.RData")


# Data analyses -----------------------------------------------------------

rm(list=ls())
load("matrices_Joaquin_abril20b.RData")
library("rstan")
Sys.setenv(LOCAL_CPPFLAGS = '-march=native')
rstan_options(auto_write = TRUE)
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
          "rhosq", "etasq", "delta", "sigmaee")


fit <- stan(file = 'model_joaquin_may20.stan', 
            data = stan_dat,
            pars = pars, 
            iter = 10, 
            chains = 1)
save.image("fit_10000.RData")


# Analyze posteriors ------------------------------------------------------
rm(list=ls()) 
load("data_Joaquin_may20.RData")
load("matrices_Joaquin_abril20b.RData")
library(sfsmisc)
library(ggplot2)
library(rstan)
library(ggpubr)
load("resJ.RData")

# Rhat Neff ---------------------------------------------------------------


fit_summary <- summary(fit,  probs = c(0.025, 0.05, 0.5, 0.95, 0.975))$summary
op <- par(mfrow = c(1,2))
hist(fit_summary[,10], main = "R-hat")
hist(fit_summary[,9], main = "n-eff" )
par(op)


N <- dim(fit_summary)[[1]]
for (n in 1:N) {
  rhat <- fit_summary[,10][n]
  if (rhat > 1.1 || is.infinite(rhat) || is.nan(rhat)) {
    print(sprintf('Rhat for parameter %s is %s!',
                  rownames(fit_summary)[n], rhat))
  }
}

max(na.omit(fit_summary[,10]))
draws <- rstan::extract(fit)
class(fit)

tmp = as.data.frame(fit_summary)
min(tmp$n_eff[which(!is.na(tmp$n_eff))])
mean(tmp$n_eff[which(!is.na(tmp$n_eff))])
max(tmp$Rhat[which(!is.na(tmp$Rhat))])

# Probability of detection ------------------------------------------------

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

prob_det = read.csv("detection_probability_per_spp.csv")

df = data.frame(x = as.character(prob_det$tag), y = prob_det$probs,
                L = plogis(prob_det$mean.2.5.), U = plogis(prob_det$mean.97.5.))
plot_det1 = ggplot(df, aes(x = x, y = y)) +
  geom_linerange(aes(ymin = df$L, ymax = df$U), color = "black")+
  geom_point(size = 4) +
  theme_classic()+
  ylab("Prob.det.")+
  xlab("")+
  theme(axis.text.x = element_text(colour = "black", size = 6, angle = 75), 
        axis.text.y = element_text(colour = "black", size = 12))+
  theme(axis.line.x = element_line(color="black"),
        axis.line.y = element_line(color="black"))+
  theme(axis.title=element_text(size=16))+  
  geom_hline(yintercept = 0, color = "grey", size = 1.2, linetype = "dashed")+
  theme(axis.text.x = element_text(margin = margin(t = 20)))
plot_det1
ggsave("det_prob.png", dpi = 300, scale = 1.5)


# Trait effects on probability of detection -------------------------------
tr_det = colnames(TT_det)
zps <- fit_summary[grepl("zp", rownames(fit_summary)),]
#calculo las f
tmp = draws$zp
f = rep(NA, ncol(tmp))
for(i in 1:ncol(tmp))
{
  tmp1 = density(tmp[,i], from = min(tmp[,i])*1.1, to = max(tmp[,i]*1.1))
  if(min(tmp[,i])*1.1 > 0)
  {tmp1 = density(tmp[,i], from = 0, to = max(tmp[,i]*1.1))}
  dd = data.frame(x = tmp1$x, y = tmp1$y)
  if(mean(tmp[,i]) > 0)
  {
    f[i] = integrate.xy(dd$x, dd$y, a = 0, b = max(dd$x))
  }else{
    f[i] = integrate.xy(dd$x, dd$y, a = min(dd$x), b = 0)
  }
}


df <- data.frame(x =  tr_det,
                 fz = zps[,1],
                 L = zps[,4],
                 U = zps[,8], 
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


#ggarrange(plot_det1, plot_det2,
#          labels = c("A", "B"),  
#          ncol = 2, nrow = 1,  
#          font.label = list(size = 12, color = "black", face = "bold"),
#          vjust = 1)

head(df2)
df2 = cbind(df2, zps[-1,9], zps[-1,10])
colnames(df2) = c("Traits", "Mean", "q0025", "q0975", "f", "Neff", "Rhat")
write.csv(df2, "traits_prob_det.csv", row.names = F)


# Rho para la probabilidad de deteccion -----------------------------------

rhop <- fit_summary[grepl("rhop", rownames(fit_summary)),]
rhop = as.vector(rhop)

plot(density(draws$rhop), main = "rho detection probability", xlab = "")
rho_pd = data.frame(x = "rho_prob_detec", 
                    mean = rhop[1], q0025 = rhop[4],
                    q0975 = rhop[8], Neff = rhop[9], 
                    Rhat = rhop[10])
write.csv(rho_pd, "rho_pd.csv")



# Covariate effects (betas) -----------------------------------------------

bs <- fit_summary[grepl("betas", rownames(fit_summary)),]
post_bs = draws$betas
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
  tmp = post_bs[,init:fin]
  f = rep(NA, ncol(tmp))
  for(i in 1:ncol(tmp))
  {
    tmp1 = density(tmp[,i], from = min(tmp[,i])*1.1, to = max(tmp[,i]*1.1))
    if(min(tmp[,i])*1.1 > 0)
    {tmp1 = density(tmp[,i], from = 0, to = max(tmp[,i]*1.1))}
    dd = data.frame(x = tmp1$x, y = tmp1$y)
    if(mean(tmp[,i]) > 0)
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
  df = data.frame(x = as.character(tag), y = tmp_b$mean, L = tmp_b$X2.5., U = tmp_b$X97.5., f= f)
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
                              mean = tmp_b[,1], 
                              q0025 = tmp_b[,4], 
                              q0975 = tmp_b[,8],
                              f = f,
                              Neff = tmp_b[,9], 
                              Rhat = tmp_b[,10])
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
zs <- fit_summary[grepl("z", rownames(fit_summary)),]
tmp = which(grepl("zp", rownames(zs)) == TRUE)
zs = zs[-tmp,]#quito las zp de las probabilidades de detecci?n
#View(zs)
post_zs = draws$z
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
  tmp = post_zs[,init:fin]
  f = rep(NA, ncol(tmp))
  for(i in 1:ncol(tmp))
  {
    tmp1 = density(tmp[,i], from = min(tmp[,i])*1.1, to = max(tmp[,i]*1.1))
    if(min(tmp[,i])*1.1 > 0)
    {tmp1 = density(tmp[,i], from = 0, to = max(tmp[,i]*1.1))}
    dd = data.frame(x = tmp1$x, y = tmp1$y)
    if(mean(tmp[,i]) > 0)
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
  df = data.frame(x = as.character(tr2), y = tmp_z$mean, L = tmp_z$X2.5., U = tmp_z$X97.5., f= f)
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
     colnames(df2) = c("Trait", "Mean", "q0025", "q0975", "f")
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
                ground = log(uruguay_traits$Ground_stratum+1), 
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
  xlab("Ground use (log)")+
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
View(OUT)
OUT = data.frame(OUT)
write.csv(OUT, "trait_effects_on_betas.csv", row.names= F)



# rho  presencia ----------------------------------------------------------


rho = fit_summary[grepl("rho", rownames(fit_summary)),]
rho = as.vector(rho)

plot(density(draws$rho), main = "rho presence", xlab = "")
rho_presence = data.frame(x = "rho_presence", 
                          mean = rho[1], q0025 = rho[4],
                          q0975 = rho[8], Neff = rho[9], 
                          Rhat = rho[10])
write.csv(rho_presence, "rho_presence.csv")


# Autocorrelaci?n espacial ------------------------------------------------

op = par(mfrow=c(2,2))
plot(density(draws$etasq), main = "eta")
plot(density(draws$rhosq), main = "rho_spatial")
plot(density(draws$delta), main = "delta")

par(op)


#Hago un plot con la autocorrelaci?n
x = seq(from = 0, to = max(D), by = 0.01)

curve(median(draws$etasq)*exp(-median(draws$rhosq)*x^2), from = 0, to = 1,
      xlab = "distance 100 km", ylab = "covariance", ylim = c(0,0.3), col = "black",
      lwd = 2)



nn = 1000
N = length(draws$etasq)
tmp = sample(1:N, nn, replace = T)
tmp1 = draws$etasq[tmp]
tmp2 = draws$rhosq[tmp]
tmp_y = array(NA, c(length(x), nn))

for(i in 1:nn)
{
  tmp_y[,i] = tmp1[i]*exp(-tmp2[i]*x^2)
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

nn = 1000
post_bs = draws$betas
pnames = c("Height", "brush", "tree", "land_use", "Intercept", "Landscape")
np = length(pnames)
N = nrow(post_bs)
id = sample(1:N, nn, replace = T)#iteraciones que he sampleado

resp = vector("list", np)
for(p in 1:np)
{
  init = ((p-1)*n.s)+1
  fin = init+(n.s-1)
  resp[[p]]=  post_bs[id,init:fin]#guardo las betas para las simulaciones
}

sigma = draws$sigmae[id]#sd de la normal del intercepto (centrado en establecimieto)
sigma2 = draws$sigmaee[id]#sd del efecto a?o

det_prob = plogis(draws$ps[id,])

#par?metros de autocorrelaci?n especial

rhosq_ = draws$rhosq[id]
etasq_ = draws$etasq[id]
delta_ = draws$delta[id]


X_ = X[1:n.su, ]
est_ = est[1:n.su]
visits_ = visits[1:n.su]
year_ = years[1:n.su]
# Simulate responses ------------------------------------------------------
load("matrices_Joaquin_abril20b.RData")

#Voy a guardar la ocupancia de cada especie su tasa de ocupaci?n (naive), es decir 
#su presencia y detecci?n


OCC = array(NA, c(length(id), n.s))


for(i in 1:length(id))# length(id))
{
  print(i/length(id))
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


# New images --------------------------------------------------------------
rm(list=ls()) 
load("data_Joaquin_may20.RData")
load("matrices_Joaquin_abril20b.RData")
library(sfsmisc)
library(ggplot2)
library(rstan)
library(ggpubr)
load("resJ.RData")

fit_summary <- summary(fit,  probs = c(0.025, 0.05, 0.5, 0.95, 0.975))$summary
draws <- rstan::extract(fit)

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

bs <- fit_summary[grepl("betas", rownames(fit_summary)),]
post_bs = draws$betas
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
  tmp = post_bs[,init:fin]
  f = rep(NA, ncol(tmp))
  for(i in 1:ncol(tmp))
  {
    tmp1 = density(tmp[,i], from = min(tmp[,i])*1.1, to = max(tmp[,i]*1.1))
    if(min(tmp[,i])*1.1 > 0)
    {tmp1 = density(tmp[,i], from = 0, to = max(tmp[,i]*1.1))}
    dd = data.frame(x = tmp1$x, y = tmp1$y)
    if(mean(tmp[,i]) > 0)
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
  tmp_L = tmp_b$X2.5.[id_sorted]
  tmp_U = tmp_b$X97.5.[id_sorted]
  tmp_f = f[id_sorted]
  names2 = tolower(names_sorted)
  df = data.frame(x= as.character(names2), y = tmp_mean, L = tmp_L, U = tmp_U, f= tmp_f, vuln = as.factor(VUL))
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



# Trait effects on responses ----------------------------------------------
plots = vector("list",6)

#Las zs est?n organizadas p1t1, p1t2,...,p1tt --- ppt1,...,pptt
betas_ = c("Height", "brush", "tree", "land_use", "intercept", "landscape")
traits = c("Intercept", "Size", "Ground", "Insect", "Greg")


betas = read.csv("betas.csv")
uruguay_traits$Greg = as.character(uruguay_traits$Greg)
uruguay_traits$Greg[which(uruguay_traits$Greg == "SI")] = "YES"
df = data.frame(b_height= betas$mean[which(betas$param == "Height")],
                b_brush = betas$mean[which(betas$param == "brush")],
                b_use = betas$mean[which(betas$param == "land_use")],
                b_lands = betas$mean[which(betas$param == "Landscape")],
                b_tree = betas$mean[which(betas$param == "tree")],
                ground = log(uruguay_traits$Ground_stratum+1), 
                size= log(uruguay_traits$body_size), 
                greg = uruguay_traits$Greg)

#Pasture height
#Sampleo las iteraciones con las que voy a hacer las curvas de respuesta
post_zs = draws$z
niter = 100000
N = nrow(post_zs)
id = sample(1:N, niter, replace = T)#iteraciones que he sampleado
tmp_z = post_zs[id,]

#Matriz de los traits raw
size_raw = log(uruguay_traits$body_size)
ground_raw = uruguay_traits$Ground_stratum
insect_raw = log(uruguay_traits$Insect+1)
TT_raw = cbind(1, size_raw, ground_raw, insect_raw, TT_pres[,5])

#Pasture height ~body  size############################

nt = ncol(TT_pres)
p = 1
init = (p-1)*nt+1
fin = init + (nt-1)
tmp_z2 = tmp_z[,init:fin]

id_col = which(colnames(TT_pres) ==  "Size")
tmp_x = seq(from = min(TT_raw[,id_col]), to= max(TT_raw[,id_col]), length = nrow(TT_pres))
tmp_x = cbind(rep(1, length(tmp_x)), tmp_x)
y_pred = array(NA, c(niter, nrow(tmp_x)))
for(i in 1:niter)
{
  tmp =  c(tmp_z2[i,1], tmp_z2[i,id_col])
  b_est = tmp[2]/sd(TT_raw[,id_col])
  inter = tmp[1] - mean(TT_raw[,id_col]) * b_est
  par_est = c(inter, b_est)
  tmp_TT = cbind(TT_raw[,1], TT_raw[,id_col])
  y_pred [i,] = par_est %*% t(tmp_x)
}

mean_ = colMeans(y_pred)
L_ = apply(y_pred, 2, function(x) quantile(x, probs = 0.025))
U_ = apply(y_pred, 2, function(x) quantile(x, probs = 0.975))
tmp_df = data.frame(x = TT_raw[,id_col], y = df$b_height, mean = mean_, L = L_, U  = U_ , x2 = tmp_x[,2])
plots[[1]] = ggplot(tmp_df, aes(x = x, y = y, mean = mean_, L = L_, U = U_, x2 = x2)) +
  geom_point(size = 3, color = "darkgrey") +
  theme_classic()+
  ylab("Pasture height")+
  xlab("Body size (log)")+
  theme(axis.text.x = element_text(colour = "black", size = 14), 
        axis.text.y = element_text(colour = "black", size = 14))+
  theme(axis.line.x = element_line(color="black"),
        axis.line.y = element_line(color="black"))+
  theme(axis.title=element_text(size=16))+
  geom_line(aes(x=x2, y=mean), color = "darkblue", size = 1.5) +
  geom_line(aes(x=x2, y=L), color = "darkblue", size = 1, linetype = "dashed") +
  geom_line(aes(x=x2, y=U), color = "darkblue", size = 1, linetype = "dashed")


# Pasture height ground use ------------------------------------------------

nt = ncol(TT_pres)
p = 1
init = (p-1)*nt+1
fin = init + (nt-1)
tmp_z2 = tmp_z[,init:fin]

id_col = which(colnames(TT_pres) ==  "Ground")
tmp_x = seq(from = min(TT_raw[,id_col]), to= max(TT_raw[,id_col]), length = nrow(TT_pres))
tmp_x = cbind(rep(1, length(tmp_x)), tmp_x)
y_pred = array(NA, c(niter, nrow(tmp_x)))
for(i in 1:niter)
{
  tmp =  c(tmp_z2[i,1], tmp_z2[i,id_col])
  b_est = tmp[2]/sd(TT_raw[,id_col])
  inter = tmp[1] - mean(TT_raw[,id_col]) * b_est
  par_est = c(inter, b_est)
  tmp_TT = cbind(TT_raw[,1], TT_raw[,id_col])
  y_pred [i,] = par_est %*% t(tmp_x)
}

mean_ = colMeans(y_pred)
L_ = apply(y_pred, 2, function(x) quantile(x, probs = 0.025))
U_ = apply(y_pred, 2, function(x) quantile(x, probs = 0.975))
tmp_df = data.frame(x = TT_raw[,id_col], y = df$b_height, mean = mean_, L = L_, U  = U_ , x2 = tmp_x[,2])
plots[[2]] = ggplot(tmp_df, aes(x = x, y = y, mean = mean_, L = L_, U = U_, x2 = x2)) +
  geom_point(size = 3, color = "darkgrey") +
  theme_classic()+
  ylab("Pasture height")+
  xlab("Ground use (%)")+
  theme(axis.text.x = element_text(colour = "black", size = 14), 
        axis.text.y = element_text(colour = "black", size = 14))+
  theme(axis.line.x = element_line(color="black"),
        axis.line.y = element_line(color="black"))+
  theme(axis.title=element_text(size=16))+
  geom_line(aes(x=x2, y=mean), color = "darkblue", size = 1.5) +
  geom_line(aes(x=x2, y=L), color = "darkblue", size = 1, linetype = "dashed") +
  geom_line(aes(x=x2, y=U), color = "darkblue", size = 1, linetype = "dashed")



# Brush cover -------------------------------------------------------------

nt = ncol(TT_pres)
p = 2
init = (p-1)*nt+1
fin = init + (nt-1)
tmp_z2 = tmp_z[,init:fin]

id_col = which(colnames(TT_pres) ==  "Size")
tmp_x = seq(from = min(TT_raw[,id_col]), to= max(TT_raw[,id_col]), length = nrow(TT_pres))
tmp_x = cbind(rep(1, length(tmp_x)), tmp_x)
y_pred = array(NA, c(niter, nrow(tmp_x)))
for(i in 1:niter)
{
  tmp =  c(tmp_z2[i,1], tmp_z2[i,id_col])
  b_est = tmp[2]/sd(TT_raw[,id_col])
  inter = tmp[1] - mean(TT_raw[,id_col]) * b_est
  par_est = c(inter, b_est)
  tmp_TT = cbind(TT_raw[,1], TT_raw[,id_col])
  y_pred [i,] = par_est %*% t(tmp_x)
}

mean_ = colMeans(y_pred)
L_ = apply(y_pred, 2, function(x) quantile(x, probs = 0.025))
U_ = apply(y_pred, 2, function(x) quantile(x, probs = 0.975))
tmp_df = data.frame(x = TT_raw[,id_col], y = df$b_brush, mean = mean_, L = L_, U  = U_ , x2 = tmp_x[,2])
plots[[3]] = ggplot(tmp_df, aes(x = x, y = y, mean = mean_, L = L_, U = U_, x2 = x2)) +
  geom_point(size = 3, color = "darkgrey") +
  theme_classic()+
  ylab("Brush cover")+
  xlab("Size(log)")+
  theme(axis.text.x = element_text(colour = "black", size = 14), 
        axis.text.y = element_text(colour = "black", size = 14))+
  theme(axis.line.x = element_line(color="black"),
        axis.line.y = element_line(color="black"))+
  theme(axis.title=element_text(size=16))+
  geom_line(aes(x=x2, y=mean), color = "darkblue", size = 1.5) +
  geom_line(aes(x=x2, y=L), color = "darkblue", size = 1, linetype = "dashed") +
  geom_line(aes(x=x2, y=U), color = "darkblue", size = 1, linetype = "dashed")

# Tree cover --------------------------------------------------------------

nt = ncol(TT_pres)
p = 3
init = (p-1)*nt+1
fin = init + (nt-1)
tmp_z2 = tmp_z[,init:fin]

id_col = which(colnames(TT_pres) ==  "Size")
tmp_x = seq(from = min(TT_raw[,id_col]), to= max(TT_raw[,id_col]), length = nrow(TT_pres))
tmp_x = cbind(rep(1, length(tmp_x)), tmp_x)
y_pred = array(NA, c(niter, nrow(tmp_x)))
for(i in 1:niter)
{
  tmp =  c(tmp_z2[i,1], tmp_z2[i,id_col])
  b_est = tmp[2]/sd(TT_raw[,id_col])
  inter = tmp[1] - mean(TT_raw[,id_col]) * b_est
  par_est = c(inter, b_est)
  tmp_TT = cbind(TT_raw[,1], TT_raw[,id_col])
  y_pred [i,] = par_est %*% t(tmp_x)
}

mean_ = colMeans(y_pred)
L_ = apply(y_pred, 2, function(x) quantile(x, probs = 0.025))
U_ = apply(y_pred, 2, function(x) quantile(x, probs = 0.975))
tmp_df = data.frame(x = TT_raw[,id_col], y = df$b_tree, mean = mean_, L = L_, U  = U_ , x2 = tmp_x[,2])
plots[[4]] = ggplot(tmp_df, aes(x = x, y = y, mean = mean_, L = L_, U = U_, x2 = x2)) +
  geom_point(size = 3, color = "darkgrey") +
  theme_classic()+
  ylab("Tree cover")+
  xlab("Size(log)")+
  theme(axis.text.x = element_text(colour = "black", size = 14), 
        axis.text.y = element_text(colour = "black", size = 14))+
  theme(axis.line.x = element_line(color="black"),
        axis.line.y = element_line(color="black"))+
  theme(axis.title=element_text(size=16))+
  geom_line(aes(x=x2, y=mean), color = "darkblue", size = 1.5) +
  geom_line(aes(x=x2, y=L), color = "darkblue", size = 1, linetype = "dashed") +
  geom_line(aes(x=x2, y=U), color = "darkblue", size = 1, linetype = "dashed")


# Tree cover Ground -------------------------------------------------------

id_col = which(colnames(TT_pres) ==  "Ground")
tmp_x = seq(from = min(TT_raw[,id_col]), to= max(TT_raw[,id_col]), length = nrow(TT_pres))
tmp_x = cbind(rep(1, length(tmp_x)), tmp_x)
y_pred = array(NA, c(niter, nrow(tmp_x)))
for(i in 1:niter)
{
  tmp =  c(tmp_z2[i,1], tmp_z2[i,id_col])
  b_est = tmp[2]/sd(TT_raw[,id_col])
  inter = tmp[1] - mean(TT_raw[,id_col]) * b_est
  par_est = c(inter, b_est)
  tmp_TT = cbind(TT_raw[,1], TT_raw[,id_col])
  y_pred [i,] = par_est %*% t(tmp_x)
}

mean_ = colMeans(y_pred)
L_ = apply(y_pred, 2, function(x) quantile(x, probs = 0.025))
U_ = apply(y_pred, 2, function(x) quantile(x, probs = 0.975))
tmp_df = data.frame(x = TT_raw[,id_col], y = df$b_tree, mean = mean_, L = L_, U  = U_ , x2 = tmp_x[,2])
plots[[5]] = ggplot(tmp_df, aes(x = x, y = y, mean = mean_, L = L_, U = U_, x2 = x2)) +
  geom_point(size = 3, color = "darkgrey") +
  theme_classic()+
  ylab("Tree cover")+
  xlab("Size(log)")+
  theme(axis.text.x = element_text(colour = "black", size = 14), 
        axis.text.y = element_text(colour = "black", size = 14))+
  theme(axis.line.x = element_line(color="black"),
        axis.line.y = element_line(color="black"))+
  theme(axis.title=element_text(size=16))+
  geom_line(aes(x=x2, y=mean), color = "darkblue", size = 1.5) +
  geom_line(aes(x=x2, y=L), color = "darkblue", size = 1, linetype = "dashed") +
  geom_line(aes(x=x2, y=U), color = "darkblue", size = 1, linetype = "dashed")


plots[[6]] = ggplot(df, aes(x=greg, y=b_use)) + 
  geom_boxplot()+
  theme_classic()+
  ylab("Antrop.")+
  xlab("Gregarious")+
  theme(axis.text.x = element_text(colour = "black", size = 14), 
        axis.text.y = element_text(colour = "black", size = 14))+
  theme(axis.line.x = element_line(color="black"),
        axis.line.y = element_line(color="black"))+
  theme(axis.title=element_text(size=16))

plots[[1]]
ggsave("Altura_pasto_size_scatter.png", dpi = 300, scale = 1)
plots[[2]]
ggsave("Altura_pasto_ground_scatter.png", dpi = 300, scale = 1)
plots[[3]]
ggsave("Brush_size.png", dpi = 300, scale = 1)
plots[[4]]
ggsave("Tree_size.png", dpi = 300, scale = 1)
plots[[5]]
ggsave("Tree_ground.png", dpi = 300, scale = 1)
plots[[6]]
ggsave("Antr_greg.png", dpi = 300, scale = 1)

plots2 = vector("list", 6)
plots2[[1]] = plots[[1]]+ coord_fixed(2)+ ylab("Pasture height")
plots2[[2]] = plots[[2]]+ coord_fixed(40)+ ylab("Pasture height")
plots2[[3]] = plots[[3]]+ coord_fixed(1.7)+ ylab("Brush cover")
plots2[[4]] = plots[[4]]+ coord_fixed(1.5)+ ylab("Tree cover")
plots2[[5]] = plots[[5]]+ coord_fixed(30)+ ylab("Tree cover")
plots2[[6]] = plots[[6]]+ coord_fixed(1)+ ylab("Antrop.")


plots2 = vector("list", 6)
plots2[[1]] = plots[[1]]+  ylab("Pasture height")
plots2[[2]] = plots[[2]]+  ylab("Pasture height")
plots2[[3]] = plots[[3]]+  ylab("Brush cover")
plots2[[4]] = plots[[4]]+  ylab("Tree cover")
plots2[[5]] = plots[[5]]+  ylab("Tree cover")
plots2[[6]] = plots[[6]]+ ylab("Antrop.")


ggarrange(plots2[[1]], plots2[[2]], plots2[[3]],
          plots2[[4]], plots2[[5]], plots2[[6]],
          labels = c("A", "B", "C", "D", "E", "F", "G"),  
          ncol = 2, nrow = 3, common.legend = T, 
          font.label = list(size = 10, color = "black", face = "bold"),
          vjust = 1)
ggsave("Scatter_trait_effects.png", dpi = 300, scale = 2)

save.image("Analyses_Joaquin.RData")

