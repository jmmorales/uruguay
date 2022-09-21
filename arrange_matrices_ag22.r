#Remove raptors from previous data -------------------------------------------------
rm(list = ls())
load("data_Joaquin_may20.RData")

# fixing some shit
tmpx = X_data$X[285]
tmpy = X_data$Y[285]
X_data$Y[285] = tmpx
X_data$X[285] = tmpy

spp22 = read.csv("bsp_list_ag22.csv")
colnames(phy)
tmp = which(colnames(phy) %in% spp22$Full_name)
length(tmp) == nrow(spp22)
phy2 = phy[tmp, tmp]
nrow(phy2)
ncol(phy2)

head(uruguay_traits)
tmp = which(uruguay_traits$spp %in% spp22$Full_name)
uruguay_traits2 = uruguay_traits[tmp,]
nrow(uruguay_traits2)

head(Y_data)
tmp = which(Y_data$Especie %in% spp22$Full_name)
Y_data2 = Y_data[tmp,]
length(unique(Y_data2$Especie))

length(unique(X_data$TransectaID))
length(unique(Y_data2$TransectaID))

tr_Y = unique(Y_data2$TransectaID)
tr_X = unique(X_data$TransectaID)

tr_X[which(!tr_X %in% tr_Y)]

tmp = which(tr_X %in% tr_Y)
X_data2 = X_data[tmp,]

length(unique(X_data2$TransectaID))
length(unique(Y_data2$TransectaID))


tag2 = rep(NA, nrow(spp22))
for(i in 1:length(tag2))
{
  tmp = strsplit(spp22$Full_name[i], "_")
  gen_ = toupper(substring(tmp[[1]][1], 1, 4))
  spp_ = toupper(substring(tmp[[1]][2], 1, 4))
  tag2[i] = paste(gen_, spp_, sep = "")
}
length(unique(tag2)) == nrow(spp22)

Y_data = Y_data2
X_data = X_data2
uruguay_traits = uruguay_traits2
phy = phy2
spp = spp22$Full_name
tag = tag2

save(Y_data, X_data, phy, uruguay_traits, spp, tag, file = "data_Joaquin_sep22.RData")



# Arrange matrices for data_analyses --------------------------------------


# Load data and libraries -------------------------------------------------
#rm(list = ls())
library(raster)
load("data_Joaquin_sep22.RData")

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
                          scale((uruguay_traits$Insect)),Greg))
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

# Matriz de distancias entre establecimientos -----------------------------------

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
D = D/10
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
     ny, years, file  = "matrices_Joaquin_sep22.RData")


