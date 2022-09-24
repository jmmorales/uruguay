

df = read.csv("betas.csv")

logit_reg_pred = function(x, a, b, increase){
  b_est = b / sd(x)
  a_est = a - mean(x) * b_est
  return(plogis(a_est + b_est * (mean(x) + increase)))
}

go = c(
  "Ammodramus_humeralis",
  "Anthus_correndera",
  "Anthus_furcatus",
  "Anthus_hellmayri",
  "Athene_cunicularia",
  "Emberizoides_ypiranganus",
  "Embernagra_platensis",
  "Gallinago_paraguaiae",
  "Nothura_maculosa",
  "Rhynchotus_rufescens",
  "Sicalis_luteola",
  "Sturnella_superciliaris",
  "Tyrannus_savana",
  "Vanellus_chilensis",
  "Xanthopsar_flavus")

# Endangered Species (global y nacional)
en = c(
  "Xanthopsar_flavus",
  "Xolmis_dominicanus",
  "Spartonoica_maluroides",
  "Limnoctites_rectirostris"
  #"Amblyramphus_holosericeus",
  #"Emberizoides_ypiranganus"
)


res_grass = data.frame(sp = spp, from = NA, to = NA)

for(i in 1:length(spp)){
  id = which(spp == spp[i])
  a = df$mean[which(df$param=="Intercept" & df$spp == spp[id])]
  b = df$mean[which(df$param=="Height" & df$spp == spp[id])]
  res_grass$from[i] = logit_reg_pred(X_data$Altura_p, a, b, 0)
  res_grass$to[i] = logit_reg_pred(X_data$Altura_p, a, b, mean(X_data$Altura_p) * 1.1)
}

res_tussock = data.frame(sp = spp, from = NA, to = NA)

for(i in 1:length(spp)){
  id = which(spp == spp[i])
  a = df$mean[which(df$param=="Intercept" & df$spp == spp[id])]
  b = df$mean[which(df$param=="brush" & df$spp == spp[id])]
  res_tussock$from[i] = logit_reg_pred(X_data$Cob_arb, a, b, 0)
  res_tussock$to[i] = logit_reg_pred(X_data$Cob_arb, a, b, mean(X_data$Cob_paj) * 1.1)
}

which(spp=="Tyrannus_savana")

res_tree = data.frame(sp = spp, from = NA, to = NA)

for(i in 1:length(spp)){
  id = which(spp == spp[i])
  a = df$mean[which(df$param=="Intercept" & df$spp == spp[id])]
  b = df$mean[which(df$param=="tree" & df$spp == spp[id])]
  res_tree$from[i] = logit_reg_pred(X_data$Cob_arb, a, b, 0)
  res_tree$to[i] = logit_reg_pred(X_data$Cob_arb, a, b, mean(X_data$Cob_arb) * 1.1)
}



res_lan = data.frame(sp = spp, from = NA, to = NA)

for(i in 1:length(spp)){
  id = which(spp == spp[i])
  a = df$mean[which(df$param=="Intercept" & df$spp == spp[id])]
  b = df$mean[which(df$param=="Landscape" & df$spp == spp[id])]
  res_lan$from[i] = logit_reg_pred(X_data$Uso_suelo, a, b, 0)
  res_lan$to[i] = logit_reg_pred(X_data$Uso_suelo, a, b, 1)
}


df$mean[which(df$param == "brush" & df$f >= 0.95)]

Vanellus_chilensis
Theristicus_caerulescens
Gallinago_paraguaiae
Colaptes_melanochloros
Colaptes_campestris
Xolmis_cinereus
Molothrus_bonariensis
Xolmis_irupero
Turdus_amaurochalinus
Mimus_saturninus
Falco_sparverius
Furnarius_rufus
Satrapa_icterophrys
Pseudoleistes_virescens
Paroaria_coronata
Athene_cunicularia
Pitangus_sulphuratus
Sturnella_superciliaris
Sicalis_luteola
Sporophila_caerulescens
Embernagra_platensis
Phacellodomus_striaticollis






b_est = fit$coefficients[2] / sd(x)
a_est = fit$coefficients[1] - mean(x) * b_est


df$mean[which(df$param=="Height" & df$spp == spp[id])]/sd(X_data$Altura_p)

plogis(df$mean[which(df$param=="Intercept" & df$spp == spp[id])])

plogis(df$mean[which(df$param=="Intercept" & df$spp == spp[id])]-
         mean(X_data$Altura_p) + 
         df$mean[which(df$param=="Height" & df$spp == spp[id])]/sd(X_data$Altura_p) *
         (mean(X_data$Altura_p) +10))


p_scaled = scale(X_data$Altura_p)
plot(X_data$Altura_p, plogis(df[which(df$param=="Intercept" & df$spp == "Vanellus_chilensis"),4] +
       df[which(df$param=="Height" & df$spp == "Vanellus_chilensis"),4] * p_scaled ))

curve(plogis(df[which(df$param=="Intercept" & df$spp == "Vanellus_chilensis"),4] +
               df[which(df$param=="Height" & df$spp == "Vanellus_chilensis"),4] * x), 
      xlim = c(-1, 11))
