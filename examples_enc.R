library(ggplot2)
library(ggridges)
library(cmdstanr)
library(forcats)
library(tidyverse)

load("data_Joaquin_sep22.RData")
load("matrices_Joaquin_sep22.rdata")
#spp = read.csv("bsp_list_ag22.csv")
load("fitc.RData")


# Endangered Species (global y nacional)
en = c(
  "Xanthopsar_flavus",
  "Xolmis_dominicanus",
  "Spartonoica_maluroides",
  "Limnoctites_rectirostris"
  #"Amblyramphus_holosericeus",
  #"Emberizoides_ypiranganus"
)

spx = en

for(i in 1:length(en)){
  tmp = strsplit(en[i], "_")
  spx[i] = paste0(tmp[[1]][1], " ", tmp[[1]][2])
}




logit_reg_pred = function(x, a, b, increase){
  b_est = b / sd(x)
  a_est = a - mean(x) * b_est
  return(plogis(a_est + b_est * (mean(x) + increase)))
}

m_bs = fit$summary("betas")
bm = fit$draws('betas')
pos_bm = NULL
for(i in 1:dim(bm)[2]){
  tmp = as.vector(bm[,i,])
  dim(tmp) = c(dim(bm)[1], dim(bm)[3])
  pos_bm = rbind(pos_bm, tmp)
}

str(pos_bm)

ns = length(spp)
mb1 = m_bs$mean[1:ns]
mb2 = m_bs$mean[(ns+1) : (ns * 2)]
mb3 = m_bs$mean[((ns*2)+1) : (ns * 3)]
mb4 = m_bs$mean[(ns*3+1) : (ns * 4)]
mb5 = m_bs$mean[(ns*4+1) : (ns * 5)]
mb6 = m_bs$mean[(ns*5+1) : (ns * 6)]

b1 = pos_bm[, 1:ns]
b2 = pos_bm[, (ns+1) : (ns * 2)]
b3 = pos_bm[, ((ns*2)+1) : (ns * 3)]
b4 = pos_bm[, (ns*3+1) : (ns * 4)]
b5 = pos_bm[, (ns*4+1) : (ns * 5)]
b6 = pos_bm[, (ns*5+1) : (ns * 6)]

rm(pos_bm)

tmp = which(spp %in% en)

d = NULL
for(i in 1:length(tmp)){
  d = rbind(d, cbind(b1[,tmp[i]], 
                     b2[,tmp[i]],
                     b3[,tmp[i]],
                     b4[,tmp[i]],
                     b5[,tmp[i]],
                     b6[,tmp[i]],
                     rep(spp[tmp[i]], nrow(b1))
  )
  )
}

d = data.frame(b1 = as.numeric(d[,1]), 
               b2 = as.numeric(d[,2]),
               b3 = as.numeric(d[,3]),
               b4 = as.numeric(d[,4]),
               b5 = as.numeric(d[,5]),
               b6 = as.numeric(d[,6]),
               sp = d[,7])

# predictions for valleys
# grass height

rs = numeric(nrow(bm) * ncol(bm))
rs1 = numeric(nrow(bm) * ncol(bm))
res = NULL


ini = (log(mean(X_data$Altura_p)) - mean(log(X_data$Altura_p)))/sd(log(X_data$Altura_p))
incr = (log( mean(X_data$Altura_p) * 2) - mean(log(X_data$Altura_p)) ) / sd(log(X_data$Altura_p))

for(i in 1:length(tmp)){
  a = d$b5[which(d$sp == spp[tmp[i]]) ]
  b = d$b1[which(d$sp == spp[tmp[i]]) ]
  for(j in 1:length(a)){
    rs[j] = plogis(a[j] + b[j] * ini)  #logit_reg_pred(log(X_data$Altura_p), a[j], b[j], 0)
    rs1[j] = plogis(a[j] + b[j] * incr)
      
      # logit_reg_pred(log(X_data$Altura_p), a[j], b[j], 
      #                       mean(log(X_data$Altura_p)) * 1.1)
  }
  res = rbind(res, cbind(c(rs1,rs), rep((rs1-rs),2),
                         c(rep("base", length(a)), 
                           rep("+10%", length(a))),
                         rep(spx[i], length(a) * 2)))
}

res_grass = data.frame(x = as.numeric(res[,1]),
                       dif = as.numeric(res[,2]),
                       inc = res[,3],
                       sp = res[,4])

res_grass %>% 
  mutate(sp = fct_reorder(.f = sp, .x = x, .fun = mean)) %>%
  ggplot(aes(y = sp)) +
  geom_density_ridges(
    aes(x = x, fill = paste(sp, inc)), 
    alpha = 0.8, color = "white", from = 0, to = 0.25
  ) +
  ylab("") +
  xlab("Pr of occurrence") +
  ggtitle("Grass height") +
  theme(axis.text.y =element_text(face="italic"))+
  scale_fill_cyclical(
    breaks = c("+10%", "base"),
    #labels = c(`1980 Indy` = "Indy", `1980 Unionist` = "Unionist"),
    values = c("#fde725ff", "#35b479ff"),
    # name = "Option", guide = "legend"
  ) 

ggsave("grass_changec2_en.pdf",
       dpi = 320,
       width = 4, height = 2.5)

#+  theme_minimal()

# ggplot(res_grass, aes(x = x, y = as.factor(sp), fill = sp)) +
#   geom_density_ridges() +
#   guides(fill = FALSE)
# 

#-----------------

rs = numeric(nrow(bm) * ncol(bm))
rs1 = numeric(nrow(bm) * ncol(bm))
res = NULL

for(i in 1:length(tmp)){
  a = d$b5[which(d$sp == spp[tmp[i]]) ]
  b = d$b2[which(d$sp == spp[tmp[i]]) ]
  for(j in 1:length(a)){
    rs[j] = plogis(a[j])
    rs1[j] = plogis( a[j] + b[j])
  }
  res = rbind(res, cbind(c(rs1,rs),
                         c(rep("base", length(a)), rep("+10%", length(a))), 
                         rep(spx[i], length(a) * 2)))
}

res_tussock = data.frame(x = as.numeric(res[,1]), 
                         inc = res[,2], 
                         sp = res[,3])

res_tussock %>% 
  mutate(sp = fct_reorder(.f = sp, .x = x, .fun = mean)) %>%
  ggplot(aes(y = sp)) +
  geom_density_ridges(
    aes(x = x, fill = paste(sp, inc)), 
    alpha = 0.8, color = "white", from = 0, to = 0.25
  ) +
  ylab("") +
  xlab("Pr of occurrence") +
  ggtitle("Tussock") +
  theme(axis.text.y =element_text(face="italic"))+
  scale_fill_cyclical(
    breaks = c("base", "+10%"),
    #labels = c(`1980 Indy` = "Indy", `1980 Unionist` = "Unionist"),
    values = c("#fde725ff", "#35b479ff"),
    # name = "Option", guide = "legend"
  ) 
ggsave("tussock_changec_en.pdf", 
       dpi = 320,
       width = 4, height = 2.5)

mean(res_tussock$x[which(res_tussock$sp == "Limnoctites rectirostris" 
                         & res_tussock$inc == "base")])
mean(res_tussock$x[which(res_tussock$sp == "Limnoctites rectirostris" 
                         & res_tussock$inc != "base")])


#------

ini = (sqrt(mean(X_data$Cob_arb)) - mean(sqrt(X_data$Cob_arb)))/sd(sqrt(X_data$Cob_arb))
incr = (sqrt( mean(X_data$Cob_arb) * 2) - mean(sqrt(X_data$Cob_arb)) ) / sd(sqrt(X_data$Cob_arb))

rs = numeric(nrow(bm) * ncol(bm))
rs1 = numeric(nrow(bm) * ncol(bm))
res = NULL

for(i in 1:length(tmp)){
  a = d$b5[which(d$sp == spp[tmp[i]]) ]
  b = d$b3[which(d$sp == spp[tmp[i]]) ]
  for(j in 1:length(a)){
    
    rs[j] = plogis(a[j] + b[j] * ini)  #logit_reg_pred(log(X_data$Altura_p), a[j], b[j], 0)
    rs1[j] = plogis(a[j] + b[j] * incr)
    
    # rs[j] = logit_reg_pred(sqrt(X_data$Cob_arb), a[j], b[j], 0)
    # rs1[j] = logit_reg_pred(sqrt(X_data$Cob_arb), a[j], b[j], 
    #                         mean(sqrt(X_data$Cob_arb)) * 1.1)
  }
  res = rbind(res, cbind(c(rs1,rs),
                         c(rep("base", length(a)), rep("+10%", length(a))), 
                         rep(spx[i], length(a) * 2)))
}

res_tree = data.frame(x = as.numeric(res[,1]), 
                      inc = res[,2], 
                      sp = res[,3])

res_tree %>% 
  mutate(sp = fct_reorder(.f = sp, .x = x, .fun = mean)) %>%
  ggplot(aes(y = sp)) +
  geom_density_ridges(
    aes(x = x, fill = paste(sp, inc)), 
    alpha = 0.8, color = "white", from = 0, to = 0.25
  ) +
  ylab("") +
  xlab("Pr of occurrence") +
  ggtitle("Tree cover") +
  theme(axis.text.y =element_text(face="italic"))+
  scale_fill_cyclical(
    breaks = c("base", "+10%"),
    #labels = c(`1980 Indy` = "Indy", `1980 Unionist` = "Unionist"),
    values = c("#fde725ff", "#35b479ff"),
    # name = "Option", guide = "legend"
  ) 

ggsave("tree_changec2_en.pdf", 
       dpi = 320,
       width = 4, height = 2.5)

mean(res_tree$x[which(res_tree$sp == "Xanthopsar flavus" & res_tree$inc == "base")])
mean(res_tree$x[which(res_tree$sp == "Xanthopsar flavus" & res_tree$inc != "base")])
mean(res_tree$x[which(res_tree$sp == "Xolmis dominicanus" & res_tree$inc == "base")])
mean(res_tree$x[which(res_tree$sp == "Xolmis dominicanus" & res_tree$inc != "base")])




#------

# grassland type
rs = numeric(nrow(bm) * ncol(bm))
rs1 = numeric(nrow(bm) * ncol(bm))
res = NULL

for(i in 1:length(tmp)){
  a = d$b5[which(d$sp == spp[tmp[i]]) ]
  b = d$b4[which(d$sp == spp[tmp[i]]) ]
  for(j in 1:length(a)){
    rs[j] = plogis(a[j])
    rs1[j] = plogis(a[j] + b[j])
  }
  res = rbind(res, cbind(c(rs1,rs),
                         c(rep("artificial", 
                               length(a)), 
                           rep("natural", length(a))), 
                         rep(spx[i], length(a) * 2)))
}

res_type = data.frame(x = as.numeric(res[,1]), 
                      inc = res[,2], 
                      sp = res[,3])

res_type %>% 
  mutate(sp = fct_reorder(.f = sp, .x = x, .fun = mean)) %>%
  ggplot(aes(y = sp)) +
  geom_density_ridges(
    aes(x = x, fill = paste(sp, inc)), 
    alpha = 0.8, color = "white", from = 0, to = 0.25
  ) +
  ylab("") +
  xlab("Pr of occurrence") +
  ggtitle("Pasture type") +
  theme(axis.text.y =element_text(face="italic"))+
  scale_fill_cyclical(
    breaks = c("base", "+10%"),
    #labels = c(`1980 Indy` = "Indy", `1980 Unionist` = "Unionist"),
    values = c("#fde725ff", "#35b479ff"),
    # name = "Option", guide = "legend"
  ) 
ggsave("type_changec_en.pdf",
       dpi = 320,
       width = 4, height = 2.5)
