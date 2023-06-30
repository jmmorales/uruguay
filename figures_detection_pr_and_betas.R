library(ggridges)
library(cmdstanr)
library(sfsmisc)
library(ggplot2)
library(ggpubr)
library(coda)
# library(extrafont)
# font_import()
# loadfonts()

load("data_Joaquin_sep22.RData")
load("matrices_Joaquin_sep22.rdata")
load("fitc.RData")

pars <- c("Omega", "tau", "betas", "rho", "z",
          "taup", "ps", "zp", "rhop", "sigmae",
          "rhosq", "etasq", "sigmaee", "pa")


fit_summary = fit$summary(pars)

# Vulnerable spp ----------------------------------------------------------

vuln = c("Limnoctites_rectirostris", 
         "Spartonoica_maluroides", 
         "Xanthopsar_flavus", 
         "Xolmis_dominicanus")
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

VUL = rep("a", length(spp))
GE = rep(0, length(spp))
GE[which(spp %in% go)] = 1
VUL[which(spp %in% go)] = "b"
VUL[which(spp %in% vuln)] = "c"


bs <- fit_summary[grepl("betas", fit_summary$variable),]
post_bs = fit$draws('betas')
b = as.vector(post_bs[,,])
dim(b) = c(dim(post_bs)[1] * dim(post_bs)[2], dim(post_bs)[3])


##Los betas vienen ordenados as? b1sp1...b1spj,,, bpsp1...bpspj
np = ncol(X)+ ncol(XE)
pnames = c("Height", 
           "brush",
           "tree", 
           "land_use", 
           "Intercept",
           "Landscape")
np = length(pnames)
ylabs = c("Grass height", 
          "Tussock", 
          "Tree cover", 
          "Pasture type", 
          "Intercept", 
          "Sierras")


colorsf <- c("a" = "lightgrey", "b" = "darkgrey", "c" = "black")
colors <- c("a" = "white", "b" = "white", "c" = "black")
shapes <- c("a" = 16, "b" = 17, "c" = 18)


# Probability of detection ------------------------------------------------

pdet <- fit_summary[grepl("ps", fit_summary$variable),]
post_det = fit$draws('ps')
det = as.vector(post_det[,,])
dim(det) = c(dim(post_det)[1] * dim(post_det)[2], dim(post_det)[3])

idx = sort(pdet$mean, index.return = TRUE)$ix
spx = spp

for(i in 1:length(idx)){
  aca = strsplit(spp[i], "_")
  spx[i] = paste0(aca[[1]][1], " ", aca[[1]][2])
}

cri = HPDinterval(as.mcmc(plogis(det)))
cri1 = HPDinterval(as.mcmc(plogis(det)), prob = 0.8)

prob_det = data.frame(spp = as.character(spx[idx]),
                      vuln = (VUL[idx]),
                      probs = round(plogis(pdet$mean[idx]), digits = 3),
                      #mean = pdet$mean,
                      Lo = round((cri[idx,1]), digits = 3), 
                      Hi = round((cri[idx,2]), digits = 3),
                      L = cri1[idx,1],
                      U = cri1[idx,2],
                      ess_bulk = round(pdet$ess_bulk[idx], digits = 0),
                      ess_tail = round(pdet$ess_tail[idx], digits = 0),
                      Rhat = round(pdet$rhat[idx], digits = 3))

# write.csv(prob_det, "detection_per_spp_c.csv", row.names = F)

plot_det = 
  ggplot(prob_det, aes(x = factor(spp, levels = spp), y = probs,
                 color = vuln, fill = vuln, shape = vuln)) +
  geom_hline(yintercept = 0, color = "grey", size = 1) +
  geom_linerange(aes(ymin =L , ymax = U, size = 0.8), color = "grey", alpha = 0.5)+
  geom_linerange(aes(ymin =Lo , ymax = Hi), color = "black")+
  #geom_point(aes(shape=vuln, color=f), size = 4)+
  geom_point(size = 3, stroke = 1) +
  scale_color_manual(values = c("#453781FF","#35b479ff", "#fde725ff")) +
  scale_fill_manual(values = c("lightgrey", "darkgrey", "black")) +
  theme_minimal() +
  ylab("Detection probability") +
  xlab("") +
  theme(axis.text.y =element_text(face="italic")) +
  coord_flip()+
  theme(legend.position = "none") 

plot_det

ggsave("detection.pdf", width = 8, height = 15, dpi = 320,
       scale = 1.4, bg = "white", units = "cm")
#ggsave("detectionc.png", width = 4, height = 6)

ps = fit$draws('ps')
pos_ps = NULL
for(i in 1:dim(ps)[2]){
  tmp = as.vector(ps[,i,])
  dim(tmp) = c(dim(ps)[1], dim(ps)[3])
  pos_ps = rbind(pos_ps, tmp)
}

d = NULL
for(i in 1:ncol(pos_ps)) {
  d = rbind(d, cbind(pos_ps[,i],
                     rep(spx[i], nrow(pos_ps))
  )
  )
}

d = data.frame(x = plogis(as.numeric(d[,1])), sp = d[,2])
ggplot(d, aes(x = x, y = as.factor(sp), fill = sp)) +
  geom_density_ridges() +
  guides(fill = FALSE)


# plot betas -----------------------------------------------------------------
plot_list = vector("list", length(ylabs))

for(p in 1:length(ylabs))
{
  init = ((p-1)*n.s)+1
  fin = init+(n.s-1)
  tmp_b = data.frame(bs[init:fin,])
  tmp_bs = b[,init:fin]
  cri = HPDinterval(as.mcmc(tmp_bs), prob = 0.95)
  cri1 = HPDinterval(as.mcmc(tmp_bs), prob = 0.8)
  idx = sort(tmp_b$mean, index.retur = TRUE)$ix
  df = data.frame(spp= as.character(spx[idx]), 
                  y = tmp_b$mean[idx], 
                  L = cri[idx,1],
                  U = cri[idx,2],
                  L1 = cri1[idx,1],
                  U1 = cri1[idx,2],
                  vuln = (VUL[idx]))
  
  plot_b = ggplot(df, aes(x = factor(spp, levels = spp), y = y,
                          color = vuln, fill = vuln, shape = vuln)) +
    geom_hline(yintercept = 0, color = "grey", size = 1)+
    geom_linerange(aes(ymin =L1 , ymax = U1, size = 0.8), color = "grey", alpha = 0.5)+
    geom_linerange(aes(ymin =L , ymax = U), color = "black")+
    #geom_point(aes(shape=vuln, color=f), size = 4)+
    geom_point(size = 3, stroke = 1) +
    scale_color_manual(values = c("#453781FF","#35b479ff", "#fde725ff")) +
    scale_fill_manual(values = c("lightgrey", "darkgrey", "black")) +
    theme_minimal()+
    ylab(as.character(ylabs[[p]])) +
    xlab("") +
    theme(axis.text.y =element_text(face="italic")) +
    coord_flip() +
    theme(legend.position = "none") 
  
  plot_list[[p]] =  plot_b  
}


plot_list[[1]] 
ggsave("grassc.pdf", width = 8, height = 15, dpi = 320,
       scale = 1.4, bg = "white", units = "cm")

plot_list[[2]] 
ggsave("tussockc.pdf", width = 8, height = 15, dpi = 320,
       scale = 1.4, bg = "white", units = "cm")

plot_list[[3]] 
ggsave("treec.pdf", width = 8, height = 15, dpi = 320,
       scale = 1.4, bg = "white", units = "cm")

plot_list[[4]] 
ggsave("typec.pdf", width = 8, height = 15, dpi = 320,
       scale = 1.4, bg = "white", units = "cm")

# ggarrange(plot_list[[1]], plot_list[[2]], plot_list[[3]],plot_list[[4]],
#           labels = c("A", "B", "C", "D"),  
#           ncol = 2, nrow = 2, common.legend = F, 
#           font.label = list(size = 12, color = "black", face = "bold"),
#           vjust = 1)
# ggsave("Cov_2.png", dpi = 300, scale = 2)
# 