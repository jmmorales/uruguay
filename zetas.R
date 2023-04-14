library(coda)
library(ggplot2)
library(cmdstanr)
library(ggpubr)
library(dplyr)

load("data_Joaquin_sep22.RData")
load("matrices_Joaquin_sep22.rdata")
load("fitc_summary.RData")
load("fitc.RData")

betas = read.csv("betasc.csv")

bms = fit_summary$mean[grepl("betas", fit_summary$variable)]
bms = matrix(bms, 69,6)

uruguay_traits$Greg = as.character(uruguay_traits$Greg)
uruguay_traits$Greg[which(uruguay_traits$Greg == "SI")] = "YES"
df = data.frame(
  b_interecpt = betas$mean[which(betas$param == "Intercept")],
  b_height= betas$mean[which(betas$param == "Height")], 
  f_height = betas$f[which(betas$param == "Height")],
  b_tussock = betas$mean[which(betas$param == "brush")],
  f_tussock = betas$f[which(betas$param == "brush")],
  b_tree = betas$mean[which(betas$param == "tree")],
  f_tree = betas$f[which(betas$param == "tree")],
  b_art_pasture = betas$mean[which(betas$param == "land_use")],
  f_art_pasture = betas$f[which(betas$param == "land_use")],
  a_hill = bms[,5],
  b_hill = betas$mean[which(betas$param == "Landscape")],
  size= log(uruguay_traits$body_size),
  ground = (uruguay_traits$Ground_stratum),
  insect = (uruguay_traits$Insect),
  greg = uruguay_traits$Greg
)

# Trait effects en los betas (Z) ----------------------------------------------
zs <- fit_summary[grepl("z", fit_summary$variable),]
tmp = which(grepl("zp",  zs$variable) == TRUE)
zs = zs[-tmp,] #quito las zp de las probabilidades de detección

nbs = rep(c("grass height", "tussock", "tree cover", "pasture type",
      "valley", "sierras"), each = 5)
n = rep(c("intercept", "size", "ground", "insectivory", "gregarious"), each = 6)

zs$coeff = nbs
zs$pred = n

# get posteriors and rearrange
post_zs = fit$draws('z')
z = as.vector(post_zs[,,])
dim(z) = c(dim(post_zs)[1] * dim(post_zs)[2], dim(post_zs)[3])

# calculate f and means
ff = numeric(ncol(z))
mz = numeric(ncol(z))
for(i in 1:ncol(z)){
  tmp = mean(z[,i])
  ff[i] = sum(sign(z[,i]) == sign(tmp))/nrow(z)
  mz[i] = tmp
}

zs$f = ff

write.csv(zs, file = "zetas.csv")


# size effects
scaled_size = seq(min(TT_pres[,2]), max(TT_pres[,2]), by = 0.1)

s_grass = matrix(NA, nrow(z), length(scaled_size))
s_tuss = matrix(NA, nrow(z), length(scaled_size))
s_tree = matrix(NA, nrow(z), length(scaled_size))
s_type = matrix(NA, nrow(z), length(scaled_size))
s_a = matrix(NA, nrow(z), length(scaled_size))
s_b = matrix(NA, nrow(z), length(scaled_size))

for(i in 1:nrow(z)){
  Z = matrix(z[i,], 5, 6)
  s_grass[i,] = Z[1,1] + scaled_size * Z[2,1]
  s_tuss[i,] = Z[1,2] + scaled_size * Z[2,2]
  s_tree[i,] = Z[1,3] + scaled_size * Z[2,3]
  s_type[i,] = Z[1,4] + scaled_size * Z[2,4]
  s_a[i,] = Z[1,5] + scaled_size * Z[2,5]
  s_b[i,] = Z[1,6] + scaled_size * Z[2,6]
}

cris = HPDinterval(as.mcmc(s_grass))
crits = HPDinterval(as.mcmc(s_tuss))
critr = HPDinterval(as.mcmc(s_tree))
crity= HPDinterval(as.mcmc(s_type))
cria= HPDinterval(as.mcmc(s_a))
crib = HPDinterval(as.mcmc(s_b))

size = exp(scaled_size*sd(log(uruguay_traits$body_size)) + mean(log(uruguay_traits$body_size)))


plot_list = vector("list", 4)

dt = data.frame(x = (size), 
                y = colMeans(s_grass), 
                ylow = cris[,1], 
                yup = cris[,2], 
                variable = 1)

pl = ggplot(dt, aes(x = x, y = y)) +
  geom_ribbon(aes(ymin = ylow, ymax = yup, 
                  fill = as.factor(variable)), 
              alpha = 0.5, fill = "grey70") +
  geom_line(color = "black", size = 2, aes(alpha=0.5)) +
  geom_hline(yintercept = 0) +
  xlab("Bird size") +
  ylab("Grass height coefficient") +
  theme_minimal() +
  ylim(-1.5,3.3) +
  geom_point(data = df, aes(x = exp(size), 
                            y = b_height,
                            color = 1-f_height), size = 2) +
  scale_colour_viridis_c(option = "D") +
#  scale_color_gradient(low = "lightgrey", high = "black") +
  scale_x_log10() +
  theme(legend.position = "none") 

plot_list[[1]] = pl

# tussock --------------------------------
dt = data.frame(x = size, 
                y = colMeans(s_tuss), 
                ylow = crits[,1], 
                yup = crits[,2], 
                variable = 1)

pl = ggplot(dt, aes(x = x, y = y)) +
  geom_ribbon(aes(ymin = ylow, ymax = yup, 
                  fill = as.factor(variable)), 
              alpha = 0.5, fill = "grey70") +
  geom_line(color = "black", size = 2, aes(alpha=0.5)) +
  geom_hline(yintercept = 0) +
  xlab("Bird size") +
  ylab("Tussock coefficient") +
  theme_minimal() +
  ylim(-1.5,3.3) +
  geom_point(data = df, aes(x = exp(size), 
                            y = b_tussock,
                            color = 1-f_tussock
                            ), size = 2) +
  scale_colour_viridis_c(option = "D") +
  scale_x_log10() +
  theme(legend.position = "none") 

plot_list[[2]] = pl

# tree --------------------------------------------------------
dt = data.frame(x = size, 
                y = colMeans(s_tree), 
                ylow = critr[,1], 
                yup = critr[,2], 
                variable = 1)

pl = ggplot(dt, aes(x = x, y = y)) +
  geom_ribbon(aes(ymin = ylow, ymax = yup, 
                  fill = as.factor(variable)), 
              alpha = 0.5, fill = "grey70") +
  geom_line(color = "black", size = 2, aes(alpha=0.5)) +
  geom_hline(yintercept = 0) +
  xlab("Bird size") +
  ylab("Tree coefficient") +
  theme_minimal() +
  ylim(-1.5,3.3) +
  geom_point(data = df, aes(x = exp(size), 
                            y = b_tree,
                            color = 1-f_tree
                            ), size = 2) +
  scale_colour_viridis_c(option = "D") +
  scale_x_log10() +
  theme(legend.position = "none") 

plot_list[[3]] = pl

# pasture type ----------------------------------------
dt = data.frame(x = size, 
                y = colMeans(s_type), 
                ylow = crity[,1], 
                yup = crity[,2], 
                variable = 1)

pl = ggplot(dt, aes(x = x, y = y)) +
  geom_ribbon(aes(ymin = ylow, ymax = yup, 
                  fill = as.factor(variable)), 
              alpha = 0.5, fill = "grey70") +
  geom_line(color = "black", size = 2, aes(alpha=0.5)) +
  geom_hline(yintercept = 0) +
  xlab("Bird size") +
  ylab("Pasture type coefficient") +
  theme_minimal() +
  ylim(-1.5,1) +
  geom_point(data = df, aes(x = exp(size), 
                            y = b_art_pasture,
                            color = 1- f_art_pasture
                            ), size = 2) +
  scale_colour_viridis_c(option = "D") +
  scale_x_log10() +
  theme(legend.position = "none") 

plot_list[[4]] = pl

ggarrange(plot_list[[1]], 
          plot_list[[2]], 
          plot_list[[3]],
          plot_list[[4]],
          labels = c("A", "B", "C", "D"),  
          ncol = 2, nrow = 2, common.legend = F, 
          font.label = list(size = 12, color = "black", face = "bold"),
          vjust = 1)

#ggsave("size.png", width = 3, height = 3, dpi = 320, scale = 2)

size_plots = plot_list

#--------------
#--------------

pl_list = vector("list", 2)
# a ----------------------
dt = data.frame(x = size, 
                y = (colMeans(s_a)), 
                ylow = (cria[,1]), 
                yup = (cria[,2]), 
                variable = 1)

pl = ggplot(dt, aes(x = x, y = y)) +
  geom_ribbon(aes(ymin = ylow, ymax = yup, 
                  fill = as.factor(variable)), 
              alpha = 0.5, fill = "grey70") +
  geom_line(color = "black", size = 2, aes(alpha=0.5)) +
  xlab("Bird size") +
  ylab("Valley") +
  theme_minimal() +
  #ylim(0,1) +
  geom_point(data = df, aes(x = exp(size), 
                            y = (a_hill), 
                            alpha = 0.5), size = 2) +
  scale_x_log10() +
  theme(legend.position = "none") 

pl_list[[1]] = pl


# b ----------------------
dt = data.frame(x = size, 
                y = (colMeans(s_b)), 
                ylow = (crib[,1]), 
                yup = (crib[,2]), 
                variable = 1)

pl = ggplot(dt, aes(x = x, y = y)) +
  geom_ribbon(aes(ymin = ylow, ymax = yup, 
                  fill = as.factor(variable)), 
              alpha = 0.5, fill = "grey70") +
  geom_line(color = "black", size = 2, aes(alpha=0.5)) +
  geom_hline(yintercept = 0) +
  xlab("Bird size") +
  ylab("Sierras") +
  theme_minimal() +
  ylim(-2,1) +
  geom_point(data = df, aes(x = exp(size), 
                            y = (b_hill), 
                            alpha = 0.5), size = 2) +
  scale_x_log10() +
  theme(legend.position = "none") 

pl_list[[2]] = pl

ggarrange(pl_list[[1]], 
          pl_list[[2]], 
          labels = c("A", "B"),  
          ncol = 2, nrow = 1, common.legend = F, 
          font.label = list(size = 12, color = "black", face = "bold"),
          vjust = 1)

#ggsave("size_sierras.png", width = 3, height = 1.7, dpi = 300, scale = 2)

size_sierras = pl_list
#-------------------
#-------------------

scaled_ground = seq(min(TT_pres[,3]), max(TT_pres[,3]), by = 0.1)

g_grass = matrix(NA, nrow(z), length(scaled_ground))
g_tuss = matrix(NA, nrow(z), length(scaled_ground))
g_tree = matrix(NA, nrow(z), length(scaled_ground))
g_type = matrix(NA, nrow(z), length(scaled_ground))
g_a = matrix(NA, nrow(z), length(scaled_ground))
g_b = matrix(NA, nrow(z), length(scaled_ground))

for(i in 1:nrow(z)){
  Z = matrix(z[i,], 5, 6)
  g_grass[i,] = Z[1,1] + scaled_ground * Z[3,1]
  g_tuss[i,] = Z[1,2] + scaled_ground * Z[3,2]
  g_tree[i,] = Z[1,3] + scaled_ground * Z[3,3]
  g_type[i,] = Z[1,4] + scaled_ground * Z[3,4]
  g_a[i,] = Z[1,5] + scaled_ground * Z[3,5]
  g_b[i,] = Z[1,6] + scaled_ground * Z[3,6]
}

gris = HPDinterval(as.mcmc(g_grass))
grits = HPDinterval(as.mcmc(g_tuss))
gritr = HPDinterval(as.mcmc(g_tree))
grity= HPDinterval(as.mcmc(g_type))
gria= HPDinterval(as.mcmc(g_a))
grib = HPDinterval(as.mcmc(g_b))

ground = scaled_ground * sd(uruguay_traits$Ground_stratum) + 
  mean(uruguay_traits$Ground_stratum)

  
plot_list = vector("list", 4)

dt = data.frame(x = ground, 
                y = colMeans(g_grass), 
                ylow = gris[,1], 
                yup = gris[,2], 
                variable = 1)

pl = ggplot(dt, aes(x = x, y = y)) +
  geom_ribbon(aes(ymin = ylow, ymax = yup, 
                  fill = as.factor(variable)), 
              alpha = 0.5, fill = "grey70") +
  geom_line(color = "black", size = 2, aes(alpha=0.5)) +
  geom_hline(yintercept = 0) +
  xlab("Ground use") +
  ylab("Grass height coefficient") +
  theme_minimal() +
  ylim(-1.5,1) +
  
  geom_point(data = df, aes(x = ground, 
                            y = b_height,
                            color = 1 - f_height),
                            size = 2) +
  scale_colour_viridis_c(option = "D") +
  theme(legend.position = "none") 

plot_list[[1]] = pl

# tussock --------------------------------
dt = data.frame(x = ground, 
                y = colMeans(g_tuss), 
                ylow = grits[,1], 
                yup = grits[,2], 
                variable = 1)

pl = ggplot(dt, aes(x = x, y = y)) +
  geom_ribbon(aes(ymin = ylow, ymax = yup, 
                  fill = as.factor(variable)), 
              alpha = 0.5, fill = "grey70") +
  geom_line(color = "black", size = 2, aes(alpha=0.5)) +
  geom_hline(yintercept = 0) +
  xlab("Ground use") +
  ylab("Tussock coefficient") +
  theme_minimal() +
  ylim(-1.5,3.3) +
  geom_point(data = df, aes(x = ground, 
                            y = b_tussock, 
                            alpha = 0.5), size = 2) +
  theme(legend.position = "none") 

plot_list[[2]] = pl

# tree --------------------------------------------------------
dt = data.frame(x = ground, 
                y = colMeans(g_tree), 
                ylow = gritr[,1], 
                yup = gritr[,2], 
                variable = 1)

pl = ggplot(dt, aes(x = x, y = y)) +
  geom_ribbon(aes(ymin = ylow, ymax = yup, 
                  fill = as.factor(variable)), 
              alpha = 0.5, fill = "grey70") +
  geom_line(color = "black", size = 2, aes(alpha=0.5)) +
  geom_hline(yintercept = 0) +
  xlab("Ground use") +
  ylab("Tree coefficient") +
  theme_minimal() +
  ylim(-1.5,3.3) +
  geom_point(data = df, aes(x = ground, 
                            y = b_tree, 
                            alpha = 0.5), size = 2) +
  theme(legend.position = "none") 

plot_list[[3]] = pl

# pasture type ----------------------------------------
dt = data.frame(x = ground, 
                y = colMeans(g_type), 
                ylow = grity[,1], 
                yup = grity[,2], 
                variable = 1)

pl = ggplot(dt, aes(x = x, y = y)) +
  geom_ribbon(aes(ymin = ylow, ymax = yup, 
                  fill = as.factor(variable)), 
              alpha = 0.5, fill = "grey70") +
  geom_line(color = "black", size = 2, aes(alpha=0.5)) +
  geom_hline(yintercept = 0) +
  xlab("Ground use") +
  ylab("Pasture type coefficient") +
  theme_minimal() +
  ylim(-1.5,3.3) +
  geom_point(data = df, aes(x = (ground), 
                            y = b_art_pasture, 
                            alpha = 0.5), size = 2) +
  theme(legend.position = "none") 

plot_list[[4]] = pl

ground_plots = plot_list

ggarrange(plot_list[[1]], 
          plot_list[[2]], 
          plot_list[[3]],
          plot_list[[4]],
          labels = c("A", "B", "C", "D"),  
          ncol = 2, nrow = 2, common.legend = F, 
          font.label = list(size = 12, color = "black", face = "bold"),
          vjust = 1)

# ggsave("ground.png", width = 3, height = 3, dpi = 300, scale = 2)

#--------------
#--------------

pl_list = vector("list", 2)
# a ----------------------
dt = data.frame(x = ground, 
                y = (colMeans(g_a)), 
                ylow = (gria[,1]), 
                yup = (gria[,2]), 
                variable = 1)

pl = ggplot(dt, aes(x = x, y = y)) +
  geom_ribbon(aes(ymin = ylow, ymax = yup, 
                  fill = as.factor(variable)), 
              alpha = 0.5, fill = "grey70") +
  geom_hline(yintercept = 0) +
  geom_line(color = "black", size = 2, aes(alpha=0.5)) +
  xlab("Ground use") +
  ylab("Valley") +
  theme_minimal() +
#  ylim(0,1) +
  geom_point(data = df, aes(x = (ground), 
                            y = (a_hill), 
                            alpha = 0.5), size = 2) +
  theme(legend.position = "none") 

pl_list[[1]] = pl


# b ----------------------
dt = data.frame(x = ground, 
                y = (colMeans(g_b)), 
                ylow = (grib[,1]), 
                yup = (grib[,2]), 
                variable = 1)

pl = ggplot(dt, aes(x = x, y = y)) +
  geom_ribbon(aes(ymin = ylow, ymax = yup, 
                  fill = as.factor(variable)), 
              alpha = 0.5, fill = "grey70") +
  geom_line(color = "black", size = 2, aes(alpha=0.5)) +
  geom_hline(yintercept = 0) +
  xlab("Ground use") +
  ylab("Sierras") +
  theme_minimal() +
  ylim(-2,1) +
  geom_point(data = df, aes(x = (ground), 
                            y = (b_hill), 
                            alpha = 0.5), size = 2) +
  theme(legend.position = "none") 

pl_list[[2]] = pl

ggarrange(pl_list[[1]], 
          pl_list[[2]], 
          labels = c("A", "B"),  
          ncol = 2, nrow = 1, common.legend = F, 
          font.label = list(size = 12, color = "black", face = "bold"),
          vjust = 1)

#ggsave("ground_sierras.png", width = 3, height = 1.7, dpi = 300, scale = 2)

gound_sierras = pl_list
#-------------------
#-------------------
  
  pg = ggplot(data = df, aes(x = greg, 
                            y = b_art_pasture) ) +
  
  geom_violin(trim=FALSE) +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=1, alpha = 0.5,
               aes(color = 1 - f_art_pasture)) +
 # scale_colour_viridis_c(option = "D") +
  geom_hline(yintercept = 0) +
    xlab("Gregarious") +
    ylab("Pasture type") +
  theme_minimal()

#------


ggarrange(size_plots[[1]], 
          size_plots[[2]], 
          size_plots[[3]], 
          size_plots[[4]], 
          plot_list[[1]], 
          pg,
          #size_sierras[[1]],
          #size_sierras[[2]],
          #pl_list[[2]],
          labels = paste0("(", letters[1:6],")"),
          hjust = -8,
          ncol = 3, nrow = 2, common.legend = F, 
          font.label = list(size = 12, color = "black", face = "bold"),
          vjust = 2)

#--------------
# Trait effects on probability of detection -------------------------------
tr_det = colnames(TT_det)
zps <- fit_summary[grepl("zp", fit_summary$variable),]

post_zs = fit$draws('zp')
z = as.vector(post_zs[,,])
dim(z) = c(dim(post_zs)[1] * dim(post_zs)[2], dim(post_zs)[3])

# calculate f and means
ff = numeric(ncol(z))
mz = numeric(ncol(z))
for(i in 1:ncol(z)){
  tmp = mean(z[,i])
  ff[i] = sum(sign(z[,i]) == sign(tmp))/nrow(z)
  mz[i] = tmp
}

zps$f = ff


df <- data.frame(x =  tr_det,
                 fz = zps$mean,
                 L = zps$q5,
                 U = zps$q95, 
                 f = zps$f)
df2 = df[-1,]

plot_det2 = ggplot(df2, aes(x = x, y = fz)) +
  geom_linerange(aes(ymin =L , ymax = U), color = "black")+
  geom_point(size = 4) +
#  scale_colour_gradient(low = "lightgrey", high = "black",limits = c(0.5,1.1))+
  theme_minimal()+
  ylab("Coefficient for detection probability")+
  xlab("")+
  theme(axis.text.x = element_text(colour = "black", size = 14, angle = 45), 
        axis.text.y = element_text(colour = "black", size = 12))+
  theme(axis.line.x = element_line(color="black"),
        axis.line.y = element_line(color="black"))+
  theme(axis.title=element_text(size=14))+  
  geom_hline(yintercept = 0, color = "grey", size = 1.2, linetype = "dashed")+
  theme(axis.text.x = element_text(margin = margin(t = 10)))

plot_det2

ggsave("Traits_prob_det.pdf", width = 8, height = 8, dpi = 320,
       scale = 1.4, bg = "white", units = "cm")
#ggsave("Traits_prob_det.png", dpi = 300, scale = 1.5)


#ggarrange(plot_det1, plot_det2,
#          labels = c("A", "B"),  
#          ncol = 2, nrow = 1,  
#          font.label = list(size = 12, color = "black", face = "bold"),
#          vjust = 1)

head(df2)
df2 = cbind(df2, zps[-1,9], zps[-1,10])
colnames(df2) = c("Traits", "Mean", "q0025", "q0975", "f", "Neff", "Rhat")
write.csv(df2, "traits_prob_det.csv", row.names = F)


