library(ggplot2)
library(cowplot)
library(cmdstanr)
library(sfsmisc)
library(ggpubr)
library(coda)


# Load data and analyses --------------------------------------------------

load("data_Joaquin_sep22.RData")
load("matrices_Joaquin_sep22.rdata")

#load("resJ.RData")
load("Occ_PPC.RData") 
load("fit_summary.RData")
load("fit.RData")

# fit_summary <- summary(fit,  probs = c(0.025, 0.05, 0.5, 0.95, 0.975))$summary
# draws <- rstan::extract(fit)
# Histogram of species detectability --------------------------------------

df = read.csv("detection_probability_per_spp.csv")

 ggplot(df, aes(x=probs)) + 
  geom_histogram(color="black", fill="grey", bins = 30)+
  theme_classic()+
  ylab("Count")+
  xlab("Probability of detection")+
  theme(axis.text.x = element_text(colour = "black", size = 14), 
        axis.text.y = element_text(colour = "black", size = 14))+
  theme(axis.line.x = element_line(color="black"),
        axis.line.y = element_line(color="black"))+
  theme(axis.title=element_text(size=16))

 ggsave("hist_detec.png",  width = 9, height = 6, dpi = 300)
 

# Posterior predictive check ----------------------------------------------

 Occ_obs = rep(NA, n.s)
 id = which(grepl("Visita", colnames(Y_data)) == TRUE)
 for(s in 1:n.s)
 {
   tmp = Y_data[which(Y_data$Especie == spp[s]), id]
   tmp2 = apply(tmp, 1, function(x) sum(na.omit(x)))
   Occ_obs[s] = length(which(tmp2 >0))/n.su
 }
 
 #Hago la base de datos y el plot
 
 L = HPDinterval(as.mcmc(OCC), prob  = 0.95)[,1]
 U = HPDinterval(as.mcmc(OCC), prob  = 0.95)[,2]
 df = data.frame(spp = as.character(spp), 
                 tags = as.character(tag), 
                 obs = Occ_obs, L = L, U = U)

 plots = vector("list", 2)
 
 plots[[1]]= ggplot(df, aes(x = tag, y = obs))+
   theme_classic() +
   geom_errorbar(aes(ymin=L, ymax=U), width=.05, size = 1, col = "darkgrey")+
   geom_point(size = 4, alpha=0.3)+
   ylab("Occupancy")+
   xlab("")+
   theme(axis.text.x = element_text(colour = "black", size = 6, angle = 75), 
         axis.text.y = element_text(colour = "black", size = 12))+
   theme(axis.line.x = element_line(color="black"),
         axis.line.y = element_line(color="black"))+
   theme(axis.title=element_text(size=16))+  
   theme(axis.text.x = element_text(margin = margin(t = 20)))
 
 
 mean_pred = apply(OCC, 2, mean)
 df2 = data.frame(x = Occ_obs, y = mean_pred)
 
plots[[2]]= ggplot(df2, aes(x = x, y =y)) +
   geom_point(size = 4, alpha=0.3) +
   theme_classic()+
   ylab("predicted occupancy")+
   xlab("observed occupancy")+
   theme(axis.text.x = element_text(colour = "black", size = 14), 
         axis.text.y = element_text(colour = "black", size = 14))+
   theme(axis.line.x = element_line(color="black"),
         axis.line.y = element_line(color="black"))+
   theme(axis.title=element_text(size=16))+
   geom_line(aes(x=x, y=x), color = "darkgray", size = 1) #, linetype = "dashed") 

plot1 =plots[[1]]+geom_text(x=5, y=0.57, label="A", size = 6)
 plot2 = plots[[2]]+geom_text(x=0.05, y=0.57, label="B", size = 6)+
   coord_fixed(1)

 plot_grid(plot1, plot2)#, rel_widths = c(0.6, 0.4))
tmp = plot_grid(plot1, plot2) #, rel_widths = c(0.6, 0.4))
tmp
ggsave("PPC.png", width = 9, height = 6, dpi = 300, scale = 1)


# Pannel of covariates ----------------------------------------------------
# Sort species according to their occupancy (commonest to rarest) -------------------------------

year = sort(unique(Y_data$A.o))
ny = length(year)
Occ_obs = array(NA, c(n.s, ny))
id = which(grepl("Visita", colnames(Y_data)) == TRUE)
for(y in 1:ny)
{
   xxx = Y_data[which(Y_data$A.o == year[y]),]
   nsites = length(unique(xxx$TransectaID))
   for(s in 1:n.s)
   {
      tmp = Y_data[which(Y_data$Especie == spp[s] & Y_data$A.o == year[y]), id]
      if(nrow(tmp) > 0)
      {
         tmp2 = apply(tmp, 1, function(x) sum(na.omit(x)))
         Occ_obs[s, y] = length(which(tmp2 >0))/nsites
      }else{
         Occ_obs[s, y] = 0  
      }
   }
}
Occ_obs2 = rowMeans(Occ_obs)
Occ_sorted = sort(Occ_obs2, decreasing = T)
length(Occ_sorted) == n.s

#create a vector identifying spp

id_sorted = c()
names_sorted = c()
names_sorted2 = c()
for(s in 1:n.s)
{
   tmp = which(Occ_obs2 == Occ_sorted[s])
   id_sorted = c(id_sorted, tmp)
   names_sorted = c(names_sorted, tag[tmp])
   names_sorted2 = c(names_sorted2, spp[tmp])
}
id_sorted = unique(id_sorted)
names_sorted = unique(names_sorted)
names_sorted2 = unique(names_sorted2)

bs <- fit_summary[grepl("betas", fit_summary$variable),]
post_bs = fit$draws('betas')

##Los betas vienen ordenados as? b1sp1...b1spj,,, bpsp1...bpspj
np = ncol(X)+ ncol(XE)
pnames = c("Height", "brush", "tree", "land_use", "Intercept", "Landscape")
np = length(pnames)
ylabs = c("Grass height (cm)", "Tussock (prop)", "Tree cover (prop)", 
          "Artificial pasture", "Intercept", "Hill")
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
      #tmp1 = density(tmp[,,i], from = min(tmp[,,i])*1.1, to = max(tmp[,,i]*1.1))
      tmp1 = as.vector(tmp[,,i])
      f[i] = length(which(tmp1 < 0 )) / length(tmp1)
      if(mean(tmp1) > 0 ) f[i] = length(which(tmp1 > 0 )) / length(tmp1)
   }
   #reordeno en funci?n de su abundancia
   tmp_mean = tmp_b$mean[id_sorted]
   tmp_L = tmp_b$q5[id_sorted]
   tmp_U = tmp_b$q95[id_sorted]
   tmp_f = f[id_sorted]
   names2 = tolower(names_sorted)
   df = data.frame(x= 1:n.s, y = tmp_mean, L = tmp_L, U = tmp_U, f= tmp_f)
   x =  ggplot(df, aes(x = x, y = y, color=f)) +
      geom_linerange(aes(ymin =L , ymax = U), color = "darkgrey")+
      scale_x_continuous(breaks= seq(1,n.s, by = 3))+
     ylim(-2.7,2.55)+
   
      #geom_point(aes(shape=vuln, color=f), size = 4)+
      geom_point(size = 2)+
      theme(axis.text.x = element_text(colour = "black", size = 9, angle =90 ), 
            axis.text.y = element_text(colour = "black", size = 12))+
      scale_colour_gradient(low = "lightgrey", high = "black",limits = c(0.5,1))+
      theme_classic()+
      ylab(as.character(ylabs[[p]]))+
      xlab("")+
      theme(axis.line.x = element_line(color="black"),
            axis.line.y = element_line(color="black"))+
      theme(axis.title=element_text(size=12))+  
      geom_hline(yintercept = 0, color = "black", size = 0.8, linetype = "dashed")+
      theme(axis.text.x = element_text(margin = margin(t = 0)))+
      coord_flip()+
    #  theme(legend.position = "none") +
      theme(panel.grid.minor.y = element_blank())+
      grids(axis= "y",linetype = "dotted", color = "darkgrey", size = 0.1) 
   plot_list[[p]] =  x
}
plot_list[[4]]


plot1 = plot_list[[1]]+ theme(legend.position='none')+
      theme(plot.margin = unit(c(2,0,2,0), "points"))+theme(aspect.ratio = 0.8)
plot2 = plot_list[[2]]+ theme(legend.position='none')+
       theme(plot.margin = unit(c(2,0,2,0), "points"))+theme(aspect.ratio = 0.8)
plot3 = plot_list[[3]]+ theme(legend.position='none')+
        theme(plot.margin = unit(c(2,0,2,0), "points"))+theme(aspect.ratio = 0.8)
plot4 = plot_list[[4]]+ theme(legend.position='none')+
       theme(plot.margin = unit(c(2,0,2,0), "points"))+theme(aspect.ratio = 0.8)

tmp = ggdraw() +
   draw_plot(plot3, x = 0, y = 0, width = .5, height = .5) +
   draw_plot(plot4, x = .5, y = 0, width = .5, height = .5) +
   draw_plot(plot1, x = 0, y = 0.5, width = 0.5, height = 0.5) +
   draw_plot(plot2, x = 0.5, y = 0.5, width = 0.5, height = 0.5) +
   draw_plot_label(label = c("A", "B", "C", "D"), size = 13,
                   x = c(0.04, 0.54, 0.04, 0.54), y = c(1, 1, 0.5, 0.5))

plot_grid(plot_list[[1]],plot_list[[2]],plot_list[[3]],plot_list[[4]],
          labels = c("A", "B", "C", "D"),
          scale = 0.9,
          label_x = 0.8,
          label_y = 0.9,
          align = 'v')


tmp
ggsave("Covariates_v3.png", dpi = 300, scale = 1)

ggarrange(plot_list[[1]],plot_list[[2]],plot_list[[3]],plot_list[[4]],
          labels = c("A", "B", "C", "D"),  
          ncol = 4, nrow = 1, common.legend = T, legend = "none",
          font.label = list(size = 12, color = "black", face = "bold"),
          vjust = 1.3, hjust = -16)

ggsave("Covariates_local.png", height = 6, width = 9, dpi = 300, scale = 1)

plot_list[[6]]

# Create table of species responses ---------------------------------------

df = read.csv("betas.csv")
tmp_p = unique(df$param)
tmp_p = c("Height", "brush", "tree", "land_use", "Landscape") 
OUT = NULL
for(s in 1:n.s)
{
   id = which(names_sorted2 == spp[s])
   tmp = df[which(df$spp == spp[s]),]
   out = rep(NA, length(tmp_p))
   for(p in 1:length(tmp_p))
   {
      tmp2 = tmp[which(tmp$param == tmp_p[p]),]
      xxx = paste(round(tmp2$mean,2), " [", round(tmp2$q5,2), ",", round(tmp2$q95,2), "]", sep = "")
      out[p] = xxx
      }
 OUT = rbind(OUT,c(spp[s], tag[s], id, round(Occ_sorted[id],4), out))  
}
colnames(OUT) = c("Spp", "Tag", "Nr", "Occ", "Height", "Tussock", "Tree", "Artif.", "Highl.")
write.csv(OUT, "TableB2_betas.csv")



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
                ground = scale(uruguay_traits$Ground_stratum), 
                size= log(uruguay_traits$body_size), 
                greg = uruguay_traits$Greg)

#Pasture height
#Sampleo las iteraciones con las que voy a hacer las curvas de respuesta
post_zs = fit$draws('z')


niter = nrow(post_zs) * ncol(post_zs)
# N = nrow(post_zs)
# id = sample(1:N, niter, replace = T)#iteraciones que he sampleado
# tmp_z = post_zs[id,]

#Matriz de los traits raw
size_raw = log(uruguay_traits$body_size)
ground_raw = uruguay_traits$Ground_stratum
insect_raw = scale(uruguay_traits$Insect)
TT_raw = cbind(1, size_raw, ground_raw, insect_raw, TT_pres[,5])

#Pasture height ~body  size############################

nt = ncol(TT_pres)
p = 1
init = (p-1)*nt+1
fin = init + (nt-1)
tmp_z2 = as.vector(post_zs[,,init:fin])
dim(tmp_z2) = c(dim(post_zs)[1] * dim(post_zs)[2], nt)

id_col = which(colnames(TT_pres) == "Size")
tmp_x = seq(from = min(TT_raw[,id_col]), to = max(TT_raw[,id_col]), length = nrow(TT_pres))
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
plots[[1]] = ggplot(tmp_df, 
                    aes(x, y, mean = mean_, L = L_, U = U_, x2 = x2)) +
   geom_point(size = 3, alpha=0.5, color = "darkgrey") +
   theme_classic()+
  ylim(-1.5,2)+
   ylab("Pasture height")+
   xlab("Body size (log)")+
   theme(axis.text.x = element_text(colour = "black", size = 14), 
         axis.text.y = element_text(colour = "black", size = 14))+
   theme(axis.line.x = element_line(color="black"),
         axis.line.y = element_line(color="black"))+
   theme(axis.title=element_text(size=16))+
   geom_line(aes(x=x2, y=mean), color = "darkgray", size = 1.5) +
   geom_line(aes(x=x2, y=L), color = "darkgray", size = 1, linetype = "solid") +
   geom_line(aes(x=x2, y=U), color = "darkgray", size = 1, linetype = "solid")


# Pasture height ground use ------------------------------------------------

nt = ncol(TT_pres)
p = 1
init = (p-1)*nt+1
fin = init + (nt-1)
tmp_z2 = as.vector(post_zs[,,init:fin])
dim(tmp_z2) = c(dim(post_zs)[1] * dim(post_zs)[2], nt)

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
   geom_point(size = 3, alpha=0.5, color = "darkgrey") +
   theme_classic()+
  ylim(-1.5,2)+
   ylab("Pasture height")+
   xlab("Ground use (%)")+
   theme(axis.text.x = element_text(colour = "black", size = 14), 
         axis.text.y = element_text(colour = "black", size = 14))+
   theme(axis.line.x = element_line(color="black"),
         axis.line.y = element_line(color="black"))+
   theme(axis.title=element_text(size=16))+
   geom_line(aes(x=x2, y=mean), color = "darkgray", size = 1.5) +
   geom_line(aes(x=x2, y=L), color = "darkgray", size = 1, linetype = "dashed") +
   geom_line(aes(x=x2, y=U), color = "darkgray", size = 1, linetype = "dashed")



# Brush cover -------------------------------------------------------------

nt = ncol(TT_pres)
p = 2
init = (p-1)*nt+1
fin = init + (nt-1)
tmp_z2 = as.vector(post_zs[,,init:fin])
dim(tmp_z2) = c(dim(post_zs)[1] * dim(post_zs)[2], nt)

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
   geom_point(size = 3, alpha = 0.5, color = "darkgrey") +
   theme_classic()+
  ylim(-1.5,2)+
   ylab("Brush cover")+
   xlab("Size(log)")+
   theme(axis.text.x = element_text(colour = "black", size = 14), 
         axis.text.y = element_text(colour = "black", size = 14))+
   theme(axis.line.x = element_line(color="black"),
         axis.line.y = element_line(color="black"))+
   theme(axis.title=element_text(size=16))+
   geom_line(aes(x=x2, y=mean), color = "darkgray", size = 1.5) +
   geom_line(aes(x=x2, y=L), color = "darkgray", size = 1, linetype = "dashed") +
   geom_line(aes(x=x2, y=U), color = "darkgray", size = 1, linetype = "dashed")

# Tree cover --------------------------------------------------------------

nt = ncol(TT_pres)
p = 3
init = (p-1)*nt+1
fin = init + (nt-1)
tmp_z2 = as.vector(post_zs[,,init:fin])
dim(tmp_z2) = c(dim(post_zs)[1] * dim(post_zs)[2], nt)


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
   geom_point(size = 3, alpha = 0.5, color = "darkgrey") +
   theme_classic()+
  ylim(-1.5,2)+
   ylab("Tree cover")+
   xlab("Size(log)")+
   theme(axis.text.x = element_text(colour = "black", size = 14), 
         axis.text.y = element_text(colour = "black", size = 14))+
   theme(axis.line.x = element_line(color="black"),
         axis.line.y = element_line(color="black"))+
   theme(axis.title=element_text(size=16))+
   geom_line(aes(x=x2, y=mean), color = "darkgray", size = 1.5) +
   geom_line(aes(x=x2, y=L), color = "darkgray", size = 1, linetype = "dashed") +
   geom_line(aes(x=x2, y=U), color = "darkgray", size = 1, linetype = "dashed")


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
   geom_point(size = 3, alpha = 0.5, color = "darkgrey") +
   theme_classic()+
   ylim(-1.5,2)+
   ylab("Tree cover")+
   xlab("Ground use (%)")+
   theme(axis.text.x = element_text(colour = "black", size = 14), 
         axis.text.y = element_text(colour = "black", size = 14))+
   theme(axis.line.x = element_line(color="black"),
         axis.line.y = element_line(color="black"))+
   theme(axis.title=element_text(size=16))+
   geom_line(aes(x=x2, y=mean), color = "darkgray", size = 1.5) +
   geom_line(aes(x=x2, y=L), color = "darkgray", size = 1, linetype = "dashed") +
   geom_line(aes(x=x2, y=U), color = "darkgray", size = 1, linetype = "dashed")


plots[[6]] = ggplot(df, aes(x=greg, y=b_use)) + 
   geom_boxplot()+
   theme_classic()+
  ylim(-1.5,2)+
   ylab("Antrop.")+
   xlab("Gregarious")+
   theme(axis.text.x = element_text(colour = "black", size = 14), 
         axis.text.y = element_text(colour = "black", size = 14))+
   theme(axis.line.x = element_line(color="black"),
         axis.line.y = element_line(color="black"))+
   theme(axis.title=element_text(size=16))

# combine plots ------
plot1 = plots[[1]]+  ylab("Height (cm)")
plot2 = plots[[2]]+  ylab("Height (cm)")
plot3= plots[[3]]+  ylab("Tussock (prop)")
plot4 = plots[[4]]+  ylab("Tree cover (prop)")
plot5 = plots[[5]]+  ylab("Tree cover (prop)")
plot6= plots[[6]]+ ylab("Artif. pasture")


plot_grid(plot1, plot2, plot3, plot4, plot5, plot6,
          labels = c("A", "B", "C", "D", "E", "F"),
          scale = 0.9,
          label_x = 0.8,
          label_y = 0.9,
          align = 'v')

ggsave("Scatter_trait_effects.png", 
       # device = "pdf",
       # width = 10, height = 10, units = "cm",
       dpi = 300, scale = 1)

# Create table of traits effects ------------------------------------------

df = read.csv("trait_effects_on_betas.csv")
names(df)[1] = "param" 
tmp = which(df$param == "Intercept")
df = df[-tmp,]
tmp_p = unique(df$param)
labs = c("Height (cm)", "Tussock (prop)", "Tree cover (prop)", 
          "Artif. pasture", "Highland")
tmp_t = unique(df$Trait)
lab_tr = c("Body size", "Ground use", "Insectivory", "Gregarism")

OUT  = NULL
for(t in 1:length(tmp_t))
{
   tmp = df[which(df$Trait == tmp_t[t]),]
   out = rep(NA, length(tmp_p))
   for(p in 1:length(tmp_p))
   {
      tmp2 = tmp[which(tmp$param == tmp_p[p]),]
      out[p]= paste(round(tmp2$Mean,2), " [", round(tmp2$q5,2), ",", round(tmp2$q95,2), "]", sep = "")
   }
   OUT = rbind(OUT, c(as.character(lab_tr[t]), out))
}

colnames(OUT) = c("Trait", labs)
write.csv(OUT, "TableB3_trait_effects.csv", row.names = F)
