#########################
# Plotting model output #
#########################

rm(list = ls())

library(ggplot2)
library(ggpubr)

# Plotting settings ####
mytheme <-   theme(panel.grid.minor = element_blank(),
                   panel.grid.major = element_blank(),
                   panel.background = element_rect(colour = "black", fill = "white", 
                                                   size = 1.5),
                   axis.text = element_text(size = 15, colour ="black"),
                   axis.text.x = element_text(margin = unit(c(0.5,0.5,0.5,0.5), "cm")),
                   axis.text.y = element_text(margin = unit(c(0.5,0.5,0.5,0.5), "cm")),                   
                   axis.title = element_text(size = 17, colour ="black"),
                   axis.title.y = element_text(margin = unit(c(-0.2,-0.2,-0.2,-0.2), "cm")),
                   axis.title.x = element_text(margin = unit(c(-0.2,-0.2,-0.2,-0.2), "cm")),
                   axis.ticks.length = unit(-0.25, "cm"),
                   plot.margin = unit(c(1, 1.5, 0.5, 0.5), "lines"),
                   legend.position = c(0.15, 0.75),
                   legend.title = element_blank(),
                   legend.text = element_text(size=15, colour ="black"),
                   legend.key = element_rect(fill = "white"),
                   plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
                   strip.background = element_rect(color = "black", fill = "white"),
                   strip.text = element_text(size = 15, color = "black")
)

#################################
# Relative proteome investments #
#################################

# Temperature dependent data
del <- read.csv("data_output_schaum.csv", header = T)
del$cue <- 1 - (del$resp/del$phot)
del$vol <- 4/3*pi*(del$size/2)^3

n <- nrow(del)
datap <- data.frame(nt = rep(del$temp, 6),
                    panel = c(rep("a", n*2), rep("b", n*2), rep("c", n*2)),
                    proteome = c(rep("photosystems", n), rep("rubisco", n), rep("ribosomes", n),
                                 rep("repair", n), rep("glycolysis", n), rep("lipid degradation", n)),
                    value = c(del$ichl, del$irub, del$irib, del$icha, 
                              del$icbd, del$ilpd)
)
datap$proteome <- factor(datap$proteome, levels = c("photosystems", "rubisco", "ribosomes", 
                                                    "repair", "glycolysis", "lipid degradation"),)

# Relative proteome investments
ggplot(data = datap, aes(x = nt, y = value, color = proteome)) +
  geom_line(size = 1.0) +
  geom_point(size = 3) +
  facet_grid(panel ~ .) +
  xlim(3,40) +
  scale_color_manual(values = c("#018571", "#a6611a", "#ca0020", "#92c5de", "forestgreen", "deeppink3")) +
  scale_y_continuous(breaks = c(0,0.2,0.4), limits = c(0.0,0.4)) +
  ylab("Relative investment") +
  xlab(expression(paste("Temperature  (", degree,"C)"))) +
  mytheme +
  theme(strip.background = element_blank(),
        strip.text = element_blank(),
        legend.position = "none")

# Temperature independent data (validated for the Schaum dataset)
delto <- read.csv("data_output_schaum_tradeoff.csv", header = T)
delto$cue <- 1 - (delto$resp/delto$phot)
delto$vol <- 4/3*pi*(delto$size/2)^3

n <- nrow(delto)
datap <- data.frame(nt = rep(delto$temp, 6),
                    panel = c(rep("a", n*2), rep("b", n*2), rep("c", n*2)),
                    proteome = c(rep("photosystems", n), rep("rubisco", n), rep("ribosomes", n),
                                 rep("repair", n), rep("glycolysis", n), rep("lipid degradation", n)),
                    value = c(delto$ichl, delto$irub, delto$irib, delto$icha, 
                              delto$icbd, delto$ilpd)
)
datap$proteome <- factor(datap$proteome, levels = c("photosystems", "rubisco", "ribosomes", 
                                                    "repair", "glycolysis", "lipid degradation"),)

# Relative proteome investments
pall <- ggplot(data = datap, aes(x = nt, y = value, color = proteome)) +
  geom_line(size = 1.0) +
  geom_point(size = 3) +
  facet_grid(panel ~ .) +
  xlim(3,40) +
  scale_color_manual(values = c("#018571", "#a6611a", "#ca0020", "#92c5de", "forestgreen", "deeppink3")) +
  scale_y_continuous(breaks = c(0,0.2,0.4), limits = c(0.0,0.4)) +
  #scale_x_continuous(breaks = c(0.2, 0.4, 0.6, 0.8, 1.0, 1.2), limits = c(0.1,1.3)) +
  ylab("Relative investment") +
  xlab(expression(paste("Temperature  (", degree,"C)"))) +
  mytheme +
  theme(strip.background = element_blank(),
        strip.text = element_blank(),
        legend.position = "none")

# Sensitivity analyses for the lipid degradation pathway
s <- read.csv("sensi_main.csv", header = T)
sbase <- subset(s, sensi == "a")

sld <- ggplot(data = s, aes(x = temp, y = value, alpha = sensi)) +
  geom_line(size = 1, color = "deeppink3") +
  scale_alpha_manual(values = c(1,0.8,0.5,0.2)) +
  geom_point(data = sbase, aes(x = temp, y = value), col = "deeppink3", size = 3) +
  facet_grid(panel ~ .) +
  xlim(3,40) +
  ylab("Relative investment") +
  xlab(expression(paste("Temperature  (", degree,"C)"))) +
  mytheme +
  theme(legend.position = c(0.1,0.5),
        strip.background = element_blank(),
        strip.text = element_blank())

tiff('fig3.tiff', res = 600, height = 3200, width = 6500)
ggarrange(pall, sld,
          labels = c("a", "b"), font.label = list(size = 22), ncol = 2, nrow = 1)
dev.off()

# Proteome investments for the Liang and O'Donnell datasets

# Temperature independent data (validated for the Liang dataset)
dliang <- read.csv("data_output_liang_tradeoff.csv", header = T)
dliang$cue <- 1 - (dliang$resp/dliang$phot)
dliang$vol <- 4/3*pi*(dliang$size/2)^3

n <- nrow(dliang)
datap <- data.frame(nt = rep(dliang$temp, 6),
                    panel = c(rep("a", n*2), rep("b", n*2), rep("c", n*2)),
                    proteome = c(rep("photosystems", n), rep("rubisco", n), rep("ribosomes", n),
                                 rep("repair", n), rep("glycolysis", n), rep("lipid degradation", n)),
                    value = c(dliang$ichl, dliang$irub, dliang$irib, dliang$icha, 
                              dliang$icbd, dliang$ilpd)
)
datap$proteome <- factor(datap$proteome, levels = c("photosystems", "rubisco", "ribosomes", 
                                                    "repair", "glycolysis", "lipid degradation"),)

# Relative proteome investments
pliang <- ggplot(data = datap, aes(x = nt, y = value, color = proteome)) +
  geom_line(size = 1.0) +
  geom_point(size = 3) +
  facet_grid(panel ~ .) +
  xlim(3,31) +
  scale_color_manual(values = c("#018571", "#a6611a", "#ca0020", "#92c5de", "forestgreen", "deeppink3")) +
  scale_y_continuous(breaks = c(0,0.2,0.4), limits = c(0.0,0.4)) +
  ylab("Relative investment") +
  xlab(expression(paste("Temperature  (", degree,"C)"))) +
  mytheme +
  theme(strip.background = element_blank(),
        strip.text = element_blank(),
        legend.position = "none")

# Temperature independent data (validated for the O'Donnell dataset)
don <- read.csv("data_output_elena_tradeoff.csv", header = T)
don$cue <- 1 - (don$resp/don$phot)
don$vol <- 4/3*pi*(don$size/2)^3

n <- nrow(don)
datap <- data.frame(nt = rep(don$temp, 6),
                    panel = c(rep("a", n*2), rep("b", n*2), rep("c", n*2)),
                    proteome = c(rep("photosystems", n), rep("rubisco", n), rep("ribosomes", n),
                                 rep("repair", n), rep("glycolysis", n), rep("lipid degradation", n)),
                    value = c(don$ichl, don$irub, don$irib, don$icha, 
                              don$icbd, don$ilpd)
)
datap$proteome <- factor(datap$proteome, levels = c("photosystems", "rubisco", "ribosomes", 
                                                    "repair", "glycolysis", "lipid degradation"),)

# Relative proteome investments
pdon <- ggplot(data = datap, aes(x = nt, y = value, color = proteome)) +
  geom_line(size = 1.0) +
  geom_point(size = 3) +
  facet_grid(panel ~ .) +
  xlim(3,31) +
  scale_color_manual(values = c("#018571", "#a6611a", "#ca0020", "#92c5de", "forestgreen", "deeppink3")) +
  scale_y_continuous(breaks = c(0,0.2,0.4), limits = c(0.0,0.4)) +
  ylab("Relative investment") +
  xlab(expression(paste("Temperature  (", degree,"C)"))) +
  mytheme +
  theme(strip.background = element_blank(),
        strip.text = element_blank(),
        legend.position = "none")

tiff('figS4.tiff', res = 600, height = 3200, width = 6500)
ggarrange(pliang, pdon,
          labels = c("a", "b"), font.label = list(size = 22), ncol = 2, nrow = 1)
dev.off()

# Proteome investments when changing Kre



###############################################################
# Plot changes in cell size, transporters and membrane lipids #
###############################################################

# kresp is temperature dependent, M is temperature dependent, Ea_ktr is halved (intermediary)
base_kresp_vol <- del$vol
cmli_T <- del$clip
ctr_half <- del$ctr

# kresp is temperature independent
delto <- read.csv("data_output_schaum_tradeoff.csv", header = T)
delto$vol <- 4/3*pi*(delto$size/2)^3
kresp_vol <- delto$vol

# M is temperature independent
delto <- read.csv("data_output_schaum_mlip.csv", header = T)
delto$vol <- 4/3*pi*(delto$size/2)^3
kmlip_vol <- delto$vol
cmli_noT <- delto$clip

# Mmax is 1.4 times higher than Mmin (intermediary)
delto <- read.csv("output_mlip_intermediary.csv", header = T)
cmli_in <- delto$clip

# ktr is temperature dependent
delto <- read.csv("data_output_schaum_noTtr.csv", header = T)
delto$vol <- 4/3*pi*(delto$size/2)^3
ktr_vol <- delto$vol
ctr_noT <- delto$ctr

# ktr is temperature independent
delT <- read.csv("data_output_schaum_tr.csv", header = T)
delT$vol <- 4/3*pi*(delT$size/2)^3
base_ktr_vol <- delT$vol
ctr_T <- delT$ctr

# Merging datasets
basevol <- data.frame(temp = rep(del$temp, 6),
                      vol = c(base_kresp_vol, kresp_vol, 
                              base_kresp_vol, kmlip_vol,
                              base_ktr_vol, ktr_vol),
                      run = c(rep("baseline", 43),
                              rep("sens", 43),
                              rep("baseline", 43),
                              rep("sens", 43),
                              rep("baseline", 43),
                              rep("sens", 43)),
                      par = c(rep("respiration", 43*2),
                                rep("lipids", 43*2),
                                rep("transporters", 43*2)))

basevol$par  <- factor(basevol$par,  levels = c("transporters", "lipids", "respiration"))
basevol$run  <- factor(basevol$run,  labels = c("T-dependent", "T-independent"))
basevol <- subset(basevol, temp < 39) 

# Plotting changes in cell size

s2 <- ggplot(data = basevol, aes(x = temp, y = log(vol), col = run, linetype = run)) +
  geom_line(size = 1) +
  facet_wrap(~ par) +
  theme(legend.position = "none") +
  scale_linetype_manual(values = c("solid", "dashed")) +
  scale_color_manual(values = c("#de2d26", "#2b8cbe"))+
  mytheme +
  xlab(expression(paste("Temperature  (", degree,"C)"))) +
  ylab("")+
  ylab(expression(paste("Log (", mu,m^3,")", sep=""))) +
  theme(legend.position = c(0.16,0.79),
        legend.text = element_text(size=12))

tiff('fig4.tiff', res = 600, height = 2800/1.5, width = 4500)
ggarrange(s2, font.label = list(size = 22), ncol = 1, nrow = 1)
dev.off()

# PLotting the concentration of transporters and membrane lipids per biovolume

memb <- data.frame(temp = rep(del$temp, 6),
                   value = c(ctr_T, ctr_noT, ctr_half, cmli_T, cmli_noT, cmli_in),
                   par = c(rep("nutrient transporters", 3*length(ctr_T)),rep("membrane lipids", 3*length(ctr_T))),
                   run = c(rep("T-dependent", length(ctr_T)),rep("T-independent", length(ctr_T)), rep("Intermediary", length(ctr_T)),
                           rep("T-dependent", length(ctr_T)),rep("T-independent", length(ctr_T)), rep("Intermediary", length(ctr_T))))
memb$par <- factor(memb$par, levels = c("nutrient transporters", "membrane lipids"))
memb$run <- factor(memb$run, levels = c("T-dependent", "Intermediary", "T-independent"))

s3 <- ggplot(data = memb, aes(x = temp, y = value, col = par, linetype = run, alpha = run)) +
  geom_line(size = 1.3) +
  facet_wrap(~ par, scales = "free") +
  theme(legend.position = "none") +
  scale_linetype_manual(values = c(1,2,5)) +
  scale_alpha_manual(values = c(1,0.7,0.4)) +
  scale_color_manual(values = c("slateblue4", "goldenrod4")) +
  mytheme +
  ylab(expression(paste("molecules ", mu,m^-3, sep=""))) +
  xlab(expression(paste("Temperature  (", degree,"C)"))) +
  theme(legend.position = c(0.26,0.85),
        legend.text = element_text(size=12)) +
  guides(color = "none",
         linetype = guide_legend(override.aes = list(size = 0.7)))

tiff('figS5.tiff', res = 600, height = 2800, width = 5000)
ggarrange(s3, font.label = list(size = 22), ncol = 1, nrow = 1)
dev.off()

###########################################
# Photosystems: functional versus damaged #
###########################################

# Quasi-steady state approximation
del <- read.csv("data_output_schaum.csv", header = T)
# Optimizing for the damages pool, Kre = 10
del2 <- read.csv("data_output_schaum_phot.csv", header = T)
# Optimizing for the damages pool, Kre = 1000
del3 <- read.csv("data_output_schaum_phot3.csv", header = T)

pho <- data.frame(temp = rep(del$temp, 6),
                  value = c(del$cp, del$cup, del2$cp, del2$cup, del3$cp, del3$cup),
                  par = c(rep("functional", length(del$cp)),
                          rep("damaged", length(del$cup)),
                          rep("functional", length(del$cp)),
                          rep("damaged", length(del$cup)),
                          rep("functional", length(del$cp)),
                          rep("damaged", length(del$cup))),
                  run = c(rep("quasi stedy-state", length(del$cp)*2), rep("optimizing damaged pool, Kre = 10", length(del$cp)*2),
                          rep("optimizing damaged pool, Kre = 1000", length(del$cp)*2)))

pho$run <- factor(pho$run, levels = c("quasi stedy-state", "optimizing damaged pool, Kre = 10", "optimizing damaged pool, Kre = 1000"))

phot <- ggplot(data = pho, aes(x = temp, y = value, col = par, linetype = par, size = par)) +
  geom_line(size = 1) +
  facet_wrap(~ run) +
  xlim(3,39) +
  scale_linetype_manual(values = c("solid", "dashed")) +
  scale_color_manual(values = c("#de2d26", "#2b8cbe")) +
  ylab(expression(paste("photosystems ", mu,m^-3, sep=""))) +
  xlab(expression(paste("Temperature  (", degree,"C)"))) +
  mytheme +
  theme(legend.position = c(0.5, 0.33))

tiff('figS2.tiff', res = 600, height = 2200, width = 2400*3)
ggarrange(phot, font.label = list(size = 22), ncol = 1, nrow = 1)
dev.off()

#############
# Ribosomes #
#############

pri <- ggplot(data = del, aes(x = temp, y = cri)) +
  geom_line(size = 1, col = "#ca0020") +
  ylab(expression(paste("ribosomes ", mu,m^-3, sep=""))) +
  xlab(expression(paste("Temperature  (", degree,"C)"))) +
  xlim(3,39) +
  mytheme

tiff('figS3.tiff', res = 600, height = 2200, width = 2800)
ggarrange(pri, font.label = list(size = 22), ncol = 1, nrow = 1)
dev.off()

#################################################
# Extra sensitivity analyses: cother, fdr, Mmax #
#################################################

# Baseline dataset
b <- read.csv("output_baseline.csv", header = T)
base <- b$mu
bases <- b$size
basel <- b$ilpd

# Sensitivity tests
scd <- read.csv("output_cother_double.csv", header = T)
sch <- read.csv("output_cother_halved.csv", header = T)
smd <- read.csv("output_Mmax_doubled.csv", header = T)
smh <- read.csv("output_Mmax_halved.csv", header = T)
sfd <- read.csv("output_fdr_double.csv", header = T)
sfh <- read.csv("output_fdr_halved.csv", header = T)
n <- length(del$temp)

senmu <- data.frame(temp = rep(del$temp, 9),
                    value = c(scd$mu, sch$mu, base, smd$mu, smh$mu, base, sfd$mu, sfh$mu, base),
                    sensi = c(rep("cother", 3*n), rep("Mmax", 3*n), rep("fdr", 3*n)),
                    run = c(rep("doubled", n), rep("halved", n), rep("baseline", n),
                            rep("doubled", n), rep("halved", n), rep("baseline", n),
                            rep("doubled", n), rep("halved", n), rep("baseline", n)),
                    par = rep("growth rate", n*9))

sensize <- data.frame(temp = rep(del$temp, 9),
                      value = c(scd$size, sch$size, bases, smd$size, smh$size, bases, sfd$size, sfh$size, bases),
                      sensi = c(rep("cother", 3*n), rep("Mmax", 3*n), rep("fdr", 3*n)),
                      run = c(rep("doubled", n), rep("halved", n), rep("baseline", n),
                              rep("doubled", n), rep("halved", n), rep("baseline", n),
                              rep("doubled", n), rep("halved", n), rep("baseline", n)),
                      par = rep("cell size", n*9))

senld <- data.frame(temp = rep(del$temp, 9),
                    value = c(scd$ilpd, sch$ilpd, basel, smd$ilpd, smh$ilpd, basel, sfd$ilpd, sfh$ilpd, basel),
                    sensi = c(rep("cother", 3*n), rep("Mmax", 3*n), rep("fdr", 3*n)),
                    run = c(rep("doubled", n), rep("halved", n), rep("baseline", n),
                            rep("doubled", n), rep("halved", n), rep("baseline", n),
                            rep("doubled", n), rep("halved", n), rep("baseline", n)),
                    par = rep("lipid degradation", n*9))

sensi <- rbind(senmu, sensize, senld)

sensi$par <- factor(sensi$par, levels = c("growth rate",
                                          "lipid degradation",
                                          "cell size"))

psensi <- ggplot(data = sensi, aes(x = temp, y = value, linetype = run, col = run)) +
  geom_line(size = 1) +
  scale_color_manual(values = c("black", "#92c5de", "#ca0020")) +
  scale_linetype_manual(values = c("solid", "dashed", "dashed")) +
  facet_grid(rows = vars(par),
             cols = vars(sensi),
             scales = "free") +
  xlab(expression(paste("Temperature  (", degree,"C)"))) +
  #scale_y_continuous(breaks = seq(0, 75, 20)) +
  ylab("") +
  mytheme +
  theme(legend.position = c(0.11, 0.21))

tiff('figS6.tiff', res = 600, height = 4000, width = 5000)
ggarrange(psensi, font.label = list(size = 22), ncol = 1, nrow = 1)  
dev.off()





