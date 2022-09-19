#############################
# Proteome model validation #
# Suzana G Leles            #
# 11 March 2022             #
#############################

rm(list = ls())

library(ggplot2)
library(ggpubr)
library(viridis)
library(scales)

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

# Import model output
data_output <- read.csv("data_output_schaum.csv", header = T)
data_output$cue <- 1 - (data_output$resp/data_output$phot)
data_output$vol <- 4/3*pi*(data_output$size/2)^3

data_output_zoe <- read.csv("data_output_liang.csv", header = T)
data_output_zoe$cue <- 1 - (data_output_zoe$resp/data_output_zoe$phot)
data_output_zoe$vol <- 4/3*pi*(data_output_zoe$size/2)^3

data_output_elena <- read.csv("data_output_elena.csv", header = T)
data_output_elena$cue <- 1 - (data_output_elena$resp/data_output_elena$phot)
data_output_elena$vol <- 4/3*pi*(data_output_elena$size/2)^3

# Empirical data from Schaum et al 2018 Nat Comm
data_growth <- read.csv("data_growth.csv", header = T)
data_growth$temp <- as.double(data_growth$temp)
data_growth_val <- subset(data_growth, lineage == "t22")

# Empirical data from Schaum et al 2018 Nat Comm
data_cue <- read.table("data_schaum_cue.txt", header = T)
data_cue <- subset(data_cue, treatment == "t22")
temp_size <- c(15.0, 20.0, 22.0, 26.0, 30.0, 32.0, 35.0)
size <- c(11.5, NA, 10.6, 10.4, 10.2, 9.7, 9.3)  # cell diameter (um)
vol <- (4/3*pi*(size/2)^3)
CN <- c(7.32, 7.6, 6.94, 7.31, 7.38, 8.24, 7.86) # units of molC/molN
NC <- (1/CN) 
ChlC <- c(0.03, 0.04, 0.08, 0.1, 0.05, 0.02, 0.02) # mg Chl/ mg C
Cvol <- c(1.83, 2.77, 2.87, 3.32, 6.41, 7.42, 9.5) # fgC/um3
chlvol <- Cvol * ChlC
data_schaum2 <- data.frame(temp_size, size, NC, ChlC, Cvol, chlvol, vol)

# Empirical data from O'Donnell et al 2021 L&O
data_elen <- data.frame(temp = c(10,10,16,16,26,26,31,31),
                        chlvol = c(0.0042,0.0052,0.0043,0.0059,0.0085,0.0125,0.007,0.01),
                        NC = c(0.19, 0.22, 0.18, 0.21, 0.16, 0.17, 0.16, 0.17))

# Empirical data from O'Donnell et al 2018
data_elen_mu <- data.frame(temp = c(3.4, 3.4, 9.6, 9.6, 15.8, 15.8, 20.3, 20.3, 25, 25, 26.2, 26.2, 28.7, 31, 31, 32.4),
                           mu = c(0.02, 0.2, 0.33, 0.43, 0.96, 0.65, 1.11, 0.88, 1.2, 0.9, 1.1, 0.9, 0.96, 0.7, 0.43, 0.02))

# Empirical data from Liang et al 2019
liang_mu <- c(0.63, 0.66, 0.78, 0.88, 0.91, 1.02, 1.07, 1.1, 1.46, 1.42, 1.57, 1.64, 1.68, 1.5, 1.4, 1.33, 0)
liang_muT <- c(15.2, 15.2, 17.5, 17.5, 17.5, 20, 20, 20, 22, 22, 25, 25, 25, 28.4, 28.4, 28.4, 30)
data_liang_mu <- data.frame(mu = liang_mu, temp = liang_muT)
vol <- c(1696, 1360, 969, 452, 556, 937)
data_liang <- data.frame(mu = c(0.78, 1.17, 0.56, 1.06, 1.63, 1.41),
                         vol = c(1696, 1360, 969, 452, 556, 937),
                         chlvol = c(1.4, 1.8, 2.6, 2.3, 2.8, 2.9),
                         temp = c(12, 16, 22, 19, 25, 28),
                         source = rep("Liang et al 2019", 6))

# Compiling into a single dataframe

# Growth rates
data_mu <- data.frame(mu2 = c(data_growth_val$mu, data_liang_mu$mu),
                      temp = c(data_growth_val$temp, data_liang_mu$temp),
                      source = c(rep("Schaum et al 2018", nrow(data_growth_val)), rep("Liang et al 2019", nrow(data_liang_mu))))
data_mu$source <- factor(data_mu$source, levels = c("Schaum et al 2018","Liang et al 2019"))

# Cell volume
data_size <- data.frame(vol = c(data_schaum2$vol, data_liang$vol),
                        temp = c(data_schaum2$temp_size, data_liang$temp),
                        source = c(rep("Schaum et al 2018", 7), rep("Liang et al 2019", 6)))
data_size$source <- factor(data_size$source, levels = c("Schaum et al 2018","Liang et al 2019"))

# Chlorophyll per volume
data_chlvol <- data.frame(chlvol =c(data_schaum2$chlvol, data_liang$chlvol),
                          temp = c(data_schaum2$temp_size, data_liang$temp),
                          source = c(rep("Schaum et al 2018", 7), rep("Liang et al 2019", 6)))
data_chlvol$source <- factor(data_chlvol$source, levels = c("Schaum et al 2018","Liang et al 2019"))

# Nitrogen to carbon cellular quota
data_nc <- data.frame(NC = c(data_schaum2$NC, data_elen$NC),
                      temp = c(data_schaum2$temp_size, data_elen$temp),
                      source = c(rep("Schaum et al 2018", 7), rep("O'Donnell et al 2021", 8)))
data_nc$source <- factor(data_nc$source, levels = c("Schaum et al 2018","O'Donnell et al 2021"))

# Plotting

pgrowth <- ggplot(data_output, aes(x = temp, y = mu)) +
  geom_line(size = 1.0) +
  geom_line(data = data_output_zoe, aes(x = temp, y = mu), color = "black", size = 1.0, linetype = "dashed") +
  geom_point(data = data_mu, aes(x = temp, y = mu2, shape = source), fill = "darkgray", size = 3) +
  scale_shape_manual(values = c(21, 24)) +
  scale_fill_viridis(discrete=TRUE) +
  scale_y_continuous(breaks = c(0.0,1.0,2.0), limits = c(0.0,2), labels = label_number(accuracy = 0.1)) +
  ylab(expression(paste("Growth rate (", day^-1,")", sep=""))) +
  xlab(expression(paste("Temperature  (", degree,"C)"))) +
  mytheme +
  theme(legend.position = c(0.25, 0.87))

data_cue$phot2 <- data_cue$phot * 0.000125*10^3 #(fmoles C/um3/day)
data_output$phot2 <- (data_output$phot*60*24*10^15)/(6.02 * 1e23) # (fmoles C/um3/day)

pru <- ggplot(data = data_output, aes(x = temp, y = phot2)) +
  geom_line(size = 1.0) +
  geom_point(data = data_cue, aes(x = temp, y = phot2), fill = "darkgray", size = 3, pch = 21) +
  scale_shape_manual(values = c(25)) +
  ylab(expression(paste("(fmol C", mu,  m^-3, day^-1,")", sep=""))) +
  xlab(expression(paste("Temperature  (", degree,"C)"))) +
  mytheme +
  theme(legend.position = c(0.15, 0.25)) +
  annotate("text", x = 12, y = 2.2, label= "Schaum et al 2018", size = 5.5)

psize <- ggplot(data = data_output, aes(x = temp, y = vol/1000)) +
  geom_line(size = 1.0) +
  geom_line(data = data_output_zoe, aes(x = temp, y = (vol+150)/1000), color = "black", size = 1.0, linetype = "dashed") +
  scale_y_continuous(breaks = c(0.0,2.5,5.0), limits = c(0.0,5), labels = label_number(accuracy = 0.1)) +
  xlim(4,42) +
  geom_point(data = data_size, aes(x = temp, y = vol/1000, shape = source), size = 3, fill = "darkgray") + 
  scale_shape_manual(values = c(21, 24)) +
  ylab(expression(paste("Cell volume (\u00D7 ", 10^3, mu,m^3,")", sep=""))) +
  xlab(expression(paste("Temperature  (", degree,"C)"))) +
  mytheme +
  theme(legend.position = c(0.35, 0.83))

pcue <- ggplot(data = data_output, aes(x = temp, y = cue)) +
  geom_line(size = 1.0) +
  scale_shape_manual(values = c(25)) +
  geom_point(data = data_cue, aes(x = temp, y = cue), fill = "darkgray", size = 3, pch = 21) +
  scale_y_continuous(breaks = c(0.0,0.5,1.0), limits = c(0.0,1)) +
  ylab("Carbon use efficiency") +
  xlab(expression(paste("Temperature  (", degree,"C)"))) +
  mytheme +
  theme(legend.position = c(0.15, 0.25)) +
  annotate("text", x = 12, y = 0.93, label= "Schaum et al 2018", size = 5.5)


pchl <- ggplot(data = data_output, aes(x = temp, y = cp* 140 / 6.02e23 * 893.5 * 10^15)) +
  geom_line(size = 1.0) +
  geom_line(data = data_output_zoe, aes(x = temp, y = cp* 140 / 6.02e23 * 893.5 * 10^15), color = "black", size = 1.0, linetype = "dashed") +
  geom_point(data = data_chlvol, aes(x = temp, y = chlvol, shape = source), fill = "darkgray", size = 3) +
  scale_shape_manual(values = c(21, 24)) +
  scale_y_continuous(breaks = c(0.0,2.0,4.0), limits = c(0.0,4), labels = label_number(accuracy = 0.1)) +
  ylab(expression(paste("Chl-a (fg  ", mu,m^-3,")", sep=""))) +
  xlab(expression(paste("Temperature  (", degree,"C)"))) +
  mytheme +
  theme(legend.position = c(0.25, 0.83))

pnc <- ggplot(data = data_output, aes(x = temp, y = NC2)) +
  geom_line(size = 1.0) +
  geom_line(data = data_output_elena, aes(x = temp, y = NC), color = "black", size = 1.0, linetype = "dotted") +
  ylim(0.05, 0.35) +
  geom_point(data = data_nc, aes(x = temp, y = NC, shape = source), size = 3, fill = "darkgray") +
  scale_shape_manual(values = c(21, 23)) +
  ylab("N:C ratio") +
  xlab(expression(paste("Temperature  (", degree,"C)"))) +
  mytheme +
  theme(legend.position = c(0.25, 0.83))


tiff('fig2.tiff', res = 600, height = 4900, width = 11000)
ggarrange(pgrowth, pru, pcue, pchl, pnc, psize, 
          labels = c("a", "b", "c", "d", "e", "f"), font.label = list(size = 22), ncol = 3, nrow = 2)
dev.off()

###########################
# EXTRA VALIDATION PLOTS #
##########################

# Thermal growth curves

model <- read.csv("model_val_growth.csv", header = T)

model$treatment <- factor(model$treatment, 
                                labels = c("Liang14", "Liang24", "ODonnell16", "ODonnell31", 
                                           "SchaumAnc", "Schaum22", "Schaum26", "Schaum32", "SchaumFS"))

data_growth$treatment <- factor(data_growth$treatment, 
                                labels = c("Liang14", "Liang24", "ODonnell16", "ODonnell31", 
                                           "SchaumAnc", "Schaum22", "Schaum26", "Schaum32", "SchaumFS"))

p1 <- ggplot(data_growth, aes(x = temp, y = mu, shape = ref)) +
  geom_point(size = 3.0, fill = "darkgray") +
  scale_shape_manual(values = c(24, 23, 21)) +
  scale_linetype_manual(values = c("dashed", "dotted", "solid")) +
  geom_line(data = model, aes(x = temp, y = mu, linetype = ref), size = 1.0) +
  ylim(0, 2.5) +
  ylab(expression(paste("Growth rate (", day^-1,")", sep=""))) +
  xlab(expression(paste("Temperature  (", degree,"C)"))) +
  mytheme +
  theme(legend.position = c(.45,.2)) +
  facet_wrap(~ treatment)

tiff('figS1.png', res = 600, height = 4900, width = 4900)
ggarrange(p1)
dev.off()
