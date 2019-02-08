### Thatch effects ###
rm(list = ls())
## paper idea: lasting effects of thatch... the ghost of competition past is more important than present competition for native annual forbs

#### Load Libraries ####

library(plyr)
library(dplyr)
library(sjstats)
library(Rmisc)
library(ggplot2)
library(lme4)
library(lmerTest)
library(reshape2)

#### Load R data files ####

# Bootstrapped parameters
load("Data/Analysis-Output/Bootstrapped-params/BS-F.Rdata")
load("Data/Analysis-Output/Bootstrapped-params/BS-g.Rdata")
load("Data/Analysis-Output/Bootstrapped-params/BS-m.Rdata")

# Simulated Lambda
load("Data/Analysis-Output/Bootstrapped-params/BS-L.Rdata")

#### Process Data ####

###
# Grass
###
grass <- read.csv("Data/Post-Processing/grass-cover.csv")
names(grass)[3] <- "cover.g"
grass <- filter(grass, Year == 2017, Treat.Code == "C")
treat <- read.csv("Data/Marina-Treatment-30.csv")

###
# Germ and mortality data
###

dem <- read.csv("Data/Post-Processing/dem-data-17.csv")
dem <- merge(dem, grass, by = c("Year", "Plot", "Subplot", "Treat.Code"), all.x = T)
dem <- dem[, -c(8,9,12)] # get rid of extra columns
dem <- filter(dem, Year == 2017, Treat.Code == "C")
dem$Subplot <- factor(dem$Subplot, levels = c("No Grass", "Grass", "Thatch"))
dem$p.mort <- dem$tot.mort/dem$germ.proj
dem$p.germ <- dem$germ.tot/dem$viable

###
# Flowering/seed set data 
###

flo.seed <- read.csv("Data/Post-Processing/final-flo-seed.csv")
flo.seed <- filter(flo.seed, Year == 2017, Treat.Code == "C")
flo.seed$Subplot <- factor(flo.seed$Subplot, levels = c("No Grass", "Grass", "Thatch"))
flo.seed$n.seed.ind <- flo.seed$n.seed.ind * flo.seed$p.viable # adjust for viability
flo.seed <- flo.seed[,-16]
flo.seed$Subplot <- factor(flo.seed$Subplot, levels = c("No Grass", "Grass", "Thatch"))

###
# Seed survival
###
sb <- read.csv("Data/Post-Processing/seed-carryover-plot.csv")[,c(1:9)]
sb <- filter(sb, Year == 2017, Treat.Code == "C")
names(sb)[1] <- "Species"
sb <- sb[,-c(5,6)]

sb.spp <- ddply(sb, .(Species), summarize, p.sd.surv = mean(p.surv))

#### Prep: Traits ####
trait.w <- read.csv("Data/Post-Processing/final-traits-w.csv")
trait.w <- trait.w[,-c(5,7,9)]

sla.13 <- read.csv("Data/Post-Processing/final-sla-13c.csv")[,c(1,3,5)]

PCA.G <- prcomp(trait.w[,c(18, 20)], scale = T) # Greenhouse trait
PCA.F <- prcomp(trait.w[, c(7, 8)], scale = T) # Field trait
PCA.s13 <- prcomp(sla.13[, c(2, 3)], scale = T) # Full SLA + D13C field


trait.w$PC.G <- PCA.G$x[,1] 
trait.w$PC.F <- PCA.F$x[,1] 
trait.w$PC.s13 <- PCA.s13$x[c(2,3,5,6,7,9),1] 

# merge datasets with traits
dem <- merge(dem, trait.w[,c(1,23:25)], by = "Species")
sb <- merge(sb, trait.w[,c(1,23:25)], by = "Species")
flo.seed <- merge(flo.seed, trait.w[,c(1,23:25)], by = "Species")

rm(PCA.F, PCA.G, PCA.s13, sla.13, grass)

# merge datasets with grass cover
#dem <- merge(dem, grass, all.x = T)

#### Prep: Plot-level Lambda ####

full <- merge(dem, flo.seed[,-c(6:13)], by = c("Year","Plot","Treat.Code", "Subplot","Species", "PC.F", "PC.G", "PC.s13"), all = T)

full <- merge(full, sb, by = c("Year","Plot","Treat.Code","Species","PC.F", "PC.G", "PC.s13"), all.x = T)

full$p.germ <- full$germ.tot/full$viable
full$L.sb <- full$p.surv*(1-full$p.germ)
full$L.sa <- ifelse(full$germ.proj > 0 , ((full$germ.proj - full$tot.mort)/full$germ.proj), 0)
full$L.seeds <- ifelse(full$n.seed.ind >= 0, full$L.sa * full$n.seed.ind, 0)

full$L <- ifelse(full$L.sa == 0,  full$L.sb, full$L.sa*full$L.seeds) #5 missing 


#### M1: Germination ####
m1.t <- glmer(cbind(germ.tot, viable-germ.tot) ~ Subplot * PC.F + (1 + Subplot|Plot:Species), family = binomial, data = dem[dem$Subplot != "Grass",])
plot(fitted(m1.t), resid(m1.t))
overdisp((m1.t))
summary(m1.t) 

# germination is lower in thatch across species, with avoiders more strongly affected than tolerators

#### M2: Mortality ####
m2.t <- glmer(cbind(tot.mort, germ.proj-tot.mort) ~ Subplot * PC.F + (1|Plot:Species), family = binomial, data = dem, glmerControl(calc.derivs = F))
plot(fitted(m2.t), resid(m2.t))
overdisp(m2.t)
summary(m2.t) # thatch increases mortality, still higher mortality in Tolerators; no interaction between current grass/drought tolerance or thatch/drought tolerance

ranef(m2.t)
#### M3: Seed set ####

# random effect structure too complex with slope for sublot; missing 20 levels because of low germ/high mort, with the nestedness structure of the previous models there are too many levels of random effects for estimates
m3.t <- lmer(log(n.seed.ind + 1) ~ Subplot * PC.F + (1|Plot:Species), flo.seed)
plot(fitted(m3.t), resid(m3.t))
qqnorm(resid(m3.t))
qqline(resid(m3.t), col = 2, lwd = 2, lty = 2)
summary(m3.t)


# Though grass lowers seed set, the majority of effects of grass occurs through thatch (by lowering light levels?); these species are adapted to high light/lower moisture environments; thatch lowers germination in all species, increases mortality in all species and lowers seed set in all species but moreso in drought avoiders than drought tolerators; I can see the mechanism for germination (light levels?) but how is thatch increasing mortality? light still? it cant be watering in this type of year. 

#### M4: Plot Lambda ####
hist(log(full$L + .3))

m4 <- lmer(log(L + .3) ~ Subplot * PC.F + (1|Plot:Species), full)
plot(fitted(m4), resid(m4))
qqnorm(resid(m4))
qqline(resid(m4), col = 2, lwd = 2, lty = 2)
summary(m4)

#### M5: Grass Cover ####
hist(sqrt(grass$cover.g))
m5 <- lm(sqrt(cover.g) ~ Subplot, grass)
plot(fitted(m5), resid(m5))
qqnorm(resid(m5))
qqline(resid(m5), col = 2, lwd = 2, lty = 2)
summary(m5)

#### M6: Grass seed set ####



#### Bootstraps for simulated Lambda ####
# Dataframe for parameter boots
df.bt <- expand.grid(Subplot = c("Grass", "No Grass", "Thatch"), PC.F = unique(trait.w$PC.F)) 


###
# Germination
###
# 
# BS.g <- bootMer(m1.t, function(x) predict(x, newdata = df.bt[df.bt$Subplot != "Grass",], type = 'response', allow.new.levels = T, re.form = NA), nsim = 1000, ncpus = 2, verbose = T) #26/1000 didnt converge
# 
# save(BS.g, file = "Data/Analysis-Output/Bootstrapped-params/BS-g.Rdata")

###
# Mortality
###
# BS.m <- bootMer(m2.t, function(x) predict(x, newdata = df.bt, type = 'response', allow.new.levels = T, re.form = NA), nsim = 1000, ncpus = 2, verbose = T)
# 
# save(BS.m, file = "Data/Analysis-Output/Bootstrapped-params/BS-m.Rdata")

##
# Seed Set
##

# BS.F <- bootMer(m3.t, function(x) predict(x, newdata = df.bt, type = 'response', allow.new.levels = T, re.form = NA), nsim = 1000, ncpus = 2, verbose = T) 
# 
# save(BS.F, file = "Data/Analysis-Output/Bootstrapped-params/BS-F.Rdata")

#### Prep: Boots seed avg ####

BS <- list(BS.m, BS.F, BS.g)
# Create dataframe of parameter means and sds
params <- c("p.mort.b", "seed.set.b", "p.mort.sd",  "seed.set.sd")
pm <- as.data.frame(matrix(NA, 18, length(params)))
names(pm) = params
pm <- cbind(df.bt, pm)

for(i in 1:2){
    pm[, i+2] = apply(BS[[i]]$t, 2, function(x) mean(x))
    pm[, i+4] = apply(BS[[i]]$t, 2, function(x) sd(x))
}

# grass germ = no grass germ
pm.gm <- as.data.frame(apply(BS[[3]]$t, 2, function(x) mean(x)))
pm.gsd <- as.data.frame(apply(BS[[3]]$t, 2, function(x) sd(x)))
pm.g <- cbind(pm.gm, pm.gsd)
names(pm.g) <- c("p.germ.b", "p.germ.sd")
pm.g <- merge(pm.g, df.bt, by = "row.names")
pm <- merge(pm, pm.g, by = c("Subplot", "PC.F"), all = T)
pm[pm$Subplot == "Grass",]$p.germ.b <- pm[pm$Subplot == "No Grass",]$p.germ.b
pm[pm$Subplot == "Grass",]$p.germ.sd <- pm[pm$Subplot == "No Grass",]$p.germ.sd
pm <- pm[,-7]

# merge with seed survival averages
pm <- merge(pm, trait.w[,c(1,24)], by = "PC.F")
pm <- merge(pm, sb.spp, by = "Species")

# Revalue levels for creating scenarios and matching later
pm$Species <- revalue(pm$Species, c("Agoseris heterophylla" =  "AGHE", "Calycadenia pauciflora" = "CAPA", "Clarkia purpurea" = "CLPU", "Hemizonia congesta" = "HECO", "Lasthenia californica" = "LACA", "Plantago erecta" = "PLER"))
pm$Subplot <- revalue(pm$Subplot, c("No Grass" = "N", "Grass" = "G", "Thatch" = "T"))
pm$Scenario <- paste(pm$Subplot, pm$Species, sep = ".")
Scenario <- pm$Scenario
later <- pm[,c(1,2,3,11)]

# Create empty Lambda dataframe to fill
L.sim <- as.data.frame(Scenario)
sims = 10001
L.sim[,c(2:sims)] <- NA

# Other dfs for storage and later merger
g.sim <- L.sim
m.sim <- L.sim
s.sim <- L.sim
F.sim <- L.sim

# Calculate L from random draws of parameter distributions
for(j in 2:sims){
  for(i in 1:length(Scenario)){
    g.sim[i,j] <- rbinom(n = 1, size = 100, prob = pm[i, "p.germ.b"])/100 #germ
    s.sim[i,j] <- rbinom(n = 1, size = 100, prob = pm[i, "p.sd.surv"])/100 #seed surv
    #m.sim[i,j] <- rbinom(n = 1, size = 100, prob = pm[i, "p.mort.b"])/100
    m.sim[i,j] <- rnorm(n = 1, mean = pm[i, "p.mort.b"], sd = pm[i, "p.mort.sd"])
    F.sim[i,j] <- exp(rnorm(n = 1, mean = pm[i, "seed.set.b"], sd = pm[i, "seed.set.sd"])) - 1 #fecund
    L.sim[i,j] = s.sim[i,j]*(1 - g.sim[i,j]) + g.sim[i,j]*(1 - m.sim[i,j])*F.sim[i,j]
  }
}

# melt all dataframes
sim.list <- list(L.sim, g.sim, s.sim, m.sim, F.sim)
for(i in 1:5){ 
  sim.list[[i]] <- melt(sim.list[[i]])
}

# merge all dataframes
names(sim.list[[1]])[3] <- "L"
names(sim.list[[2]])[3] <- "g"
names(sim.list[[3]])[3] <- "s"
names(sim.list[[4]])[3] <- "m"
names(sim.list[[5]])[3] <- "Fe"

L.sim <- merge(sim.list[[1]], sim.list[[2]], by = c("Scenario", "variable"))
L.sim <- merge(L.sim, sim.list[[3]], by = c("Scenario", "variable"))
L.sim <- merge(L.sim, sim.list[[4]], by = c("Scenario", "variable"))
L.sim <- merge(L.sim, sim.list[[5]], by = c("Scenario", "variable"))
L.sim <- merge(L.sim, later, by = "Scenario")

B.L <- L.sim[,c(8:10,3:7)]

B.L <- ddply(B.L, .(Species, Subplot, PC.F), summarize, L = mean(L), g = mean(g), sb = mean(s), m = mean(m), Fe = mean(Fe))
save(B.L, file = "Data/Analysis-Output/Bootstrapped-params//BS-L.Rdata")
rm(g.sim, s.sim, m.sim, F.sim, params, pm, i, j, later, Scenario, sim.list, sims, BS, L.sim)


#### M7: Boot Lambda ####

###
# Lambda with seed survival averages
###
B.L$Subplot <- factor(B.L$Subplot, levels = c("N", "G", "T"))
B.L$Subplot <- revalue(B.L$Subplot, c("N" = "No Grass", "G" = "Grass", "T" = "Thatch"))
hist(B.L$L)
hist(log(B.L$L))
hist(log(B.L$L + 1))
B.L$log.L <- log(B.L$L + 1)
m5.t <- lm(log.L ~ Subplot * PC.F, B.L)
plot(fitted(m5.t), resid(m5.t))
qqnorm(resid(m5.t))
qqline(resid(m5.t), col = 2, lwd = 2, lty = 2) 
summary(m5.t)
hist(resid(m5.t))



 #### Plot Germination ####

germ.t.sum <- summarySE(dem, measurevar = "p.germ", groupvars = c("Subplot","Species","PC.F"), na.rm = T)

ggplot(germ.t.sum, aes(x = PC.F, y = p.germ, col = Subplot, group = Subplot)) +
  geom_point() +
  theme_classic() +
  theme(legend.title = element_blank(), 
        axis.text = element_text(size = 10), 
        axis.title = element_text(size = 15), 
        strip.text = element_text(size = 15),
        legend.text = element_text(size = 11),
        legend.position = c(.85, .84),
        legend.key.size = unit(1.4, 'lines'),
        axis.line = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  labs(y = "Germination Rate", x = "Drought Tolerance") +
  geom_smooth(method = "lm", se = F) + 
  ylim(.1,1) +
  scale_color_viridis_d() 

#### Plot Mortality ####

dem.sum <- summarySE(dem, measurevar = "p.mort", groupvars = c("Subplot","Species","PC.F"), na.rm = T)

ggplot(dem.sum, aes(x = PC.F, y = p.mort, col = Subplot, group = Subplot)) +
  geom_point() +
  theme_classic() +
  theme(legend.title = element_blank(), 
        axis.text = element_text(size = 10), 
        axis.title = element_text(size = 15), 
        strip.text = element_text(size = 15),
        legend.text = element_text(size = 11),
        legend.position = c(.15, .78),
        legend.key.size = unit(1.4, 'lines'),
        axis.line = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  labs(y = "Moratlity Rate", x = "Drought Tolerance") +
  geom_smooth(method = "lm", se = F) + 
  scale_color_viridis_d()


#### Plot Seed set ####

flo.seed.sum <- summarySE(flo.seed, measurevar = "n.seed.ind", groupvars = c("Subplot","Species","PC.F"), na.rm = T)

ggplot(flo.seed.sum, aes(x = PC.F, y = n.seed.ind, col = Subplot, group = Subplot)) +
  geom_point() +
  theme_classic() +
  theme(legend.title = element_blank(), 
        axis.text = element_text(size = 10), 
        axis.title = element_text(size = 15), 
        strip.text = element_text(size = 15),
        legend.text = element_text(size = 11),
        legend.position = c(.85, .78),
        legend.key.size = unit(1.4, 'lines'),
        axis.line = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  labs(y = "Seed Set", x = "Drought Tolerance") +
  geom_smooth(method = "lm", se = F) + 
  scale_color_viridis_d()

#### Plot Lambda (plot) ####

full.t.sum <- summarySE(full, measurevar = "L", groupvars = c("Subplot","Species","PC.F"), na.rm = T)

ggplot(full.t.sum, aes(y = L, x = PC.F, col = Subplot, group = Subplot)) +
  geom_smooth(method = "lm", se = F) + 
  geom_point() +
  theme_classic() +
  theme(legend.title = element_blank(), 
        axis.text = element_text(size = 10), 
        axis.title = element_text(size = 15), 
        strip.text = element_text(size = 15),
        legend.text = element_text(size = 11),
        legend.position = c(.85, .78),
        legend.key.size = unit(1.4, 'lines'),
        axis.line = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  labs(y = "Population Growth Rate", x = "Drought Tolerance") +
  scale_color_viridis_d()

#### Plot Lambda (BS) ####
ggplot(B.L, aes(y = L, x = PC.F, col = Subplot, group = Subplot)) +
  geom_smooth(method = "lm", se = F) + 
  geom_point() +
  theme_classic() +
  theme(legend.title = element_blank(), 
        axis.text = element_text(size = 10), 
        axis.title = element_text(size = 15), 
        strip.text = element_text(size = 15),
        legend.text = element_text(size = 11),
        legend.position = c(.85, .78),
        legend.key.size = unit(1.4, 'lines'),
        axis.line = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  labs(y = "Population Growth Rate", x = "Drought Tolerance") +
  scale_color_viridis_d()
