### Thatch effects ###
rm(list = ls())
#### Load Libraries ####

library(plyr)
library(dplyr)
library(sjstats)
library(Rmisc)
library(ggplot2)
library(lme4)
library(lmerTest)
library(reshape2)
library(emmeans)
library(multcomp)

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
#flo.seed <- flo.seed[,-16]
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
trait.w$strat <- ifelse(trait.w$PC.F > 0, "ST", "SA")

# merge datasets with traits
dem <- merge(dem, trait.w[,c(1,26)], by = "Species")
sb <- merge(sb, trait.w[,c(1,26)], by = "Species")
flo.seed <- merge(flo.seed, trait.w[,c(1,26)], by = "Species")

rm(PCA.F, PCA.G, PCA.s13, sla.13, grass)

# merge datasets with grass cover
#dem <- merge(dem, grass, all.x = T)


#### Prep: Plot-level Lambda ####

# full <- merge(dem, flo.seed[,-c(6:13)], by = c("Year","Plot","Treat.Code", "Subplot", "Species", "strat"), all = T)
# 
# full <- merge(full, sb, by = c("Year","Plot","Treat.Code","Species","strat"), all.x = T)
# 
# full$p.germ <- full$germ.tot/full$viable
# 
# # calculate L according: L = s*(1-g) + g*(1-m)*F
# full$L.sb <- full$p.surv*(1-full$p.germ)
# full$L.sa <- full$p.germ * (1 -  full$p.mort)
# full$L.seeds <- full$L.sa * full$n.seed.ind
# full$L <- full$L.sb + full$L.seeds
# full$L <- ifelse(full$L.sa == 0, full$L.sb, full$L)
# full$L <- ifelse(full$p.germ == 0, full$L.sb, full$L)

full <- merge(dem, flo.seed[,-c(6:9,11:12)], by = c("Year","Plot","Treat.Code", "Subplot","Species", "strat"), all = T)

# replace NAs in n.seed.inf with species averages
flo.seed.sum <- summarySE(flo.seed, measurevar = "avg.seed.inf", groupvars = "Species", na.rm = T)

for(j in unique(full[is.na(full$avg.seed.inf),]$Species)) {
    full[is.na(full$avg.seed.inf) & full$Species == j,]$avg.seed.inf <- flo.seed.sum[flo.seed.sum$Species == j,]$avg.seed.inf
  }

# replace NAs in avg.flo with species averages
flo.seed.sum <- summarySE(flo.seed, measurevar = "avg.flo", groupvars = "Species", na.rm = T)

for(j in unique(full[is.na(full$avg.flo),]$Species)) {
    full[is.na(full$avg.flo) & full$Species == j,]$avg.flo <- flo.seed.sum[flo.seed.sum$Species == j,]$avg.flo
  }

# Put viability estimates back in
flo.seed.sum <- summarySE(flo.seed, measurevar = "p.viable", groupvars = "Species", na.rm = T)

for(j in unique(full[is.na(full$p.viable),]$Species)) {
    full[is.na(full$p.viable) & full$Species == j,]$p.viable <- flo.seed.sum[flo.seed.sum$Species == j,]$p.viable
  }

# Multiply together to get n.seed.ind
full$n.seed.ind <- full$avg.seed.inf * full$avg.flo * full$p.viable

# merge with seed bank data
full <- merge(full, sb, by = c("Year","Plot","Treat.Code","Species", "strat"), all = T)

# calculate L according: L = s*(1-g) + g*(1-m)*F
full$L.sb <- full$p.surv*(1-full$p.germ)
full$L.sa <- full$p.germ * (1 -  full$p.mort)
full$L.seeds <- full$L.sa * full$n.seed.ind
full$L <- full$L.sb + full$L.seeds
full$L <- ifelse(full$L.sa == 0, full$L.sb, full$L)
full$L <- ifelse(is.na(full$L.sa) == T, full$L.sb, full$L)

#### M1: Germination ####
dem$Subplot2 <- ifelse(dem$Subplot == "Grass", "No Grass", as.character(dem$Subplot))
m1.t <- glmer(cbind(germ.tot, viable-germ.tot) ~ Subplot2 * strat + (1|Plot:Species), family = binomial, data = dem, glmerControl(calc.derivs = F))
plot(fitted(m1.t), resid(m1.t))
summary(m1.t) 
overdisp(m1.t)

fixef(m1.t)
# contrasts
SA.N <- c(1, 0, 0, 0, 0, 0)
SA.G <- c(1, 1, 0, 0, 0, 0)
SA.T <- c(1, 0, 1, 0, 0, 0)

ST.N <- c(1, 0, 0, 1, 0, 0)
ST.G <- c(1, 1, 0, 1, 1, 0)
ST.T <- c(1, 0, 1, 1, 0, 1)

K <- rbind("ST.G - ST.N" = ST.G - ST.N,
           "ST.T - ST.N" = ST.T - ST.N,
           "ST.T - ST.G" = ST.T - ST.G,
           "SA.G - SA.N" = SA.G - SA.N,
           "SA.T - SA.N" = SA.T - SA.N,
           "SA.T - SA.G" = SA.T - SA.G)

summary(glht(m1.t, linfct = K), test = adjusted("BH"))

contrast(emmeans(m1.t, ~ Subplot | strat), interaction = "revpairwise", adjust = "FDR")

contrast(emmeans(m1.t, ~ strat:Subplot), interaction = "revpairwise", adjust = "FDR")

contrast(emmeans(m1.t, ~ strat | Subplot), interaction = "revpairwise", adjust = "FDR")

#### M2: Mortality ####
m2.t <- glmer(cbind(tot.mort, germ.proj-tot.mort) ~ Subplot * strat + (1|Plot:Species), family = binomial, data = dem, glmerControl(calc.derivs = F))
plot(fitted(m2.t), resid(m2.t))
summary(m2.t) # thatch increases mortality, still higher mortality in Tolerators; no interaction between current grass/drought tolerance or thatch/drought tolerance

summary(glht(m2.t, linfct = K), test = adjusted("BH"))

contrast(emmeans(m2.t, ~ Subplot | strat), interaction = "revpairwise", adjust = "FDR")

contrast(emmeans(m2.t, ~ strat:Subplot), interaction = "revpairwise", adjust = "FDR")

contrast(emmeans(m2.t, ~ strat | Subplot), interaction = "revpairwise", adjust = "FDR")

#### M3: Seed set ####

# random effect structure too complex with slope for sublot; missing 20 levels because of low germ/high mort, with the nestedness structure of the previous models there are too many levels of random effects for estimates

m3.t <- lmer(log(n.seed.ind + 1) ~ Subplot * strat + (1|Plot:Species), flo.seed)
plot(fitted(m3.t), resid(m3.t))
qqnorm(resid(m3.t))
qqline(resid(m3.t), col = 2, lwd = 2, lty = 2)
summary(m3.t)

summary(glht(m3.t, linfct = K), test = adjusted("BH"))

contrast(emmeans(m3.t, ~ Subplot | strat), interaction = "revpairwise", adjust = "FDR")

contrast(emmeans(m3.t, ~ strat:Subplot), interaction = "revpairwise", adjust = "FDR")

contrast(emmeans(m3.t, ~ strat | Subplot), interaction = "revpairwise", adjust = "FDR")

# Though grass lowers seed set, the majority of effects of grass occurs through thatch (by lowering light levels?); these species are adapted to high light/lower moisture environments; thatch lowers germination in all species, increases mortality in all species and lowers seed set in all species but moreso in drought avoiders than drought tolerators; I can see the mechanism for germination (light levels?) but how is thatch increasing mortality? light still? it cant be watering in this type of year. 

#### M4: Plot Lambda ####
hist(log(full$L + .5))

m4.t <- lmer(log(L + .7) ~ Subplot * strat + (1|Plot:Species), full)
plot(fitted(m4.t), resid(m4.t))
qqnorm(resid(m4.t))
qqline(resid(m4.t), col = 2, lwd = 2, lty = 2)
summary(m4.t)

summary(glht(m4.t, linfct = K), test = adjusted("BH"))

contrast(emmeans(m4.t, ~ Subplot | strat), interaction = "revpairwise", adjust = "FDR")

contrast(emmeans(m4.t, ~ strat:Subplot), interaction = "revpairwise", adjust = "FDR")

contrast(emmeans(m4.t, ~ strat | Subplot), interaction = "revpairwise", adjust = "FDR")

# Is it due to live grass?
m4.t <- lmer(log(L + .3) ~ Subplot * strat + (1|Plot:Species), full[full$Subplot != "Thatch",])
plot(fitted(m4.t), resid(m4.t))
qqnorm(resid(m4.t))
qqline(resid(m4.t), col = 2, lwd = 2, lty = 2)
summary(m4.t)

fixef(m4.t)

# contrasts
SA.N <- c(1, 0, 0, 0)
SA.G <- c(1, 1, 0, 0)

ST.N <- c(1, 0, 1, 0)
ST.G <- c(1, 1, 1, 1)

K2 <- rbind("ST.G - ST.N" = ST.G - ST.N,
           "SA.G - SA.N" = SA.G - SA.N)

summary(glht(m4.t, linfct = K2), test = adjusted("BH"))

contrast(emmeans(m4.t, ~ strat:Subplot), interaction = "revpairwise", adjust = "FDR")

#### M5: Grass Cover ####
hist(sqrt(grass$cover.g))
m5 <- lm(sqrt(cover.g) ~ Subplot, grass)
plot(fitted(m5), resid(m5))
qqnorm(resid(m5))
qqline(resid(m5), col = 2, lwd = 2, lty = 2)
summary(m5)

#### Plot Germination ####

germ.t.sum <- summarySE(dem, measurevar = "p.germ", groupvars = c("Subplot","strat"), na.rm = T)
germ.t.sum$strat <- revalue(germ.t.sum$strat, c("SA" = "Stress Avoider", "ST" = "Stress Tolerator"))
germ.t.sum$Subplot <- factor(germ.t.sum$Subplot, levels = c("No Grass", "Grass", "Thatch"))
germ.t.sum$Subplot <- revalue(germ.t.sum$Subplot, c("Thatch" = "Grass + Thatch"))

ggplot(germ.t.sum, aes(x = strat, y = p.germ, fill = Subplot)) +
  geom_bar(position = "dodge", stat = "identity") +
  geom_errorbar(aes(ymin = p.germ - se, ymax = p.germ + se, width = 0.2), position = position_dodge(0.9)) + 
  theme_classic() +
  theme(legend.title = element_blank(), 
        axis.text.y = element_text(size = 10), 
        axis.text.x = element_text(size = 13), 
        axis.title = element_text(size = 15), 
        axis.title.x = element_blank(),
        strip.text = element_text(size = 15),
        legend.text = element_text(size = 11),
        legend.position = c(.8, .84),
        legend.key.size = unit(1.4, 'lines'),
        axis.line = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  labs(y = "Germination Rate", x = "Subplot") +
  ylim(0,.85) +
  scale_fill_viridis_d()

#### Plot Mortality ####

dem.sum <- summarySE(dem, measurevar = "p.mort", groupvars = c("Subplot","strat"), na.rm = T)

dem.sum$Subplot <- factor(dem.sum$Subplot, levels = c("No Grass", "Grass", "Thatch"))
dem.sum$Subplot <- revalue(dem.sum$Subplot, c("Thatch" = "Grass + Thatch"))
dem.sum$strat <- revalue(dem.sum$strat, c("SA" = "Stress Avoider", "ST" = "Stress Tolerator"))

ggplot(dem.sum, aes(x = strat, y = p.mort, fill = Subplot)) +
  geom_bar(position = "dodge", stat = "identity") +
  geom_errorbar(aes(ymin = p.mort - se, ymax = p.mort + se, width = 0.2), position = position_dodge(0.9)) + 
  theme_classic() +
  theme(legend.title = element_blank(), 
        axis.text.y = element_text(size = 10), 
        axis.text.x = element_text(size = 13), 
        axis.title = element_text(size = 15), 
        axis.title.x = element_blank(),
        strip.text = element_text(size = 15),
        legend.text = element_text(size = 11),
        legend.position = c(.2, .8),
        legend.key.size = unit(1.4, 'lines'),
        axis.line = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  labs(y = "Mortality Rate", x = "Subplot") +
  scale_fill_viridis_d()


#### Plot Seed set ####

flo.seed.sum <- summarySE(flo.seed, measurevar = "n.seed.ind", groupvars = c("Subplot","strat"), na.rm = T)
flo.seed.sum$Subplot <- factor(flo.seed.sum$Subplot, levels = c("No Grass", "Grass", "Thatch"))
flo.seed.sum$Subplot <- revalue(flo.seed.sum$Subplot, c("Thatch" = "Grass + Thatch"))
flo.seed.sum$strat <- revalue(flo.seed.sum$strat, c("SA" = "Stress Avoider", "ST" = "Stress Tolerator"))

ggplot(flo.seed.sum, aes(x = strat, y = n.seed.ind, fill = Subplot)) +
  geom_bar(position = "dodge", stat = "identity") +
  geom_errorbar(aes(ymin = n.seed.ind - se, ymax = n.seed.ind + se, width = 0.2), position = position_dodge(0.9)) + 
  theme_classic() +
  theme(legend.title = element_blank(), 
        axis.text.y = element_text(size = 10), 
        axis.text.x = element_text(size = 13), 
        axis.title = element_text(size = 15), 
        axis.title.x = element_blank(),
        strip.text = element_text(size = 15),
        legend.text = element_text(size = 11),
        legend.position = c(.8, .84),
        legend.key.size = unit(1.4, 'lines'),
        axis.line = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  labs(y = "Seeds per individual", x = "Subplot") +
  scale_fill_viridis_d()


#### Plot Lambda (plot) ####

full.t.sum <- summarySE(full, measurevar = "L", groupvars = c("Subplot","strat"), na.rm = T)

full.t.sum$strat <- revalue(full.t.sum$strat, c("SA" = "Stress Avoider", "ST" = "Stress Tolerator"))
full.t.sum$Subplot <- factor(full.t.sum$Subplot, levels = c("No Grass", "Grass", "Thatch"))
full.t.sum$Subplot <- revalue(full.t.sum$Subplot, c("Thatch" = "Grass + Thatch"))

# ggplot(full.t.sum, aes(x = Subplot, y = L, fill = strat)) +
#   geom_bar(position = "dodge", stat = "identity") +
#   geom_errorbar(aes(ymin = L - se, ymax = L + se, width = .2), position = position_dodge(.9)) +
#   theme_classic() +
#   theme(legend.title = element_blank(), 
#         axis.text.y = element_text(size = 10), 
#         axis.text.x = element_text(size = 13), 
#         axis.title = element_text(size = 15), 
#         axis.title.x = element_blank(),
#         strip.text = element_text(size = 15),
#         legend.text = element_text(size = 11),
#         legend.position = c(.8, .84),
#         legend.key.size = unit(1.4, 'lines'),
#         axis.line = element_blank(),
#         panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
#   labs(y = "Lambda", x = "Subplot") +
#   scale_fill_viridis_d()

ggplot(full.t.sum, aes(x = strat, y = L, fill = Subplot)) +
  geom_bar(position = "dodge", stat = "identity") +
  geom_errorbar(aes(ymin = L - se, ymax = L + se, width = 0.2), position = position_dodge(.9)) + 
  theme_classic() +
  theme(legend.title = element_blank(), 
        axis.text.y = element_text(size = 10), 
        axis.text.x = element_text(size = 13), 
        axis.title = element_text(size = 15), 
        axis.title.x = element_blank(),
        strip.text = element_text(size = 15),
        legend.text = element_text(size = 11),
        legend.position = c(.8, .84),
        legend.key.size = unit(1.4, 'lines'),
        axis.line = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  labs(y = "Lambda", x = "Subplot") +
  scale_fill_viridis_d()
