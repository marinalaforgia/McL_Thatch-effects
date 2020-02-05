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
library(glmmTMB)
library(arm)
library(tidyverse)
library(car)

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

biplot(PCA.F)
summary(PCA.F)

trait.w$PC.G <- PCA.G$x[,1] 
trait.w$PC.F <- PCA.F$x[,1] 
trait.w$PC.s13 <- PCA.s13$x[c(2,3,5,6,7,9),1] 
trait.w$strat <- ifelse(trait.w$PC.F > 0, "ST", "SA")

# merge datasets with traits
dem <- merge(dem, trait.w[,c(1,26)], by = "Species")
sb <- merge(sb, trait.w[,c(1,26)], by = "Species")
flo.seed <- merge(flo.seed, trait.w[,c(1,26)], by = "Species")

rm(PCA.F, PCA.G, PCA.s13, sla.13)

# merge datasets with grass cover
#dem <- merge(dem, grass, all.x = T)


#### Prep: Plot-level Lambda ####

full <- merge(dem, flo.seed[,-c(6:13)], by = c("Year","Plot","Treat.Code", "Subplot", "Species", "strat"), all = T)

full <- merge(full, sb, by = c("Year","Plot","Treat.Code","Species","strat"), all.x = T)

full$p.germ <- full$germ.tot/full$viable

# calculate L according: L = s*(1-g) + g*(1-m)*F
full$L.sb <- full$p.surv*(1-full$p.germ)
full$L.sa <- full$p.germ * (1 -  full$p.mort)
full$L.seeds <- full$L.sa * full$n.seed.ind
full$L <- full$L.sb + full$L.seeds
full$L <- ifelse(full$L.sa == 0, full$L.sb, full$L)
full$L <- ifelse(full$p.germ == 0, full$L.sb, full$L)

# full <- merge(dem, flo.seed[,-c(6:9,11:12)], by = c("Year","Plot","Treat.Code", "Subplot","Species", "strat"), all = T)
# 
# # replace NAs in n.seed.inf with species averages
# flo.seed.sum <- summarySE(flo.seed, measurevar = "avg.seed.inf", groupvars = "Species", na.rm = T)
# 
# for(j in unique(full[is.na(full$avg.seed.inf),]$Species)) {
#     full[is.na(full$avg.seed.inf) & full$Species == j,]$avg.seed.inf <- flo.seed.sum[flo.seed.sum$Species == j,]$avg.seed.inf
#   }
# 
# # replace NAs in avg.flo with species averages
# flo.seed.sum <- summarySE(flo.seed, measurevar = "avg.flo", groupvars = "Species", na.rm = T)
# 
# for(j in unique(full[is.na(full$avg.flo),]$Species)) {
#     full[is.na(full$avg.flo) & full$Species == j,]$avg.flo <- flo.seed.sum[flo.seed.sum$Species == j,]$avg.flo
#   }
# 
# # Put viability estimates back in
# flo.seed.sum <- summarySE(flo.seed, measurevar = "p.viable", groupvars = "Species", na.rm = T)
# 
# for(j in unique(full[is.na(full$p.viable),]$Species)) {
#     full[is.na(full$p.viable) & full$Species == j,]$p.viable <- flo.seed.sum[flo.seed.sum$Species == j,]$p.viable
#   }
# 
# # Multiply together to get n.seed.ind
# full$n.seed.ind <- full$avg.seed.inf * full$avg.flo * full$p.viable
# 
# # merge with seed bank data
# full <- merge(full, sb, by = c("Year","Plot","Treat.Code","Species", "strat"), all = T)
# 
# # calculate L according: L = s*(1-g) + g*(1-m)*F
# full$L.sb <- full$p.surv*(1-full$p.germ)
# full$L.sa <- full$p.germ * (1 -  full$p.mort)
# full$L.seeds <- full$L.sa * full$n.seed.ind
# full$L <- full$L.sb + full$L.seeds
# full$L <- ifelse(full$L.sa == 0, full$L.sb, full$L)
# full$L <- ifelse(is.na(full$L.sa) == T, full$L.sb, full$L)

### Calculate elasticity ###
# full$elas.g <- abs(full$p.germ*((1-full$p.mort)*full$n.seed.ind - full$p.surv)/full$L)
# full$elas.F <- abs((full$p.germ*(1-full$p.mort)*full$n.seed.ind)/full$L)

#### M1: Germination ####
dem$Subplot2 <- ifelse(dem$Subplot == "Grass", "No Grass", as.character(dem$Subplot))
# m1.t <- glmer(cbind(germ.tot, viable-germ.tot) ~ Subplot2 * strat + (1|Plot:Species), family = binomial, data = dem, glmerControl(calc.derivs = F))
# plot(fitted(m1.t), resid(m1.t))
# summary(m1.t) 
# (OD <- overdisp(m1.t)$ratio)
# (new.p <- as.data.frame(2*pnorm(-abs(fixef(m1.t)/(se.fixef(m1.t)*sqrt(OD))))))
m1.t1 <- glmmTMB(cbind(germ.tot, viable-germ.tot) ~ Subplot2 + strat + (1|Plot:Species), family = betabinomial(link = "logit"), dem) 
m1.t2 <- glmmTMB(cbind(germ.tot, viable-germ.tot) ~ Subplot2 * strat + (1|Plot:Species), family = betabinomial(link = "logit"), dem) 

anova(m1.t1, m1.t2)

summary(m1.t2)

# contrasts
SA.N <- c(1, 0, 0, 0)
SA.T <- c(1, 1, 0, 0)

ST.N <- c(1, 0, 1, 0)
ST.T <- c(1, 1, 1, 1)

K <- rbind("ST.T - ST.N" = ST.T - ST.N,
           "SA.T - SA.N" = SA.T - SA.N)

# workaround to get glmmTMB to work with glht; https://rdrr.io/cran/glmmTMB/f/vignettes/model_evaluation.rmd
glht_glmmTMB <- function (model, ..., component="cond") {
    glht(model, ...,
         coef. = function(x) fixef(x)[[component]],
         vcov. = function(x) vcov(x)[[component]],
         df = NULL)
}
modelparm.glmmTMB <- function (model, coef. = function(x) fixef(x)[[component]],
                               vcov. = function(x) vcov(x)[[component]],
                               df = NULL, component="cond", ...) {
    multcomp:::modelparm.default(model, coef. = coef., vcov. = vcov.,
                        df = df, ...)
}

summary(glht(m1.t2, linfct = K), test = adjusted("BH"))  

# contrast(emmeans(m1.t, ~ Subplot | strat), interaction = "revpairwise", adjust = "FDR")
# 
# contrast(emmeans(m1.t, ~ strat:Subplot), interaction = "revpairwise", adjust = "FDR")
# 
# contrast(emmeans(m1.t, ~ strat | Subplot), interaction = "revpairwise", adjust = "FDR")

#### M2: Mortality ####
#dem$Subplot <- factor(dem$Subplot, levels = c("Grass", "No Grass", "Thatch"))
m2.t1 <- glmer(cbind(tot.mort, germ.proj-tot.mort) ~ Subplot + strat + (1|Plot:Species), family = binomial, data = dem, glmerControl(calc.derivs = F))
m2.t2 <- glmer(cbind(tot.mort, germ.proj-tot.mort) ~ Subplot * strat + (1|Plot:Species), family = binomial, data = dem, glmerControl(calc.derivs = F))

anova(m2.t1, m2.t2)

plot(fitted(m2.t1), resid(m2.t1))
summary(m2.t1)

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

summary(glht(m2.t, linfct = K), test = adjusted("BH"))

# contrast(emmeans(m2.t, ~ Subplot | strat), interaction = "revpairwise", adjust = "FDR")
# 
# contrast(emmeans(m2.t, ~ strat:Subplot), interaction = "revpairwise", adjust = "FDR")
# 
# contrast(emmeans(m2.t, ~ strat | Subplot), interaction = "revpairwise", adjust = "FDR")

#### M3: Seed set ####
m3.t1 <- lmer(log(n.seed.ind + 1) ~ Subplot + strat + (1|Plot:Species), flo.seed)
m3.t2 <- lmer(log(n.seed.ind + 1) ~ Subplot * strat + (1|Plot:Species), flo.seed)
anova(m3.t1, m3.t2)

plot(fitted(m3.t2), resid(m3.t2))
qqnorm(resid(m3.t2))
qqline(resid(m3.t2), col = 2, lwd = 2, lty = 2)
summary(m3.t2)

summary(glht(m3.t2, linfct = K), test = adjusted("BH"))

# contrast(emmeans(m3.t, ~ Subplot | strat), interaction = "revpairwise", adjust = "FDR")
# 
# contrast(emmeans(m3.t, ~ strat:Subplot), interaction = "revpairwise", adjust = "FDR")
# 
# contrast(emmeans(m3.t, ~ strat | Subplot), interaction = "revpairwise", adjust = "FDR")

#### M4: Plot Lambda ####
hist(log(full$L + .7))
m4.t1 <- lmer(log(L + .7) ~ Subplot + strat + (1|Plot:Species), full)
m4.t2 <- lmer(log(L + .7) ~ Subplot * strat + (1|Plot:Species), full)
anova(m4.t1, m4.t2)

plot(fitted(m4.t2), resid(m4.t2))
qqnorm(resid(m4.t2))
qqline(resid(m4.t2), col = 2, lwd = 2, lty = 2)
summary(m4.t2)

summary(glht(m4.t2, linfct = K), test = adjusted("BH"))

contrast(emmeans(m4.t2, ~ Subplot | strat), interaction = "revpairwise", adjust = "FDR")

contrast(emmeans(m4.t2, ~ strat:Subplot), interaction = "revpairwise", adjust = "FDR")

contrast(emmeans(m4.t2, ~ strat | Subplot), interaction = "revpairwise", adjust = "FDR")

#### M5: Grass Cover ####
hist(sqrt(grass$cover.g))
m5 <- lm(sqrt(cover.g) ~ Subplot, grass)
plot(fitted(m5), resid(m5))
qqnorm(resid(m5))
qqline(resid(m5), col = 2, lwd = 2, lty = 2)
summary(m5)

#### Plot RGR v WUE ####
trait.w$Species.short <- revalue(trait.w$Species, c("Agoseris heterophylla" =  "AGHE", "Calycadenia pauciflora" = "CAPA", "Clarkia purpurea" = "CLPU", "Hemizonia congesta" = "HECO", "Lasthenia californica" = "LACA", "Plantago erecta" = "PLER"))
trait.w$strat <- revalue(trait.w$strat, c("SA" = "Acquisitive", "ST" = "Conservative"))

plot.trait <- ggplot(trait.w, aes(x = D13C.F, y = RGR.la.F, col = strat)) +
  theme_classic() +
  geom_smooth(method = "lm", se = F, col = "black", size = .7) +
  geom_point(size = 1.2) +
  geom_label(aes(label = Species.short), size = 2, nudge_x = .7, label.padding = unit(0.15, "lines"), col = "black") +
  theme(axis.text = element_text(size = 8),
        axis.title = element_text(size = 10),
        legend.position = c(0.78, 0.97),
        legend.text = element_text(size = 8),
        legend.title = element_blank(),
        legend.spacing.y = unit(1.0, 'cm'),
        legend.key.size = unit(0.5, 'cm'),
        legend.margin=margin(t = 0, unit='cm')) +
  scale_x_reverse(lim = c(25, 19.5)) +
  scale_color_manual(values = c("#009E73", "#E69F00")) +
  labs(x = "Carbon isotope discrimination\n(∆, \u2030)", y = expression(atop("Relative growth rate", paste("(", cm^{2}, ")"%.%"(", cm^{2}, ")"^{-1}%.%day^{-1})))) 

ggsave(plot.trait, filename = "Figures/trait.tiff", width = 3, height = 2.5, units = "in", dpi = 600)

#### Plot Germination ####

germ.t.sum <- summarySE(dem, measurevar = "p.germ", groupvars = c("Subplot2", "strat"), na.rm = T)
germ.t.sum$strat <- revalue(germ.t.sum$strat, c("SA" = "Acquisitive", "ST" = "Conservative"))
germ.t.sum$Subplot2 <- revalue(germ.t.sum$Subplot, c("Thatch" = "Litter addition", "No Grass" = "No grass"))
germ.t.sum$Subplot2 <- factor(germ.t.sum$Subplot2, levels = c("No grass", "Litter addition"))

plot.germ <- ggplot(germ.t.sum, aes(x = strat, y = p.germ, fill = Subplot2)) +
  geom_bar(position = "dodge", stat = "identity") +
  geom_errorbar(aes(ymin = p.germ - se, ymax = p.germ + se, width = 0.2), position = position_dodge(0.9)) + 
  theme_classic() +
  theme(legend.title = element_blank(), 
        axis.text = element_text(size = 10, color = "black"), 
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 12), 
        axis.title.x = element_blank(),
        legend.text = element_text(size = 9),
        legend.position = c(.77, .84),
        legend.key.size = unit(1.2, 'lines')) +
  labs(y = "Germination", x = "Subplot") +
  ylim(0,.85) +
  scale_fill_viridis_d()

ggsave(plot.germ, filename = "Figures/Germ.tiff", width = 3, height = 2.7, units = "in", dpi = 600)

#### Plot Mortality ####

dem.sum <- summarySE(dem, measurevar = "p.mort", groupvars = c("Subplot","strat"), na.rm = T)

dem.sum$Subplot <- factor(dem.sum$Subplot, levels = c("No Grass", "Grass", "Thatch"))
dem.sum$Subplot <- revalue(dem.sum$Subplot, c("Thatch" = "Litter addition", "No Grass" = "No grass"))
dem.sum$strat <- revalue(dem.sum$strat,  c("SA" = "Acquisitive", "ST" = "Conservative"))

plot.mort <- ggplot(dem.sum, aes(x = strat, y = p.mort, fill = Subplot)) +
  geom_bar(position = "dodge", stat = "identity") +
  geom_errorbar(aes(ymin = p.mort - se, ymax = p.mort + se, width = 0.2), position = position_dodge(0.9)) + 
  theme_classic() +
  theme(legend.title = element_blank(), 
        axis.text.y = element_text(size = 10), 
        axis.text.x = element_text(size = 10, color = "black"), 
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 12), 
        axis.title.x = element_blank(),
        legend.text = element_text(size = 9),
        legend.position = c(.26, .8),
        legend.key.size = unit(1.2, 'lines')) +
  labs(y = "Mortality", x = "Subplot") +
  scale_fill_viridis_d()

ggsave(plot.mort, filename = "Figures/Mort.tiff", width = 3, height = 2.7, units = "in", dpi = 600)

#### Plot Seed set ####

flo.seed.sum <- summarySE(flo.seed, measurevar = "n.seed.ind", groupvars = c("Subplot","strat"), na.rm = T)
flo.seed.sum$Subplot <- factor(flo.seed.sum$Subplot, levels = c("No Grass", "Grass", "Thatch"))
flo.seed.sum$Subplot <- revalue(flo.seed.sum$Subplot, c("Thatch" = "Litter addition", "No Grass" = "No grass"))
flo.seed.sum$strat <- revalue(flo.seed.sum$strat,  c("SA" = "Acquisitive", "ST" = "Conservative"))

plot.seed <- ggplot(flo.seed.sum, aes(x = strat, y = n.seed.ind, fill = Subplot)) +
  geom_bar(position = "dodge", stat = "identity") +
  geom_errorbar(aes(ymin = n.seed.ind - se, ymax = n.seed.ind + se, width = 0.2), position = position_dodge(0.9)) + 
  theme_classic() +
  theme(legend.title = element_blank(), 
        axis.text.y = element_text(size = 10), 
        axis.text.x = element_text(size = 10, color = "black"), 
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 12), 
        axis.title.x = element_blank(),
        legend.text = element_text(size = 9),
        legend.position = c(.75, .8),
        legend.key.size = unit(1.2, 'lines')) +
  labs(y = "Seeds per individual", x = "Subplot") +
  scale_fill_viridis_d()

ggsave(plot.seed, filename = "Figures/Seed.tiff", width = 3, height = 2.7, units = "in", dpi = 600)

#### Plot Lambda (plot) ####

full.t.sum <- summarySE(full, measurevar = "L", groupvars = c("Subplot","strat"), na.rm = T)

full.t.sum$strat <- revalue(full.t.sum$strat, c("SA" = "Acquisitive", "ST" = "Conservative"))
full.t.sum$Subplot <- factor(full.t.sum$Subplot, levels = c("No Grass", "Grass", "Thatch"))
full.t.sum$Subplot <- revalue(full.t.sum$Subplot, c("Thatch" = "Litter addition", "No Grass" = "No grass"))

# ggplot(full.t.sum, aes(x = Subplot, y = L, fill = strat)) +
#   geom_bar(position = "dodge", stat = "identity") +
#   geom_errorbar(aes(ymin = L - se, ymax = L + se, width = .2), position = position_dodge(.9)) +
#   theme_classic() +
#   theme(legend.title = element_blank(), 
#         axis.text.y = element_text(size = 10), 
#         axis.text.x = element_text(size = 13, color = "black"), 
#         axis.title = element_text(size = 15), 
#         axis.title.x = element_blank(),
#         legend.text = element_text(size = 11),
#         legend.position = c(.8, .84),
#         legend.key.size = unit(1.4, 'lines'),
#         axis.line = element_blank(),
#         panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
#   labs(y = "Lambda", x = "Subplot") +
#   scale_fill_viridis_d()

plot.lambda <- ggplot(full.t.sum, aes(x = strat, y = L, fill = Subplot)) +
  geom_bar(position = "dodge", stat = "identity") +
  geom_errorbar(aes(ymin = L - se, ymax = L + se, width = 0.2), position = position_dodge(.9)) + 
  theme_classic() +
  theme(legend.title = element_blank(), 
        axis.text.y = element_text(size = 10), 
        axis.text.x = element_text(size = 10, color = "black"), 
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 12), 
        axis.title.x = element_blank(),
        legend.text = element_text(size = 9),
        legend.position = c(.75, .8),
        legend.key.size = unit(1.2, 'lines')) +
  labs(y = "Lambda", x = "Subplot") +
  scale_fill_viridis_d()

ggsave(plot.lambda, filename = "Figures/Lambda.tiff", width = 3, height = 2.7, units = "in", dpi = 600)

#### Spp: Lambda ####

full.sum.spp <- summarySE(full, measurevar = "L", groupvars = c("Subplot", "strat", "Species"), na.rm = T)

full.sum.spp$Subplot <- factor(full.sum.spp$Subplot, levels = c("No Grass", "Grass", "Thatch"))
full.sum.spp$Subplot <- revalue(full.sum.spp$Subplot, c("Thatch" = "Litter addition", "No Grass" = "No grass"))

# Avoiders
plot.lambda.SA <- ggplot(full.sum.spp[full.sum.spp$strat == "SA",], aes(x = Species, y = L, fill = Subplot)) +
  geom_bar(position = "dodge", stat = "identity") +
  geom_errorbar(aes(ymin = L - se, ymax = L + se, width = 0.2), position = position_dodge(.9)) + 
  theme_classic() +
  theme(legend.title = element_blank(), 
        axis.text.y = element_text(size = 10), 
        axis.text.x = element_text(size = 12, color = "black"), 
        axis.ticks.x = element_blank(),
        axis.title = element_text(size = 15), 
        axis.title.x = element_blank(),
        strip.text = element_text(size = 15),
        legend.text = element_text(size = 12),
        legend.position = c(.84, .76),
        legend.key.size = unit(1.4, 'lines')) +
  labs(y = "Lambda", x = "Subplot") +
  ylim(0,40) +
  scale_fill_viridis_d()

ggsave(plot.lambda.SA, filename = "Figures/lambda-sa.tiff", width = 6, height = 3, units = "in", dpi = 600)

# Tolerator
plot.lambda.ST <- ggplot(full.sum.spp[full.sum.spp$strat == "ST",], aes(x = Species, y = L, fill = Subplot)) +
  geom_bar(position = "dodge", stat = "identity") +
  geom_errorbar(aes(ymin = L - se, ymax = L + se, width = 0.2), position = position_dodge(.9)) + 
  theme_classic() +
  theme(legend.title = element_blank(), 
        legend.position = "none",
        axis.text.y = element_text(size = 10), 
        axis.text.x = element_text(size = 12, color = "black"), 
        axis.ticks.x = element_blank(),
        axis.title = element_text(size = 15), 
        axis.title.x = element_blank(),
        strip.text = element_text(size = 15),
        legend.text = element_text(size = 12),
        legend.key.size = unit(1.4, 'lines')) +
  labs(y = "Lambda", x = "Subplot") +
  #scale_y_continuous(breaks=c(0,3,6,9)) +
  ylim(0,40) +
  scale_fill_viridis_d()

ggsave(plot.lambda.ST, filename = "Figures/lambda-st.tiff", width = 6, height = 3, units = "in", dpi = 600)

#### Spp: Germ ####
germ.sum.spp <- summarySE(dem, measurevar = "p.germ", groupvars = c("Subplot2", "strat", "Species"), na.rm = T)
germ.sum.spp$Subplot <- factor(germ.sum.spp$Subplot, levels = c("No Grass", "Grass", "Thatch"))
germ.sum.spp$Subplot <- revalue(germ.sum.spp$Subplot, c("Thatch" = "Litter addition", "No Grass" = "No grass"))

# Avoiders
plot.germ.SA <- ggplot(germ.sum.spp[germ.sum.spp$strat == "SA",], aes(x = Species, y = p.germ, fill = Subplot)) +
  geom_bar(position = "dodge", stat = "identity") +
  geom_errorbar(aes(ymin = p.germ - se, ymax = p.germ + se, width = 0.2), position = position_dodge(0.9)) + 
  theme_classic() +
  theme(legend.title = element_blank(), 
        axis.text.y = element_text(size = 10), 
        axis.text.x = element_text(size = 12, color = "black"), 
        axis.ticks.x = element_blank(),
        axis.title = element_text(size = 15), 
        axis.title.x = element_blank(),
        strip.text = element_text(size = 15),
        legend.text = element_text(size = 12),
        legend.position = c(.86, .88),
        legend.key.size = unit(1, 'lines')) +
  labs(y = "Proportion Germinated", x = "Subplot") +
  ylim(0,1) +
  scale_fill_viridis_d()

ggsave(plot.germ.SA, filename = "Figures/germ-sa.tiff", width = 6, height = 3, units = "in", dpi = 600)

# Tolerators
plot.germ.ST <- ggplot(germ.sum.spp[germ.sum.spp$strat == "ST",], aes(x = Species, y = p.germ, fill = Subplot)) +
  geom_bar(position = "dodge", stat = "identity") +
  geom_errorbar(aes(ymin = p.germ - se, ymax = p.germ + se, width = 0.2), position = position_dodge(0.9)) + 
  theme_classic() +
  theme(legend.title = element_blank(), 
        axis.text.y = element_text(size = 10), 
        axis.text.x = element_text(size = 12, color = "black"), 
        axis.ticks.x = element_blank(),
        axis.title = element_text(size = 15), 
        axis.title.x = element_blank(),
        strip.text = element_text(size = 15),
        legend.text = element_text(size = 12),
        legend.position = "none",
        legend.key.size = unit(1.2, 'lines')) +
  labs(y = "Proportion Germinated", x = "Subplot") +
  ylim(0,1) +
  scale_fill_viridis_d()

ggsave(plot.germ.ST, filename = "Figures/germ-st.tiff", width = 6, height = 3, units = "in", dpi = 600)

#### Spp: Mort ####
mort.sum.spp <- summarySE(dem, measurevar = "p.mort", groupvars = c("Subplot","strat", "Species"), na.rm = T)
mort.sum.spp$Subplot <- factor(mort.sum.spp$Subplot, levels = c("No Grass", "Grass", "Thatch"))
mort.sum.spp$Subplot <- revalue(mort.sum.spp$Subplot, c("Thatch" = "Litter addition", "No Grass" = "No grass"))

# Avoiders
plot.mort.SA <- ggplot(mort.sum.spp[mort.sum.spp$strat == "SA",], aes(x = Species, y = p.mort, fill = Subplot)) +
  geom_bar(position = "dodge", stat = "identity") +
  geom_errorbar(aes(ymin = p.mort - se, ymax = p.mort + se, width = 0.2), position = position_dodge(0.9)) + 
  theme_classic() +
  theme(legend.title = element_blank(), 
        axis.text.y = element_text(size = 10), 
        axis.text.x = element_text(size = 12, color = "black"),
        axis.ticks.x = element_blank(),
        axis.title = element_text(size = 15), 
        axis.title.x = element_blank(),
        strip.text = element_text(size = 15),
        legend.text = element_text(size = 12),
        legend.position = c(.84, .8),
        legend.key.size = unit(1.2, 'lines')) +
  labs(y = "Mortality", x = "Subplot") +
  ylim(0,1) +
  scale_fill_viridis_d()

ggsave(plot.mort.SA, filename = "Figures/mort-sa.tiff", width = 6, height = 3, units = "in", dpi = 600)

# Tolerators
plot.mort.ST <- ggplot(mort.sum.spp[mort.sum.spp$strat == "ST",], aes(x = Species, y = p.mort, fill = Subplot)) +
  geom_bar(position = "dodge", stat = "identity") +
  geom_errorbar(aes(ymin = p.mort - se, ymax = p.mort + se, width = 0.2), position = position_dodge(0.9)) + 
  theme_classic() +
  theme(legend.title = element_blank(), 
        axis.text.y = element_text(size = 10), 
        axis.text.x = element_text(size = 12, color = "black"), 
        axis.ticks.x = element_blank(),
        axis.title = element_text(size = 15), 
        axis.title.x = element_blank(),
        strip.text = element_text(size = 15),
        legend.text = element_text(size = 12),
        legend.position = "none",
        legend.key.size = unit(1.2, 'lines')) +
  labs(y = "Mortality", x = "Subplot") +
  ylim(0,1) +
  scale_fill_viridis_d()

ggsave(plot.mort.ST, filename = "Figures/mort-st.tiff", width = 6, height = 3, units = "in", dpi = 600)

#### Spp: Seed set ####

flo.seed.spp <- summarySE(flo.seed, measurevar = "n.seed.ind", groupvars = c("Subplot","strat", "Species"), na.rm = T)
flo.seed.spp$Subplot <- factor(flo.seed.spp$Subplot, levels = c("No Grass", "Grass", "Thatch"))
flo.seed.spp$Subplot <- revalue(flo.seed.spp$Subplot, c("Thatch" = "Litter addition", "No Grass" = "No grass"))

# Avoider
plot.seed.SA <- ggplot(flo.seed.spp[flo.seed.spp$strat == "SA",], aes(x = Species, y = n.seed.ind, fill = Subplot)) +
  geom_bar(position = "dodge", stat = "identity") +
  geom_errorbar(aes(ymin = n.seed.ind - se, ymax = n.seed.ind + se, width = 0.2), position = position_dodge(0.9)) + 
  theme_classic() +
  theme(legend.title = element_blank(), 
        axis.text.y = element_text(size = 10), 
        axis.text.x = element_text(size = 12, color = "black"),
        axis.ticks.x = element_blank(),
        axis.title = element_text(size = 15), 
        axis.title.x = element_blank(),
        strip.text = element_text(size = 15),
        legend.text = element_text(size = 12),
        legend.position = c(.84, .8),
        legend.key.size = unit(1.2, 'lines')) +
  labs(y = "Seeds per individual", x = "Subplot") +
  ylim(0,60) +
  scale_fill_viridis_d()

ggsave(plot.seed.SA, filename = "Figures/seed-sa.tiff", width = 6, height = 3, units = "in", dpi = 600)

# Tolerator
plot.seed.ST <- ggplot(flo.seed.spp[flo.seed.spp$strat == "ST",], aes(x = Species, y = n.seed.ind, fill = Subplot)) +
  geom_bar(position = "dodge", stat = "identity") +
  geom_errorbar(aes(ymin = n.seed.ind - se, ymax = n.seed.ind + se, width = 0.2), position = position_dodge(0.9)) + 
  theme_classic() +
  theme(legend.title = element_blank(), 
        axis.text.y = element_text(size = 10), 
        axis.text.x = element_text(size = 12, color = "black"), 
        axis.ticks.x = element_blank(),
        axis.title = element_text(size = 15), 
        axis.title.x = element_blank(),
        strip.text = element_text(size = 15),
        legend.text = element_text(size = 12),
        legend.position = "none",
        legend.key.size = unit(1.2, 'lines')) +
  labs(y = "Seeds per individual", x = "Subplot") +
  ylim(0,60) +
  scale_fill_viridis_d()

ggsave(plot.seed.ST, filename = "Figures/seed-st.tiff", width = 6, height = 3, units = "in", dpi = 600)

#### Elasticity Table ####
# full.elas <- gather(full[,c(2, 4:6, 23, 24)], key = "stage", value = "elas", -c(Plot, Species, Subplot, strat), na.rm = T)
# 
# ggplot(full.elas, aes(x = strat, y = elas, fill = stage)) +
#   geom_boxplot() +
#   facet_wrap(~Subplot)
# 
# hist((full.elas$elas + 1)^(1/2))
# 
# m.elas <- lmer(sqrt(elas + 1) ~ stage + (1|Plot:Species), data = full.elas[full.elas$strat == "SA" & full.elas$Subplot == "Thatch",]) # SA sig more elastic wrt fecundity... I think
# plot(fitted(m.elas), resid(m.elas))
# qqnorm(resid(m.elas))
# qqline(resid(m.elas), col = 2, lwd = 2, lty = 2)
# summary(m.elas)
# 
# full.elas.sum <- ddply(full, .(Subplot, strat), summarize, se.eg = se(elas.g, na.rm = T), se.eF = se(elas.F, na.rm = T), elas.g = mean(elas.g, na.rm = T), elas.F = mean(elas.F, na.rm = T))
# 
# full.elas.sum
# 
# ggplot(full.elas.sum, aes(x = Subplot, y = elas.F)) +
#   geom_point() +
#   geom_errorbar(aes(ymin = elas.F - se.eF, ymax = elas.F + se.eF), width = 0.02) +
#   facet_wrap(~strat)
# 
# ggplot(full.elas.sum, aes(x = Subplot, y = elas.g)) +
#   geom_point() +
#   geom_errorbar(aes(ymin = elas.g - se.eg, ymax = elas.g + se.eg), width = 0.02) +
#   facet_wrap(~strat)
# 
