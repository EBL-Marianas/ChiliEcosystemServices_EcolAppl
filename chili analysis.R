### Define functions

ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}


### Load packages
packages <- c("plyr", "lme4", "coxme", "multcomp", "lsmeans", "pscl")
ipak(pkg = packages)


### Load data
load("chili.RData")


### Germination analysis
# Cox Regression
germination$treatment.sp <- relevel(germination$treatment.sp, ref="mech")
cr.germ.mod <- coxme(Surv(days.since.sown, germ01) ~ treatment.sp + (1|date.collected), data = germination)
summary(cr.germ.mod)
summary(glht(cr.germ.mod, linfct = mcp(treatment.sp = "Tukey")))

# GLMM
germination.short$treatment.sp <- relevel(germination.short$treatment.sp, ref="mech")
gl.germ.mod <- glmer(cbind(germ,right.censored) ~ treatment.sp + (1|date.collected),
                data = germination.short,
                family = "binomial")
summary(gl.germ.mod)
summary(glht(gl.germ.mod, linfct = mcp(treatment.sp = "Tukey")))


### Removal analysis
# Cox Regression
cr.mod <- coxme(Surv(day, death) ~ condition * location  + (1|pileID), data = removal)
summary(cr.mod)
# Make this digestible for glht
cr.mod <- coxme(Surv(day, death) ~ cl + (1|pileID), data = removal)
summary(cr.mod)
summary(glht(cr.mod, rbind("GP.Far vs WF.Near" =  c(0, 0, 0, 0, -1))))
summary(glht(cr.mod, rbind("GP.Far vs GP.Near" =  c(0, 0, 1, 0, 0))))
removal$cl <- relevel(removal$cl, ref = "GP.Near") # relevel to make last comparison
cr.mod <- coxme(Surv(day, death) ~ cl + (1|pileID), data = removal)
summary(glht(cr.mod, rbind("GP.Near vs WF.Near" =  c(0, 0, 0, 0, -1))))

# GLMM
gl.mod <- glmer(cbind(removed,not.removed) ~ condition * location + (1|pileID),
                data = removal.short,
                family = "binomial")
summary(gl.mod)
# Make this digestible for glht
gl.mod <- glmer(cbind(removed,not.removed) ~ cl - 1 + (1|pileID),
                data = removal.short,
                family = "binomial")
summary(glht(gl.mod, rbind("GP.Far vs WF.Near" =  c(1, 0, 0, 0, 0, -1))))
summary(glht(gl.mod, rbind("GP.Near vs WF.Near" =  c(0, 0, 0, 1, 0, -1))))
summary(glht(gl.mod, rbind("GP.Far vs GP.Near" =  c(1, 0, 0, -1, 0, 0))))


### Abundance analysis
glm.nbin <- glm.nb(count ~ island, data = abundance)
summary(glht(glm.nbin, linfct = mcp(island = "Tukey")))
