library(MCMCglmm)
library(lubridate)
library(ggplot2)
library(gridExtra)
df <- read.csv("dataset for submitting.csv", sep = ",")
KP <- df[,-c(5:6,8,15:16, 18, 21:24)]
View(df)
View(KP)
#--------
#converting to Julian date
date <- dmy(KP$date)
KP$j_date <- yday(date)
#------------------
#---- Makin dataframes specific for continent and island and deleting C. dubius

SMA.KP <- lm(log(weight)~log(wing), data=KP, na.action = na.omit)## Scale Mass Index, step 1 ##
KP$smi <- KP$weight*((mean(KP$wing, na.rm = T))/(KP$wing))^SMA.KP$coefficients[2]## Scale Mass Index, step 2 ##


#prevalences per population: 
#POOLED TOGETHER
prev <- table(KP$location, KP$pooled_0)
prev[1,2]*100/(prev[1,2]+prev[1,1])       # prev KP  Cadiz    28.85%
prev[2,2]*100/(prev[2,2]+prev[2,1])       # prev KP  Canary   11.11%
prev[3,2]*100/(prev[3,2]+prev[3,1])       # prev KP  Donana   23.08%
prev[4,2]*100/(prev[4,2]+prev[4,1])       # prev KP  Maio     19.32%


#CAMPY
prev1 <- table(KP$location, KP$campy_0)
prev1[1,2]*100/(prev1[1,2]+prev1[1,1])       # prev KP  Cadiz    1.92%
prev1[2,2]*100/(prev1[2,2]+prev1[2,1])       # prev KP  Canary   0%
prev1[3,2]*100/(prev1[3,2]+prev1[3,1])       # prev KP  Donana   1.92%
prev1[4,2]*100/(prev1[4,2]+prev1[4,1])       # prev KP  Maio     2.27%


#CHLAMY
prev2 <- table(KP$location, KP$chlamy_0)
prev2[1,2]*100/(prev2[1,2]+prev2[1,1])       # prev KP  Cadiz    7.69%
prev2[2,2]*100/(prev2[2,2]+prev2[2,1])       # prev KP  Canary   0%
prev2[3,2]*100/(prev2[3,2]+prev2[3,1])       # prev KP  Donana   1.92%
prev2[4,2]*100/(prev2[4,2]+prev2[4,1])       # prev KP  Maio     5.68%


#SALMO
prev3 <- table(KP$location, KP$salmo_0)
prev3[1,2]*100/(prev3[1,2]+prev3[1,1])       # prev KP  Cadiz    19.23%
prev3[2,2]*100/(prev3[2,2]+prev3[2,1])       # prev KP  Canary   11.11%
prev3[3,2]*100/(prev3[3,2]+prev3[3,1])       # prev KP  Donana   19.23%
prev3[4,2]*100/(prev3[4,2]+prev3[4,1])       # prev KP  Maio     11.36%




##### calculating Body weight means and SDs of males feamles island and mainland ######
mean(KP.cont$smi[KP$sex=="M"], na.rm = T)
sd(KP.cont$smi[KP$sex=="M"], na.rm = T)

mean(KP.cont$smi[KP$sex=="F"], na.rm = T)
sd(KP.cont$smi[KP$sex=="F"], na.rm = T)

mean(KP.isl$smi[KP.isl$sex=="M"], na.rm = T)
sd(KP.isl$smi[KP.isl$sex=="M"], na.rm = T)

mean(KP.isl$smi[KP.isl$sex=="F"], na.rm = T)
sd(KP.isl$smi[KP.isl$sex=="F"], na.rm = T)
#############################################################################




prior3 <- list(R = list(V = 1, nu = 0.002, n=1, fix=1),
               G = list(G1=list(V = 1*0.02, nu =7),
                        G1=list(V = 1*0.02, nu =7)))
#  --------
#  --------  MCMC  ----- CAMPY
mcmc1 <- MCMCglmm(campy_0 ~ lands_str*sex,
                  random = ~site + j_date,
                  data = KP,
                  prior = prior3,
                  family = "categorical", #binomial distribution in MCMCglmm
                  nitt=1000000, thin=600, burnin=1500)
summary(mcmc1)
plot(mcmc1)
mcmc1b <- MCMCglmm(campy_0 ~ lands_str*sex,
                  random = ~site + j_date,
                  data = KP,
                  prior = prior3,
                  family = "categorical", #binomial distribution in MCMCglmm
                  nitt=1000000, thin=600, burnin=1500)
summary(mcmc1b)
plot(mcmc1b)
diag(autocorr(mcmc1$VCV)[2, , ])
diag(autocorr(mcmc1b$VCV)[2, , ])
gelman.diag(mcmc.list(mcmc1$Sol, mcmc1b$Sol)) #diag of convergence for males
gelman.plot(mcmc.list(mcmc1$Sol, mcmc1b$Sol))


#  --------
#  --------  MCMC  ----- CHLAMY
mcmc2 <- MCMCglmm(chlamy_0 ~ lands_str*sex,
                  random = ~site + j_date,
                  data = KP,
                  prior = prior3,
                  family = "categorical", #binomial distribution in MCMCglmm
                  nitt=1000000, thin=600, burnin=1500)
summary(mcmc2)
plot(mcmc2)
mcmc2b <- MCMCglmm(chlamy_0 ~ lands_str*sex,
                   random = ~site + j_date,
                   data = KP,
                   prior = prior3,
                   family = "categorical", #binomial distribution in MCMCglmm
                   nitt=1000000, thin=600, burnin=1500)
summary(mcmc2b)
plot(mcmc2b)
diag(autocorr(mcmc2$VCV)[2, , ])
diag(autocorr(mcmc2b$VCV)[2, , ])
gelman.diag(mcmc.list(mcmc2$Sol, mcmc2b$Sol)) #diag of convergence for males
gelman.plot(mcmc.list(mcmc2$Sol, mcmc2b$Sol))


#  --------
#  --------  MCMC  ----- SALMO
mcmc3 <- MCMCglmm(salmo_0 ~ lands_str*sex,
                  random = ~site + j_date,
                  data = KP,
                  prior = prior3,
                  family = "categorical", #binomial distribution in MCMCglmm
                  nitt=1000000, thin=600, burnin=1500)
summary(mcmc3)
plot(mcmc3)
mcmc3b <- MCMCglmm(salmo_0 ~ lands_str*sex,
                   random = ~site + j_date,
                   data = KP,
                   prior = prior3,
                   family = "categorical", #binomial distribution in MCMCglmm
                   nitt=1000000, thin=600, burnin=1500)
summary(mcmc3b)
plot(mcmc3b)

diag(autocorr(mcmc3$VCV)[2, , ])
diag(autocorr(mcmc3b$VCV)[2, , ])
gelman.diag(mcmc.list(mcmc3$Sol, mcmc3b$Sol)) #diag of convergence for males
gelman.plot(mcmc.list(mcmc3$Sol, mcmc3b$Sol))


#  --------
#  --------  MCMC  ----- POOLED
mcmc4 <- MCMCglmm(pooled_0 ~ lands_str*sex,
                  random = ~site + j_date,
                  data = KP,
                  prior = prior3,
                  family = "categorical", #binomial distribution in MCMCglmm
                  nitt=1000000, thin=600, burnin=1500)
summary(mcmc4)
plot(mcmc4)
mcmc4b <- MCMCglmm(pooled_0 ~ lands_str*sex,
                   random = ~site + j_date,
                   data = KP,
                   prior = prior3,
                   family = "categorical", #binomial distribution in MCMCglmm
                   nitt=1000000, thin=600, burnin=1500)
summary(mcmc4b)
plot(mcmc4b)

diag(autocorr(mcmc4$VCV)[2, , ])
diag(autocorr(mcmc4b$VCV)[2, , ])
gelman.diag(mcmc.list(mcmc4$Sol, mcmc4b$Sol)) #diag of convergence for males
gelman.plot(mcmc.list(mcmc4$Sol, mcmc4b$Sol))







#  --------
#  --------  MCMC  ----- BODY CONDITION
prior2 <- list(R = list(V = 1, nu = 0.002),
               G = list(G1=list(V = 1*0.02, nu =7),
                        G1=list(V = 1*0.02, nu =7)))
mcmc5 <- MCMCglmm(smi ~ pooled_0*lands_str + pooled_0*sex + lands_str*sex,
                  random = ~site + j_date,
                  data = KP,
                  prior = prior2,
                  family = "gaussian", #binomial distribution in MCMCglmm
                  nitt=1000000, thin=600, burnin=1500)
summary(mcmc5)
plot(mcmc5)
mcmc5b <- MCMCglmm(smi ~ pooled_0*lands_str + pooled_0*sex + lands_str*sex,
                   random = ~site + j_date,
                   data = KP,
                   prior = prior2,
                   family = "gaussian", #binomial distribution in MCMCglmm
                   nitt=1000000, thin=600, burnin=1500)
summary(mcmc5b)
plot(mcmc5b)

diag(autocorr(mcmc5$VCV)[2, , ])
diag(autocorr(mcmc5b$VCV)[2, , ])
gelman.diag(mcmc.list(mcmc5$Sol, mcmc5b$Sol)) #diag of convergence for males
gelman.plot(mcmc.list(mcmc5$Sol, mcmc5b$Sol))





# PLOTS #########################################################

KP.salmo <- KP
KP.pooled <- KP

#### SEX in Salmonella infection
KP.salmo.f <- KP.salmo[KP.salmo$sex=="F",]
KP.salmo.m <- KP.salmo[KP.salmo$sex=="M",]
KP.salmo.f$salmo_0 <- KP.salmo.f$salmo_0*100/length((KP.salmo.f$salmo_0))
KP.salmo.m$salmo_0 <- KP.salmo.m$salmo_0*100/length((KP.salmo.m$salmo_0))
KP.salmo.sex <- rbind(KP.salmo.f, KP.salmo.m)

plot1 <- ggplot(data = KP.salmo.sex, aes(x = sex, y = salmo_0, fill=sex)) + 
  stat_summary(fun.y=sum,geom="bar",colour="black", width = .8)+
  scale_y_continuous(limits=c(0, 50))+
  theme_classic()+
  scale_fill_manual(values=c("darkorange1", "cyan3"))+
  labs(title = "Salmonella infection")+
  ylab("Percentage of birds infected")+
  theme(axis.text=element_text(size=14),
        axis.title.y = element_text(size=14),
        axis.title.x = element_text(size=14))+
  theme(legend.position="none")
plot1


#### SEX in Pooled infection
KP.pooled.f <- KP.pooled[KP.pooled$sex=="F",]
KP.pooled.m <- KP.pooled[KP.pooled$sex=="M",]
KP.pooled.f$pooled_0 <- KP.pooled.f$pooled_0*100/length((KP.pooled.f$pooled_0))
KP.pooled.m$pooled_0 <- KP.pooled.m$pooled_0*100/length((KP.pooled.m$pooled_0))
KP.pooled.sex <- rbind(KP.pooled.f, KP.pooled.m)


plot1.1 <- ggplot(data = KP.pooled.sex, aes(x = sex, y = pooled_0, fill=sex)) + 
  stat_summary(fun.y=sum,geom="bar",colour="black", width = .8)+
  scale_y_continuous(limits=c(0, 50))+
  theme_classic()+
  scale_fill_manual(values=c("darkorange1", "cyan3"))+
  labs(title ="Pooled bacteria infection") +
  ylab("Percentage of birds infected")+
  theme(axis.text=element_text(size=14),
        axis.title.y = element_text(size=14),
        axis.title.x = element_text(size=14))+
theme(legend.position="none")
plot1.1


grid.arrange(plot1, plot1.1,  nrow = 1)



#### BODY CONDITION and sex
plot2 <- ggplot(data = KP, aes(x = lands_str, y = smi, fill=sex)) +
  geom_boxplot()+
  theme_classic()+
  ylab("Body condition")+
  xlab(" ")+
  scale_fill_manual(values=c("darkorange1", "cyan3"))+
  theme(axis.text=element_text(size=14),
        axis.title.y = element_text(size=14),
        axis.title.x = element_text(size=14))+
  theme(legend.position="none")
plot2
  
plot(weight ~ sex, data = KP[KP$lands_str=="continent",])
t.test(weight ~ sex, data = KP[KP$lands_str=="continent",])
plot(smi ~ sex, data = KP[KP$lands_str=="continent",])
t.test(smi ~ sex, data = KP[KP$lands_str=="continent",])

### ISLAND vs CONTINENT #########################################################
# FULL

KP.plot <- KP
KP.plot$pooled_0[KP.plot$lands_str=="continent"] <- KP.plot$pooled_0[KP.plot$lands_str=="continent"]*100/length(KP.plot$pooled_0[KP.plot$lands_str=="continent"]) #FULL
KP.plot$campy_0[KP.plot$lands_str=="continent"] <- KP.plot$campy_0[KP.plot$lands_str=="continent"]*100/length(KP.plot$campy_0[KP.plot$lands_str=="continent"])    #CAMPY
KP.plot$chlamy_0[KP.plot$lands_str=="continent"] <- KP.plot$chlamy_0[KP.plot$lands_str=="continent"]*100/length(KP.plot$chlamy_0[KP.plot$lands_str=="continent"]) #CHLAMI
KP.plot$salmo_0[KP.plot$lands_str=="continent"] <- KP.plot$salmo_0[KP.plot$lands_str=="continent"]*100/length(KP.plot$salmo_0[KP.plot$lands_str=="continent"])    #SALMO

KP.plot$pooled_0[KP.plot$lands_str=="island"] <- KP.plot$pooled_0[KP.plot$lands_str=="island"]*100/length(KP.plot$pooled_0[KP.plot$lands_str=="island"]) #FULL
KP.plot$campy_0[KP.plot$lands_str=="island"] <- KP.plot$campy_0[KP.plot$lands_str=="island"]*100/length(KP.plot$campy_0[KP.plot$lands_str=="island"])    #CAMPY
KP.plot$chlamy_0[KP.plot$lands_str=="island"] <- KP.plot$chlamy_0[KP.plot$lands_str=="island"]*100/length(KP.plot$chlamy_0[KP.plot$lands_str=="island"]) #CHLAMI
KP.plot$salmo_0[KP.plot$lands_str=="island"] <- KP.plot$salmo_0[KP.plot$lands_str=="island"]*100/length(KP.plot$salmo_0[KP.plot$lands_str=="island"])    #SALMO


sum(KP.plot$pooled_0[KP.plot$lands_str=="continent"])

plot2 <- ggplot(data = KP.plot, aes(x = lands_str, y = pooled_0, fill=lands_str)) + 
  stat_summary(fun.y=sum,geom="bar",colour="black", width = .8)+
  scale_y_continuous(limits=c(0, 50))+
  theme_classic()+
  labs(title = "Pooled infeciton")+
  xlab(" ") + 
  ylab("Percentage of birds infected")+
  scale_fill_manual(values=c("goldenrod3", "skyblue4"))+
  theme(axis.text=element_text(size=14),
        axis.title.y = element_text(size=14),
        axis.title.x = element_text(size=14))+
  theme(legend.position="none")
plot2

# CAMPY

plot2.1 <- ggplot(data = KP.plot, aes(x = lands_str, y = campy_0, fill=lands_str)) + 
  stat_summary(fun.y=sum,geom="bar",colour="black", width = .8)+
  scale_y_continuous(limits=c(0, 50))+
  theme_classic()+
  labs(title = "Campylobacter")+
  xlab(" ") + 
  ylab("Percentage of birds infected")+
  scale_fill_manual(values=c("goldenrod3", "skyblue4"))+
  theme(axis.text=element_text(size=14),
        axis.title.y = element_text(size=14),
        axis.title.x = element_text(size=14))+
  theme(legend.position="none")
plot2.1

# CHLAMY

plot2.2 <- ggplot(data = KP.plot, aes(x = lands_str, y = chlamy_0, fill=lands_str)) + 
  stat_summary(fun.y=sum,geom="bar",colour="black", width = .8)+
  scale_y_continuous(limits=c(0, 50))+
  theme_classic()+
  labs(title = "Chlamydia")+
  xlab(" ") + 
  ylab("Percentage of birds infected")+
  scale_fill_manual(values=c("goldenrod3", "skyblue4"))+
  theme(axis.text=element_text(size=14),
        axis.title.y = element_text(size=14),
        axis.title.x = element_text(size=14))+
  theme(legend.position="none")
plot2.2

# SALMO

plot2.3 <- ggplot(data = KP.plot, aes(x = lands_str, y = salmo_0, fill=lands_str)) + 
  stat_summary(fun.y=sum,geom="bar",colour="black", width = .8)+
  scale_y_continuous(limits=c(0, 50))+
  theme_classic()+
  labs(title = "Salmonella")+
  xlab(" ") + 
  ylab("Percentage of birds infected")+
  scale_fill_manual(values=c("goldenrod3", "skyblue4"))+
  theme(axis.text=element_text(size=14),
        axis.title.y = element_text(size=14),
        axis.title.x = element_text(size=14))+
  theme(legend.position="none")
plot2.3

grid.arrange(plot2.1, plot2.2, plot2.3,plot2,  nrow = 1)

