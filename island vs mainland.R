library(MCMCglmm)

KP <- read.csv("submit.csv")
View(KP)

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


# THE END #
