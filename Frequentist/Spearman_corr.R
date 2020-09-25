###################****************************Hybrid*******************************##############

setwd("C:/Users/Volpato/OneDrive/Sorgo MTME/Analises_prior/SpearmanRankCorr")
PBV<-read.csv("Hybrid_PBV.csv",header=T)
head(PBV)

corr_flor <- cor.test(x=PBV$Flor_Bayes, y=PBV$Flor_Freq, method = 's')
corr_flor

corr_AP <- cor.test(x=PBV$AP_Bayes, y=PBV$AP_Freq, method = 's')
corr_AP

corr_prod <- cor.test(x=PBV$Prod_Bayes, y=PBV$Prod_Freq, method = 's')
corr_prod

library(boot)

#Flor
foo <- function(data, indices, cor.type){
  dt<-data[indices,]
  cor(dt[,2], dt[,3], method=cor.type)
}

set.seed(777)
myBootstrap <- boot(PBV, foo, R=1000, cor.type='s')

head(myBootstrap$t)
myBootstrap$t0
colMeans(myBootstrap$t)-myBootstrap$t0
apply(myBootstrap$t,2,sd)
plot(myBootstrap, index=1)

boot.ci(myBootstrap, conf = c(0.90, 0.95),
        type = c("norm", "basic", "perc", "bca"))

#The 95% CI for the normal bootstrap is obtained by calculating:
with(myBootstrap, 2*t0 - mean(t) + qnorm(c(0.025, 0.975)) %o% sqrt(var(t)[1,1]))

#The p-value is thus obtained:
#H0: p = 0 
#Ha: p different of 0
with(myBootstrap, pnorm(abs((2*t0 - mean(t) - 1) / sqrt(var(t)[1,1])), lower.tail=F)*2)


#AP
foo <- function(data, indices, cor.type){
  dt<-data[indices,]
  cor(dt[,4], dt[,5], method=cor.type)
}

set.seed(777)
myBootstrap <- boot(PBV, foo, R=1000, cor.type='s')

head(myBootstrap$t)
myBootstrap$t0
colMeans(myBootstrap$t)-myBootstrap$t0
apply(myBootstrap$t,2,sd)
plot(myBootstrap, index=1)

boot.ci(myBootstrap, conf = c(0.90, 0.95),
        type = c("norm", "basic", "perc", "bca"))

#The 95% CI for the normal bootstrap is obtained by calculating:
with(myBootstrap, 2*t0 - mean(t) + qnorm(c(0.025, 0.975)) %o% sqrt(var(t)[1,1]))

#The p-value is thus obtained:
#H0: p = 0 
#Ha: p different of 0
with(myBootstrap, pnorm(abs((2*t0 - mean(t) - 1) / sqrt(var(t)[1,1])), lower.tail=F)*2)


#Prod
foo <- function(data, indices, cor.type){
  dt<-data[indices,]
  cor(dt[,6], dt[,7], method=cor.type)
}

set.seed(777)
myBootstrap <- boot(PBV, foo, R=1000, cor.type='s')

head(myBootstrap$t)
myBootstrap$t0
colMeans(myBootstrap$t)-myBootstrap$t0
apply(myBootstrap$t,2,sd)
plot(myBootstrap, index=1)

boot.ci(myBootstrap, conf = c(0.90, 0.95),
        type = c("norm", "basic", "perc", "bca"))

#The 95% CI for the normal bootstrap is obtained by calculating:
with(myBootstrap, 2*t0 - mean(t) + qnorm(c(0.025, 0.975)) %o% sqrt(var(t)[1,1]))

#The p-value is thus obtained:
#H0: p = 0 
#Ha: p different of 0
with(myBootstrap, pnorm(abs((2*t0 - mean(t) - 1) / sqrt(var(t)[1,1])), lower.tail=F)*2)

###################****************************Restored line*******************************##############

setwd("C:/Users/Volpato/OneDrive/Sorgo MTME/Analises_prior/SpearmanRankCorr")
PBV<-read.csv("Pai_PBV.csv",header=T)
head(PBV)

corr_flor <- cor.test(x=PBV$FL_Bayes, y=PBV$FL_Freq, method = 's')
corr_flor

corr_AP <- cor.test(x=PBV$PH_Bayes, y=PBV$PH_Freq, method = 's')
corr_AP

corr_prod <- cor.test(x=PBV$GY_Bayes, y=PBV$GY_Freq, method = 's')
corr_prod

library(boot)

#Flor
foo <- function(data, indices, cor.type){
  dt<-data[indices,]
  cor(dt[,2], dt[,3], method=cor.type)
}

set.seed(777)
myBootstrap <- boot(PBV, foo, R=1000, cor.type='s')

head(myBootstrap$t)
myBootstrap$t0
colMeans(myBootstrap$t)-myBootstrap$t0
apply(myBootstrap$t,2,sd)
plot(myBootstrap, index=1)

boot.ci(myBootstrap, conf = c(0.90, 0.95),
        type = c("norm", "basic", "perc", "bca"))

#The 95% CI for the normal bootstrap is obtained by calculating:
with(myBootstrap, 2*t0 - mean(t) + qnorm(c(0.025, 0.975)) %o% sqrt(var(t)[1,1]))

#The p-value is thus obtained:
#H0: p = 0 
#Ha: p different of 0
with(myBootstrap, pnorm(abs((2*t0 - mean(t) - 1) / sqrt(var(t)[1,1])), lower.tail=F)*2)


#AP
foo <- function(data, indices, cor.type){
  dt<-data[indices,]
  cor(dt[,4], dt[,5], method=cor.type)
}

set.seed(777)
myBootstrap <- boot(PBV, foo, R=1000, cor.type='s')

head(myBootstrap$t)
myBootstrap$t0
colMeans(myBootstrap$t)-myBootstrap$t0
apply(myBootstrap$t,2,sd)
plot(myBootstrap, index=1)

boot.ci(myBootstrap, conf = c(0.90, 0.95),
        type = c("norm", "basic", "perc", "bca"))

#The 95% CI for the normal bootstrap is obtained by calculating:
with(myBootstrap, 2*t0 - mean(t) + qnorm(c(0.025, 0.975)) %o% sqrt(var(t)[1,1]))

#The p-value is thus obtained:
#H0: p = 0 
#Ha: p different of 0
with(myBootstrap, pnorm(abs((2*t0 - mean(t) - 1) / sqrt(var(t)[1,1])), lower.tail=F)*2)


#Prod
foo <- function(data, indices, cor.type){
  dt<-data[indices,]
  cor(dt[,6], dt[,7], method=cor.type)
}

set.seed(777)
myBootstrap <- boot(PBV, foo, R=1000, cor.type='s')

head(myBootstrap$t)
myBootstrap$t0
colMeans(myBootstrap$t)-myBootstrap$t0
apply(myBootstrap$t,2,sd)
plot(myBootstrap, index=1)

boot.ci(myBootstrap, conf = c(0.90, 0.95),
        type = c("norm", "basic", "perc", "bca"))

#The 95% CI for the normal bootstrap is obtained by calculating:
with(myBootstrap, 2*t0 - mean(t) + qnorm(c(0.025, 0.975)) %o% sqrt(var(t)[1,1]))

#The p-value is thus obtained:
#H0: p = 0 
#Ha: p different of 0
with(myBootstrap, pnorm(abs((2*t0 - mean(t) - 1) / sqrt(var(t)[1,1])), lower.tail=F)*2)


####################################

###################****************************Male sterile*******************************##############

setwd("C:/Users/Volpato/OneDrive/Sorgo MTME/Analises_prior/SpearmanRankCorr")
PBV<-read.csv("Mae_PBV.csv",header=T)
head(PBV)

corr_flor <- cor.test(x=PBV$FL_Bayes, y=PBV$FL_Freq, method = 's')
corr_flor

corr_AP <- cor.test(x=PBV$PH_Bayes, y=PBV$PH_Freq, method = 's')
corr_AP

corr_prod <- cor.test(x=PBV$GY_Bayes, y=PBV$GY_Freq, method = 's')
corr_prod

library(boot)

#Flor
foo <- function(data, indices, cor.type){
  dt<-data[indices,]
  cor(dt[,2], dt[,3], method=cor.type)
}

set.seed(777)
myBootstrap <- boot(PBV, foo, R=1000, cor.type='s')

head(myBootstrap$t)
myBootstrap$t0
colMeans(myBootstrap$t)-myBootstrap$t0
apply(myBootstrap$t,2,sd)
plot(myBootstrap, index=1)

boot.ci(myBootstrap, conf = c(0.90, 0.95),
        type = c("norm", "basic", "perc", "bca"))

#The 95% CI for the normal bootstrap is obtained by calculating:
with(myBootstrap, 2*t0 - mean(t) + qnorm(c(0.025, 0.975)) %o% sqrt(var(t)[1,1]))

#The p-value is thus obtained:
#H0: p = 0 
#Ha: p different of 0
with(myBootstrap, pnorm(abs((2*t0 - mean(t) - 1) / sqrt(var(t)[1,1])), lower.tail=F)*2)


#AP
foo <- function(data, indices, cor.type){
  dt<-data[indices,]
  cor(dt[,4], dt[,5], method=cor.type)
}

set.seed(777)
myBootstrap <- boot(PBV, foo, R=1000, cor.type='s')

head(myBootstrap$t)
myBootstrap$t0
colMeans(myBootstrap$t)-myBootstrap$t0
apply(myBootstrap$t,2,sd)
plot(myBootstrap, index=1)

boot.ci(myBootstrap, conf = c(0.90, 0.95),
        type = c("norm", "basic", "perc", "bca"))

#The 95% CI for the normal bootstrap is obtained by calculating:
with(myBootstrap, 2*t0 - mean(t) + qnorm(c(0.025, 0.975)) %o% sqrt(var(t)[1,1]))

#The p-value is thus obtained:
#H0: p = 0 
#Ha: p different of 0
with(myBootstrap, pnorm(abs((2*t0 - mean(t) - 1) / sqrt(var(t)[1,1])), lower.tail=F)*2)


#Prod
foo <- function(data, indices, cor.type){
  dt<-data[indices,]
  cor(dt[,6], dt[,7], method=cor.type)
}

set.seed(777)
myBootstrap <- boot(PBV, foo, R=1000, cor.type='s')

head(myBootstrap$t)
myBootstrap$t0
colMeans(myBootstrap$t)-myBootstrap$t0
apply(myBootstrap$t,2,sd)
plot(myBootstrap, index=1)

boot.ci(myBootstrap, conf = c(0.90, 0.95),
        type = c("norm", "basic", "perc", "bca"))

#The 95% CI for the normal bootstrap is obtained by calculating:
with(myBootstrap, 2*t0 - mean(t) + qnorm(c(0.025, 0.975)) %o% sqrt(var(t)[1,1]))

#The p-value is thus obtained:
#H0: p = 0 
#Ha: p different of 0
with(myBootstrap, pnorm(abs((2*t0 - mean(t) - 1) / sqrt(var(t)[1,1])), lower.tail=F)*2)
