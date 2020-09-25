### MTME FREQUENTISTA E BAYESIANO ###
setwd("C:/Users/Volpato/Desktop/MTME")
dadosFen = read.csv("C:\\Users\\Volpato\\Desktop\\MTME\\dados3.csv", header = T)

dadosFen$prog = factor(dadosFen$prog)
dadosFen$rep = factor(dadosFen$rep)
dadosFen$int = factor(dadosFen$int)
dadosFen$env = factor(dadosFen$env)

library(asreml)

###Mod Unicaracterístico###

modA <- asreml(fixed=DM ~ rep,
               random = ~ prog + int,
               data=dadosFen)

varComp = summary(modA)$varcomp
varComp

(l1 = modA$logl)

(K1 = length(modA$gammas))

LRT1 = -2* (summary(modA))$loglik
LRT1

##AIC##

(AIC = -2*l1 + 2*K1)
 

BLUPs=(summary(modA,all=TRUE)$coef.random)
BLUPs 

write.table(BLUPs,"BLUPsDM.txt",row.names=T,quote=F)

modA <- asreml(fixed=SW ~ rep,
               random = ~ prog + int,
               data=dadosFen)

varComp = summary(modA)$varcomp
varComp


(l1 = modA$logl)

(K1 = length(modA$gammas))

LRT1 = -2* (summary(modA))$loglik
LRT1

##AIC##

(AIC = -2*l1 + 2*K1)


BLUPs=(summary(modA,all=TRUE)$coef.random)
BLUPs 

write.table(BLUPs,"BLUPsSW.txt",row.names=T,quote=F)

modA <- asreml(fixed=SY ~ rep,
               random = ~ prog + int,
               data=dadosFen)

varComp = summary(modA)$varcomp
varComp


(l1 = modA$logl)

(K1 = length(modA$gammas))

LRT1 = -2* (summary(modA))$loglik
LRT1

##AIC##

(AIC = -2*l1 + 2*K1)

BLUPs=(summary(modA,all=TRUE)$coef.random)
BLUPs 

write.table(BLUPs,"BLUPsSY.txt",row.names=T,quote=F)

###Mod Mult-Trait Multi-Env - ASReml ###
modUS <- asreml(fixed=cbind(DM,SW,SY) ~ trait:rep - 1,
                random = ~ us(trait):prog + us(trait):int,
                rcov = ~ units:us(trait),
                maxiter=1000,
                data=dadosFen)

varComp = summary(modUS)$varcomp
varComp

(l1 = modUS$logl) ###logK

(K1 = length(modUS$gammas))  ###número de parâmetros

##AIC##

(AIC = -2*l1 + 2*K1)

##Obtenção dos BLUPS das progênies

BLUPs=(summary(modUS,all=TRUE)$coef.random)
BLUPs ###Para obter os valores de str error (SEP) de progênies para o cálculo da acurácia

write.table(BLUPs,"BLUPsMT.txt",row.names=T,quote=F)
