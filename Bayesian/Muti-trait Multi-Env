#############################modBAY########BAYESIANO###

setwd("C:/Users/Volpato/Desktop/MTME")
dadosFen = read.csv("C:\\Users\\Volpato\\Desktop\\MTME\\dados3.csv", header = T)


dadosFen$prog = factor(dadosFen$prog)
dadosFen$rep = factor(dadosFen$rep)
dadosFen$int = factor(dadosFen$int)
dadosFen$env = factor(dadosFen$env)

library(MCMCglmm)

Niter=1000000  #1000000
Burnin=500000  #500000
Thin=5
nefit=(Niter-Burnin)/Thin

### Modelo Unicaracteristico - BST ###

prioriST <- list(R = list(V=1, nu=0.002), 
              G = list(G1 = list(V=1, nu=0.002),
                       G2 = list(V=1, nu=0.002)))
prioriST

modBU = MCMCglmm(fixed = DM ~ rep,
                  random = ~ prog + int, data=dadosFen,
                  prior=prioriST, family = "gaussian", DIC=TRUE, verbose=FALSE, nitt=Niter, thin=Thin, burnin=Burnin,pr=TRUE,singular.ok=TRUE)

summary(modBU)

heidel.diag(modBU$VCV)

modBU2 = MCMCglmm(fixed = SW ~ rep,
                 random = ~ prog + int, data=dadosFen,
                 prior=prioriST, family = "gaussian", DIC=TRUE, verbose=FALSE, nitt=Niter, thin=Thin, burnin=Burnin,pr=TRUE,singular.ok=TRUE)

summary(modBU2)

heidel.diag(modBU2$VCV)

modBU3 = MCMCglmm(fixed = SY ~ rep,
                  random = ~ prog + int, data=dadosFen,
                  prior=prioriST, family = "gaussian", DIC=TRUE, verbose=FALSE, nitt=Niter, thin=Thin, burnin=Burnin,pr=TRUE,singular.ok=TRUE)

summary(modBU3)

heidel.diag(modBU3$VCV)

### a1 e bi shoud be modificated ###
a1=modBU3$Sol
b1=modBU3$VCV
a1
b1

#herdabilidade do MCMC para DM UNICACTERISTICO
heritDMu <- b1[,"prog"]/(b1[, "prog"] + b1[, "int"]/2 + b1[, "units"]/6) #herdabilidade do MCMC
effectiveSize(heritDMu)
mean(heritDMu)
median(heritDMu)
HPDinterval(heritDMu) #Display 95% credible interval

#herdabilidade do MCMC para SW UNICACTERISTICO
heritSWu <- b1[,"prog"]/(b1[, "prog"] + b1[, "int"]/2 + b1[, "units"]/6) #herdabilidade do MCMC
effectiveSize(heritSWu)
mean(heritSWu)
median(heritSWu)
HPDinterval(heritSWu) #Display 95% credible interval

#herdabilidade do MCMC para SY UNICACTERISTICO
heritSYu <- b1[,"prog"]/(b1[, "prog"] + b1[, "int"]/2 + b1[, "units"]/6) #herdabilidade do MCMC
effectiveSize(heritSYu)
mean(heritSYu)
median(heritSYu)
HPDinterval(heritSYu) #Display 95% credible interval

##Breeding values BST - DM##

a1=modBU$Sol
b1=modBU$VCV
a1
b1

dim(b1)
Prog_est= modBU$Sol[,7:209]  #valores genéticos de cada prog
Prog_est

trait_VG=modBU$VCV[,1]
trait_VG

trait_INT=modBU$VCV[,2]
trait_INT

trait_RES=modBU$VCV[,3]
trait_RES

DIC=rep(modBU$DIC,nefit)
DIC

result_modBU = cbind(Prog_est,trait_VG,trait_INT,trait_RES,DIC)
result_modBU

write.table(result_modBU,"result_modBAYSTdm.txt",row.names=F,quote=F)

library(boa)
saida_modBU=boa.stats(result_modBU,c(0.025,0.975),100000)
write.table(saida_modBU,"saida_modBAYSTdm.txt",row.names=T,quote=F)

##Breeding values BST - SW##

a1=modBU2$Sol
b1=modBU2$VCV
a1
b1

dim(b1)
Prog_est= modBU2$Sol[,7:209]  #valores genéticos de cada prog
Prog_est

trait_VG=modBU2$VCV[,1]
trait_VG

trait_INT=modBU2$VCV[,2]
trait_INT

trait_RES=modBU2$VCV[,3]
trait_RES

DIC=rep(modBU2$DIC,nefit)
DIC

result_modBU2 = cbind(Prog_est,trait_VG,trait_INT,trait_RES,DIC)
result_modBU2

write.table(result_modBU2,"result_modBAYSTsw.txt",row.names=F,quote=F)

library(boa)
saida_modBU2=boa.stats(result_modBU2,c(0.025,0.975),100000)
write.table(saida_modBU2,"saida_modBAYSTsw.txt",row.names=T,quote=F)

##Breeding values BST - SY##

a1=modBU3$Sol
b1=modBU3$VCV
a1
b1

dim(b1)
Prog_est= modBU3$Sol[,7:209]  #valores genéticos de cada prog
Prog_est

trait_VG=modBU3$VCV[,1]
trait_VG

trait_INT=modBU3$VCV[,2]
trait_INT

trait_RES=modBU3$VCV[,3]
trait_RES

DIC=rep(modBU3$DIC,nefit)
DIC

result_modBU3 = cbind(Prog_est,trait_VG,trait_INT,trait_RES,DIC)
result_modBU3

write.table(result_modBU3,"result_modBAYSTsy.txt",row.names=F,quote=F)

library(boa)
saida_modBU3=boa.stats(result_modBU3,c(0.025,0.975),100000)
write.table(saida_modBU3,"saida_modBAYSTsy.txt",row.names=T,quote=F)

### TESTE DE SIGNIFICANCIA DIC DO Modelo Unicaracteristico - BST ###

prioriST1 <- list(R = list(V=1, nu=0.002), 
                 G = list(G1 = list(V=1, nu=0.002)))
prioriST1

modBUtest = MCMCglmm(fixed = SY ~ rep,
                 random = ~  int, data=dadosFen,
                 prior=prioriST1, family = "gaussian", DIC=TRUE, verbose=FALSE, nitt=Niter, thin=Thin, burnin=Burnin,pr=TRUE,singular.ok=TRUE)

summary(modBUtest)

modBUtest2 = MCMCglmm(fixed = SY ~ rep,
                  random = ~ prog, data=dadosFen,
                  prior=prioriST1, family = "gaussian", DIC=TRUE, verbose=FALSE, nitt=Niter, thin=Thin, burnin=Burnin,pr=TRUE,singular.ok=TRUE)

summary(modBUtest2)


### Modelo BAYES MTME ######

priori=list(R=list(V=diag(3)/2, nu=2), 
G=list(G1=list(V=diag(3)/2, nu=2), 
 G2=list(V=diag(3)/2, nu=2)))

priori

modBAY = MCMCglmm(fixed = cbind(DM,SW,SY) ~ trait:rep - 1,
                  random = ~ us(trait):prog + us(trait):int,
                  rcov = ~ us(trait):units, data=dadosFen,
                  prior=priori, family = c("gaussian", "gaussian","gaussian"), DIC=TRUE, verbose=FALSE, nitt=Niter, thin=Thin, burnin=Burnin,pr=TRUE,singular.ok=TRUE)

summary(modBAY)

###Verificar se o número de interações foi suficiente ###
#plot(modBAY$Sol) ##Distribuições das interações fixa (Não usual, pois fazendo um numero alto de interações não há problemas)
#plot(modBAY$VCV) ##Distribuições das interações aleatórias

###Teste de significancia DIC ###

priori1 =list(R=list(V=diag(3)/2, nu=2), 
            G=list(G1=list(V=diag(3)/2, nu=2)))
priori1


modBAY1 = MCMCglmm(fixed = cbind(DM,SW,SY) ~ trait:rep - 1,
                  random = ~ us(trait):int,
                  rcov = ~ us(trait):units, data=dadosFen,
                  prior=priori1, family = c("gaussian", "gaussian","gaussian"), DIC=TRUE, verbose=FALSE, nitt=Niter, thin=Thin, burnin=Burnin,pr=TRUE,singular.ok=TRUE)

summary(modBAY1)


modBAY2 = MCMCglmm(fixed = cbind(DM,SW,SY) ~ trait:rep - 1,
                   random = ~ us(trait):prog,
                   rcov = ~ us(trait):units, data=dadosFen,
                   prior=priori1, family = c("gaussian", "gaussian","gaussian"), DIC=TRUE, verbose=FALSE, nitt=Niter, thin=Thin, burnin=Burnin,pr=TRUE,singular.ok=TRUE)

summary(modBAY2)


####################

heidel.diag(modBAY$VCV) ##diagnostic tests of convergence (The p-values must exceed 0:05)

a1=modBAY$Sol
b1=modBAY$VCV
a1
b1

#herdabilidade do MCMC para DM
heritDM <- b1[,"traitDM:traitDM.prog"]/(b1[, "traitDM:traitDM.prog"] + b1[, "traitDM:traitDM.int"]/2 + b1[, "traitDM:traitDM.units"]/6) #herdabilidade do MCMC
effectiveSize(heritDM)
mean(heritDM)
median(heritDM)
HPDinterval(heritDM) #Display 95% credible interval
plot(heritDM,main="VarDM")

#herdabilidade do MCMC para SW
heritSW <- b1[,"traitSW:traitSW.prog"]/(b1[, "traitSW:traitSW.prog"] + b1[, "traitSW:traitSW.int"]/2 + b1[, "traitSW:traitSW.units"]/6) #herdabilidade do MCMC
effectiveSize(heritSW)
mean(heritSW)
median(heritSW)
HPDinterval(heritSW) #Display 95% credible interval
plot(heritSW,main="VarSW")

#herdabilidade do MCMC para SY
heritSY <- b1[,"traitSY:traitSY.prog"]/(b1[, "traitSY:traitSY.prog"] + b1[, "traitSY:traitSY.int"]/2 + b1[, "traitSY:traitSY.units"]/6) #herdabilidade do MCMC
effectiveSize(heritSY)
mean(heritSY)
median(heritSY)
HPDinterval(heritSY) #Display 95% credible interval
plot(heritSY,main="VarSY")


### Plot herdability ST and MTME Bayes ###
##Plot somente dos ST##
#install.packages("bayesplot")
library(bayesplot)

###Plot das H2prog BST models ###
btest1<-matrix(c(heritDMu),,1,byrow=T)
btest2<-matrix(c(heritSWu),,1,byrow=T)
btest3<-matrix(c(heritSYu),,1,byrow=T)
BST<- cbind(btest1[,1], btest2[,1], btest3[,1])
BST
colnames(BST)<-c("DM","SW","SY")
BST
color_scheme_set("darkgray")
mcmc_areas(BST,pars = c("DM","SW","SY"), prob = 0.95,point_est = "mean")  +
  ggplot2::labs(
    x="Heritability",
    title = "Posterior distributions - DM, SW and SY",
    subtitle = "with mean and 95% intervals")

##Plot dos componentes BST ##
b1
color_scheme_set("darkgray")
mcmc_areas(b1,pars = c("prog","int","units"), prob = 0.95,point_est = "mean")  +
  ggplot2::labs(
    x="Estimates of variance components",
    title = "Posterior distributions - SY",
    subtitle = "with mean and 95% intervals")

##Plot dos componentes BMTME ##
#install.packages("gridExtra")
b1
PlotBMTME1<- cbind(b1[,1], b1[,5],b1[,9])
colnames(PlotBMTME1)<-c("DM","SW","SY")
color_scheme_set("darkgray")
mcmc_pairs(PlotBMTME1,pars = c("DM","SW","SY"), off_diag_args = list(size=1.5)) 

## Plot H2 prog BMTME ##
btest1<-matrix(c(heritDM),,1,byrow=T)
btest2<-matrix(c(heritSW),,1,byrow=T)
btest3<-matrix(c(heritSY),,1,byrow=T)
PlotDM<- cbind(btest1[,1], btest2[,1],btest3[,1])
colnames(PlotDM)<-c("DM","SW","SY")
PlotDM
color_scheme_set("darkgray")
mcmc_pairs(PlotDM,pars = c("DM","SW","SY"), off_diag_args = list(size=1.5)) 

## Plot ST + MTME ##

btest1<-matrix(c(heritDMu),,1,byrow=T)
btest2<-matrix(c(heritDM),,1,byrow=T)
PlotDM<- cbind(btest1[,1], btest2[,1])
colnames(PlotDM)<-c("BST","BMTME")
PlotDM

color_scheme_set("darkgray")
mcmc_areas(PlotDM,pars = c("BST","BMTME"), prob = 0.95,point_est = "mean")  +
  ggplot2::labs(
    x="Heritability",
    title = "Posterior distributions - DM",
    subtitle = "with mean and 95% intervals")


##genetic correlation:
corr.genSWDM<-modBAY$VCV[,"traitSW:traitDM.prog"]/
  sqrt(modBAY$VCV[,"traitSW:traitSW.prog"]*modBAY$VCV[,"traitDM:traitDM.prog"]) 
mean(corr.genSWDM)

corr.genSYDM<-modBAY$VCV[,"traitSY:traitDM.prog"]/
  sqrt(modBAY$VCV[,"traitSY:traitSY.prog"]*modBAY$VCV[,"traitDM:traitDM.prog"]) 
mean(corr.genSYDM)

corr.genSYSW<-modBAY$VCV[,"traitSY:traitSW.prog"]/
  sqrt(modBAY$VCV[,"traitSY:traitSY.prog"]*modBAY$VCV[,"traitSW:traitSW.prog"]) 
mean(corr.genSYSW)

##Breeding values BMTME ##

dim(b1)
Prog_est= modBAY$Sol[,19:627]  #valores genéticos de cada prog
Prog_est

trait_VG_DMDM = modBAY$VCV[,1] #var gen para DM
trait_COV_SWDM = modBAY$VCV[,2] #COV gen para SW:DM
trait_COV_SYDM = modBAY$VCV[,3] #COV gen para SY:DM
trait_VG_SWSW = modBAY$VCV[,5] #var gen para SW
trait_COV_SYSW = modBAY$VCV[,6] #COV gen para SY:SW
trait_VG_SYSY = modBAY$VCV[,9] #var gen para SY

trait_INT_DMDM = modBAY$VCV[,10] #var int para DM
trait_INT_SWDM = modBAY$VCV[,11] #COV int para SW:DM
trait_INT_SYDM = modBAY$VCV[,12] #COV int para SY:DM
trait_INT_SWSW = modBAY$VCV[,14] #var int para SW
trait_INT_SYSW = modBAY$VCV[,15] #COV int para SY:SW
trait_INT_SYSY = modBAY$VCV[,18] #var int para SY

trait_RES_DMDM = modBAY$VCV[,19] #var res para DM
trait_RES_SWDM = modBAY$VCV[,20] #COV res para SW:DM
trait_RES_SYDM = modBAY$VCV[,21] #COV res para SY:DM
trait_RES_SWSW = modBAY$VCV[,23] #var res para SW
trait_RES_SYSW = modBAY$VCV[,24] #COV res para SY:SW
trait_RES_SYSY = modBAY$VCV[,27] #var res para SY

DIC=rep(modBAY$DIC,nefit)
DIC

result_modBAY=cbind(Prog_est,trait_VG_DMDM,trait_COV_SWDM,trait_COV_SYDM,trait_VG_SWSW,trait_COV_SYSW,trait_VG_SYSY,trait_INT_DMDM,trait_INT_SWDM,trait_INT_SYDM,trait_INT_SWSW,trait_INT_SYSW,trait_INT_SYSY,trait_RES_DMDM,trait_RES_SWDM,trait_RES_SYDM,trait_RES_SWSW,trait_RES_SYSW,trait_RES_SYSY, DIC)

result_modBAY

write.table(result_modBAY,"result_modBAYMTME.txt",row.names=F,quote=F)

library(boa)
saida_modBAY=boa.stats(result_modBAY,c(0.025,0.975),100000)
write.table(saida_modBAY,"saida_modBAYMTME.txt",row.names=T,quote=F)


## Correlações do BST ##

dadosFencor = read.csv("C:\\Users\\Volpato\\Desktop\\MTME\\dadostestcor.csv", header = T)

dadosFencor$prog = factor(dadosFencor$prog)
dadosFencor$rep = factor(dadosFencor$rep)
dadosFencor$int = factor(dadosFencor$int)
dadosFencor$env = factor(dadosFencor$env)

library(MCMCglmm)

Niter=1000000  #1000000
Burnin=500000  #500000
Thin=5
nefit=(Niter-Burnin)/Thin

### Modelo Unicaracteristico - BST ###

prioriST <- list(R = list(V=1, nu=0.002), 
                 G = list(G1 = list(V=1, nu=0.002),
                          G2 = list(V=1, nu=0.002)))
prioriST

modBUc = MCMCglmm(fixed = DMSW ~ rep,
                 random = ~ prog + int, data=dadosFencor,
                 prior=prioriST, family = "gaussian", DIC=TRUE, verbose=FALSE, nitt=Niter, thin=Thin, burnin=Burnin,pr=TRUE,singular.ok=TRUE)

summary(modBUc)

heidel.diag(modBUc$VCV)

modBU2c = MCMCglmm(fixed = DMSY ~ rep,
                  random = ~ prog + int, data=dadosFencor,
                  prior=prioriST, family = "gaussian", DIC=TRUE, verbose=FALSE, nitt=Niter, thin=Thin, burnin=Burnin,pr=TRUE,singular.ok=TRUE)

summary(modBU2c)

heidel.diag(modBU2c$VCV)

modBU3c = MCMCglmm(fixed = SWSY ~ rep,
                  random = ~ prog + int, data=dadosFencor,
                  prior=prioriST, family = "gaussian", DIC=TRUE, verbose=FALSE, nitt=Niter, thin=Thin, burnin=Burnin,pr=TRUE,singular.ok=TRUE)

summary(modBU3c)

heidel.diag(modBU3c$VCV)
