#########LME4##########
###Script com efeito do Pop e Pop Interação###

##Variável NDF##
dados=read.table("dir")
dados
dados$Rep=factor(dados$Rep)
dados$Local=factor(dados$Local)
dados$Fix=factor(dados$Fix)
dados$Prog=factor(dados$Prog)
dados$Pop=factor(dados$Pop)
dados$Inter=factor(dados$Inter)
require(lme4)
mod=lmer(NDF~Rep+Fix+(1|Prog)+(1|Pop)+(1|Prog:Local)+(1|Pop:Local),data=dados,REML=TRUE)  
summary(mod)
ranef(mod)
REMLcrit(mod)
AIC(mod)
