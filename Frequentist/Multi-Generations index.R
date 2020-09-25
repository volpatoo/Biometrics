####### R O T I N A S   M U L T I G E R A Ç O E S #######
##Objetivo: substituir a matrix de corr (A) do modelo 180 do Selegen
##obtendo as correlações através dos conjunto de dados. 
##Para tal, temos diferentes estruturas de VCOV a serem testadas. 

### Carregandp os dados Fenotipicos

setwd("")
dadosFen2 = read.csv(")

### Todos os efeitos do modelo devem ser fatores

dadosFen2$PLOT = factor(dadosFen2$PLOT)
dadosFen2$GENO = factor(dadosFen2$GENO)
dadosFen2$REP = factor(dadosFen2$REP)
dadosFen2$SITE = factor(dadosFen2$SITE)
dadosFen2$CHECK = factor(dadosFen2$CHECK)


### Matrix de efeitos da interação
## SÃo USADOS P/ O INDICE MULTIGERAÇÕES

## Matrix A (paper): VCOV de cada genotipo em diferentes geracoes

corr = c(1,1,1,1,1.5,1.5,1,1.5,1.75) # De acordo com o paper do Deon
matA = matrix(corr, nrow = 3, ncol = 3, byrow = T)
matA

## Matrix G: assumindo covariancia 0 entre familias: é uma identidade
matG = diag(nlevels(dadosFen2$GENO))
matG
dim(matG)

## Matrix da interação
matINT = matA %x% matG
dim(matINT)
matINT

## Inversa
matINV = solve(matINT)
rownames(matINV) = levels(dadosFen2$INTER)
colnames(matINV) = levels(dadosFen2$INTER)

matINV

#######################################################
##MODELO 180 SELEGEN## 
#O modelo ajustado é o seguinte: y = Xg + Za + Wb + Ti + e, em que y é o
#vetor de dados, g é o vetor dos efeitos de geração (assumidos como fixos) somados à
#média geral, a é o vetor dos efeitos genéticos aditivos de famílias (assumidos como
#aleatórios), b é o vetor dos efeitos de blocos (assumidos como aleatórios), i é vetor
#dos efeitos da interação famílias x gerações (aleatórios) e e é o vetor de erros ou
#resíduos (aleatórios). As letras maiúsculas representam as matrizes de incidência
#para os referidos efeitos.

library(asreml)

str(dadosFen)

dadosFen

dadosFen = read.csv("C:\\Users\\Volpato\\Desktop\\Artigo 3\\Script Lorena\\DadosFenVolpatoTESTE.csv", header = T)

dadosFen$GENERATION = factor(dadosFen$GENERATION)
dadosFen$GENO = factor(dadosFen$GENO)
dadosFen$REP = factor(dadosFen$REP)
dadosFen$SITE = factor(dadosFen$SITE)
dadosFen$BLOC= factor(dadosFen$BLOC)

modA <- asreml(fixed=YIELD ~ 1 + GENERATION,
                random = ~ GENO:GENERATION + GENO + BLOC,
                data=dadosFen)

summary(modA) ## RESULTADOS IDENTICOS AO MODELO 180 DO SELEGEN ##

##Obtenção dos BLUPS das progênies (valores genotípicos iguais do SELEGEN)
BLUPS = (predict(modA, classify="GENO", sed=T))$predictions$pvals[,1:2]

library(dplyr)
BLUPS = BLUPS %>%
  arrange(-predicted.value)
BLUPS


## VIA METODOLOGIA BAYESIANA ##

Niter=1000  #10000
Burnin=500  #5000
Thin=5
nefit=(Niter-Burnin)/Thin

priori=list(R=list(V=1, nu=0.002), 
            G=list(G1=list(V=1, nu=0.002), 
                   G2=list(V=1, nu=0.002),
                   G3=list(V=1, nu=0.002)))

priori

modBAY1 = MCMCglmm(fixed=YIELD ~ 1 + GENERATION,
                   random = ~ GENO:GENERATION + GENO + BLOC,
                   data=dadosFen,
                   prior=priori, family = c("gaussian"), DIC=TRUE, verbose=FALSE, nitt=Niter, thin=Thin, burnin=Burnin,pr=TRUE,singular.ok=TRUE)

summary(modBAY1) ##RESUTLADOS SEMELHANTES AO REML/BLUP

##########################################################

##MODELO 180 SELEGEN com efeitos de locais:progenies## 
#O modelo ajustado é o seguinte: y = Xg + Za + Wb + Ti + Lp + e

library(asreml)
dadosFen = read.csv("C:\\Users\\Volpato\\Desktop\\Artigo 3\\Script Lorena\\DadosFenVolpatoTESTEBMG.csv", header = T)
dadosFen$GENERATION = factor(dadosFen$GENERATION)
dadosFen$GENO = factor(dadosFen$GENO)
dadosFen$REP = factor(dadosFen$REP)
dadosFen$SITE = factor(dadosFen$SITE)
dadosFen$BLOC= factor(dadosFen$BLOC)

str(dadosFen)

dadosFen

modA1 <- asreml(fixed=YIELD ~ 1 + GENERATION,
               random = ~ GENO:GENERATION + GENO + GENO:SITE + BLOC,
               data=dadosFen)

summary(modA1) 

## VIA METODOLOGIA BAYESIANA ##

Niter=10000  #10000
Burnin=5000  #5000
Thin=5
nefit=(Niter-Burnin)/Thin

priori=list(R=list(V=1, nu=0.002), 
            G=list(G1=list(V=1, nu=0.002), 
                   G2=list(V=1, nu=0.002),
                   G3=list(V=1, nu=0.002),
                   G4=list(V=1, nu=0.002)))

priori

modBAY1 = MCMCglmm(fixed=NDF ~ 1 + GENERATION,
                   random = ~ GENO:GENERATION + GENO + GENO:SITE +BLOC,
                   data=dadosFen,
                   prior=priori, family = c("gaussian"), DIC=TRUE, verbose=FALSE, nitt=Niter, thin=Thin, burnin=Burnin,pr=TRUE,singular.ok=TRUE)

summary(modBAY1) ##RESUTLADOS SEMELHANTES AO REML/BLUP


##############################################################

## MODELOS MULTIGERAÇÕES ##
###MODELOS PARA SE OBTER OS BLUPS### 
###Ou seja, como se faz no Selegen mod 180, se computa primeiro
###os componentes de variância com o modelo referido acima
###e depois o Selegen aplica os pesos obtidos de acordo com os parametros
###obtidos (h2, Nrep, etc)

###Dúvida a ser esclarecida: É correto obter os BLUPs incluindo o efeito "GENO" no modelo?
###Ou os BLUPs devem serem obtidos a partir da estrutura de VCOV sugerida (GENO:GENERATION ou GENO:INTER)? (acredito que é mais coerente fazer o modelo sem o efeito GENO.


### Indice Multigeracoes (segundo teoria do Índice Multigerações)
###Verificar convergência do modelo incluindo o efeito GENO:SITE

library(asreml)
modIM <- asreml(fixed=YIELD ~ 1 + REP,
                random = ~ giv(INTER) + GENERATION:REP,
                ginverse=list(INTER = matINV), data=dadosFen)

varComp = summary(modIM)$varcomp
varComp

##Obtenção dos BLUPS das progênies (BLUPS SEMELHANTES, MAS NÃO IGUAIS, PORÉM ORDENAMENTO SEMELHANTES)
BLUPS = (predict(modIM, classify="INTER", sed=T))$predictions$pvals[,1:2]

library(dplyr)
BLUPS = BLUPS %>%
  arrange(-predicted.value)
BLUPS
###########################################

##FREQUENTISTA##

modus <- asreml(fixed=NDF ~ 1 + REP,
                random = ~us(GENERATION):GENO + GENO:SITE + GENERATION:REP,
                rcov = ~ us(GENERATION):REP:GENO,
                data=dadosFen)
varComp = summary(modSC)$varcomp
varComp


##BAYESIANO##
###################################################################

dadosFen = read.csv("C:\\Users\\Volpato\\Desktop\\Artigo 3\\Script Lorena\\DadosFenVolpatoTESTEBMG.csv", header = T)

dadosFen$GENERATION = factor(dadosFen$GENERATION)
dadosFen$GENO = factor(dadosFen$GENO)
dadosFen$REP = factor(dadosFen$REP)
dadosFen$SITE = factor(dadosFen$SITE)
dadosFen$BLOC = factor(dadosFen$BLOC)
str(dadosFen)

## Modelo Bayesiano MultGerações ##

library(MCMCglmm)

Niter=1000000 
Burnin=500000 
Thin=5
nefit=(Niter-Burnin)/Thin

priori2=list(R=list(V=diag(3)/2, nu=2), 
             G=list(G1=list(V=diag(3)/2, nu=2), 
                    G2=list(V=1, nu=0.002),
                    G3=list(V=1, nu=0.002)))

priori2

modBAY1 = MCMCglmm(fixed=NDF ~ 1 + REP,
                   random = ~ us(GENERATION):GENO + SITE:GENO + GENERATION:REP,
                   rcov = ~ us(GENERATION):units, data=dadosFen,
                   prior=priori2, family = c("gaussian"), DIC=TRUE, verbose=FALSE, nitt=Niter, thin=Thin, burnin=Burnin,pr=TRUE,singular.ok=TRUE)

summary(modBAY1)

##MODELO PARA SE OBTER A CORRELAÇÃO DAS 3 GERAÇÕES##
##Na verdade é um modelo de Simetria composta (SC), ou seja, VCOV na diagonal constante e fora igual a 0
##Esse modelo foi escolhido pois não consigo enxergar como obter a correlação dos genótipos nas três gerações,
##e não par a par como sai direto no modelo US.

modSC <- asreml(fixed=NDF ~ 1 + REP,
                random = ~cor(GENERATION):GENO + BLOC,
                rcov = ~ cor(GENERATION):REP:GENO,
                data=dadosFen)
varComp = summary(modSC)$varcomp
varComp

## Correlação genética:
(rgger = modSC1$G.param$`GENERATION:GENO`$GENERATION$initial) 

###VALE A PENA OBSERVAR OS RESULTADOS DO SEGUINTE MODELO, E VERIFICAR SUA UTILIDADE (SE É MAIS CORRETO QUE OS ANTERIORES)
##Cada geração uma variável
##Dados para váriavel PROD (em sacos por ha). Obs:demais análises fiz considerando produção por plot (mais correto).
##Modelo multigeneration and multienvironment (MGME)

##Frequentista##

dadosFen = read.csv("C:\\Users\\Volpato\\Desktop\\Artigo 3\\Script Lorena\\DadosFenVolpato2MTME.csv", header = T)

dadosFen$GENO = factor(dadosFen$GENO)
dadosFen$REP = factor(dadosFen$REP)
dadosFen$SITE = factor(dadosFen$SITE)
dadosFen$INTER = factor(dadosFen$INTER)

library(asreml)

modUS2 = asreml(fixed = cbind(S1,S2,S3) ~ trait:REP - 1,
                  random = ~ us(trait):GENO + us(trait):INTER + GENO,
                  rcov = ~ us(trait):units, data=dadosFen)
                 
varComp = summary(modUS2)$varcomp
varComp

##Bayesiano##

library(MCMCglmm)

Niter=10000  #10000
Burnin=5000  #5000
Thin=5
nefit=(Niter-Burnin)/Thin

priori=list(R=list(V=diag(3)/2, nu=2), 
            G=list(G1=list(V=diag(3)/2, nu=2), 
                   G2=list(V=diag(3)/2, nu=2),
                   G3=list(V=1, nu=0.002)))

priori

modBAY = MCMCglmm(fixed = cbind(S1,S2,S3) ~ trait:REP - 1,
                  random = ~ us(trait):GENO + us(trait):INTER + GENO,
                  rcov = ~ us(trait):units, data=dadosFen,
                  prior=priori, family = c("gaussian", "gaussian","gaussian"), DIC=TRUE, verbose=FALSE, nitt=Niter, thin=Thin, burnin=Burnin,pr=TRUE,singular.ok=TRUE)

summary(modBAY)
