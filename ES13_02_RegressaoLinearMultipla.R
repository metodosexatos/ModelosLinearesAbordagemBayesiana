#********************************************************#
# www.metodosexatos.com                                  #
# Prof.Ms. André Santos | andre@metodosexatos.com.br     #
# Curso Modelos Lineares com Abordagem Bayesiana         #
# Modelo Linear de Regressão Múltipla (MLRM)             #
# Data: 02/08/2020                                       #
#********************************************************#

# I. Entrada dos dados

# Opção A .....................................................................

# Diretórios e Arquivos:
# getwd() # Qual o diretório que o script está apontando
# list.files() # Quais arquivos estão contidos no diretório
# setwd("C:/Users/andre/Downloads")
# Leitura de uma base externa
# dados <- read.csv2(file = "fosforo.csv")

# Opção B .....................................................................

dados <- data.frame(disp=c(64,60,71,61,54,77,81,93,93,51,76,96,77,93,95,54,99),
                    inorg=c(0.4,0.4,3.1,0.6,4.7,1.7,9.4,10.1,11.6,12.6,10.9,
                            23.1,23.1,21.6,23.1,1.9,29.9),
                    org=c(53,23,19,34,24,65,44,31,29,58,37,46,50,44,56,36,51))

k <- ncol(dados)
#str(dados)
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# 1. Diagramas de dispersão

plot(cbind(dados,log(dados)))   # escala original vs escala logarítmica
plot(dados)                     # dados na escala original
plot(log(dados))                # dados na escala logarítmica
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# 2. Análise de regressão bayesiana com variáveis padronizadas

# a) Variáveis explicativas padronizadas
dados$z1 <- (dados$inorg-mean(dados$inorg))/sd(dados$inorg)
dados$z2 <- (dados$org-mean(dados$org))/sd(dados$org)

# b) Modelo linear com abordagem clássica (estatística frequentista)
modelo <- lm(disp~z1+z2,data=dados)

# c) Extração do erro padrão, coeficientes e matriz de covariância
se <- summary(modelo)$sigma
est <- coef(modelo)
v.beta <- summary(modelo)$cov.unscaled

# d) Tamanho da amostra da distribuição posterior
m <- 1000

# e) Precisão da distribuição Gamma e Variância
prec <- rgamma(m,(nrow(dados)-k)/2,se^2*(nrow(dados)-k)/2)
sig2 <- 1/prec

# f) Números aleatórios de uma normal multivariada e coeficientes do MLRM
library(MASS)
beta.conj <- numeric()
for (i in 1:m) {
  beta.conj <- rbind(beta.conj,mvrnorm(1,est,sig2[i]*v.beta))
  
}

# g) Primeiras linhas da amostra de tamanho m da dist. posterior conjunta.
cbind(beta.conj,sig2)[1:3,]
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# 3. Representação gráfica da amostra da distribuição posterior

plot(data.frame(cbind(beta.conj,sig2)), cex=0.4)
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# 4. Predição para novos valores das variáveis explicativas (x)

# a) Quantidade das variáveis explicativas que queremos fazer a predição:
x1p <- 15
x2p <- 40

# b) Padronização das variáveis explicativas
z1p <- (x1p-mean(dados$inorg))/sd(dados$inorg)
z2p <- (x2p-mean(dados$org))/sd(dados$org)

# c) Média da distribuição preditiva para cada uma das m amostras
mup <- beta.conj%*%c(1,z1p,z2p)

# d) Retirada da amostra da distribuição posterior para cada m amostras
ytil <- numeric()
for(i in 1:length(mup)){
  ytil <- c(ytil,rnorm(1,mup[i],sqrt(sig2[i])))
}

# e) Primeiros valores da amostra preditiva
ytil[1:3]
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# 5. Representação gráfica e sumário da distribuição posterior preditiva

# a) Histograma
hist(ytil,main="")

# b) sumário
quantile(ytil,prob=c(.025,.25,.5,.75,.975))
