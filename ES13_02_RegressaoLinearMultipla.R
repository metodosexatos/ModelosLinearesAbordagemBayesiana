#********************************************************#
# www.metodosexatos.com                                  #
# Prof.Ms. Andr� Santos | andre@metodosexatos.com.br     #
# Curso Modelos Lineares com Abordagem Bayesiana         #
# Modelo Linear de Regress�o M�ltipla (MLRM)             #
# Data: 02/08/2020                                       #
#********************************************************#

# I. Entrada dos dados

# Op��o A .....................................................................

# Diret�rios e Arquivos:
# getwd() # Qual o diret�rio que o script est� apontando
# list.files() # Quais arquivos est�o contidos no diret�rio
# setwd("C:/Users/andre/Downloads")
# Leitura de uma base externa
# dados <- read.csv2(file = "fosforo.csv")

# Op��o B .....................................................................

dados <- data.frame(disp=c(64,60,71,61,54,77,81,93,93,51,76,96,77,93,95,54,99),
                    inorg=c(0.4,0.4,3.1,0.6,4.7,1.7,9.4,10.1,11.6,12.6,10.9,
                            23.1,23.1,21.6,23.1,1.9,29.9),
                    org=c(53,23,19,34,24,65,44,31,29,58,37,46,50,44,56,36,51))

k <- ncol(dados)
#str(dados)
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# 1. Diagramas de dispers�o

plot(cbind(dados,log(dados)))   # escala original vs escala logar�tmica
plot(dados)                     # dados na escala original
plot(log(dados))                # dados na escala logar�tmica
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# 2. An�lise de regress�o bayesiana com vari�veis padronizadas

# a) Vari�veis explicativas padronizadas
dados$z1 <- (dados$inorg-mean(dados$inorg))/sd(dados$inorg)
dados$z2 <- (dados$org-mean(dados$org))/sd(dados$org)

# b) Modelo linear com abordagem cl�ssica (estat�stica frequentista)
modelo <- lm(disp~z1+z2,data=dados)

# c) Extra��o do erro padr�o, coeficientes e matriz de covari�ncia
se <- summary(modelo)$sigma
est <- coef(modelo)
v.beta <- summary(modelo)$cov.unscaled

# d) Tamanho da amostra da distribui��o posterior
m <- 1000

# e) Precis�o da distribui��o Gamma e Vari�ncia
prec <- rgamma(m,(nrow(dados)-k)/2,se^2*(nrow(dados)-k)/2)
sig2 <- 1/prec

# f) N�meros aleat�rios de uma normal multivariada e coeficientes do MLRM
library(MASS)
beta.conj <- numeric()
for (i in 1:m) {
  beta.conj <- rbind(beta.conj,mvrnorm(1,est,sig2[i]*v.beta))
  
}

# g) Primeiras linhas da amostra de tamanho m da dist. posterior conjunta.
cbind(beta.conj,sig2)[1:3,]
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# 3. Representa��o gr�fica da amostra da distribui��o posterior

plot(data.frame(cbind(beta.conj,sig2)), cex=0.4)
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# 4. Predi��o para novos valores das vari�veis explicativas (x)

# a) Quantidade das vari�veis explicativas que queremos fazer a predi��o:
x1p <- 15
x2p <- 40

# b) Padroniza��o das vari�veis explicativas
z1p <- (x1p-mean(dados$inorg))/sd(dados$inorg)
z2p <- (x2p-mean(dados$org))/sd(dados$org)

# c) M�dia da distribui��o preditiva para cada uma das m amostras
mup <- beta.conj%*%c(1,z1p,z2p)

# d) Retirada da amostra da distribui��o posterior para cada m amostras
ytil <- numeric()
for(i in 1:length(mup)){
  ytil <- c(ytil,rnorm(1,mup[i],sqrt(sig2[i])))
}

# e) Primeiros valores da amostra preditiva
ytil[1:3]
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# 5. Representa��o gr�fica e sum�rio da distribui��o posterior preditiva

# a) Histograma
hist(ytil,main="")

# b) sum�rio
quantile(ytil,prob=c(.025,.25,.5,.75,.975))
