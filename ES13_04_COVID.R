#********************************************************#
# www.metodosexatos.com                                  #
# Prof.Ms. Andr� Santos | andre@metodosexatos.com.br     #
# Curso Modelos Lineares com Abordagem Bayesiana         #
# Modelo Linear de Regress�o M�ltipla (MLRM)             #
# Data: 05/08/2020                                       #
#********************************************************#

# I. Entrada dos dados

# Diret�rios e Arquivos:
# getwd() # Qual o diret�rio que o script est� apontando
# list.files() # Quais arquivos est�o contidos no diret�rio
# setwd("C:/Users/andre/Downloads")
# Leitura de uma base externa
covid <- read.csv2(file = "coronavirus.csv")
str(covid)
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# 1. An�lise preliminar (sele��o da vari�vel resposta)

covid$taxa <- covid$obitos/covid$habitantes*10^6
covid$ltaxa <- log(covid$taxa)
head(covid)
attach(covid) # base de dados anexada

# sum�rio para os n�veis das vari�veis de interesse

# vari�vel original
data.frame(media=tapply(taxa, regiao, mean),
           mediana=tapply(taxa, regiao, median), var=tapply(taxa,regiao,var))

# logar�tmo da vari�vel taxa
data.frame(media=tapply(ltaxa, regiao, mean),
           mediana=tapply(ltaxa, regiao, median), var=tapply(ltaxa,regiao,var))

# an�lise visual (boxplot)
par(mfrow = c(1,2))
boxplot(taxa~regiao, main = "Dados originais")
boxplot(ltaxa~regiao, main = "Logar�tmo da taxa")

#detach(base) # base de dados desanexada
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# 2. Estat�sticas b�sicas para distribui��es posteriores

fnsq <- function(x) sum((x-mean(x))^2)   # fun��o para SQR
ygbar <- tapply(ltaxa, regiao, mean)     # m�dias
sqg <- tapply(ltaxa, regiao, fnsq)       # soma dos quadrados dos res�duos(SQR)
ng <- tapply(ltaxa, regiao, length)      # n�mero de amostras
nG <- sum(ng)                            # n�mero total de amostras
G <- length(ng)                          # n�mero de n�veis do fator
se2G <- sum(sqg)/(nG-G)                  # vari�ncia das estimativas
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# 3. Amostra simulada das distribui��es posteriores marginais

Vmu <- diag(1/ng)                          # matriz diagonal
m <- 3000                                  # amostras
tauG <- rgamma(m, (nG-G)/2, se2G*(nG-G)/2) # amostra da precis�o
varG <- 1/tauG                             # vari�ncia
sigma <- sqrt(1/tauG)                      # desvio padr�o

# retirada da amostra:
library(MASS)           # biblioteca para retirada de n� aleat�rio (mvrnorm)
muGpost <- numeric()    # amostra da posterior dos efeitos (m�dia g)
for (i in 1:m) {
  muGpost <- rbind(muGpost,mvrnorm(1,as.numeric(ygbar), varG[i]*Vmu))
}
emuGpost <- exp(muGpost)
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# 4. Sum�rio estat�stico para a amostra da posterior do desvio padr�o (sigma)

# Intervalo de credibilidade (ICr95% -> [2.5%; 97.5%]):
exp(c(quantile(sigma,prob=c(0.025,0.5,0.975)),media=mean(sigma),dp=sd(sigma)))

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# 5. Sum�rio estat�stico para os efeitos dos n�veis do fator taxa

for (i in 1:ncol(emuGpost)) {
  print(c(quantile(emuGpost[,i], prob=c(0.025, 0.5, 0.975)),
          media=mean(emuGpost[,i]), dp=sd(emuGpost[,i])))
}
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# 6. Gr�fico com as marginais posteriores para os efeitos

par(mfrow = c(1,1))
plot(density(emuGpost[,1]),xlab = "mu",
     xlim=c(0,7000),ylab = "p(mu|D)", main = "")   # centro-oeste
lines(density(emuGpost[,2]), col="red")           # nordeste
lines(density(emuGpost[,3]), col="blue")          # norte
lines(density(emuGpost[,4]), col="green")         # sudeste
lines(density(emuGpost[,5]), col="purple")        # sul

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# 7. Probabilidades entre n�veis

# que mu3 (norte) seja maior que mu4 (sudeste)
sum(emuGpost[,3] > emuGpost[,4])/m
# que mu1 (centro-oeste) seja maior que mu5 (sul)
sum(emuGpost[,1] > emuGpost[,5])/m










