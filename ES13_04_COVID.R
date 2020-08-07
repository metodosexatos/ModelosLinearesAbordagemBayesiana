#********************************************************#
# www.metodosexatos.com                                  #
# Prof.Ms. André Santos | andre@metodosexatos.com.br     #
# Curso Modelos Lineares com Abordagem Bayesiana         #
# Modelo Linear de Regressão Múltipla (MLRM)             #
# Data: 05/08/2020                                       #
#********************************************************#

# I. Entrada dos dados

# Diretórios e Arquivos:
# getwd() # Qual o diretório que o script está apontando
# list.files() # Quais arquivos estão contidos no diretório
# setwd("C:/Users/andre/Downloads")
# Leitura de uma base externa
covid <- read.csv2(file = "coronavirus.csv")
str(covid)
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# 1. Análise preliminar (seleção da variável resposta)

covid$taxa <- covid$obitos/covid$habitantes*10^6
covid$ltaxa <- log(covid$taxa)
head(covid)
attach(covid) # base de dados anexada

# sumário para os níveis das variáveis de interesse

# variável original
data.frame(media=tapply(taxa, regiao, mean),
           mediana=tapply(taxa, regiao, median), var=tapply(taxa,regiao,var))

# logarítmo da variável taxa
data.frame(media=tapply(ltaxa, regiao, mean),
           mediana=tapply(ltaxa, regiao, median), var=tapply(ltaxa,regiao,var))

# análise visual (boxplot)
par(mfrow = c(1,2))
boxplot(taxa~regiao, main = "Dados originais")
boxplot(ltaxa~regiao, main = "Logarítmo da taxa")

#detach(base) # base de dados desanexada
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# 2. Estatísticas básicas para distribuições posteriores

fnsq <- function(x) sum((x-mean(x))^2)   # função para SQR
ygbar <- tapply(ltaxa, regiao, mean)     # médias
sqg <- tapply(ltaxa, regiao, fnsq)       # soma dos quadrados dos resíduos(SQR)
ng <- tapply(ltaxa, regiao, length)      # número de amostras
nG <- sum(ng)                            # número total de amostras
G <- length(ng)                          # número de níveis do fator
se2G <- sum(sqg)/(nG-G)                  # variância das estimativas
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# 3. Amostra simulada das distribuições posteriores marginais

Vmu <- diag(1/ng)                          # matriz diagonal
m <- 3000                                  # amostras
tauG <- rgamma(m, (nG-G)/2, se2G*(nG-G)/2) # amostra da precisão
varG <- 1/tauG                             # variância
sigma <- sqrt(1/tauG)                      # desvio padrão

# retirada da amostra:
library(MASS)           # biblioteca para retirada de nº aleatório (mvrnorm)
muGpost <- numeric()    # amostra da posterior dos efeitos (média g)
for (i in 1:m) {
  muGpost <- rbind(muGpost,mvrnorm(1,as.numeric(ygbar), varG[i]*Vmu))
}
emuGpost <- exp(muGpost)
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# 4. Sumário estatístico para a amostra da posterior do desvio padrão (sigma)

# Intervalo de credibilidade (ICr95% -> [2.5%; 97.5%]):
exp(c(quantile(sigma,prob=c(0.025,0.5,0.975)),media=mean(sigma),dp=sd(sigma)))

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# 5. Sumário estatístico para os efeitos dos níveis do fator taxa

for (i in 1:ncol(emuGpost)) {
  print(c(quantile(emuGpost[,i], prob=c(0.025, 0.5, 0.975)),
          media=mean(emuGpost[,i]), dp=sd(emuGpost[,i])))
}
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# 6. Gráfico com as marginais posteriores para os efeitos

par(mfrow = c(1,1))
plot(density(emuGpost[,1]),xlab = "mu",
     xlim=c(0,7000),ylab = "p(mu|D)", main = "")   # centro-oeste
lines(density(emuGpost[,2]), col="red")           # nordeste
lines(density(emuGpost[,3]), col="blue")          # norte
lines(density(emuGpost[,4]), col="green")         # sudeste
lines(density(emuGpost[,5]), col="purple")        # sul

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# 7. Probabilidades entre níveis

# que mu3 (norte) seja maior que mu4 (sudeste)
sum(emuGpost[,3] > emuGpost[,4])/m
# que mu1 (centro-oeste) seja maior que mu5 (sul)
sum(emuGpost[,1] > emuGpost[,5])/m










