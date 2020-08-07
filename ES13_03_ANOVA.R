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
martelo <- read.csv2(file = "martelo.csv")
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# 1. Análise preliminar (seleção da variável resposta)

attach(martelo) # base de dados anexada

# sumário para os níveis das variáveis de interesse

# variável original
data.frame(media=tapply(cpue, area, mean),
           mediana=tapply(cpue, area, median), var=tapply(cpue,area,var))

# logarítmo da variável cpue
data.frame(media=tapply(lcpue, area, mean),
           mediana=tapply(lcpue, area, median), var=tapply(lcpue,area,var))

# análise visual (boxplot)
par(mfrow = c(1,2))
boxplot(cpue~area, main = "Dados originais")
boxplot(lcpue~area, main = "Logarítmo de CPUE")

#detach(martelo) # base de dados desanexada
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# 2. Estatísticas básicas para distribuições posteriores

fnsq <- function(x) sum((x-mean(x))^2)   # função para SQR
ygbar <- tapply(lcpue, area, mean)       # médias
sqg <- tapply(lcpue, area, fnsq)         # soma dos quadrados dos resíduos(SQR)
ng <- tapply(lcpue, area, length)        # número de amostras
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

# 5. Sumário estatístico para os efeitos dos níveis do fator CPUE

for (i in 1:ncol(emuGpost)) {
  print(c(quantile(emuGpost[,i], prob=c(0.025, 0.5, 0.975)),
          media=mean(emuGpost[,i]), dp=sd(emuGpost[,i])))
}
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# 6. Gráfico com as marginais posteriores para os efeitos

par(mfrow = c(1,1))
plot(density(emuGpost[,1]),xlab ="mu",xlim=c(30,300),
     ylab ="p(mu|D)",main ="")                                  # 1º nível
lines(density(emuGpost[,2]),col="red")                          # 2º nível  
lines(density(emuGpost[,3]),col="blue")                         # 3º nível
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# 7. Probabilidade que mu3 (3º nível) seja maior que mu1 (1º nível)

sum(emuGpost[,3] > emuGpost[,1])/m
