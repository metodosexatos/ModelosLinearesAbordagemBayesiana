#********************************************************#
# www.metodosexatos.com                                  #
# Prof.Ms. Andr� Santos | andre@metodosexatos.com.br     #
# Curso Modelos Lineares com Abordagem Bayesiana         #
# Regress�o Linear Simples                               #
# Data: 01/08/2020                                       #
#********************************************************#

# 1.Entrada dos dados

x <- c(1.9, 3.1, 3.3, 4.8, 5.3, 6.1, 6.4, 7.6, 9.8, 12.4)
y <- c(2, 1, 5, 5, 20, 20, 23, 10, 30, 25)
#--------------------------------------------------------

# 2.Sum�rio do modelo de regress�o linear simples

modelo <- lm(y~x)
summary(modelo)
#--------------------------------------------------------

# 3. Par�metros de escala das distribui�es posteriores

# erro padr�o dos res�duos:
se <- summary(modelo)$sigma
se

# coeficientes:
est <- coef(modelo)
est
#--------------------------------------------------------

# 4. Simula��es posteriores conjuntas (opcional)

# matriz VB:
v.beta <- summary(modelo)$cov.unscaled
v.beta

# tamanho da amostra:
m <- 5000

# Precis�o de uma distribui��o gama:
prec <- rgamma(m, (length(x)-2)/2, se^2*(length(x)-2)/2)

# Valores de sigma� para uma amostra gama-inversa:
sig2 <- 1/prec

# Amostra da distribui��o posterior conjunta:
library(MASS)
set.seed(123)
beta.conj <- numeric()
for(i in 1:m){
  beta.conj <- rbind(beta.conj, mvrnorm(1,est, sig2[i]*v.beta))
}

# Os 3 primeiros vetores (beta 0, beta 1):
beta.conj[1:3,]
#--------------------------------------------------------

# 5.Determinar a distribui��o preditiva

# Medidas preliminares:
xp <- 11
n <- length(x)
xbar <- mean(x)
sxx <- mean(x*x)-xbar^2

# Escala e m�dia da distribui��o de student:
escala <- sqrt(se^2*((n+1)/n+(xp-xbar)^2/(n*sxx)))
media <- est[1]+est[2]*xp
c(media, escala)
#--------------------------------------------------------

# 6. Extrair informa��es de interesse (ICr95%)

qt(c(0.025,0.975), n-2)*escala+media