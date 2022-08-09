library(tidyverse)
library(ggplot2)
set.seed(250)
#10a

RNGkind()
RNGkind(kind="Wichmann-Hill")
RNGkind()

#Stimulating Binomial variable using Inversion method
N = 1000
n = 10
p = 0.4
u <- runif(N)
bins <- .bincode(u, breaks = c(0, pbinom(0:10, 10, 0.4)), right = F, include.lowest = T)

Inversion.bin <- numeric()
for(i in 1:N){
  Inversion.bin[i] <- bins[i]-1
}
#Simulating a binomial random variable by simulating corresponding Bernoulli random variables by inversion method
Bernoulli.bin <- numeric()
for (i in 1:N){
  v <- runif(n)
  Bernoulli.bin[i] <- sum(v < p)
}
#Simulating of a binomial random variable with rbinom
Binomial.rbinom <- rbinom(N, n, p)

#Plot 
Inversion.bin <- data.frame(Inversion.bin)
Inversion.bin$method <- rep("Inverse CDF", 1000)
colnames(Inversion.bin) <- c("rand_num", "method")

Bernoulli.bin <- data.frame(Bernoulli.bin)
Bernoulli.bin$method <- rep(" Bernoulli", 1000)
colnames(Bernoulli.bin) <- c("rand_num", "method")

Binomial.rbinom <- data.frame(Binomial.rbinom)
Binomial.rbinom$method <- rep("rbinom", 1000)
colnames(Binomial.rbinom) <- c("rand_num", "method")

Dataframe <- rbind(Inversion.bin, Bernoulli.bin, Binomial.rbinom)

 ggplot(Dataframe, aes(x = rand_num, fill = method)) +
  geom_histogram( binwidth=.5, position="dodge") 

#Switch the random number generator back to its default
RNGkind(kind = "default", normal.kind = NULL)

#10b

f <- function(x){
  ((2*pi)^(-1/2))*exp(-(x^2)/2)
}

g <- function(x){
  (pi*(1 + x^2))^(-1)
}

#Determining the best value of the constant c, such that f(x) <= cg(x)
x <- -1500:1500
c <- max(f(x)/g(x))

#10000 standard normal random variables using accept-reject method, generating Cauchy distributed random variables using inversion method
N <- 10000
j <- 0
Rand.num <- numeric()
while(length(Rand.num) != N){
  w <- runif(1) #step 1
  cauchy <- tan((w-(1/2))*pi) 
  U <- runif(1) #step 2
  if(U*c*g(cauchy) <= f(cauchy)){
    Rand.num[j] <- cauchy
    j <- j + 1
  }
}

#Histogram 
k <- rnorm(N)
Dataframe <- data.frame(Rand.num, k)
 ggplot(Dataframe) +
  geom_histogram(aes( x = Rand.num, y = ..density.., colour = Rand.num), colour ="white") +
  geom_density(aes(x = k), colour = "blue")

#QQ-plot
 ggplot(data = Dataframe, mapping = aes(sample = Rand.num)) +
  stat_qq()


