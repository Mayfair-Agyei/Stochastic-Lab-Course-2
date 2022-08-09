install.packages("bootstrap")
library(tidyverse)
library(bootstrap)


SHHS <- read.delim("R/Stochastic Lab Course/shhs2.txt", sep = '\t'  )

#9a
# Monte Carlo samples and significance level
M <- 1000        
alpha <- 0.95 

# Weibull distribution parameters
lambda <- 13
k <- 1          

x = rweibull(100, shape = k, scale = lambda)
sigma = sd(x)
x.med = median(x)
c(sigma, x.med)

#Bootstrap Percentile Confidence interval and Coverage Probability for n=100, R=1000

n <- 100        
R <- 1000       

#Two-sided bootstrap percentile confidence intervals for sigma and x.med
Ts.Median <- 0
Ts.SD <- 0


ConIntMedian.left <- 0
ConIntMedian.right <- 0
ConIntSd.left <- 0
ConIntSd.right <- 0

for (j in 1:M) {
  Weibull.Samples <- rweibull(n, k, lambda)
  Med <- 0
  SD <- 0
  for (i in 1:R) {
    Bootstrap.Values <- sample(Weibull.Samples, n, replace = T)
    Med[i] <- median(Bootstrap.Values)
    SD[i] <- sd(Bootstrap.Values)
  }
  
  Med <- sort(Med)
  SD <- sort(SD)
  
  ConIntMedian.left[j] <- Med[floor(R * (1-alpha))/2]
  ConIntMedian.right[j] <- Med[floor(R * (1-(1-alpha)/2))]
  ConIntSd.left[j] <- SD[floor(R * (1-alpha))/2]
  ConIntSd.right[j] <- SD[floor(R * (1-(1-alpha)/2))]
  
  ConInts <- data.frame(ConIntMedian.left, ConIntMedian.right, ConIntSd.left, ConIntSd.right)
  
}

ConInts <- ConInts %>% mutate(Ts.Median = x.med >= ConIntMedian.left & x.med <= ConIntMedian.right)

ConInts <- ConInts %>% mutate(Ts.SD = sigma >= ConIntSd.left & sigma <= ConIntSd.right)


#Estimating Coverage probability and Average interval length for n=100, R=1000
Coverage.Prob <- c(sum(ConInts$Ts.Median)/ M, sum(ConInts$Ts.SD)/M)


AverageInt.length <- c(sum(ConInts$ConIntMedian.right - ConInts$ConIntMedian.left)/nrow(ConInts), 
                             sum(ConInts$ConIntSd.right - ConInts$ConIntSd.left)/nrow(ConInts))

c(Coverage.Prob, AverageInt.length)

# Bootstrap Percentile Confidence interval and Coverage Probability for n=R=1000
n <- 1000     
R <- 1000      

#Two-sided bootstrap percentile confidence intervals for sigma and x.med
Ts.Median <- 0
Ts.SD <- 0


ConIntMedian.left <- 0
ConIntMedian.right <- 0
ConIntSd.left <- 0
ConIntSd.right <- 0

for (j in 1:M) {
  Weibull.Samples <- rweibull(n, k, lambda)
  Med <- 0
  SD <- 0
  for (i in 1:R) {
    Bootstrap.Values <- sample(Weibull.Samples, n, replace = T)
    Med[i] <- median(Bootstrap.Values)
    SD[i] <- sd(Bootstrap.Values)
  }
  
  Med <- sort(Med)
  SD <- sort(SD)
  
  ConIntMedian.left[j] <- Med[floor(R * (1-alpha))/2]
  ConIntMedian.right[j] <- Med[floor(R * (1-(1-alpha)/2))]
  ConIntSd.left[j] <- SD[floor(R * (1-alpha))/2]
  ConIntSd.right[j] <- SD[floor(R * (1-(1-alpha)/2))]
  
  ConInts <- data.frame(ConIntMedian.left, ConIntMedian.right, ConIntSd.left, ConIntSd.right)
  
}

ConInts <- ConInts %>% mutate(Ts.Median = x.med >= ConIntMedian.left & x.med <= ConIntMedian.right)

ConInts <- ConInts %>% mutate(Ts.SD = sigma >= ConIntSd.left & sigma <= ConIntSd.right)

#Estimating Coverage probability and Average interval length for n=R=1000
Coverage.Prob <- c(sum(ConInts$Ts.Median)/ M, sum(ConInts$Ts.SD)/M)


AverageInt.length <- c(sum(ConInts$ConIntMedian.right - ConInts$ConIntMedian.left)/nrow(ConInts), 
                             sum(ConInts$ConIntSd.right - ConInts$ConIntSd.left)/nrow(ConInts))

c(Coverage.Prob, AverageInt.length)

# Bootstrap Percentile Confidence interval and Coverage Probability for n=100 R=5000
n <- 100        
R <- 5000     

#Two-sided bootstrap percentile confidence intervals for sigma and x_med
Ts.Median <- 0
Ts.SD <- 0


ConIntMedian.left <- 0
ConIntMedian.right <- 0
ConIntSd.left <- 0
ConIntSd.right <- 0

for (j in 1:M) {
  Weibull.Samples <- rweibull(n, k, lambda)
  Med <- 0
  SD <- 0
  for (i in 1:R) {
    Bootstrap.Values <- sample(Weibull.Samples, n, replace = T)
    Med[i] <- median(Bootstrap.Values)
    SD[i] <- sd(Bootstrap.Values)
  }
  
  Med <- sort(Med)
  SD <- sort(SD)
  
  ConIntMedian.left[j] <- Med[floor(R * (1-alpha))/2]
  ConIntMedian.right[j] <- Med[floor(R * (1-(1-alpha)/2))]
  ConIntSd.left[j] <- SD[floor(R * (1-alpha))/2]
  ConIntSd.right[j] <- SD[floor(R * (1-(1-alpha)/2))]
  
  ConInts <- data.frame(ConIntMedian.left, ConIntMedian.right, ConIntSd.left, ConIntSd.right)
  
}

ConInts <- ConInts %>% mutate(Ts.Median = x.med >= ConIntMedian.left & x.med <= ConIntMedian.right)

ConInts <- ConInts %>% mutate(Ts.SD = sigma >= ConIntSd.left & sigma <= ConIntSd.right)

#Estimating Coverage probability and Average interval length for n=100, R=5000
Coverage.Prob <- c(sum(ConInts$Ts.Median)/ M, sum(ConInts$Ts.SD)/M)

AverageInt.length <- c(sum(ConInts$ConIntMedian.right - ConInts$ConIntMedian.left)/nrow(ConInts), 
                             sum(ConInts$ConIntSd.right - ConInts$ConIntSd.left)/nrow(ConInts))
c(Coverage.Prob, AverageInt.length)


#9aii
#Buliding Bootstrap Accelerated Bias-Corrected Interval

ConIntMedian.left <- rep(0, M)
ConIntMedian.right <- rep(0, M)
zHat.0Median <- 0
aHat.Median <- 0

for (j in 1:M) {
  Weibull.Samples <- rweibull(100, k, lambda)
  A <- bcanon(Weibull.Samples, R, theta = median, alpha = c(0.025, 0.975))
  zHat.0Median[j] <- A$z0
  aHat.Median[j] <- A$acc
  ConIntMedian.left[j] <- A$confpoints[1,2]
  ConIntMedian.right[j] <- A$confpoints[2,2]
}

DataFramebca.PointMedian <- data.frame(ConIntMedian.left, ConIntMedian.right)

DataFramebca.PointMedian <- DataFramebca.PointMedian %>% mutate(Ts.Median = x.med >= ConIntMedian.left & x.med <= ConIntMedian.right)
#Using the Standard Deviation
ConIntSd.left <- 0
ConIntSd.right <- 0
zHat.0SD <- 0
aHat.SD <- 0
for (j in 1:M) {
  Weibull <- rweibull(100, k, lambda)
  A <- bcanon(Weibull, R, theta = sd, alpha = c(0.025, 0.975))
  zHat.0SD[j] <- A$z0
  aHat.SD[j] <- A$acc
  ConIntSd.left[j] <- A$confpoints[1,2]
  ConIntSd.right[j] <- A$confpoints[2,2]
}

DataFramebca.PointSD <- data.frame(ConIntSd.left, ConIntSd.right)

DataFramebca.PointSD <- DataFramebca.PointSD %>% 
  mutate(Ts.SD = sigma >= ConIntSd.left & sigma <= ConIntSd.right)


Estimate <-  data.frame(zHat.0SD, aHat.SD, zHat.0Median, aHat.Median)
Statistics <- data.frame(DataFramebca.PointMedian, DataFramebca.PointSD) 

#Estimating Coverage probability and Average interval length 
Coverage.Prob <- c(sum(Statistics$Ts.Median)/M, sum(Statistics$Ts.SD)/M)
AverageInt.length <- c(sum(Statistics$ConIntMedian.right - Statistics$ConIntMedian.left)/nrow(Statistics), 
                             sum(Statistics$ConIntSd.right - Statistics$ConIntSd.left)/nrow(Statistics))

c(mean(zHat.0SD),mean(aHat.SD), mean(zHat.0Median), mean(aHat.Median))

#9b

#Histogram and Empirical distribution of the variable rdi4p
ggplot(SHHS, aes(x = rdi4p, y = ..density..)) +
  geom_histogram(color = 'white') +
  geom_density(aes(x = rdi4p), colour = "blue") + 
  labs(title = "histogram and empirical distribution of the variable rdi4p",
           x = "Sample",y = "Relative Frequency", colour = "Method") + theme_gray()


#Buildng Bootstrap percentile (two-sided) confidence intervals for the standard deviation and median
rdi4p.median <- median(SHHS$rdi4p)
rdi4p.sigma <- sd(SHHS$rdi4p)
alpha <- 0.95

n <- length(SHHS$rdi4p) 
R <- 1000 

Med <- 0
SD <- 0
for (i in 1:R) {
  Bootstrap.Values <- sample(SHHS$rdi4p, n, replace = T)
  Med[i] <- median(Bootstrap.Values)
  SD[i] <- sd(Bootstrap.Values)
}

Med <- sort(Med)
SD <- sort(SD)
ConIntMedian.left <- Med[floor(R *(1-alpha))/2]
ConIntMedian.right <- Med[floor(R *(1-(1-alpha)/2))]
ConIntSd.left <- SD[floor(R *(1-alpha))/2]
ConIntSd.right <- SD[floor(R *(1-(1-alpha)/2))]

boostrap_percentile_CI <- data.frame(ConIntMedian.left, ConIntMedian.right, ConIntSd.left, ConIntSd.right)
Ts.Median <- rdi4p.median >= ConIntMedian.left & rdi4p.median <= ConIntMedian.right
Ts.SD <- rdi4p.sigma >= ConIntSd.left & rdi4p.sigma <= ConIntSd.right

#Bootstrap accelerated bias-corrected confidence intervals 

# Median
rdi4p.median.rem <- SHHS$rdi4p [! SHHS$rdi4p %in% median(SHHS$rdi4p)]
A <- bcanon(rdi4p.median.rem, R, theta = median, alpha = c(0.025, 0.975))
z.Median <- A$z0
a.0Median <- A$acc
ConIntMedian.left <- A$confpoints[1,2]
ConIntMedian.right <- A$confpoints[2,2]

Ts.MedianA = rdi4p.median >= ConIntMedian.left & rdi4p.median <= ConIntMedian.right
# Standard deviation
A <- bcanon(SHHS$rdi4p, R, theta = sd, alpha = c(0.025, 0.975))
z.SD <- A$z0
a.0SD <- A$acc
ConIntSd.left <- A$confpoints[1,2]
ConIntSd.right <- A$confpoints[2,2]

Ts.SDA = rdi4p.sigma >= ConIntSd.left & rdi4p.sigma <= ConIntSd.right

BCa_CI <- data.frame(ConIntMedian.left, ConIntMedian.right, ConIntSd.left, ConIntSd.right, z.Median, a.0Median, z.SD, a.0SD)








