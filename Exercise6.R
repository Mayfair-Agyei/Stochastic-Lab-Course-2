library(tidyverse)
library(haven)
library(ggplot2)
library(NonpModelCheck)
library(dplyr)
install.packages("splines2")
install.packages("nlme")
library(splines)
library(nlme)


#reading data stemcells.txt
txtfile<-  file.path(getwd(),"R", "Stochastic Lab Course", "stemcell.txt")
Stemcells<- read_table(txtfile, col_names = "number")

Spline.Regression <- function(Y, X, k, deg){
  # Y: response
  # X: covariate
  # k: number of knots
  # deg: degree of spline
  
  m <- deg + 1 
   
  knots = seq(min(X), max(X), length.out = k+2)
  knots = c(min(X) - (m-1):1 , knots, max(X) + 1:(m-1))
  
 
  N <- spline.des(knots, X, ord = m)$design
  
  
  beta <- lm(Y ~ N - 1)$coefficients 
  
  
  f.hat <- function(x){
    n <- spline.des(knots, x, ord = m,  outer.ok = TRUE)$design
    as.numeric(n %*% beta)
  }
  
  return(f.hat)
}
 #Estimate f with different splines and a degree of 2
Y <- Stemcells$number
X <- seq(10, 1440, 10)

Fit.EstmateKnot1 <- Spline.Regression(Y, X, 3, deg = 2)
Fit.EstmateKnot2 <- Spline.Regression(Y, X, 8, deg = 2)
Fit.EstmateKnot3 <- Spline.Regression(Y, X, 13, deg = 2)
Fit.EstmateKnot4 <- Spline.Regression(Y, X, 17, deg = 2)

Time = seq(10, 1440, 1)
Data.FrameKnot1 <- data.frame(x = Time, y = Fit.EstmateKnot1(Time))
Data.FrameKnot2 <- data.frame(x = Time, y = Fit.EstmateKnot2(Time))
Data.FrameKnot3 <- data.frame(x = Time, y = Fit.EstmateKnot3(Time))
Data.FrameKnot4 <- data.frame(x = Time, y = Fit.EstmateKnot4(Time))

#plot
ggplot(data = data.frame(x = X, y = Y), aes(x = x, y = y)) +
  geom_point(alpha = 0.3) +
  geom_line(data = Data.FrameKnot1, aes(color = "3")) +
  geom_line(data = Data.FrameKnot2, aes(color = "8")) +
  geom_line(data = Data.FrameKnot3, aes(color = "13")) +
  geom_line(data = Data.FrameKnot4, aes(color = "17")) +
  scale_colour_manual(name = "Number of \nKnots", values = c("red", "blue", "yellow", "green")) +
  xlab("Time in minutes") +
  ylab("Order Parameter") +
  labs(title = "Estimate with different number of knots and degree of 2") +
  theme_gray(base_size = 12) + #font size
  theme(legend.key = element_rect(fill = "white", colour = "black")) 

#Estimate f using a degree from 1 to 4 and number of knots to 4
Fit.EstmateDeg1 <- Spline.Regression(Y, X, 4, deg = 1)
Fit.EstmateDeg2 <- Spline.Regression(Y, X, 4, deg = 2)
Fit.EstmateDeg3 <- Spline.Regression(Y, X, 4, deg = 3)
Fit.EstmateDeg4 <- Spline.Regression(Y, X, 4, deg = 4)

#data frame
Data.FrameDeg1 <- data.frame(x = Time, y = Fit.EstmateDeg1(Time))
Data.FrameDeg2 <- data.frame(x = Time, y = Fit.EstmateDeg2(Time))
Data.FrameDeg3 <- data.frame(x = Time, y = Fit.EstmateDeg3(Time))
Data.FrameDeg4 <- data.frame(x = Time, y = Fit.EstmateDeg4(Time))

#plot
ggplot(data = data.frame(x = X, y = Y), aes(x = x, y = y)) +
  geom_point(alpha = 0.3) +
  geom_line(data = Data.FrameDeg1, aes(color = "1")) +
  geom_line(data = Data.FrameDeg2, aes(color = "2")) +
  geom_line(data = Data.FrameDeg3, aes(color = "3")) +
  geom_line(data = Data.FrameDeg4, aes(color = "4")) +
  scale_colour_manual(name = "Degrees", values = c("magenta", "purple", "yellow", "green")) +
  xlab("Time in minutes") + ylab("Order Parameter") + 
  labs(title = "Estimate with number of knots as 4 and degree from 1 to 4")+
  theme_gray(base_size = 12) +
  theme(legend.key = element_rect(fill = "white", colour = "black")) 
  
  # function that estimates the optimal number of equidistant knots using Generalised Cross Validation (GCV)
  GCV <- function(k, deg){
    Spline.fit <- Spline.Regression(Y, X, k, deg = deg)
    Spline.fit <- Spline.fit(X)
    
    a <- norm(Y - Spline.fit, type = "2")
    n <- length(Y)
    
    res <- (a**2)/(1 - k/n)**2
    
    return(res)
  }
  #Vectorize GCV
  GCV <- Vectorize(GCV, vectorize.args = c("k"))  

  #Optimize over k for degrees 1 to 4
  max.k <- 50
  
  Optk.deg1 <- which(GCV(1:max.k, deg = 1) == min(GCV(1:max.k, deg = 1)))
  Optk.deg2 <- which(GCV(1:max.k, deg = 2) == min(GCV(1:max.k, deg = 2)))
  Optk.deg3 <- which(GCV(1:max.k, deg = 3) == min(GCV(1:max.k, deg = 3)))
  Optk.deg4 <- which(GCV(1:max.k, deg = 4) == min(GCV(1:max.k, deg = 4)))
  
  DataFrame.Optk <- data.frame(deg1 = Optk.deg1, deg2 = Optk.deg2,
                         deg3 = Optk.deg3, deg4 = Optk.deg4)
  row.names(DataFrame.Optk) <- "Optimal GCV knot number"
  DataFrame.Optk
  
  # Calculate fits with GCV knots number, for degrees 1 to 4
  Fit.GCV1 <- Spline.Regression(Y, X, k = Optk.deg1, deg = 1)
  Fit.GCV2 <- Spline.Regression(Y, X, k = Optk.deg2, deg = 2)
  Fit.GCV3 <- Spline.Regression(Y, X, k = Optk.deg3, deg = 3)
  Fit.GCV4 <- Spline.Regression(Y, X, k = Optk.deg4, deg = 4)
  
  
  #Dataframes
  Data.FrameGCV1 <- data.frame(x = Time, y =  Fit.GCV1(Time))
  Data.FrameGCV2 <- data.frame(x = Time, y =  Fit.GCV2(Time))
  Data.FrameGCV3 <- data.frame(x = Time, y =  Fit.GCV3(Time))
  Data.FrameGCV4 <- data.frame(x = Time, y =  Fit.GCV4(Time))
  
  #plot
   ggplot(data = data.frame(x = X, y = Y), aes(x = x, y = y)) +
    geom_point(alpha = 0.3) +
    geom_line(data = Data.FrameGCV1, aes(color = "1")) +
    geom_line(data = Data.FrameGCV2, aes(color = "2")) +
    geom_line(data = Data.FrameGCV3, aes(color = "3")) +
    geom_line(data = Data.FrameGCV4, aes(color = "4")) +
    scale_colour_manual(name = "Degree",  values = c("yellow", "blue", "red", "green")) +
    xlab("Time in minutes") +
    ylab("Order Parameter") +  labs(title = "") +
    theme_gray(base_size = 12) +
    theme(legend.key = element_rect(fill = "white", colour = "gray19"))
 
   #6c
  #Update function for regression splines so that it takes into account that 
  #the errors follow an autoregressive process of order one
  
  Spline.Autoregressive <- function(Y, X, k, deg){
    #Y: response variable
    #X: covariate
    #k: number of knots
    #deg: degree of spline
    
    m = deg + 1 
    
    # Determine knot points 
    knots = seq(min(X), max(X), length.out = k+2)
    knots = c(min(X) - (m-1):1 , knots, max(X) + 1:(m-1)) 
    
    N = spline.des(knots, X, ord = m)$design
    
    a = 0.55
    n = length(Y)
    R = toeplitz(a**(0:(n-1)))
    R.inverse = solve(R)
    NRN = t(N) %*% R.inverse %*% N
    NRY = t(N) %*% R.inverse %*% Y
    M = solve(NRN) %*% NRY
    
    #Regression spline
    f.hat = function(x){
    
      n = spline.des(knots, x, ord = m,  outer.ok = TRUE)$design
      as.numeric(n %*% M)
    }
    return(f.hat)
  }
  
  #Update function for GCV so that it takes into account that the errors follow an autoregressive process of order one
  
  GCV.Autoregressive <- function(k, deg){
    Spline.fit = Spline.Autoregressive(Y, X, k, deg = deg)
    Spline.fit = Spline.fit(X)
    
    a = 0.55
    n = length(Y)
    R = toeplitz(a**(0:(n-1)))
    R.inverse <- solve(R)
    
    aux = t(Y - Spline.fit) %*% R.inverse %*% (Y - Spline.fit)
    n = length(Y)
    
    res = aux/(1 - k/n)**2
    
    return(res)
  }
  
  GCV.Autoregressive <- Vectorize(GCV.Autoregressive, vectorize.args = c("k"))
  
  # Optimise over k for degrees 1 to 4 for the updated model
  Optk.Autoregressive1 <- which(GCV.Autoregressive(1:max.k, deg = 1) == min(GCV.Autoregressive(1:max.k, deg = 1)))
  Optk.Autoregressive2 <- which(GCV.Autoregressive(1:max.k, deg = 2) == min(GCV.Autoregressive(1:max.k, deg = 2)))
  Optk.Autoregressive3 <- which(GCV.Autoregressive(1:max.k, deg = 3) == min(GCV.Autoregressive(1:max.k, deg = 3)))
  Optk.Autoregressive4 <- which(GCV.Autoregressive(1:max.k, deg = 4) == min(GCV.Autoregressive(1:max.k, deg = 4)))
  
  
  DFOptk.Autoregressive <- data.frame(deg1 = Optk.Autoregressive1, deg2 = Optk.Autoregressive2,
                              deg3 = Optk.Autoregressive3, deg4 = Optk.Autoregressive4)
  row.names(DFOptk.Autoregressive) <- "Optimal GCV knot number (autoregressive)"
  DFOptk.Autoregressive
  
  # Calculating estimators with knot numbers from the updated GCV for degrees 1 to 4
  Fit.Autoregressive1 <- Spline.Autoregressive(Y, X, k = Optk.Autoregressive1, deg = 1)
  Fit.Autoregressive2 <- Spline.Autoregressive(Y, X, k = Optk.Autoregressive2, deg = 2)
  Fit.Autoregressive3 <- Spline.Autoregressive(Y, X, k = Optk.Autoregressive3, deg = 3)
  Fit.Autoregressive4 <- Spline.Autoregressive(Y, X, k = Optk.Autoregressive4, deg = 4)
  
 
  #dataframes for ggplot
  Data.FrameAutoregressive1 <- data.frame(x = Time, y = Fit.Autoregressive1(Time))
  Data.FrameAutoregressive2 <- data.frame(x = Time, y = Fit.Autoregressive2(Time))
  Data.FrameAutoregressive3 <- data.frame(x = Time, y = Fit.Autoregressive3(Time))
  Data.FrameAutoregressive4 <- data.frame(x = Time, y = Fit.Autoregressive4(Time))
  
  #Plot
  ggplot(data = data.frame(x = X, y = Y), aes(x = x, y = y)) +
    geom_point(alpha = 0.3) +
    geom_line(data = Data.FrameAutoregressive1, aes(color = "1")) +
    geom_line(data = Data.FrameAutoregressive2, aes(color = "2")) +
    geom_line(data = Data.FrameAutoregressive3, aes(color = "3")) +
    geom_line(data = Data.FrameAutoregressive4, aes(color = "4")) +
    scale_colour_manual(name = "Degree", values = c("yellow", "blue", "red", "green")) +
    xlab("Time in minutes") +
    ylab("Order Parameter") +
    theme_gray(base_size = 12) +
    theme(legend.key = element_rect(fill = "white", colour = "black"))  
  
  # Fitting parametric model of 4th degree
  FourthDeg.model <- gls(Y ~ X + I(X**2) + I(X**3) + I(X**4), correlation = corAR1(0.55))
  
  Polynomial.fit4 <- function(x){ as.numeric(FourthDeg.model$coefficients %*% x**(0:4)) }
  Polynomial.fit4 <- Vectorize(Polynomial.fit4)
  
  ggplot(data = data.frame(x = X, y = Y), aes(x = x, y = y)) +
    geom_point(alpha = 0.3) +
    geom_line(data = Data.FrameAutoregressive1, aes(color = "Degree 1")) +
    geom_line(data = Data.FrameAutoregressive2, aes(color = "Degree 2")) +
    geom_line(data = Data.FrameAutoregressive3, aes(color = "Degree 3")) +
    geom_line(data = Data.FrameAutoregressive4, aes(color = "Degree 4")) +
    stat_function(fun = Polynomial.fit4, aes(color = "Parametric Fit \n(degree 4)")) +
    scale_colour_manual(name = "", values = c("yellow", "blue", "black", "green", "red")) +
    xlab("Time in minutes") +
    ylab("Order Parameter") +
    theme_gray(base_size = 12) + #font size
    theme(legend.key = element_rect(fill = "white", colour = "black"))  #legend keys editing
  
  # For 3rd and 5th degree 
  ThirdDeg.model <- gls(Y ~ X + I(X**2) + I(X**3), correlation = corAR1(0.55))
  FifthDeg.model <- gls(Y ~ X + I(X**2) + I(X**3) + I(X**4) + I(X**5), correlation = corAR1(0.55))
  
  Polynomial.fit3 <- function(x){ as.numeric(ThirdDeg.model$coefficients %*% x**(0:3)) }
  Polynomial.fit5 <- function(x){ as.numeric(FifthDeg.model$coefficients %*% x**(0:5)) }
  Polynomial.fit3 <- Vectorize(Polynomial.fit3)
  Polynomial.fit5 <- Vectorize(Polynomial.fit5)
  
  ggplot(data = data.frame(x = X, y = Y), aes(x = x, y = y)) +
    geom_point(alpha = 0.3) +
    stat_function(fun = Polynomial.fit4, aes(color = "4")) +
    stat_function(fun = Polynomial.fit3, aes(color = "3")) +
    stat_function(fun = Polynomial.fit5, aes(color = "5")) + 
    scale_colour_manual(name = "Parametric fit \ndegree ", values = c("green", "red","blue")) +
    xlab("Time in minutes") +
    ylab("Order Parameter") +
    theme_gray(base_size = 12) +
    theme(legend.key = element_rect(fill = "white", colour = "black"))
  
  
  
  
