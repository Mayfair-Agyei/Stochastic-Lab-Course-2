library(tidyverse)
library(haven)
library(ggplot2)
install.packages("NonpModelCheck")
library(NonpModelCheck)
library(dplyr)


#Childrenfinal3 <- read.delim("R/Stochastic Lab Course/ChildrenfinalTib.txt", sep = "\t") not separating for some reason 
#so I had to change my format for reading the file into R.

# reading the data children.dta using read_dta command
dtafile <-  file.path(getwd(), "R","Stochastic Lab Course","childrenfinal.dta")
childrenfinal <- read_dta(dtafile)

childrenfinal_new <- childrenfinal %>% select(-matches("^[svm][0-9]"))

childrenfinal_new1<- childrenfinal_new %>% mutate_if(is.double, as.double)

Childrenfinal3 <- select(childrenfinal_new1, hypage,ruralfacto,female,zstunt,zweight,zwast,adm2)


#Kernel Functions
uniform <- function(x) (abs(x)< 1) * 0.5
gaussian <- function(x)1/sqrt(2*pi)* exp(-(x^2)/2)
triangular <- function(x)(abs(x)< 1)* (1 -abs(x))
epanechnikov <- function(x)(abs(x)< 1) *(0.75 * (1 - x^2))

KerNames <- list("uniform", "triangular", "epanechnikov", "gaussian")
KerFunc <- list(uniform,triangular, epanechnikov, gaussian)                                                 
                                      


Local.Polynomial<- function(Y, X, bw, deg, kernel = "epanechnikov"){
  # Y = response variable
  # X = covariate
  # bw = bandwidth
  # deg = polynomial degree
  # kernel = epanechnikov 
  
  kernel <- KerFunc[[which(KerNames == kernel)]]
  Y <- as.matrix(Y)
  X <- as.matrix(X)
  d <- ncol(X)
  n <- nrow(X)
  X_mat <- matrix(0, d*n, deg+1)
  
  A.hat <- function(x, derivative){
    # derivative = number of derivatives
    equation1 <- (X - x)/bw
    equation2 <- as.vector(t(X - x))
    
    V <- apply(equation1, MARGIN = 1, FUN = kernel)
    
    for (i in 0:deg){
      X_mat[ ,i+1] <- equation2^i
    }
    res <- lm(Y ~ X_mat - 1, weights = V) # weighted least squares (exclude intercept)
    res <- res$coefficients
    res <- factorial(derivative)*res[derivative + 1]
    return(res)
  }
  
  VA.hat <- Vectorize(A.hat, vectorize.args = "x")
  return(VA.hat)
}


  # Estimate of f with 4 different bandwidth and polynomial degree of 1
Fit.Estimate1 <-Local.Polynomial(Childrenfinal3$zwast, Childrenfinal3$hypage, bw = 2, deg = 1)
Fit.Estimate2 <-Local.Polynomial(Childrenfinal3$zwast, Childrenfinal3$hypage, bw = 9, deg = 1)
Fit.Estimate3 <-Local.Polynomial(Childrenfinal3$zwast, Childrenfinal3$hypage, bw = 15, deg = 1)
Fit.Estimate4 <-Local.Polynomial(Childrenfinal3$zwast, Childrenfinal3$hypage, bw = 24, deg = 1)
  
  # Fit of Domain of hypage
  mininum = min(Childrenfinal3$hypage)
  maximum = max(Childrenfinal3$hypage)
  
 
Domain.Fit1 <- Fit.Estimate1(x = 0:maximum, derivative = 0)
Domain.Fit2 <- Fit.Estimate2(x = 0:maximum, derivative = 0)
Domain.Fit3 <- Fit.Estimate3(x = 0:maximum, derivative = 0)
Domain.Fit4 <- Fit.Estimate4(x = 0:maximum, derivative = 0)
  
  #dataframes for ggplot
  Data.Frame1 <- data.frame(a = 0:maximum, b = Domain.Fit1)
  Data.Frame2 <- data.frame(a = 0:maximum, b = Domain.Fit2)
  Data.Frame3 <- data.frame(a = 0:maximum, b = Domain.Fit3)
  Data.Frame4 <- data.frame(a = 0:maximum, b = Domain.Fit4)
  
 ggplot(Childrenfinal3, aes(x = hypage, y = zwast)) +
    geom_point(color = "black") +
    geom_line(data = Data.Frame1, aes(x = a, y = b, color = "2"), size = 1) +
    geom_line(data = Data.Frame2, aes(x = a, y = b, color = "9"), size = 1) +
    geom_line(data = Data.Frame3, aes(x = a, y = b, color = "15"), size = 1) +
    geom_line(data = Data.Frame4, aes(x = a, y = b, color = "24"), linetype = "dashed", size = 1) +
    labs(color = "Bandwidth") +
     labs(title = "",x = "hypage",y = "zwast", colour = "BANDWIDTH") + theme_gray()
 
 #Estimate of f with the 4 kernels from exercise 4 to a degree of 1
 Fit.Estimate5 <- Local.Polynomial(Childrenfinal3$zwast, Childrenfinal3$hypage, bw = 15, deg = 1, kernel = "epanechnikov")
 Fit.Estimate6 <- Local.Polynomial(Childrenfinal3$zwast, Childrenfinal3$hypage, bw = 15, deg = 1, kernel = "uniform")
 Fit.Estimate7 <- Local.Polynomial(Childrenfinal3$zwast, Childrenfinal3$hypage, bw = 15, deg = 1, kernel = "triangular")
 Fit.Estimate8 <- Local.Polynomial(Childrenfinal3$zwast, Childrenfinal3$hypage, bw = 15, deg = 1, kernel = "gaussian")
 
 #Fitting domain of hypage = 0:59
 Domain.Fit5 <- Fit.Estimate5(x = 0:59, derivative = 0)
 Domain.Fit6 <- Fit.Estimate6(x = 0:59, derivative = 0)
 Domain.Fit7 <- Fit.Estimate7(x = 0:59, derivative = 0)
 Domain.Fit8 <- Fit.Estimate8(x = 0:59, derivative = 0)
 
 #dataframes for plot
 Data.Frame5 <- data.frame(a = 0:59, b = Domain.Fit5)
 Data.Frame6 <- data.frame(a = 0:59, b = Domain.Fit6)
 Data.Frame7 <- data.frame(a = 0:59, b = Domain.Fit7)
 Data.Frame8 <- data.frame(a = 0:59, b = Domain.Fit8)
 
 ggplot(Childrenfinal3, aes(x = hypage, y = zwast)) +
   geom_point(color = "black") +
   geom_line(data = Data.Frame5, aes(x = a, y = b, color = "epanechnikov"), size = 1) +
   geom_line(data = Data.Frame6, aes(x = a, y = b, color = "uniform"), size = 1) +
   geom_line(data = Data.Frame7, aes(x = a, y = b, color = "triangular"), size = 1) +
   geom_line(data = Data.Frame8, aes(x = a, y = b, color = "gaussian"), linetype = "dashed", size = 1) +
   labs(color = "Kernel") + labs(title = "",x = "hypage",y = "zwast", colour = "KERNEL") + theme_gray()
 

  #function to calculate the optimal bandwidth using Generalised Cross Validation (GCV)
 
 GCV <- function(Y, X, bw, deg){ 
   # Y,X,bw,deg are defined as in Local Polynomial
   X_values <- unique(sort(X))
   n <- length(Y)
   
   X_mat_list <- list()
   for (i in X_values){
     index = which(X_values == i)
     X_mat_list[[index]] = matrix(0, n, 4+1)
     aux = (X - i)
     for (j in 0:4){
       X_mat_list[[index]][, j+1] = aux**j
     }
   }
   sum_square = rep(0, length(X_values))
   W_trace = matrix(0, length(X_values)) 
   
   
   Polynomial.fit = Local.Polynomial(Y, X, bw, deg)
   
  
   for (i in X_values){
     index = which(X_values == i)
     aux = (X - i)
     
     X_mat = X_mat_list[[index]][, 1:(deg+1)]
     
     V = diag(epanechnikov(aux/bw))
     
     weight_vector = solve(t(X_mat) %*% V %*% X_mat) %*% t(X_mat) %*% V
     weight_vector = weight_vector[1, ]
     
     W_trace[index] = sum((X == i) * weight_vector)
     
     sum_square[index] = sum((Y[(X == i)] - Polynomial.fit(i, derivative = 0))**2)
   }
   res = sum(sum_square)/(1 - sum(W_trace)/n)**2
   return(res)
   
 }
 
 # Optimal bandwidth with GCV to estimate f using polynomial degree 1 to 4.
 Y <- Childrenfinal3$zwast
 X <- Childrenfinal3$hypage
 
 GCV1 <- Vectorize(function(bw){GCV(Y, X, bw, deg =  1)}) 
 GCV2 <- Vectorize(function(bw){GCV(Y, X, bw, deg =  2)})
 GCV3 <- Vectorize(function(bw){GCV(Y, X, bw, deg =  3)})
 GCV4 <- Vectorize(function(bw){GCV(Y, X, bw, deg =  4)})
 
 GCV.Opt1 <- optimize(GCV1, interval = c(2, 11))$minimum
 GCV.Opt2 <- optimize(GCV1, interval = c(2, 11))$minimum
 GCV.Opt3 <- optimize(GCV3, interval = c(2, 11))$minimum
 GCV.Opt4 <- optimize(GCV4, interval = c(2, 11))$minimum
 

 Fit.Estimate9 <- Local.Polynomial(childrenfinal3$zwast, childrenfinal3$hypage, bw = GCV.Opt1, deg = 1, kernel = "epanechnikov")
 Fit.Estimate10 <- Local.Polynomial(childrenfinal3$zwast, childrenfinal3$hypage, bw = GCV.Opt2, deg = 2, kernel = "epanechnikov")
 Fit.Estimate11 <- Local.Polynomial(childrenfinal3$zwast, childrenfinal3$hypage, bw = GCV.Opt3, deg = 3, kernel = "epanechnikov")
 Fit.Estimate12 <- Local.Polynomial(childrenfinal3$zwast, childrenfinal3$hypage, bw = GCV.Opt4, deg = 4, kernel = "epanechnikov")
 
 #Fitting on the domain of hypage = 0:59
 Domain.Fit9 <- Fit.Estimate9(x = 0:maximum, derivative = 0)
 Domain.Fit10 <- Fit.Estimate10(x = 0:maximum, derivative = 0)
 Domain.Fit11 <- Fit.Estimate11(x = 0:maximum, derivative = 0)
 Domain.Fit12 <- Fit.Estimate12(x = 0:maximum, derivative = 0)
 
 #dataframes for ggplot
 Data.Frame9 <- data.frame(a = 0:maximum, b = Domain.Fit9)
 Data.Frame10 <- data.frame(a = 0:maximum, b = Domain.Fit10)
 Data.Frame11 <- data.frame(a = 0:maximum, b = Domain.Fit11)
 Data.Frame12 <- data.frame(a = 0:maximum, b = Domain.Fit12)
 
 #Then the plot
 ggplot(children_subset, aes(x = hypage, y = zwast)) +
   geom_point(color = "gray31") +
   geom_line(data = Domain.Fit9, aes(x = a, y = b, color = "1")) +
   geom_line(data = Domain.Fit10, aes(x = a, y = b, color = "2")) +
   geom_line(data = Domain.Fit11, aes(x = a, y = b, color = "3")) +
   geom_line(data = Domain.Fit12, aes(x = a, y = b, color = "4")) +
   labs(color = "polynomial \n degrees")
 
 Derivative1 <- localpoly.reg(X, Y, bandwidth = 4.99996, degree.pol = 1, deriv = 1)
 Derivative2 <- localpoly.reg(X, Y, bandwidth = 10.9999, degree.pol = 2, deriv = 1)
 Derivative3 <- localpoly.reg(X, Y, bandwidth = 10.9999, degree.pol = 3, deriv = 1)
 Derivative4 <- localpoly.reg(X, Y, bandwidth = 8.40399, degree.pol = 4, deriv = 1)
 #dataframes for ggplot
 Data.Frame13 <- data.frame(a = unique(Derivative1$x), b = unique(Derivative1$predict))
 Data.Frame14 <- data.frame(a = unique(Derivative2$x), b = unique(Derivative2$predict))
 Data.Frame15 <- data.frame(a = unique(Derivative3$x), b = unique(Derivative3$predict))
 Data.Frame16 <- data.frame(a = unique(Derivative4$x), b = unique(Derivative4$predict))
 
 #Plot all four derivative fits on one plot
  ggplot(Childrenfinal3, aes(x = hypage, y = zwast)) +
   geom_point(color = "blue") +
   geom_line(data = Data.Frame13, aes(x = a, y = b, color = "1")) +
   geom_line(data = Data.Frame14, aes(x = a, y = b, color = "2")) +
   geom_line(data = Data.Frame15, aes(x = a, y = b, color = "3")) +
   geom_line(data = Data.Frame16, aes(x = a, y = b, color = "4")) +
   labs(color = "polynomial degrees")
 
 
 
   
   