library(tidyverse)
library(readr)
library(ggplot2)

StudentsPerformance <- read_csv("R/Stochastic Lab Course/StudentsPerformance.csv")
math.score <-StudentsPerformance$`math score`
reading.score <-StudentsPerformance$`reading score`
writing.score <-StudentsPerformance$`writing score`

#4a 

# Kernel Functions
uniform <- function(x) (abs(x)< 1) * 0.5
gaussian <- function(x)1/sqrt(2*pi)* exp(-(x^2)/2)
triangular <- function(x)(abs(x)< 1)* (1 -abs(x))
epanechnikov <- function(x)(abs(x)< 1) *(0.75 * (1 - x^2))

KerNames <- list("uniform", "triangular", "epanechnikov", "gaussian")
KerFunc <- list(uniform,triangular, epanechnikov, gaussian)

#Kernel Density Estimator
KerDesEst <- function(sample, bandwith, kernel){
  K<- KerFunc[[which(KerNames == kernel)]]
  a <- function(x){
    t <- (x - sample)/bandwith
    (1/(bandwith*length(sample)))*(sum(K(t)))
  }
  return(Vectorize(a))
}

# Plot of Epanechnikov kernel with different bandwidths
ggplot(StudentsPerformance, aes(x =  `math score`)) + 
  stat_function(fun = KerDesEst(math.score, 2,"epanechnikov"),aes (colour = "2")) + 
  stat_function(fun = KerDesEst(math.score, 9,"epanechnikov"),aes (colour = "9")) + 
  stat_function(fun = KerDesEst(math.score, 15,"epanechnikov"),aes (colour = "15")) + 
  stat_function(fun = KerDesEst(math.score, 24,"epanechnikov"),aes (colour = "24")) +
  theme(legend.position = c (0.96, 0.95), legend.justification = c("right", "bottom"), legend.title = element_text()) +
  labs(colour = "Bandwidth") + labs(title = "Epanechnikov Kernel with different bandwidths (Kernel Density Estimators)",
                                    x = "Mathematics Score", y = "Density", colour = "Bandwidth") + theme_gray()

# Plot of different Kernel functions with bandwidth of 15
ggplot(StudentsPerformance, aes(x =  `math score`)) + 
  stat_function(fun = KerDesEst(math.score, 15,"uniform"),aes (colour = "uniform")) + 
  stat_function(fun = KerDesEst(math.score, 15,"epanechnikov"),aes (colour = "epanechnikov")) + 
  stat_function(fun = KerDesEst(math.score, 15,"gaussian"),aes (colour = "gaussian")) + 
  stat_function(fun = KerDesEst(math.score, 15,"triangular"),aes (colour = "triangular")) +
  theme(legend.position = c (0.96, 0.95), legend.justification = c("right", "bottom"), legend.title = element_text()) +
  labs(colour = "Kernel") + labs(title = "Different Kernel functions with bandwidth of 15 (Kernel Density Estimators)",
                                    x = "Mathematics Score", y = "Density", colour = "Kernel Function") + theme_gray()

#4b

# Cross Validation Criterion
# Cross Validation Function
Cross_Validation <- function(sample, bandwidth){
  n <- length(sample)
  
  f_Hat <- KDE(sample, bandwidth, "gaussian")
  
  f_Hat.sqr <- function(x){f_Hat(x)**2} 
  
  A <- integrate(f_Hat.sqr, lower = Inf, upper = Inf) 
  
  B <- (sum(gaussian(outer(sample, sample, FUN = "-")/bandwidth)) - n * gaussian (0)) * 2/(n *(n - 1)* bandwidth)
  
  C <- A$value - B
  
  return(C)
}
CV.math <- Vectorize(function(bandwidth){Cross_Validation(math.score, bandwidth)})
CV.reading <- Vectorize(function(bandwidth){Cross_Validation(reading.score, bandwidth)})
CV.writing <-Vectorize(function(bandwidth){Cross_Validation(writing.score, bandwidth)})

 # Optimising the objective function
CVMath.Opt <- optimize(CV.math, interval = c(1, 25))$minimum
CVReading.Opt <- optimize(CV.reading, interval = c(1,15))$minimum
CVWriting.opt <- optimize(CV.writing, interval = c(1, 50))$minimum
 
# built-in bw.ucv
bw.ucv.Math.opt <- bw.ucv(math.score, lower = 1, upper = 25)
bw.ucv.Reading.opt <- bw.ucv(reading.score, lower = 1, upper = 15)
bw.ucv.Writing.opt <- bw.ucv(writing.score, lower = 1, upper = 50)

#built-in bw.bcv
bw.bcv.Math.opt <- bw.bcv(math.score, lower = 1, upper = 25)
bw.bcv.Reading.opt <- bw.bcv(reading.score, lower = 1, upper = 15)
bw.bcv.Writing.opt <- bw.bcv(writing.score, lower = 1, upper = 50)


#4c
# Students who took part in the  test preparation course
# Cross Validation
Students.Participation <- StudentsPerformance %>% dplyr::filter(`test preparation course` == "completed")

math.score.participated <- Students.Participation$`math score`
reading.score.participated <- Students.Participation$`reading score`
writing.score.participated <- Students.Participation$`writing score`

CVMath.participated <- Vectorize(function(bandwidth){Cross_Validation(math.score.participated, bandwidth)})
CVReading.participated <- Vectorize(function(bandwidth){Cross_Validation(reading.score.participated, bandwidth)})
CVWriting.participated <- Vectorize(function(bandwidth){Cross_Validation(writing.score.participated, bandwidth)})

CVMath.Opt.participated <- optimize(CVMath.participated, interval = c(1, 25))$minimum
CVReading.Opt.participated <- optimize(CVMath.participated, interval = c(1, 15))$minimum
CVWriting.Opt.participated <- optimize(CVMath.participated, interval = c(1, 50))$minimum
 
# Kernel Density Estimators
KDEMath.participated <- KerDesEst(Students.Participation$`math score`, CVMath.Opt.participated, "gaussian") 
KDEReading.participated <- KerDesEst(Students.Participation$`reading score`, CVReading.Opt.participated, "gaussian") 
KDEWriting.participated <- KerDesEst(Students.Participation$`writing score`, CVMath.Opt.participated, "gaussian") 

#Students who did not take part in the test preparation course
# Cross validation
Students.NoParticipation <- StudentsPerformance %>% dplyr::filter(`test preparation course` == "none")

math.score.Notparticipated <- Students.NoParticipation$`math score`
reading.score.Notparticipated <- Students.NoParticipation$`reading score`
writing.score.Notparticipated <- Students.NoParticipation$`writing score`

CVMath.Notparticipated <- Vectorize(function(bandwidth){Cross_Validation(math.score.Notparticipated, bandwidth)})
CVReading.Notparticipated <- Vectorize(function(bandwidth){Cross_Validation(reading.score.Notparticipated, bandwidth)})
CVWriting.Notparticipated <- Vectorize(function(bandwidth){Cross_Validation(writing.score.Notparticipated, bandwidth)})

CVMath.Opt.Notparticipated <- optimize(CVMath.Notparticipated, interval = c(1, 25))$minimum
CVReading.Opt.Notparticipated <- optimize(CVMath.Notparticipated, interval = c(1, 15))$minimum
CVWriting.Opt.Notparticipated <- optimize(CVMath.Notparticipated, interval = c(1, 50))$minimum


# Kernel Density Estimators
KDEMath.Notparticipated <- KerDesEst(Students.NoParticipation$`math score`, CVMath.Opt.Notparticipated, "gaussian") 
KDEReading.Notparticipated <- KerDesEst(Students.NoParticipation$`reading score`, CVReading.Opt.Notparticipated, "gaussian") 
KDEWriting.Notparticipated <- KerDesEst(Students.NoParticipation$`writing score`, CVMath.Opt.Notparticipated, "gaussian") 

# Plot comparison betwween students who took the test preparation course and students who did not.
ggplot(StudentsPerformance, aes(x = `math score`)) + 
  stat_function(fun = KDEMath.participated, aes(color = "completed")) +
  stat_function(fun = KDEMath.Notparticipated, aes(color = "none")) +
  scale_colour_manual(name="Test preparation course", values = c("blue", "green")) +
  theme(legend.position = c(0.96, 0.95),legend.justification = c("right", "bottom"), legend.title = element_text()) + 
  labs(title = "", x = "Math Score",y = "density", colour = "Preparatory Course") + theme_gray()

ggplot(StudentsPerformance, aes(x = `reading score`)) + 
  stat_function(fun = KDEReading.participated, aes(color = "completed")) +
  stat_function(fun = KDEReading.Notparticipated, aes(color = "none")) +
  scale_colour_manual(name="Test preparation course", values = c("blue", "green")) +
  theme(legend.position = c(0.96, 0.95),legend.justification = c("right", "bottom"), legend.title = element_text()) + 
  labs(title = "", x = "Reading Score",y = "density", colour = "Preparatory Course") + theme_gray()

ggplot(StudentsPerformance, aes(x = `writing score`)) + 
  stat_function(fun = KDEWriting.participated, aes(color = "completed")) +
  stat_function(fun = KDEWriting.Notparticipated, aes(color = "none")) +
  scale_colour_manual(name="Test preparation course", values = c("blue", "green")) +
  theme(legend.position = c(0.96, 0.95),legend.justification = c("right", "bottom"), legend.title = element_text()) + 
  labs(title = "", x = "Writing Score",y = "density", colour = "Preparatory Course") + theme_gray()

