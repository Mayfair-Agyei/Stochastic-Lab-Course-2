library(tidyverse)
install.packages("fitdistrplus")
library(fitdistrplus)
csvfile <- file.path(getwd(),"R", "Stochastic Lab Course", "student-mat.csv")
studentMat <- read.csv(csvfile)
G1 <- data.frame(studentMat$G1, rep('G1', nrow(studentMat)))
colnames(G1) <- c("grades", "types")
G2 <- data.frame(studentMat$G2, rep('G2', nrow(studentMat)))
colnames(G2) <- c("grades", "types")
G3 <- data.frame(studentMat$G3, rep('G3', nrow(studentMat)))
colnames(G3) <- c("grades", "types")
All.groups <- rbind(G1, G2, G3)
#Test for Normal distribution
ggplot(data = All.groups, mapping = aes(sample = grades)) + 
  geom_density(aes(x = grades), fill = "light blue") +
  facet_grid(. ~types) 
labs(title = "Density plot of student grades",caption = "Grades data of students for mathematics course"
     ,x = "Grades",y = "Density") + theme_gray()

ggplot(data = All.groups, mapping = aes(sample = grades)) + 
  stat_qq(distribution = stats::qnorm, dparams = list(mean = mean(All.groups$grades), sd = sd(All.groups$grades))) +
  geom_abline(alpha = 0.25) +
  facet_grid(. ~types) + 
  labs(title = "QQ-Plot plot of student grades",caption = "Grades data of students for mathematics course"
       ,x = "Theoretical",y = "Sample") + theme_gray()

shapiro.test(studentMat$G1)
shapiro.test(studentMat$G2)
shapiro.test(studentMat$G3)

#Test for a Poisson distribution
ggplot(data = All.groups, mapping = aes(sample = grades)) + 
  stat_qq(distribution = stats::qpois, dparams = list(lambda = mean(All.groups$grades))) +
  geom_abline(alpha = 0.25) +
  facet_grid(. ~types) + labs(title = "QQ-Plot plot of student grades in Mathematics",caption = "Grades data of students for mathematics course"
                              ,x = "Theoretical",y = "Sample") + theme_gray()
par(mfrow=c(1,1))
n=length(All.groups$grades)
(x=table(All.groups$grades))
k=as.numeric(names(x))
plot(k,log(x)+lfactorial(k))                                                                



#
Model1 <- glm(G1~.-G2-G3, family=poisson, data=studentMat)
summary(Model1)
# Pearson's Residual for Model1
PearsonResidual <- residuals(Model1, "pearson")
PearsonResidual <- data.frame(PearsonResidual)

# Visualization for Pearson Residual for Model1
ggplot(data = PearsonResidual, mapping = aes(sample = PearsonResidual)) + 
  stat_qq(distribution = stats::qnorm, dparams = list(mean = mean(PearsonResidual$PearsonResidual),
                                                      sd = sd(PearsonResidual$PearsonResidual))) +
  geom_abline(alpha = 0.25) + labs(title = "Plot of Pearson Residuals"
                                   ,x = "Theoretical",y = "Sample") + theme_gray()

# Function for Anscombe Residual
anscombe <- function(y, mu){
  (3 *(y ** (2/3)- mu ** (2/3)))/2*(mu ** (1/6))
}

# Anscombe Residual for Model1
AnscombeResidual <- anscombe(studentMat$G1, Model1$fitted.values)
AnscombeResidual <- data.frame(AnscombeResidual)

# Visualization for Anscombe Residual for Model1
ggplot(data = AnscombeResidual, mapping = aes(sample = AnscombeResidual)) + 
  stat_qq(distribution = stats::qnorm, dparams = list(mean = mean(AnscombeResidual$AnscombeResidual), 
                                                      sd = sd(AnscombeResidual$AnscombeResidual))) +
  geom_abline(alpha = 0.25) + labs(title = "Plot of Anscombe Residuals"
                                   ,x = "Theoretical",y = "Sample") + theme_gray()

#Residual Analysis for Model1
par(mfrow = c(2, 2))
plot(Model1)

#Model2
Model2 <- glm(formula = G1 ~ sex + Fedu + studytime + failures + schoolsup + famsup + goout , 
               family = poisson, data = studentMat) 
summary(Model2)

#Analysis of Deviance
anova(Model1, Model2, test = "Chisq")

#Model3
Model3 <- glm(formula = G1 ~ sex + Fedu + studytime + failures + schoolsup + famsup + Walc , family = poisson, data = studentMat) 
summary(Model3)

#Residual Analysis for Model3
par(mfrow = c(2, 2))
plot(Model3)


