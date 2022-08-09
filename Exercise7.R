library(tidyverse)
library(ggplot2)
library(Matrix)
install.packages("pls")
library(pls)


#7a
Users <- read_csv("R/Stochastic Lab Course/users.csv")
head(Users)
Likes <- read_csv("R/Stochastic Lab Course/likes.csv")
head(Likes)
Users.Likes <- read_csv("R/Stochastic Lab Course/users-likes.csv")
head(Users.Likes)


#Adding 2 extra columns
Users.row <- match(Users.Likes$userid, Users$userid)

Likes.row <- match(Users.Likes$likeid, Likes$likeid)

Users.Likes <- cbind(Users.Likes, Users.row, Likes.row)

#Building the Matrix
M <- sparseMatrix(i = Users.Likes$Users.row, j = Users.Likes$Likes.row, x = 1)
N <- M 

rownames(M) <- Users$userid
colnames(M) <- Likes$name

repeat {
  i <- sum(dim(M))
  M <- M[rowSums(M) >= 80, colSums(M) >= 150]
  if (sum(dim(M)) == i) break
}
Users <- Users[match(rownames(M),Users$userid), ]

Likes <- Likes[match(colnames(M),Likes$likeid), ]

Users.Likes.M <- as.matrix(M)

#7b
#Creating train and test sets from the data
set.seed(1122)

n <- nrow(Users.Likes.M)
Train.indices <- sample(1:n, size = round((2/3)*n ))
Train <- Users.Likes.M[Train.indices, ]
Test <- Users.Likes.M[-Train.indices, ]

Train.age <- Users$age[Train.indices]
Test.age <- Users$age[-Train.indices]


#fitting data in the regression
model <- plsr(Train.age ~ Train, ncomp = 50)

#prediction on the test set
Predict.age <- predict(model, newdata = Test)


corr <- 0
for (i in 1:50) {
  corr[i] <- cor(Predict.age[,1,i], Test.age, method = "pearson")
}

#Plot
 plot(corr)

#opt
d.opt <- which(corr == max(corr))
d.opt
#plot
Predicted.age <-Predict.age[ , ,d.opt]
Dataframe <- data.frame(Predicted.age, Test.age)

 ggplot(Dataframe, aes(x = Predicted.age, y =  Test.age)) +
  geom_point() +
  geom_abline( aes(intercept=0, slope=1,color='red'), size = 1) +
  scale_colour_manual(name = " ", values = c("red")) +
  theme_gray(base_size = 15) +
  expand_limits(x = 0, y = 0)


#7c
a <- tail(sort(model$coefficients[ ,1 , d.opt]), 6)
b <- head(sort(model$coefficients[ ,1 , d.opt]), 6)

#7d 
# Rerunning the analysis
Users <- read_csv("R/Stochastic Lab Course/users.csv")

Likes <- read_csv("R/Stochastic Lab Course/likes.csv")

Users.Likes <- read_csv("R/Stochastic Lab Course/users-likes.csv")


Users.row <- match(Users.Likes$userid, Users$userid)

Likes.row <- match(Users.Likes$likeid, Likes$likeid)

Users.Likes <- cbind(Users.Likes, Users.row, Likes.row)

#Matrix
N<- sparseMatrix(i = Users.Likes$Users.row, j = Users.Likes$Likes.row, x = 1)

rownames(N) <- Users$userid
colnames(N) <- Likes$name

repeat {
  i <- sum(dim(N))
  N <- N[rowSums(N) >= 60, colSums(N) >= 120]
  if (sum(dim(N)) == i) break
}

Users <- Users[match(rownames(N),Users$userid), ]

Likes <- Likes[match(colnames(N),Likes$likeid), ]

Users.Likes.N <- as.matrix(N)


set.seed(1122)
n <- nrow(Users.Likes.M)
Train.indices <- sample(1:n, size = round((2/3)*n ))
Train <- Users.Likes.M[Train.indices, ]
Test <- Users.Likes.M[-Train.indices, ]

Train.age2 <- Users$age[Train.indices]
Test.age2 <- Users$age[-Train.indices]


#fitting data in the regression
model2 <- plsr(Train.age ~ Train, ncomp = 50)

#prediction on the test set
Predict.age2 <- predict(model, newdata = Test)


corr <- 0
for (i in 1:50) {
  corr[i] <- cor(Predict.age2[,1,i], Test.age, method = "pearson")
}

#Plot
plot(corr)

d.opt <- which(corr == max(corr))
d.opt
Predicted.age2 <-Predict.age2[ , ,d.opt]
Dataframe <- data.frame(Predicted.age, Test.age)


