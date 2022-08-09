install.packages("JoSAE")
install.packages("conflicted")
library(tidyverse)
library(rgdal)
library(maptools)
library(Matrix)
library(mapdata)
library(nlme)
library(ggrepel)
library(maps)
library(nlme)
library(JoSAE)
library(conflicted)

#Loading the Data
data(landsat)
Satellite.data <- select(landsat, -outlier)


head(landsat)
tail(landsat)
summary(landsat)

#3a
#Creat groupedData
Corn.group <- groupedData(HACorn ~ PixelsCorn | CountyName, data = landsat)
Soybeans.group <-  groupedData(HASoybeans ~ PixelsSoybeans | CountyName, data = landsat)

# Linear Model for Corn and soya group
LinMod.corn <- lmList(Corn.group)
LinMod.corn
 plot(LinMod.corn)

LinMod.soybeans <- lmList(Soybeans.group)
LinMod.soybeans
 plot(LinMod.soybeans)


#3b

#Fitting linear mixed model for both crops so that segments share the same countywide random effect

# Linear Mixed Model for Corn and Soybeans group
LinMixMod.corn= lme(HACorn ~ PixelsCorn, data = Corn.group, random = ~ 1)
LinMixMod.corn
summary(LinMixMod.corn)
Corn.beta <- LinMixMod.corn$coefficients$fixed
Corn.beta
plot(LinMixMod.corn,HACorn~fitted(.)|CountyName,abline=c(0,1))

LinMixMod.soybeans <- lme(HASoybeans ~ PixelsSoybeans, data = Soybeans.group, random = ~ 1)
LinMixMod.soybeans
summary(LinMixMod.soybeans)
Soybeans.beta <- LinMixMod.soybeans$coefficients$fixed
Soybeans.beta
plot(LinMixMod.soybeans,HASoybeans~fitted(.)|CountyName,abline=c(0,1))


#3c

#Population for explanatory variables
Corn.popmean <- unique(Satellite.data$MeanPixelsCorn)
Soybeans.popmean <- unique(Satellite.data$MeanPixelsSoybeans)

# Mean of over observed segments only
Segments.mean = aggregate(Satellite.data[3:6], by = list(Satellite.data$CountyName), mean)


#Number of observations by county
Numobserved.county <- plyr::count(Satellite.data, "CountyName")
County.names <- Numobserved.county[, 1]
Num.observed <- Numobserved.county[, 2]

#Calculating parameters of Corn and Soybeans
Corn.varestimate <- VarCorr(LinMixMod.corn)
Soybeans.varestimate <- VarCorr(LinMixMod.soybeans)
Corn.sigmaestimate <- as.numeric(Corn.varestimate[2])
Soybeans.sigmaestimate <- as.numeric(Soybeans.varestimate[2])
Corn.randsigma <- as.numeric(Corn.varestimate[1])
Soybeans.randsigma <- as.numeric(Soybeans.varestimate[1])
 
#The covariance matrix V.hat of beta.hat
Corn.covmatrix<- list()
Soybeans.covmatrix <- list()
for (i in 1:12){
  Corn.covmatrix[[i]] = matrix(Corn.randsigma, Num.observed[i], Num.observed[i])
  Soybeans.covmatrix[[i]] = matrix(Soybeans.randsigma, Num.observed[i], Num.observed[i])
}

#Parameters for the estimates
Corn.idmatrix <-  diag(x = Corn.sigmaestimate, nrow = sum(Num.observed), ncol = sum(Num.observed))
soybeans.idmatrix <- diag(x = Soybeans.sigmaestimate, nrow = sum(Num.observed), ncol = sum(Num.observed))
Corn.BlockCVM <- bdiag(Corn.covmatrix)
Soybeans.BlockCVM <- bdiag(Soybeans.covmatrix)
V.Corn <- Corn.idmatrix + Corn.BlockCVM
V.Soybeans <- soybeans.idmatrix + Soybeans.BlockCVM

corn.a<- cbind(1, Satellite.data$PixelsCorn)
soybeans.a<- cbind(1, Satellite.data$PixelsSoybeans)
VHat.Corn <- solve(t(corn.a) %*% solve(V.Corn) %*% corn.a)
VHat.Soybeans <- solve(t(soybeans.a) %*% solve(V.Soybeans) %*% soybeans.a)

corn.gamma <- Corn.randsigma/(Corn.randsigma + Corn.sigmaestimate/Num.observed)
soybeans.gamma <- Soybeans.randsigma/(Soybeans.randsigma + Soybeans.sigmaestimate/Num.observed)

#Predictors

#regression predictor:
Corn.RegPred <- cbind(1, Corn.popmean) %*% Corn.beta
Soybeans.RegPred <- cbind(1, Soybeans.popmean) %*% Soybeans.beta

#adjusted survey predictor:
Corn.ASP <- cbind(1, Corn.popmean) %*% Corn.beta +
  (Segments.mean$HACorn - (cbind(1, Segments.mean$PixelsCorn) %*% Corn.beta))
Soybeans.ASP <- cbind(1, Soybeans.popmean) %*% Soybeans.beta +
  (Segments.mean$HASoybeans - (cbind(1, Segments.mean$PixelsSoybeans) %*% Soybeans.beta))

#Empirical BLUP (EBLUP)
corn.EBLUP <- cbind(1, Corn.popmean) %*% Corn.beta +
  corn.gamma*(Segments.mean$HACorn - (cbind(1, Segments.mean$PixelsCorn) %*% Corn.beta))
soybeans.EBLUP <- cbind(1, Soybeans.popmean) %*% Soybeans.beta +
  soybeans.gamma*(Segments.mean$HASoybeans - (cbind(1, Segments.mean$PixelsSoybeans) %*% Soybeans.beta))

#Survey predictor:
Corn.SP <- Segments.mean$HACorn
Soybeans.SP <- Segments.mean$HASoybeans

#Creating a dataframe for predictors
Dataframe.pred <- data.frame(County = County.names,Corn.RegPred,Soybeans.RegPred,
                      Corn.ASP,Soybeans.ASP,
                      corn.EBLUP,soybeans.EBLUP,
                      Corn.SP,Soybeans.SP)
view(Dataframe.pred)

#Mean Squared Error(mse) for predictors

#Parameters:
#Let denote predictors by p i.e.: 
#p=0(regression predictor); 
#p=1(adjusted survey predictor);
#p=gamma(Empirical BLUP) 
#p=3(survey predictor)
#crop= corn or soybeans

MSE.hat = function(p, crop){ 
  if (length(p) == 1){
    p = rep(p, 12)
  }
  if (crop == "corn"){
    sigmaestimate = Corn.sigmaestimate
    randsigma = Corn.randsigma 
    gamma = corn.gamma
    popmean = Corn.popmean
    Segments.mean = Segments.mean$PixelsCorn
    VHat = VHat.Corn} 
  else {
    sigmaestimate = Soybeans.sigmaestimate
    randsigma = Soybeans.randsigma
    gamma = soybeans.gamma
    popmean = Soybeans.popmean
    Segments.mean = Segments.mean$PixelsSoybeans
    VHat = VHat.Soybeans
  }
  
  res = rep(0, 12)
  for (i in 1:12){
    if (p[1] == 2){
      aux = (cbind(1, popmean[i]) - cbind(1, Segments.mean[i]))
      
      res[i] = sigmaestimate/Num.observed[i] + aux %*% VHat %*% t(aux)
    } else {
      aux1 = (cbind(1, popmean[i]) - p[i]*cbind(1, Segments.mean[i]))
      aux2 = cbind(1, Segments.mean[i])
      
      equation1 = (1 - p[i])^2*randsigma + p[i]^2*sigmaestimate/Num.observed[i]
      equation2 = 2*(p[i] - gamma[i])*aux1 %*% VHat %*% t(aux2)
      equation3 = aux1 %*% VHat %*% t(aux1)
      
      res[i] = equation1 + equation2 + equation3
    }
  }
  return(res)
}

Corn.RegPredMSE<- MSE.hat(0, crop = "corn")
Soybeans.RegPredMSE <- MSE.hat(0, crop = "soybeans")
Corn.ASPMSE <- MSE.hat(1, crop = "corn")
Soybeans.ASPMSE <- MSE.hat(1, crop = "soybeans")
corn.EBLUPMSE <- MSE.hat(corn.gamma, crop = "corn")
soybeans.EBLUPMSE <- MSE.hat(soybeans.gamma, crop = "soybeans")
Corn.SPMSE <- MSE.hat(2, crop = "corn")
Soybeans.SPMSE <- MSE.hat(2, crop = "soybeans")

# Creating a dataframe for the predictors MSE
Dataframe.MSE <- data.frame(County = County.names,
                     Corn.RegPred = Corn.RegPredMSE,
                     Soybeans.RegPred = Soybeans.RegPredMSE,
                     Corn.ASP = Corn.ASPMSE,
                     Soybeans.ASP = Soybeans.ASPMSE,
                     corn.EBLUP = corn.EBLUPMSE,
                     soybeans.EBLUP = soybeans.EBLUPMSE,
                     Corn.SP = Corn.SPMSE,
                     Soybeans.SP = Soybeans.SPMSE)

view(Dataframe.MSE)



#3d

#Estimating total county field size for both crops
corn.segmentsCounty <- unique(Satellite.data$MeanPixelsCorn)
soybeans.segmentsCounty<- unique(Satellite.data$MeanPixelsSoybeans)

#Estimate the total county field size
corntotal.EBLUP <- corn.EBLUP * corn.segmentsCounty
Soybeanstotal.EBLUP <- soybeans.EBLUP * soybeans.segmentsCounty
Corntotal.SP <- Corn.SP * corn.segmentsCounty
Soybeanstotal.SP <- Soybeans.SP * soybeans.segmentsCounty

#Creating a dataframe for total BLUP and total Survey Predictor
Dataframe.total <- data.frame(County = County.names, 
                       corn.EBLUP = corntotal.EBLUP,
                       soybeans.EBLUP = Soybeanstotal.EBLUP,
                       Corn.SP = Corntotal.SP,
                       Soybeans.SP = Soybeanstotal.SP)
view(Dataframe.total)

#Plotting the results on the map of Iowa
states <- map_data("state")
iowa <- subset(states, region == "iowa")
counties <- map_data("county")
iowa.county <- subset(counties, region == "iowa")




iowacounty.polygon <- dplyr::select(iowa.county, long, lat, subregion)
centroids <- aggregate(iowacounty.polygon[,1:2], by=list(iowacounty.polygon$subregion), FUN = mean)
centroids$County <- unique(iowacounty.polygon[,3])
centroids <- dplyr::filter(centroids, County %in% tolower(as.character(County.names)))


#Creating a dataframe to categorize the crops

#Corn 
corn.category <- data.frame(total.eblup = "blup:", blup = round(Dataframe.total$corn.EBLUP, 0),
                       total.survey = "survey:", survey = round(Dataframe.total$Corn.SP, 0))
corn.category <- data.frame(blup = paste("", corn.category$total.eblup, "", corn.category$blup),
                       survey = paste(corn.category$total.survey, corn.category$survey))
corn.category <- data.frame(Total.crop = paste("", corn.category$blup, "\n", corn.category$survey))


#Soybeans 
soybeans.category <- data.frame(total.eblup = "blup:", blup = round(Dataframe.total$soybeans.EBLUP, 0),
                      total.survey = "survey:", survey = round(Dataframe.total$Soybeans.SP, 0))
soybeans.category <-  data.frame(blup = paste("", soybeans.category$total.eblup, "", soybeans.category$blup),
                       survey = paste(soybeans.category$total_survey, soybeans.category$survey))
soybeans.category <- data.frame(Total.crop = paste("", soybeans.category$blup, "\n", soybeans.category$survey))

iowa.county$fill_value = 0
iowa.county$fill_value[iowa.county$subregion %in% tolower(as.character(County.names))] = 1


#plots

ggplot(data = iowa.county, aes(x = long, y = lat), fill = ) + 
  geom_polygon(aes(x = long, y = lat, group = group), 
               fill = "lightblue", color = "purple") +
  geom_polygon(data = iowa.county, aes(x = long, y = lat, group = group, fill = factor(fill_value)), 
               color = "white", show.legend = FALSE) +
  geom_polygon(data = iowa, aes(x = long, y = lat, group = group), 
               color = "black", fill = NA) +
  geom_label_repel(data = centroids, aes(long, lat, label = corn.category$Total.crop),
                   size = 3.5, alpha = 0.7, point.padding = 1.5,
                   min.segment.length = 0, segment.size = 0.6) +
  ggtitle("Estimated Total County Field Size for Corn in Iowa") +
  scale_fill_manual(values = c("blue","red")) 


ggplot(data = iowa.county, aes(x = long, y = lat), fill = ) + 
  geom_polygon(aes(x = long, y = lat, group = group), 
               fill = "lightblue", color = "purple") +
  geom_polygon(data = iowa.county, aes(x = long, y = lat, group = group, fill = factor(fill_value)), 
               color = "white", show.legend = FALSE) +
  geom_polygon(data = iowa, aes(x = long, y = lat, group = group), 
               color = "black", fill = NA) +
  geom_label_repel(data = centroids, aes(long, lat, label = soybeans.category$Total.crop),
                   size = 3.5, alpha = 0.7, point.padding = 1.5,
                   min.segment.length = 0, segment.size = 0.6) +
  ggtitle("Estimated Total County Field Size for Soybeans in Iowa") +
  scale_fill_manual(values = c("green","purple")) 

