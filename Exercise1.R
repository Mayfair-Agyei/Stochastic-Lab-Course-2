
# installing tidyverse
install.packages("tidyverse")
library(tidyverse)
library(haven)
library(dplyr)
library(ggplot2)

#1a
# reading the data children.dta using read_dta command
  dtafile <- file.path(getwd(), "R","Stochastic Lab Course","childrenfinal.dta")
childrenfinal <- read_dta(dtafile)

# data cleaning. Removing variables starting with s,v and m followed by numbers
childrenfinal_new <- childrenfinal %>% select(-matches("^[svm][0-9]"))

childrenfinal_new1<- childrenfinal_new %>% mutate_if(is.double, as.double)


#1b
ChildrenfinalTib <- select(childrenfinal_new1, hypage,ruralfacto,female,zstunt,zweight,zwast,adm2)

#scatter plot hypage against zstunt. Add a smooth line

ggplot(ChildrenfinalTib,aes(x=hypage, y=zstunt)) + geom_point()
ggplot(ChildrenfinalTib,aes(x=hypage, y=zstunt)) + geom_point() + geom_smooth(method = "gam", se=F, fullrange=F)
 # scatter plot for females and males
ggplot(ChildrenfinalTib, aes(x = hypage, y = zstunt,  colour = factor(female))) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "auto",se = F) +
  scale_colour_manual(labels = c("male", "female"), values = c("red", "green")) +
  guides(colour = guide_legend(title="Gender")) + theme_gray()
#scatter plot for rural and urban children
ggplot(ChildrenfinalTib, aes(x = hypage, y = zstunt,  colour = factor(ruralfacto))) +
  geom_point(alpha=0.4) +
  geom_smooth(method = "auto",se = F) +
  scale_colour_manual(labels = c("urban", "rural"), values = c("orange", "purple")) +
  guides(colour = guide_legend(title="Area")) + theme_gray()

#1c
# map of kenya
install.packages("maptools")
install.packages("rgdal")
install.packages("maps")
install.packages("mapdata")
library(maptools)
library(rgdal)
library(maps)
library(mapdata)
install.packages("raster")
library(raster)
KenyaMap <- raster::getData("GADM", country = "KE", level = 1)
KenyaMap1<-spTransform(KenyaMap, CRS("+init=epsg:32537"))
colnames(ChildrenfinalTib)[7] <- "NAME_1" 
ChildrenfinalTib<- ChildrenfinalTib[order(ChildrenfinalTib$NAME_1),]
KenyaMap1@data<- KenyaMap1@data[order(KenyaMap1@data$NAME_1),]
childrenFinal2 <- ChildrenfinalTib %>% group_by(NAME_1) %>%summarise(mean = mean(zstunt), n = n())
# adding missing county
childrenFinal2[nrow(childrenFinal2) + 1,] <- NA
childrenFinal2$NAME_1[47] <- "Isiolo"
childrenFinal2<- childrenFinal2[order(childrenFinal2$NAME_1),]

KenyaMap1@data$id <- rownames(KenyaMap1@data)
KenyaMap1@data <- mutate(KenyaMap1@data, zstunt.mean= childrenFinal2$mean)
KenyaDF <- fortify(KenyaMap1)
KenyaDF <- full_join(KenyaDF,KenyaMap1@data, by="id")

centroidsDF <- as.data.frame(coordinates(KenyaMap1))
names(centroidsDF) <- c("long", "lat")
childrenFinal2<- childrenFinal2[order(childrenFinal2$NAME_1),]
centroidsDF$NAME_1 <- KenyaMap1@data$NAME_1
centroidsDF$zstunt.mean <- childrenFinal2$mean

ggplot(data = KenyaDF, aes(x = long, y = lat, group = group, fill = zstunt.mean)) + 
  geom_polygon(color = "black", size = 0.25) +
  geom_text(data = centroidsDF, aes(x = long, y = lat, label = NAME_1, group = NULL), size = 3, position=position_jitter(width=1,height=1)) +
  scale_fill_distiller(name="Mean zstunt by county", palette = "Set2") +
  theme(aspect.ratio = 1) + theme_classic()

#1d
write.csv(ChildrenfinalTib,"R/Stochastic Lab Course/ChildrenfinalTib.txt")


