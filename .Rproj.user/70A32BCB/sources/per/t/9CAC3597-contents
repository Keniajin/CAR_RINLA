## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, comment=NA, tidy = T)
## load the packages needed for this module
library(dplyr) ## package for data manipulation
library(ggplot2)
library(rgeos) ## Geomegry Engine- Open Source (GEOS)
library(maptools) ##  R package with useful map tools
library(broom) ## for converting a map to a data frame
library(rgdal) ## Geospatial Data Analysis Library (GDAL)"
library(rgeos) # # "Geomegry Engine- Open Source (GEOS)" -- check your shape file
library(spdep)
library(data.table)
library(haven) ## read in stata files
library("PerformanceAnalytics") ## EDA
library(INLA) ## inla

## remove the data in the environment
rm(list = ls())


## ------------------------------------------------------------------------
zim_shp <- maptools::readShapePoly("zim_shape/ZWE_adm2.shp",
                                        IDvar="ID_2")
plot(zim_shp, border="red", axes=TRUE, las=1 )


## ------------------------------------------------------------------------
zim_OGR <-  rgdal::readOGR("zim_shape/ZWE_adm2.shp")
plot(zim_OGR, border="blue", axes=TRUE, las=1)


## ------------------------------------------------------------------------
zim_valid <- data.table(valid=gIsValid(zim_OGR,byid=TRUE))
sum(zim_valid$valid==F)


## ------------------------------------------------------------------------
#zim_child_data <- haven::read_dta("data/Data2.dta")
zim_child_data <- readstata13::read.dta13("data/Data2.dta")
#dplyr::glimpse(zim_child_data , width = 5)


## ------------------------------------------------------------------------
## EDA
eda_data <- zim_child_data %>% 
  select(Stunting, Employed, b19, Education,  BMI)
PerformanceAnalytics::chart.Correlation(eda_data, histogram=TRUE, pch=19)


## ------------------------------------------------------------------------
zim_child_model <- zim_child_data %>%  
  select(id_2 , name_1 , name_2 , Stunting, Employed, b19, Education, b4, v025, BMI)
zim_child_model <- zim_child_model[complete.cases(zim_child_model),]

glimpse(zim_child_model)


## ------------------------------------------------------------------------
table(zim_shp@data$ID_2 %in% zim_child_model$id_2)
table(zim_shp@data$NAME_2 %in% zim_child_model$name_2)


## ------------------------------------------------------------------------
## Define the formula of the model
form_fit <- Stunting ~ Employed +b19  +  as.factor(Education) + b4+ v025+ BMI

## fit the model
lin_mod <- glm(form_fit, 
               data=zim_child_data, 
               family=gaussian(link="identity"))


## ------------------------------------------------------------------------
summary(lin_mod)$coefficient
AIC(lin_mod)  # AIC => 17272
BIC(lin_mod)  # BIC => 17333


## ------------------------------------------------------------------------
par(mfrow = c(2,2)) # display a unique layout for all graphs
plot(lin_mod)
par(mfrow = c(1,1)) 


## ---- tidy=FALSE---------------------------------------------------------
form_fit <-  Stunting ~1+ Employed +b19  + 
  as.factor(Education) + b4+ v025+ BMI

model_linear <- inla(form_fit,
                     family="gaussian",
                     data=zim_child_model,
                    control.compute=list(dic=TRUE, waic=TRUE))


## ------------------------------------------------------------------------
summary(model_linear)
summary(lin_mod)


## ------------------------------------------------------------------------
round(model_linear$summary.fixed[,1:5],3)


## ------------------------------------------------------------------------
par(mfrow = c(2,2)) # display a unique layout for all graphs
## Posterior density plot for the intercepts
plot(model_linear$marginals.fixed[[1]],
     type="l",
     main="",
     ylab="",
     xlab=expression(beta[0]))
plot(model_linear$marginals.fixed[[2]],
     type="l",
     main="",
     ylab="",
     xlab=expression(beta[1]))
par(mfrow = c(1,1))


## ------------------------------------------------------------------------
zim_nb<- spdep::poly2nb(zim_shp)

nb2INLA("zim_shape/zim_inla.graph", zim_nb)
zim_adj <- paste(getwd(),"/zim_shape/zim_inla.graph", sep="")


## ------------------------------------------------------------------------
H <- inla.read.graph(filename="zim_shape/zim_inla.graph")
image(inla.graph2matrix(H),xlab="",ylab="")


## ------------------------------------------------------------------------
formula_inla <- Stunting ~ 1 + as.factor(Employed) +
  b19  + as.factor(Education) + b4+ v025+ BMI+
  
  ##our CAR specification
  f(id_2, model="bym",graph=zim_adj, scale.model=TRUE,
    ## spefiying the priors for the unstri and str 
                     hyper=list(prec.unstruct=list(prior="loggamma",param=c(1,0.001)), 
                                prec.spatial=list(prior="loggamma",param=c(1,0.001))))



## ---- tidy=FALSE---------------------------------------------------------
model_inla <- inla(formula_inla,family="gaussian",
                     data=zim_child_model,
                     control.compute=list(dic=TRUE))
 



## ------------------------------------------------------------------------
summary(model_inla) #inla --> DIC 17218.91(pd=15.09) ##linear model--> DIC=17272.83(pd=9.444)


## ------------------------------------------------------------------------
round(model_inla$summary.fixed,3)


## ------------------------------------------------------------------------
round(head(model_inla$summary.random$id_2),3)


## ------------------------------------------------------------------------
csi <- model_inla$marginals.random$id_2[1:60]


## ------------------------------------------------------------------------
# *** Code for posterior probablity
a <- 0
prob_csi <- lapply(csi, function(x) {1 - inla.pmarginal(a, x)})
## for each location estimate the probability in continous
cont_prob_cs1 <- data.frame(maps_cont_prob_csi=unlist(prob_csi)) %>% 
  tibble::rownames_to_column("ID_2") %>% 
  mutate(ID_2=gsub("index.","",ID_2))

maps_cont_prob_csi <- cont_prob_cs1


## ---- eval=F-------------------------------------------------------------
# for each location estimate the probability in groups
prob_csi_cutoff <- c(0,0.2,0.4,0.6,0.8,1) ## can change accordingly
cat_prob_csi <- cut(unlist(prob_csi),
                    breaks=prob_csi_cutoff, 
                    include.lowest=TRUE)


maps_cat_prob_csi <- data.frame(ID_2=unique(zim_child_model$id_2), ## check whether it joins well
                                cat_prob_csi=cat_prob_csi)
maps_cat_prob_csi$ID_2 <- as.character(maps_cat_prob_csi$ID_2)


## ---- tidy=F ,, message=FALSE, warning=FALSE-----------------------------
zim_shp_df <- broom::tidy(zim_shp, region = "ID_2")
zim_shp_df <- zim_shp_df %>% 
 # left_join(maps_cat_prob_csi, by=c("id"="ID_2")) %>% 
  left_join(maps_cont_prob_csi, by=c("id"="ID_2"))

glimpse(zim_shp_df)


## ------------------------------------------------------------------------
p2 <- ggplot() + 
  geom_polygon(data = zim_shp_df, aes(x = long, y = lat,
                                      group = group,
                                      fill = maps_cont_prob_csi), 
               colour = "white") + theme_void() +
  ggtitle("RINLA Fit") + labs(fill = "P of  high HAZ") +
  scale_fill_continuous(high = "#fff7ec", low = "#7F0000")
p2


## ---- echo=FALSE ,  fig.height=10, eval=FALSE----------------------------
p1 <- p1 + scale_fill_continuous(high = "#fff7ec", low = "#7F0000") +
  ggtitle("Open Bugs Fit") + labs(fill = "Average HAZ")
gridExtra::grid.arrange(p1,p2, nrow=2 ) 

