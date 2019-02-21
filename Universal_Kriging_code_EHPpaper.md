## Universal Kriging

* Here we use the R package *geoR* to perform Universal Kriging on a subset of the US 2011 air quality data; the downscaler model is fitted to PM2.5 concentration on the log scale.
*  The spatial process of PM2.5 concentration on a given day is modeled as having an exponential covariance function with covariance parameters sigma.sq, tau.sq and phi.
* The covariance parameters are assumed to be constant across time.
* The time-invariant covariance parameters are estimated by taking the mean across time of the covariance parameters estimated for each day.
* Daily estimates of the covariance parameters are estimated via REML using the geoR package.

* The R workspace contains:

* Data:

- train_dat: log PM2.5 concentration and covariates at 50 training sites from Jan 1st-7th, 2011
- test_dat:  log PM2.5 concentration and covariates at 50 testing sites from Jan 1st-7th, 2011


```r

load(file="DataExample.RData")
library(sp)
library(rgdal)
library(geoR)


names.var.train <- names(train_dat)
# no.days with data
n.days.train <- length(unique(train_dat$t))

# Defining the Lambert conformal projection to transform longitude and latitude
LCC <- "+proj=lcc +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0
+y_0=0 +ellps=GRS80 +datum=NAD83 +units=km +no_defs"


# Matrix with the estimated daily covariance parameters
cov.pars.daily <- matrix(NA,n.days.train,3)

# Vectors that will host all predicted values at test day and lower and upper bounds of 95%
# predictive interval
all.pred.uk <- NULL
all.low.bd.pred.uk <- NULL
all.upp.bd.pred.uk <- NULL


# This first loop is to estimate the covariance parameters of the exponential covariance
for(j in 1:n.days.train){
    # only using the data for the j-th day
    index.train.day <- which(train_dat$t==unique(train_dat$t)[j])
    coord.train.day <- cbind(train_dat$Longitude[index.train.day],train_dat$Latitude[index.train.day])
    # transforming latitude and longitude using the Lambert conformal projection
    sp.train.day <- SpatialPoints(coords=coord.train.day,proj4string=CRS("+proj=longlat +ellps=WGS84"))
    sp.train.trans <- spTransform(sp.train.day,CRS(LCC))
    sp.train.day.trans.coord <- as.matrix(cbind(sp.train.trans$coords.x1,sp.train.trans$coords.x2))
    
    pm.train.day <- train_dat$PM[index.train.day]
    cmaq.train.day <- train_dat$CMAQ1[index.train.day]
    
    # geodata object with observed daily PM concentration and daily CMAQ output
    pm.ukrig.df <- data.frame(cbind(sp.train.day.trans.coord[,1],sp.train.day.trans.coord[,2],pm.train.day,cmaq.train.day))
    pm.ukrig.geo <- as.geodata(pm.ukrig.df,coords.col=c(1,2),data.col=c(3,4))
    
    # estimation of the covariance parameters via REML
    pm.ukrig.reml <- likfit(geodata=pm.ukrig.geo,coord=pm.ukrig.geo$coords,data=pm.ukrig.geo$data[,1],trend = ~pm.ukrig.geo$data[, 2],cov.model="exponential",ini=c(0.5,2000.0), nugget=0.01, fix.nug = FALSE, lik.met="REML",lambda=1)
    
    cov.pars.daily[j,] <- c(pm.ukrig.reml$sigmasq,pm.ukrig.reml$phi,pm.ukrig.reml$tausq)
}

cov.pars.uk <- as.numeric(apply(cov.pars.daily,2,mean))


for(j in 1:n.days.train){
    # only using the data for the j-th day
    index.train.day <- which(train_dat$t==unique(train_dat$t)[j])
    index.test.day <- which(test_dat$t==unique(test_dat$t)[j])
    coord.train.day <- cbind(train_dat$Longitude[index.train.day],train_dat$Latitude[index.train.day])
    coord.test.day <- cbind(test_dat$Longitude[index.test.day],test_dat$Latitude[index.test.day])
    # transforming latitude and longitude using the Lambert conformal projection
    sp.train.day <- SpatialPoints(coords=coord.train.day,proj4string=CRS("+proj=longlat +ellps=WGS84"))
    sp.test.day <- SpatialPoints(coords=coord.test.day,proj4string=CRS("+proj=longlat +ellps=WGS84"))
    sp.train.trans <- spTransform(sp.train.day,CRS(LCC))
    sp.test.trans <- spTransform(sp.test.day,CRS(LCC))
    sp.train.day.trans.coord <- as.matrix(cbind(sp.train.trans$coords.x1,sp.train.trans$coords.x2))
    sp.test.day.trans.coord <- as.matrix(cbind(sp.test.trans$coords.x1,sp.test.trans$coords.x2))


    pm.train.day <- train_dat$PM[index.train.day]
    pm.test.day <- test_dat$PM[index.test.day]
    cmaq.train.day <- train_dat$CMAQ1[index.train.day]
    cmaq.test.day <- test_dat$CMAQ1[index.test.day]

    # geodata object with observed daily PM concentration and daily CMAQ output
    pm.ukrig.train.df <- data.frame(cbind(sp.train.day.trans.coord[,1],sp.train.day.trans.coord[,2],pm.train.day,cmaq.train.day))
    pm.ukrig.train.geo <- as.geodata(pm.ukrig.train.df,coords.col=c(1,2),data.col=c(3,4))
    pm.ukrig.test.df <- data.frame(cbind(sp.test.day.trans.coord[,1],sp.test.day.trans.coord[,2],cmaq.test.day))
    pm.ukrig.test.geo <- as.geodata(pm.ukrig.test.df,coords.col=c(1,2),data.col=c(3))
    
    # Specifying all the options for Universal Kriging
    kc.uk.control <- krige.control(type.krige="ok",trend.d=~pm.ukrig.train.geo$data[,2],trend.l=~pm.ukrig.test.geo$data,cov.model="exponential",cov.pars=c(cov.pars.uk[1],cov.pars.uk[2]),nugget=cov.pars.uk[3],lambda=1)

    # Predicting at test sites via Universal Kriging
    pred.uk.day <- krige.conv(pm.ukrig.train.geo, coords=pm.ukrig.train.geo$coords,data=pm.ukrig.train.geo$data[,1],locations=pm.ukrig.test.geo$coords,krige=kc.uk.control)

    # Predicted values and boundaries of the 95% prediction intervals
    #new.pred.uk.day <- as.numeric(pred.uk.day$predict)
    #new.low.bd.pred.uk.day <- new.pred.uk.day+qnorm(0.025)*sqrt(as.numeric(pred.uk.day$krige.var))
    #new.upp.bd.pred.uk.day <- new.pred.uk.day+qnorm(0.975)*sqrt(as.numeric(pred.uk.day$krige.var))
    
    new.pred.uk.day <- exp(as.numeric(pred.uk.day$predict))
    new.low.bd.pred.uk.day <- new.pred.uk.day+qnorm(0.025)*sqrt(as.numeric(pred.uk.day$krige.var))*new.pred.uk.day
    new.upp.bd.pred.uk.day <- new.pred.uk.day+qnorm(0.975)*sqrt(as.numeric(pred.uk.day$krige.var))*new.pred.uk.day
    
    # Concatenating all predictions
    all.pred.uk <- c(all.pred.uk,new.pred.uk.day)
    all.low.bd.pred.uk <- c(all.low.bd.pred.uk,new.low.bd.pred.uk.day)
    all.upp.bd.pred.uk <- c(all.upp.bd.pred.uk,new.upp.bd.pred.uk.day)
    
}

rmse.test = sqrt(mean((all.pred.uk - exp(test_dat$PM))^2))
mad.test = mean(abs(all.pred.uk  - exp(test_dat$PM)))
cor.test = cor(all.pred.uk,exp(test_dat$PM))
cov.test = mean(exp(test_dat$PM) < all.upp.bd.pred.uk & exp(test_dat$PM) > all.low.bd.pred.uk)

output = list(rmse.test,mad.test,cor.test)
print(output)
```





