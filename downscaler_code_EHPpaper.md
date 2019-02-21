## Downscaler model

* Here we use the R package *spBayes* to fit the downscaler model to a subset of the US 2011 air quality data; the downscaler model is fitted to PM2.5 concentration on the log scale.
* The predictions are backtransformed to the original scale for predictive performance assessment.
* Parameters of the downscaler model are estimated by running an MCMC algorithm for n.samples number of iterations. Posterior inference is based on the MCMC samples post burn-in (burn.in number of iterations).
* The downscaler model employs an exponential covariance function to model the spatial dependence in observed PM2.5 concentration measured at monitors with covariance parameters sigma.sq, phi and tau.sq.
* A flat priors is placed on the regression coefficients (e.g. intercept and coefficient of CMAQ), an Inverse Gamma prior is placed on sigma.sq and tau.sq, and a Uniform prior is specified for phi.
* The R workspace contains:

* Data:

- train_dat: log PM2.5 concentration and covariates at 50 training sites from Jan 1st-7th, 2011
- test_dat:  log PM2.5 concentration and covariates at 50 testing sites from Jan 1st-7th, 2011


```r

load(file="DataExample.RData")
library(sp)
library(rgdal)
library(spBayes)


names.var.train <- names(train_dat)
# no.days with data
n.days.train <- length(unique(train_dat$t))

# Defining the Lambert conformal projection to transform longitude and latitude
LCC <- "+proj=lcc +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0
+y_0=0 +ellps=GRS80 +datum=NAD83 +units=km +no_defs"

# Defining the number of MCMC samples
n.samples <- 100000
# Determining the number of iterations to discard for burn-in
burn.in <- (3*n.samples)/4
# Initial values for the MCMC algorithm:
starting.values <- list("phi"=0.05,"sigma.sq"=0.05,"tau.sq"=0.005)

# Tuning parameters for proposal distributions in Metropolis algorithm:
tuning.pars <- list("phi"=0.005,"sigma.sq"=0.1,"tau.sq"=0.05)

# Prior distributions: Uniform for range parameter (phi), Inverse Gamma
# for spatial variance (sigma.sq), Inverse Gamma for non-spatial variance (tau.sq)
priors.pars <- list("beta.Flat","phi.Unif"=c(0.0001,0.1),"sigma.sq.IG"=c(2,0.05),"tau.sq.IG"=c(2,0.005))

# Thinning factor for the MCMC samples
n.thin <- 2

# Vectors that will host all predicted values at test day and lower and upper bounds of 95%
# predictive interval
all.pred.down <- NULL
all.low.bd.pred.down <- NULL
all.upp.bd.pred.down <- NULL



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
    cmaq.test.day <-matrix(cbind(rep(1,length(pm.test.day)),test_dat$CMAQ1[index.test.day]),nrow=length(pm.test.day),ncol=2)
    

    # Fitting the downscaler model using the spBayes package
    coords.train.spbayes <- as.matrix(sp.train.day.trans.coord,nrow=dim(sp.train.day.trans.coord)[1],2)
    down.model.day <- spLM(pm.train.day~cmaq.train.day,coords=coords.train.spbayes,starting=starting.values,tuning=tuning.pars,priors=priors.pars,cov.model="exponential",n.samples=n.samples,verbose=FALSE, n.report=500)
    
    # Predicting at test sites; using a thinning of n.thin
    pred.down.day <- spPredict(down.model.day, pred.coords=sp.test.day.trans.coord, pred.covars=cmaq.test.day, start=burn.in, thin=n.thin)

    # Posterior predictive median and boundaries of the 95% predictive intervals on the original scale
    new.pred.down.day <- as.numeric(apply(t(exp(pred.down.day$p.y.predictive.samples)),2,quantile,0.5))
    new.low.bd.pred.down.day <- as.numeric(apply(t(exp(pred.down.day$p.y.predictive.samples)),2,quantile,0.025))
    new.upp.bd.pred.down.day <- as.numeric(apply(t(exp(pred.down.day$p.y.predictive.samples)),2,quantile,0.975))
    
    # Concatenating all predictions
    all.pred.down <- c(all.pred.down,new.pred.down.day)
    all.low.bd.pred.down <- c(all.low.bd.pred.down,new.low.bd.pred.down.day)
    all.upp.bd.pred.down <- c(all.upp.bd.pred.down,new.upp.bd.pred.down.day)
    
}

rmse.test = sqrt(mean((all.pred.down - exp(test_dat$PM))^2))
mad.test = mean(abs(all.pred.down  - exp(test_dat$PM)))
cor.test = cor(all.pred.down,exp(test_dat$PM))
cov.test = mean(exp(test_dat$PM) < all.upp.bd.pred.down & exp(test_dat$PM) > all.low.bd.pred.down)

output = list(rmse.test,mad.test,cor.test)
print(output)
```




