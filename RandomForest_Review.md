##  Random Forests

* Here we use the r package *randomForest* to fit random forest to a subset of the US 2011 air quality data; random forest is fitted to PM2.5 concentration at the log scale.
* We use the default value for $m = \left \lfloor{p/3}\right \rfloor$, the number of covariates selected at each split.
* The final model selected has n=500 trees, which minimizes the out-of-bag errors.
* The R workspace contains:

    + train_dat: log PM2.5 concentration and covariates at 50 training sites from Jan 1st-7th, 2011
    + test_dat:  log PM2.5 concentration and covariates at 50 testing sites from Jan 1st-7th, 2011


```r
load(file="DataExample.RData") # contains the data and names of the selected covariates
library(randomForest)
# selected covariates and no cmaq
namesx <- colnames(train_dat)[5:15]
cat(namesx)
varnames = c("Longitude","Latitude","t",namesx,"CMAQ1")
ntrees = 500 # fix the number of trees

# fit random forest
rf.PM=randomForest(formula(paste0("PM~",paste0(varnames,collapse = "+"))),data=train_dat,importance=TRUE,ntree=ntrees)
# prediction
yhat.rf = predict(rf.PM ,newdata = test_dat ,predict.all = T)
  
# exponentiate the data back to the original data scale
expyhat.rf = exp(yhat.rf$individual)
ytest   = exp(test_dat$PM)

# compute prediction quantiles
LU = apply(expyhat.rf,1,quantile,prob=c(0.025,0.975))

oob.errors = rf.PM$mse[ntrees]
rmse.test = sqrt(mean((rowMeans(expyhat.rf) - ytest)^2))
mad.test = mean(abs(rowMeans(expyhat.rf)  - ytest))
cor.test = cor(rowMeans(expyhat.rf),ytest)
cov.test = mean(ytest<LU[2,] & ytest >LU[1,])

# collect output
output = data.frame(oob = oob.errors,rmse = rmse.test, mad = mad.test,cor = cor.test, cov = cov.test)
print(output)
```