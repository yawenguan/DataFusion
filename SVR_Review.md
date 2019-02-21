##Support Vector Regression
* Here we use the r package *e1071* to fit support vector regression to a subset of the US 2011 air quality data; support vector regression is fitted to PM2.5 concentration at the log scale.
* The model is trained on PM2.5 concentration data, with Longtitude, Latitude, time, the selected covaraites and CMAQ data as inputs. 
* The final selected model was tuned on full data where gamma=0.3.
 
```r

load(file="DataExample.RData")
library("e1071")

## fit svm
m   <- svm(PM~.,data=train_dat, kernel = "radial", gamma=0.3)
## prediction
Yhat = predict(m ,newdata = test_dat, predict.all = T)

## plot on regular scale
pred <- exp(Yhat)
test <- exp(test_dat$PM)
plot(pred,test)
abline(a=0,b=1)

rmse.test = sqrt(mean((pred - test)^2))
mad.test = mean(abs(pred - test))
cor.test = cor(pred,test)

## collect results
output = list(rmse.test,mad.test,cor.test)
print(output)
```
