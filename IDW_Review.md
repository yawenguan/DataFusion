##  Inverse Distance Weighting


```r
###############################################
# INPUT
#  s   := training data spatial locations (n_obs x 2)
#  Y   := training data response (vector of length n_obs)
#  sp  := test set locations (n_pred x 2)
#  pow := Power of the kernel smoothing function
#
# OUTPUT:
# Yhat := vector of predictions at the n_pred site
###############################################

load(file="DataExample.RData")

idw <- function(s,Y,sp,pow){
   library(fields)
   s    <- s[!is.na(Y),]
   Y    <- Y[!is.na(Y)]
   w    <- rdist(sp,s)^(-pow)
   Yhat <- w%*%Y/rowSums(w)
return(Yhat)}

# Create an empty variable for predictions
Yp <- c()

# Perform IDW each day with power=2
for(day in unique(train_dat$t)){
  sub_train_dat <- subset(train_dat, t==day)
  sub_test_dat <- subset(test_dat, t==day)
  Y = sub_train_dat$PM
  s = cbind(sub_train_dat$Longitude,sub_train_dat$Latitude)
  sp = cbind(sub_test_dat$Longitude,sub_test_dat$Latitude)
  Yp <- c(Yp,idw(s,Y,sp,pow=2))
}

# plot
pred <- exp(Yp)
test <- exp(test_dat$PM)
plot(pred,test)
abline(a=0,b=1)

rmse.test = sqrt(mean((pred - test)^2))
mad.test = mean(abs(pred - test))
cor.test = cor(pred,test)


output = list(rmse.test,mad.test,cor.test)
print(output)
```
