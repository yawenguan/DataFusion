##  Neural Networks

* We fit the multilayer neural network using the deep learning package ``keras".
* The network is trained on PM2.5 concentration data, with Longtitude, Latitude, time, the selected covaraites and CMAQ data as inputs. 
* The final selected model has 2 hidden layers with 2000 and 100 neurons for the first and second layer, respectively. The final tuning parameter values are 200 epochs, 0.4 drop-out rate, 0.003 learning rate and 1024 minibatch size.



```r
load(file="DataExample.RData")
library(keras)

# transform the data back to the original scale
train_dat[,"PM"] <- exp(train_dat[,"PM"])
train_dat[,"CMAQ1"] <- exp(train_dat[,"CMAQ1"])

test_dat[,"PM"] <- exp(test_dat[,"PM"])
test_dat[,"CMAQ1"] <- exp(test_dat[,"CMAQ1"])

namesx = colnames(train_dat)[5:15]
varnames = c("Longitude","Latitude","t",namesx,"CMAQ1")

# Scale the variables to be in [0,1]
maxs    <- apply(train_dat[,colnames(train_dat)%in%varnames], 2, max)
mins    <- apply(train_dat[,colnames(train_dat)%in%varnames], 2, min)

train_dat[,colnames(train_dat)%in%varnames] <- scale(train_dat[,colnames(train_dat)%in%varnames], center = mins, scale = maxs - mins)
test_dat[,colnames(test_dat)%in%varnames]   <- scale(test_dat[,colnames(test_dat)%in%varnames], center = mins, scale = maxs - mins)

# Make a formula to go into the neural nets package
names <- c("PM", varnames)
test_names <- c("Site","t")
  
# fit NN with covariates = cmaq + covs
Y_test    <- test_dat[,"PM"] 
X_test    <- test_dat[,colnames(test_dat)%in%varnames]
Y_train   <- train_dat[,"PM"]
X_train   <- train_dat[,colnames(train_dat)%in%varnames]
  
model <- keras_model_sequential()
model %>%
  layer_dense(units = 2000,activation="relu",input_shape = c(ncol(X_train))) %>%
  layer_dropout(rate = 0.4) %>%
  layer_dense(units = 100,activation="relu") %>%
  layer_dense(units = 1)

model %>% compile(
  loss = 'mean_squared_error',
  optimizer = optimizer_adam(lr=0.003, decay = 0.0005),
  metrics = c('mse')
)


# Fit model to data
batch_size <- 1024
epochs <- 200
history <- model %>% fit(
  as.matrix(X_train), c(Y_train),
  batch_size = batch_size,
  epochs = epochs,
  verbose = 0,
  validation_split = 0.1
)

# plot(history)
pred <- data.frame(y = predict(model, as.matrix(X_test)))
pred <- pred[[1]] 
test <- Y_test 

plot(pred,test)
abline(a=0,b=1)
rmse.test = sqrt(mean((pred - test)^2))
mad.test = mean(abs(pred - test))
cor.test = cor(pred,test)

output = list(rmse.test,mad.test,cor.test)
print(output)
```