## functions to predict response in test set from the mean of the training set.
# TSAW

#-------------------------
fit_mean_training <- function(x, y){
  
  fit <- mean(y)
  return(fit)
}

#-------------------------
# prediction function will predict new y values for a given fit and test_x.
predict_mean_test <- function(fit,test_x){
  
  predictions <- rep(fit,length=length(test_x))
  
  return(predictions)
}


