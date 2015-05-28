## collection of helper functions for performing crossvalidated model fits.
# TSAW.

#------------------------
# crossvalidation with each fold fit in parallel using parallel toolbox.

cross_validation_parallel <- function(x,y,fit_model,predict_model,n_folds=10,cores=4,seed,...){
  # set up for parallel sampling:
  library(doMC)
  options(cores=cores)
  registerDoMC()
  
  set.seed(seed)
  
  folds <- sample(rep(seq_len(n_folds), length.out=NROW(x)))
  
  predictions <- rep(NA, NROW(x))
  
  # do fits.  
  fold_predictions <- foreach(fold = 1:n_folds) %dopar% {
    training_ids <- which(folds!=fold)
    ifelse(NCOL(x)==1, training_x <- x[training_ids], training_x <- x[training_ids,])
    training_y <- y[training_ids]
    
    ifelse(NCOL(x)==1, test_x <- x[which(folds==fold)], test_x <- x[which(folds==fold),])
    
    model <- fit_model(training_x,training_y,...)
    preds <- predict_model(model,test_x,...)
    return(preds)
  }
  
  # fill in predictions.
  for (fold in 1:n_folds) {
    predictions[which(folds==fold)] <- fold_predictions[[fold]]
  }
  
  return(predictions)
}

#------------------------
# crossvalidation (serial, for testing).

cross_validation_serial <- function(x,y,fit_model,predict_model,n_folds=10,seed,...){
  
  set.seed(seed)
  
  folds <- sample(rep(seq_len(n_folds), length.out=NROW(x)))
  
  predictions <- rep(NA, NROW(x))
  
  # prediction generation function:
  fun <- function(fold){
    training_ids <- which(folds!=fold)
    ifelse(NCOL(x)==1, training_x <- x[training_ids], training_x <- x[training_ids,])
    training_y <- y[training_ids]
    
    ifelse(NCOL(x)==1, test_x <- x[which(folds==fold)], test_x <- x[which(folds==fold),])
    
    model <- fit_model(training_x,training_y,...)
    return(predict_model(model,test_x,...))
  }
  
  # do fits.
  fold_predictions <- sapply(1:n_folds,fun)
  
  # fill in predictions.
  for (fold in 1:n_folds) {
    predictions[which(folds==fold)] <- fold_predictions[[fold]]
  }
  
  return(predictions)
}

# #---------------------
# # fit model function for testing:
# test_fit_model <- function(x,y){
#   return(fm <- glm(y ~ x,family=binomial()))
# }
# 
# #---------------------
# # predict model function for testing:
# test_predict_model <- function(model,test_data){
#   df <- data.frame(x = test_data)
#   return(predict(model,df,type="response")) # predict method is default for lm. must take data frame with named columns.
# }
# 
# #-------------------
# # test:
# preds1 <- cross_validation_parallel(x,y,test_fit_model,test_predict_model,seed=12345)
# plot(y,preds1)
# 
# preds2 <- cross_validation_serial(x,y,test_fit_model,test_predict_model,seed=12345)
# plot(y,preds2)
# 
# all(preds1==preds2)

#------------------------
# evaluate classifier results using ROC.

output_roc <- function(y,y_hat){
  library(ROCR)
  # performance predictions:
  pred <- prediction(y_hat,y)
  a_roc <- performance(pred,measure='auc')
  roc_curve <- performance(pred,measure='tpr',x.measure='fpr')
  
  # output a list of useful metrics:
  roc <- list(auc = round(a_roc@y.values[[1]],digits=4),
              x = roc_curve@x.values[[1]],
              y = roc_curve@y.values[[1]])
  return(roc)
}
