library(osqp)
library(Matrix)
library(dplyr)
library(progress)
library(nnls)
library(caret)

source('./util.r')
source('./generate.r')
source('./real_optim.r')

opt_pernmf = function(res_opt, gen_data, train_data, test_data){
  Yi_list = list()
  Yi_test_list = list()
  username <- sort(unique(train_data$user))
  N <- length(username)
  
  user_matrix_list_test <- lapply(split(test_data, test_data$user), function(x){
    d <- x[, -c(1:23, 4493)]
    d <- as.matrix(d)
    return(d)
  })
  
  user_data_list <- lapply(split(train_data, train_data$user), function(x){
    d <- x[,c(1:22)]
    return(d)
  })
  
  user_data_list_test <- lapply(split(test_data, test_data$user), function(x){
    d <- x[,c(1:22)]
    return(d)
  })
  
  user_matrix_list_test <- user_matrix_list_test[as.character(username)]
  user_data_list <- user_data_list[as.character(username)]
  user_data_list_test <- user_data_list_test[as.character(username)]
  recon_train <- list()
  recon_test <- list()
  
  for (i in 1:N){
    dat = data.frame(user_data_list[[i]], H =  t(res_opt$updated_var_list$Hi_list[[i]]),
                     G = t(res_opt$updated_var_list$Gi_list[[i]]))
    Yi_list[[i]] = dat
    
    W_all <- cbind(res_opt$updated_var_list$Wc, res_opt$updated_var_list$Wi_list[[i]])
    Y_test <- t(user_matrix_list_test[[i]][,rownames(res_opt$updated_var_list$Yi_list[[i]])])
    fit_res <- matrix(0, nrow = gen_data$kc + gen_data$k, ncol = ncol(Y_test))
    
    for(j in 1:ncol(Y_test)){
      fit <- nnls(W_all, Y_test[,j])
      fit_res[,j] = fit$x
    }
    
    dat_test <- data.frame(user_data_list_test[[i]], t(fit_res))
    colnames(dat_test) <- colnames(dat)
    Yi_test_list[[i]] <- dat_test
  }
  
  train_opt = do.call(rbind, Yi_list)
  test_opt = do.call(rbind, Yi_test_list)
  return(list(train_opt = train_opt, test_opt = test_opt, recon_train = recon_train, recon_test = recon_test))
}

xgb_res_calc <- function(path){
  load(path)
  
  set.seed(1)
  
  res2 <- res_opt
  
  res_nnls <- opt_pernmf(res2, gen_data, train_set_user_1, test_set_user_1)
  
  
  xgb_fit_pernmf <- caret::train(like~., data = res_nnls$train_opt, method = 'xgbTree', 
                                 trControl = trainControl(allowParallel = T), nthread = 1, verbose = F, verbosity = 0)
  
  save(res2, res_nnls, xgb_fit_pernmf, gen_data, file = './real/real_result.Rda')
  
  return(list(train_rmse = pred_res(xgb_fit_pernmf, res_nnls$train_opt),
              train_recon = res_nnls$recon_train %>% unlist %>% mean,
              test_rmse = pred_res(xgb_fit_pernmf, res_nnls$test_opt),
              test_recon =  res_nnls$recon_test %>% unlist %>% mean))
}