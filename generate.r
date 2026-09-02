library(osqp)
library(Matrix)
library(dplyr)
library(progress)

source('./util.r')

generate_data_unif = function(data, kc, k){
  set.seed(1)
  
  username <- sort(unique(data$user))
  nvec <- table(data$user)
  
  nonzero <- data[, -c(1:23, 4493)]
  nonzero <- nonzero[, colSums(nonzero) != 0]
  data_zero <- data.frame(user = data$user, nonzero)
  p <- ncol(nonzero)
  
  user_matrix_list <- lapply(split(data_zero, data_zero$user), function(x){
    d <- x[,-1]
    d <- as.matrix(d)
    return(d)
  })

  Yi_list <- list()
  for (i in username){
    Yi_list[[as.character(i)]] <- t(user_matrix_list[[as.character(i)]])
  }
  
  Yi_list <- Yi_list[as.character(username)]
  
  Wc <- matrix(runif(p*kc), nrow=p)
  Wi_list <- list()
  for (i in username){
    Wi_list[[as.character(i)]] <- matrix(runif(p*k), nrow=p)
  }
  Hi_list <- list()
  for (i in username){
    Hi_list[[as.character(i)]] <- matrix(runif(kc*nvec[[as.character(i)]]), nrow=kc)
  }
  Gi_list = list()
  for (i in username){
    Gi_list[[as.character(i)]] <- matrix(runif(k*nvec[[as.character(i)]]), nrow=k)
  }
  
  Wi_list <- Wi_list[as.character(username)]
  Hi_list <- Hi_list[as.character(username)]
  Gi_list <- Gi_list[as.character(username)]
  Yi_list <- Yi_list[as.character(username)]

  
  
  var_list = list(Wc = Wc,
                  Wi_list = Wi_list,
                  Hi_list = Hi_list,
                  Gi_list = Gi_list,
                  Yi_list = Yi_list
  )

  
  return(list(var_list = var_list, p = p, kc = kc, k = k, N = length(username), nvec = nvec))
}

generate_simul_grid <- function(
    N = 30,
    ni = 300,
    p = 900,
    Kc = 4,
    Ki = 4,
    Wc_row_sp = 0.25,
    overall_row_sp = 0,
    noise_lv = 0,
    seed){
  set.seed(seed)
  
  nvec <- rep(ni, N)
  
  nrow_nz_Wc <- p - round(p * Wc_row_sp)
  
  Wc_true <- matrix(0, nrow = p, ncol = Kc)
  if (Kc > 0 && nrow_nz_Wc > 0) {
    Wc_true[1:nrow_nz_Wc, ] <- runif(Kc * nrow_nz_Wc, min = 0, max = 1)
  }
  
  Wi_true_list <- list()
  Hi_true_list <- list()
  Gi_true_list <- list()
  
  raw_Yi_list <- list()
  sparse_Yi_list <- list()
  Yi_list <- list()
  
  for (i in 1:N) {
    nrow_nz_Wi <- round(p * Wc_row_sp)
    Wi_true_list[[i]] <- matrix(0, nrow = p, ncol = Ki)
    
    Hi_true <- if (Kc > 0) {
      matrix(runif(Kc * ni, min = 0, max = 1), nrow = Kc, ncol = ni)
    } else {
      matrix(0, nrow = 0, ncol = ni)
    }
    Hi_true_list[[i]] <- Hi_true
    
    Gi_true <- if (Ki > 0) {
      matrix(runif(Ki * ni, min = 0, max = 1), nrow = Ki, ncol = ni)
    } else {
      matrix(0, nrow = 0, ncol = ni)
    }
    Gi_true_list[[i]] <- Gi_true
    
    raw_Yi <- safe_reconstruct(
      Wc = Wc_true,
      Wi = Wi_true,
      Hi = Hi_true,
      Gi = Gi_true
    )
    raw_Yi_list[[i]] <- raw_Yi
    
    vec_len <- length(raw_Yi)
    zero_mask <- rep(1, vec_len)
    
    if (overall_row_sp > 0) {
      zero_ind <- sample(
        1:vec_len,
        size = round(vec_len * overall_row_sp),
        replace = FALSE
      )
      zero_mask[zero_ind] <- 0
    }
    
    sparse_Yi <- raw_Yi * matrix(
      zero_mask,
      nrow = nrow(raw_Yi),
      ncol = ncol(raw_Yi)
    )
    sparse_Yi_list[[i]] <- sparse_Yi
    
    if (noise_lv == 0) {
      Yi <- sparse_Yi
    } else {
      sig_sz <- sum(raw_Yi^2)
      
      noise_mat <- matrix(
        rnorm(vec_len),
        nrow = nrow(raw_Yi),
        ncol = ncol(raw_Yi)
      )
      
      noise_sz <- sum(noise_mat^2)
      scale <- sqrt((sig_sz / noise_lv) / noise_sz)
      noise_mat <- noise_mat * scale
      
      Yi <- pmax(sparse_Yi + noise_mat, 0)
    }
    
    Yi_list[[i]] <- Yi
  }
  
  Wc_init <- make_basis_matrix(p, Kc)
  
  Wi_init_list <- list()
  Hi_init_list <- list()
  Gi_init_list <- list()
  
  for (i in 1:N) {
    Wi_init_list[[i]] <- make_basis_matrix(p, Ki)
    Hi_init_list[[i]] <- make_coef_matrix(Kc, ni)
    Gi_init_list[[i]] <- make_coef_matrix(Ki, ni)
  }
  
  var_list <- list(
    Wc = Wc_init,
    Wi_list = Wi_init_list,
    Hi_list = Hi_init_list,
    Gi_list = Gi_init_list,
    Yi_list = Yi_list
  )
  
  true_list <- list(
    Wc_true = Wc_true,
    Wi_true_list = Wi_true_list,
    Hi_true_list = Hi_true_list,
    Gi_true_list = Gi_true_list,
    raw_Yi_list = raw_Yi_list,
    sparse_Yi_list = sparse_Yi_list
  )
  
  out <- list(
    var_list = var_list,
    true_list = true_list,
    p = p,
    kc = Kc,
    k = Ki,
    N = N,
    nvec = nvec,
    setting = list(
      Kc = Kc,
      Ki = Ki,
      Wc_row_sp = Wc_row_sp,
      overall_row_sp = overall_row_sp,
      noise_lv = noise_lv,
      seed = seed
    )
  )
  
  return(out)
}

generate_simul_test_grid <- function(train_gen_data,
                                     ni_test = 50,
                                     overall_row_sp = NULL,
                                     noise_lv = NULL,
                                     seed){
  set.seed(seed)
  
  p  <- train_gen_data$p
  N  <- train_gen_data$N
  Kc <- train_gen_data$kc
  Ki <- train_gen_data$k
  Kc_est <- train_gen_data$kc_est
  Ki_est <- train_gen_data$ki_est
  
  if (is.null(overall_row_sp)) {
    overall_row_sp <- train_gen_data$setting$overall_row_sp
  }
  
  if (is.null(noise_lv)) {
    noise_lv <- train_gen_data$setting$noise_lv
  }
  
  nvec_test <- rep(ni_test, N)
  
  Wc_true <- train_gen_data$true_list$Wc_true
  Wi_true_list <- train_gen_data$true_list$Wi_true_list
  
  Hi_test_true_list <- list()
  Gi_test_true_list <- list()
  
  raw_Yi_test_list <- list()
  sparse_Yi_test_list <- list()
  Yi_test_list <- list()
  
  for (i in 1:N) {
    Wi_true <- Wi_true_list[[i]]
    
    Hi_test_true <- if (Kc > 0) {
      matrix(
        runif(Kc * ni_test, min = 0, max = 1),
        nrow = Kc,
        ncol = ni_test
      )
    } else {
      matrix(0, nrow = 0, ncol = ni_test)
    }
    Hi_test_true_list[[i]] <- Hi_test_true
    
    Gi_test_true <- if (Ki > 0) {
      matrix(
        runif(Ki * ni_test, min = 0, max = 1),
        nrow = Ki,
        ncol = ni_test
      )
    } else {
      matrix(0, nrow = 0, ncol = ni_test)
    }
    Gi_test_true_list[[i]] <- Gi_test_true
    
    raw_Yi_test <- safe_reconstruct(
      Wc = Wc_true,
      Wi = Wi_true,
      Hi = Hi_test_true,
      Gi = Gi_test_true
    )
    raw_Yi_test_list[[i]] <- raw_Yi_test
    
    
    vec_len <- length(raw_Yi_test)
    zero_mask <- rep(1, vec_len)
    
    if (overall_row_sp > 0) {
      zero_ind <- sample(
        1:vec_len,
        size = round(vec_len * overall_row_sp),
        replace = FALSE
      )
      zero_mask[zero_ind] <- 0
    }
    
    sparse_Yi_test <- raw_Yi_test * matrix(
      zero_mask,
      nrow = nrow(raw_Yi_test),
      ncol = ncol(raw_Yi_test)
    )
    sparse_Yi_test_list[[i]] <- sparse_Yi_test
    
    
    if (noise_lv == 0) {
      Yi_test <- sparse_Yi_test
    } else {
      sig_sz <- sum(raw_Yi_test^2)
      
      noise_mat <- matrix(
        rnorm(vec_len),
        nrow = nrow(raw_Yi_test),
        ncol = ncol(raw_Yi_test)
      )
      
      noise_sz <- sum(noise_mat^2)
      scale <- sqrt((sig_sz / noise_lv) / noise_sz)
      noise_mat <- noise_mat * scale
      
      Yi_test <- pmax(sparse_Yi_test + noise_mat, 0)
    }
    
    Yi_test_list[[i]] <- Yi_test
  }
  
  test_data <- list(
    Yi_test_list = Yi_test_list,
    true_list = list(
      Wc_true = Wc_true,
      Wi_true_list = Wi_true_list,
      Hi_test_true_list = Hi_test_true_list,
      Gi_test_true_list = Gi_test_true_list,
      raw_Yi_test_list = raw_Yi_test_list,
      sparse_Yi_test_list = sparse_Yi_test_list
    ),
    p = p,
    kc = Kc,
    k = Ki,
    kc_est = Kc_est,
    ki_est = Ki_est,
    N = N,
    nvec = nvec_test,
    setting = list(
      Kc = Kc,
      Ki = Ki,
      Kc_est = Kc_est,
      Ki_est = Ki_est,
      Wc_row_sp = train_gen_data$setting$Wc_row_sp,
      overall_row_sp = overall_row_sp,
      noise_lv = noise_lv,
      train_seed = train_gen_data$setting$seed,
      test_seed = seed
    )
  )
  
  return(test_data)
}