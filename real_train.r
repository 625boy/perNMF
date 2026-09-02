library(osqp)
library(Matrix)
library(dplyr)
library(progress)

source('./util.r')
source('./generate.r')

update_Wc <- function(var_list, kc, ki, lambda, gamma, weight){
  new_var_list <- var_list
  Wc <- var_list$Wc
  Wi <- var_list$Wi_list
  Hi <- var_list$Hi_list
  Gi <- var_list$Gi_list
  Yi <- var_list$Yi_list
  
  N <- length(Wi)
  p <- nrow(Wc)
  
  grad <- matrix(0, p, kc)
  ck_1 <- c(0)
  ck_2 <- c(0)

  update_res <- solve_osqp_column(V, ck, alpha)
  
  new_var_list$Wc <- update_res
  return(new_var_list)
}


update_Wi <- function(var_list, kc, ki, lambda, gamma, weight){
  new_var_list <- var_list
  Wc <- var_list$Wc
  Wi <- var_list$Wi_list
  Hi <- var_list$Hi_list
  Gi <- var_list$Gi_list
  Yi <- var_list$Yi_list
  
  N <- length(Wi)
  p <- nrow(Wc)
  
  new_Wi_list <- list()
  
  for(i in 1:N){
    update_res <- solve_osqp_column(V, ck, alpha)
    
    new_Wi_list[[i]] <- update_res
  }
  
  new_var_list$Wi_list <- new_Wi_list
  return(new_var_list)
}


update_Hi <- function(var_list, kc, ki, mu_h, weight){
  new_var_list <- var_list
  Wc <- var_list$Wc
  Wi <- var_list$Wi_list
  Hi <- var_list$Hi_list
  Gi <- var_list$Gi_list
  Yi <- var_list$Yi_list
  
  N <- length(Wi)
  p <- nrow(Wc)
  
  new_Hi_list <- list()
  
  for(i in 1:N){
    update_res <- solve_prox(V, ck, alpha)
    new_Hi_list[[i]] <- update_res
  }
  
  new_var_list$Hi_list <- new_Hi_list
  return(new_var_list)
}


update_Gi <- function(var_list, kc, ki, mu_g, weight){
  new_var_list <- var_list
  Wc <- var_list$Wc
  Wi <- var_list$Wi_list
  Hi <- var_list$Hi_list
  Gi <- var_list$Gi_list
  Yi <- var_list$Yi_list
  
  N <- length(Wi)
  p <- nrow(Wc)
  
  new_Gi_list <- list()
  
  for(i in 1:N){
    update_res <- solve_prox(V, ck, alpha)
    new_Gi_list[[i]] <- update_res
  }
  
  new_var_list$Gi_list <- new_Gi_list
  return(new_var_list)
}


run_userbatch <- function(train_data, kc, k, 
                          max_iter,
                          eps,
                          lambda, gamma, mu_h, mu_g, rho,
                          user_id,
                          save_dir = "./real") {
  gen_data <- generate_data_unif(train_data, kc, k)


  res_opt <- opt_osqp(
    max_iter,
    gen_data,
    eps,
    rho,
    lambda,
    gamma,
    mu_h,
    mu_g,
    user_id = user_id
  )

  file_name <- sprintf(
    "%s/res_osqp_user%d_%d_%d_%s_%s_%s_%s.Rda",
    save_dir, user_id, kc, k, eps, rho[1], gamma, lambda
  )

  save(res_opt, gen_data, file = file_name)
  
  message(sprintf("Finished (kc=%d, k=%d) -> %s", kc, k, file_name))
  
  invisible(list(res_opt = res_opt, gen_data = gen_data))
}

train_set_list <- list(
  train_set_user_1,
  train_set_user_2,
  train_set_user_3,
  train_set_user_4
)

sbatch_args <- commandArgs(trailingOnly = T)
user_batch <- as.integer(sbatch_args[1])
kc_batch <- as.integer(sbatch_args[2])
k_batch <- as.integer(sbatch_args[3])
eps_batch <- as.numeric(sbatch_args[4])

lambda_batch <- as.integer(sbatch_args[5])
gamma_batch <- as.integer(sbatch_args[6])
mu_h_batch <- as.integer(sbatch_args[7])
mu_g_batch <- as.integer(sbatch_args[8])

rho1_batch <- as.integer(sbatch_args[9])
rho2_batch <- as.integer(sbatch_args[10])
rho3_batch <- as.integer(sbatch_args[11])
rho4_batch <- as.integer(sbatch_args[12])
rho_batch <- c(rho1_batch, rho2_batch, rho3_batch, rho4_batch)

max_iter_batch <- as.integer(sbatch_args[13])

cat("Running with user =", user_batch, ", kc =", kc_batch, ", k =", k_batch ,"\n",
"rho = ", rho_batch, "lambda = ", lambda_batch, "gamma = ", gamma_batch, "\n'",
"mu_h = ", mu_h_batch, "mu_g = ", mu_g_batch, "max_iter = ", max_iter_batch, "\n")

train_batch <- train_set_list[[user_batch]]

run_userbatch(train_data = train_batch, kc = kc_batch, k = k_batch, 
            max_iter = max_iter_batch, eps = eps_batch, lambda = lambda_batch,
            gamma = gamma_batch, mu_h = mu_h_batch, mu_g = mu_g_batch,
            rho = rho_batch, user_id = user_batch, save_dir = "./real")