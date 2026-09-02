library(osqp)
library(Matrix)
library(dplyr)
library(progress)

source('./util.r')
source('./generate.r')

update_Wc <- function(var_list, kc, ki, lambda, gamma, weight, upper){
  if (kc <= 0) return(var_list)
  
  new_var_list <- var_list
  Wc <- var_list$Wc
  Wi <- var_list$Wi_list
  Hi <- var_list$Hi_list
  Gi <- var_list$Gi_list
  Yi <- var_list$Yi_list
  
  N <- length(Wi)
  p <- nrow(Wc)
  cross_pen <- safe_cross_penalty(kc, ki, gamma)
  
  grad <- matrix(0, p, kc)
  ck_1 <- 0
  ck_2 <- 0
  
  new_var_list$Wc <- solve_osqp_column(V, ck, alpha)
  return(new_var_list)
}

update_Wc_nocol <- function(var_list, kc, ki, lambda, gamma, weight){
  if (kc <= 0) return(var_list)
  
  new_var_list <- var_list
  Wc <- var_list$Wc
  Wi <- var_list$Wi_list
  Hi <- var_list$Hi_list
  Gi <- var_list$Gi_list
  Yi <- var_list$Yi_list
  
  N <- length(Wi)
  p <- nrow(Wc)
  cross_pen <- safe_cross_penalty(kc, ki, gamma)
  
  grad <- matrix(0, p, kc)
  ck_1 <- 0
  ck_2 <- 0
    
  new_var_list$Wc <- solve_prox(V, ck, alpha)
  return(new_var_list)
}

update_Wi <- function(var_list, kc, ki, lambda, gamma, weight, upper){
  if (ki <= 0) return(var_list)
  
  new_var_list <- var_list
  Wc <- var_list$Wc
  Wi <- var_list$Wi_list
  Hi <- var_list$Hi_list
  Gi <- var_list$Gi_list
  Yi <- var_list$Yi_list
  
  p <- nrow(Wi[[1]])
  cross_pen <- safe_cross_penalty(kc, ki, gamma)
  
  new_Wi_list <- list()
  
  for (i in 1:length(Wi)){
    new_Wi_list[[i]] <- solve_osqp_column(V, ck, alpha)
  }
  
  new_var_list$Wi_list <- new_Wi_list
  return(new_var_list)
}

update_Wi_nocol <- function(var_list, kc, ki, lambda, gamma, weight){
  if (ki <= 0) return(var_list)
  
  new_var_list <- var_list
  Wc <- var_list$Wc
  Wi <- var_list$Wi_list
  Hi <- var_list$Hi_list
  Gi <- var_list$Gi_list
  Yi <- var_list$Yi_list
  
  p <- nrow(Wi[[1]])
  cross_pen <- safe_cross_penalty(kc, ki, gamma)
  
  new_Wi_list <- list()
  
  for (i in 1:length(Wi)){
    new_Wi_list[[i]] <- solve_prox(V, ck, alpha)
  }
  
  new_var_list$Wi_list <- new_Wi_list
  return(new_var_list)
}

update_Hi <- function(var_list, kc, ki, mu_h, weight){
  if (kc <= 0) return(var_list)
  
  new_var_list <- var_list
  Wc <- var_list$Wc
  Wi <- var_list$Wi_list
  Hi <- var_list$Hi_list
  Gi <- var_list$Gi_list
  Yi <- var_list$Yi_list
  
  new_Hi_list <- list()
  
  for (i in 1:length(Wi)){
    new_Hi_list[[i]] <- solve_prox(V, ck, alpha)
  }
  
  new_var_list$Hi_list <- new_Hi_list
  return(new_var_list)
}

update_Gi <- function(var_list, kc, ki, mu_g, weight){
  if (ki <= 0) return(var_list)
  
  new_var_list <- var_list
  Wc <- var_list$Wc
  Wi <- var_list$Wi_list
  Hi <- var_list$Hi_list
  Gi <- var_list$Gi_list
  Yi <- var_list$Yi_list
  
  new_Gi_list <- list()
  
  for (i in 1:length(Wi)){
    new_Gi_list[[i]] <- solve_prox(V, ck, alpha)
  }
  
  new_var_list$Gi_list <- new_Gi_list
  new_var_list
}

run_simulbatch <- function(Kc_ind, Ki_ind,
                           Wc_row_sp_ind, overall_ind, noise_ind,
                           max_iter,
                           eps,
                           lambda, gamma, mu_h, mu_g, rho,
                           nocol,
                           seed,
                           upper,
                           save_dir = "./simul",
                           param) {
  set.seed(seed)

  Wc_row_sp_vec <- c(0.25, 0.5, 0.75)
  overall_row_sp_vec <- c(0, 0.25, 0.5)
  noise_lv_vec <- c(0, 0.25, 0.5, 1)
  

  gen_data <- generate_simul_grid(
    N = 10,
    ni = 50,
    p = 900,
    Kc = Kc_ind,
    Ki = Ki_ind,
    Wc_row_sp = Wc_row_sp_vec[Wc_row_sp_ind],
    overall_row_sp = overall_row_sp_vec[overall_ind],
    noise_lv = noise_lv_vec[noise_ind],
    seed = seed
  )
  
  if (nocol == 0) {
    
    res_opt <- opt_osqp(
      max_iter,
      gen_data,
      eps,
      rho,
      lambda,
      gamma,
      mu_h,
      mu_g,
      upper = upper,
      param
    )
    
    file_name <- sprintf(
      "%s/res_simul_%d_%d_%d_%d_%d_upper%s_seed%s.Rda",
      save_dir,
      Kc_ind, Ki_ind, Wc_row_sp_ind, overall_ind, noise_ind, upper, seed
    )
    
  } else if (nocol == 1) {
    
    res_opt <- opt_osqp_nocol(
      max_iter,
      gen_data,
      eps,
      rho,
      lambda,
      gamma,
      mu_h,
      mu_g,
      param
    )
    
    file_name <- sprintf(
      "%s/res_simul_%d_%d_%d_%d_%d_nocol_seed%s.Rda",
      save_dir,
      Kc_ind, Ki_ind, Wc_row_sp_ind, overall_ind, noise_ind, seed
    )
    
  } else {
    stop("nocol must be either 0 or 1.")
  }

  recon_error <- compute_simul_recon_error(res_opt, gen_data)
  
  save(res_opt, gen_data, recon_error, file = file_name)
  
  message(sprintf(
    "Finished simulation: Kc_ind=%d, Ki_ind=%d, Wc_sp=%d, overall=%d, noise=%d, nocol=%d -> %s",
    Kc_ind, Ki_ind, Wc_row_sp_ind, overall_ind, noise_ind, nocol, file_name
  ))
  
  invisible(list(
    res_opt = res_opt,
    gen_data = gen_data,
    recon_error = recon_error
  ))
}

sbatch_args <- commandArgs(trailingOnly = T)
kc_batch <- as.integer(sbatch_args[1])
k_batch <- as.integer(sbatch_args[2])
wc_ind_batch <- as.integer(sbatch_args[3])
overall_ind_batch <- as.integer(sbatch_args[4])
noise_ind_batch <- as.integer(sbatch_args[5])
eps_batch <- as.numeric(sbatch_args[6])

lambda_batch <- as.integer(sbatch_args[7])
gamma_batch <- as.integer(sbatch_args[8])
mu_h_batch <- as.integer(sbatch_args[9])
mu_g_batch <- as.integer(sbatch_args[10])

rho1_batch <- as.integer(sbatch_args[11])
rho2_batch <- as.integer(sbatch_args[12])
rho3_batch <- as.integer(sbatch_args[13])
rho4_batch <- as.integer(sbatch_args[14])
rho_batch <- c(rho1_batch, rho2_batch, rho3_batch, rho4_batch)

max_iter_batch <- as.integer(sbatch_args[15])
nocol_batch <- as.integer(sbatch_args[16])
seed_batch <- as.integer(sbatch_args[17])
upper_batch <- as.numeric(sbatch_args[18])

cat("\n\n\n kc =", kc_batch, ", k =", k_batch , ", nocol = ", nocol_batch, "\n",
"rho = ", rho_batch, "lambda = ", lambda_batch, "gamma = ", gamma_batch, "upper = ", upper_batch, "\n'",
"mu_h = ", mu_h_batch, "mu_g = ", mu_g_batch, "max_iter = ", max_iter_batch, "seed = ", seed_batch, "\n",
"wc_ind = ", wc_ind_batch, "overall_ind = ", overall_ind_batch, "noise_ind = ", noise_ind_batch, "\n")

set.seed(seed_batch)

run_simulbatch(Kc_ind = kc_batch, Ki_ind = k_batch, Wc_row_sp_ind = wc_ind_batch, overall_ind = overall_ind_batch,
                noise_ind = noise_ind_batch, max_iter = max_iter_batch, eps = eps_batch,
                lambda = lambda_batch, gamma = gamma_batch, mu_h = mu_h_batch, mu_g = mu_g_batch, rho = rho_batch,
                nocol = nocol_batch, seed = seed_batch, upper = upper_batch, save_dir = "./simul", param = sbatch_args)