library(osqp)
library(Matrix)
library(dplyr)
library(progress)
library(nnls)

source('./util.r')
source('./generate.r')
source('./simul_optim.r')

compute_simul_test_recon_error <- function(res_opt, gen_data, test_data){
  
  Wc_hat <- res_opt$updated_var_list$Wc
  Wi_hat_list <- res_opt$updated_var_list$Wi_list
  
  N  <- gen_data$N
  Kc <- gen_data$kc_est
  Ki <- gen_data$ki_est
  total_rank <- Kc + Ki
  
  Yi_test_list <- test_data$Yi_test_list
  
  recon_test_obs <- numeric(N)
  recon_test_sparse <- numeric(N)
  recon_test_raw <- numeric(N)
  
  H_test_hat_list <- list()
  G_test_hat_list <- list()
  Yhat_test_list <- list()
  
  for (i in 1:N) {
    
    Y_test_obs <- Yi_test_list[[i]]
    ni_test <- ncol(Y_test_obs)
    
    fit_res <- matrix(0, nrow = total_rank, ncol = ni_test)
    
    if (total_rank > 0) {
      osqp_res <- solve_osqp(
        P = P,
        q = q,
        A = A,
        l = l,
        u = u,
        pars = list(
          verbose = FALSE,
          eps_abs = 1e-5,
          eps_rel = 1e-5,
          adaptive_rho = FALSE
        )
      )
      
      fit_res <- matrix(
        osqp_res$x,
        nrow = total_rank,
        ncol = ni_test
      )
    }

    idx_h <- if (Kc > 0) 1:Kc else integer(0)
    idx_g <- if (Ki > 0) (Kc + 1):(Kc + Ki) else integer(0)
    
    H_test_hat <- if (length(idx_h) > 0) {
      fit_res[idx_h, , drop = FALSE]
    } else {
      matrix(0, nrow = 0, ncol = ni_test)
    }
    
    G_test_hat <- if (length(idx_g) > 0) {
      fit_res[idx_g, , drop = FALSE]
    } else {
      matrix(0, nrow = 0, ncol = ni_test)
    }
    
    H_test_hat_list[[i]] <- H_test_hat
    G_test_hat_list[[i]] <- G_test_hat
    
    Yhat_test <- safe_reconstruct(
      Wc = Wc_hat,
      Wi = Wi_hat_list[[i]],
      Hi = H_test_hat,
      Gi = G_test_hat
    )
    
    Yhat_test_list[[i]] <- Yhat_test
    
    Y_test_sparse <- test_data$true_list$sparse_Yi_test_list[[i]]
    Y_test_raw <- test_data$true_list$raw_Yi_test_list[[i]]
    
    recon_test_obs[i] <- norm(Y_test_obs - Yhat_test, "F")^2 / prod(dim(Y_test_obs))
    recon_test_sparse[i] <- norm(Y_test_sparse - Yhat_test, "F")^2 / prod(dim(Y_test_sparse))
    recon_test_raw[i] <- norm(Y_test_raw - Yhat_test, "F")^2 / prod(dim(Y_test_raw))
  }
  
  out <- list(
    recon_test_obs = recon_test_obs,
    recon_test_sparse = recon_test_sparse,
    recon_test_raw = recon_test_raw,
    
    mean_recon_test_obs = mean(recon_test_obs),
    mean_recon_test_sparse = mean(recon_test_sparse),
    mean_recon_test_raw = mean(recon_test_raw),
    
    H_test_hat_list = H_test_hat_list,
    G_test_hat_list = G_test_hat_list,
    Yhat_test_list = Yhat_test_list
  )
  
  return(out)
}