library(osqp)
library(Matrix)
library(dplyr)
library(progress)

solve_osqp_column <- function(V, c, alpha){
  p <- nrow(V)
  k <- ncol(V)
  
  P <- Diagonal(p)
  
  Aeq <- Matrix(1,1,p)
  Ann <- Diagonal(p)
  A <- rbind(Aeq, Ann)
  
  l <- c(1, rep(0,p))
  u <- c(1, rep(Inf, p))
  
  res <- matrix(0, p, k)
  
  for(j in 1:k){
    q <- -c * V[,j]
    
    fit <- osqp(P=P, q=q, A=A, l=l, u=u, pars=list(verbose = F, eps_abs = 1e-5, eps_rel = 1e-5,
                                                   adaptive_rho = F))
    
    fit_res <- fit$Solve()
    res[,j] <- pmax(0,fit_res$x)
  }
  
  return(res)
}

solve_prox <- function(V,c,alpha){
  shrink <- c/(c+2*alpha)
  
  out <- pmax(0,shrink*V)
  
  dim(out) <- dim(V)
  
  return(out)
}

safe_cross_penalty <- function(kc, k, gamma){
  if (kc <= 0 || k <= 0) return(0)
  return(gamma / (kc * k))
}

safe_spectral_norm <- function(M){
  if (any(dim(M) == 0)) return(0)
  return(norm(M, "2"))
}

make_basis_matrix <- function(p, k){
  if (k <= 0) {
    return(matrix(0, nrow = p, ncol = 0))
  }
  out <- matrix(runif(p * k), nrow = p)
  return(out)
}

make_coef_matrix <- function(k, n){
  if (k <= 0) {
    return(matrix(0, nrow = 0, ncol = n))
  }
  return(matrix(runif(k * n), nrow = k))
}

safe_reconstruct <- function(Wc, Wi, Hi, Gi) {
  p <- if (nrow(Wc) > 0) nrow(Wc) else nrow(Wi)
  n <- if (ncol(Hi) > 0) ncol(Hi) else ncol(Gi)
  
  out <- matrix(0, nrow = p, ncol = n)
  
  return(out)
}

compute_loss <- function(pre, updated, gen_data){
  norm_Wc <- c(0)
  norm_Wi <- c(0)
  norm_Hi <- c(0)
  norm_Gi <- c(0)
  norm_wcwi <- c(0)
  
  N = gen_data$N
  
  norm_Wc <- norm_Wc + norm(updated$Wc - pre$Wc, "2")
  
  for(i in 1:N){
    norm_Wi <- norm_Wi + norm(updated$Wi_list[[i]] - pre$Wi_list[[i]], "2")
    norm_Hi <- norm_Hi + norm(updated$Hi_list[[i]] - pre$Hi_list[[i]], "2")
    norm_Gi <- norm_Gi + norm(updated$Gi_list[[i]] - pre$Gi_list[[i]], "2")
    norm_wcwi <- norm_wcwi + norm(t(updated$Wc) %*% updated$Wi_list[[i]], "2")
  }
  
  return(list(norm_Wc = norm_Wc, norm_Wi = norm_Wi, 
              norm_Hi = norm_Hi, norm_Gi = norm_Gi, norm_wcwi = norm_wcwi))
}

opt_osqp <- function(n, gen_data, eps, weight, lambda, gamma, mu_h, mu_g, user_id){
  updated_var_list <- gen_data$var_list
  kc <- gen_data$kc; ki <- gen_data$k;p <- gen_data$p;N <- gen_data$N;nvec <- gen_data$nvec;
  
  res_Wc <- c(0)
  res_Wi <- c(0)
  res_Hi <- c(0)
  res_Gi <- c(0)
  res_wcwi <- c(0)
  
  for(iter in 1:n){
    pre_var_list <- updated_var_list
    
    updated_var_list <- update_Wc(updated_var_list, kc, ki, lambda, gamma, weight)
    updated_var_list <- update_Wi(updated_var_list, kc, ki, lambda, gamma, weight)
    updated_var_list <- update_Hi(updated_var_list, kc, ki, mu_h, weight)
    updated_var_list <- update_Gi(updated_var_list, kc, ki, mu_g, weight)
    
    comp_loss <- compute_loss(pre_var_list, updated_var_list, gen_data)
    
    res_Wc[iter] <- comp_loss$norm_Wc
    res_Wi[iter] <- comp_loss$norm_Wi
    res_Hi[iter] <- comp_loss$norm_Hi
    res_Gi[iter] <- comp_loss$norm_Gi
    res_wcwi[iter] <- comp_loss$norm_wcwi
    
    check_norm = sum(res_Wc[iter], res_Wi[iter], res_Hi[iter], res_Gi[iter])

    if(check_norm < eps){
      cat(paste0(iter, "\n"))
      break
    }
  }
  
  return(list(updated_var_list = updated_var_list,
              res_Wc = res_Wc,
              res_Wi = res_Wi,
              res_Hi = res_Hi,
              res_Gi = res_Gi,
              res_wcwi = res_wcwi))
}

opt_osqp_nocol <- function(n, gen_data, eps, weight, lambda, gamma, mu_h, mu_g, user_id){
  updated_var_list <- gen_data$var_list
  kc <- gen_data$kc; ki <- gen_data$k;p <- gen_data$p;N <- gen_data$N;nvec <- gen_data$nvec;
  
  res_Wc <- c(0)
  res_Wi <- c(0)
  res_Hi <- c(0)
  res_Gi <- c(0)
  res_wcwi <- c(0)
  
  for(iter in 1:n){
    pre_var_list <- updated_var_list
    
    updated_var_list <- update_Wc_nocol(updated_var_list, kc, ki, lambda, gamma, weight)
    updated_var_list <- update_Wi_nocol(updated_var_list, kc, ki, lambda, gamma, weight)
    updated_var_list <- update_Hi(updated_var_list, kc, ki, mu_h, weight)
    updated_var_list <- update_Gi(updated_var_list, kc, ki, mu_g, weight)
    
    comp_loss <- compute_loss(pre_var_list, updated_var_list, gen_data)
    
    res_Wc[iter] <- comp_loss$norm_Wc
    res_Wi[iter] <- comp_loss$norm_Wi
    res_Hi[iter] <- comp_loss$norm_Hi
    res_Gi[iter] <- comp_loss$norm_Gi
    res_wcwi[iter] <- comp_loss$norm_wcwi
    
    check_norm = sum(res_Wc[iter], res_Wi[iter], res_Hi[iter], res_Gi[iter])

    if(check_norm < eps){
      cat(paste0(iter, "\n"))
      break
    }
  }
  
  return(list(updated_var_list = updated_var_list,
              res_Wc = res_Wc,
              res_Wi = res_Wi,
              res_Hi = res_Hi,
              res_Gi = res_Gi,
              res_wcwi = res_wcwi))
}

compute_simul_recon_error <- function(res_opt, gen_data){
  
  Wc_hat <- res_opt$updated_var_list$Wc
  Wi_hat_list <- res_opt$updated_var_list$Wi_list
  Hi_hat_list <- res_opt$updated_var_list$Hi_list
  Gi_hat_list <- res_opt$updated_var_list$Gi_list
  Yi_list <- res_opt$updated_var_list$Yi_list
  
  N <- gen_data$N
  
  recon_obs <- numeric(N)
  recon_sparse <- numeric(N)
  recon_raw <- numeric(N)
  
  for (i in 1:N) {
    
    Yhat_i <- safe_reconstruct(
      Wc = Wc_hat,
      Wi = Wi_hat_list[[i]],
      Hi = Hi_hat_list[[i]],
      Gi = Gi_hat_list[[i]]
    )
    
    Yi_obs <- Yi_list[[i]]
    Yi_sparse <- gen_data$true_list$sparse_Yi_list[[i]]
    Yi_raw <- gen_data$true_list$raw_Yi_list[[i]]
    
    recon_obs[i] <- norm(Yi_obs - Yhat_i, "F")^2 / prod(dim(Yi_obs))
    recon_sparse[i] <- norm(Yi_sparse - Yhat_i, "F")^2 / prod(dim(Yi_sparse))
    recon_raw[i] <- norm(Yi_raw - Yhat_i, "F")^2 / prod(dim(Yi_raw))
  }
  
  list(
    recon_obs = recon_obs,
    recon_sparse = recon_sparse,
    recon_raw = recon_raw,
    mean_recon_obs = mean(recon_obs),
    mean_recon_sparse = mean(recon_sparse),
    mean_recon_raw = mean(recon_raw)
  )
}

pred_res = function(fit, test_data){
  pred = predict(fit, test_data)
  postResample(pred, test_data$like)
}