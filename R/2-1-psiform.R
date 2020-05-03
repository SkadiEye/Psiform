###########################################################
### ttsDE

#' A tail-based test statistics (TTS) for differential expression detection.
#'
#' @param ... One or multiple matrices of expression levels from multiple platforms.
#' @param y Subject labels.
#' @param G Pathway information.
#' @param dim_u Number of dimension of u.
#' @param weight Type of weight on the pathway. Possible values include "sqrt-degree",
#'   "none" or "degree". The default is "none".
#' @param scale Whether to scale input matrices.
#' @param lambda Hyper-parameters in SCAD, TLP, and L1 of phenotype differences.
#' @param a Hyper-parameter in SCAD. Default is 3.7.
#' @param param_init If not NULL, should be a list including c, beta, u, v and sigma_sq.
#' @param n_iter Number of maximum iterations.
#' @param epsilon A small number to avoid zero-denominators.
#' @param conv_thresh Threshold for convergence.
#' @param method Type of TLP.
#'
#' @return Returns a \code{Psiform} object.
#'
#' @importFrom CVXR Variable
#' @importFrom CVXR Minimize
#' @importFrom CVXR solve
#' @importFrom CVXR Problem
#' @importFrom pracma gramSchmidt
#' @importFrom CVXR max_elemwise
#' @importFrom stats na.omit
#' @importFrom methods new
#' @importFrom stats rnorm
#' @importFrom stats runif
#' @importFrom stats var
#'
#' @export
psiform <- function(..., y, G = NULL, dim_u = 5, weight = "none", scale = FALSE, lambda = NULL, a = 3.7,
                    param_init = NULL, n_iter = 100, epsilon = 10**-8, conv_thresh = 10**-6, method = 1) {

  mat_list <- list(...)
  d <- length(mat_list)
  n <- length(y)
  p <- numeric(d)
  for(i in 1:d) p[i] <- dim(mat_list[[i]])[2]
  p_total <- sum(p)
  k <- length(unique(y))
  if(is.null(lambda)) lambda <- list(SCAD = rep(0, d), TLP = rep(0, d), L1 = 0)
  if(d > 1 && length(lambda$SCAD) == 1) lambda$SCAD <- rep(lambda$SCAD, d)
  if(d > 1 && length(lambda$TLP) == 1) lambda$TLP <- rep(lambda$TLP, d)

  y_label <- levels(factor(y))
  ind <- as.numeric(factor(y))
  x <- matrix(NA, n, 0)
  for(i in 1:d) x <- cbind(x, mat_list[[i]])
  if(scale) {
    x <- scale(x)
    x_center <- attr(x, "scaled:center")
    x_scale <- attr(x, "scaled:scale")
  } else {
    x_center <- rep(0, dim(x)[2])
    x_scale <- rep(1, dim(x)[2])
  }
  MTM <- diag(table(ind))
  M <- t(sapply(ind, function(x) {y <- numeric(k); y[x] <- 1; y}))
  if(!is.null(colnames(x))) {
    x_names <- colnames(x)
  } else {
    x_names <- paste0("V", 1:p_total)
  }

  if(is.null(param_init)) {

    c_est <- stats::runif(k)
    c_est <- c_est - mean(c_est[ind])
    c_est <- c_est/sqrt(sum(c_est[ind]**2))
    beta_est <- stats::rnorm(p_total)
    u_est <- matrix(stats::rnorm(dim_u*n), n, dim_u)
    v_est <- matrix(stats::rnorm(p_total*dim_u), p_total, dim_u)
    sigma_sq_est <- apply(x, 2, stats::var)

    #### Initialization
    # n_iter0 <- 5
    # for(i in 1:n_iter0) {
    #   c_est <- c(solve(MTM*sum(beta_est**2/sigma_sq_est)) %*%
    #                (t(M) %*% ((x - u_est %*% t(v_est)) %*% (beta_est/sigma_sq_est))))
    #   u_est <- ((x - c_est[ind] %*% t(beta_est)) / (rep(1, n) %*% t(sigma_sq_est))) %*%
    #     v_est %*% solve(t(v_est) %*% (v_est / (sigma_sq_est %*% t(rep(1, dim_u)))))
    #
    #   ## Gram–Schmidt process
    #   c_est <- c_est - mean(c_est[ind])
    #   c_est <- c_est / sqrt(sum((c_est[ind])**2))
    #   u_est <- pracma::gramSchmidt(cbind(c_est[ind], u_est))$Q[, -1]
    #
    #   beta_est <- c(t(x - u_est %*% t(v_est)) %*% c_est[ind])
    #   v_est <- t(x - c_est[ind] %*% t(beta_est)) %*% u_est
    #   sigma_sq_est <- colMeans((x - c_est[ind] %*% t(beta_est) - u_est %*% t(v_est))**2)
    # }
  } else {

    c_est <- param_init$c
    beta_est <- param_init$beta
    u_est <- param_init$u
    v_est <- param_init$v
    sigma_sq_est <- param_init$sigma_sq
  }

  if(!is.null(G)) {
    degree_G <- sapply(1:p_total, function(x) sum(c(G) == x))
    if(weight == "sqrt-degree") W <- sqrt(degree_G)
    if(weight == "none") W <- rep(1, p_total)
    if(weight == "degree") W <- degree_G

    if(class(G) == "list") {
      Gx <- G[[1]]
      if(length(G) > 1) {
        add_ <- p[1]
        for(i in 2:length(G)) {
          Gx <- rbind(Gx, G[[i]] + add_)
          add_ <- add_ + p[i]
        }
      }
      G <- Gx
    }
  }

  n_iter_m <- 10
  loglik_ <- rep(NA, n_iter)
  loss_ <- rep(NA, n_iter)
  lambda_1 <- c(unlist(apply(cbind(lambda$SCAD, p), 1, function(x) rep(x[1], x[2]))))
  lambda_2 <- c(unlist(apply(cbind(lambda$TLP, p), 1, function(x) rep(x[1], x[2]))))
  lambda_c <- lambda$L1
  for(i in 1:n_iter) {

    old_c <- c_est
    old_beta <- beta_est

    if(method == 1) {

      #### METHOD 1 ####
      b_mat <- abs(c_est %*% t(rep(1, k)) - rep(1, k) %*% t(c_est)) + epsilon
      Q <- diag(colSums(1/b_mat)) - 1/b_mat
      c_est <- c(solve(MTM*sum(beta_est**2/sigma_sq_est) + 2*Q*lambda_c) %*%
                   (t(M) %*% ((x - u_est %*% t(v_est)) %*% (beta_est/sigma_sq_est))))

      u_est <- ((x - c_est[ind] %*% t(beta_est)) / (rep(1, n) %*% t(sigma_sq_est))) %*% v_est %*%
        solve(t(v_est) %*% (v_est / (sigma_sq_est %*% t(rep(1, dim_u)))))

      ## Gram–Schmidt process
      c_est <- c_est - mean(c_est[ind])
      c_est <- c_est / sqrt(sum((c_est[ind])**2))
      u_est <- pracma::gramSchmidt(cbind(c_est[ind], u_est))$Q[, -1]

      b_ <- CVXR::Variable(p_total)
      x_res <- x - u_est %*% t(v_est)
      b_w <- c(t(x_res) %*% c_est[ind])
      b_l1_weight <- lambda_1*((abs(beta_est) <= lambda_1) +
                                 (a*lambda_1 - abs(beta_est))*(abs(beta_est) > lambda_1)*(abs(beta_est) < a*lambda_1)/(a-1))

      if(!is.null(G)) {

        # browser()
        bg1_weight <- (1 + (abs(beta_est[G[, 1]])/W[G[, 1]] > lambda_1[G[, 1]]))/W[G[, 1]]*sign(beta_est[G[, 1]])
        bg2_weight <- (1 + (abs(beta_est[G[, 2]])/W[G[, 2]] > lambda_1[G[, 1]]))/W[G[, 2]]*sign(beta_est[G[, 2]])
        obj <- CVXR::Minimize((sum(b_^2*(1/sigma_sq_est)) - 2*sum(b_*(b_w/sigma_sq_est)))/2 + sum(b_l1_weight*abs(b_)) +
                                sum(2*lambda_2[G[, 1]]/lambda_1[G[, 1]]*
                                      CVXR::max_elemwise(abs(b_[G[, 1]])/W[G[, 1]] + CVXR::max_elemwise(abs(b_[G[, 2]])/W[G[, 2]] - lambda_1[G[, 1]], 0),
                                                   abs(b_[G[, 2]])/W[G[, 2]] + CVXR::max_elemwise(abs(b_[G[, 1]])/W[G[, 1]] - lambda_1[G[, 1]], 0))) -
                                (sum(lambda_2[G[, 1]]/lambda_1[G[, 1]]*bg1_weight*b_[G[, 1]]/W[G[, 1]]) +
                                   sum(lambda_2[G[, 1]]/lambda_1[G[, 1]]*bg2_weight*b_[G[, 2]]/W[G[, 2]])))
      } else {

        obj <- CVXR::Minimize((sum(b_^2*(1/sigma_sq_est)) - 2*sum(b_*(b_w/sigma_sq_est)))/2 + sum(b_l1_weight*abs(b_)))
      }

      sol <- CVXR::solve(CVXR::Problem(obj))
      if(length(sol$getValue(b_)) > 1)
        beta_est <- c(sol$getValue(b_))
      else
        break

      v_est <- t(x - c_est[ind] %*% t(beta_est)) %*% u_est
      sigma_sq_est <- colMeans((x - c_est[ind] %*% t(beta_est) - u_est %*% t(v_est))**2)

      loglik_[i] <- -sum(log(sigma_sq_est))*n/2
      if(!is.null(G)) {
        loss_[i] <- -loglik_[i] +
          sum(ifelse(abs(beta_est) > a*lambda_1, (a+1)*lambda_1**2/2,
                     ifelse(abs(beta_est) > lambda_1,
                            lambda_1*abs(beta_est) - (abs(beta_est) - lambda_1)**2/2/(a-1),
                            lambda_1*abs(beta_est)))) +
          sum(lambda_2[G[, 1]]*abs(sapply(abs(beta_est[G[, 1]])/W[G[, 1]]/lambda_1[G[, 1]], function(x) min(x, 1)) -
                                     sapply(abs(beta_est[G[, 2]])/W[G[, 2]]/lambda_1[G[, 2]], function(x) min(x, 1)))) +
          lambda_c*sum(abs((c_est %*% t(rep(1, k))) - (rep(1, k) - c_est)))
      } else {
        loss_[i] <- -loglik_[i] +
          sum(ifelse(abs(beta_est) > a*lambda_1, (a+1)*lambda_1**2/2,
                     ifelse(abs(beta_est) > lambda_1,
                            lambda_1*abs(beta_est) - (abs(beta_est) - lambda_1)**2/2/(a-1),
                            lambda_1*abs(beta_est)))) +
          lambda_c*sum(abs((c_est %*% t(rep(1, k))) - (rep(1, k) - c_est)))
      }
    } else {

      #### METHOD 2 ####
      for(j in 1:n_iter_m) {

        b_mat <- abs(c_est %*% t(rep(1, k)) - rep(1, k) %*% t(c_est)) + epsilon
        Q <- diag(colSums(1/b_mat)) - 1/b_mat
        c_est <- c(solve(MTM*sum(beta_est**2/sigma_sq_est) + 2*Q*lambda_c) %*%
                     (t(M) %*% ((x - u_est %*% t(v_est)) %*% (beta_est/sigma_sq_est))))
      }

      u_est <- ((x - c_est[ind] %*% t(beta_est)) / (rep(1, n) %*% t(sigma_sq_est))) %*% v_est %*%
        solve(t(v_est) %*% (v_est / (sigma_sq_est %*% t(rep(1, dim_u)))))

      ## Gram–Schmidt process
      c_est <- c_est - mean(c_est[ind])
      c_est <- c_est / sqrt(sum((c_est[ind])**2))
      u_est <- pracma::gramSchmidt(cbind(c_est[ind], u_est))$Q[, -1]

      # browser()
      b_l1_weight <- lambda_1*((abs(beta_est) <= lambda_1) +
                                 (a*lambda_1 - abs(beta_est))*(abs(beta_est) > lambda_1)*(abs(beta_est) < a*lambda_1)/(a-1))
      ind_G <- apply(cbind(G, lambda_1[G[, 1]]), 1, function(z) abs(beta_est[z[1]]/W[z[1]] - beta_est[z[2]/W[z[2]]]) <= z[3])
      G0 <- G[ind_G, ]

      for(j in 1:n_iter_m) {

        b0 <- abs(beta_est + epsilon)
        c0 <- abs(beta_est %*% t(rep(1, p_total)) - rep(1, p_total) %*% t(beta_est)) + epsilon
        if(!is.null(G)) {

          Q0 <- matrix(0, p_total, p_total)
          for(sss in 1:dim(G0)[1])
            Q0[c(G0[sss, 1], G0[sss, 2]), c(G0[sss, 1], G0[sss, 2])] <-
              Q0[c(G0[sss, 1], G0[sss, 2]), c(G0[sss, 1], G0[sss, 2])] +
              matrix(c(1/W[G0[sss, 1]]**2, -1/W[G0[sss, 1]]/W[G0[sss, 2]],
                       -1/W[G0[sss, 1]]/W[G0[sss, 2]], 1/W[G0[sss, 2]]**2), 2, 2)*(1/c0[G0[sss, 1], G0[sss, 2]])
          beta_est <- c(solve(diag(1/sigma_sq_est) + 2*diag(b_l1_weight/b0) + 2*Q0*lambda_2) %*%
                          ((t(x - u_est %*% t(v_est)) %*% c_est[ind])/sigma_sq_est))
        } else {

          beta_est <- c(solve(diag(1/sigma_sq_est) + 2*diag(b_l1_weight/b0)) %*%
                          ((t(x - u_est %*% t(v_est)) %*% c_est[ind])/sigma_sq_est))
        }
      }

      v_est <- t(x - c_est[ind] %*% t(beta_est)) %*% u_est
      sigma_sq_est <- colMeans((x - c_est[ind] %*% t(beta_est) - u_est %*% t(v_est))**2)

      loglik_[i] <- -sum(log(sigma_sq_est))*n/2
      if(!is.null(G)) {

        loss_[i] <- -loglik_[i] +
          sum(ifelse(abs(beta_est) > a*lambda_1, (a+1)*lambda_1**2/2,
                     ifelse(abs(beta_est) > lambda_1,
                            lambda_1*abs(beta_est) - (abs(beta_est) - lambda_1)**2/2/(a-1),
                            lambda_1*abs(beta_est)))) +
          sum(lambda_2[G[, 1]]*abs(sapply(abs(beta_est[G[, 1]]/W[G[, 1]]/lambda_1[G[, 1]] -
                                                beta_est[G[, 2]]/W[G[, 2]]/lambda_1[G[, 2]]), function(x) min(x, 1)))) +
          lambda_c*sum(abs((c_est %*% t(rep(1, k))) - (rep(1, k) - c_est)))
      } else {

        loss_[i] <- -loglik_[i] +
          sum(ifelse(abs(beta_est) > a*lambda_1, (a+1)*lambda_1**2/2,
                     ifelse(abs(beta_est) > lambda_1,
                            lambda_1*abs(beta_est) - (abs(beta_est) - lambda_1)**2/2/(a-1),
                            lambda_1*abs(beta_est)))) +
          lambda_c*sum(abs((c_est %*% t(rep(1, k))) - (rep(1, k) - c_est)))
      }
    }

    print(c(c_est, loglik_[i], loss_[i]))
    print(Sys.time())
    print(c(mean(abs(c_est - old_c))/mean(abs(old_c) + 1), mean(abs(beta_est - old_beta))/mean(abs(old_beta) + 1)))
    if(max(mean(abs(c_est - old_c))/mean(abs(old_c) + 1),
           mean(abs(beta_est - old_beta))/mean(abs(old_beta) + 1)) < conv_thresh)
      break
  }
  names(c_est) <- y_label
  names(beta_est) <- x_names

  return(methods::new("Psiform", c = c_est, beta = beta_est, u = u_est, v = v_est, sigma_sq = sigma_sq_est,
                      param = list(is.null.G = is.null(G), dim_u = dim_u, weight = weight,
                                   scale = scale, lambda = lambda, a = a, n_iter = n_iter,
                                   epsilon = epsilon, conv_thresh = conv_thresh, method = method,
                                   x_scale = x_scale, x_center = x_center, log_lik = stats::na.omit(loglik_),
                                   loss = na.omit(loss_))))
}
