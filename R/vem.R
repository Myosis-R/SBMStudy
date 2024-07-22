# programmation de l'algorithme variationnel Expectation-Maximization
# pour l'inférence des params d'un SBM avec K fixé

ve_f = function(N,K,pi_v,gamma_m,adj_v,Z_v){
  Z_v_new = Z_v
  for(i in 1:N){
    for(k in 1:K){
      Z_v_new[(i-1)*K+k] = pi_v[k] * prod(((rep(gamma_m[k,],N)**rep(adj_v[((i-1)*N+1):(i*N)],each=K)*(1-rep(gamma_m[k,],N))**rep(1-adj_v[((i-1)*N+1):(i*N)],each=K))**Z_v)[-c(((i-1)*K+1):(i*K))])
    }
    Z_v_new[((i-1)*K+1):(i*K)] = Z_v_new[((i-1)*K+1):(i*K)]/sum(Z_v_new[((i-1)*K+1):(i*K)])
  }
  return(Z_v_new)
}

#' Implementation of the Variationnal Expectation-Maximization algorithm for SBM
#'
#' @param pi_v Probability vector to belong to a given class.
#' @param gamma_m Probability matrix for an edge between two nodes knowing their class.
#' @param adj_m Adjacency matrix.
#' @returns list of SBM parameters infered and nodes' class assignment
#' @export
vem_f = function(pi_v,gamma_m,adj_m){
  R = 100
  adj_v = as.vector(adj_m)
  K = length(pi_v)
  N = dim(adj_m)[1]
  Z_m = matrix(rep(pi_v_init,N),nrow=K) # dim KxN
  stop = 0

  for(r in 1:R){
    Z_v = as.vector(Z_m)
    ve_f_ = purrr::partial(ve_f,N=N,K=K,pi_v=pi_v,gamma_m=gamma_m,adj_v=adj_v)
    FP = FixedPoint::FixedPoint(ve_f_,Z_v)
    if(FP$fpevals==1) return(list(tau=Z_m,gamma=gamma_m,pi=pi_v))
    Z_m = matrix(FP$FixedPoint,nrow=K)
    pi_v = apply(Z_m,1,mean)
    for(k in 1:K){
      for(l in 1:K){
        temp_m = matrix(rep(Z_m[k,],N),nrow = N)*t(matrix(rep(Z_m[l,],N),nrow = N))
        gamma_m[k,l] = sum(temp_m*adj_m)/sum(temp_m*(1-diag(N)))
      }
    }
  }
  return(list(tau=Z_m,gamma=gamma_m,pi=pi_v))
}


