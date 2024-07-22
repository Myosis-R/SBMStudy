# programmation de l'algorithme de Lin Zhang 2019
# adaptant la validation croisée aux graphes

method_col = function(adj_m,train_idx,directed=TRUE){
  sbm_lib <- sbm::estimateSimpleSBM(adj_m[train_idx,train_idx],
                               "bernoulli",
                               directed = directed,
                               estimOptions = list(verbosity = 0, plot = FALSE))
  return(list(Z=sbm_lib$memberships,G=sbm_lib$connectParam$mean,P=sbm_lib$blockProp))
}

method_gc = function(adj_m,train_idx,directed=TRUE){
  N = length(train_idx)
  Z_init = rmultinom(N,size = 1,prob = c(0.2,0.2,0.2,0.2,0.2))
  Z = exICL_f(N,Z_init,adj_m[train_idx,train_idx])
  theta = Z_to_theta(Z,adj_m[train_idx,train_idx])
  return(list(Z=vec_to_num(Z),G=theta$gamma,P=theta$pi))
}

cross_val = function(adj,test_idx,method,directed=TRUE){
  N = dim(adj)[1]
  N_test = length(test_idx)
  train_idx = (1:N)[-test_idx]
  N_train = N - N_test

  sbm_lib <- method(adj,train_idx)

  Z_idx = sbm_lib$Z
  gamma = sbm_lib$G
  gamma[gamma == 0] = 1e-6 # avoids log(0)
  gamma[gamma == 1] = 1-1e-6
  pi = sbm_lib$P

  K = length(pi)
  idxs = cbind(rep(1:K,each=N_train),rep(Z_idx,K))
  g_temp = matrix(gamma[idxs],ncol=K)
  gt_temp = matrix(t(gamma)[idxs],ncol=K)

  Z_test = sapply(test_idx,
         \(i){
           which.max(
             t(adj[train_idx,i,drop=FALSE])%*%log(g_temp)+
               t(1-adj[train_idx,i,drop=FALSE])%*%log(1-g_temp)+
               adj[i,train_idx,drop=FALSE]%*%log(gt_temp)+
               (1-adj[i,train_idx,drop=FALSE])%*%log(1-gt_temp)+
               pi
           )
         })

  adj_test = matrix(gamma[cbind(rep(Z_test,each=N_test),rep(Z_test,N_test))],N_test,N_test,byrow = TRUE)
  adj_test[as.logical(diag(N_test))] = 0
  err = sum((adj[test_idx,test_idx]-adj_test)**2)/(N_test*(N_test-1))
  return(err)
}

#' Implementation of the vertices-fold algorithm for SBM
#'
#' @param adj Adjacency matrix.
#' @param v Number of folds
#' @param method method to use to infer sbm method_col/method_gc corresponding to VEM/exICL.
#' @returns Average vfold error.
#' @export
vfold = function(adj,v=5,method=method_col){
  N = dim(adj)[1]
  rand_nodes = sample(1:N)
  v_interval = as.integer(seq(1,N+1,length.out=(v+1)))
  v_err = sapply(1:v,
         \(i){
           cross_val(adj,rand_nodes[v_interval[i]:(v_interval[i+1]-1)],method)
         }
  )
  return(mean(v_err))
}


