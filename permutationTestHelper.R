# BH iter + permutation test

# Example:
  # permutation.number(k = 10, q = 0.2, K = 100, prob.error = 0.1, delta = 0.5)
permutation.number = function(k, q = 0.2, K, prob.error = 0.2, delta = 0.5, reuse = T, M.min = 0.1){ # M.min = 0.1
  if(reuse){
    mk = 2 / q * (K / k) * (log(1 / prob.error) + log(K)) / delta^2 * M.min # new
  }else{
    mk = 2 / q * (K / k) * (log(1 / prob.error) + log(K)) / delta^2 * M.min # new
  } 
  return(mk)
}


# BH iter + permutation test
# Input:
  # permutation.pval.func: a function takes in the index of the hypotheses
# Example:
  # p.true = rep(0.001, 100)
  # test = BH.permutation(p.true = p.true)
BH.permutation = function(permutation.pval.func, K, q = 0.2, prob.error = 0.2, delta = 0.5, trace = F, adaptive = T, n.permute = NULL){
  if(adaptive){
    # initialization
    R0 = R = K
    R.set = seq(1, K)
    m0 = 0
    m = ceiling(permutation.number(k = R, q = q, K = K, prob.error = prob.error, delta = delta))
    pval = rep(1, K); names(pval) = seq(1, K)
    threshold = q
    count = rep(0, K)
    for(i in 1 : K){
      if(trace){print(paste("iteration:", " ", i, ";", " ", "number of rejections:", " ", length(R.set), sep = ""))}
      if(m0 < m){
        pval[R.set] = (1 + (permutation.pval.func(R.set, M = (m - m0)) * (m - m0 + 1) - 1) + (pval[R.set] * (m0 + 1) - 1)) / (m + 1)
      }
      count[R.set] = count[R.set] + (m - m0)
      R = sum(pval[R.set] <= threshold)
      R.set = as.numeric(names(pval[pval <= threshold]))
      if(R == R0){break}
      if(R == 0){break} # if no more rejections, break
      threshold = q * R / K
      m0 = m
      m = ceiling(permutation.number(k = R, q = q, K = K, prob.error = prob.error, delta = delta))
      R0 = R
    }
    
    return(list(rejection.set = R.set, number.permutation = count, pval = pval))
  }
  # standard
  if(is.null(n.permute)){n.permute = 1 / (q / K) * log(K)}
  pval = permutation.pval.func(seq(1, K), M = n.permute)
  R.set = which(pval <= BH.threshold(pval = pval, q = q))
  return(list(rejection.set = R.set, number.permutation = rep(n.permute, K), pval = pval))
}


BH.threshold = function(pval, q = 0.1){
  # preprocess
  m = length(pval)
  pval = sort(pval, decreasing = FALSE)
  
  # compute the threshold of the BH procedure
  FDR.hat = m * pval/seq(1, m)
  pval.index = which(FDR.hat <= q)
  if(length(pval.index) == 0){return(-1e6)}
  threshold = pval[max(pval.index)]
  
  return(threshold)  
}


BH.infty = function(p.true, q){
  threshold = BH.threshold(pval = p.true, q = q)
  return(list(rejection.set = which(p.true <= threshold)))
}

