# helper functions for q(c)
qc.mFDR = function(q, c, c0 = 0, prior.mean = 0, tau.lower = NULL, sigma = 1){
  c.lower = min(c)
  mu = c(prior.mean, prior.mean)
  Sigma = matrix(sigma^2, nrow = 2, ncol = 2); Sigma[2, 2] = Sigma[2, 2] + 1
  pic0 = pnorm(c0, mean = prior.mean, sd = sigma)
  pic = pnorm(c, mean = prior.mean, sd = sigma)
  if(is.null(tau.lower)){
    tau.seq = seq(-10 + prior.mean, 10 + prior.mean, length.out = 1000)
    FDR = sapply(tau.seq, function(x){pmvnorm(lower = c(-Inf, x), upper = c(c0, Inf), mean = mu, sigma = Sigma)}) / (1 - pnorm(tau.seq, mean = prior.mean, sd = sqrt(1 + sigma^2)))
    tau.lower = max(min(tau.seq[FDR <= pic0 * q]), qnorm(1 - q, 0, 1))
  }
  mFDRc = sapply(c, function(x){pmvnorm(lower = c(-Inf, tau.lower), upper = c(x, Inf), mean = mu, sigma = Sigma)})
  mFDRc0 = pmvnorm(lower = c(-Inf, tau.lower), upper = c(c0, Inf),
                   mean = mu, sigma = Sigma)
  return(q * pic0 / pic * mFDRc / mFDRc0[1])
}

qc.free = function(q, c, c0 = 0, K){
  c.lower = min(c)
  c.upper = 400
  ratio.lower = (pnorm(qnorm(q) + c0 - c))  / pnorm(qnorm(q) + c0 - c0)
  ratio.upper = pnorm(qnorm(q / K) + c0 - c)  / pnorm(qnorm(q / K) + c0 - c0)
  return(pmax(ratio.upper, ratio.lower) * q)
}

q.star = function(hypothesis.seq, hypothesis.seq.prior = 0, q.curve = 1, m){
  q.curve.naive = c(1, q.curve)[as.numeric(cut(hypothesis.seq, breaks = c(-Inf, hypothesis.seq.prior, Inf)))]
  q.curve.star = pmin(1, apply(sapply(seq(1, length(q.curve)), function(index){qc.free(q = q.curve[index], c = hypothesis.seq, c0 = hypothesis.seq.prior[index], K = m)}), 1, min))
  return(list(q.curve.naive = q.curve.naive,
              q.curve.star = q.curve.star))
}

BH.threshold = function(pval, q = 0.1){
  m = length(pval)
  pval = sort(pval, decreasing = FALSE)
  
  FDR.hat = m * pval/seq(1, m)
  pval.index = which(FDR.hat < q)
  if(length(pval.index) == 0){return(-1e6)}
  threshold = pval[max(pval.index)]
  return(threshold)
}
