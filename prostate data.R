# BH + permutation tests
# prostate data
source("/Users/zijungao/Desktop/Research/Qingyuan/iterative BH/code/qcHelper.R")
source("/Users/zijungao/Desktop/Research/Qingyuan/iterative BH/code/permutationTestHelper.R")

# load sda library
library("sda")

# helper functions
# permutation p-value based on the t-test statistics
# the function uses the dataset singh2002.
# index: the index of the hypotheses being tested. Default is NULL, and all 6033 hypotheses in the prostate data will be used.
# n: number of permutations
# Example: 
# set.seed(318); permutation.pval.func(n = 10)[1:10]
permutation.pval.func = function(index = NULL, M = 100){
  if(is.null(index)){index = seq(1, dim(singh2002$x)[2])}
  t.stats.permute = t(replicate(M, {
    sapply(index, function(i){W.permute = sample(dim(singh2002$x)[1]); t.test(x = singh2002$x[singh2002$y[W.permute] == "cancer", i], y = singh2002$x[singh2002$y[W.permute] == "healthy", i])$statistic})}, simplify = "array")) # independently permute for each gene
  # observed test statistics
  t.stats = sapply(index, function(i){t.test(x = singh2002$x[singh2002$y == "cancer", i], y = singh2002$x[singh2002$y == "healthy", i])$statistic}) 
  p.val = sapply(1 : length(index), function(i){(sum(t.stats.permute[, i] >= t.stats[i]) + 1) / (1 + M)})
  return(p.val)
}


# load prostate cancer dataset
data(singh2002) # 102 units (rows), 6033 genes (columns)
# hyperparamters
q = 0.1 # FDR level

# parametric p-values
# two sample t-test
t.stats = sapply(seq(1, dim(singh2002$x)[2]), 
                 function(i){t.test(x = singh2002$x[singh2002$y == "cancer", i], 
                                    y = singh2002$x[singh2002$y == "healthy", i])$statistic}) # compute two sample t-statistics
z.value = qnorm(pt(t.stats, df = 100)) # compute z-values
p.value = 1 - pnorm(z.value) # compute p-values
threshold = BH.threshold(pval = p.value, q = q) # BH, 0.000438
result.infty = list()
result.infty$pval = p.value
result.infty$rejection.set = which(p.value <= threshold) # q = 0.1: 28 rejections
hist(z.value, 
     main = "parametric p-values", xlab = "z-values",
     xlim = c(-4, 4), breaks = 50)

# permutation p-values
# standard BH
M.seq = seq(400, 2000, by = 400) 
result.standard = list(); result.standard$M = result.standard$rejection.set = result.standard$pval = result.standard$time = list()
set.seed(318)
for(i in seq(1, length(M.seq))){
  M = M.seq[i]
  result.temp = list()
  time.start = proc.time()
  result.temp$pval = permutation.pval.func(M = M)
  time.end = proc.time()
  
  print(M)
  result.temp$time = (time.end - time.start)[1] / 3600 # in hours
  result.standard$M[[i]] = M
  result.standard$rejection.set[[i]] = BH.infty(result.temp$pval, q = q)
  result.standard$pval[[i]] =  result.temp$pval
  result.standard$time[[i]] =  result.temp$time
}


# BH with acceleration
eps.seq =  seq(0.05, 0.25, by = 0.05) # seq(0.05, 0.25, by = 0.025) # 0.2
delta.seq = seq(0.1, 0.5, by = 0.1) # seq(0.1, 0.5, by = 0.05) # 1/3; 0.2
parameter.seq = data.frame(eps = c(eps.seq, rep(0.2, length(delta.seq))), delta = c(rep(1/3, length(eps.seq)), delta.seq))
result.accelerate = list(); result.accelerate$eps = result.accelerate$delta = result.accelerate$rejection.set = result.accelerate$number.permutation = result.accelerate$pval = result.accelerate$time = list()
set.seed(318)
for(i in 1 : dim(parameter.seq)[1]){
  eps = parameter.seq$eps[i]; delta = parameter.seq$delta[i]
  time.start = proc.time()
  BH.permutation.result = BH.permutation(permutation.pval.func = permutation.pval.func, K = dim(singh2002$x)[2], q = 0.1, prob.error = eps, delta = delta, trace = T)
  time.end = proc.time()

  print(i); print(eps); print(delta)
  result.accelerate$eps[[i]] = eps; result.accelerate$delta[[i]] = delta
  result.accelerate$time[[i]] = time
  result.accelerate$pval[[i]] = BH.permutation.result$pval
  result.accelerate$rejection.set[[i]] = BH.permutation.result$rejection.set
  result.accelerate$number.permutation[[i]] = BH.permutation.result$number.permutation
}

# results
print(result)