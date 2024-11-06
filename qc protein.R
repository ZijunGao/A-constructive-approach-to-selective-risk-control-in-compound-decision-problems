# FDR curve with the protein dataset.
rm(list = ls())
library(ggplot2)
library(gridExtra)
source("qcHelper.R")

# Load data
data = readRDS("protein.rds")

# Construct test statistics
# Use negative control proteins to construct the reference distribution
sd.diff = mad((data$trt.129C - data$cnt.127C)[data$FP == 1])
mean.diff = mean((data$trt.129C - data$cnt.127C)[data$FP == 1])
test.stats = (data$trt.129C - data$cnt.127C - mean.diff )[data$TP == 1] / sd.diff # larger test stats suggests stronger evidence of rejecting the protein

# Compute p-values for different nulls
m = length(test.stats)
hypothesis.seq = seq(-3.5, 1.5, length = 400)

hypothesis.seq.prior.list = list(0,
                                 c(-0.5, 0),
                                 c(-1, -0.5, 0),
                                 c(-1.5,-1,-0.5,0)) 
q.curve.list = list(0.2,
                    c(0.3, 0.2), 
                    c(0.3, 0.2, 0.1),
                    c(0.5, 0.3, 0.2, 0.1)) 

g = list()
result = list()
for(k in 1 : length(q.curve.list)){

  hypothesis.seq.prior = hypothesis.seq.prior.list[[k]]
  q.curve = q.curve.list[[k]]
  
  p.value = sapply(hypothesis.seq.prior, function(x){pnorm(-test.stats - x, 0, 1)})

  # BH with simultaneously control of multiple FDRs
  p.val.aggregate = pmin(1, apply(p.value %*% (diag(length(q.curve)) / q.curve), 1, max))
  threshold = BH.threshold(pval = p.val.aggregate, q = 1)
  R = sum(p.val.aggregate <= threshold)
  
  # The above rejection set controls the FDRs corresponding to hypothesis.seq at q.curve.star for free.
  result[[k]] = q.star(hypothesis.seq = hypothesis.seq, hypothesis.seq.prior = hypothesis.seq.prior, q.curve = q.curve, m = m)
  # Create a data frame for ggplot
  plot.data = data.frame(hypothesis.seq = hypothesis.seq,
                         q.curve.naive = result[[k]]$q.curve.naive,
                         q.curve.star = result[[k]]$q.curve.star)
  
  # Plot
  g[[k]] = ggplot(plot.data, aes(x = hypothesis.seq)) +
    geom_ribbon(aes(ymin = q.curve.naive, ymax = q.curve.star), fill = "lightpink", alpha = 0.5) +
    geom_line(aes(y = q.curve.naive, color = "q(c)"), size = 1) +
    geom_line(aes(y = q.curve.star, color = "q*(c)"), size = 1) +
    scale_color_manual(name = "", values = c( "q(c)" = "black", "q*(c)" = "red")) +
    labs(x = "Effect size", y = "FDR level") +
    ggtitle(paste("R = ", R, sep= "")) + 
    xlim(-3.5, 1.5) + 
    theme_minimal() +
    theme(text = element_text(size = 16), 
         plot.title = element_text(hjust = 0.5),
         axis.title = element_text(size = 18), 
         axis.text = element_text(size = 14), 
         legend.title = element_text(size = 16),
         legend.text = element_text(size = 14), 
         legend.position = "bottom")
  
  print(g[[k]])
}

grid.arrange(g[[1]], g[[2]], g[[3]], g[[4]], ncol = 2, nrow = 2)


