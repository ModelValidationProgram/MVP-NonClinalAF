

setwd("~/Documents/GitHub/MVP-NonClinalAF/")
library(LEA)
data("tutorial")
write.lfmm(tutorial.R, "genotypes.lfmm")
write.geno(tutorial.R, "genotypes.geno")
write.env(tutorial.C, "gradients.env")

head(tutorial.R)
head(tutorial.C)


# Simulate non-null effect sizes for 10 target loci #individuals
n = 100
#loci
L = 1000
# Environmental variable
X = as.matrix(rnorm(n))
# effect sizes
B = rep(0, L)
target = sample(1:L, 10) 
B[target] = runif(10, -10, 10)

U = t(tcrossprod(as.matrix(c(-1,0.5,1.5)), X)) + matrix(rnorm(3*n), ncol = 3)

V <- matrix(rnorm(3*L), ncol = 3)

str(U)
str(V)

Y <- tcrossprod(as.matrix(X), B) + tcrossprod(U, V) +
  matrix(rnorm(n*L, sd = .5), nrow = n) 
Y <- matrix(as.numeric(Y > 0), ncol = L)
str(Y)
table(Y)

write.lfmm(Y, "genotypes.lfmm")
pc = pca("genotypes.lfmm", scale = TRUE)
tw = tracy.widom(pc)
plot(tw$percentage)
a <- tw$percentage[1:(length(tw$percentage)-2)]
b <- tw$percentage[2:(length(tw$percentage)-1)]
K <- max(which(a > b*1.5))

mod <- lfmm2(input = Y, env = X, K = 3)
pv <- lfmm2.test(object = mod, input = Y,
                 env = X,
                 linear = TRUE)
plot(-log10(pv$pvalues), col = "grey", cex = .6, pch = 19) 
points(target, -log10(pv$pvalues[target]), col = "red")
