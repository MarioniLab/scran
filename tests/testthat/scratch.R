set.seed(20000)
p1 <- runif(1000)
p2 <- runif(1000)
p3 <- runif(1000)

some.na <- sample(length(p1), length(p1)/2)
p1.na <- p1
p1.na[some.na] <- NA
p2.na <- p2
p2.na[-some.na] <- NA

library(scran)
method <- "simes"; weights<- NULL
ref <- combinePValues(p1[2], p3[2], method=method, weights=weights[c(1,3)])
ref2 <- combinePValues(p2[2], p3[2], method=method, weights=weights[c(2,3)])
out <- combinePValues(p1.na[2], p2.na[2], p3[2], method=method, weights=weights)


