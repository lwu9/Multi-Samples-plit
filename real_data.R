########## real data #########
# data <- read.table("/Users/lwu9/Documents/lab_proj/SP/real data/b09data.txt")
# source("/Users/lwu9/Documents/lab_proj/SP/functions2.r") # in laptop
# setwd("/Users/lwu9/Documents/lab_proj/SP/real data/")

data <- read.table("/home/lwu9/sampleSplit/b09data.txt")
source("/home/lwu9/sampleSplit/functions2.r") 

for (j in 1:dim(data)[2]) {
  data <- data[!data[,j]=='.',]
}
data <- data[complete.cases(data), ]
head(data)
a <- data$V3-1
survive <- data$V4
disease.free <- data$V6
 # y <- rep(0, dim(data)[1]); y[(survive+disease.free)==0] <- 1
y <- rep(0, dim(data)[1]); y[data$V7 >= 36] <- 1 # time to first event 
age <- data$V8
PR <- as.numeric(as.character(data$V10))
ER <- as.numeric(as.character(data$V9))
tumor.size <- as.numeric(as.character(data$V11))
num.node <- as.numeric(as.character(data$V12))
x <- as.matrix(data.frame(age, log1PR=log(1+PR)))#, log1ER=log(1+ER), tumor.size, log1node=log(1+num.node)))
set.seed(1)
dlist <- listdtr(y, a, x, stage.x = rep(1, dim(x)[2]), maxlen = 3L)
yrec <- predict(dlist, x, 1)
plot(dlist) ## real_data

ind <- (x[, 1] <= 50 & x[, 2] <= 2.398)
simon.trt <- rep(1, length(y)); simon.trt[ind] <- 0
dlist[[1]]$rule[1,5] <- 50; dlist[[1]]$rule[1,2]  <- "LL"; dlist[[1]]$rule[1,2]  <- "LL"; dlist[[1]]$rule[1,6]  <- 2.398
dlist[[1]]$rule <- dlist[[1]]$rule[-2,]
plot(dlist) ## gail_simon

list.rec <- as.numeric(as.character(yrec))

mean(simon.trt==list.rec)


n <- length(y); nB <- 600; rs <- 1:nB
registerDoParallel(4) 
pvalues <- foreach(seed=rs, .combine=rbind, .errorhandling = 'remove',
                         .packages=c('listdtr','mvtnorm','randomForest','ramify')) %dopar% {
                           set.seed(seed)
                           ind.in <- sample(1:n, size=n/2)
                           return(pValue.oneSplit.krr(ind.in, x, a, y, Switch=F, delta_m=T, true.diff=0))
                         }
stopImplicitCluster()
print(dim(pvalues)); B <- dim(pvalues)[1]
# pvalues <- pValue.multiSplit.krr(x, a, y, nB=nB, Switch=F, delta_m = T)
save(pvalues, file="real_data.RData")

gamma.min <- 0.01; nRep <- 1; alpha <- 0.05
pvalues <- unlist(pvalues[,1])
p.values <- c(pvalues[1])
for (b in 2:588) {
    pValues <- pvalues[1:b]
    Q_inf <- optimize(Q.gamma, c(gamma.min, 1), tol = 0.00001, pvalues = pValues)$objective
    p.value <- min(1, (1-log(gamma.min))*Q_inf) ## adaptive way to find gamma
    # p.value <- Q.gamma(0.5, pValues)
    p.values <-c(p.values, p.value)
}
cat(p.values,"\n")