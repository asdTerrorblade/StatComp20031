## -----------------------------------------------------------------------------
set.seed(1)

n <- 500
u <- runif(n)
x <- 2/sqrt(1-u) # Pareto(2,2) random numbers

# Graph the density histogram of the sample
hist(x, prob = TRUE, breaks = 50, main = expression(f(x) == 8/x^3)) 

# superimpose the Pareto(2,2) density function curve
y <- seq(2,20, 0.1)
lines(y, 8/y^3, col = 'red') 

## -----------------------------------------------------------------------------

n <- 10000
# generate three iid uniform distribution random number
u1 <- runif(n, min = -1, max = 1)
u2 <- runif(n, min = -1, max = 1)
u3 <- runif(n, min = -1, max = 1)

x <- seq(1:n) # initialize
# use the algorithm to generate variates from fe
for(i in 1:n)
{
  if((abs(u3[i]) >= abs(u2[i]))&(abs(u3[i]) >= abs(u1[i])))
    {x[i] <- u2[i]}
  else
    {x[i] <- u3[i]}
}

# construct the histogram
hist(x, probability = TRUE, main = expression(f(x)==3/4 * (1-x^2)))

# superimpose the density function curve
y <- seq(-1, 1, 0.1)
lines(y, 3/4 * (1 - y^2), col = 'red')

## -----------------------------------------------------------------------------
set.seed(1)

# generate random observations
n <- 1000
g <- rgamma(n, shape = 4, rate = 2)
y <- rexp(n,rate = g)

# graph the density histogram of the sample
hist(y, probability = TRUE, main = expression(f(y)==64/(1+y)^5), breaks = 50)

#superimpose the Pareto density curve
x <- seq(0, 5, 0.1)
lines(x, 64/(2 + x)^5, col = 'red')

## -----------------------------------------------------------------------------
set.seed(123)

n <- 1e4
u <- runif(n, min = 0, max = pi/3)
theta.hat <- mean(sin(u)) * pi / 3 # estimated value
theta = cos(0) - cos(pi / 3) # exact value
print(c(theta, theta.hat))


## -----------------------------------------------------------------------------
set.seed(1234)
n <- 1e4

# MC method function
MC <- function(number, antithetic) {
  u <- runif(number/2)
  if (!antithetic) v <- runif(number/2) else v <- 1 - u
  u <- c(u, v)
  population <- seq(1, number, 1)
  for (i in 1:number) {
    population[i] <- exp(u[i])
  }
  mean(population)
}

# function to compute empirical variance
variance <- function(Antithetic = FALSE){
  N <- 1e4
  population <- seq(1, N, 1)
  for (i in 1:length(population)){
    population[i] <- MC(n, antithetic = Antithetic)
  }
  var(population)
}

# estimate theta with the simple Monte Carlo method
theta.hat1 <- MC(n, antithetic = FALSE)

# estimate theta with the antithetic variate method
theta.hat2 <- MC(n, antithetic = TRUE)

# estimate variance of theta with the simple Monte Carlo method
theta.hat1.variance <- variance()

# estimate variance of theta with the antithetic variate method
theta.hat2.variance <- variance(Antithetic = TRUE)

e <- exp(1)
reduction <- (10*e-3*e^2-5)/(4*e-e^2-3) # theoretical the percent reduction value
reduction.hat <- theta.hat2.variance/theta.hat1.variance # empirical the percent reduction value
print(c(reduction,reduction.hat))

## -----------------------------------------------------------------------------
x <- seq(1,10,0.1)
w <- 2
f <- x^2*exp(-x^2/2)/(sqrt(2*pi))
f1 <- exp(-(x-1)^2/2)/(sqrt(2*pi))
f2 <- exp(-x/2)/(sqrt(2*pi))

plot(x, f/f1, type = 'l', main = '', ylab = '', lwd = w,ylim = c(0,2), xlim = c(1,2))
lines(x, f/f2, lty = 2, lwd = w)

## -----------------------------------------------------------------------------
n <- 1000
k <- 5
M <- 50
estimator <-numeric(k)
estimators <- numeric(M)
f <- function(x) {
  exp(-x) / (1+x^2) * (x > 0) * (x < 1)
}
set.seed(123)

for(i in 1:M){
  for (j in 1:k) {
    g <- function(x) {
      exp(-x)/(exp(-(j-1)/5)-exp(-j/5))
    }
    u <- runif(n/k)
    x <- -log(exp(-(j-1)/5) - (u * exp(-(j-1)/5) - u * exp(-j/5)))
    estimator[j] <- mean(f(x)/g(x))
   } 
  estimators[i] <- sum(estimator)
}
print(c(mean(estimators), sd(estimators)))


## -----------------------------------------------------------------------------
n <- 20
m <- 1000
alpha <- 0.025
UCL <- numeric(m)
set.seed(123)

for (i in 1:m){
  x <- rlnorm(n, meanlog = 0, sdlog = 2)
  lx <- log(x)
  theta.hat <- mean(lx)
  sigma.hat <- sd(lx)
  left <- theta.hat + qt(alpha, df = n-1) * sigma.hat / sqrt(n)
  right <- theta.hat - qt(alpha, df = n-1) * sigma.hat / sqrt(n)
  UCL[i] <- left < 0 & 0 < right
}
print(mean(UCL))

## -----------------------------------------------------------------------------
n <- 20
m <- 1000
UCL <- 0
set.seed(123)

for (i in 1:m) {
  x <- rchisq(n,2)
  left <- mean(x) - sd(x) * qt(0.975, n - 1)/(sqrt(n))
  right <- mean(x) + sd(x) * qt(0.975, n - 1)/(sqrt(n))
  if (mean(left)<2 & 2<mean(right)){
    UCL = UCL + 1
  }
}
print(UCL/m)

## -----------------------------------------------------------------------------
set.seed(123)
sk <- function(x) {
#the sample skewness coeff.
x.mean <- mean(x)
moment.3 <- mean((x - x.mean)^3)
moment.2 <- mean((x - x.mean)^2)
return( moment.3 / moment.2^1.5 )
}


alpha <- 0.1
n <- 30
m <- 2500
parameters <- seq(1,10)
N <- length(parameters)
beta.pwr <- numeric(N)
t.pwr <- numeric(N)
cv <- qnorm(1-alpha/2, 0, sqrt(6*(n-2)/((n+1)*(n+3))))
set.seed(123)


for (j in 1:N) {
  parameter <- parameters[j]
  beta.sktests <- numeric(m)
  t.sktests <- numeric(m)
  for (i in 1:m) {
    x <- rbeta(n, parameter, parameter)
    y <- rt(n, parameter)
    beta.sktests[i] <- as.integer(abs(sk(x))>=cv)
    t.sktests[i] <- as.integer(abs(sk(y))>=cv)
  }
  beta.pwr[j] <- mean(beta.sktests)
  t.pwr[j] <- mean(t.sktests)
}
plot(parameters, beta.pwr, type = 'b', ylim = c(0,1))
lines(parameters, t.pwr)

## -----------------------------------------------------------------------------
set.seed(123)
#count 5 test function
c5t <- function(x, y) {
  x.st <- x - mean(x)
  y.st <- y - mean(y)
  outx <- sum(x.st >  max(y.st)) + sum(x.st < min(y.st))
  outy <- sum(y.st > max(x.st)) + sum(y.st < min(x.st))
  return(as.integer(max(c(outx, outy)) > 5))
}

#F test function
Ft <- function(x, y) {
  Fvalue <- (sd(x)/sd(y))^2
  return(qf(0.055/2, length(x)-1,length(y)-1) > Fvalue|  Fvalue > qf(1-0.055/2,length(x)-1,length(y)-1))
}

m <- c(10, 100, 1000)
for (i in seq(1,3,1)) {
power.c5t <- mean(replicate(1000, expr = {
  x <- rnorm(m[i], 0, 1)
  y <- rnorm(m[i], 0, 1.5)
  c5t(x,y)
}))

power.Ft <- mean(replicate(1000, expr = {
  x <- rnorm(m[i], 0, 1)
  y <- rnorm(m[i], 0, 1.5)
  Ft(x,y)
}))

print(c(power.c5t,power.Ft))
}

## -----------------------------------------------------------------------------
set.seed(123)
# mardia test  function
mt.s <- function(x) {
  r <- nrow(x)
  var.inver <- solve(var(x))
  stastic <- sum((x%*%var.inver%*%t(x))^3)/r^2
}

n <- c(10, 20, 30, 50, 100, 500)
p.reject <- numeric(length(n))
m <- 100

for (i in 1:length(n)) {
  mts <- numeric(m)
  for (j in 1:m) {
    x <- as.matrix(rnorm(n[i]))
    mts[j] <- as.integer(abs(mt.s(x)) >= qchisq(0.95,1))
  }
  p.reject[i] <- mean(mts)
}
p.reject

## -----------------------------------------------------------------------------
set.seed(123)
mt.s <- function(x) {
  r <- nrow(x)
  var.inver <- solve(var(x))
  stastic <- sum((x%*%var.inver%*%t(x))^3)/r^2
}

alpha <- .1
n <- 30
m <- 2500
epsilon <- c(seq(0, .15, .01), seq(.15, 1, .05))
N <- length(epsilon)
pwr <- numeric(N)
#critical value for the skewness test

for (j in 1:N) { #for each epsilon
  e <- epsilon[j]
  sktests <- numeric(m)
  for (i in 1:m) { #for each replicate
    sigma <- sample(c(1, 10), replace = TRUE,size = n, prob = c(1-e, e))
    x <- matrix(rnorm(n, 0, sigma))
    sktests[i] <- as.integer(mt.s(x) >= qchisq(0.9,1))
  }
  pwr[j] <- mean(sktests)
}
#plot power vs epsilon
plot(epsilon, pwr, type = "b",xlab = bquote(epsilon), ylim = c(0,1))

## -----------------------------------------------------------------------------
set.seed(1234)
library(bootstrap)
n <- nrow(law)
theta.hat <- cor(law$LSAT, law$GPA)
theta <- numeric(n)
for (i in 1:n) {
  LSAT <- law$LSAT[-i]
  GPA <- law$GPA[-i]
  theta[i] <- cor(LSAT, GPA)
}
theta.bar <- mean(theta)
bias <- (n-1) * (theta.bar - theta.hat)
sde <- (n-1)/(sqrt(n)) * sd(theta)
print(c(bias, sde))

## -----------------------------------------------------------------------------
set.seed(1234)
library(boot)
boot.theta <- function(x,i) mean(x[i])
R <- aircondit$hours
de <- boot(data = R, statistic = boot.theta, R = 999)
ci <- boot.ci(de,type = c('norm','basic','perc','bca'))
ci

## -----------------------------------------------------------------------------
set.seed(1234)
library(bootstrap)
n <- nrow(scor)
lambda <- eigen(cov(scor))$values
theta <- lambda[1] / sum(lambda)
theta.hat <- numeric(n)

for (i in 1:n) {
  x <- scor[-i,]
  lambda.hat <- eigen(cov(x))$values
  theta.hat[i] <- lambda.hat[1] / sum(lambda.hat)
}

bias <- (n-1) * (mean(theta.hat) - theta)
sde <- (n-1)/(sqrt(n)) * sd(theta.hat)

print(c(bias, sde))

## ----message = FALSE----------------------------------------------------------
set.seed(1)
library(DAAG)
library(boot)
attach(ironslag)
n <- length(magnetic)
e1 <- e2 <- e3 <- e4 <- numeric((n-1)*n)
number <- 1
for (i in 1:n) {
  for (j in 1:n) {
    if (j != i) {
      y <- magnetic[c(-i,-j)]
      x <- chemical[c(-i,-j)]
      
      J1 <- lm(y ~ x)
      yhat1 <- J1$coef[1] + J1$coef[2] * chemical[i]
      yhat2 <- J1$coef[1] + J1$coef[2] * chemical[j]
      e1[number] <- magnetic[i] - yhat1
      e1[number + 1] <- magnetic[j] - yhat2
      
      J2 <- lm(y ~ x + I(x^2))
      yhat1 <- (J2$coef[1] + J2$coef[2] * chemical[i] + J2$coef[3] * chemical[i]^2)
      yhat2 <- (J2$coef[1] + J2$coef[2] * chemical[j] + J2$coef[3] * chemical[j]^2)
      e2[number] <- magnetic[i] - yhat1
      e2[number + 1] <- magnetic[j] - yhat2
      
      J3 <- lm(log(y) ~ x)
      logyhat1 <- J3$coef[1] + J3$coef[2] * chemical[i]
      yhat1 <- exp(logyhat1)
      logyhat2 <- J3$coef[1] + J3$coef[2] * chemical[j]
      yhat2 <- exp(logyhat2)
      e3[number] <- magnetic[i] - yhat1
      e3[number + 1] <- magnetic[j] - yhat2
      
      J4 <- lm(log(y) ~ log(x))
      logyhat1 <- J4$coef[1] + J4$coef[2] * log(chemical[i])
      yhat1 <- exp(logyhat1)
      logyhat2 <- J4$coef[1] + J4$coef[2] * log(chemical[j])
      yhat2 <- exp(logyhat2)
      e4[number] <- magnetic[i] - yhat1
      e4[number + 1] <- magnetic[j] - yhat2

      number = number + 1
    }
  }
}

c(mean(e1^2), mean(e2^2), mean(e3^2), mean(e4^2))


L2 <- lm(magnetic ~ chemical + I(chemical^2))




## -----------------------------------------------------------------------------
set.seed(123)


CFT <- function(x,y) {
  xbar <- x - mean(x)
  ybar <- y - mean(y)
  X.out <- sum(xbar > max(ybar)) + sum(xbar < min(ybar))
  Y.out <- sum(ybar > max(xbar)) + sum(ybar < min(xbar))
  
  return(as.integer(max(c(X.out, Y.out)) > 5))
}

CFTP <- function(x,y,r) {
  z <- c(x,y)
  output <- numeric(r)
  n <- length(z)
  for (i in 1:r) {
    index <- sample(1:n, n, replace = FALSE)
    index1 <- index[1 : (n/2)]
    index2 <- index[(n/2+1) : n]
    output[i] <- CFT(z[index1], z[index2])
  }
  mean(output)
}

n1 <- 20
n2 <- 50
mean1 <- 0
mean2 <- 0
sigm <- 1
R <- 100

count5test <- mean(replicate(R,expr = {
  x <- rnorm(n1, mean1, sigm)
  y <- rnorm(n2, mean2, sigm)
  CFT(x,y)
}))

count5test.p <- mean(replicate(R, expr = {
  x <- rnorm(n1, mean1, sigm)
  y <- rnorm(n2, mean2, sigm)
  CFTP(x, y, 100)
}) )

print(c(count5test, count5test.p))

## -----------------------------------------------------------------------------
library(RANN)
library(boot)
library(energy)
library(Ball)
library(mvtnorm)
set.seed(123)

Tn <- function(z, ix, sizes,k) {
  n1 <- sizes[1]; n2 <- sizes[2]; n <- n1 + n2
  if(is.vector(z)) z <- data.frame(z,0);
  z <- z[ix, ];
  NN <- nn2(data=z, k=k+1) # what's the first column?
  block1 <- NN$nn.idx[1:n1,-1] 
  block2 <- NN$nn.idx[(n1+1):n,-1] 
  i1 <- sum(block1 < n1 + .5); i2 <- sum(block2 > n1+.5) 
  (i1 + i2) / (k * n)
}

m <- 100; k<-3; p<-2; mu <- 0.3; 
n1 <- n2 <- 20; R<-999; n <- n1+n2; N = c(n1,n2)
eqdist.nn <- function(z,sizes,k){
  boot.obj <- boot(data=z,statistic=Tn,R=R,
  sim = "permutation", sizes = sizes,k=k)
  ts <- c(boot.obj$t0,boot.obj$t)
  p.value <- mean(ts>=ts[1])
  list(statistic=ts[1],p.value=p.value)
}

p.values <- matrix(NA,m,3)

#Unequal variances and equal expectations  
for(i in 1:m){
  x <- rmvnorm(20, mean = c(0,0), sigma = matrix(c(1,0,0,1),ncol = 2))
  y <- rmvnorm(20, mean = c(0,0), sigma = matrix(c(3,0,0,3),ncol = 2))
  z <- rbind(x,y)
  p.values[i,1] <- eqdist.nn(z,N,k)$p.value
  p.values[i,2] <- eqdist.etest(z,sizes=N,R=R)$p.value
  p.values[i,3] <- bd.test(x=x,y=y,num.permutations=999,seed=i*12345)$p.value
}

alpha <- 0.1; 
pow <- colMeans(p.values<alpha)

pow

#Unequal variances and unequal expectations
for(i in 1:m){
  x <- rmvnorm(20, mean = c(0,0), sigma = matrix(c(1,0,0,1),ncol = 2))
  y <- rmvnorm(20, mean = c(0.6,0.6), sigma = matrix(c(2,0,0,2),ncol = 2))
  z <- rbind(x,y)
  p.values[i,1] <- eqdist.nn(z,N,k)$p.value
  p.values[i,2] <- eqdist.etest(z,sizes=N,R=R)$p.value
  p.values[i,3] <- bd.test(x=x,y=y,num.permutations=999,seed=i*12345)$p.value
}

alpha <- 0.1; 
pow <- colMeans(p.values<alpha)

pow
#t distribution with 1 df and bimodel distribution  
for(i in 1:m){
  x <- rmvt(20, sigma = diag(2), df = 1)
  y1 <- rt(20,3); y2 <- rt(20,5); y <- cbind(y1, y2)
  z <- rbind(x,y)
  p.values[i,1] <- eqdist.nn(z,N,k)$p.value
  p.values[i,2] <- eqdist.etest(z,sizes=N,R=R)$p.value
  p.values[i,3] <- bd.test(x=x,y=y,num.permutations=999,seed=i*12345)$p.value
}

alpha <- 0.1; 
pow <- colMeans(p.values<alpha)

pow

#Unbalanced samples
for(i in 1:m){
  x <- rmvnorm(30, mean = c(0,0), sigma = matrix(c(1,0,0,1),ncol = 2))
  y <- rmvnorm(10, mean = c(1,1), sigma = matrix(c(2,0,0,2),ncol = 2))
  z <- rbind(x,y) 
  p.values[i,1] <- eqdist.nn(z,N,k)$p.value
  p.values[i,2] <- eqdist.etest(z,sizes=N,R=R)$p.value
  p.values[i,3] <- bd.test(x=x,y=y,num.permutations=999,seed=i*12345)$p.value
}

alpha <- 0.1; 
pow <- colMeans(p.values<alpha)

pow

## -----------------------------------------------------------------------------
set.seed(1)


Lf <- function(x) {
  exp(-abs(x))
}

rwm <- function(sigma) {
  x <- numeric(2000)
  x[1] <- 25
  u <- runif(2000)
  k <- 0
  for (i in 2:2000) {
    y <- rnorm(1, x[i-1], sigma)
    if (u[i] <= (Lf(y) / Lf(x[i-1]))) x[i] = y
    else {
      x[i] = x[i-1]
      k = k + 1
    }
  }
  list(x = x, k = k)
}

sigma <- c(0.05, 0.5, 2, 16)
rw1 <- rwm(sigma[1])
rw2 <- rwm(sigma[2])
rw3 <- rwm(sigma[3])
rw4 <- rwm(sigma[4])

#number of candidate points rejected
    ac.rate <- data.frame(sigma=sigma,no.reject=c((2000-rw1$k)/2000, (2000-rw2$k)/2000, (2000-rw3$k)/2000, (2000-rw4$k)/2000))



    rw <- cbind(rw1$x, rw2$x, rw3$x,  rw4$x)
    

## -----------------------------------------------------------------------------


    Lf <- function(x) {
       exp(-abs(x))
      }
    Gelman.Rubin <- function(psi) {
        # psi[i,j] is the statistic psi(X[i,1:j])
        # for chain in i-th row of X
        psi <- as.matrix(psi)
        n <- ncol(psi)
        k <- nrow(psi)

        psi.means <- rowMeans(psi)     #row means
        B <- n * var(psi.means)        #between variance est.
        psi.w <- apply(psi, 1, "var")  #within variances
        W <- mean(psi.w)               #within est.
        v.hat <- W*(n-1)/n + (B/n)     #upper variance est.
        r.hat <- v.hat / W             #G-R statistic
        return(r.hat)
        }

    normal.chain <- function(sigma, N, X1) {
        #generates a Metropolis chain for Normal(0,1)
        #with Normal(X[t], sigma) proposal distribution
        #and starting value X1
        x <- rep(0, N)
        x[1] <- X1
        u <- runif(N)

        for (i in 2:N) {
            xt <- x[i-1]
            y <- rnorm(1, xt, sigma)     #candidate point
            r1 <- Lf(y) * dnorm(xt, y, sigma)
            r2 <- Lf(xt) * dnorm(y, xt, sigma)
            r <- r1 / r2
            if (u[i] <= r) x[i] <- y else
                 x[i] <- xt
            }
        return(x)
        }

    sigma <- 1     #parameter of proposal distribution
    k <- 4          #number of chains to generate
    n <- 15000      #length of chains
    b <- 1000       #burn-in length

    #choose overdispersed initial values
    x0 <- c(-10, -5, 5, 10)

    #generate the chains
    set.seed(12345)
    X <- matrix(0, nrow=k, ncol=n)
    for (i in 1:k)
        X[i, ] <- normal.chain(sigma, n, x0[i])

    #compute diagnostic statistics
    psi <- t(apply(X, 1, cumsum))
    for (i in 1:nrow(psi))
        psi[i,] <- psi[i,] / (1:ncol(psi))

    #plot psi for the four chains
   
    for (i in 1:k)
        plot(psi[i, (b+1):n], type="l",
            xlab=i, ylab=bquote(psi))
    
    
    
    
    
    rhat <- rep(0, n)
    for (j in (b+1):n)
        rhat[j] <- Gelman.Rubin(psi[,1:j])
    

## -----------------------------------------------------------------------------
set.seed(1)
S <- function(a,k) {
  constant = sqrt((a^2)*k/(k+1-a^2))
  pt(constant, df = k)
}

root <- function(k) {
  x <- uniroot(function(a){(S(a,k)-S(a,k-1))}, c(0.05,sqrt(k)-0.05))
  x$root
}

x <- c(4:25,100,500,1000)
y <- numeric(25)
for (i in 4:25) {
  y[i-3] <- root(i)
}
y[23] <- root(100)
y[24] <- root(500)
y[25] <- root(1000)

data.frame(degree = x, intersection = y)

## -----------------------------------------------------------------------------
set.seed(1)
library(nloptr)

nA <- 444
nB <- 132
nOO <- 361
nAB <- 63
f <- function(x, y) {
  r1 <- 1-sum(y)
  nAA <- nA*y[1]^2/(y[1]^2+2*y[1]*r1)
  nBB <- nB*y[2]^2/(y[2]^2+2*y[2]*r1)
  r <- 1-sum(x)
  return(-2*nAA*log(x[1])-2*nBB*log(x[2])-2*nOO*log(r)-(nA-nAA)*log(2*x[1]*r)-(nB-nBB)*log(2*x[2]*r)-nAB*log(2*x[1]*x[2]))
}

g <- function(x,y) {
  return(sum(x)-0.999999)
}

opts <- list("algorithm"="NLOPT_LN_COBYLA",
             "xtol_rel"=1.0e-8)
mle <- NULL
r <- rbind(matrix(0,1,2),c(0.5,0.5))
i <- 2
while (sum(abs(r[i,] - r[i-1,])) > 1e-6) {
      res <- nloptr( x0=c(0.2,0.2),eval_f=f,lb = c(0,0), ub = c(1,1), eval_g_ineq = g, opts = opts, y=r[i,] )
      i <- i+1
      r <- rbind(r,res$solution)
      mle <- c(mle,f(x=r[i,],y=r[i-1,]))
}

#the result of EM algorithm
r

#the max likelihood values
mle


## -----------------------------------------------------------------------------
set.seed(1)
attach(mtcars)

formulas <- list(
  mpg ~ disp,
  mpg ~ I(1 / disp),
  mpg ~ disp + wt,
  mpg ~ I(1 / disp) + wt
)

# lapply

lapply(formulas, function(x) lm(x, mtcars))

# loops

out <- vector("list", length(formulas))
for (i in seq_along(formulas)) {
  out[[i]] <- lm(formulas[[i]], mtcars)
}

out


## -----------------------------------------------------------------------------
set.seed(1)
trials <- replicate(
  100,
  t.test(rpois(10,10), rpois(7,10)),
  simplify = FALSE
)

sapply(trials,  function(x) x[["p.value"]])


## -----------------------------------------------------------------------------
set.seed(1)
attach(mtcars)

#lapply()
lapply(mtcars, is.numeric)

#vapply() + Map()
f <- function(x) {x}
Map(f,vapply(mtcars, is.numeric, logical(1)))

## -----------------------------------------------------------------------------
library(Rcpp)
library(microbenchmark)
# R function
lap_f = function(x) exp(-abs(x))

rw.Metropolis = function(sigma, x0, N){
 x = numeric(N)
 x[1] = x0
 u = runif(N)
 k = 0
 for (i in 2:N) {
  y = rnorm(1, x[i-1], sigma)
  if (u[i] <= (lap_f(y) / lap_f(x[i-1]))) x[i] = y 
  else {
  x[i] = x[i-1]
  k = k+1
  }
 }
 return(list(x = x, k = k))
}

# Cpp function

cppFunction('NumericVector rwMCpp(double sigma,double x0,int N) {
  NumericVector x(N);
  x[0]=x0;
  NumericVector u(N);
  u=runif(N);
  for(int i=1;i<N;i++){
    NumericVector y(1);
    y=rnorm(1, x[i-1], sigma);
    if(u[i] <= ((exp(-abs(y[0]))) / (exp(-abs(x[i-1]))))){
      x[i] = y[0];
    }
    else{
      x[i] = x[i-1];
    }
  }
  return x;
}')


# Compare by qqplot
N = 2000
sigma = c(.05, .5, 2, 16)
x0 = 25
for (i in 1:length(sigma)) {
    x = rw.Metropolis(sigma[i],x0,N)$x
    y = rwMCpp(sigma[i],x0,N)
    qqplot(x,y)
}


# compare the time

for (i in 1:length(sigma)) {
  time <- microbenchmark(rw.Metropolis(sigma[i],x0,N),rwMCpp(sigma[i],x0,N))
  print(summary(time)[,c(1,3,5,6)])
}


