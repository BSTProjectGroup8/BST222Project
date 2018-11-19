#########################################
####            Setup                ####
#########################################

lib <- c("DescTools", "statmod", "pwr", "SuppDists")
install.packages(lib)
lapply(lib, require, character.only = TRUE)

#########################################
####          Hypothesis Test        ####
#########################################

# 1. Wald test
WaldStat <- function(z, O11, N){
  n <- N/2
  p1 <- O11/n
  p2 <- (z-O11)/n
  return(log(p1/(1-p1))-log(p2/(1-p2)))
}

# 2. Chi-square
ChiStat <- function(z,O11,N){
  n <- N/2
  x <- matrix(c(O11, z-O11, n-O11, n-z+O11),2,2)
  return(chisq.test(x)$statistic)
}

# 3. LRT
LRTstat <- function(z,O11,N){
  n <- N/2
  x <- matrix(c(O11, z-O11, n-O11, n-z+O11),2,2)
  return(GTest(x)$statistic)
}

# 4. Exact
Exactstat <- function(z,O11,N){
  n <- N/2
  x <- matrix(c(O11, z-O11, n-O11, n-z+O11),2,2)
  return(fisher.test(x)$estimate)
}

alpha <- 0.05
N <- 30   # the total number of subjects; 15 to case & 15 to control
z <- 10    # the total number of outcome events in the study
O11 <- 6    # the value in the first cell.

WaldStat(z,O11,N)
ChiStat(z,O11,N)
LRTstat(z,O11,N)
Exactstat(z,O11,N)

#########################################
####            Power                ####
#########################################

P <- data.frame(TreaEven=c(0:z), contrEven=z-c(0:z))
## P is a data frame containing all possible 2 by 2 table with given z
# TreaEven column contains the possible value for cell 11
# contrEven column contains the possible value for cell 21
# H1 column contains the log odds ratio under H1

P$H1 <- log((P$TreaEven/(N/2))/(1-P$TreaEven/(N/2)))-log((P$contrEven/(N/2))/(1-P$contrEven/(N/2)))

# Wald test
WaldPower <- function(z,O11,N,H1){
  mu <- H1
  sigma <- 1/O11 + 1/(z-O11) + 1/(N/2-O11) + 1/(N/2-z+O11)
  return(pnorm(alpha/2, mu, sigma)+1-pnorm(1-alpha/2,mu,sigma))
}

WaldPower(z,O11,N,H1=P$H1)

# Chi-squate
ChiPower <- function(z,O11,N,P){
  chi <- ((O11-P[,1])^2)/P[,1]+((O11-P[,1])^2)/(15-P[,1])+((z-O11-P[,2])^2)/(P[,2])+((z-O11-P[,2])^2)/(15-P[,2])
  return(pwr.chisq.test(sqrt(chi/N), N, 1, sig.level = alpha))
}
ChiPower(z,O11,N,P)

# Exact test
s = sort(rexp(100))
range01 <- function(x, na.rm = FALSE)
{(x-min(x))/(max(x)-min(x))}
range01(s) # 100 numbers rangling between 0 and 1

p1 <- range01(s) #p1 ~ any proportion possible between 0 and 1
p2 <- 1-p1

power.fisher.test(p1, p2, 15, 15, 100, alpha = 0.05, alternative = "two.sided")
dghyper(seq(0,10,1), z, N/2, N)

# G-test to be written

#########################################
####      Confidence Interval        ####
#########################################

## CI coverage

## CI width