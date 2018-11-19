####Simulation on power function(only row fixed)
set.seed(123)
p1 <- seq(0.01, 0.99, 0.01)
p2 <- sample(p1, length(p1), replace=F)
n1 <- 15
n2 <- 15
alpha <- 0.005
B <- 1000

WaldStat <- function(x){
  z <- sum(x[,1])
  O11 <- x[1,1]
  N <- sum(x)
  n <- N/2
  p1 <- O11/n
  p2 <- (z-O11)/n
  return(log(p1/(1-p1))-log(p2/(1-p2)))
}

chisq.stat <- function(x){
  return(chisq.test(x)$statistic)
}

WaldRej <- function(x, n1, n2, alpha){
  O11 <- x[,1]
  O21 <- x[,2]
  p1 <- O11/n1
  p2 <- O21/n2
  stat <- log(p1/(1-p1))-log(p2/(1-p2))
  sigma <- 1/O11 + 1/O21 + 1/(n1-O11) + 1/(n2-O21)
  return(abs((stat<qnorm(alpha/2, 0, sigma))|(stat>qnorm(1-alpha/2, 0, sigma))))
}

WaldRejAd <- function(x, n1, n2, alpha){
  O11 <- x[,1]+0.5
  O21 <- x[,2]+0.5
  p1 <- O11/(n1+1)
  p2 <- O21/(n2+1)
  stat <- log(p1/(1-p1))-log(p2/(1-p2))
  sigma <- sqrt(1/O11 + 1/O21 + 1/(n1+1-O11) + 1/(n2+1-O21))
  return(abs(((stat-0.5)<qnorm(alpha/2, 0, sigma))|((stat+0.5)>qnorm(1-alpha/2, 0, sigma))))
}

ChiqRej <- function(x, n1, n2, alpha){
  O11 <- x[,1]
  O21 <- x[,2]
  X <- array(0, dim=c(2,2,length(O11)))
X[1,1,] <- O11
X[2,1,] <- O21
X[1,2,] <- n1-O11
X[2,2,] <- n2-O21
  stat <- apply(X, 3, chisq.stat)
  return(abs((stat<qchisq(alpha/2,1))|(stat>qchisq(1-alpha/2,1))))
}

prop.stat <- function(x){
  return(prop.test(x)$p.value)
} 

PropRej <- function(x, n1, n2, alpha){
  O11 <- x[,1]
  O21 <- x[,2]
  X <- array(0, dim=c(2,2,length(O11)))
X[1,1,] <- O11
X[2,1,] <- O21
X[1,2,] <- n1-O11
X[2,2,] <- n2-O21
  stat <- apply(X, 3, prop.stat)
  return(abs((stat<qchisq(alpha/2,1))|(stat>qchisq(1-alpha/2,1))))
}

G.stat <- function(x){
  library(DescTools)
  return(GTest(x)$statistic)
}

GRej <- function(x, n1, n2, alpha){
  O11 <- x[,1]
  O21 <- x[,2]
  X <- array(0, dim=c(2,2,length(O11)))
X[1,1,] <- O11
X[2,1,] <- O21
X[1,2,] <- n1-O11
X[2,2,] <- n2-O21
stat <- apply(X, 3, G.stat)
  return(abs((stat<qchisq(alpha/2,1))|(stat>qchisq(1-alpha/2,1))))  
}

###generate 2X2 table
O11 <- apply(matrix(p1,nrow=1), 2, rbinom, n=B, size=n1)
O21 <- apply(matrix(p2,nrow=1), 2, rbinom, n=B, size=n2)
X <- array(0, dim=c(B,2,length(p1)))
for(i in 1:length(p1)){
  X[,,i] <- cbind(O11[,i],O21[,i])
}
Wald <- na.omit(apply(X,3,WaldRej, n1=n1, n2=n2, alpha=alpha))
WaldAd <- na.omit(apply(X,3,WaldRejAd, n1=n1, n2=n2, alpha=alpha))
Chisq <- na.omit(apply(X,3,ChiqRej, n1=n1, n2=n2, alpha=alpha))
Score <- na.omit(apply(X,3,PropRej, n1=n1, n2=n2, alpha=alpha))
LRT <- na.omit(apply(X,3,GRej, n1=n1, n2=n2, alpha=alpha))
H1 <- log(p1/(1-p1))-log(p2/(1-p2))
WaldPow <- apply(Wald,2,mean)
WaldAdPow <- apply(WaldAd,2,mean)
ChisqPow <- apply(Chisq,2,mean)
ScorePow <- apply(Score,2,mean)
LRTPow <- apply(LRT,2,mean)

ord=order(H1)
plotdata <- data.frame(log_OR=rep(sort(H1),5), Power=c(WaldPow[ord],
 WaldAdPow[ord], ChisqPow[ord], ScorePow[ord], LRTPow[ord]),
 type=rep(c("Wald", "Wald Correction", "Chi-Square", "Score", "LRT"), each=length(H1)))
library(ggplot2)
library(tidyverse)
####Power function plot
plotdata %>% 
ggplot(aes(x=log_OR, y=Power, color=type))+
  geom_smooth(se=F)+
  labs(title = "Comparison of Power when row margin fixed(alpha=0.005)")+
  geom_hline(yintercept = alpha, lty=2, col=2)+
scale_colour_manual(values = c('orange','red','blue','green','purple'))

########simulation on CI
####there are wald CI, fisher CI, adWald
WaldCI <- function(x, n1, n2){
  O11 <- x[,1]
  O21 <- x[,2]
  p1 <- O11/n1
  p2 <- O21/n2
  stat <- log(p1/(1-p1))-log(p2/(1-p2))
  sigma <- sqrt(1/O11 + 1/O21 + 1/(n1-O11) + 1/(n2-O21))
  return(c(conf.inf=c(stat-1.96*sigma, stat+1.96*sigma)))
}

WaldAdCI <- function(x, n1, n2){
  O11 <- x[,1]+0.5
  O21 <- x[,2]+0.5
  p1 <- O11/(n1+1)
  p2 <- O21/(n2+1)
  stat <- log(p1/(1-p1))-log(p2/(1-p2))
  sigma <- sqrt(1/O11 + 1/O21 + 1/(n1+1-O11) + 1/(n2+1-O21))
  return(c(conf.inf=c(stat-1.96*sigma, stat+1.96*sigma)))
}

FisherCI <- function(x){
  return(log(fisher.test(x)$conf.int))
}

ExactCI <- function(x, n1, n2){
  X <- array(0, dim=c(2,2,nrow(x)))
  X[1,1,] <- x[,1]
  X[2,1,] <- x[,2]
  X[1,2,] <- n1-x[,1]
  X[2,2,] <- n2-x[,2]
  return(apply(X, 3, FisherCI))
}

score <- function(x, n1, n2){
  library(ratesci)
  O11 <- x[1]
  O21 <- x[2] 
  return(log(scoreci(O11,n1,O21,n2,contrast="OR")$estimates[c(1,3)]))
}

ScoreCI <- function(x, n1, n2){
  return(apply(x,1,score, n1=n1, n2=n2))
}

waldCI <- apply(X, 3, WaldCI, n1=n1, n2=n2)
waldAdCI <- apply(X, 3, WaldAdCI, n1=n1, n2=n2)
exactCI <- apply(X, 3, ExactCI, n1=n1, n2=n2)
scoreCI <- apply(X, 3, ScoreCI, n1=n1, n2=n2)

waldCIAr <- array(0, dim=c(2, length(p1), B))
waldAdCIAr <- array(0, dim=c(2, length(p1), B))
exactCIAr <- array(0, dim=c(2, length(p1), B))
scoreCIAr <- array(0, dim=c(2, length(p1), B))

for(i in 1:B){
 waldCIAr[,,i] <- waldCI[c(i,(1000+i)),]
 waldAdCIAr[,,i] <- waldAdCI[c(i,(1000+i)),]
 exactCIAr[,,i] <- exactCI[c(((i-1)*2+1):(2*i)),]
 scoreCIAr[,,i] <- scoreCI[c(((i-1)*2+1):(2*i)),]
}
logOR <- log(p1/(1-p1))-log(p2/(1-p2))

Coverage <- function(x, logOR=logOR){
  lower <- x[1,]
  upper <- x[2,]
  return(abs((lower<=logOR)&(upper>=logOR)))
}
waldCIcontain <- apply(waldCIAr,3,Coverage,logOR=logOR)
waldCIcontain[is.na(waldCIcontain)] <- 0
waldAdCIcontain <- apply(waldAdCIAr,3,Coverage,logOR=logOR)

exactCIcontain <- apply(exactCIAr,3,Coverage,logOR=logOR)

scoreCIcontain <- apply(scoreCIAr,3,Coverage,logOR=logOR)

waldCoProb <- apply(waldCIcontain, 1, mean)
waldAdCoProb <- apply(waldAdCIcontain, 1, mean)
exactCoProb <- apply(exactCIcontain, 1, mean)
scoreCoProb <- apply(scoreCIcontain, 1, mean)

ord=order(H1)
plotdata1 <- data.frame(log_OR=rep(sort(H1),4), Coverage_Probability=c(waldCoProb[ord],
 waldAdCoProb[ord], exactCoProb[ord], scoreCoProb[ord]),
 type=rep(c("Wald", "Wald Correction", "Exact", "Score"), each=length(H1)))
####Power function plot
plotdata1 %>% 
ggplot(aes(x=log_OR, y=Coverage_Probability, color=type))+
  geom_line()+
  labs(title = "Comparison of Coverage Probability for different 95% confidence intervals")+
  geom_hline(yintercept = 0.95, lty=2, col=2)+
scale_colour_manual(values = c('orange','red','blue','green'))

Width <- function(x){
  lower <- x[1,]
  upper <- x[2,]
  return(abs(upper-lower))
} 

Mean <- function(x){
  y <- x[x<Inf]
  return(mean(y))
}

waldCIwid <- apply(waldCIAr,3,Width)
waldCIwid[is.na(waldCIwid)] <- 0
waldAdCIwid <- apply(waldAdCIAr,3,Width)
exactCIwid <- apply(exactCIAr,3,Width)
scoreCIwid <- apply(scoreCIAr,3,Width)

waldAveWid <- apply(waldCIwid,1,Mean)
waldAdAveWid <- apply(waldAdCIwid,1,Mean)
exactAveWid <- apply(exactCIwid,1,Mean)
scoreAveWid <- apply(scoreCIwid,1,Mean)

ord=order(H1)
plotdata2 <- data.frame(log_OR=rep(sort(H1),4), Average_Width=c(waldAveWid[ord],
 waldAdAveWid[ord], exactAveWid[ord], scoreAveWid[ord]),
 type=rep(c("Wald", "Wald Correction", "Exact", "Score"), each=length(H1)))
####Power function plot
plotdata2%>% 
ggplot(aes(x=log_OR, y=Average_Width, color=type))+
  geom_smooth(se=F)+
  labs(title = "Comparison of Average Width for different 95% confidence intervals")+
  geom_hline(yintercept = 0, lty=2, col=2)+
scale_colour_manual(values = c('orange','red','blue','green'))



