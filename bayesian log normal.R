library(readxl)
library(lubridate)
library(MASS)
library(vcd)
library(mltools)
library(installr)
library(survival)
library("writexl")
library(stringr)
library(pec)
library(My.stepwise)
library(devtools)
library(ggbiplot)
library(factoextra)
library(reliaR)
library(actuar)
library(fitdistrplus)
library(survival)
library(data.table)
library(goft)
library(psych)
require(GGally)
require(VGAM)
library(model4you)
library(eha)
library(glmnet)
library(flexsurv)
library(BayesSurvival)
library(rstanarm)
library(R2jags)
library(coda)

df<- read_excel('data_surv.xlsx')
cens<-matrix(c(df$lenght_stay,rep(NA,length(df$lenght_stay))),nrow = length(df$lenght_stay), ncol = 2)
df$lenght_stay[df$state == 0] <- NA
is.censored <- as.numeric(is.na(df$lenght_stay))

X<- model.matrix(~DM+Hemoglobin+Neutrofils+fibrinogen+Ddimer+LRA+INR+pH+pCO2+FiO2+FC,data=df)

d.jags <- list(n = nrow(df), time = df$lenght_stay, cens = cens, X = X, is.censored = is.censored, Nbetas = ncol(X))
i.jags <- function(){ list(beta = rnorm(ncol(X),0,50), tau = runif(1)) }
p.jags <- c("beta", "tau")

modelo1<-function(){
  for(i in 1:n){
    is.censored[i] ~dinterval(time[i],cens[i,1])
    time[i]~dlnorm(mu[i], tau)
    mu[i]<-inprod(beta[],X[i,])	
  }
  
  
  for(l in 1:Nbetas){ beta[l]~dnorm(0,50)}
  tau~dunif(0,10)
  
}

set.seed(123)
mod1 <-  jags(d.jags, i.jags, p.jags, n.chains=3, n.iter=50000,n.burnin=10000, model.file=modelo1)

beta<-mod1$BUGSoutput$mean$beta
beta<-as.array(beta)
ypred<-X%*%beta

scale<-as.numeric(mod1$BUGSoutput$mean$tau)

my_resid<- (log(df$lenght_stay)-ypred)/scale
my_resid<- data.frame(lp=ypred,my_resid=my_resid)

hist(my_resid$my_resid,freq = T)


pv1<-round(ks.test(my_resid$my_resid,pnorm,0,sd(my_resid$my_resid))$p.value,2)
pv1

ggplot(data=my_resid, aes(x=my_resid)) +
  geom_histogram(aes(y=..density..),breaks=seq(-0.3,0.3,by=0.1),
                 colour="Black", fill="#2171B5")+
  xlab("Residuals") +
  ylab("Density") +
  theme_minimal() +
  geom_function(aes(colour=paste("Normal","-",pv1)), fun = dnorm, 
                args = list(0, sd(my_resid$my_resid)), size=1.2)+
  scale_color_manual(values=c("red2"))+
  xlim(-1,1) +
  labs(colour="Distribuição")

ggplot(data=my_resid, aes(x=lp,y=my_resid)) +
  geom_point(size=2, colour="#2171B5")+
  theme_minimal()+
  xlab("Linear Predictor") +
  ylab("residuals") 

ggplot(data=my_resid, aes(sample=my_resid)) +
  theme_minimal()+
  stat_qq()+
  stat_qq_line()+
  xlab("Theoretical") +
  ylab("Sample")

var_cont<-c("Hemoglobin","Neutrofils","fibrinogen","Ddimer","INR","pH","pCO2","FiO2","FC")


ln_curve <- function(t){
  beta<-mod1$BUGSoutput$mean$beta
  beta<-as.array(beta)
  y_hat=mean(X%*%beta)
  scale=log(as.numeric(mod1$BUGSoutput$mean$tau))
  S_i = 1-pnorm((log(t)-y_hat)/scale)
  return(S_i)
}

ln_curve_factor <- function(t){
  x_cont=colMeans(df[,var_cont])
  x_cat=c(0,0)
  b0<-mod1$BUGSoutput$mean$beta[1]
  b_cat<-mod1$BUGSoutput$mean$beta[c(2,7)]
  b_cont<-mod1$BUGSoutput$mean$beta[-c(1,2,7)]
  y_hat=as.numeric(b0+b_cat%*%x_cat+b_cont%*%x_cont)
  scale=log(as.numeric(mod1$BUGSoutput$mean$tau))
  S_i = 1-pnorm((log(t)-y_hat)/scale)
  return(S_i)
}


ln_curve_factor_1 <- function(t){
  x_cont=colMeans(df[,var_cont])
  x_cat=c(1,0)
  b0<-mod1$BUGSoutput$mean$beta[1]
  b_cat<-mod1$BUGSoutput$mean$beta[c(2,7)]
  b_cont<-mod1$BUGSoutput$mean$beta[-c(1,2,7)]
  y_hat=as.numeric(b0+b_cat%*%x_cat+b_cont%*%x_cont)
  scale=log(as.numeric(mod1$BUGSoutput$mean$tau))
  S_i = 1-pnorm((log(t)-y_hat)/scale)
  return(S_i)
}

ln_curve_factor_2 <- function(t){
  x_cont=colMeans(df[,var_cont])
  x_cat=c(0,1)
  b0<-mod1$BUGSoutput$mean$beta[1]
  b_cat<-mod1$BUGSoutput$mean$beta[c(2,7)]
  b_cont<-mod1$BUGSoutput$mean$beta[-c(1,2,7)]
  y_hat=as.numeric(b0+b_cat%*%x_cat+b_cont%*%x_cont)
  scale=log(as.numeric(mod1$BUGSoutput$mean$tau))
  S_i = 1-pnorm((log(t)-y_hat)/scale)
  return(S_i)
}

ln_curve_factor_3 <- function(t){
  x_cont=colMeans(df[,var_cont])
  x_cat=c(1,1)
  b0<-mod1$BUGSoutput$mean$beta[1]
  b_cat<-mod1$BUGSoutput$mean$beta[c(2,7)]
  b_cont<-mod1$BUGSoutput$mean$beta[-c(1,2,7)]
  y_hat=as.numeric(b0+b_cat%*%x_cat+b_cont%*%x_cont)
  scale=log(as.numeric(mod1$BUGSoutput$mean$tau))
  S_i = 1-pnorm((log(t)-y_hat)/scale)
  return(S_i)
}


curve(ln_curve_factor(x),from=0,to=60,col=2, lwd=2)
curve(ln_curve_factor_1(x), add = TRUE, col = 3, lwd=2)
curve(ln_curve_factor_2(x), add = TRUE, col = "yellow", lwd=2)
curve(ln_curve_factor_3(x), add = TRUE, col = "blue", lwd=2)
curve(ln_curve, add = TRUE, col = "black", lwd=2)
legend(x = "topright",          # Position
       legend = c("DM=0|LRA=0", "DM=1|LRA=0","DM=0|LRA=1","DM=1|LRA=1","Mean"),  # Legend texts
       lty = c(1, 1,1,1,1),           # Line types
       col = c(2, 3,"yellow","blue","black"),
       lwd = 2) 
