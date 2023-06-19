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
library("rjags")
library(R2jags)
library(coda)

df<- read_excel('data_surv.xlsx')

cens<-matrix(c(df$lenght_stay,rep(NA,length(df$lenght_stay))),nrow = length(df$lenght_stay), ncol = 2)
df$lenght_stay[df$state == 0] <- NA
is.censored <- as.numeric(is.na(df$lenght_stay))

X<- model.matrix(~DM  + Hemoglobin + Neutrofils + fibrinogen + Ddimer + 
                   LRA + INR + pH + pCO2 + HCO3 + FiO2 + FC + Troponine,data=df)


d.jags <- list(n = nrow(df), time = log(df$lenght_stay), cens = cens, X = X, is.censored = is.censored, Nbetas = ncol(X))
i.jags <- function(){ list(beta = rnorm(ncol(X),0,50), tau = runif(1)) }
p.jags <- c("beta", "tau")



modelo4<-function(){
  for(i in 1:n){
    is.censored[i] ~dinterval(time[i],cens[i,1])
    time[i]~dlogis(mu[i],tau)
    mu[i]<-inprod(beta[],X[i,])	
  }
  
  
  for(l in 1:Nbetas){ beta[l]~dnorm(0,50)}
  tau~dunif(0,10)
  
}

set.seed(123)
mod4 <-  jags(d.jags, i.jags, p.jags, n.chains=3, n.iter=50000,n.burnin=10000, model.file=modelo4)

mod4
ll_curve <- function(t){
  beta<-mod4$BUGSoutput$mean$beta
  beta<-as.array(beta)
  y_hat=mean(X%*%beta)
  scale=as.numeric(mod4$BUGSoutput$mean$tau)
  S_i = (1+exp((log(t)-y_hat)/scale))^(-1)
  return(S_i)
}


var_cont<-c("Hemoglobin","Neutrofils","fibrinogen","Ddimer","INR","pH","pCO2","HCO3","FiO2","FC","Troponine")



beta<-mod4$BUGSoutput$mean$beta
beta<-as.array(beta)
ypred<-X%*%beta
scale_1<-as.numeric(mod4$BUGSoutput$mean$tau)
res<- (log(df$lenght_stay)-ypred)/scale_1
res
exp(ypred
res<- data.frame(lp=ypred,res=res)
res

pllogis=fitdist(res$res, "logis")
pllogis
pv2<-round(ks.test (res$res, "plogis",location = pllogis$estimate[1],scale = pllogis$estimate[2])$p.value,2)

pv2



ggplot(data=res, aes(x=res)) +
  geom_histogram(aes(y=..density..), breaks=seq(-4,4,by=0.7), 
                 colour="Black", fill="#2171B5") + 
  ylim(0,0.75)+
  xlim(-4,4) +
  xlab("residuals") +
  ylab("Density") +
  theme_minimal() +
  geom_function(aes(colour=paste("Logistic","-",pv2)), fun = dlogis, 
                args = list(pllogis$estimate[1], pllogis$estimate[2]), size=1.2)+
  scale_color_manual(values=c("red2"))+
  labs(colour="Distribuição")

ggplot(data=res, aes(x=lp,y=res)) +
  geom_point(size=2, colour="#2171B5")+
  theme_minimal()+
  xlab("Linear predictor") +
  ylab("residuals") 


ggplot(data=res, aes(sample=res)) +
  theme_minimal()+
  stat_qq(distribution = qlogis)+
  stat_qq_line(distribution = qlogis)+
  xlab("Theoretical") +
  ylab("Sample")



ll_curve_factor <- function(t){
  x_cont=colMeans(df[,var_cont])
  x_cat=c(0,0)
  b0<-mod4$BUGSoutput$mean$beta[1]
  b_cat<-as.array(mod4$BUGSoutput$mean$beta[c(2,7)])
  b_cont<-as.array(mod4$BUGSoutput$mean$beta[-c(1,2,7)])
  y_hat=as.numeric(b0+b_cat%*%x_cat+b_cont%*%x_cont)
  scale=as.numeric(mod4$BUGSoutput$mean$tau)
  S_i = (1+exp((log(t)-y_hat)/scale))^(-1)
  return(S_i)
}


ll_curve_factor_1 <- function(t){
  x_cont=colMeans(df[,var_cont])
  x_cat=c(1,0)
  b0<-mod4$BUGSoutput$mean$beta[1]
  b_cat<-mod4$BUGSoutput$mean$beta[c(2,7)]
  b_cont<-mod4$BUGSoutput$mean$beta[-c(1,2,7)]
  y_hat=as.numeric(b0+b_cat%*%x_cat+b_cont%*%x_cont)
  scale=as.numeric(mod4$BUGSoutput$mean$tau)
  S_i = (1+exp((log(t)-y_hat)/scale))^(-1)
  return(S_i)
}


ll_curve_factor_2 <- function(t){
  x_cont=colMeans(df[,var_cont])
  x_cat=c(0,1)
  b0<-mod4$BUGSoutput$mean$beta[1]
  b_cat<-as.array(mod4$BUGSoutput$mean$beta[c(2,7)])
  b_cont<-as.array(mod4$BUGSoutput$mean$beta[-c(1,2,7)])
  y_hat=as.numeric(b0+b_cat%*%x_cat+b_cont%*%x_cont)
  scale=as.numeric(mod4$BUGSoutput$mean$tau)
  S_i = (1+exp((log(t)-y_hat)/scale))^(-1)
  return(S_i)
}

ll_curve_factor_3 <- function(t){
  x_cont=colMeans(df[,var_cont])
  x_cat=c(1,1)
  b0<-mod4$BUGSoutput$mean$beta[1]
  b_cat<-as.array(mod4$BUGSoutput$mean$beta[c(2,7)])
  b_cont<-as.array(mod4$BUGSoutput$mean$beta[-c(1,2,7)])
  y_hat=as.numeric(b0+b_cat%*%x_cat+b_cont%*%x_cont)
  scale=as.numeric(mod4$BUGSoutput$mean$tau)
  S_i = (1+exp((log(t)-y_hat)/scale))^(-1)
  return(S_i)
}



curve(ll_curve_factor(x),from=0,to=60,col=2, lwd=2)
curve(ll_curve_factor_1(x), add = TRUE, col = 3, lwd=2)
curve(ll_curve_factor_2(x), add = TRUE, col = "yellow", lwd=2)
curve(ll_curve_factor_3(x), add = TRUE, col = "blue", lwd=2)
curve(ll_curve, add = TRUE, col = "black", lwd=2)
legend(x = "topright",          # Position
       legend = c("DM=0|LRA=0", "DM=1|LRA=0","DM=0|LRA=1","DM=1|LRA=1","Mean"),  # Legend texts
       lty = c(1, 1,1,1,1),           # Line types
       col = c(2, 3,"yellow","blue","black"),
       lwd = 2) 
