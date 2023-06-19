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
library(BayesPostEst)

df<- read_excel('data_surv.xlsx')

cens<-matrix(c(df$lenght_stay,rep(NA,length(df$lenght_stay))),nrow = length(df$lenght_stay), ncol = 2)
df$lenght_stay[df$state == 0] <- NA
is.censored <- as.numeric(is.na(df$lenght_stay))

X<- model.matrix(~DM  + Hemoglobin + Neutrofils + fibrinogen + Ddimer + 
                   LRA + INR + pH + pCO2 + HCO3 + FiO2 + FC + Troponine,data=df)

d.jags <- list(n = nrow(df), time = log(df$lenght_stay), cens = cens, X = X, is.censored = is.censored, Nbetas = ncol(X))
i.jags <- function(){ list(beta = rnorm(ncol(X),0,50), tau =rgamma(0.001,0.001)) }
p.jags <- c("beta", "tau","sigma")

modelo4<-function(){
  for(i in 1:n){
    is.censored[i] ~dinterval(time[i],cens[i,1])
    time[i]~dlogis(mu[i],sqrt(tau))
    mu[i]<-inprod(beta[],X[i,])	
    residuos[i]<- (time[i]- mu[i])/sigma
    S_0[i] <-1/(1+ exp((log(0)-mu[i])*sqrt(tau)))
    S_1[i] <- 1/(1+ exp((log(2)-mu[i])*sqrt(tau)))
    S_2[i] <- 1/(1+ exp((log(4)-mu[i])*sqrt(tau)))
    S_3[i] <- 1/(1+ exp((log(6)-mu[i])*sqrt(tau)))
    S_4[i] <- 1/(1+ exp((log(8)-mu[i])*sqrt(tau)))
    S_5[i] <- 1/(1+ exp((log(10)-mu[i])*sqrt(tau)))
    S_6[i] <- 1/(1+ exp((log(15)-mu[i])*sqrt(tau)))
    S_7[i] <- 1/(1+ exp((log(20)-mu[i])*sqrt(tau)))
    S_8[i] <- 1/(1+ exp((log(25)-mu[i])*sqrt(tau)))
    S_9[i] <- 1/(1+ exp((log(30)-mu[i])*sqrt(tau)))
    S_10[i] <-1/(1+ exp((log(40)-mu[i])*sqrt(tau)))
    S_11[i] <- 1/(1+ exp((log(50)-mu[i])*sqrt(tau)))
  }
  
  
  for(l in 1:Nbetas){ beta[l]~dnorm(0,55)}
  tau ~ dgamma(0.001,0.001)
  sigma <- sqrt(1/tau)
  
}
mod4 <-  jags(d.jags, i.jags, p.jags, n.chains=3, n.iter=50000,n.burnin=10000, model.file=modelo4)
mcmcTab(mod4)

par(mar = rep(2, 4))
traceplot(mod4,ask=FALSE,mfrow = c(4, 5))



S0<-mean(mod4$BUGSoutput$mean$S_0)
S1<-mean(mod4$BUGSoutput$mean$S_1)
S2<-mean(mod4$BUGSoutput$mean$S_2)
S3<-mean(mod4$BUGSoutput$mean$S_3)
S4<-mean(mod4$BUGSoutput$mean$S_4)
S5<-mean(mod4$BUGSoutput$mean$S_5)
S6<-mean(mod4$BUGSoutput$mean$S_6)
S7<-mean(mod4$BUGSoutput$mean$S_7)
S8<-mean(mod4$BUGSoutput$mean$S_8)
S9<-mean(mod4$BUGSoutput$mean$S_9)
S10<-mean(mod4$BUGSoutput$mean$S_10)
S11<-mean(mod4$BUGSoutput$mean$S_11)

t<-c(0,2,4,6,8,10,15,20,25,30,40,50)

S<-c(S0,S1,S2,S3,S4,S5,S6,S7,S8,S9,S10,S11)
surv_ll<-c()
surv_ll<-data.frame(t=t,S=S)




beta<-mod4$BUGSoutput$mean$beta
beta<-as.array(beta)
ypred<-X%*%beta
res<-mod4$BUGSoutput$mean$residuos
res<- data.frame(i=ypred,res=res)

pllogis=fitdist(res$res, "logis")
pllogis
pv2<-round(ks.test (res$res, "plogis",location = pllogis$estimate[1],scale = pllogis$estimate[2])$p.value,2)

pv2




ggplot(data=res, aes(x=res)) +
  geom_histogram(aes(y=..density..), breaks=seq(-3,3,by=0.6), 
                 colour="Black", fill="#2171B5") + 
  ylim(0,1.75)+
  xlim(-3.5,3.5) +
  xlab("residuals") +
  ylab("Density") +
  theme_minimal() +
  geom_function(aes(colour=paste("Logistic","- p = ",pv2)), fun = dlogis, 
                args = list(pllogis$estimate[1], pllogis$estimate[2]), size=1.2)+
  scale_color_manual(values=c("red2"))+
  labs(colour="Distribution")

ggplot(data=res, aes(x=ypred,y=res)) +
  geom_point(size=2, colour="#2171B5")+
  theme_minimal()+
  xlab("Predicted values") +
  ylab("residuals") 


ggplot(data=res, aes(sample=res)) +
  theme_minimal()+
  stat_qq(distribution = qlogis)+
  stat_qq_line(distribution = qlogis)+
  xlab("Theoretical quantiles") +
  ylab("Sample quantiles")











X_1<- model.matrix(~DM+Hemoglobin+Neutrofils+fibrinogen+Ddimer+LRA+INR+pH+pCO2+FiO2+FC,data=df)

d.jags <- list(n = nrow(df), time = df$lenght_stay, cens = cens, X = X_1, is.censored = is.censored, Nbetas = ncol(X_1))
i.jags <- function(){ list(beta = rnorm(ncol(X_1),0,50), tau = rgamma(0.001,0.001)) }
p.jags <- c("beta", "tau","sigma")

modelo1<-function(){
  for(i in 1:n){
    is.censored[i] ~dinterval(time[i],cens[i,1])
    time[i]~dlnorm(mu[i], tau)
    mu[i]<-inprod(beta[],X[i,])	
    residuos[i]<- (log(time[i])- mu[i])/sigma
    S_0[i] <- 1- phi((log(0)-mu[i])*sqrt(tau))
    S_1[i] <- 1- phi((log(2)-mu[i])*sqrt(tau))
    S_2[i] <- 1- phi((log(4)-mu[i])*sqrt(tau))
    S_3[i] <- 1- phi((log(6)-mu[i])*sqrt(tau))
    S_4[i] <- 1- phi((log(8)-mu[i])*sqrt(tau))
    S_5[i] <- 1- phi((log(10)-mu[i])*sqrt(tau))
    S_6[i] <- 1- phi((log(15)-mu[i])*sqrt(tau))
    S_7[i] <- 1- phi((log(20)-mu[i])*sqrt(tau))
    S_8[i] <- 1- phi((log(25)-mu[i])*sqrt(tau))
    S_9[i] <- 1- phi((log(30)-mu[i])*sqrt(tau))
    S_10[i] <- 1- phi((log(40)-mu[i])*sqrt(tau))
    S_11[i] <- 1- phi((log(50)-mu[i])*sqrt(tau))
  }
  
  
  for(l in 1:Nbetas){ beta[l]~dnorm(0,55)}
  tau ~ dgamma(0.001,0.001)
  sigma <- sqrt(1/tau)
  
}



mod1 <-  jags(d.jags, i.jags, p.jags, n.chains=3, n.iter=50000,n.burnin=10000, model.file=modelo1)
mod1
mcmcTab(mod1)
par(mar = rep(2, 4))
traceplot(mod1,ask=FALSE,mfrow = c(3, 5))


S0<-mean(mod1$BUGSoutput$mean$S_0)
S1<-mean(mod1$BUGSoutput$mean$S_1)
S2<-mean(mod1$BUGSoutput$mean$S_2)
S3<-mean(mod1$BUGSoutput$mean$S_3)
S4<-mean(mod1$BUGSoutput$mean$S_4)
S5<-mean(mod1$BUGSoutput$mean$S_5)
S6<-mean(mod1$BUGSoutput$mean$S_6)
S7<-mean(mod1$BUGSoutput$mean$S_7)
S8<-mean(mod1$BUGSoutput$mean$S_8)
S9<-mean(mod1$BUGSoutput$mean$S_9)
S10<-mean(mod1$BUGSoutput$mean$S_10)
S11<-mean(mod1$BUGSoutput$mean$S_11)

t<-c(0,2,4,6,8,10,15,20,25,30,40,50)

S<-c(S0,S1,S2,S3,S4,S5,S6,S7,S8,S9,S10,S11)
surv<-c()
surv<-data.frame(t=t,S=S)

km<-survfit(Surv(df$lenght_stay,df$state)~1)
plot(surv,type = "l",lwd=2,col=2)
lines(km,conf.int=F,col ="blue",lwd=2)

beta<-mod1$BUGSoutput$mean$beta
beta<-as.array(beta)
ypred<-X_1%*%beta
res<-mod1$BUGSoutput$mean$residuos
my_resid<- data.frame(i=ypred,my_resid=res)


pv1<-round(ks.test(my_resid$my_resid,pnorm,0,1)$p.value,3)
pv1

ggplot(data=my_resid, aes(x=my_resid)) +
  geom_histogram(aes(y=..density..),breaks=seq(-3,3,by=0.7),
                 colour="Black", fill="#2171B5")+
  xlab("Residuals") +
  ylab("Density") +
  theme_minimal() +
  geom_function(aes(colour=paste("Normal","- p = ",pv1)), fun = dnorm, 
                args = list(0, sd(my_resid$my_resid)), size=1.2)+
  scale_color_manual(values=c("red2"))+
  xlim(-3,3) +
  labs(colour="Distribution")

ggplot(data=my_resid, aes(x=i,y=my_resid)) +
  geom_point(size=2, colour="#2171B5")+
  theme_minimal()+
  xlab("Predicted values") +
  ylab("residuals") 

ggplot(data=my_resid, aes(sample=my_resid)) +
  theme_minimal()+
  stat_qq()+
  stat_qq_line()+
  xlab("Theoretical quantiles") +
  ylab("Sample quantiles")








km<-survfit(Surv(df$lenght_stay,df$state)~1)
plot(surv,type = "l",lwd=2,col=2,xlab="LOS",ylab="Survival probability")
lines(surv_ll,col ="blue",lwd=2)
lines(km,conf.int=F,col ="green",lwd=2)
legend(x = "topright",          # Position
       legend = c("log-normal", "log-logistic","km"),  # Legend texts
       lty = c(1, 1,1),           # Line types
       col = c(2,"blue",3),
       lwd = 2) 




fit<- survreg(formula = Surv(lenght_stay, state) ~ DM + Hemoglobin + 
                Neutrofils + fibrinogen + Ddimer + LRA + INR + pH + pCO2 + 
                FiO2 + FC , data = df, dist = "lognormal")

ln_curve <- function(t){
  y_hat=mean(fit$linear.predictors)
  scale=fit$scale
  S_i = 1-pnorm((log(t)-y_hat)/scale)
  return(S_i)
}

fit_ll<-survreg(formula = Surv(lenght_stay, state) ~
                  DM +  Hemoglobin + Neutrofils + fibrinogen + Ddimer + 
                  LRA + INR + pH + pCO2 + HCO3 + FiO2 + FC + Troponine, data = df, 
                dist = "loglogistic", maxiter = 100)


ll_curve <- function(t){
  y_hat=mean(fit_ll$linear.predictors)
  scale=fit_ll$scale
  S_i = (1+exp((log(t)-y_hat)/scale))^(-1)
  return(S_i)
}


km<-survfit(Surv(df$lenght_stay,df$state)~1)
plot(surv_ll,type = "l",lwd=2,col=2,xlab="LOS",ylab="Survival probability",ylim=c(0,1))
curve(ln_curve(x), add = TRUE, col = "blue", lwd=2)
lines(km,conf.int=F,col ="green",lwd=2)
legend(x = "topright",          # Position
       legend = c("Bayesian", "Frequentist","KM"),  # Legend texts
       lty = c(1, 1,1),           # Line types
       col = c(2,"blue",3),
       lwd = 2) 

