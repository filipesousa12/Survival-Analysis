
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
library(survminer)

df<- read_excel('data_surv.xlsx')





obj.cox<-coxph(formula = Surv(lenght_stay, state) ~
             Hypertension + DM + Hemoglobin + Neutrofils + fibrinogen + 
             Ddimer + LRA + TP + INR + FA + pH + pCO2 + HCO3 + FiO2 + 
             FC + Troponine, data = df, iter.max = 50)
summary(obj.cox)
AIC(obj.cox)


cox_fit <- survfit(obj.cox)
km<-survfit(Surv(df$lenght_stay,df$state)~1)



fit<- survreg(formula = Surv(lenght_stay, state) ~ DM + Hemoglobin + 
                Neutrofils + fibrinogen + Ddimer + LRA + INR + pH + pCO2 + 
                FiO2 + FC , data = df, dist = "lognormal")


summary(fit)
AIC(fit)


ln_curve <- function(t){
  y_hat=mean(fit$linear.predictors)
  scale=fit$scale
  S_i = 1-pnorm((log(t)-y_hat)/scale)
  return(S_i)
}








my_resid<- (log(df$lenght_stay)-fit$linear.predictors)/fit$scale
my_resid<- data.frame(lp=fit$linear.predictors,my_resid=my_resid)
my_resid
plot(my_resid$lp~my_resid$my_resid)



mean(my_resid$my_resid)
sd<-sd(my_resid$my_resid)
shapiro.test(my_resid$my_resid)
pv1<-round(ks.test(my_resid$my_resid,pnorm,0,sd(my_resid$my_resid))$p.value,2)
pv1

ggplot(data=my_resid, aes(x=my_resid)) +
  geom_histogram(aes(y=..density..),breaks=seq(-3,3,by=0.5),
                 colour="Black", fill="#2171B5")+
  xlab("Residuals") +
  ylab("Density") +
  theme_minimal() +
  geom_function(aes(colour=paste("Normal","- p = ",pv1)), fun = dnorm, 
                args = list(0, sd), size=1.2)+
  scale_color_manual(values=c("red2"))+
  ylim(0,0.6)+
  xlim(-3,3) +
  labs(colour="Distribution")

ggplot(data=my_resid, aes(x=lp,y=my_resid)) +
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





fit_ll<-survreg(formula = Surv(lenght_stay, state) ~
                  DM +  Hemoglobin + Neutrofils + fibrinogen + Ddimer + 
                  LRA + INR + pH + pCO2 + HCO3 + FiO2 + FC + Troponine, data = df, 
                dist = "loglogistic", maxiter = 100)

summary(fit_ll)
AIC(fit_ll)

ll_curve <- function(t){
  y_hat=mean(fit_ll$linear.predictors)
  scale=fit_ll$scale
  S_i = (1+exp((log(t)-y_hat)/scale))^(-1)
  return(S_i)
}







plot(cox_fit,conf.int=F,col =3, lwd=2)
curve(ll_curve(x),add=TRUE,col="red", lwd=2)
curve(ln_curve(x), add = TRUE, col = "blue", lwd=2)
legend(x = "topright",          # Position
       legend = c("Cox","Log-Logistic","Log-Normal"),  # Legend texts
       lty = c(1,1,1),           # Line types
       col = c(3,"red","blue"),
       lwd = 2) 


res<-(log(df$lenght_stay)-fit_ll$linear.predictor)/fit_ll$scale
res<- data.frame(lp=fit_ll$linear.predictor,res=res)
l<-exp(res)
pllogis=fitdist(res$res, "logis")
pllogis
pv2<-round(ks.test (res$res, "plogis",location = pllogis$estimate[1],scale = pllogis$estimate[2])$p.value,2)

pv2
hist(res$res,freq=F)

ggplot(data=res, aes(x=res)) +
  geom_histogram(aes(y=..density..), breaks=seq(-4,4,by=0.7), 
                 colour="Black", fill="#2171B5") + 
  ylim(0,0.75)+
  xlim(-4,4) +
  xlab("residuals") +
  ylab("Density") +
  theme_minimal() +
  geom_function(aes(colour=paste("Logistic","- p = ",pv2)), fun = dlogis, 
                args = list(pllogis$estimate[1], pllogis$estimate[2]), size=1.2)+
  scale_color_manual(values=c("red2"))+
  labs(colour="Distribution")

ggplot(data=res, aes(x=lp,y=res)) +
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

