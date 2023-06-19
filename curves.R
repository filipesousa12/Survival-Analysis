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

df<- read_excel('data_surv.xlsx')
fit<- survreg(formula = Surv(lenght_stay, state) ~ DM + Hemoglobin + 
                Neutrofils + fibrinogen + Ddimer + LRA + INR + pH + pCO2 + 
                FiO2 + FC , data = df, dist = "lognormal")

ln_curve <- function(t){
  y_hat=mean(fit$linear.predictors)
  scale=fit$scale
  S_i = 1-pnorm((log(t)-y_hat)/scale)
  return(S_i)
}

var_cont<-c("Hemoglobin","Neutrofils","fibrinogen","Ddimer","INR","pH","pCO2","FiO2","FC")
var_cat<-c("DM","LRA")


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

var_cont_l<-c("Hemoglobin","Neutrofils","fibrinogen","Ddimer","INR","pH","pCO2","HCO3","FiO2","FC","Troponine")

fit_ll$coefficients


ll_curve_factor <- function(t){
  x_cont=colMeans(df[,var_cont_l])
  x_cat=c(0,0)
  b0<-fit_ll$coefficients[1]
  b_cat<-fit_ll$coefficients[c(2,7)]
  b_cont<-fit_ll$coefficients[-c(1,2,7)]
  y_hat=as.numeric(b0+b_cat%*%x_cat+b_cont%*%x_cont)
  scale=fit$scale
  S_i = (1+exp((log(t)-y_hat)/scale))^(-1)
  return(S_i)
}

ll_curve_factor_1 <- function(t){
  x_cont=colMeans(df[,var_cont_l])
  x_cat=c(1,0)
  b0<-fit_ll$coefficients[1]
  b_cat<-fit_ll$coefficients[c(2,7)]
  b_cont<-fit_ll$coefficients[-c(1,2,7)]
  y_hat=as.numeric(b0+b_cat%*%x_cat+b_cont%*%x_cont)
  scale=fit$scale
  S_i = (1+exp((log(t)-y_hat)/scale))^(-1)
  return(S_i)
}

ll_curve_factor_2 <- function(t){
  x_cont=colMeans(df[,var_cont_l])
  x_cat=c(0,1)
  b0<-fit_ll$coefficients[1]
  b_cat<-fit_ll$coefficients[c(2,7)]
  b_cont<-fit_ll$coefficients[-c(1,2,7)]
  y_hat=as.numeric(b0+b_cat%*%x_cat+b_cont%*%x_cont)
  scale=fit$scale
  S_i = (1+exp((log(t)-y_hat)/scale))^(-1)
  return(S_i)
}

ll_curve_factor_3 <- function(t){
  x_cont=colMeans(df[,var_cont_l])
  x_cat=c(1,1)
  b0<-fit_ll$coefficients[1]
  b_cat<-fit_ll$coefficients[c(2,7)]
  b_cont<-fit_ll$coefficients[-c(1,2,7)]
  y_hat=as.numeric(b0+b_cat%*%x_cat+b_cont%*%x_cont)
  scale=fit$scale
  S_i = (1+exp((log(t)-y_hat)/scale))^(-1)
  return(S_i)
}

ln_curve_factor <- function(t){
  x_cont=colMeans(df[,var_cont])
  x_cat=c(0,0)
  b0<-fit$coefficients[1]
  b_cat<-fit$coefficients[c(2,7)]
  b_cont<-fit$coefficients[-c(1,2,7)]
  y_hat=as.numeric(b0+b_cat%*%x_cat+b_cont%*%x_cont)
  scale=fit$scale
  S_i = 1-pnorm((log(t)-y_hat)/scale)
  return(S_i)
}


ln_curve_factor_1 <- function(t){
  x_cont=colMeans(df[,var_cont])
  x_cat=c(1,0)
  b0<-fit$coefficients[1]
  b_cat<-fit$coefficients[c(2,7)]
  b_cont<-fit$coefficients[-c(1,2,7)]
  y_hat=as.numeric(b0+b_cat%*%x_cat+b_cont%*%x_cont)
  scale=fit$scale
  S_i = 1-pnorm((log(t)-y_hat)/scale)
  return(S_i)
}

ln_curve_factor_2 <- function(t){
  x_cont=colMeans(df[,var_cont])
  x_cat=c(0,1)
  b0<-fit$coefficients[1]
  b_cat<-fit$coefficients[c(2,7)]
  b_cont<-fit$coefficients[-c(1,2,7)]
  y_hat=as.numeric(b0+b_cat%*%x_cat+b_cont%*%x_cont)
  scale=fit$scale
  S_i = 1-pnorm((log(t)-y_hat)/scale)
  return(S_i)
}

ln_curve_factor_3 <- function(t){
  x_cont=colMeans(df[,var_cont])
  x_cat=c(1,1)
  b0<-fit$coefficients[1]
  b_cat<-fit$coefficients[c(2,7)]
  b_cont<-fit$coefficients[-c(1,2,7)]
  y_hat=as.numeric(b0+b_cat%*%x_cat+b_cont%*%x_cont)
  scale=fit$scale
  S_i = 1-pnorm((log(t)-y_hat)/scale)
  return(S_i)
}


curve(ll_curve_factor_3(x),from=0,to=60,col="black", lwd=2,xlab = "LOS",ylab="Survival probability")
curve(ln_curve_factor_3(x),from=0,to=60,col="blue", lwd=2,add=T)
legend(x = "topright",          # Position
       legend = c("Log-Logistic","Log-Normal"),  # Legend texts
       lty = c(1, 1),           # Line types
       col = c("black","blue"),
       lwd = 2) 
