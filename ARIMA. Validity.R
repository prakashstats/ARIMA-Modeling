rm(list=ls())
#Importing two year 2015-16 SBI share price data from Excel to R
setwd("C:/Users/PRAKASH CHANDRA/Desktop/Project")
data1<-read.csv(file.choose())
attach(data1);View(data1)
plot.ts(Close,ylab="Stock_price",col="Red",main="SBIN data plotting Jan 2015 to Dec 2016")

# // Testing the data stationary 
library(tseries)                      # Apply Dicky Fuler test for stationarity
adf.test(Close,alternative="stationary")
y=diff(Close)       
adf.test(y,alternative="stationary")
plot(y,ylab="Stock_price",xlab="Time",col="blue",type="l",main="Plotting of Stock_price with 1st difference")
# Now we apply ARIMA model for forecasting

# // Finding the value of Auto correlation
acf1<-function(p){
n<<-length(y)
m<<-mean(y)
y1<<-c(0);y2<<-c(0);j=1;l=1
A<<-c(0);B<<-c(0)
for(i in (p+1):n){
 A[j]<<-y[i]-m
 B[j]<<-y[i-p]-m
 j=j+1
}
for(k in 1:n){
 y1[l]=(y[k]-m)^2
 l=l+1
}
z1<<-sqrt((sum(y1))/n)
y2=(B^2)
z2<<-sqrt((sum(y2))/(n-p))
C=sum(A*B)
num<<-C/(n-p)
den<<-z1*z2
acf2<<-num/den
return(acf2)
}
acf(y,plot=FALSE)


# //Finding the value of Partial Autocorrelation
pacf1<-function(p){
q4=matrix(c(acf1(0),acf1(1),acf1(2),acf1(3),acf1(1),acf1(0),acf1(1),acf1(2),acf1(2),acf1(1),acf1(0),acf1(1)),4,4,T)
ac4<-matrix(c(acf1(1),acf1(2),acf1(3),acf1(4)),4,1,T)
q1=q4[1,1]
q2=matrix(0,2,2)
q3=matrix(0,3,3)
ac1=ac4[1,1]
ac2=matrix(c(0),2,1,T)
ac3=matrix(c(0),3,1,T)
for(i in 1:2){
for(j in 1:2){
q2[i,j]=q4[i,j]
ac2[i,1]=ac4[i,1]
}}
for(i in 1:3){
for(j in 1:3){
q3[i,j]=q4[i,j]
ac3[i,1]=ac4[i,1]
}}
if(p==1){((q=q1) && (ac=ac1))}
else if (p==2){((q=q2) && (ac=ac2))}
else if (p==3){((q=q3) && (ac=ac3))}
else if (p==4){((q=q4) && (ac=ac4))}
phi=solve(q)%*%ac       # Finding the Cofficient of AR model  
return(phi)             # we use yule's Walker Equation
}
acf(y)
pacf(y)   # Use two inbuilt command in R for finding Order of P & Q in ARIMA model 
pacf1(1)  # From the Correlogram it is clear that the value of P & Q is 1
acf1(1)


# finding the cofficient of MA model

f<-function(a,b,c){
X=(-b+sqrt(b^2-4*a*c))/(2*a)
Y=(-b-sqrt(b^2-4*a*c))/(2*a)
z=c(X,Y)
return(z)
}
MA<-function(q){
 a1=f(acf1(1),1,acf1(1))
 q1=a1[which(a1>=-1 && a1<=1)]
 a2=f(acf1(2),1,(acf1(2)+(acf1(2)*q1^2)-q1^2))
  q2=a2[which(a2>=-1 && a2<=1)]
 a3<-f(acf1(3),1,((acf1(3)*(1+q1^2+q2^2))-2*q1*q2))
 q3<-a3[which(a3>=-1 && a3<=1)]
  if(q==1){ ma=q1}
  else if (q==2){ma=c(q1,q2)}
  else if (q==3){ma=c(q1,q2,q3)}
 return(ma)
}
# Now we apply ARIMA model
# take maximum value of p & q is 3.
Pacf1=as.numeric(pacf1(1))
Pacf2=as.numeric(pacf1(2))
Pacf3=as.numeric(pacf1(3))
ma=MA(3)
#Importing 5 month data of 2017 for validation of the modal
data<-read.csv(file.choose())
attach(data)
close1=c(Close[length(Close)],close)         # small close is the data of 2017
actual<-diff(close1)                         # Capital Close is the data of 2015-16
ARIMA<-function(p,q){
# ARMA(0,1)
P01=c(0)             # In every modal we take 1st forcasting error is 0.
P01[1]=ma[1]*0 
for(t in 2:(length(actual))){
P01[t]=(ma[1]*(actual[t-1]-P01[t-1])) 
}
# ARMA(0,2)
P02=c(0)        
P02[1]=ma[1]*0 + ma[2]*0
P02[2]=ma[1]*(actual[1]-P02[1])+(ma[2]*0)
for(t in 3:(length(actual))){
P02[t]=ma[1]*(actual[t-1]-P02[t-1])+(ma[2]*(actual[t-2]-P02[t-2])) 
}
# ARMA(0,3)
P03=c(0)        
P03[1]=(ma[1]*0)+(ma[2]*0)+(ma[3]*0)
P03[2]=(ma[1]*(actual[1]-P03[1]))+(ma[2]*0)+(ma[3]*0)
P03[3]=(ma[1]*(actual[1]-P03[1]))+(ma[2]*actual[2]-P03[2])+(ma[3]*0)
for(t in 4:(length(actual))){
P03[t]=(ma[1]*(actual[t-1]-P03[t-1]))+(ma[2]*(actual[t-2]-P03[t-2]))+(ma[3]*(actual[t-3]-P03[t-3])) 
}
# ARMA(1,0)
P10=c(0)         
P10[1]=(Pacf1[1]*y[length(y)])  
for(t in 2:(length(actual))){
P10[t]=(Pacf1[1]*P10[t-1])
}
# ARMA(1,1)
P11=c(0)         
P11[1]=(Pacf1[1]*y[length(y)])+ (ma[1]*0)  
for(t in 2:(length(actual))){
P11[t]=(Pacf1[1]*P11[t-1])+ (ma[1]*(actual[t-1]-P11[t-1])) 
}
# ARMA(1,2)
P12=c(0)        
P12[1]=(Pacf1[1]*y[length(y)]) + (ma[1]*0)+(ma[2]*0)
P12[2]=(Pacf1[1]*P12[1])+ (ma[1]*(actual[1]-P12[1]))+(ma[2]*0)
for(t in 3:(length(actual))){
P12[t]=(Pacf1[1]*P12[t-1])+(ma[1]*(actual[t-1]-P12[t-1]))+(ma[2]*(actual[t-2]-P12[t-2])) 
}
# ARMA(1,3)
P13=c(0)        
P13[1]=(Pacf1[1]*y[length(y)])+ (ma[1]*0)+(ma[2]*0)+(ma[3]*0)
P13[2]=(Pacf1[1]*P13[1])+ (ma[1]*(actual[1]-P13[1]))+(ma[2]*0)+(ma[3]*0)
P13[3]=(Pacf1[1]*P13[2]) + (ma[1]*(actual[1]-P13[1]))+(ma[2]*(actual[2]-P13[2]))+(ma[3]*0)
for(t in 4:(length(actual))){
P13[t]=(Pacf1[1]*P13[t-1])+(ma[1]*(actual[t-1]-P13[t-1]))+(ma[2]*(actual[t-2]-P13[t-2]))+(ma[2]*(actual[t-3]-P13[t-3])) 
}
# ARMA(2,0)
P20=c(0)        
P20[1]=(Pacf2[1]*y[length(y)])+ (Pacf2[2]*y[length(y)-1]) 
P20[2]=(Pacf2[1]*P20[1])+ (Pacf2[2]*y[length(y)]) 
for(t in 3:(length(actual))){
P20[t]=(Pacf2[1]*P20[t-1])+ (Pacf2[2]*P20[t-2]) 
}
# ARMA(2,1)
P21=c(0)        
P21[1]=(Pacf2[1]*y[length(y)])+ (Pacf2[2]*y[length(y)-1]) + (ma[1]*0)
P21[2]=(Pacf2[1]*P21[1])+ (Pacf2[2]*y[length(y)]) + (ma[1]*(actual[1]-P21[1]))
for(t in 3:(length(actual))){
P21[t]=(Pacf2[1]*P21[t-1])+ (Pacf2[2]*P21[t-2])+(ma[1]*(actual[t-1]-P21[t-1])) 
}
# ARMA(2,2)
P22=c(0)        
P22[1]=(Pacf2[1]*y[length(y)])+ (Pacf2[2]*y[length(y)-1]) + (ma[1]*0)+(ma[2]*0)
P22[2]=(Pacf2[1]*P22[1])+ (Pacf2[2]*y[length(y)]) + (ma[1]*(actual[1]-P22[1]))+(ma[2]*0)
for(t in 3:(length(actual))){
P22[t]=(Pacf2[1]*P22[t-1])+ (Pacf2[2]*P22[t-2])+(ma[1]*(actual[t-1]-P22[t-1]))+(ma[2]*(actual[t-2]-P22[t-2])) 
}
# ARMA(2,3)
P23=c(0)        
P23[1]=(Pacf2[1]*y[length(y)])+ (Pacf2[2]*y[length(y)-1]) + (ma[1]*0)+(ma[2]*0)+(ma[3]*0)
P23[2]=(Pacf2[1]*P23[1])+ (Pacf2[2]*y[length(y)]) + (ma[1]*(actual[1]-P23[1]))+(ma[2]*0)+(ma[3]*0)
P23[3]=(Pacf2[1]*P23[2])+ (Pacf2[2]*P23[1]) + (ma[1]*(actual[1]-P23[1]))+(ma[2]*(actual[2]-P23[2]))+(ma[3]*0)
for(t in 4:(length(actual))){
P23[t]=(Pacf2[1]*P23[t-1])+ (Pacf2[2]*P23[t-2])+(ma[1]*(actual[t-1]-P23[t-1]))+(ma[2]*(actual[t-2]-P23[t-2]))+(ma[2]*(actual[t-3]-P23[t-3])) 
}
# ARMA(3,0)
P30=c(0)        
P30[1]=(Pacf3[1]*y[length(y)])+ (Pacf3[2]*y[length(y)-1])+(Pacf3[3]*y[length(y)-2]) 
P30[2]=(Pacf3[1]*P30[1])+ (Pacf3[2]*y[length(y)])+(Pacf3[3]*y[length(y)-1]) 
P30[3]=(Pacf3[1]*P30[2])+ (Pacf3[2]*P30[1])+(Pacf3[3]*y[length(y)]) 
for(t in 4:(length(actual))){
P30[t]=(Pacf3[1]*P30[t-1])+ (Pacf3[2]*P30[t-2])+(Pacf3[3]*P30[t-3])
}
# ARMA(3,1)
P31=c(0)        
P31[1]=(Pacf3[1]*y[length(y)])+ (Pacf3[2]*y[length(y)-1])+(Pacf3[3]*y[length(y)-2]) + (ma[1]*0)
P31[2]=(Pacf3[1]*P31[1])+ (Pacf3[2]*y[length(y)])+(Pacf3[3]*y[length(y)-1]) + (ma[1]*(actual[1]-P31[1]))
P31[3]=(Pacf3[1]*P31[2])+ (Pacf3[2]*P31[1])+(Pacf3[3]*y[length(y)]) + (ma[1]*(actual[1]-P31[1]))
for(t in 4:(length(actual))){
P31[t]=(Pacf3[1]*P31[t-1])+ (Pacf3[2]*P31[t-2])+(Pacf3[3]*P31[t-3])+(ma[1]*(actual[t-1]-P31[t-1])) 
}
# ARMA(3,2)
P32=c(0)        
P32[1]=(Pacf3[1]*y[length(y)])+ (Pacf3[2]*y[length(y)-1])+(Pacf3[3]*y[length(y)-2]) + (ma[1]*0)+(ma[2]*0)
P32[2]=(Pacf3[1]*P32[1])+ (Pacf3[2]*y[length(y)])+(Pacf3[3]*y[length(y)-1]) + (ma[1]*(actual[1]-P32[1]))+(ma[2]*0)
P32[3]=(Pacf3[1]*P32[2])+ (Pacf3[2]*P32[1])+(Pacf3[3]*y[length(y)]) + (ma[1]*(actual[1]-P32[1]))+(ma[2]*(actual[2]-P32[2]))
for(t in 4:(length(actual))){
P32[t]=(Pacf3[1]*P32[t-1])+ (Pacf3[2]*P32[t-2])+(Pacf3[3]*P32[t-3])+(ma[1]*(actual[t-1]-P32[t-1]))+(ma[2]*(actual[t-2]-P32[t-2])) 
}
# ARMA(3,3)
P33=c(0)        
P33[1]=(Pacf3[1]*y[length(y)])+ (Pacf3[2]*y[length(y)-1])+(Pacf3[3]*y[length(y)-2]) + (ma[1]*0)+(ma[2]*0)+(ma[3]*0)
P33[2]=(Pacf3[1]*P33[1])+ (Pacf3[2]*y[length(y)])+(Pacf3[3]*y[length(y)-1]) + (ma[1]*(actual[1]-P33[1]))+(ma[2]*0)+(ma[3]*0)
P33[3]=(Pacf3[1]*P33[2])+ (Pacf3[2]*P33[1])+(Pacf3[3]*y[length(y)]) + (ma[1]*(actual[1]-P33[1]))+(ma[2]*(actual[2]-P33[2]))+(ma[3]*0)
for(t in 4:(length(actual))){
P33[t]=(Pacf3[1]*P33[t-1])+ (Pacf3[2]*P33[t-2])+(Pacf3[3]*P33[t-3])+(ma[1]*(actual[t-1]-P33[t-1]))+(ma[2]*(actual[t-2]-P33[t-2]))+(ma[3]*(actual[t-3]-P33[t-3])) 
}
if((p==0) && (q==1)) { ARMA=P01}
else if ((p==0) && (q==2)){ ARMA=P02}
else if ((p==0) && (q==3)){ ARMA=P03}
else if ((p==1) && (q==0)){ ARMA=P10}
else if ((p==1) && (q==1)){ ARMA=P11}
else if ((p==1) && (q==2)){ ARMA=P12}
else if ((p==1) && (q==3)){ ARMA=P13}
else if ((p==2) && (q==0)){ ARMA=P20}
else if ((p==2) && (q==1)){ ARMA=P21}
else if ((p==2) && (q==2)){ ARMA=P22}
else if ((p==2) && (q==3)){ ARMA=P23}
else if ((p==3) && (q==0)){ ARMA=P30}
else if ((p==3) && (q==1)){ ARMA=P31}
else if ((p==3) && (q==2)){ ARMA=P32}
else if ((p==3) && (q==3)){ ARMA=P33}
return(ARMA)
}
predict=ARIMA(1,1)     
Predict=c(0)         # Predict is the predicted value w.r.to close value
Predict[1]=Close[length(Close)]+Predict[1]
for(i in 2:length(predict)){
Predict[i]=close[i-1]+predict[i]
}
Predict

error=close-Predict
error_percent=(abs(error)*100)/close
compare=cbind(close,Predict,error,error_percent)
compare
Actual=close
plot(Actual,type="l",xlab="Time",ylab="Stock_Price",col="red",lty=1,main="Model Valadity")
lines(Predict,type="l",col="blue",lty=2)
legend("topleft", legend=c("Actual","Predict"),lty=1:2,col=c("red","blue"))








