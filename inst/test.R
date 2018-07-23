#library(devtools)
#install_github("mcfanda/cpower")
library(cpower)
### examples
##### d.cont ###
cont<-c(-3,-1,1,3)
means<-c(10,12,10,12)
d.contr(cont,means = means,sd=2,scale = "g")
### different scaling
d.contr(cont,means = means,sd=2,scale = "z")

### raw data
y<-rep(means,1000)+rnorm(4000,0,2)
x<-factor(rep(1:4,1000))
d.contr(cont,y = y,x = x)

###### ci.contr ######
cont<-c(-3,-1,1,3)
means<-c(10,12,10,12)
d<-d.contr(cont,means = means,sd=2,scale = "g")
ci.contr(cont,d=d,scale = "g",n=100)

###### sims

cont<-c(-3,-1,1,3)
means<-c(10,12,10,12)
one.sample(cont,means,sd=1,n=10)

## check actual power with 10 simulations
pwr<-power.contr.simulate(10,cont,10,means,sd=1)
mean(pwr$`Pr(>|t|)`>=.05)
summary(pwr)

#### eta2.contr
### sim some data
cont<-c(-3,-1,1,3)
means<-c(10,12,10,12)
y<-rep(means,1000)+rnorm(4000,0,2)
x<-factor(rep(1:4,1000))
### compute eta-squared
mod<-lm(y~x)
eta2.contr(cont,xname = "x",model = mod)

#### check from the d coefficient

d<-d.contr(cont,x=x,y=y)
eta2.contr.d(cont = cont,d=d,scale = "g")



#### contr.custom

### sim some data
cont<-c(-3,-1,1,3)
means<-c(10,12,10,12)
y<-rep(means,1000)+rnorm(4000,0,2)
x<-factor(rep(1:4,1000))
### assign contrast weights
contrasts(x)<-contr.custom(cont)
contrasts(x)
cor(contr.custom(cont))
summary(lm(y~x))


### test.contr

cont<-c(-3,-1,1,3)
means<-c(10,12,10,12)
y<-rep(means,1000)+rnorm(4000,0,10)
x<-rep(1:4,1000)
dat<-as.data.frame(cbind(y,x))
dat$x<-factor(dat$x)

test.contr(data = dat,cont=cont,yname = "y",xname = "x")

# check the contrast value
observed<-tapply(dat$y,dat$x,mean)
observed%*%cont

# check the t-test and p-value
contrasts(dat$x)<-contr.custom(cont)
summary(lm(y~x,data=dat))



### some try ###########
library(cpower)
power.t.test(n=10,delta=1)
power.contrast.t(cont=c(1,-1),d=1,n=10)

power.t.test(power=.80,delta=1)
power.contrast.t(cont=c(1,-1),d=1,power=.80)


####### check if theoretical power works with simulated data #########
cont<-c(1,-1/2,-1/2)
m<-c(30,25,21)
sp<-15
n=12
res<-power.contr.simulate(5000,n=n,cont = cont,means = m,sd=sp)
summary.simulations(res)
(dg<-d.contr(cont,m,sp))
power.contrast.t(cont=cont,d=dg,n=n)

####### check if expected N works with simulated data #########

(pp<-power.contrast.t(cont=cont,d=dg,power=.70))
n<-round(pp$n)
res<-power.contr.simulate(5000,n=n,cont = cont,means = m,sd=sp)
summary.simulations(res)



