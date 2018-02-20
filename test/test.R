### some try ###########
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

