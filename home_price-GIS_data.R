setwd("C:/Users/ZHOUHANGXU/Desktop/baltimore")
library(spBayes)
library(MBA)
library(fields)

###[1]data process
#Y:PRICE1
#X:NROOM, AGE, LOTSE1, CITCOU
DATA <- read.csv("baltim.csv",header=T)
DATA= data.frame(DATA)
DATA$PRICE1 <- log(DATA$PRICE)
DATA$LOTSZ1 <- log(DATA$LOTSZ)
COORDS <- as.matrix(DATA[,c('X','Y')]) #coordinates


###[2] hierarchical model
##model comparison
n.sample <- 1000
burn.in <- floor(0.75*n.sample)

#model1:ols
ols <- lm(PRICE1~NROOM+AGE+LOTSZ1+CITCOU, data = DATA)
cand.1 <- bayesLMRef(ols, n.sample) #Bayesian Linear Regression

#model2:空間だけ
cand.2 <- spLM(PRICE~1, data=DATA, coords=COORDS, 
               starting=list("phi"=3/200, "sigma.sq"=0.08,"tau.sq"=0.02),
               tuning=list("phi"=0.1,"sigma.sq"=0.05,"tau.sq"=0.05),
               priors=list("phi.Unif"=c(3/1500, 3/50), "sigma.sq.IG"=c(2, 0.08),"tau.sq.IG"=c(2, 0.02)),
               cov.model="exponential",n.samples=n.sample)

#model3: すべて
cand.3 <- spLM(PRICE1~NROOM+AGE+LOTSZ1+CITCOU, data=DATA, coords=COORDS, 
               starting=list("phi"=3/200, "sigma.sq"=0.08,"tau.sq"=0.02),
               tuning=list("phi"=0.1,"sigma.sq"=0.05,"tau.sq"=0.05),
               priors=list("phi.Unif"=c(3/1500, 3/50), "sigma.sq.IG"=c(2, 0.08),"tau.sq.IG"=c(2, 0.02)),
               cov.model="exponential",n.samples=n.sample)

cand.1.DIC <- spDiag(cand.1, start=burn.in, verbose=FALSE)
cand.2 <- spRecover(cand.2, start=burn.in, verbose=FALSE)
cand.2.DIC <- spDiag(cand.2, verbose=FALSE)
cand.3 <- spRecover(cand.3, start=burn.in, verbose=FALSE)
cand.3.DIC <- spDiag(cand.3, start=burn.in, verbose=FALSE)

cand.1.DIC
cand.2.DIC
cand.3.DIC


##model3
#事後分布
beta.samples = cand.3$p.beta.recover.samples #beta
theta.samples = cand.3$p.theta.recover.samples #sigma^2, tau^2 and phi
w.samples = cand.3$p.w.recover.samples #spatial

#サンプルパス
#beta
par(mfrow = (c(2, 3)))
plot(beta.samples[,'(Intercept)'], density=FALSE, main='path.intercept')
plot(beta.samples[,'NROOM'], density=FALSE, main='path.NROOM')
plot(beta.samples[,'AGE'], density=FALSE, main='path.AGE')
plot(beta.samples[,'LOTSZ1'], density=FALSE, main='path.LOTSZ1')
plot(beta.samples[,'CITCOU'], density=FALSE, main='CITCOU')

#sigma^2, tau^2, phi
par(mfrow = (c(1, 3)))
plot(theta.samples[,'sigma.sq'], density=FALSE, main='path.sigma.sq')
plot(theta.samples[,'tau.sq'], density=FALSE, main='path.tau.sq')
plot(theta.samples[,'phi'], density=FALSE, main='path.phi')


#事後密度
#beta
dens.i = density(beta.samples[,'(Intercept)'])
dens.n = density(beta.samples[,'NROOM'])
dens.a = density(beta.samples[,'AGE'])
dens.l = density(beta.samples[,'LOTSZ1'])
dens.c = density(beta.samples[,'CITCOU'])

par(mfrow = (c(2, 3)))
plot(dens.i, main='density.intercept')
plot(dens.n, main='density.NROOM')
plot(dens.a, main='density.AGE')
plot(dens.l, main='density.LOTSZ1')
plot(dens.c, main='CITCOU')

##sigma^2, tau^2, phi
dens.s = density(theta.samples[,'sigma.sq'])
dens.t = density(theta.samples[,'tau.sq'])
dens.p = density(theta.samples[,'phi'])

par(mfrow = (c(1,3)))
plot(dens.s, main='density.sigma.sq')
plot(dens.t, main='density.tau.sq')
plot(dens.p, main='density.phi')


#自相関関数
#beta
par(mfrow = (c(2, 3)))
ts.i = ts(beta.samples[,'(Intercept)'])
acf.i <- acf(ts.i)
ts.n = ts(beta.samples[,'NROOM'])
acf.n <- acf(ts.n)
ts.a = ts(beta.samples[,'AGE'])
acf.a <- acf(ts.a)
ts.l = (beta.samples[,'LOTSZ1'])
acf.l <- acf(ts.l)
ts.c = (beta.samples[,'CITCOU'])
acf.c <- acf(ts.c)

par(mfrow = (c(2, 3)))
plot(acf.i, main='acf.intercept')
plot(acf.n, main='acf.NROOM')
plot(acf.a, main='acf.AGE')
plot(acf.l, main='acf.LOTSZ1')
plot(acf.c, main='acf.CITCOU')

###sigma^2, tau^2, phi
par(mfrow = (c(1, 3)))
ts.s = ts(theta.samples[,'sigma.sq'])
acf.s <- acf(ts.s)
ts.t = ts(theta.samples[,'tau.sq'])
acf.t <- acf(ts.t)
ts.p = ts(theta.samples[,'phi'])
acf.p <- acf(ts.p)

par(mfrow = (c(1, 3)))
plot(acf.s, main='acf.sigma.sq')
plot(acf.t, main='acf.tau.sq')
plot(acf.c, main='acf.phi')

###[3]
x.res=100 #格点の数
y.res=100
par(mfrow = (c(1, 3)))

##y
surf <- mba.surf(cbind(COORDS, DATA$PRICE1), no.X=x.res, no.Y=y.res, extend=FALSE)$xyz.est #二维表面
z.lim <- range(surf[[3]], na.rm=TRUE) #表面数据的取值范围
image.plot(surf, xaxs="r", yaxs="r", main="被説明変数")

##平均空間効果
#空間効果の事後平均と標準偏差
w.hat.mu <- apply(w.samples,1,mean)
w.hat.sd <- apply(w.samples,1,sd)

surf <- mba.surf(cbind(COORDS, w.hat.mu), no.X=x.res, no.Y=y.res, extend=FALSE)$xyz.est
image.plot(surf, xaxs="r", yaxs="r", main="空間成分の事後平均")

##残差
r1 <- DATA$PRICE1 - w.hat.mu
surf <- mba.surf(cbind(COORDS, r1), no.X=x.res, no.Y=y.res, extend=FALSE)$xyz.est
image.plot(surf, xaxs="r", yaxs="r", main="残差")

