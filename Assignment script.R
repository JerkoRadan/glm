# Assigment 1
# Jerko Radan Cruz
# r0773435

# ----------------- Q1 ----------------------------------

par(mfrow=c(1,2))
x <- seq(1/10,5,0.2)
plot(x, (1/x)-1+log(x), main="KL [from 1/10 to 5 every 0.2]", ylab="KL(fg,fe)", xlab= "lambda / theta", type="b")

x1 <- seq(1/100,100,0.005)
plot(x1, (1/x1)-1+log(x1), main="KL [from 1/100 to 100 every 0.005]", ylab="KL(fg,fe)", xlab= "lambda / theta", type="l")


# ----------------- Q2: Part a ----------------------------------

library(MASS)
library(MuMIn)


fulldata <- read.table("dataHW.txt", sep="", header = T, na.strings="")

# Construct regression models and variable selection methods

# Create a personal subsample
studentnumber = 773435
set.seed(studentnumber)
rownumbers = sample(1:nrow(fulldata),600, replace=F)
mydata= fulldata[rownumbers,]

hist(mydata$Y,main="Histogram for number of crashes") # Poisson (discrete response)
hist(mydata$X2,main="Histogram for X2 (daily traffic)")

glm.full=glm(Y~.,data=mydata, family=poisson())
summary(glm.full)
plot(glm.full)


bm_aic<- stepAIC(glm.full, k=2, direction='both', scope=list(upper=~.^2, lower=~1))
bm_bic<- stepAIC(glm.full, k=log(nrow(mydata)), direction='both', scope=list(upper=~.^2, lower=~1))
bm_hqic<- stepAIC(glm.full, k=log(log(nrow(mydata))), direction='both', scope=list(upper=~.^2, lower=~1))


AIC(bm_aic)
BIC(bm_bic)
# const <- 2+nrow(mydata)+nrow(mydata)*log(2*pi) # Not needed, stepAIC is calculating the correct value for AIC, BIC and HQIC

summary(bm_aic)
bm_bic
bm_hqic

# Calculating for all submodels
glm.full2=glm(Y~.^2,data=mydata, family=poisson())
summary(glm.full2)

options(na.action = "na.fail")
ms.aic <- dredge(glm.full2, rank=AIC)
ms.bic <- dredge(glm.full2, rank=BIC)
#ms.hqic <- dredge(glm.full2, rank=HQIC) #Doesn't work

options(na.action = "na.omit")



#------------------ Q2: Part b -------------------------

#install.packages('fic')
library(fic)

table(mydata$X1)
table(mydata$X2)
table(mydata$X3)
table(mydata$X4) # 559 observations present a value of 12
table(mydata$X5) # 244 observations present a value of 8

hist(mydata$X2)
median(mydata$X2)
mean(mydata$X2)
quantile(mydata$X2)
max(mydata$X2)

names(coef(glm.full2))

glm.full2<-glm(Y~.^2 -X1:X3 -X2:X3 -X3:X4 -X3:X5,data=mydata, family=poisson())
inds0 = c(1,0,1,rep(0, c(length(coef(glm.full2))-3)))
combs <- all_inds(glm.full2, inds0)
combs2 <- with(combs, combs[!(((combs$X1*combs$X2)!=1 & combs$'X1:X2'==1) |
                                ((combs$X1*combs$X4)!=1 & combs$'X1:X4'==1) |
                                ((combs$X1*combs$X5)!=1 & combs$'X1:X5'==1) |
                                ((combs$X2*combs$X4)!=1 & combs$'X2:X4'==1) |
                                ((combs$X2*combs$X5)!=1 & combs$'X2:X5'==1) |
                                ((combs$X4*combs$X5)!=1 & combs$'X4:X5'==1)),])

names(combs)

X<- c(1,median(mydata$X1[mydata$X2==2325]),2325,median(mydata$X3[mydata$X2==2325]),median(mydata$X4[mydata$X2==2325]),median(mydata$X5[mydata$X2==2325]),
      median(mydata$X1[mydata$X2==2325])*2325, median(mydata$X1[mydata$X2==2325])*median(mydata$X4[mydata$X2==2325]),median(mydata$X1[mydata$X2==2325])*median(mydata$X5[mydata$X2==2325]),
      2325*median(mydata$X4[mydata$X2==2325]), 2325*median(mydata$X5[mydata$X2==2325]),
      median(mydata$X4[mydata$X2==2325])*median(mydata$X5[mydata$X2==2325]))

X2<- c(1,median(mydata$X1[mydata$X2==32050]),32050,median(mydata$X3[mydata$X2==32050]),median(mydata$X4[mydata$X2==32050]),median(mydata$X5[mydata$X2==32050]),
       median(mydata$X1[mydata$X2==32050])*32050, median(mydata$X1[mydata$X2==32050])*median(mydata$X4[mydata$X2==32050]),median(mydata$X1[mydata$X2==32050])*median(mydata$X5[mydata$X2==32050]),
       32050*median(mydata$X4[mydata$X2==32050]), 32050*median(mydata$X5[mydata$X2==32050]),
       median(mydata$X4[mydata$X2==32050])*median(mydata$X5[mydata$X2==32050]))

X3<- c(1,median(mydata$X1[mydata$X2==12780]),12780,median(mydata$X3[mydata$X2==12780]),median(mydata$X4[mydata$X2==12780]),median(mydata$X5[mydata$X2==12780]),
       median(mydata$X1[mydata$X2==12780])*12780, median(mydata$X1[mydata$X2==12780])*median(mydata$X4[mydata$X2==12780]),median(mydata$X1[mydata$X2==12780])*median(mydata$X5[mydata$X2==12780]),
       12780*median(mydata$X4[mydata$X2==12780]), 12780*median(mydata$X5[mydata$X2==12780]),
       median(mydata$X4[mydata$X2==12780])*median(mydata$X5[mydata$X2==12780]))

focus2 <- function(par, X) exp(X %*% par)

fic1.exp <- fic(wide=glm.full2, inds=combs2, inds0=inds0, focus=focus2, X=X)
fic2.exp <- fic(wide=glm.full2, inds=combs2, inds0=inds0, focus=focus2, X=X2)
fic3.exp <- fic(wide=glm.full2, inds=combs2, inds0=inds0, focus=focus2, X=X3)

summary(fic1.exp, adj=T)
top3_fic1 <- fic1.exp[order(fic1.exp$rmse.adj),]
top3_fic1

summary(fic2.exp, adj=T)
top3_fic2 <- fic2.exp[order(fic2.exp$rmse.adj),]
top3_fic2

summary(fic3.exp, adj=T)
top3_fic3 <- fic3.exp[order(fic3.exp$rmse.adj),]
top3_fic3

