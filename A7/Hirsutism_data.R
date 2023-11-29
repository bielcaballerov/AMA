hirs <- read.table("hirsutism.dat",header=T, sep="\t",fill=TRUE)

summary(hirs)
attach(hirs)

boxplot(hirs[,2:5])

par(mfrow=c(2,2))
boxplot(hirs[,2]~Treatment,ylim=c(0,30), main=names(hirs)[2], xlab="Treatment")
boxplot(hirs[,3]~Treatment,ylim=c(0,30), main=names(hirs)[3], xlab="Treatment")
boxplot(hirs[,4]~Treatment,ylim=c(0,30), main=names(hirs)[4], xlab="Treatment")
boxplot(hirs[,5]~Treatment,ylim=c(0,30), main=names(hirs)[5], xlab="Treatment")
par(mfrow=c(1,1))

par(mfrow=c(2,2))
boxplot(hirs[Treatment==0,2:5],ylim=c(0,30), main="Treatment 0")
boxplot(hirs[Treatment==1,2:5],ylim=c(0,30), main="Treatment 1")
boxplot(hirs[Treatment==2,2:5],ylim=c(0,30), main="Treatment 2")
boxplot(hirs[Treatment==3,2:5],ylim=c(0,30), main="Treatment 3")
par(mfrow=c(1,1))
