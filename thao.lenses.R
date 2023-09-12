setwd("~/Documents/Thao")
require(plyr)
require(ggplot2)
require(data.table)
library(reshape2)

mal.dat <- read.csv("Photophores_Malacosteus_niger.csv")
mal.dat <- subset(mal.dat,foam!=is.na(foam))
ech.dat	 <- read.csv("Echiostoma_barbatum.csv")
ech.dat <- subset(ech.dat,foam!=is.na(foam))
hist(ech.dat$sl,col=factor(ech.dat$sex))

mal.dat[mal.dat$sex=="UF","sex"]<- factor("F")

mal.dat$sex <- gsub("UF","f",mal.dat$sex)
mal.dat$sex <- gsub("UM","M",mal.dat$sex)
mal.dat$sex <- gsub("DM","M",mal.dat$sex)
mal.dat$sex <- gsub("DF","f",mal.dat$sex)
mal.dat$sex <- gsub("M","m",mal.dat$sex)
mal.dat$MCZ.Number <- as.character(mal.dat$MCZ.Number)


po.mass <- function(dat){p <- ggplot(aes(x=log(mass),y=log(po)), data=dat)+geom_point(aes(color=sex));print(p)}
lens.mass <- function(dat){p <- ggplot(aes(x=log(mass),y=log(vol)), data=dat)+geom_point(aes(color=sex))+geom_smooth(method="lm",se=F,aes(color=as.factor(sex)));print(p)}
eye.mass <- function(dat){p <- ggplot(aes(x=log(mass),y=log(eye)), data=dat)+geom_point(aes(color=sex));print(p)}

pdf("mal.prelim.pdf")
qplot(data=mal.dat,x=mass,y=log(vol),color=factor(sex))+geom_smooth(method="lm",se=F,aes(colour=factor(sex)))
qplot(data=mal.dat,x=mass,y=log(eye),color=factor(sex))+geom_smooth(method="lm",se=F,aes(colour=factor(sex)))
dev.off()
eye.mass.mal <- eye.mass(mal.dat)
lens.mass.mal <- lens.mass(mal.dat)+geom_smooth(method="lm",se=F,aes(colour=sex))
po.mass.mal <- po.jl(mal.dat)+geom_text(aes(label=as.factor(MCZ.Number)))



prelim.dat <- read.csv("sboa.vol.foam1and2.csv")
prelim.dat <- data.frame(prelim.dat,log.unwrap=log(prelim.dat$unwrap.vol/10^9),log.wrap=log(prelim.dat$wrap.vol/10^9),log.sl=log(prelim.dat$SL),unwrap.r=prelim.dat$unwrap.vol/10^9/(4/3*pi))

qplot(sl,data=foam.dat,geom="histogram",binwidth=10)+facet_grid(sex~.)

qplot(data=prelim.dat,x=SL,y=log(unwrap.r),color=factor(Sex))+geom_smooth(method="lm",se=F,aes(colour=factor(Sex)))

+geom_text(aes(label=lense))

wrap.lm.m <- lm(log.unwrap~log.sl,data=subset(prelim.dat,Sex=="m"))
wrap.lm.f <- lm(log.unwrap~log.sl,data=subset(prelim.dat,Sex=="f"))
coefficients(wrap.lm.f)
coefficients(wrap.lm.m)

summary(aov(log.unwrap~log.sl*Sex,data=prelim.dat))
