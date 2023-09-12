library(tidyverse)
library(vroom)
library(lamW)





se <- function(x,na.rm=T){if(na.rm==T)
{x <- x[is.na(x)==F]}else{x <- x}
  sd(x)/sqrt(length(x))}

phot.dat <- read.csv("Photostomias_guernei.csv")
qplot(log(mass),log(po),data=phot.dat,col=sex)+geom_label(aes(label=MCZ.Number),size=2)
qplot(log(jl),log(mass),data=phot.dat,col=sex)+geom_label(aes(label=MCZ.Number),size=2)
qplot(sl,log(volume),data=phot.dat,col=sex)+geom_label(aes(label=MCZ.Number),size=2)

phot.mass.miss <- phot.dat[is.na(phot.dat$mass) & is.na(phot.dat$jl)==F,]
phot.mass.miss2 <- phot.dat[is.na(phot.dat$sl) & is.na(phot.dat$jl)==F,]
phot.dat[is.na(phot.dat$mass) & is.na(phot.dat$jl)==F,"mass"] <- exp(predict(lm(log(mass)~log(jl),phot.dat),newdata=data.frame(jl=phot.mass.miss$jl)))
phot.dat[is.na(phot.dat$sl) & is.na(phot.dat$jl)==F,"sl"] <- predict(lm(sl~jl,phot.dat),newdata=data.frame(jl=phot.mass.miss2$jl))

qplot(log(sl),log(mass),data=phot.dat,col=sex)+geom_point()
qplot(log(jl),log(mass),data=phot.dat,col=sex)+geom_label(aes(label=MCZ.Number),size=2)
qplot(log(mass,base=10),log(volume,base = 10),data=phot.dat,col=sex)

summary(aov(log(volume,base=10)~log(mass,base=10)+sex,phot.dat))

ech.dat	 <- read.csv("Echiostoma_barbatum2.csv")

ech.dat %>% 
  filter(sex %in% c("M","F")) %>% 
  ggplot(aes(log(sl),log(volume),col=sex))+geom_point()
qplot(log(sl,base=10),log(po,base=10),data=ech.dat,col=sex)+geom_label(aes(label=MCZ),size=2)
qplot(log(sl),log(mass),data=ech.dat,col=sex)+geom_label(aes(label=MCZ),size=2)
qplot(log(sl),log(volume),data=ech.dat,col=sex)+geom_label(aes(label=MCZ),size=2)

summary(aov(log(volume)~log(sl)+sex,ech.dat))


mal.dat.thao	 <- read.csv("Malacosteus_niger.csv")


mal.dat.thao %>%
  filter(sex %in% c("M","F")) %>% 
  ggplot(aes(log(sl),log(volume),col=sex))+geom_point()
  
mal.dat.thao %>%
  mutate(log_vol=log(volume),log_sl=log(sl),log_po=log(po)) %>% 
  filter(sex %in% c("M","F")) %>% 
  ggplot(aes(log_sl,log_po,col=sex))+geom_point()+geom_smooth(method="lm")

mal.dat.thao %>%
  mutate(log_vol=log(volume),log_sl=log(sl),log_po=log(po)) %>% 
  filter(sex %in% c("M","F")) %>% 
  ggplot(aes(log_sl,log_vol,col=sex))+geom_point()+geom_smooth(method="lm")

summary(aov(log(volume,base=10)~log(sl)+sex,mal.dat.thao))

summary(aov(log(volume,base=10)~log(sl)+sex,mal.dat.thao))

#setwd("cpk.data/mniger")
mal.files <- list.files("cpk.data/mniger",full.names =T )
mal.files <- mal.files[!grepl("POs",mal.files)]

cpk.mal.dat <- mal.dat.thao %>% select(mcz,sex,mass,sl,volume,po,jl) %>% rename(po_t=po,jl_t=jl) %>%  full_join(
mal.files %>% vroom(id = "file") %>% 
  mutate(file=gsub(".xls","",file)) %>% 
  rename(dim=`...1`) %>% 
  left_join(tibble(dim=1:4,name=c("eye","ao","po","po.l"))) %>% 
  select(-dim) %>% 
  mutate(mcz=gsub("(.*)_.*","\\1",basename(file)) %>% as.integer(),
         jl=gsub("(.*)_(.*).*","\\2",basename(file)) %>% as.numeric()) %>% 
  mutate(value=ifelse(name=="po.l",Length,Area)) %>% 
  select(-c(Area:Length)) %>% 
  pivot_wider(id_cols=c(file,mcz,jl),names_from = name,values_from = value),
  
) %>% mutate(po=ifelse(is.na(po),po_t,po))

  
cpk.mal.dat %>% 
  na.omit() %>% 
  filter(sex%in% c("M","F")) %>% 
  ggplot(aes(log(sl),log(po),col=sex))+geom_point()+geom_smooth(method="lm")+xlab(length.exp)+ylab(po.exp)
  

   
  

# data.names <- c("mcz","sex","sl","jl","mass","po","eye","volume")
# colnames(ech.dat)[1] <- colnames(phot.dat)[1]<- "mcz"
# ech.dat$ao <- NA
#all.dat <- rbind(data.frame(sp="E. barbatum",ech.dat[,data.names]),data.frame(sp="P. guernei",phot.dat[,data.names]),data.frame(sp="M. niger",mal.dat[,data.names]))

po.exp <-  expression(paste("log PO area. "," (",mm^3,")",sep=" "))
vol.exp <- expression(paste("log lens vol. "," (",mu,m^3,")",sep=" "))
length.exp <- expression(paste("log "," ",L[s]," (mm)",sep=" "))



# log.vol.p <- ggplot(all.dat, aes(x=sl,y=volume/10^9,col=sex))+geom_smooth(method="lm",se=F,alpha=0.5)+geom_point(size=4,alpha=0.5)+ylab(vol.exp)+theme_bw(20)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+xlab(length.exp)+scale_x_continuous(trans='log2',breaks=seq(0,30,10))+scale_y_continuous(trans='log2',breaks=seq(0.1,1,0.1))+facet_grid(.~sp)


#setwd("~/Google Drive/deep sea fishes - eyes")

#add flux
flux.dat <- read.csv("mes.case.flux.dat2.csv") %>% mutate(po_area=pi*(PO/2)^2)
flux.lm <- glm(log(flux)~po_area,flux.dat,family=gaussian)

flux.dat$pred <- predict(flux.lm) %>% exp
flux.dat %>% ggplot(aes(po_area,flux))+geom_point()+geom_point(aes(y=pred),col="red")

cpk.mal.dat <- cpk.mal.dat %>%na.omit() %>% 
  filter(sex%in% c("M","F")) %>% 
  mutate(flux=predict(flux.lm,newdata = data.frame(po_area=po)) %>% exp)


cpk.mal.dat %>%  
  ggplot(aes(log(sl),log(flux),col=sex))+geom_point()+geom_smooth(method="lm")

anova(lm(log(flux)~log(sl)+sex,cpk.mal.dat))

  
 
 


vol.lm <- lm(log.volume~log.sl+sex,mal.dat)
for(i in 1:100){try(vol.seg <- segmented(vol.lm,seg.Z=~log.sl,psi = 2.2));print(i)}

summary(aov(vol.lm))
AICc(vol.lm)
AICc(vol.seg)

#warrant (2000) model
#N=[EA^2/16r^2]e^-(alph*r)

#reduces to r*ln(r)=ln(sqrt(PA^2/N))
#so r=e^w(ln(sqrt(PA^2/N)))
#https://www.r-bloggers.com/getting-that-x-with-the-glog-function-and-lamberts-w/
#matmatica where a=alpha (attenuation); b is (E*A^2)/(16*N)
#https://www.wolframalpha.com/input/?i=solve+r%5E2*e%5E(a*r)%3Db+for+r
#install.packages("gsl")
library(gsl)
#lambert_W0() is w()


#where r is distance, atten is attentuation (0.05/m)
r.dist <- function(E,A,atten=0.05,N=5){
  if(is.na(A)){return(NA)}else{
  b <- (E*A^2)/(16*N)
r <- (2*lambertW0((atten*sqrt(b))/2))/atten
return(r)}
}

r.dist(E=10^10,A=0.0073) #Warrants (2000) estimation for pupil 7.3 mm and 10^10 flux



#diam funciton
diam <- function(vol){
  r <- (vol/(4/3*pi))^(1/3)
  return(r*2)
}

cpk.mal.dat$diam <- diam(cpk.mal.dat$vol/10^18)*5

#calculate
#cpk.mal.dat$pupil <- (sqrt(all.dat$eye/pi)*2)/1000




cpk.mal.dat$detect <- sapply(cpk.mal.dat$diam,function(x) r.dist(E=10^10,A=x))

cpk.mal.dat %>% 
  ggplot(aes(sl,detect,col=sex))+geom_point()+geom_smooth(method="lm")
anova(lm(detect~sl+sex,cpk.mal.dat))

sls <- seq(8,15,0.1)

fem.flux <- predict(lm(log(flux)~log(sl),cpk.mal.dat %>% filter(sex=="F")),newdata = data.frame(sl=sls))
male.flux <- predict(lm(log(flux)~log(sl),cpk.mal.dat %>% filter(sex=="M")),newdata = data.frame(sl=sls))

fem.diam <- predict(lm(log(diam)~log(sl),cpk.mal.dat %>% filter(sex=="F")),newdata = data.frame(sl=sls)) 
male.diam <- predict(lm(log(diam)~log(sl),cpk.mal.dat %>% filter(sex=="M")),newdata = data.frame(sl=sls)) 

mod_dat <- 
  tibble(sl=sls,fem_flux=fem.flux %>% exp,fem_diam=fem.diam %>% exp,male_flux=male.flux %>% exp,male_diam=male.diam %>% exp) %>% na.omit

mod_dat2 <- mod_dat %>% 
  group_by(sl) %>% 
  mutate(male_det=r.dist(E=fem_flux,A=male_diam),fem_det=r.dist(E=male_flux,A=fem_diam))  %>% 
  mutate(male_det2=r.dist(E=fem_flux,A=fem_diam)) %>% 
  mutate(det_ratio=male_det/fem_det,det_ratio2=male_det2/fem_det) %>% 
  mutate(improve=det_ratio-det_ratio2)

mod_dat2 %>% 
  ungroup() %>% 
  select(sl,det_ratio) %>% 
  mutate(male_sl=sl,female_sl=sl) %>% select(-sl) %>% expand.grid() %>% view
  ggplot(aes(male_sl,female_sl,fill=det_ratio),size=0.2)+geom_tile()

mod_dat3 <- expand_grid(male_sl=sls,female_sl=sls) %>% 
  left_join(
    mod_dat2 %>% select(sl,fem_det) %>% rename(female_sl=sl)
  ) %>% 
  left_join(
    mod_dat2 %>% select(sl,male_det) %>% rename(male_sl=sl)
  ) %>% 
  mutate(female_gap=fem_det-male_det,male_gap=male_det-fem_det,det_ratio=male_det/fem_det)
  
library(viridis)


bar_1 <- mod_dat3 %>% filter(round(det_ratio,1)==1) %>% mutate(female_sl=predict(lm(female_sl~poly(male_sl,2))))
bar_1.5 <- mod_dat3 %>% filter(plyr::round_any(det_ratio,1.5)==1.5) %>% mutate(female_sl=predict(lm(female_sl~poly(male_sl,2))))
bar_3 <- mod_dat3 %>% filter(plyr::round_any(det_ratio,3)==3) %>% mutate(female_sl=predict(lm(female_sl~poly(male_sl,2))))

bar_05 <- mod_dat3 %>% filter(plyr::round_any(det_ratio,0.5)==0.5) %>% mutate(female_sl=predict(lm(female_sl~poly(male_sl,3))))
pal <- wes_palette("Zissou1",n=21,type="continuous")

male_bars <- mod_dat3 %>% filter(round(male_gap,1)%in%seq(0,-125,-25)) %>% mutate(female_sl=predict(lm(female_sl~poly(male_sl,3))))
male_gap_p<- mod_dat3 %>% 
  ggplot(aes(male_sl,female_sl,fill=male_gap))+geom_tile()+scale_fill_viridis(option="H")

p2 <- mod_dat3 %>% 
  ggplot(aes(male_sl,female_sl,fill=female_gap))+geom_tile()+scale_fill_viridis(option="H")#+geom_line(data=bar_1)+geom_line(data=bar_05)+geom_line(data=bar_3)

p3 <- mod_dat3 %>% 
  ggplot(aes(male_sl,female_sl,fill=det_ratio))+geom_tile()+scale_fill_viridis(option="H")#+geom_line(data=bar_1)+geom_line(data=bar_05)+geom_line(data=bar_3)

library(cowplot)

plot_grid(p1,p2,p3,ncol=3)
  

wdmod_dat2 %>%
  pivot_longer(det_ratio:det_ratio2) %>% 
  ggplot(aes(sl,value,col=name))+geom_point()
  
mod_dat2 %>%
  ggplot(aes(sl,improve))+geom_point()

fem.det <- r.dist(E=10^10,A=x)
all.dat.sum <-ddply(all.dat,.(sp,sex),summarize,m=mean(detect,na.rm=T),se=se(detect))


detect.p <- ggplot(dat=subset(all.dat.sum,sex!="?" & sex!="" & sp!="E. barbatum"))+geom_point(aes(x=sex,y=m))+facet_grid(.~sp)+geom_errorbar(aes(x=sex,ymin=m-se,ymax=m+se))+theme_classic(15)

pdf("detect.p.pdf")
print(detect.p)
dev.off()

all.dat <- subset(all.dat,sp!="E. barbatum")
log.detect.p <- ggplot(all.dat, aes(x=sl,y=detect,col=sex))+geom_smooth(method="lm",se=F,alpha=0.5)+geom_point(size=4,alpha=0.5)+ylab("log det. distance (m)")+theme_bw(20)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+xlab(length.exp)+scale_x_continuous(trans='log2',breaks=seq(0,30,10))+scale_y_continuous(trans='log2')+facet_grid(.~sp,scales="free_y")

summary(aov(lm(detect~sex*sl,subset(all.dat,sp=="P. guernei" &sex !="?"  & sex !=""))))

pdf("log.detect.p")
print(log.detect.p)
dev.off()

