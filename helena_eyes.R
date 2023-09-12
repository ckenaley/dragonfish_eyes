library(tidyverse)
library(vroom)
library(lamW)
library(ggsci)
library(viridis)
library(metR)





se <- function(x,na.rm=T){if(na.rm==T)
{x <- x[is.na(x)==F]}else{x <- x}
  sd(x)/sqrt(length(x))}

phot.dat <- read.csv("Photostomias_guernei2.csv")
qplot(log(mass),log(po),data=phot.dat,col=sex)


phot.dat %>% write_csv("Photostomias_guernei2.csv")


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

summary(aov(log(volume)~log(sl)+sex,mal.dat.thao))

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

  
phot.files <- list.files("cpk.data/pguernei",full.names =T )
phot.files <- phot.files[!grepl("POs",phot.files)]
cpk.phot.dat <- phot.dat %>% select(mcz,sex,mass,sl,volume,po,jl) %>% rename(po_t=po,jl_t=jl) %>% mutate(mcz2=gsub("_\\d","",mcz)) %>% full_join(
  phot.files %>% vroom(id = "file") %>% 
    mutate(file=gsub(".txt","",file)) %>% 
    rename(dim=`...1`) %>% 
    left_join(tibble(dim=1:4,name=c("eye","ao","po","po.l"))) %>% 
    select(-dim) %>% 
    mutate(mcz2=gsub("MCZ(\\d+)_.*","\\1",basename(file)),
           jl=gsub(".+guernei_(.+)JL.*","\\1",basename(file)) %>% as.numeric()) %>% 
    mutate(value=ifelse(name=="po.l",Length,Area)) %>% 
  select(-c(Area:Length)) %>% 
    pivot_wider(id_cols=c(file,mcz2,jl),names_from = name,values_from = value),
  
) %>% mutate(po=ifelse(is.na(po),po_t,po))

cpk.phot.dat %>% view


cpk.phot.dat %>% 
  na.omit() %>% 
  filter(!is.na(file)) %>% 
  filter(sex%in% c("M","F")) %>% 
  ggplot(aes(log(volume),log(eye),col=sex))+geom_point()+geom_smooth(method="lm")


mal_phot_p <- cpk.mal.dat %>% 
  na.omit() %>% 
  filter(sex%in% c("M","F")) %>% 
  ggplot(aes(log(sl),log(po),shape=sex))+geom_point()+geom_smooth(method="lm")+xlab("Standard Length (log cm)")+ylab(po.exp)+theme_classic(15) 
  ggsave("mal_po_size.pdf",mal_phot_p)
  
  mal_eye_p <- cpk.mal.dat %>% 
    na.omit() %>% 
    filter(sex%in% c("M","F")) %>% 
    ggplot(aes(log(sl),log(volume),shape=sex))+geom_point()+geom_smooth(method="lm")+xlab("Standard Length (log cm)")+ylab(vol.exp)+theme_classic(15) 
  ggsave("mal_eye_size.pdf",mal_eye_p)
  
all.dat <- rbind(
cpk.mal.dat %>% 
  select(-setdiff(colnames(cpk.mal.dat),colnames(phot.dat))) %>% 
  mutate(genus="Malacosteus"),
phot.dat %>% 
  select(-setdiff(colnames(phot.dat),colnames(cpk.mal.dat)))%>% 
  mutate(genus="Photostomias")
)


po.exp <-  expression(paste("log PO area. "," (",mm^3,")",sep=" "))
vol.exp <- expression(paste("lens vol. "," (",mu,m^3,")",sep=" "))
length.exp <- expression(paste("log "," ",L[s]," (mm)",sep=" "))



# log.vol.p <- ggplot(all.dat, aes(x=sl,y=volume/10^9,col=sex))+geom_smooth(method="lm",se=F,alpha=0.5)+geom_point(size=4,alpha=0.5)+ylab(vol.exp)+theme_bw(20)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+xlab(length.exp)+scale_x_continuous(trans='log2',breaks=seq(0,30,10))+scale_y_continuous(trans='log2',breaks=seq(0.1,1,0.1))+facet_grid(.~sp)


#setwd("~/Google Drive/deep sea fishes - eyes")

#add flux
flux.dat <- read.csv("mes.case.flux.dat2.csv") %>% mutate(po_area=pi*(PO/2)^2)
flux.lm <- glm(log(flux)~po_area,flux.dat,family=gaussian)

flux.dat$pred <- predict(flux.lm) %>% exp
flux.dat %>% ggplot(aes(po_area,flux))+geom_point()+geom_point(aes(y=pred),col="red")

all.dat <- all.dat  %>% 
  filter(sex%in% c("M","F")) %>% 
  mutate(flux=predict(flux.lm,newdata = data.frame(po_area=po)) %>% exp)


all.dat %>%  
  ggplot(aes(log(sl),log(flux),col=sex))+geom_point()+geom_smooth(method="lm")+facet_wrap(~genus)

anova(lm(log(flux)~log(sl)*sex+genus,all.dat))

flux.exp <- expression(paste('flux (',photons%*%s^{-1},')'))
bquote(X-Axis^superscript)
all.dat %>%  
  ggplot(aes(sl,flux,col=sex))+geom_point()+geom_smooth(method="lm")+facet_wrap(~genus,ncol=1,scales = "free")+theme_classic(15) +
  scale_y_continuous(trans='log10',labels = function(x) round(x,-10))+ylab(flux.exp)+scale_color_aaas(na.translate = F)+xlab("standard length (cm)")

all.dat %>%  
  ggplot(aes(sl,volume,col=sex))+geom_point()+geom_smooth(method="lm")+facet_wrap(~genus,ncol=1,scales = "free")+theme_classic(15) +
  scale_y_continuous(trans='log10',labels = function(x) round(x,-6))+ylab(vol.exp)+scale_color_aaas(na.translate = F)+xlab("standard length (cm)")
  
 
#warrant (2000) model
#N=[EA^2/16r^2]e^-(alph*r)

#reduces to r*ln(r)=ln(sqrt(PA^2/N))
#so r=e^w(ln(sqrt(PA^2/N)))
#https://www.r-bloggers.com/getting-that-x-with-the-glog-function-and-lamberts-w/
#matmatica where a=alpha (attenuation); b is (E*A^2)/(16*N)
#https://www.wolframalpha.com/input/?i=solve+r%5E2*e%5E(a*r)%3Db+for+r
#install.packages("gsl")
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

all.dat$diam <- diam(all.dat$vol/10^18)*5

#calculate
#cpk.mal.dat$pupil <- (sqrt(all.dat$eye/pi)*2)/1000




all.dat$detect <- sapply(all.dat$diam,function(x) r.dist(E=10^10,A=x))

all.dat %>% 
  ggplot(aes(sl,detect,col=sex))+geom_point()+geom_smooth(method="lm")+facet_wrap(.~genus)
anova(lm(detect~sl*sex+genus,all.dat))

sls <- seq(6,15,0.1)

mal.fem.flux <- predict(lm(log(flux)~log(sl),all.dat %>% filter(sex=="F",genus=="Malacosteus")),newdata = data.frame(sl=sls))
mal.male.flux <- predict(lm(log(flux)~log(sl),all.dat %>% filter(sex=="M",genus=="Malacosteus")),newdata = data.frame(sl=sls))

mal.fem.diam <- predict(lm(log(diam)~log(sl),all.dat %>% filter(sex=="F",genus=="Malacosteus")),newdata = data.frame(sl=sls)) 
mal.male.diam <- predict(lm(log(diam)~log(sl),all.dat %>% filter(sex=="M",genus=="Malacosteus")),newdata = data.frame(sl=sls)) 

mal.mod_dat <- 
  tibble(sl=sls,fem_flux=mal.fem.flux %>% exp,fem_diam=mal.fem.diam %>% exp,male_flux=mal.male.flux %>% exp,male_diam=mal.male.diam %>% exp) %>% na.omit

mal.mod_dat2 <- mal.mod_dat %>% 
  group_by(sl) %>% 
  mutate(male_det=r.dist(E=fem_flux,A=male_diam),fem_det=r.dist(E=male_flux,A=fem_diam))  %>% 
  mutate(male_det2=r.dist(E=fem_flux,A=fem_diam)) %>% 
  mutate(det_ratio=male_det/fem_det,det_ratio2=male_det2/fem_det) %>% 
  mutate(improve=det_ratio-det_ratio2)

mal.mod_dat2 %>% 
  ungroup() %>% 
  select(sl,det_ratio) %>% 
  mutate(male_sl=sl,female_sl=sl) %>% select(-sl) %>% expand.grid() %>% 
  ggplot(aes(male_sl,female_sl,fill=det_ratio),size=0.2)+geom_tile()

mal.mod_dat3 <- expand_grid(male_sl=sls,female_sl=sls) %>% 
  left_join(
    mal.mod_dat2 %>% select(sl,fem_det) %>% rename(female_sl=sl)
  ) %>% 
  left_join(
    mal.mod_dat2 %>% select(sl,male_det) %>% rename(male_sl=sl)
  ) %>% 
  mutate(female_gap=fem_det-male_det,male_gap=male_det-fem_det,det_ratio=male_det/fem_det)



mal.male_gap_p<- mal.mod_dat3 %>% 
  ggplot(aes(male_sl,female_sl,fill=male_gap))+geom_tile(alpha=0.85)+scale_fill_viridis("male detection gap (m)",option="H")+geom_contour2(aes(z=male_gap,label = ..level..),col="black",binwidth = 25)+theme_minimal(15)+theme(legend.position = "bottom",legend.key.width= unit(0.5, 'cm'))+coord_fixed()+xlab("Male Standard Length (cm)")+ylab("Female Standard Length (cm)")

ggsave("malacosteus_male_gap.pdf",mal.male_gap_p)


mal.female_gap_p<- mal.mod_dat3 %>% 
  ggplot(aes(male_sl,female_sl,fill=female_gap))+geom_tile(alpha=0.85)+scale_fill_viridis("female gap (m)",option="H")+geom_contour2(aes(z=female_gap,label = ..level..),col="black",binwidth = 25)+theme_minimal(9)+theme(legend.position = "bottom",axis.title.y = element_blank())+coord_fixed()

mal.det_ratio_p<- mal.mod_dat3 %>% 
  ggplot(aes(male_sl,female_sl,fill=det_ratio))+geom_tile(alpha=0.85)+scale_fill_viridis("detectio ratio",option="H")+geom_contour2(aes(z=det_ratio,label = ..level..),col="black",binwidth = .25)+theme_minimal(9)+theme(legend.position = "bottom",axis.title.y = element_blank())+coord_fixed()

mal.mod_dat3 %>% 
  filter(male_sl==8.0) %>% 
  ggplot(aes(female_sl,fem_det))+geom_point()

mal.mod_dat3 %>% 
  filter(female_sl==10.0|male_sl==8.0) %>% 
  ggplot(aes(male_sl,male_det))+geom_point(aes(female_sl,fem_det),col="red")


#photostomias


phot.fem.flux <- predict(lm(log(flux)~log(sl),all.dat %>% filter(sex=="F",genus=="Photostomias")),newdata = data.frame(sl=sls))
phot.male.flux <- predict(lm(log(flux)~log(sl),all.dat %>% filter(sex=="M",genus=="Photostomias")),newdata = data.frame(sl=sls))

phot.fem.diam <- predict(lm(log(diam)~log(sl),all.dat %>% filter(sex=="F",genus=="Photostomias")),newdata = data.frame(sl=sls)) 
phot.male.diam <- predict(lm(log(diam)~log(sl),all.dat %>% filter(sex=="M",genus=="Photostomias")),newdata = data.frame(sl=sls)) 

phot.mod_dat <- 
  tibble(sl=sls,fem_flux=phot.fem.flux %>% exp,fem_diam=phot.fem.diam %>% exp,male_flux=phot.male.flux %>% exp,male_diam=phot.male.diam %>% exp) %>% na.omit

phot.mod_dat2 <- phot.mod_dat %>% 
  group_by(sl) %>% 
  mutate(male_det=r.dist(E=fem_flux,A=male_diam),fem_det=r.dist(E=male_flux,A=fem_diam))  %>% 
  mutate(male_det2=r.dist(E=fem_flux,A=fem_diam)) %>% 
  mutate(det_ratio=male_det/fem_det,det_ratio2=male_det2/fem_det) %>% 
  mutate(improve=det_ratio-det_ratio2)

phot.mod_dat2 %>% 
  ungroup() %>% 
  select(sl,det_ratio) %>% 
  mutate(male_sl=sl,female_sl=sl) %>% select(-sl) %>% expand.grid() %>% 
  ggplot(aes(male_sl,female_sl,fill=det_ratio),size=0.2)+geom_tile()

phot.mod_dat3 <- expand_grid(male_sl=sls,female_sl=sls) %>% 
  left_join(
    phot.mod_dat2 %>% select(sl,fem_det) %>% rename(female_sl=sl)
  ) %>% 
  left_join(
    phot.mod_dat2 %>% select(sl,male_det) %>% rename(male_sl=sl)
  ) %>% 
  mutate(female_gap=fem_det-male_det,male_gap=male_det-fem_det,det_ratio=male_det/fem_det)


phot.male_gap_p<- phot.mod_dat3 %>% 
  ggplot(aes(male_sl,female_sl,fill=male_gap))+geom_tile(alpha=0.85)+scale_fill_viridis("male detection gap (m)",option="H")+geom_contour2(aes(z=male_gap,label = ..level..),col="black",binwidth = 5)+theme_minimal(15)+theme(legend.position = "bottom",legend.key.width= unit(0.5, 'cm'),legend.text = element_text(size=7))+coord_fixed()+xlab("Male Standard Length (cm)")+ylab("Female Standard Length (cm)")

ggsave("photostomias_male_gap.pdf",phot.male_gap_p)


phot.female_gap_p<- phot.mod_dat3 %>% 
  ggplot(aes(male_sl,female_sl,fill=female_gap))+geom_tile(alpha=0.85)+scale_fill_viridis("female gap (m)",option="H")+geom_contour2(aes(z=female_gap,label = ..level..),col="black",binwidth = 5)+theme_minimal(9)+theme(legend.position = "bottom",axis.title.y = element_blank())+coord_fixed()

phot.det_ratio_p<- phot.mod_dat3 %>% 
  ggplot(aes(male_sl,female_sl,fill=det_ratio))+geom_tile(alpha=0.85)+scale_fill_viridis("detectio ratio",option="H")+geom_contour2(aes(z=det_ratio,label = ..level..),col="black",binwidth = .5)+theme_minimal(9)+theme(legend.position = "bottom",axis.title.y = element_blank())+coord_fixed()

phot.mod_dat3 %>% 
  filter(male_sl==8.0) %>% 
  ggplot(aes(female_sl,fem_det))+geom_point()

phot.mod_dat3 %>% 
  filter(female_sl==10.0|male_sl==8.0) %>% 
  ggplot(aes(male_sl,male_det))+geom_point(aes(female_sl,fem_det),col="red")


library(cowplot)

mal_det_p <- plot_grid(mal.male_gap_p,mal.female_gap_p,mal.det_ratio_p,ncol=2,axis="t")

ggsave("malacosteus_detection.pdf",mal_det_p,)
  

mod_dat2 %>%
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

