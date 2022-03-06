library(brms)
library(mvtnorm)
library(lme4)
library(bayesplot)
library(ggplot2)
library(MplusAutomation)
library(relimp)
library(data.table)
library(MASS)
study<-2
studys<-c("study1_aftericc","study2","study3WLSMV","study4")
#wd<-"D:/study/current/组会/毕设/模拟数据/普通多层/"
wd<-"F:/舒方/current/本科毕设/重跑模拟/普通多层/"
#wd1<-paste(wd,"changevar/",studys[study],sep="")
wd1<-paste(wd,"every/",studys[study],sep="")
startrepi=1
totrepi=100
# Generating the group sizes. 
Lev2N <- c(20,30,40,50) #
Lev1N <- c(30,60,150)#
iccs<-c(0.05,0.1,0.2)#
iccLs<-c("005","01","02")#
if(study>=3){
  iccs<-c(0.1,0.2)#,
  iccLs<-c("01","02")#,,
}

esti<-c("ml","nonif","l_3sd_s","l_1sd_s","cor_s","r_1sd_s","r_3sd_s",#
"l_3sd_m","l_1sd_m","cor_m","r_1sd_m","r_3sd_m",
"l_3sd_w","l_1sd_w","cor_w","r_1sd_w","r_3sd_w")

if(study==4){
  esti[1]<-"wlsmv"
}

skew<-c("nor")
if(study==2){
  skew<-c("pos","neg")
}

#ref:Foldnes, N., & Gronneberg, S. (2017). 
#    The Asymptotic Covariance Matrix and its Use in Simulation Studies.
fleishman1978_abcd <- function (var, skewness, kurtosis) { 
  system.function <- function (x, skewness, kurtosis) { 
    b.=x[1L]
    c.=x[2L]
    d.=x[3L] 
    eq1 <- b.*b. + 6*b.*d. + 2*c.*c. + 15*d.*d. - var 
    eq2 <- 2*c.*(b.*b. + 24*b.*d. + 105*d.*d. + 2) - skewness 
    eq3 <- 24*(b.*d. + c.*c.*(1 + b.*b. + 28*b.*d.) 
               + d.*d.*(12 + 48*b.*d. + 141*c.*c. + 225*d.*d.)) - kurtosis 
    eq <- c(eq1,eq2,eq3) 
    sum (eq*eq) ## SS
  }
  out <- nlminb(start=c(1,0,0), objective=system.function, 
                scale=10, 
                control=list (trace =0), 
                skewness=skewness, kurtosis=kurtosis) 
  #if (out $ covergence!=0) warning ("no convergence") 
  b.<-out $ par[1L]
  c. <- out $ par [2L];d.<-out $ par[3L]
  a. <- -c.
  print(out$ objective)
  if (out$ objective > 10^{-5}) { 
    warning("! ! !Not valid combination of skew/kurtosis") 
  }
  return (c(a., b., c., d.)) 
}

effects<-c("ef")#,"zeroef"
gamma_01s = c(1.5)#,0
gamma_00 = 0  
gamma_10 = 2.3  
gamma_20 = -0.25 
gamma_30 = -0.75
gamma_11 = 0.12  
gamma_21 = -0.05  
gamma_31 = -0.75  

probw = 0.5 
var_w=probw*(1-probw)
x1_mean=0 
x1_sd=1 
var_x1=x1_sd^2
x2_pro1=0.5 
var_x2=x2_pro1*(1-x2_pro1)
x3_pro1=0.75
var_x3=x3_pro1*(1-x3_pro1)

#随机斜率的残差方差
s00<-0.5

set.seed(12345)

#for(ef in 1:length(effects)){
ef=1
gamma_01<-gamma_01s[ef]
effect<-effects[ef]

#先验信息的分布设置：
#方差分别为真值的0.1,0.2,0.5
varper<-c(0.1,0.2,0.5)
#均值分别偏离真实值3个、1个标准差，以及无偏。标准差对应改变。
errsd<-c(-3,-1,0,1,3)
#只对回归系数设置先验信息分布
setprior<-c("b01","b11","b21","b31","b10","b20","b30")
pror_true<-c(gamma_01,gamma_11,gamma_21,gamma_31,gamma_10,gamma_20,gamma_30)
pror_trues<-matrix(rep(pror_true,length(varper)*length(errsd)),
                   length(varper)*length(errsd),byrow = T)
prior_sd<-sqrt(abs(rep(varper,length(setprior))*rep(pror_true,each=length(varper))))
prior_sds<-matrix(rep(matrix(prior_sd,length(setprior)),each=length(errsd)),
                  ncol=length(setprior))
prior_vars<-prior_sds^2
prior_err<-matrix(rep(errsd,length(varper)*length(setprior)),ncol=length(setprior))
prior_means<-pror_trues+prior_err*prior_sds

#==============画先验分布图===============
strengths<-rep(c("10%Variance","20%Variance","50%Variance"),each=5)
mags<-rep(c("-3SD","-1SD","0SD","+1SD","+3SD"),3)
setwd("D:/study/current/毕设/中文修改/图片/prior exam")
for(i in 1:15){#i<-1
  stre<-strengths[i]
  title=mag=mags[i]
  x<-rnorm(1000000,prior_means[i],prior_sds[i])
  tiff(file=paste0(i,".tiff"), width=6, height=6, units="cm", res=300)
  pic<-ggplot(,aes(x = x))+geom_density()+
    xlab("b01")+
    ggtitle(title)+
    theme(plot.title = element_text(size=16, hjust = 0.5),
          axis.title = element_text(size=10),
          axis.text = element_text(size=5))+
    geom_vline(xintercept=gamma_01,linetype=2)+
    xlim(-4,8)
  print(pic)
  dev.off()
}
#=============画先验分布图结束===============

#通过改变残差来确保组内方差为希望得到的取值
if(study<=2){
  var_tot<-13.801   #参考McNeish (2016)
  var_ybs<-var_tot*iccs
  var_yws<-var_tot-var_ybs
  yj_resids<-var_ybs-gamma_01^2*var_w
  yij_resids<-var_yws-(gamma_10^2*var_x1+gamma_11^2*(var_w+var_x1)+(s00+var_x1)+
                         gamma_20^2*var_x2+gamma_21^2*(var_w+var_x2)+(s00+var_x2)+
                         gamma_30^2*var_x3+gamma_31^2*(var_w+var_x3)+(s00+var_x3))
}

if(study>=3){
  yij_resids=rep(1,length(iccs))
  if(study==3){
    var_yws<-yij_resids+(gamma_10^2*var_x1+gamma_11^2*(var_w+var_x1)+(s00+var_x1)+
                           gamma_20^2*var_x2+gamma_21^2*(var_w+var_x2)+(s00+var_x2)+
                           gamma_30^2*var_x3+gamma_31^2*(var_w+var_x3)+(s00+var_x3))
  }
  if(study==4){
    var_yws<-yij_resids+gamma_10^2*var_x1+gamma_20^2*var_x2+gamma_30^2*var_x3
  }
  var_ybs=var_yws*iccs/(1-iccs)
  yj_resids<-var_ybs-gamma_01^2*var_w
}

setwd(wd)
source('functions_every.R')
setwd(wd1)
if(study<4){
  temp<-createMultiInps(study)
  cat(temp, file = paste(wd1,"/createModels.txt", sep = ''), append = F)
  createModels("createModels.txt")
}
if(study==4){
  temp<-createStudy4()
  cat(temp, file = paste(wd1,"/createStudy4.txt", sep = ''), append = F)
  createModels("createStudy4.txt")
}

#生成数据
for (ic in 1:length(iccs)){#ic=1
  icc<-iccs[ic]
  yj_resid<-yj_resids[ic]
  if(study==4){
    gamma_11 = gamma_21 = gamma_31 = 0
    s00<-0
  }
  theta<-diag(c(yj_resid,rep(s00,3)))
  yij_resid<-yij_resids[ic]
  sigma=sqrt(yij_resid)
  
  #for skew=-2(-c,-a)/2(c,a),kurtosis=7:
  fleish_b<-fleishman1978_abcd(yj_resid,2,7)
  fleish_w<-fleishman1978_abcd(yij_resid,2,7)
  
  for (i in 1:length(Lev2N)){#i=j=repi=1
    for (j in 1:length(Lev1N)){
      for (s in skew){#s="pos"
        with_sample_cross<-rep(Lev1N[j], Lev2N[i])
        path0=paste(wd1,"/",effect,"/",s,"/",iccLs[ic],"/",Lev2N[i],"/",Lev1N[j],"/",sep="")
        datname0=paste("data ",effect,"-",s,"-",icc,"-",Lev2N[i],"-",Lev1N[j],sep="")
        c1<-seq(startrepi:totrepi)
        datalist<-paste(datname0,"_",c1,".dat",sep="")
        write.table(datalist,paste(path0,datname0,"_list.dat",sep=""),
                    row.names=F,col.names=F,quote=F)
        com_iccs<-rep(NA,totrepi)
        for (repi in startrepi:totrepi){
          com_iccs[repi]<-simi(i,j,study)
        }
        #确认ICC大致取值
        print(c(icc,mean(com_iccs)))
        write.table(com_iccs,paste(path0,datname0," ",startrepi,"-",totrepi,"_icc.dat",sep=""),
                    row.names=F,col.names=F,quote=F)
      }
    }
  }
}

runModels(wd1,recursive=TRUE, showOutput=FALSE)
