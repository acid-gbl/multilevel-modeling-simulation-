simi<-function(i,j,study){
  u_js<-mvrnorm(Lev2N[i],c(0,0,0,0),theta)
  if(study==2){
    a<-fleish_b[1]
    b<-fleish_b[2]
    c<-fleish_b[3]
    d<-fleish_b[4]
    x <- rnorm(Lev2N[i],mean = 0,sd =1)
    if(s=="pos"){
      temp<-a + b*x + c*x^2 + d*x^3
    }
    if(s=="neg"){
      temp<-c + b*x + a*x^2 + d*x^3
    }
    u_js[,1]<-temp
  }
  w1js = rbinom(Lev2N[i],1,probw)
  #btw var
  temp<-gamma_01*w1js + u_js[,1]
  variance_b<-var(temp)
  
  data_holder_cross_list <-  rep(list(list()), Lev2N[i])
  variance_W<-rep(NA,Lev2N[i])
  
  for (GroupID in 1:Lev2N[i]){
    GroupSize=with_sample_cross[GroupID]
    u0j <- u_js[GroupID,1]
    u1j <- u_js[GroupID,2]
    u2j <- u_js[GroupID,3]
    u3j <- u_js[GroupID,4]
    w1j <- w1js[GroupID]
    
    x1 <- rnorm(n = GroupSize,mean = x1_mean, sd = x1_sd)
    x2 = rbinom(GroupSize,1,x2_pro1)
    x3 = rbinom(GroupSize,1,x3_pro1)
    
    MemID <- 1:GroupSize
    e_ij = rnorm(GroupSize,mean = 0,sd =sigma)
    
    if(study==2){
      a<-fleish_w[1]
      b<-fleish_w[2]
      c<-fleish_w[3]
      d<-fleish_w[4]
      x <- rnorm(GroupSize,mean = 0,sd =1)
      if(s=="pos"){
        e_ij<-a + b*x + c*x^2 + d*x^3
      }
      if(s=="neg"){
        e_ij<-c + b*x + a*x^2 + d*x^3
      }
      e_ij<-e_ij
    }
    ### Generating the outcome Variable
    b00=gamma_00 + gamma_01*w1j + u0j #相当于Yb
    b1j=gamma_10 + gamma_11*w1j + u1j
    b2j=gamma_20 + gamma_21*w1j + u2j
    b3j=gamma_30 + gamma_31*w1j + u3j
    Y = b00 + b1j*x1 + b2j*x2 + b3j*x3 + e_ij 
    #wthin var
    variance_W[GroupID]<-var(Y-b00)
    
    if(study>=3){
      Y01<-rep(0,length(Y))
      Y01[Y>0]<-1
      Y<-Y01
    }
    ### Placing each group's data into an element of the list    
    data_holder_cross_list[[GroupID]] = data.frame(GroupID, MemID, w1j, x1, x2, x3, 
                                                   Y=Y)
  }
  
  ### Binding all groups of the data together
  data_holder_cross_centered <- do.call(rbind, data_holder_cross_list)
  data_holder_cross_centered <- as.data.frame(data_holder_cross_centered)
  #write.csv(data_holder_cross_centered, file = "final_data.csv", row.names = FALSE)
  
  datname=paste(path0,datname0,"_",repi,".dat",sep="")
  prepareMplusData(data_holder_cross_centered, datname,inpfile =FALSE) 
  
  #icc
  icc<-variance_b/(variance_b+mean(variance_W))
  return(icc)
}
### generate Mplus inputs==========================================================
prior<-""
#可能舍去后面好些有信息的贝叶斯
if(length(esti)>2){ 
  prior<-c(paste("priorM_",setprior[1],"#esti=",sep=''),rep("999",2),prior_means[1:(length(esti)-2),1],";\n",
           paste("priorV_",setprior[1],"#esti=",sep=''),rep("999",2), prior_vars[1:(length(esti)-2),1],";\n",
           paste("priorM_",setprior[2],"#esti=",sep=''),rep("999",2), prior_means[1:(length(esti)-2),2],";\n",
           paste("priorV_",setprior[2],"#esti=",sep=''),rep("999",2), prior_vars[1:(length(esti)-2),2],";\n",
           paste("priorM_",setprior[3],"#esti=",sep=''),rep("999",2), prior_means[1:(length(esti)-2),3],";\n",
           paste("priorV_",setprior[3],"#esti=",sep=''),rep("999",2), prior_vars[1:(length(esti)-2),3],";\n",
           paste("priorM_",setprior[4],"#esti=",sep=''),rep("999",2),prior_means[1:(length(esti)-2),4],";\n",
           paste("priorV_",setprior[4],"#esti=",sep=''),rep("999",2), prior_vars[1:(length(esti)-2),4],";\n",
           paste("priorM_",setprior[5],"#esti=",sep=''),rep("999",2),prior_means[1:(length(esti)-2),5],";\n",
           paste("priorV_",setprior[5],"#esti=",sep=''),rep("999",2), prior_vars[1:(length(esti)-2),5],";\n",
           paste("priorM_",setprior[6],"#esti=",sep=''),rep("999",2),prior_means[1:(length(esti)-2),6],";\n",
           paste("priorV_",setprior[6],"#esti=",sep=''),rep("999",2), prior_vars[1:(length(esti)-2),6],";\n",
           paste("priorM_",setprior[7],"#esti=",sep=''),rep("999",2),prior_means[1:(length(esti)-2),7],";\n",
           paste("priorV_",setprior[7],"#esti=",sep=''),rep("999",2), prior_vars[1:(length(esti)-2),7],";\n")
}
#可能舍去ML或者无信息的贝叶斯,或者都舍去，或者都不舍
if(length(esti)==15|length(esti)==16|length(esti)==17){ 
  prior<-c(paste("priorM_",setprior[1],"#esti=",sep=''),rep("999",length(esti)-nrow(prior_means)),prior_means[,1],";\n",
           paste("priorV_",setprior[1],"#esti=",sep=''),rep("999",length(esti)-nrow(prior_means)), prior_vars[,1],";\n",
           paste("priorM_",setprior[2],"#esti=",sep=''),rep("999",length(esti)-nrow(prior_means)), prior_means[,2],";\n",
           paste("priorV_",setprior[2],"#esti=",sep=''),rep("999",length(esti)-nrow(prior_means)), prior_vars[,2],";\n",
           paste("priorM_",setprior[3],"#esti=",sep=''),rep("999",length(esti)-nrow(prior_means)), prior_means[,3],";\n",
           paste("priorV_",setprior[3],"#esti=",sep=''),rep("999",length(esti)-nrow(prior_means)), prior_vars[,3],";\n",
           paste("priorM_",setprior[4],"#esti=",sep=''),rep("999",length(esti)-nrow(prior_means)),prior_means[,4],";\n",
           paste("priorV_",setprior[4],"#esti=",sep=''),rep("999",length(esti)-nrow(prior_means)), prior_vars[,4],";\n",
           paste("priorM_",setprior[5],"#esti=",sep=''),rep("999",length(esti)-nrow(prior_means)),prior_means[,5],";\n",
           paste("priorV_",setprior[5],"#esti=",sep=''),rep("999",length(esti)-nrow(prior_means)), prior_vars[,5],";\n",
           paste("priorM_",setprior[6],"#esti=",sep=''),rep("999",length(esti)-nrow(prior_means)),prior_means[,6],";\n",
           paste("priorV_",setprior[6],"#esti=",sep=''),rep("999",length(esti)-nrow(prior_means)), prior_vars[,6],";\n",
           paste("priorM_",setprior[7],"#esti=",sep=''),rep("999",length(esti)-nrow(prior_means)),prior_means[,7],";\n",
           paste("priorV_",setprior[7],"#esti=",sep=''),rep("999",length(esti)-nrow(prior_means)), prior_vars[,7],";\n")
}
##create single inp for each rep
createMultiInps <- function(study){
  init <- c("[[init]]\n",
            "iterators = ef L2N L1N icc esti skew rep;\n",
            "L2N = ", Lev2N,";\n",
            "L1N = ", Lev1N,";\n",
            paste("rep = ",startrepi,":",totrepi,";\n",sep=""),
            paste("esti = 1:",length(esti),";\n",sep=""),
            paste("ef = 1:",length(effects),";\n",sep=""),
            "efL#ef = ", effects,";\n",
            "efV#ef = ", gamma_01s,";\n",
            paste("icc = 1:",length(iccs),";\n",sep=""),
            "iccV#icc = ", iccs,";\n",
            "iccL#icc = ", iccLs,";\n",
            "yij_resid#icc = ", yij_resids,";\n",
            "yj_resid#icc = ", yj_resids,";\n",
            paste("skew = 1:",length(skew),";\n",sep=""),
            prior,
            "estiID#esti= ",esti,";\n",
            "skewID#skew=",skew,";\n",
            "filename = [[efL#ef]]-[[skewID#skew]]-[[iccL#icc]]-[[L2N]]-[[L1N]]-[[estiID#esti]]-[[rep]].inp;\n",
            "outputDirectory = [[efL#ef]]/[[skewID#skew]]/[[iccL#icc]]/[[L2N]]/[[L1N]];\n",
            "[[/init]]\n",
            "\n")
  
  cate<-c(""," CATEGORICAL= Y;\n"," [Y$1*0];\n",
          " INTEGRATION=MONTECARLO;\n  LINK=PROBIT;\n")
  if(study<3){
    cate<-c(" Y*[[yij_resid#icc]];\n",""," [Y*0];\n","")
  }
  body <- c("TITLE:\n",
            "[[L2N]]_[[L1N]]_[[estiID#esti]]\n",
            "\n",
            "DATA:\n",
            " FILE=data [[efL#ef]]-[[skewID#skew]]-[[iccV#icc]]-[[L2N]]-[[L1N]]_[[rep]].dat;\n",
            " !TYPE = MONTECARLO;\n",
            "\n",
            "VARIABLE: \n",
            " NAMES = GroupID MemID w1j x1 x2 x3 Y; \n",
            " USEVAR= GroupID w1j x1 x2 x3 Y;\n",
            " CLUSTER = GroupID;\n",
            " WITHIN = x1 x2 x3;\n",
            " BETWEEN = w1j;\n",
            cate[2],
            "\n",
            "ANALYSIS:\n",
            " TYPE IS TWOLEVEL RANDOM;\n",
            "[[esti==1]]\n",
            " ESTIMATOR = MLR;\n",
            cate[4],
            "[[/esti==1]]\n",
            "[[esti!=1]]\n",
            " ESTIMATOR = BAYES;\n",
            " PROCESSORS = 2;\n",
            " BITERATIONS = 70000 (2000);\n",
            "[[/esti!=1]]\n",
            "\n")
  
  model <- c("MODEL:\n",
             " %WITHIN%\n",
             cate[1],
             " b1j| Y ON x1;\n",
             " b2j| Y ON x2;\n",
             " b3j| Y ON x3;\n",
             " \n",
             " %BETWEEN%\n",
             " Y b1j b2j b3j;\n",
             cate[3],
             " [b1j*2.3] (b10);\n",
             " [b2j*-0.25] (b20);\n",
             " [b3j*-0.75] (b30);\n",
             " Y ON w1j*[[efV#ef]](b01);\n",
             " b1j ON w1j*0.12(b11);\n",
             " b2j ON w1j*-0.05(b21);\n",
             " b3j ON w1j*-0.75(b31);\n",
             " Y*[[yj_resid#icc]];\n",
             "[[esti==1]]\n",
             " b1j-b3j*0.5;\n",
             "[[/esti==1]]\n",
             "[[esti!=1]]\n",
             " b1j-b3j*0.5;\n",
             "[[/esti!=1]]\n",
             "\n")
  
  prior <- c("[[esti>2]]\n",
             "MODEL PRIORS:\n",
             paste(setprior[1],"~N([[priorM_",setprior[1],"#esti]],[[priorV_",setprior[1],"#esti]]);\n",sep=''),
             paste(setprior[2],"~N([[priorM_",setprior[2],"#esti]],[[priorV_",setprior[2],"#esti]]);\n",sep=''),
             paste(setprior[3],"~N([[priorM_",setprior[3],"#esti]],[[priorV_",setprior[3],"#esti]]);\n",sep=''),
             paste(setprior[4],"~N([[priorM_",setprior[4],"#esti]],[[priorV_",setprior[4],"#esti]]);\n",sep=''),
             paste(setprior[5],"~N([[priorM_",setprior[5],"#esti]],[[priorV_",setprior[5],"#esti]]);\n",sep=''),
             paste(setprior[6],"~N([[priorM_",setprior[6],"#esti]],[[priorV_",setprior[6],"#esti]]);\n",sep=''),
             paste(setprior[7],"~N([[priorM_",setprior[7],"#esti]],[[priorV_",setprior[7],"#esti]]);\n",sep=''),
             "[[/esti>2]]\n",
             "\n",
             "OUTPUT:\n",
             " TECH1 TECH8 CINT;! TECH9\n")
  
  template <- c(init, body, model, prior, sep='')
  return(template)
}

##study4: categorical & no random slope
createStudy4 <- function(){
  init <- c("[[init]]\n",
            "iterators = ef rep L2N L1N icc esti skew;\n",
            paste("skew = 1:",length(skew),";\n",sep=""),
            "L2N = ", Lev2N,";\n",
            "L1N = ", Lev1N,";\n",
            paste("esti = 1:",length(esti),";\n",sep=""),
            paste("ef = 1:",length(effects),";\n",sep=""),
            paste("rep = 1:",totrepi,";\n",sep=""),
            "efL#ef = ", effects,";\n",
            "efV#ef = ", gamma_01s,";\n",
            paste("icc = 1:",length(iccs),";\n",sep=""),
            "iccV#icc = ", iccs,";\n",
            "iccL#icc = ", iccLs,";\n",
            "yij_resid#icc = ", yij_resids,";\n",
            "yj_resid#icc = ", yj_resids,";\n",
            "skewID#skew=",skew,";\n",
            paste("priorM_",setprior[1],"#esti=  999 999 ",sep=''),prior_means[,1],";\n",
            paste("priorV_",setprior[1],"#esti=  999 999 ",sep=''), prior_vars[,1],";\n",
            paste("priorM_",setprior[5],"#esti=  999 999 ",sep=''),prior_means[,5],";\n",
            paste("priorV_",setprior[5],"#esti=  999 999 ",sep=''), prior_vars[,5],";\n",
            paste("priorM_",setprior[6],"#esti=  999 999 ",sep=''),prior_means[,6],";\n",
            paste("priorV_",setprior[6],"#esti=  999 999 ",sep=''), prior_vars[,6],";\n",
            paste("priorM_",setprior[7],"#esti=  999 999 ",sep=''),prior_means[,7],";\n",
            paste("priorV_",setprior[7],"#esti=  999 999 ",sep=''), prior_vars[,7],";\n",
            "estiID#esti= ",esti,";\n",
            "filename = [[efL#ef]]-[[skewID#skew]]-[[iccL#icc]]-[[L2N]]-[[L1N]]-[[estiID#esti]]-[[rep]].inp;\n",
            "outputDirectory = [[efL#ef]]/[[skewID#skew]]/[[iccL#icc]]/[[L2N]]/[[L1N]];\n",
            "[[/init]]\n",
            "\n")
  
  body <- c("TITLE:\n",
            "[[L2N]]_[[L1N]]_[[estiID#esti]]\n",
            "\n",
            "DATA:\n",
            " FILE=data [[efL#ef]]-[[skewID#skew]]-[[iccV#icc]]-[[L2N]]-[[L1N]]_[[rep]].dat;\n",
            "!TYPE = MONTECARLO;\n",
            "\n",
            "VARIABLE: \n",
            " NAMES = GroupID MemID w1j x1 x2 x3 Y;\n",
            " USEVAR= GroupID w1j x1 x2 x3 Y;\n",
            " CLUSTER = GroupID;\n",
            " WITHIN = x1 x2 x3;\n",
            " BETWEEN = w1j;\n",
            " CATEGORICAL= Y;\n",
            "\n",
            "ANALYSIS:\n",
            " TYPE IS TWOLEVEL;\n",
            "[[esti==1]]\n",
            " ESTIMATOR = WLSMV;\n",
            " INTEGRATION=MONTECARLO;\n",
            "[[/esti==1]]\n",
            "[[esti!=1]]\n",
            " ESTIMATOR = BAYES;\n",
            " PROCESSORS = 2;\n",
            " BITERATIONS = 70000 (2000);\n",
            "[[/esti!=1]]\n",
            "\n")
  
  model <- c("MODEL:\n",
             " %WITHIN%\n",
             " Y ON x1*2.3(b10);\n",
             " Y ON x2*-0.25(b20);\n",
             " Y ON x3*-0.75(b30);\n",
             " \n",
             " %BETWEEN%\n",
             " Y;\n",
             " [Y$1*0];\n",
             " Y ON w1j*1.5(b01);\n",
             " Y*[[yj_resid#icc]];\n",
             "\n")
  
  prior <- c("[[esti>1]]\n",
             "MODEL PRIORS:\n",
             paste(setprior[1],"~N([[priorM_",setprior[1],"#esti]],[[priorV_",setprior[1],"#esti]]);\n",sep=''),
             paste(setprior[5],"~N([[priorM_",setprior[5],"#esti]],[[priorV_",setprior[5],"#esti]]);\n",sep=''),
             paste(setprior[6],"~N([[priorM_",setprior[6],"#esti]],[[priorV_",setprior[6],"#esti]]);\n",sep=''),
             paste(setprior[7],"~N([[priorM_",setprior[7],"#esti]],[[priorV_",setprior[7],"#esti]]);\n",sep=''),
             "[[/esti>1]]\n",
             "\n",
             "OUTPUT:\n",
             " TECH1 TECH8 CINT;! TECH9\n")
  
  template <- c(init, body, model, prior, sep='')
  return(template)
}