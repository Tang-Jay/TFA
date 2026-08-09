#===============================================
#                    Compas  AL                #
#===============================================
epsilonG_CI<-function(group_specific_data,target=theta_P, a=-1, b=1, step=0.001, alpha=0.95){
  # 超参数赋值
  target <- theta_P
  sex_name <- unique(group_specific_data$sex)
  age_cat_name <- unique(group_specific_data$age_cat)
  
  # print(sex_name)
  if(length(sex_name)!=1){ sex_name = ''}
  if(length(age_cat_name)!=1){ age_cat_name = ''}
  
  # 计算Loss(已知Y_h=1,选出Y=1,则L=0,即L取Y)
  Li <- group_specific_data$Y
  n <- length(Li)
  m <- 1
  L <- matrix(Li, n, m)
  Theta <- matrix(theta_P, n, m)
  
  # 计算cut
  alpha <- alpha
  cut = exp(-qchisq(alpha,m)/2)
  
  # ---------------------开始模拟------------------- #
  lb = 0
  ub = 0
  times = 0
  elrMax_epsilonG = 0
  elRatio=c(0)
  
  epsilonGs = seq(a+0.01,b,step)
  for(epsilonG in epsilonGs){
    # 估计方程
    E <-  matrix(epsilonG, n, m, byrow = TRUE)
    g = L - Theta - E 
    
    # 计算EL值
    z = g
    lam = c(lambdaChen(z))
    npi = 1/(1+lam*z)
    elr = prod(npi)
    # el = 2*sum(log(1+t(lam)%*%t(z)))
    # (-2)*log(elr)==el
    if(elr>max(elRatio)){elrMax_epsilonG=epsilonG}
    elRatio = c(elRatio,elr)
    if(elr>=cut && times==0){lb=epsilonG;times=1}
    if(elr>=cut && times==1){ub=epsilonG}
  }
  
  # setEPS()
  # postscript(paste0('African-American','-',alpha,'-',sex_name,'-',age_cat_name,'.eps'), width=7, height=5)
  # 空白画布
  par(mfrow = c(1, 1))
  plot(epsilonGs,epsilonGs,
       xlim=c(a, b),ylim=c(0,1),col='white',
       yaxt = 'n', ann = F)
  axis(2, las = 1)
  title( main=paste0('alpha =',alpha,' ',sex_name,'  ',age_cat_name),
         xlab='epsilonG',ylab='elr')
  # legend('topright',legend=c('AELCP','  ELCP'),
         # col=c(col='aquamarine3','blue1'),
         # lty=c(1,1),bty='n',lwd=1.5)
  # cut图
  abline(h=cut,col='red')
  # text(a+0.01/(b-a),0.2,label = round(cut,3),col='red')
  
  # elr图
  points(epsilonGs,elRatio[-1],type='o',lty=1,lwd=0.3,cex=1,pch=20,col='black')
  abline(v=elrMax_epsilonG,col='blue1')
  abline(v=lb,col='blue1')
  abline(v=ub,col='blue1')
  text(elrMax_epsilonG+0.03/(b-a),0.95,label = paste0('el_h=',elrMax_epsilonG),col='blue1')
  text(ub+0.02/(b-a),0.2,label = paste0('ub=',round(ub,4)),col='blue1')
  text(lb-0.02/(b-a),0.2,label = paste0('lb=',round(lb,4)),col='blue1')

  # 计算置信区域长度
  ELAL = ub-lb
  
  cat('n ',' ep_h1',' lb ',' ub ',' ELAL','\n')
  cat(n,elrMax_epsilonG,lb,ub,ELAL,'\n')

  # P(Y=1|X=aa,Y_h=1) VS P(Y=1|X=ca,Y_h=1)
  mu_aa_Y <- mean(group_specific_data$Y)
  mu_ca_Y <- target
  cat('样本epsilonG', mu_aa_Y - mu_ca_Y )
  
}
