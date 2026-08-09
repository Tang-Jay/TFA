#===============================================
#                 Chen
#===============================================
lambdaChen<-function(u){# u为p*n,p是X的维数，n是数据的个数
  p=min(dim(u))
  M=rep(0,p)
  k=0
  gama=1
  tol=1e-11
  dif=1
 
  R1=rep(0,p)
  R2=R1%*%t(R1)
  while(dif>tol && k<=300){
    # 计算R1、R2
    aa=1+t(M)%*%t(u)
    for(i in 1:p){
      R1[i]=sum(t(u[,i])/aa)
      for(j in 1:p){
        R2[i,j]=-sum(u[,i]*u[,j]/aa^2)
      }
    }
    
    delta=-solve(R2)%*%R1
    dif=c(sqrt(t(delta)%*%delta))
    sigma=gama*delta
    while(min(1+t(M+sigma)%*%t(u))<=0){
      gama=gama/2
      sigma=gama*delta
    }
    
    M=M+sigma
    gama=1/sqrt(k+1)
    k=k+1
  }
  # cat(k,'\n')
  return(M)
} 
# 检验
# lam = lambdaChen(z)
# aa=1+t(lam)%*%t(z)
# glam=rowSums(t(z)/matrix(aa,p,n,byrow = TRUE))
# print(max(abs(glam)))
#===================================================#
#           当u为p*n时,p是X的维数，n是数据的个数
#===================================================#
# lambdaChen<-function(u){ 
#   p=dim(u)[1]
#   M=rep(0,p)
#   k=0
#   gama=1
#   tol=1e-11
#   dif=1
#   R<-function(lam){R0=sum(log(1+t(lam)%*%u));return(R0)}
#   
#   R1=rep(0,p)
#   R2=R1%*%t(R1)
#   while(dif>tol && k<=300){
#     # 计算R1、R2
#     aa=1+t(M)%*%u
#     for(i in 1:p){
#       R1[i]=sum(t(u[i,])/aa)
#       for(j in 1:p){
#         R2[i,j]=-sum(u[i,]*u[j,]/aa^2)
#       }
#     }
#     
#     delta=-solve(R2)%*%R1
#     dif=c(sqrt(t(delta)%*%delta))
#     sigma=gama*delta
#     while(min(1+t(M+sigma)%*%u)<=0){
#       gama=gama/2
#       sigma=gama*delta
#     }
#     
#     M=M+sigma
#     gama=1/sqrt(k+1)
#     k=k+1
#   }
#   # cat(k,'\n')
#   return(M)
# } 

# 检验
# lam = lambdaChen(z)
# aa=1+t(lam)%*%z
# glam=rowSums(z/matrix(aa,p,n,byrow = TRUE))
# print(max(abs(glam)))
