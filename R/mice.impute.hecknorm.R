mice.impute.hecknorm <-
function(y, ry, x,JointModelEq, control, ...) {

  vecvec<-2
   if(dim(control)[2]>1){
   vecvec<-apply(control,MARGIN=2, FUN=function(control){sum(ifelse(all.equal(y[ry],as.numeric(control[!is.na(control)]))==TRUE,1,0),
        ifelse(all.equal(ry,!is.na(control))==TRUE,1,0))})}
  
  formula_s<-as.formula(paste("ry~",
                              paste(row.names(JointModelEq[JointModelEq[,paste(names(control[,vecvec==2,drop=FALSE]),
                                                                 "_var_sel",sep="")]==1,]),collapse="+"),sep=""))
  formula_o<-as.formula(paste("y~",
                              paste(row.names(JointModelEq[JointModelEq[,paste(names(control[,vecvec==2,drop=FALSE]),
                                                                 "_var_out",sep="")]==1,]),collapse="+"),sep=""))

  rm(control,vecvec,JointModelEq)

  
  # 1. Estimate the Heckman's model parameters
  res <- copulaSampleSel(list(formula_s, 
                               formula_o),data=cbind(ry,y,x),fp=TRUE)

  # 2. Draw q.star
  q.star<-rmvnorm(1, res$coefficients, res$Vb, method = "chol")

  # 3. Calculate y.star
  # 3.a. Calculate linear predictors for selection equation with q.star
  var_sel<-res$X1[!ry,]
  c1<-dim(var_sel)[2]
  betaxstar_sel<-apply(var_sel,MARGIN=1,FUN=function(var_sel){sum(q.star[1:c1]*var_sel)})

  # 3.b. Calculate linear predictors for outcome equation with q.star
  var_out<-cbind(1,x[!ry,names(as.data.frame(res$X2))[2:dim(res$X2)[2]]])
  c2<-dim(var_out)[2]
  betaxstar_out<-apply(var_out,MARGIN=1,FUN=function(var_out){sum(q.star[(c1+1):(c1+c2)]*var_out)})
  
  # 3.c. In SemiParSampleSel rho and sigma are transformed, recalculate it
  transtheta<-function(x){(exp(2*x)-1)/(1+exp(2*x))}
  q.star[length(q.star)]<-pmin(100,q.star[length(q.star)]) # transtheta(100)==1
  q.star[length(q.star)]<-pmax(-100,q.star[length(q.star)])# transtheta(-100)==-1
  rho.star<-transtheta(q.star[length(q.star)])
  sigma.star<-exp(q.star[,(length(q.star)-1)]/2)

  # 3.d. Calculate y.star
  y.star<-betaxstar_out+sigma.star*rho.star*(-dnorm(betaxstar_sel)/(1-pnorm(betaxstar_sel)))+rnorm(sum(!ry),0, sigma.star)

return(y.star)
}
