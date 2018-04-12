mice.impute.hecknorm2step <-
function(y, ry, x,JointModelEq, control, ...) {

  vecvec<-2
  
  if(dim(control)[2]>1){vecvec<-apply(control,MARGIN=2, FUN=function(control){
      sum(ifelse(all.equal(y[ry],as.numeric(control[!is.na(control)]))==TRUE,1,0),
                                  ifelse(all.equal(ry,!is.na(control))==TRUE,1,0))})}
  
  var_sel<-row.names(JointModelEq[JointModelEq[,paste(names(control[,vecvec==2,drop=FALSE]),"_var_sel",sep="")]==1,])
  formula_s<-as.formula(paste("ry~",paste(var_sel,collapse="+"),sep=""))
  
  var_out<-row.names(JointModelEq[JointModelEq[,paste(names(control[,vecvec==2,drop=FALSE]),"_var_out",sep="")]==1,])
  formula_o<-as.formula(paste("y~",paste(var_out,collapse="+"),sep=""))
  nv1<-2+length(var_sel) ; nv2<-1+nv1+length(var_out)
  
  rm(control,vecvec,JointModelEq)
  
  heckit<-heckit2fit(formula_s, formula_o,data=cbind(ry,y,x))
  predlin<-linearPredictors(heckit$probit)
  imr <- ifelse(ry==1,dnorm(predlin)/pnorm(predlin),-dnorm(predlin)/pnorm(-predlin))
  delta0<-imr *(imr + predlin)

  x<-as.matrix(cbind(1,x[,var_out],imr))
  parm <- .norm.draw2(y, ry, x, vcovcor=heckit$vcov[nv1:nv2,nv1:nv2]/(heckit$sigma^2),
                                correct=1-(pmin(0.99999,(coef(heckit)["rho"])^2))*delta0[ry] ,...)
  
  return(x[!ry, ] %*% parm$beta + rnorm(sum(!ry))* parm$sigma*sqrt(1-(pmin(0.99999,(coef(heckit)["rho"])^2))*delta0[!ry])) 
}
