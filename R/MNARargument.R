MNARargument <-
function(data, method = NULL, predictorMatrix = NULL,varMNAR, JointModelEq=NULL){
  
  
  control<-data[,varMNAR,drop=FALSE]
  if(is.null(predictorMatrix)){
    predictorMatrix<-mice(data,maxit=0)$predictorMatrix
  }
  if(is.null(method)){
    method<-mice(data,maxit=0)$method
  }
  
  for(i in 1:length(varMNAR)){
    
    if(sum((JointModelEq[,paste(varMNAR[i],"_var_out",sep="")]-
            JointModelEq[,paste(varMNAR[i],"_var_sel",sep="")])^2)==0){
      warning("It is preferable to specify different outcome and selection equations",call. = FALSE)}
    
    
### Modification of predictorMatrix #
    predictorMatrix<-cbind(predictorMatrix,pmin(1,mice(data,maxit=0)$nmis))
      colnames(predictorMatrix)[dim(predictorMatrix)[2]]<-paste("ind",varMNAR[i],sep="_")
    predictorMatrix[varMNAR[i],] <-1
    predictorMatrix[varMNAR[i],varMNAR[i]] <-0
    predictorMatrix[varMNAR[i],paste("ind",varMNAR[i],sep="_")] <-0 
    predictorMatrix<-rbind(predictorMatrix,0)
      row.names(predictorMatrix)[dim(predictorMatrix)[1]]<-paste("ind",varMNAR[i],sep="_")
    
### Modification of method #  
    method[varMNAR[i]]<-ifelse(is.numeric(control[,varMNAR[i]]) & length(levels(as.factor(control[,varMNAR[i]])))>2,"hecknorm","heckprob")
    method<-c(method,"")
      names(method)[length(method)]<-paste("ind",varMNAR[i],sep="_")

### Modification of the dataset #    
    if(length(levels(as.factor(control[,varMNAR[i]])))==2){
    data[,varMNAR[i]]<-as.factor(ifelse(as.character(data[,varMNAR[i]])==levels(as.factor(control[,varMNAR[i]]))[1],0,1))}
    data<-cbind(data,as.factor(!is.na(control[,i]))) 
      colnames(data)[dim(data)[2]]<-paste("ind",varMNAR[i],sep="_")
    }


  for(i in 1:length(varMNAR)){JointModelEq<-rbind(JointModelEq,0)
  row.names(JointModelEq)[dim(JointModelEq)[1]]<-paste("ind",varMNAR[i],sep="_")
  for(j in 1:length(varMNAR)){
    JointModelEq[paste("ind",varMNAR[i],sep="_"),paste(varMNAR[j],"_var_out",sep="")]<-1
  }
  JointModelEq[paste("ind",varMNAR[i],sep="_"),paste(varMNAR[i],"_var_out",sep="")]<-0
  }

### Modification of JointModelEq #  
  for(i in 1:dim(data)[2]){
    if(is.factor(data[,i])==TRUE){
      row.names(JointModelEq)[i]<-paste0(colnames(data)[i],unlist(lapply(levels(data[,i])[2],FUN=function(x){gsub(" ","_",x, perl=T)})))
      if(length(levels(data[,i]))>2){
        for (j in 3:(length(levels(data[,i])))){
          JointModelEq<-rbind(JointModelEq,JointModelEq[i,])
          row.names(JointModelEq)[dim(JointModelEq)[1]]<-paste0(colnames(data)[i],unlist(lapply(levels(data[,i])[j],FUN=function(x){gsub(" ","_",x, perl=T)})))
        }
      }
    }
  }
  
  return(list(data_mod=data,method=method, predictorMatrix = predictorMatrix,JointModelEq=JointModelEq, control=control))
}
