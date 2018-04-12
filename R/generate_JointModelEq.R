generate_JointModelEq <-
function(varMNAR,data){
  JointModelEq<-data.frame(matrix(rep(0,2*length(varMNAR)*dim(data)[2]),nrow=dim(data)[2]),
                          row.names=names(data))
  JointModelEq_colnames<-NULL
  for (i in 1:length(varMNAR)){
    JointModelEq_colnames<-c(JointModelEq_colnames,
                             c(paste(varMNAR[i],"_var_sel",sep=""),
                               paste(varMNAR[i],"_var_out",sep="")))}
  names(JointModelEq)<-JointModelEq_colnames
  return(JointModelEq)
}
