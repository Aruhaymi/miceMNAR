\name{mice.impute.hecknorm}
\alias{mice.impute.hecknorm}

\title{
Imputation by Heckman's model for continuous outcome with Missing Not At Random mechanism using one-step estimator
}
\description{
Imputes univariate missing continuous outcome using Heckman's model with a one-step maximum likelihood estimator.
}
\usage{
mice.impute.hecknorm(y, ry, x, JointModelEq, control, ...)
}

\arguments{
  \item{y}{
Incomplete data vector of length \code{n}.
}
  \item{ry}{
Vector of missing data pattern of length \code{n} (FALSE=missing, TRUE=observed).
}
  \item{x}{
Matrix (n x p) of complete covariates.
}
  \item{JointModelEq}{
Output \code{JointModelEq} provided by \code{\link[miceMNAR:MNARargument]{MNARargument}()} function.
}
  \item{control}{
Output \code{control} provided by \code{\link[miceMNAR:MNARargument]{MNARargument}()} function.
}
  \item{\dots}{
Other named arguments.
}
}
\details{
Imputation of continuous MNAR outcome using Heckman's model. This imputation model requires to define a selection and an outcome equation for the imputation model. The imputation algorithm consists of the following steps:
\enumerate{
\item Fit the one-step Heckman's model and compute the posterior distribution of parameters
\item Draw imputation parameters from their posterior distribution
\item Draw imputed values from their predictive distribution
}}
\value{
A vector of length \code{nmis} with imputations.
}
%\references{
%Galimard, J.E., Chevret, S., Curis, E. and Resche-Rigon, M. (2018). Heckman imputation models for binary or continuous MNAR missing outcomes and MAR missing predictors (submitted).
%}
\author{
Jacques-Emmanuel Galimard and Matthieu Resche-Rigon
}

\references{
Galimard, J.E., Chevret, S., Curis, E., and Resche-Rigon, M. (2018). Heckman imputation models for binary or continuous MNAR missing outcomes and MAR missing predictors. BMC Medical Research Methodology (In press).
}

\seealso{
\code{\link{miceMNAR}}
}

\examples{

require("GJRM")
require(mvtnorm)

# Generation of a simulated dataset with MNAR mechanism on a continuous outcome

X1 <- rnorm(500,0,1)
X2 <- rbinom(500,1,0.5)
X3 <- rnorm(500,1,0.5)
  
errors <- rmvnorm(500,mean=c(0,0),sigma=matrix(c(1,0.3,0.3,1),nrow=2,byrow=TRUE))

Y <- X1+X2+errors[,1]
Ry <- ifelse(0.66+1*X1-0.5*X2+X3+errors[,2]>0,1,0)

Y[Ry==0] <- NA
  
simul_data <- data.frame(Y,X1,X2,X3)

JointModelEq <- generate_JointModelEq(data=simul_data,varMNAR = "Y")

JointModelEq[,"Y_var_sel"] <- c(0,1,1,1)
JointModelEq[,"Y_var_out"] <- c(0,1,1,0)

arg <- MNARargument(data=simul_data,varMNAR="Y",JointModelEq=JointModelEq)

\dontrun{imputation2 <- mice(data = arg$data_mod,
                 method = arg$method,
                 predictorMatrix = arg$predictorMatrix,
                 JointModelEq=arg$JointModelEq,
                 control=arg$control,
                 maxit=1,m=5)

analysis2 <- with(imputation,lm(Y~X1+X2+X3))
result2 <- pool(analysis2)
summary(result2)}
}