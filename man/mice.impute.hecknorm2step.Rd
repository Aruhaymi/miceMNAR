\name{mice.impute.hecknorm2step}
\alias{mice.impute.hecknorm2step}

\title{
Imputation by Heckman's model for continuous Missing Not At Random ouctome using a two-step estimator
}
\description{
Imputes univariate continuous Missing Not At Random (MNAR) outcome using Heckman's model with a two-step estimator.
}
\usage{
mice.impute.hecknorm2step(y, ry, x, JointModelEq, control, ...)
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
Output \code{JointModelEq} provided by \code{\link[miceMNAR:MNARargument]{MNARargument}}() function.
}
  \item{control}{
Output \code{control} provided by \code{\link[miceMNAR:MNARargument]{MNARargument}}() function.
}
  \item{\dots}{
Other named arguments.
}
}
\details{
Imputation of continuous MNAR outcome using two-step Heckman's model (Galimard, 2016). This imputation model requires to define a selection and an outcome equation. 
The imputation algorithm consists of the following steps:
\enumerate{
\item Compute Inverse Mill's Ratio (IMR) from the selection equation probit model
\item Include IMR as predictor in the imputation model
\item Draw imputation parameters using approximate proper imputation for the linear model and adding the Heckman variance correction as detailed in Galimard et al (2016)
\item Draw imputed values from their predictive distribution
}
}
\value{
A vector of length \code{nmis} with imputations.
}
\references{
Galimard, J.-E., Chevret, S., Protopopescu, C., and Resche-Rigon, M. (2016) A multiple imputation approach for MNAR mechanisms compatible with Heckman's model. Statist. Med., 35: 2907-2920. doi:10.1002/sim.6902.
%
%Galimard, J.E., Chevret, S., Curis, E. and Resche-Rigon, M. (2018). Heckman imputation models for binary or continuous MNAR missing outcomes and MAR missing predictors (submitted).
}
\author{
Jacques-Emmanuel Galimard and Matthieu Resche-Rigon
}
\note{
\code{\link[miceMNAR:mice.impute.hecknorm]{mice.impute.hecknorm}()} should be preferred as \code{\link[miceMNAR:mice.impute.hecknorm2step]{mice.impute.hecknorm2step}()}. Accordingly, for continuous outcome, \code{\link[miceMNAR:mice.impute.hecknorm]{mice.impute.hecknorm}()} is a default. However, it is possible to modify the \code{method} argument by hand.
}

\seealso{
\code{\link{miceMNAR}}
}
\examples{

require("GJRM")
require(mvtnorm)
require(sampleSelection)

# Generation of a simulated dataset with MNAR mechanism on a continuous outcome

X1 <- rnorm(500,0,1)
X2 <- rbinom(500,1,0.5)
X3 <- rnorm(500,1,0.5)
  
errors <- rmvnorm(500,mean=c(0,0),sigma=matrix(c(1,0.3,0.3,1),nrow=2,byrow=TRUE))

Y <- X1+X2+errors[,1]
Ry <- ifelse(0.66+1*X1-0.5*X2+X3+errors[,2]>0,1,0)

Y[Ry==1] <- NA
  
simul_data <- data.frame(Y,X1,X2,X3)

JointModelEq <- generate_JointModelEq(data=simul_data,varMNAR = "Y")

JointModelEq[,"Y_var_sel"] <- c(0,1,1,1)
JointModelEq[,"Y_var_out"] <- c(0,1,1,0)

arg <- MNARargument(data=simul_data,varMNAR="Y",JointModelEq=JointModelEq)
arg$method["Y"] <- "hecknorm2step"

\dontrun{imputation3 <- mice(data = arg$data_mod,
                 method = arg$method,
                 predictorMatrix = arg$predictorMatrix,
                 JointModelEq=arg$JointModelEq,
                 control=arg$control,
                 maxit=1,m=5)

analysis3 <- with(imputation3,lm(Y~X1+X2+X3))
result3 <- pool(analysis3)
summary(result3)}
}
