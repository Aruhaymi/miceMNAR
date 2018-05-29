\name{mice.impute.heckprob}
\alias{mice.impute.heckprob}

\title{
Imputation by bivariate probit sample selection model for binary Missing Not At Random outcome
}
\description{
Imputes univariate binary Missing Not At Random (MNAR) outcome using the bivariate probit sample selection model.
}
\usage{
mice.impute.heckprob(y, ry, x, JointModelEq, control, ...)
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
Imputation of binary MNAR outcome using the bivariate probit sample selection model. This imputation model requires defining a selection and an outcome equation for the bivariate probit model. The imputation algorithm consists of the following steps:
\enumerate{
\item Fit the bivariate probit model and compute the posterior distribution of parameters
\item Draw imputation parameters from their posterior distributions
\item Draw imputed values from their predictive distribution
}
}
\value{
A vector of length \code{nmis} with imputations (0 or 1).
}
%\references{
%Galimard, J.E., Chevret, S., Curis, E. and Resche-Rigon, M. (2018). Heckman imputation models for binary or continuous MNAR missing outcomes and MAR missing predictors (submitted).
%}
\author{
Jacques-Emmanuel Galimard and Matthieu Resche-Rigon
}

\seealso{
\code{\link{miceMNAR}}
}

\examples{

# Import dataset with a suspected MNAR mechanism
require("GJRM")
require(mvtnorm)
require(pbivnorm)
data("hiv") 

# We select only one region (lusuka) and 5 variables
lusuka <- hiv[hiv$region==5,c("hiv", "age", "marital", "condom", "smoke")]

# Categorical variables have to be recoded as factor
lusuka$hiv <- as.factor(lusuka$hiv)

# Specify a selection (missing data mechanism) and an outcome equation (analyse model)

# Generate an empty matrix

JointModelEq <- generate_JointModelEq(data=lusuka,varMNAR = "hiv")

# Fill in with 1 for variable included in equations
JointModelEq[,"hiv_var_sel"] <- c(0,1,1,1,1)
JointModelEq[,"hiv_var_out"] <- c(0,1,1,1,0)

# Generation of argument for MNAR imputation model in "mice()" function
arg <- MNARargument(data=lusuka,varMNAR="hiv",JointModelEq=JointModelEq)

# Imputation using mice() function
# Values returned have to be included in the "mice()" function as argument:

\dontrun{imputation <- mice(data = arg$data_mod,
                 method = arg$method,
                 predictorMatrix = arg$predictorMatrix,
                 JointModelEq=arg$JointModelEq,
                 control=arg$control,
                 maxit=1,m=5)

# Because of missing data only on one variable, fix maxit=1

# Estimation on each imputed dataset and pooling               
analysis <- with(imputation,glm(hiv~age+condom+marital,family=binomial(link="probit")))
result <- pool(analysis)
summary(result)}
}
