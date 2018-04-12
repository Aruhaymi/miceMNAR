\name{generate_JointModelEq}
\alias{generate_JointModelEq}

\title{
Generation of an empty matrix for sample selection model
}
\description{
This function generates a matrix of 0 with rows corresponding to the names of variables included in the dataset and column corresponding to selection and outcome equations for each MNAR outcome. The user has to fill this matrix with 1 for variable included in selection and/or outcome equations of MNAR outcome imputation models.
}
\usage{
generate_JointModelEq(varMNAR, data)
}

\arguments{
  \item{varMNAR}{
The name of MNAR outcome to be imputed.
}
  \item{data}{
The dataset used for classical multiple imputation by chained equation and additional variables necessary for MNAR imputation models.
}
}
\details{
Be careful to not define the same selection and outcome equations for MNAR imputation models. A constraint of the sample selection model implies the inclusion of different sets of covariates, which may or not be nested in the selection equation and the outcome equation, to avoid collinearity issues. It has been recommended to include at least a supplementary variable in the selection equation. This variable should be known to be unlinked directly to the outcome.
}
\value{
A matrix such as:
\itemize{ 
\item Rows correspond to names of variables included in the \code{data} argument
\item Columns correspond to selection and outcome equations for each MNAR outcome imputation models
}
Each column is named as "MNAR outcome" names adding "_var_sel" or "_var_out" corresponding to variables included in selection and outcome equation respectively.
}

\author{
Jacques-Emmanuel Galimard
}

\seealso{
\code{\link{miceMNAR}}
\code{\link[GJRM:hiv]{hiv}}
}
\examples{
# Import dataset with a suspected MNAR mechanism
require("GJRM")
data("hiv") 

# We select only one region (lusuka) and 5 variables
lusuka <- hiv[hiv$region==5,c("hiv", "age", "highhiv", "condom", "smoke")]

# Categorical variables have to be recoded as factor
lusuka$hiv <- as.factor(lusuka$hiv)

# Specify a selection (missing data mechanism) and an outcome equation (analyse model)

# Generate an empty matrix
JointModelEq <- generate_JointModelEq(data=lusuka,varMNAR = "hiv")

# Fill in with 1 for variable included in equations
JointModelEq[,"hiv_var_sel"] <- c(0,1,1,1,1) 
# This indicates that age, highhiv, condom and smoke are included in the selection equation of hiv
JointModelEq[,"hiv_var_out"] <- c(0,1,1,1,0) 
# This indicates that age, highhiv and condom are included in the outcome equation of hiv

}

