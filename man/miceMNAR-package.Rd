\name{miceMNAR-package}
\alias{miceMNAR-package}
\alias{miceMNAR}
\docType{package}
\title{
\packageTitle{miceMNAR}
}
\description{
Provides imputation models and functions for binary or continuous Missing Not At Random (MNAR) outcomes through the use of the '\code{\link{mice}}' package. The \code{\link[miceMNAR:mice.impute.hecknorm]{mice.impute.hecknorm}()} function provides imputation model for continuous outcome based on Heckman's model also named sample selection model as described in Galimard et al (2016) <doi:10.1002/sim.6902>. The \code{\link[miceMNAR:mice.impute.heckprob]{mice.impute.heckprob}()} function provides imputation model for binary outcome based on bivariate probit model. 

As these two previous imputation models require to specify a selection and an outcome equation, \code{\link[mice:mice]{mice}()} function has to be adapted using \code{\link[miceMNAR:MNARargument]{MNARargument}()}.
}
\details{

\packageIndices{miceMNAR}
}
\author{
\packageAuthor{miceMNAR}

Maintainer: \packageMaintainer{miceMNAR}
}
\references{
%Galimard, J.E., Chevret, S., Curis, E. and Resche-Rigon, M. (2018). Heckman imputation models for binary or continuous MNAR missing outcomes and MAR missing predictors (submitted).
%
Galimard, J.-E., Chevret, S., Protopopescu, C., and Resche-Rigon, M. (2016) A multiple imputation approach for MNAR mechanisms compatible with Heckman's model. Statist. Med., 35: 2907-2920. doi:10.1002/sim.6902.
}

\section{Warning }{
This package is only validated for the imputation of MNAR outcome. However, it is implemented to impute several MNAR variables in the same process. Such implementation must be realised carefully.
}

\seealso{
\code{\link{mice}}
\code{\link[GJRM:copulaSampleSel]{copulaSampleSel}}
\code{\link[GJRM:copulaSampleSel]{SemiParBIV}}

\code{\link[sampleSelection:selection]{selection}}
}

