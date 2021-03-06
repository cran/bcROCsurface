% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bcvusMAR.R
\name{vus}
\alias{vus}
\title{Estimation methods for volume under ROC surface (VUS)}
\usage{
vus(
  method = "full",
  T,
  Dvec,
  V,
  rhoEst = NULL,
  piEst = NULL,
  ci = TRUE,
  ci.level = ifelse(ci, 0.95, NULL),
  BOOT = FALSE,
  nR = ifelse(ci, 250, NULL),
  parallel = FALSE,
  ncpus = ifelse(parallel, detectCores()/2, NULL),
  trace = TRUE
)
}
\arguments{
\item{method}{name of bias-corrected estimation method to be used for estimating the VUS in presence of verification bias. See \code{\link{ROCs}} for more details.}

\item{T}{a numeric vector containing the diagnostic test values. \code{NA} values are not admitted.}

\item{Dvec}{a n * 3  binary matrix with the three columns, corresponding to three classes of the disease status. In row i, 1 in column j indicates that the i-th subject belongs to class j, with j = 1, 2, 3. A row of \code{NA} values indicates a non-verified subject.}

\item{V}{a binary vector containing the verification status (1 verified, 0 not verified).}

\item{rhoEst}{a result of a call to \code{\link{rhoMLogit}} of \code{\link{rhoKNN}} to fit the disease model.}

\item{piEst}{a result of a call to \code{\link{psglm}} to fit the verification model.}

\item{ci}{a logical value. If TRUE (default), computes an confidence interval of VUS and tests the null hypothesis H0: VUS = 1/6.}

\item{ci.level}{an confidence level to be used for constructing the confidence interval; default 0.95.}

\item{BOOT}{a logical value. Default = \code{FALSE}. If set to \code{TRUE}, a bootstrap resampling is employed to estimate the asymptotic variance of bias-corrected VUS estimates. See \code{\link{asyVarVUS}}.}

\item{nR}{the number of bootstrap replicates, which is used for FULL or KNN estimator, or option \code{BOOT = TRUE}. Usually this will be a single positive integer.}

\item{parallel}{a logical value. If \code{TRUE}, a parallel computing is employed to the bootstrap resampling process.}

\item{ncpus}{number of processes to be used in parallel computing. Default is a half of available cores.}

\item{trace}{a logical value. If \code{TRUE}, tracing information on the progress of the estimation is produced.}
}
\value{
\code{vus} returns an object of class inheriting from "vus" class.

The function \code{\link{print.vus}} can be used to print a summary of the results.

An object of class "vus" is a list containing at least the following components:

\item{vus.fit}{the estimate of VUS.}
\item{std}{the standard error, obtained by using asymptotic theory or bootstrap resampling method.}
\item{call}{the matched call.}
\item{t.stat}{t-statistic.}
\item{p.val_norm}{p-value correspond to normal-test.}
\item{ci.norm}{the confidence interval of VUS by using normal approximation.}
\item{ci.logit}{the confidence interval of VUS via logit transform.}
\item{ci.level}{the confidence level used.}
\item{BOOT}{the value of \code{BOOT}.}
\item{nR}{the number of bootstrap replicates used.}

In addition, the name of method used to estimate VUS also is given as the attribute of \code{vus.fit}.
}
\description{
\code{vus} computes bias-corrected estimates of the volume under the ROC surface for evaluating the accuracy of a continuous diagnostic test.
}
\details{
The function implements five bias-corrected estimation methods in To Duc et al. (2016, 2018) for estimating VUS of a three-class continuous diagnostic test in presence of verification bias. The estimators are full imputation (FI), mean score imputation (MSI), inverse probability weighted (IPW), semiparametric efficient (SPE) and K nearest-neighbor (KNN), see \code{\link{ROCs}}. These esitmators work under MAR assumption.

The standard error of the estimates are obtained through the function \code{\link{asyVarVUS}}. In particular, the standard error of the FULL estimate is computed by bootrap resampling method or by Jackknife approach proposed in Guangming et al. (2013). For the bias-corrected estimates, the standard errors are computed by using asymptotic theory (with respect to FI, MSI, IPW and SPE estimator) or bootstrap resampling method (with respect to KNN estimator). A confidence interval of VUS also is given. A logit transformation is also applied for obtaining the confidence interval.

The default value of the number of bootstrap replicates is 250.

Note that, before apply the functions \code{vus}, the use of \code{\link{preDATA}} might be needed to check the monotone ordering disease classes and to create the matrix format for disease status.
}
\examples{
data(EOC)
head(EOC)


\dontrun{
# FULL data estimator
Dfull <- preDATA(EOC$D.full, EOC$CA125)
Dvec.full <- Dfull$Dvec
vus("full", T = EOC$CA125, Dvec = Dvec.full)
}

# Preparing the missing disease status
Dna <- preDATA(EOC$D, EOC$CA125)
Dfact.na <- Dna$D
Dvec.na <- Dna$Dvec
# FI estimator
rho.out <- rhoMLogit(Dfact.na ~ CA125 + CA153 + Age, data = EOC, test = TRUE)

vus("fi", T = EOC$CA125, Dvec = Dvec.na, V = EOC$V, rhoEst = rho.out)

\dontrun{
# MSI estimator
vus("msi", T = EOC$CA125, Dvec = Dvec.na, V = EOC$V, rhoEst = rho.out)

# IPW estimator
pi.out <- psglm(V ~ CA125 + CA153 + Age, data = EOC, test = TRUE)
vus("ipw", T = EOC$CA125, Dvec = Dvec.na, V = EOC$V, piEst = pi.out)

# SPE estimator
vus("spe", T = EOC$CA125, Dvec = Dvec.na, V = EOC$V, rhoEst = rho.out, piEst = pi.out)

# KNN estimator, K = 1, Mahalanobis distance
XX <- cbind(EOC$CA125, EOC$CA153, EOC$Age)
rho.maha.1nn <- rhoKNN(X = XX, Dvec = Dvec.na, V = EOC$V, K = 1, type = "mahala")
vus("knn", T = EOC$CA125, Dvec = Dvec.na, V = EOC$V, rhoEst = rho.maha.1nn)
}

}
\references{
To Duc, K., Chiogna, M. and Adimari, G. (2018)
Nonparametric estimation of ROC surfaces in presence of verification bias.
\emph{REVSTAT Statistical Journal}. Accepted.

To Duc, K., Chiogna, M. and Adimari, G. (2016)
Bias-corrected methods for estimating the receiver operating characteristic surface of continuous diagnostic tests.
\emph{Electronic Journal of Statistics}, \bold{10}, 3063-3113.

Guangming, P., Xiping, W. and Wang, Z. (2013)
Non-parameteric statistical inference for $P(X < Y < Z)$.
\emph{Sankhya A}, \bold{75}, 1, 118-138.
}
