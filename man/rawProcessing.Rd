% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rawProcessing.R
\name{rawProcessing}
\alias{rawProcessing}
\title{rawProcessing}
\usage{
rawProcessing(
  rawdata,
  stratify = FALSE,
  frequency_threshold,
  data.squashing = FALSE,
  zeroes = FALSE,
  keep_pts
)
}
\arguments{
\item{rawdata}{A data frame containing columns named: primaryid, PT (Adverse Events), prod_ai (drug), strat1 (optional), strat2 (optional), ...... stratx}

\item{stratify}{A logical scalar specifiying whether stratification is to be used (TRUE) or not (FALSE)}

\item{frequency_threshold}{the minimum frequency of each AE or drug for it to be kept in the dataset}

\item{data.squashing}{whether to conduct data squashing or not (TRUE or FALSE)}

\item{zeroes}{A logical scalar specifying if zero counts should be included.}

\item{keep_pts}{A vector of whole numbers for the number of points to leave unsquashed for each count (N). See the 'openEBGM' details section.}
}
\value{
a list including the following:
\itemize{
\item \code{N_ij} actual counts of each AE - drug pair in matrix format
\item \code{Nij} actual counts of each AE - drug pair in vector format
\item \code{E_ij} expected counts of each AE - drug pair in matrix format
\item \code{Eij} expected counts of each AE - drug pair in vector format
\item \code{processedData} dataset after deleting AE or drug with frequencies less than \code{frequency_threshold}
\item \code{drugList} a vector of drug names that are kept in \code{processedData}
\item \code{AElist} a vector of adverse event names that are kept in \code{processedData}
\item \code{N_ij_squashed} a list of squashed N - E - weight dataset for each drug
\item \code{N_ij_list} a list of original N - E - weight (= 1) dataset for each drug
\item \code{squashed} squashed \code{processedData} dataset
}
}
\description{
This \code{rawProcessing} function finds the actual (N_ij) and expected counts (E_ij) of each AE - drug pair by implementing the methodology described by DuMouchel (1999); This function outputs N and E in both matrix format and vector format.
}
\details{
An \code{primaryid} column must be included. Columns must be named as as \code{primaryid, PT, prod_ai, strat1, strat2, ... stratx}. Only variables containing 'strat' (case sensitive) will be used as stratification variables.
}
\examples{
primaryid = c(1:10)
PT = c("fatigue", "restlessness", "sleepdisorder", "fatigue", "anaemia", "blood creatinine increased", "eosinophilia", "generalised oedema", "hypertension", "nephropathy toxic")
prod_ai = c("peginterferon alfa-2a", "peginterferon alfa-2a", "peginterferon alfa-2a", "pantoprazole", "ribavirin", "cyclosporine", "cyclosporine", "cyclosporine", "ribavirin", "flunitrazepam")
strat1 = c("M", "M", "M", "F", "F", "M", "M", "M", "F", "F")
strat2 = c("Young", "Young", "Young", "Old", "Old", "Old", "Old", "Young", "Young", "Young")
dat = data.frame(primaryid = primaryid, PT = PT, prod_ai = prod_ai, strat1 = strat1, strat2 = strat2)
result = hgzips::rawProcessing(dat, stratify = TRUE, frequency_threshold = 0, data.squashing = FALSE, zeroes = TRUE, keep_pts = 190000)
}
\references{
DuMouchel W (1999). "Bayesian Data Mining in Large Frequency Tables, With an Application to the FDA Spontaneous Reporting System." The American Statistician, 53(3), 177-190.
openEBGM R - package
}
\seealso{
openEBGM
}
