#' @title Empirical Calibration Density Generator

#' @description
#' This function uses the fossil record of a clade to generate an empirical probability density function that can be used as a calibration density or calibration prior in Bayesian divergence time estimation. Relevant quantiles such as the median or the 95% quantile can also obtained to be used in analyses based on fixed ages or bounds (e.g. penalized likelihood).

#' @param ages vector of fossil ages (if ages are known exactly) or or a matrix with fossils in rows and two columns: the first with the minimum age bounds (upper stratigraphic bounds) and the second with the maximum age bounds (lower stratigraphic bounds).
#' @param p vector of probabilities indicating the quantiles to be reported.
#' @param method a character string specifying the method: \code{"StraussSadler"} (default), \code{"RobsonWhitlock"}, \code{"Beta"}, \code{"NorrisPenGap"}, \code{"NorrisGhostLin"}, or \code{"OLE"}. See [pdate] for details.
#' @param KStest if TRUE, a Kolmogorov-Smirnov test is returned testing the null hypothesis that the distribution of fossil ages is uniform.
#' @param n number of Monte Carlo replicates for generating the empirical distribution.
#' @param PDFfitting a character string specifying a standard probability density function ( \code{"lognormnal"}, \code{"gamma"}, \code{"exponential"}) to be fit to the sample. "best" returns the best-fit density among those. "skewStudent" fits a skew-Student distribution using function \code{sstdFit} in the \code{fGarch} library. NULL prevents fitting.
#' @param plot If TRUE, plots fossil ages, a histogram of the Monte Carlo replicates, the median, and the 95% quantile (and the fitted probability density function if available).
#' @param repvalues	If TRUE, all n Monte Carlo replicates are returned.
#' @param breaks Specifies the algorithm used to determine the number of histogram cells in the plot (see [hist] for alternatives, if needed).
#' @param ... additional arguments passed to \code{qdate} (see [pdate] for details)

#' @details
#' This function generates an empirical distribution representing the probability density of the age of a clade of organisms. The distribution is generated using a Monte Carlo algorithm based and methods for estimating the upper bound of a sample of fossils ages (see [pdate]). The algorithm incorporates fossil age uncertainty (if present) by sampling ages from the time interval were the fossil was found. A Kolmogorov-Smirnov test can be requested to test the null hypothesis that the distribution of fossil ages is uniform (using function [ks.test]); the bounds under the null are estimated using the StraussSadler method and if fossil ages are interval, a random sample of ages from those intervals is used in the test.

#' @return a list with the following elements:
#' \itemize{
#' \item{ages: fossil ages}
#' \item{method: the estimation method used}
#' \item{rep.values: all Monte Carlo replicates (if requested)}
#' \item{Quantiles}
#' \item{KStest: result of the Kolmogorov-Smirnov test of uniformity (if requested), as returned by [ks.test]}
#' \item{PDFfit.model: name of probability density model fit (if requested)}
#' \item{PDFfit: list including the results of the probability density fit, including estimated parameters (if requested)}
#' }
#' @import stats grDevices graphics fGarch

#' @examples
#'
#' ## Matrix of fossil ages including maximum and minimum bounds for each
#' Fages <- cbind(c(50, 38, 30, 14, 3.5), c(56, 45, 35, 14, 6))
#'
#' ## Basic execution
#'   clade.date(Fages)
#'
#' ## Alternative options
#'   clade.date(Fages, method="OLE", plot=TRUE, PDFfitting="lognormal")
#'
#' @author Santiago Claramunt, \email{sclaramunt@@rom.on.ca}

#' @references

#' Claramunt, S. & J. L. Cracraft. (2015) A new time tree reveals Earth history???s imprint on the evolution of modern birds. *Science Advances* **1**:e1501005 <https://doi.org/10.1126/sciadv.1501005>
#'
#'
#' @seealso [pdate] for point estimates when fossil ages are known exactly, quantile functions, and the random generator function.

#' @export

clade.date <- function(ages, p=c(0, 0.5, 0.95), n=10000, method="StraussSadler", KStest=FALSE, PDFfitting="best", repvalues=TRUE, plot=FALSE, ...) {
	
	ages <- as.matrix(ages)
		
	# Check format of ages
	if(ncol(ages) > 2) stop("'ages' must be a vector or a two-column matrix")

	# Check that bounds are correct for all fossils
	if(any(apply(ages, 1, diff) < 0)) {
		stop("Detected at least one fossil with older bound < younger bound") }

	# Create list for results
	RES <-list()

	RES$ages <- ages

	RES$method <- method
	
	# Obtain replciates
	rA <- rdate(ages=ages, n=n, method=method, ...)
	
	if(repvalues) { RES$rep.values <- rA }

	# Compute quantiles
	Qs <- quantile(rA, prob=p)	
			
	RES$Quantiles <- Qs
	
	# Test of Uniformity
	if(KStest && nrow(ages) > 2) {
		
		# If fossil ages are interval, take one sample from intervals
		if(ncol(ages) == 2) {
		Mages <- runif(nrow(ages), min=ages[,1], max=ages[,2])
		} else Mages <- ages	
		
		# Estimate the bounds of the distribution under the null using StraussSadler method
		range <- max(Mages) - min(Mages)
		
		N <- length(Mages)
		
		Min <- min(Mages) - range/N
		
		Max <- max(Mages) + range/N
												
		KStest <- ks.test(Mages, 'punif', min=max(c(0,Min)), max=Max)
						
		if(KStest$p.value < 0.05 && method=="StraussSadler") { warning("The data may not be uniformly distributed (Kolmogorov-Smirnov P = ", signif(KStest$p.value, digits=2), ")", call.=FALSE) }
		
		RES$KStest <- KStest
			}


	# Fit a probability distribution to the MC sample, if requested
				
	if(!is.null(PDFfitting)) {

		rA <- rA[is.finite(rA)]

		if(PDFfitting=="best") {
			
			model.names <- c("lognormal", "gamma", "exponential")
			
			PDFfit.lognormal <- MASS::fitdistr(x=rA-max(ages[,1]), densfun="lognormal")
			PDFfit.gamma <- MASS::fitdistr(x=rA-max(ages[,1]), densfun="gamma")
			PDFfit.exponential <- MASS::fitdistr(x=rA-max(ages[,1]), densfun="exponential")
			
			AICs <- c(AIC(PDFfit.lognormal), AIC(PDFfit.gamma), AIC(PDFfit.exponential))
			
			best.model <- which.min(AICs)
			
			# Assign the best model name to PDFfitting so it is plotted and reported
			
			PDFfitting <- model.names[best.model]
			
			RES$PDFfit.model <- PDFfitting
			
			RES$PDFfit <- switch(best.model, PDFfit.lognormal, PDFfit.gamma, PDFfit.exponential)

			} else if(PDFfitting=="skewStudent") {		
		
			RES$PDFfit.model <- PDFfitting

			PDFfit <- fGarch::sstdFit(rA)
				
			RES$PDFfit <- PDFfit
			} else if(PDFfitting=="skewnormal") {		
		
			PDFfit <- fGarch::snormFit(rA)
		
			RES$PDFfit.model <- PDFfitting

			RES$PDFfit <- PDFfit

			} else {		
			
			PDFfit <- MASS::fitdistr(x=rA-max(ages[,1]), densfun=PDFfitting)

			RES$PDFfit.model <- PDFfitting
				
			RES$PDFfit <- PDFfit
			}
						
		}
		
	class(RES) <- "clade.date"

	if(plot) { plot(RES)}
	
	return(RES) 
}	




