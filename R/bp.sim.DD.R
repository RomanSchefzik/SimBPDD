#'Creates control and manipulated samples based on Beta-Poisson model fits for a single-cell RNA-sequencing data set
#'
#'Creates control samples based on Beta-Poisson (BP) model fits for a single-cell RNA-sequencing data set and manipulated samples based on those control samples for simulations involving BP models
#'
#'@import BPSC
#'@name bp.sim.DD
#'@details Creates control samples based on BP model fits for a single-cell RNA-sequencing data set and manipulated samples based on those control samples for simulations involving BP models. Details regarding the design of the manipulations can be found in Schefzik (2021). Combines the functions \code{bp.sim.control} and \code{bp.sim.manipulated}.
#'@usage bp.sim.DD(DATA,N1,N2,case,degree,seedex)
#'
#'@param DATA	matrix of single-cell RNA-sequencing expression data with genes in rows and samples (cells) in columns
#'@param N1	size of the samples drawn from the control BP models
#'@param N2	size of the samples drawn from the manipulated BP models
#'@param case	case corresponding to a specific parameter manipulation in the BP model; specifically, case can be one of "DLambda", "DAlpha", "DBeta", "DAlphaBeta" or "DPZ" following the design and the respective descriptive table in Schefzik (2021)
#'@param degree	parameter to set the degree of the created difference (low to strong), see Schefzik (2021) for details and the choice of a range of possible values
#'@param seedex	seed used for sampling from the fitted BP models to ensure reproducibility
#'
#'@return A list of two:
#'\itemize{
#'\item controls:	a list of two:
#'\itemize{
#'\item samples.bp.wellfit:	matrix with simulated control single-cell RNA-sequencing expression data of dimension GxN1, where G is the number of well-fitted genes by the BP models
#'\item parameters.bp.wellfit:	matrix with G rows containing the fitting results for the corresponding BP models in the columns, namely the parameters alpha, beta, lambda1, lambda2 and p0 as in Vu et al. (2016) or Schefzik (2021) and the Monte-Carlo-method-based p-value MCpval derived to check the validity of the BP fit (here, a fit is considered to be reasonably good if MCpval>0.05)
#'}
#'\item manipulations:	matrix with simulated manipulated single-cell RNA-sequencing expression data of dimension KxN2, where K is the number of well-fitted genes by the BP models for which the desired manipulation is feasible, see Schefzik (2021) for details
#'}
#'
#'@references R. Schefzik (2021). SimBPDD: Simulating differential distributions in Beta-Poisson models, in particular for single-cell RNA-sequencing data. Annales Mathematicae et Informaticae, 53:283-298. Available at \url{https://ami.uni-eszterhazy.hu/uploads/papers/finalpdf/AMI_53_from283to298.pdf} \cr
#' T. N. Vu, Q. F. Wills, K. R. Kalari, N. Niu, L. Wang, M. Rantalainen, and Y. Pawitan (2016). Beta-Poisson model for single-cell RNA seq data analyses. Bioinformatics, 32:2128-2135.
#'
#'
#'@examples
#'N1<-500
#'N2<-500
#'degree<-1/3
#'seedex<-24
#'bp1<-bp.sim.DD(DATA.EX,N1,N2,case="DLambda",degree,seedex)
#'bp2<-bp.sim.DD(DATA.EX,N1,N2,case="DAlpha",degree,seedex)
#'bp3<-bp.sim.DD(DATA.EX,N1,N2,case="DBeta",degree,seedex)
#'bp4<-bp.sim.DD(DATA.EX,N1,N2,case="DAlphaBeta",degree,seedex)
#'bp5<-bp.sim.DD(DATA.EX,N1,N2,case="DPZ",degree,seedex)
#'
#'
#'@export
#'
bp.sim.DD<-function(DATA,N1,N2,case,degree,seedex){

  ctrl<-bp.sim.control(DATA,N1,seedex)
  manip<-bp.sim.manipulated(ctrl$parameters.bp.wellfit,N2,case,degree,seedex)

  lst<-list(ctrl,manip)
  names(lst)<-c("controls","manipulations")

  return(lst)

}




