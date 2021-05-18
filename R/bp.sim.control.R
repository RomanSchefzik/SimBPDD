#'Creates control samples for simulations involving Beta-Poisson models
#'
#'Creates control samples for simulations involving Beta-Poisson (BP) models
#'
#'@import BPSC
#'@name bp.sim.control
#'@details Creates control samples for simulations involving BP models. For a given single-cell RNA-sequencing data set, a five-parameter BP model (Vu et al., 2016) is fitted to each gene. Then, it is checked whether the BP model is actually a good fit following the procedure in Vu et al. (2016), and instances with a good fit are kept. For each of those well-fitted genes, a sample of size N1 is drawn from the respective fitted BP model.
#'@usage bp.sim.control(DATA,N1,seedex)
#'
#'@param DATA	matrix of single-cell RNA-sequencing expression data with genes in rows and samples (cells) in columns
#'@param N1	size of the samples drawn from the control BP models
#'@param seedex	seed used for sampling from the fitted BP models to ensure reproducibility
#'
#'@return A list of two:
#'\itemize{
#'\item  samples.bp.wellfit:	matrix with simulated control single-cell RNA-sequencing expression data of dimension GxN1, where G is the number of well-fitted genes by the BP models
#'\item parameters.bp.wellfit:	matrix with G rows containing the fitting results for the corresponding BP models in the columns, namely the parameters alpha, beta, lambda1, lambda2 and p0 as in Vu et al. (2016) or Schefzik (2021) and the Monte-Carlo-method-based p-value MCpval derived to check the validity of the BP fit (here, a fit is considered to be reasonably good if MCpval>0.05)
#'}
#'
#'@references R. Schefzik (2021). SimBPDD: Simulating differential distributions in Beta-Poisson models, in particular for single-cell RNA-sequencing data. Annales Mathematicae et Informaticae, 53:283-298. Available at \url{https://ami.uni-eszterhazy.hu/uploads/papers/finalpdf/AMI_53_from283to298.pdf}\cr
#' T. N. Vu, Q. F. Wills, K. R. Kalari, N. Niu, L. Wang, M. Rantalainen, and Y. Pawitan (2016). Beta-Poisson model for single-cell RNA seq data analyses. Bioinformatics, 32:2128-2135.
#'
#'
#'@examples
#'N1<-500
#'seedex<-24
#'ctrl<-bp.sim.control(DATA.EX,N1,seedex)
#'
#'
#'@export
#'
bp.sim.control<-function(DATA,N1,seedex){

###1. for a given scRNA-seq data matrix (rows: genes; columns:cells), fit the corresponding BP models
set.seed(seedex)


N.simMC<-100

Res.par1<-array(NA,dim=c(dim(DATA)[1],5))

MCpval<-rep(dim(DATA)[1],10)

for (g in 1:dim(DATA)[1]){
  Res<-BPSC::estimateBP(DATA[g,],para.num=5)
  Res.par1[g,]<-Res$par

  #Generate Monte-Carlo null distribution of the model.
  #Due to time limit, the number of simulations (sim.num) here is set by 100.
  MCnull.Res<-BPSC::getBPMCnull(Res$par,n=100,tbreak=Res$tbreak,sim.num=N.simMC)

  #Compute Monte-Carlo p-value
  MCpval[g]<-sum(MCnull.Res$X2 >= Res$X2)/length(MCnull.Res$X2)
}


Res.par<-cbind(Res.par1,MCpval)


##indices for genes for which model fits well (here based on a 5% threshold)
ind.wellfit<-which(Res.par[,6]>0.05)


##2. simulate corresponding data for genes that have been fitted well by BP model

sim.bp.wellfit<-array(data=NA,dim=c(length(ind.wellfit),N1))



for (i in 1:length(ind.wellfit)){
  set.seed(seedex)
  sim.bp.wellfit[i,]<-BPSC::rBP(n=N1,alp=Res.par[ind.wellfit[i],1],bet=Res.par[ind.wellfit[i],2],lam1=Res.par[ind.wellfit[i],3],lam2=Res.par[ind.wellfit[i],4],prob0=Res.par[ind.wellfit[i],5])
}


rownames(sim.bp.wellfit)<-paste0("Instance",1:dim(sim.bp.wellfit)[1])
colnames(sim.bp.wellfit)<-paste0("Sample",1:dim(sim.bp.wellfit)[2])

parameters.wellfit<-Res.par[ind.wellfit,]
rownames(parameters.wellfit)<-paste0("Instance",1:dim(parameters.wellfit)[1])
colnames(parameters.wellfit)<-c("alpha","beta","lambda1","lambda2","p0","MCpval")


lst<-list(sim.bp.wellfit,parameters.wellfit)
names(lst)<-c("samples.bp.wellfit","parameters.bp.wellfit")

return(lst)

}




