#'Creates manipulated samples for simulations involving Beta-Poisson models
#'
#'Creates manipulated samples based on the control samples for simulations involving Beta-Poisson (BP) models
#'
#'@import BPSC
#'@name bp.sim.manipulated
#'@details Creates manipulated samples based on the control samples for simulations involving BP models. Details regarding the design of the manipulations can be found in Schefzik (2021).
#'@usage bp.sim.manipulated(Res.par,N2,case,degree,seedex)
#'
#'@param Res.par	matrix including the fitted parameters for the well-fitted control instances that are to be manipulated (e.g. the output matrix parameters.bp.wellfit from the \code{bp.sim.control} function)
#'@param N2	size of the samples drawn from the manipulated BP models
#'@param case	case corresponding to a specific parameter manipulation in the BP model; specifically, case can be one of "DLambda", "DAlpha", "DBeta", "DAlphaBeta" or "DPZ" following the design and the respective descriptive table in Schefzik (2021)
#'@param degree	parameter to set the degree of the created difference (low to strong), see Schefzik (2021) for details and the choice of a range of possible values
#'@param seedex	seed used for sampling from the fitted BP models to ensure reproducibility
#'
#'@return matrix with simulated manipulated single-cell RNA-sequencing expression data of dimension KxN2, where K is the number of well-fitted genes by the BP models for which the desired manipulation is feasible, see Schefzik (2021) for details
#'
#'@references R. Schefzik (2021). SimBPDD: Simulating differential distributions in Beta-Poisson models, in particular for single-cell RNA-sequencing data. Annales Mathematicae et Informaticae, 53:283-298. Available at \url{https://ami.uni-eszterhazy.hu/uploads/papers/finalpdf/AMI_53_from283to298.pdf}
#'
#'@examples
#'N1<-500
#'seedex<-24
#'ctrl<-bp.sim.control(DATA.EX,N1,seedex)
#'N2<-500
#'degree<-1/3
#'man1<-bp.sim.manipulated(ctrl[[2]],N2,case="DLambda",degree,seedex)
#'man2<-bp.sim.manipulated(ctrl[[2]],N2,case="DAlpha",degree,seedex)
#'man3<-bp.sim.manipulated(ctrl[[2]],N2,case="DBeta",degree,seedex)
#'man4<-bp.sim.manipulated(ctrl[[2]],N2,case="DAlphaBeta",degree,seedex)
#'man5<-bp.sim.manipulated(ctrl[[2]],N2,case="DPZ",degree,seedex)
#'
#'
#'@export
#'
bp.sim.manipulated<-function(Res.par,N2,case,degree,seedex){
  
  ###expected values
  exp.val<-rep(NA,dim(Res.par)[1])
  for (i in 1:dim(Res.par)[1]){
    exp.val[i]<-(1-Res.par[i,5])*Res.par[i,4]*Res.par[i,3]*(Res.par[i,1]/(Res.par[i,1]+Res.par[i,2]))
  }
  
  ##variances
  variances<-rep(NA,dim(Res.par)[1])
  for (i in 1:dim(Res.par)[1]){
    var.bp4<-(Res.par[i,4]^2)*(Res.par[i,3]*(Res.par[i,1]/(Res.par[i,1]+Res.par[i,2]))+((Res.par[i,3]^2)*(Res.par[i,1]*Res.par[i,2]/(((Res.par[i,1]+Res.par[i,2])^2)*(Res.par[i,1]+Res.par[i,2]+1)))))
    exp.val.sq.bp4<-(Res.par[i,4]*Res.par[i,3]*(Res.par[i,1]/(Res.par[i,1]+Res.par[i,2])))^2
    variances[i]<-(1-Res.par[i,5])*(exp.val.sq.bp4+var.bp4)-(exp.val[i]^2)
  }
  
  
  if(case=="DLambda") {
    
    
    sim.bp.wellfit.dd<-array(data=NA,dim=c(dim(Res.par)[1],N2))
    
    for (i in 1:dim(Res.par)[1]){
      Delta<-degree
      set.seed(seedex)
      sim.bp.wellfit.dd[i,]<-BPSC::rBP(n=N2,alp=Res.par[i,1],bet=Res.par[i,2],lam1=Delta*Res.par[i,3],lam2=Res.par[i,4],prob0=Res.par[i,5])
    }
    
    rownames(sim.bp.wellfit.dd)<-paste0("Instance",1:dim(sim.bp.wellfit.dd)[1])
    colnames(sim.bp.wellfit.dd)<-paste0("Sample",1:dim(sim.bp.wellfit.dd)[2])
    
    return(sim.bp.wellfit.dd)
  }
  
  if(case=="DAlpha"){
    
    
    delta.ii<-function(g,theta){
      alpha<-Res.par[g,1]
      beta<-Res.par[g,2]
      del2<-(beta*theta)/(alpha+beta-alpha*theta)
      return(del2)
    }
    
    sim.bp.wellfit.dd<-array(data=NA,dim=c(dim(Res.par)[1],N2))
    
    for (i in 1:dim(Res.par)[1]){
      Delta<-delta.ii(i,degree)
      set.seed(seedex)
      sim.bp.wellfit.dd[i,]<-BPSC::rBP(n=N2,alp=Delta*Res.par[i,1],bet=Res.par[i,2],lam1=Res.par[i,3],lam2=Res.par[i,4],prob0=Res.par[i,5])
    }
    
    rownames(sim.bp.wellfit.dd)<-paste0("Instance",1:dim(sim.bp.wellfit.dd)[1])
    colnames(sim.bp.wellfit.dd)<-paste0("Sample",1:dim(sim.bp.wellfit.dd)[2])
    
    return(sim.bp.wellfit.dd)
  }
  
  
  if(case=="DBeta"){
    
    
    delta.iii<-function(g,theta){
      alpha<-Res.par[g,1]
      beta<-Res.par[g,2]
      del3<-(alpha+beta-alpha*theta)/(beta*theta)
      return(del3)
    }
    
    
    sim.bp.wellfit.dd<-array(data=NA,dim=c(dim(Res.par)[1],N2))
    
    for (i in 1:dim(Res.par)[1]){
      Delta<-delta.iii(i,degree)
      set.seed(seedex)
      sim.bp.wellfit.dd[i,]<-BPSC::rBP(n=N2,alp=Res.par[i,1],bet=Delta*Res.par[i,2],lam1=Res.par[i,3],lam2=Res.par[i,4],prob0=Res.par[i,5])
    }
    
    rownames(sim.bp.wellfit.dd)<-paste0("Instance",1:dim(sim.bp.wellfit.dd)[1])
    colnames(sim.bp.wellfit.dd)<-paste0("Sample",1:dim(sim.bp.wellfit.dd)[2])
    
    return(sim.bp.wellfit.dd)
    
    
  }
  
  if(case=="DAlphaBeta") {
    
    
    
    ###calculate lower limits for theta and do some quality checks
    low.limit<-function(g){
      lambda2<-Res.par[g,4]
      exp.Y<-exp.val[g]/(1-Res.par[g,5])
      lower.limit<-(1/variances[g])*(exp.val[g]*(exp.Y+lambda2-exp.val[g]))
      return(lower.limit)
    }
    
    
    lower.limits<-rep(NA,dim(Res.par)[1])
    for (n in 1:dim(Res.par)[1]){
      lower.limits[n]<-low.limit(n)
    }
    
    
    
    ###how many cases are there with lower.limit<=degree ? (such that the degrees of DD could be constructed)
    #idx.sub.iv<-which(lower.limits<=degree)
    
    
    ##create a subset such that all previous degrees of DD are possible
    # Res.par<-Res.par[idx.sub.iv,]
    
    delta.iv<-function(g,theta){
      ##different degrees of DD
      #theta.values<-c(10/11,2/3,1/2,2/5,1/3)
      exp.Y<-exp.val[g]/(1-Res.par[g,5])
      alpha<-Res.par[g,1]
      beta<-Res.par[g,2]
      lambda1<-Res.par[g,3]
      lambda2<-Res.par[g,4]
      p0<-Res.par[g,5]
      C.theta<-((1/(lambda2^2))*(((variances[g]*theta+exp.val[g]^2)/(1-p0))-exp.Y^2))-((lambda1*alpha)/(alpha+beta))
      del4<-(1/(alpha+beta))*((((lambda1^2)*alpha*beta)/(C.theta*((alpha+beta)^2)))-1)
      return(del4)
    }
    
    
    
    sim.bp.wellfit.dd<-array(data=NA,dim=c(dim(Res.par)[1],N2))
    
    for (i in 1:dim(Res.par)[1]){
      
      if(lower.limits[i]<=degree){
        Delta<-delta.iv(i,degree)
        set.seed(seedex)
        sim.bp.wellfit.dd[i,]<-BPSC::rBP(n=N2,alp=Delta*Res.par[i,1],bet=Delta*Res.par[i,2],lam1=Res.par[i,3],lam2=Res.par[i,4],prob0=Res.par[i,5])
      }else{
        sim.bp.wellfit.dd[i,]<-rep(NA,N2)
      }
      
    }
    rownames(sim.bp.wellfit.dd)<-paste0("Instance",1:dim(sim.bp.wellfit.dd)[1])
    colnames(sim.bp.wellfit.dd)<-paste0("Sample",1:dim(sim.bp.wellfit.dd)[2])
    
    return(sim.bp.wellfit.dd)
    
  }
  
  if(case=="DPZ"){
    p0.manip<-function(p0,Delta){
      if(Delta<=(1-p0)){output<-(p0+Delta)}
      else {output<-(p0-Delta)}
      return(output)
    }
    
    
    sim.bp.wellfit.dd<-array(data=NA,dim=c(dim(Res.par)[1],N2))
    
    for (i in 1:dim(Res.par)[1]){
      set.seed(seedex)
      sim.bp.wellfit.dd[i,]<-BPSC::rBP(n=N2,alp=Res.par[i,1],bet=Res.par[i,2],lam1=Res.par[i,3],lam2=Res.par[i,4],prob0=p0.manip(Res.par[i,5],degree))
    }
    
    rownames(sim.bp.wellfit.dd)<-paste0("Instance",1:dim(sim.bp.wellfit.dd)[1])
    colnames(sim.bp.wellfit.dd)<-paste0("Sample",1:dim(sim.bp.wellfit.dd)[2])
    
    return(sim.bp.wellfit.dd)
    
  }
  
}

