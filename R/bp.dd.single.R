#'Creates control and manipulated Beta-Poisson samples for pre-specified parameters
#'
#'Creates control and manipulated Beta-Poisson (BP) samples for pre-specified BP parameters
#'
#'@import BPSC
#'@name bp.dd.single
#'@details Creates control and manipulated BP samples for pre-specified BP parameters
#'@usage bp.dd.single(pars,N1,N2,case,degree,seedex)
#'
#'@param pars	vector consisting of pre-specified values for the BP parameters alpha, beta, lambda1, lambda2 and p0 (in this order!)
#'@param N1	size of the sample drawn from the control BP model
#'@param N2	size of the sample drawn from the manipulated BP model
#'@param case	case corresponding to a specific parameter manipulation in the BP model; specifically, case can be one of "DLambda", "DAlpha", "DBeta", "DAlphaBeta" or "DPZ" following the design and the respective descriptive table in Schefzik (2020)
#'@param degree	parameter to set the degree of the created difference (low to strong), see Schefzik (2020) for details and the choice of a range of possible values
#'@param seedex	seed used for sampling from the fitted BP model to ensure reproducibility
#'
#'@return A list of two:
#'\itemize{
#'\item sim.bp.ctrl: vector of simulated control BP sample
#'\item sim.bp.manip: vector of simulated manipulated BP sample
#'}
#'
#'@references R. Schefzik (2020). Simulating differential distributions in Beta-Poisson models, in particular for single-cell RNA sequencing data.
#'
#'@examples
#'#create vector consisting of pre-specified values for the BP parameters alpha, beta, lambda1,
#'#lambda2 and p0 (in this order!)
#'pars<-c(0.21,3.50,96.13,0.05218,0.02)
#'N1<-500
#'N2<-500
#'degree<-1/3
#'seedex<-24
#'dd1<-bp.dd.single(pars,N1,N2,case="DLambda",degree,seedex)
#'dd2<-bp.dd.single(pars,N1,N2,case="DAlpha",degree,seedex)
#'dd3<-bp.dd.single(pars,N1,N2,case="DBeta",degree,seedex)
#'dd4<-bp.dd.single(pars,N1,N2,case="DAlphaBeta",degree,seedex)
#'dd5<-bp.dd.single(pars,N1,N2,case="DPZ",degree,seedex)
#'
#'
#'@export
#'
bp.dd.single<-function(pars,N1,N2,case,degree,seedex){

  ####create control
  set.seed(seedex)
  sim.bp.ctrl<-BPSC::rBP(n=N1,alp=pars[1],bet=pars[2],lam1=pars[3],lam2=pars[4],prob0=pars[5])

  names(sim.bp.ctrl)<-paste0("Sample",1:N1)

  ###expected values
  exp.val<-(1-pars[5])*pars[4]*pars[3]*(pars[1]/(pars[1]+pars[2]))


  ##variances
  var.bp4<-(pars[4]^2)*(pars[3]*(pars[1]/(pars[1]+pars[2]))+((pars[3]^2)*(pars[1]*pars[2]/(((pars[1]+pars[2])^2)*(pars[1]+pars[2]+1)))))
  exp.val.sq.bp4<-(pars[4]*pars[3]*(pars[1]/(pars[1]+pars[2])))^2
  variances<-(1-pars[5])*(exp.val.sq.bp4+var.bp4)-(exp.val^2)

##create manipulation
  if(case=="DLambda") {

    set.seed(seedex)

    Delta<-degree
    sim.bp.manip<-BPSC::rBP(n=N2,alp=pars[1],bet=pars[2],lam1=Delta*pars[3],lam2=pars[4],prob0=pars[5])

    names(sim.bp.manip)<-paste0("Sample",1:N2)

    lst<-list(sim.bp.ctrl,sim.bp.manip)
    names(lst)<-c("sample.ctrl","sample.manip")

    return(lst)
  }

  if(case=="DAlpha"){


    delta.ii<-function(theta){
      alpha<-pars[1]
      beta<-pars[2]
      del2<-(beta*theta)/(alpha+beta-alpha*theta)
      return(del2)
    }


    set.seed(seedex)

    Delta<-delta.ii(degree)
    sim.bp.manip<-BPSC::rBP(n=N2,alp=Delta*pars[1],bet=pars[2],lam1=pars[3],lam2=pars[4],prob0=pars[5])

    names(sim.bp.manip)<-paste0("Sample",1:N2)

    lst<-list(sim.bp.ctrl,sim.bp.manip)
    names(lst)<-c("sample.ctrl","sample.manip")

    return(lst)
  }


  if(case=="DBeta"){


    delta.iii<-function(theta){
      alpha<-pars[1]
      beta<-pars[2]
      del3<-(alpha+beta-alpha*theta)/(beta*theta)
      return(del3)
    }



    set.seed(seedex)

    Delta<-delta.iii(degree)
    sim.bp.manip<-BPSC::rBP(n=N2,alp=pars[1],bet=Delta*pars[2],lam1=pars[3],lam2=pars[4],prob0=pars[5])

    names(sim.bp.manip)<-paste0("Sample",1:N2)

    lst<-list(sim.bp.ctrl,sim.bp.manip)
    names(lst)<-c("sample.ctrl","sample.manip")

    return(lst)


  }

  if(case=="DAlphaBeta") {


    lambda2<-pars[4]
    exp.Y<-exp.val/(1-pars[5])
    lower.limits<-(1/variances)*(exp.val*(exp.Y+lambda2-exp.val))


    if(lower.limits<=degree){

      delta.iv<-function(theta){
        exp.Y<-exp.val/(1-pars[5])
        alpha<-pars[1]
        beta<-pars[2]
        lambda1<-pars[3]
        lambda2<-pars[4]
        p0<-pars[5]
        C.theta<-((1/(lambda2^2))*(((variances*theta+exp.val^2)/(1-p0))-exp.Y^2))-((lambda1*alpha)/(alpha+beta))
        del4<-(1/(alpha+beta))*((((lambda1^2)*alpha*beta)/(C.theta*((alpha+beta)^2)))-1)
        return(del4)
      }




      set.seed(seedex)

      Delta<-delta.iv(degree)
      sim.bp.manip<-BPSC::rBP(n=N2,alp=Delta*pars[1],bet=Delta*pars[2],lam1=pars[3],lam2=pars[4],prob0=pars[5])

      names(sim.bp.manip)<-paste0("Sample",1:N2)

      lst<-list(sim.bp.ctrl,sim.bp.manip)
      names(lst)<-c("sample.ctrl","sample.manip")

      return(lst)

    } else{
      lst<-list(sim.bp.ctrl,rep(NA,N2))
      names(lst)<-c("sample.ctrl","sample.manip")
      return(lst)
    }

  }

  if(case=="DPZ"){
    p0.manip<-function(p0,Delta){
      if(Delta<=(1-p0)){output<-(p0+Delta)}
      else {output<-(p0-Delta)}
      return(output)
    }



    set.seed(seedex)

    sim.bp.manip<-BPSC::rBP(n=N2,alp=pars[1],bet=pars[2],lam1=pars[3],lam2=pars[4],prob0=p0.manip(pars[5],degree))

    names(sim.bp.manip)<-paste0("Sample",1:N2)

    lst<-list(sim.bp.ctrl,sim.bp.manip)
    names(lst)<-c("sample.ctrl","sample.manip")

    return(lst)

  }

}


