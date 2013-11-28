callCalibrationPerson <- function(seed=rep(12345,6),n=500,runpar=c(4,0.5,0.05,10,3,0.5)) {

  state <- RNGstate(); on.exit(state$reset())
  RNGkind("user")
  set.user.Random.seed(seed)
	
  out <- .Call("callCalibrationSimulation",as.integer(seed),list(n=as.integer(n),runpar=as.double(runpar)),PACKAGE="microsimulation")
  states <- c("DiseaseFree","Precursor","PreClinical","Clinical")
  curnames<-out$Key
#  curnames <- names(out)
  mat <- matrix(0,nr=10,ncol=length(states))
  colnames(mat) <- states; rownames(mat) <- seq(10,100,10)
  mat[,]<-data.matrix(as.data.frame(out$Value[match(states,out$Key)]))
	
  list(StateOccupancy=mat, TimeAtRisk=out$Value[out$Key=="TimeAtRisk"][[1]])
}
