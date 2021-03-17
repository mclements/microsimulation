#' call CalibrationPerson example
#'
#' @param seed random number seed
#' @param n number of simulations
#' @param runpar parameters
#' @param mc.cores number of cores
#' @return data-frame
#' @export
#' @rdname Examples
callCalibrationPerson <- function(seed=12345,n=500,runpar=c(4,0.5,0.05,10,3,0.5),mc.cores=1) {
  
  state <- RNGstate(); on.exit(state$reset())
  RNGkind("user")
  set.user.Random.seed(seed)
  
  f <- as.factor(rep(1:mc.cores,length.out=n))
  chunks <- split(1:n,f)
  initialSeeds <- list()
  currentSeed <- user.Random.seed()
  
  fun <- function(obj, i) if(i==1)obj else parallel::nextRNGStream(obj)
  for(i in 1:mc.cores){
    initialSeeds[[i]] <- fun(currentSeed,i)
    currentSeed <- initialSeeds[[i]]
  }
  out <- parallel::mclapply(1:mc.cores, function(i) {
    chunk <- chunks[[i]]
    set.user.Random.seed(initialSeeds[[i]])
    .Call("callCalibrationSimulation",list(n=as.integer(length(chunk)),runpar=as.double(runpar)),PACKAGE="microsimulation")
  })
  states <- c("DiseaseFree","Precursor","PreClinical","Clinical")
  out <- lapply(out, function(o){
    curnames <- names(o)
    mat <- matrix(0,nrow=10,ncol=length(states))
    colnames(mat) <- states; rownames(mat) <- seq(10,100,10)
    mat[,states[states %in% curnames]]<-data.matrix(transform(as.data.frame(o[states[states %in% curnames]])))
    list(StateOccupancy=mat, TimeAtRisk=o$TimeAtRisk)
  })
  Reduce(function(u,z) list(StateOccupancy=u$StateOccupancy + z$StateOccupancy, TimeAtRisk = u$TimeAtRisk + z$TimeAtRisk),out)
}
