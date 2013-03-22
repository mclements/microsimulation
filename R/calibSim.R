callCalibrationPerson <- function(seed=rep(12345,6),n=500,runpar=c(4,0.5,0.05,10,3,0.5)) {

	oldseed <- if (exists(".Random.seed", envir = globalenv(),
						  inherits = FALSE)) {
		get(".Random.seed", envir = globalenv(), inherits = FALSE)
	}else {
		NULL}
	RNGkind("user")
	
	out <- .Call("callCalibrationSimulation",as.integer(seed),list(n=as.integer(n),runpar=as.double(runpar)),PACKAGE="microsimulation")
	states <- c("DiseaseFree","Precursor","PreClinical","Clinical")
	curnames <- names(out)
	mat <- matrix(0,nr=10,ncol=length(states))
	colnames(mat) <- states; rownames(mat) <- seq(10,100,10)
	mat[,states[states %in% curnames]]<-data.matrix(transform(as.data.frame(out[states[states %in% curnames]])))
	
	if (!is.null(oldseed)){
		assign(".Random.seed", oldseed, envir = globalenv())
	}else {
		rm(.Random.seed, envir=globalenv())
	}
	list(StateOccupancy=mat, TimeAtRisk=out$TimeAtRisk)
}
