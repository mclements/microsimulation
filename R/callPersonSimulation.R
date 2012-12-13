callPersonSimulation <-
function(n=500L)
  .C("callPersonSimulation",as.double(1),as.integer(n),out=as.double(1:2),as.integer(2),PACKAGE="microsimulation")$out
