LdFlags <- function()
    cat(system.file("lib/libmicrosimulation.a",  package="microsimulation", mustWork=TRUE))

microsimulation.init <- function(PACKAGE="microsimulation") {
  .Call("r_create_current_stream",PACKAGE=PACKAGE)
  return(1)
}

microsimulation.exit <- function(PACKAGE="microsimulation") {
  .Call("r_remove_current_stream",PACKAGE=PACKAGE)
  return(1)
}
## http://r.789695.n4.nabble.com/How-to-construct-a-valid-seed-for-l-Ecuyer-s-method-with-given-Random-seed-td4656340.html
unsigned <- function(seed) ifelse(seed < 0, seed + 2^32, seed)
signed <- function(seed) ifelse(seed>2^31, seed-2^32, seed)

rnormPos <- function(n,mean=0,sd=1,lbound=0) {
    if (length(mean)<n) mean <- rep(mean,length=n)
    if (length(sd)<n) sd <- rep(sd,length=n)
    x <- rnorm(n,mean,sd)
    while(any(i <- which(x<lbound)))
        x[i] <- rnorm(length(i),mean[i],sd[i])
    x
}

set.user.Random.seed <- function (seed,PACKAGE="microsimulation") {
  seed <- as.double(unsigned(seed))
  if (length(seed) == 1) seed <- rep(seed,6)
  if (length(seed) == 7) seed <- seed[-1]
  .C("r_set_user_random_seed",seed = seed,PACKAGE=PACKAGE)
  return(invisible(seed))
}

advance.substream <- function (seed,n,PACKAGE="microsimulation") {
  seed <- as.double(unsigned(seed))
  if (length(seed) == 1) seed <- rep(seed,6)
  if (length(seed) == 7) seed <- seed[-1]
  .C("r_rng_advance_substream", seed = seed, n = as.integer(n), PACKAGE=PACKAGE)$seed
}

next.user.Random.substream <- function(PACKAGE="microsimulation") {
  .C("r_next_rng_substream", PACKAGE=PACKAGE)
  return(invisible(TRUE))
}

user.Random.seed <- function(PACKAGE="microsimulation") {
  c(407L,
    as.integer(signed(.C("r_get_user_random_seed", seed=as.double(rep(1,6)),
                         PACKAGE=PACKAGE)$seed)))
}

enum <- function(obj, labels, start=0) {
  if (is.logical(obj)) obj <- obj+0
  structure(factor(obj, levels=start + (0:(length(labels)-1)), labels=labels),start=start)
}

"enum<-" <- function(obj, value) {
  enum(if(is.factor(obj)) unclass(obj)-1+attr(obj,"start") else obj, value)
}

RNGstate <- function() {
  ## house-keeping for random streams (see parallel::clusterSetRNGStream)
  oldseed <- if (exists(".Random.seed", envir = .GlobalEnv,
                        inherits = FALSE))
    get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  else NULL
  oldkind <- RNGkind()
  reset <- function() {
    if (!is.null(oldseed))
      assign(".Random.seed", oldseed, envir = .GlobalEnv)
    else {
        ## clean up if we created a .Random.seed
        if (exists(".Random.seed" ,envir = .GlobalEnv, inherits = FALSE))
            rm(.Random.seed, envir = .GlobalEnv, inherits = FALSE)
            do.call(RNGkind, as.list(oldkind))
    }
  }
  list(oldseed = oldseed, reset = reset)
}

## find the concave frontier
frontier<-function(x,y)
  {
      ichull <- grDevices::chull(cbind(x,y)) # convex hull
      if (length(ichull)<2) return(ichull)
      xi <- x[ichull]
      yi <- y[ichull]          # subset to convex hull
      imin <- which(xi==min(xi))
      include <- sapply(1:length(ichull),
                        function(i)       # establish the frontier
                        i==imin || (i>imin && yi[i-1]<yi[i] && xi[i-1]<xi[i]))
      ichull[include]
}
lines.frontier <- function(x,y,pch=19,type="b",...) {
    index <- frontier(x,y)
    lines(x[index],y[index],pch=pch,type=type,...)
}


callPersonSimulation <- function(n=20,seed=rep(12345,6)) {
  state <- RNGstate(); on.exit(state$reset())
  RNGkind("user")
  set.user.Random.seed(seed)
  stateT =c("Healthy","Localised","DxLocalised","LocallyAdvanced",
    "DxLocallyAdvanced","Metastatic","DxMetastatic","Death")
  eventT =c("toDeath", "toPCDeath", "toLocalised", "toDxLocalised",
	      "toDxLocallyAdvanced",
	      "toLocallyAdvanced", "toMetastatic", "toDxMetastatic")
  out <- .Call(.callPersonSimulation,
               as.integer(rep(seed,length=6)), # magic
               list(n=as.integer(n)))
  out$state = enum(out$state,stateT)
  out$event = enum(out$event,eventT)
  as.data.frame(out)
}

callSimplePerson <- function(n=10) {
  state <- RNGstate(); on.exit(state$reset())
  RNGkind("Mersenne-Twister")
  set.seed(12345)
  stateT <- c("Healthy","Cancer","Death")
  eventT <- c("toOtherDeath", "toCancer", "toCancerDeath")
  out <- .Call(.callSimplePerson,
               parms=list(n=as.integer(n)))
  out$state = enum(out$state,stateT)
  out$event = enum(out$event,eventT)
  as.data.frame(out)
}

callSimplePerson2 <- function(n=10) {
  state <- RNGstate(); on.exit(state$reset())
  RNGkind("Mersenne-Twister")
  set.seed(12345)
  stateT <- c("Healthy","Cancer","Death")
  eventT <- c("toOtherDeath", "toCancer", "toCancerDeath")
  out <- .Call(.callSimplePerson2,
               parms=list(n=as.integer(n)))
  reader <- function(obj)
      cbind(data.frame(state=enum(obj[[1]],stateT)),
          data.frame(obj[-1]))
  out <- lapply(out,reader)
  out$events <- with(out$events, data.frame(state=state,event=enum(event,eventT),age=age,number=number))
  out$pt <- with(out$pt, data.frame(state=state,age=age,pt=pt))
  out$prev <- with(out$prev, data.frame(state=state,age=age,number=number))
  out
}


## readEventReport <- function(obj) {
##     pt <- obj$pt
##     n <- ncol(pt)
##     names(pt)[(n-1):n] <- c("age","pt")
## }

callIllnessDeath <- function(n=10L,cure=0.1,zsd=0) {
  state <- RNGstate(); on.exit(state$reset())
  RNGkind("Mersenne-Twister")
  set.seed(12345)
  stateT <- c("Healthy","Cancer")
  eventT <- c("toOtherDeath", "toCancer", "toCancerDeath")
  out <- .Call(.callIllnessDeath,
               parms=list(n=as.integer(n),cure=as.double(cure),zsd=as.double(zsd)))
  reader <- function(obj)
      cbind(data.frame(state=enum(obj[[1]],stateT)),
          data.frame(obj[-1]))
  out <- lapply(out,reader)
  out$events <- with(out$events, data.frame(state=state,event=enum(event,eventT),age=age,number=number))
  out$pt <- with(out$pt, data.frame(state=state,age=age,pt=pt))
  out$prev <- with(out$prev, data.frame(state=state,age=age,number=number))
  out
}

ICER <- function(object1, object2, ...)
    UseMethod("ICER")

## utility - not exported
assignList <- function(lst,...)
  for(i in 1:length(lst))
    assign(names(lst)[i], lst[[i]], ...)
## assignList(formals(callFhcrc),pos=1)

.testPackage <- function() {
  list(callSimplePerson(),
       callPersonSimulation(n=10),
       callSimplePerson2())
}

EventQueue <-
    setRefClass("EventQueue",
                fields = list(times = "numeric", events = "list"),
                methods = list(
                    help = function() {
                        'Reference class implementation of an event queue. Fields for the event times and the events list.'
                    },
                    insert = function(time,event) {
                        'Method to insert the event at the given time'
                        insert.ord <- findInterval(time,times)
                        times <<- append(times,time,insert.ord)
                        events <<- append(events,event,insert.ord)
                    },
                    pop = function() {
                        'Method to remove the head of the event queue and return its value'
                        head <- structure(events[[1]], time=times[1])
                        times <<- times[-1]
                        events <<- events[-1]
                        return(head)
                    },
                    empty = function() {
                        'Method to check whether there are no events in the queue'
                        length(times) == 0
                        },
                    clear = function() {
                        'Method to clear the event queue'
                        times <<- numeric()
                        events <<- list()
                    },
                    remove = function(predicate, ...) {
                        'Method to remove events that satisfy some predicate'
                        i <- sapply(events, predicate, ...)
                        stopifnot(is.logical(i))
                        i[is.na(i)] <- TRUE
                        times <<- times[!i]
                        events <<- events[!i]
                    }))

BaseDiscreteEventSimulation <-
    setRefClass("BaseDiscreteEventSimulation",
                contains = "EventQueue",
                fields = list(currentTime = "numeric",
                    previousEventTime = "numeric"),
                methods = list(
                    help = function() {
                        'Reference class implementation of an event-oriented discrete event simulation. Fields for the event times and the events list.'
                    },
                    scheduleAt = function(time, event) {
                        'Method that adds attributes for the event time and the sendingTime, and then insert the event into the event queue'
                        attr(event,"time") <- time
                        attr(event,"sendingTime") <- currentTime
                        insert(time, event)
                    },
                    init = function() {
                        'Virtual method to initialise the event queue and attributes'
                        stop("VIRTUAL!")
                    },
                    handleMessage = function(event) {
                        'Virtual method to handle the messages as they arrive'
                        stop("VIRTUAL!")
                    },
                    final = function() {
                        'Method for finalising the simulation'
                        NULL
                    },
                    now = function() currentTime,
                    reset = function(startTime = 0.0) {
                        'Method to reset the event queue'
                        clear()
                        previousEventTime <<- currentTime <<- startTime
                    },
                    run = function(startTime = 0.0) {
                        'Method to run the simulation'
                        reset(startTime)
                        init()
                        while(!empty()) {
                            event <- pop()
                            currentTime <<- attr(event, "time")
                            handleMessage(event)
                            previousEventTime <<- currentTime
                        }
                        final()
                    }))

RNGStream <- function(nextStream = TRUE, iseed = NULL) {
  stopifnot(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE) || !is.null(iseed))
  oldseed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
    get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  else NULL
  if (RNGkind()[1] != "L'Ecuyer-CMRG") RNGkind("L'Ecuyer-CMRG")
  if (!is.null(iseed)) set.seed(iseed)
  current <- if (nextStream) parallel::nextRNGStream(.Random.seed) else .Random.seed
  .Random.seed <<- startOfStream <- startOfSubStream <- current
  structure(list(resetRNGkind = function() {
      if (!is.null(oldseed))
          assign(".Random.seed", oldseed, envir = .GlobalEnv)
      else rm(.Random.seed, envir = .GlobalEnv)
  },
                 seed = function() current,
                 open = function() .Random.seed <<- current,
                 close = function() current <<- .Random.seed,
                 resetStream = function() .Random.seed <<- current <<- startOfSubStream <<- startOfStream,
                 resetSubStream = function() .Random.seed <<- current <<- startOfSubStream,
                 nextSubStream = function() .Random.seed <<- current <<- startOfSubStream <<- parallel::nextRNGSubStream(startOfSubStream),
                 nextStream = function() .Random.seed <<- current <<- startOfSubStream <<- startOfStream <<- parallel::nextRNGStream(startOfStream)),
            class="RNGStream")
  }

with.RNGStream <- function(data,expr,...) {
  data$open()
  out <- eval(substitute(expr), enclos = parent.frame(), ...)
  data$close()
  out
}
setOldClass("RNGStream")
