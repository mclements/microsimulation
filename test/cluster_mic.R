library(Rmpi)
library(snow)
library(parallel)
require(microsimulation)

mc.cores <- max(1, mpi.universe.size() - 1)
cl <- makeMPIcluster(mc.cores)
cat(sprintf("Running with %d workers\n", length(cl)))
clusterCall(cl, function() { library(microsimulation); NULL })

callFhcrc <- 
function (n = 10, screen = "noScreening", nLifeHistories = 10, 
    screeningCompliance = 0.75, seed = 12345, studyParticipation = 50/260, 
    psaThreshold = 3, mc.cores, cl) 
{
    state <- RNGstate()
    on.exit(state$reset())
    RNGkind("user")
    set.user.Random.seed(seed)
    pop1 <- data.frame(cohort = 1980:1900, pop = c(rep(17239, 
        9), 16854, 16085, 15504, 15604, 16381, 16705, 16762, 
        16853, 15487, 14623, 14066, 13568, 13361, 13161, 13234, 
        13088, 12472, 12142, 12062, 12078, 11426, 12027, 11963, 
        12435, 12955, 13013, 13125, 13065, 12249, 11103, 9637, 
        9009, 8828, 8350, 7677, 7444, 7175, 6582, 6573, 6691, 
        6651, 6641, 6268, 6691, 6511, 6857, 7304, 7308, 7859, 
        7277, 8323, 8561, 7173, 6942, 7128, 6819, 5037, 6798, 
        rep(6567, 14)))
    screenT <- c("noScreening", "randomScreen50to70", "twoYearlyScreen50to70", 
        "fourYearlyScreen50to70", "screen50", "screen60", "screen70", 
        "screenUptake", "stockholm3_goteborg", "stockholm3_risk_stratified")
    stateT <- c("Healthy", "Localised", "Metastatic")
    gradeT <- c("Gleason_le_6", "Gleason_7", "Gleason_ge_8")
    eventT <- c("toLocalised", "toMetastatic", "toClinicalDiagnosis", 
        "toCancerDeath", "toOtherDeath", "toScreen", "toBiopsy", 
        "toScreenDiagnosis", "toOrganised", "toTreatment", "toCM", 
        "toRP", "toRT", "toADT")
    diagnosisT <- c("NotDiagnosed", "ClinicalDiagnosis", "ScreenDiagnosis")
    treatmentT <- c("CM", "RP", "RT")
    psaT <- c("PSA<3", "PSA>=3")
    stopifnot(screen %in% screenT)
    stopifnot(is.na(n) || is.integer(as.integer(n)))
    stopifnot(is.integer(as.integer(nLifeHistories)))
    stopifnot(is.double(as.double(screeningCompliance)))
    screenIndex <- which(screen == screenT) - 1
    if (is.na(n)) {
        cohort <- pop1$cohort[rep.int(1:nrow(pop1), times = pop1$pop)]
        n <- length(cohort)
    }
    else cohort <- sample(pop1$cohort, n, prob = pop1$pop/sum(pop1$pop), 
        replace = TRUE)
    cohort <- sort(cohort)
    chunks <- tapply(cohort, sort((0:(n - 1))%%mc.cores), I)
    currentSeed <- user.Random.seed()
    powerFun <- function(obj, FUN, n, ...) {
        for (i in 1:n) obj <- FUN(obj, ...)
        obj
    }
    initialSeeds <- Reduce(function(seed, i) powerFun(seed, parallel::nextRNGStream, 
        10), 1:mc.cores, currentSeed, accumulate = TRUE)[-1]
    ns <- cumsum(sapply(chunks, length))
    ns <- c(0, ns[-length(ns)])
    fhcrcData$prtx$Age <- as.double(fhcrcData$prtx$Age)
    fhcrcData$prtx$DxY <- as.double(fhcrcData$prtx$DxY)
    fhcrcData$prtx$G <- fhcrcData$prtx$G - 1L
    fhcrcData$pradt$Grade <- fhcrcData$pradt$Grade - 1L
    fhcrcData$pradt$Age <- as.double(fhcrcData$pradt$Age)
    fhcrcData$pradt$DxY <- as.double(fhcrcData$pradt$DxY)
    fhcrcData$survival_local <- with(fhcrcData$survival_local, 
        data.frame(Age = as.double(AgeLow), Grade = Grade, Time = as.double(Time), 
            Survival = Survival))
    fhcrcData$survival_dist <- with(fhcrcData$survival_dist, 
        data.frame(Grade = Grade, Time = as.double(Time), Survival = Survival))
    print(system.time(out <- clusterApply(cl, 1:mc.cores, function(i) {
        chunk <- chunks[[i]]
        set.user.Random.seed(initialSeeds[[i]])
        .Call("callFhcrc", parms = list(n = as.integer(length(chunk)), 
            firstId = ns[i], screen = as.integer(screenIndex), 
            nLifeHistories = as.integer(nLifeHistories), screeningCompliance = as.double(screeningCompliance), 
            studyParticipation = as.double(studyParticipation), 
            psaThreshold = as.double(psaThreshold), cohort = as.double(chunk), 
            tables = fhcrcData), PACKAGE = "microsimulation")
    })))
    cbindList <- function(obj) if (is.list(obj)) 
        do.call("cbind", lapply(obj, cbindList))
    else data.frame(obj)
    reader <- function(obj) {
        obj <- cbindList(obj)
        out <- cbind(data.frame(state = enum(obj[[1]], stateT), 
            grade = enum(obj[[2]], gradeT), dx = enum(obj[[3]], 
                diagnosisT), psa = enum(obj[[4]], psaT), cohort = obj[[5]]), 
            data.frame(obj[, -(1:5)]))
        out
    }
    summary <- lapply(seq_along(out[[1]]$summary), function(i) do.call("rbind", 
        lapply(out, function(obj) reader(obj$summary[[i]]))))
    names(summary) <- names(out[[1]]$summary)
    states <- c("state", "grade", "dx", "psa", "cohort")
    names(summary$prev) <- c(states, "age", "count")
    names(summary$pt) <- c(states, "age", "pt")
    names(summary$events) <- c(states, "event", "age", "n")
    summary <- lapply(summary, function(obj) within(obj, year <- cohort + 
        age))
    map2df <- function(obj) as.data.frame(do.call("cbind", obj))
    lifeHistories <- do.call("rbind", lapply(out, function(obj) map2df(obj$lifeHistories)))
    parameters <- map2df(out[[1]]$parameters)
    enum(summary$events$event) <- eventT
    enum(lifeHistories$state) <- stateT
    enum(lifeHistories$dx) <- diagnosisT
    enum(lifeHistories$event) <- eventT
    enum <- list(stateT = stateT, eventT = eventT, screenT = screenT, 
        diagnosisT = diagnosisT, psaT = psaT)
    out <- list(n = n, screen = screen, enum = enum, lifeHistories = lifeHistories, 
        parameters = parameters, summary = summary)
    class(out) <- "fhcrc"
    out
}

print(test <- callFhcrc(1e6, mc.cores=mc.cores, cl=cl))

stopCluster(cl)
mpi.quit()
