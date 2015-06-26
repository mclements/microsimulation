require(ggplot2)
df <- read.csv("./implementationProfiling.csv", header=T)
df$Setting <- factor(df$Setting, levels(df$Setting)[c(3,2,1)])
refTime <- min(df$time[df$cores==1])
ggplot(data = df, aes(x = cores, y = time)) +
    geom_line(alpha=0.8, size=1, aes(group = Setting, colour = Setting)) +
    xlab("No. of cores") + ylab("Time [s]") +
    ggtitle("Microsimulation of 'no screening' with 1e7 individuals") +
    scale_x_log10(breaks = 2^seq(1,3,1), minor_breaks = NULL) +
    scale_y_log10(breaks = seq(100,500,100), minor_breaks = NULL) +
    geom_abline(intercept = log10(refTime), slope = -1, size = 0.3, linetype = "dashed")
ggsave("../doc/report/images/implementationProfiling.pdf")


df <- read.csv("./flagsProfiling.csv", header=T)
refTime <- min(df$time[df$cores==1])
ggplot(data = df, aes(x = cores, y = time)) +
    geom_line(alpha=0.8, size=1, aes(group = Setting, colour = Setting, linetype=Setting)) +
    xlab("No. of cores") + ylab("Time [s]") +
    ggtitle("Microsimulation of 'no screening' with 1e7 individuals") +
    scale_x_log10(breaks = 2^seq(1,3,1), minor_breaks = NULL) +
    scale_y_log10(breaks = seq(100,500,100), minor_breaks = NULL) +
    geom_abline(intercept = log10(refTime), slope = -1, size = 0.3, linetype = "dashed")
    ggsave("../doc/report/images/flagsProfiling.pdf")


df <- read.csv("./multinode.csv", header=T)
refTime <- min(df$time[df$cores==1])
levels(df$Setting) <- c("openMP",  "MPI", "MPI + openMP")
ggplot(data = df, aes(x = cores, y = time)) +
    geom_line(alpha=0.8, size=1, aes(group = Setting, colour = Setting)) +
    xlab("No. of cores") + ylab("Time [s]") +
    ggtitle("Microsimulation of 'no screening' with 1e7 individuals") + 
    scale_x_log10(breaks = 2^seq(1,7,1), minor_breaks = NULL) +
    scale_y_log10(breaks = c(seq(20,80,20), seq(100,500,100)), minor_breaks = NULL) +
    geom_abline(intercept = log10(refTime), slope = -1, size = 0.3, linetype = "dashed")
ggsave("../doc/report/images/multiNodeProfiling.pdf")
