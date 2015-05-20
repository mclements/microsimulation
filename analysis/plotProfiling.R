require(ggplot2)
df <- read.csv("./implementationProfiling.csv", header=T)
ggplot(data = df, aes(x = cores, y = time)) + geom_line(aes(group = Setting, colour = Setting)) + xlab("No. of cores") + ylab("Time [s]") + ggtitle("Microsimulation of 'no screening' with 1e7 individuals") 

ggsave("implementationProfiling.pdf")
