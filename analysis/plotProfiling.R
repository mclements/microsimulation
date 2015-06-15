require(ggplot2)
df <- read.csv("./implementationProfiling.csv", header=T)
df$Setting <- factor(df$Setting, levels(df$Setting)[c(3,2,1)])
ggplot(data = df, aes(x = cores, y = time, fill=Setting)) +
    geom_bar(stat="identity", position="dodge")+
    xlab("No. of cores") + ylab("Time [s]") +
    ggtitle("Microsimulation of 'no screening' with 1e7 individuals") +
  scale_fill_hue(name = "Implementation")

ggsave("../doc/report/images/implementationProfiling.pdf")


df <- read.csv("./flagsProfiling.csv", header=T)
ggplot(data = df, aes(x = cores, y = time, fill = Setting)) +
  #geom_line(aes(group = Setting, colour = Setting)) +
  geom_bar(stat="identity", position="dodge") + xlab("No. of cores") + ylab("Time [s]") + ggtitle("Microsimulation of 'no screening' with 1e7 individuals") + scale_fill_hue(name = "Optimisation flags")

## ggplot(data = df, aes(x = cores, y = time)) + geom_line(aes(group = Setting, colour = Setting)) +
##   xlab("No. of cores") + ylab("Time [s]") + ggtitle("Microsimulation of 'no screening' with 1e7 individuals")

ggsave("../doc/report/images/flagsProfiling.pdf")
