#  PROJECT 
#  =======
#  
#  Execution time profiling of the FHCRC model.
#
#  DESCRIPTION
#  ===========
#  N.B. the profiling gets messed up by knitr make the plots separately

# AUTHOR AND DATE
# Andreas Karlsson, October 2014

## @knitr Requirements
## Profiling FHCRC
require(microsimulation)
require(ggplot2)

## @knitr Profiling_FHCRC
## Profiling fhcrc implementation
## =============================================================================
n <- 1e7
Rprof(tmp <- tempfile(),memory.profiling=T)
out <- callFhcrc(n=n, screen="noScreening")

Rprof()
prof <- summaryRprof(tmp)
data <- cbind(prof$by.self, process=rownames(prof$by.self))
data$process <- factor(data$process, levels=rownames(prof$by.self))
data <- head(data,20)
ggplot(data, aes(x =process, y = self.pct )) + geom_bar(stat="identity") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    ggtitle('By self process')
ggsave("selfProcess.pdf",
    width = 10,
    height = 8,
    dpi = 1200)


data_tot <- cbind(prof$by.total, process=rownames(prof$by.total))
data_tot$process <- factor(data_tot$process, levels=rownames(prof$by.total))
data_tot_top <- head(data_tot,12)
ggplot(data_tot_top, aes(x =process, y = total.pct , fill = self.pct/total.pct)) + geom_bar(stat="identity") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle('By parent process coloured by self portion')
ggsave("parentColByPortion.pdf")
unlink(tmp)
## =============================================================================
