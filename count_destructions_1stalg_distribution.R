library(ggplot2)
library(dplyr)
library(beeswarm)

directorystr <- "extinctions"

if (exists("dftot"))
  rm("dftot")
ficheros <- Sys.glob(paste0(directorystr,"/ALL_EXTINCTIONS_halfGC*.csv"))

for (j in ficheros)
{
  dfparc <- read.csv(j,sep=",")
  if (exists("dftot"))
    dftot <- rbind(dftot,dfparc)
  else
    dftot <- dfparc
}


calc_stats <- function(dfinput, fn = mean, fnname = "mean")
{
  dfsal <- data.frame( Network=NA, giant_component=NA, krisk=NA, degree=NA, kdegree=NA,eigenc=NA, best =NA)[numeric(0), ]
  redes <- unique(dfinput$Network)
  for (j in redes)
  {
    dfparcsal <- data.frame( Network="", giant_component=0, NoOrder = 0, krisk = 0,
                             degree=0,kdegree=0, eigenc=0, best = 0)
    dfnet <- dfinput[dfinput$Network == j,]
    dfparcsal$Network <- j
    dfparcsal$giant_component <- mean(dfnet$giant_component)
    dfparcsal$NoOrder <- fn(dfnet$NoOrder)
    dfparcsal$krisk <- fn(dfnet$krisk)
    dfparcsal$kdegree <- fn(dfnet$kdegree)
    dfparcsal$degree <- fn(dfnet$degree)
    dfparcsal$eigenc <-fn(dfnet$eigenc)
    dfparcsal$best <- min(c(dfparcsal$NoOrder,dfparcsal$krisk,dfparcsal$kdegree,dfparcsal$degree,dfparcsal$eigenc))
    dfsal <- rbind(dfsal, dfparcsal)
  }
  write.csv(dfsal,paste0(directorystr,"/Destructions_1st_",fnname,".csv"))
  return(dfsal)
}

dfsalmean <- calc_stats(dftot, mean, fnname = "mean")

# auxdf <- data.frame(NetworkFile = dfsalmean$Network,giant_component = dfsalmean$giant_component)
# write.csv(auxdf,"extinctions/size_networks.csv",row.names = FALSE)

boxplot(dftot[!names(dftot) %in% c("Network","giant_component","NoOrder","best")])
 
