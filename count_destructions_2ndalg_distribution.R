library(ggplot2)
library(dplyr)

dest_method <- "GCdestroymethod"
directorystr <- "python/results"

rm("dftot")
ficheros <- Sys.glob(paste0(directorystr,"/",dest_method,"/DIAM_EXTIN_ALL*"))
for (j in ficheros)
{
  dfparc <- read.csv(j,sep="\t")
  if (exists("dftot"))
    dftot <- rbind(dftot,dfparc)
  else
    dftot <- dfparc
}
tnames <- names(dftot)
tnames[1] <- "Network"
names(dftot) <- tnames

calc_stats <- function(dfinput,index)
{
  dfsal <- data.frame( Network=NA, NoOrder=NA, MR = NA,
                       Krisk=NA,Kdegree=NA, Degree=NA, eigenc=NA)[numeric(0), ]
  redes <- unique(dfinput$Network)
  for (j in redes)
  {
    dfparcsal <- data.frame( Network="", NoOrder=0, MR = 0,
                             Krisk=0,Kdegree=0, Degree=0, eigenc=0)
    dfnet <- dfinput[dfinput$Network == j,]
    dfparcsal$Network <- j
    if (index == "median"){
      dfparcsal$NoOrder <- median(dfnet$NoOrder)
      dfparcsal$MR <- median(dfnet$MR)
      dfparcsal$Krisk <- median(dfnet$Krisk)
      dfparcsal$Kdegree <- median(dfnet$Kdegree)
      dfparcsal$Degree <- median(dfnet$Degree)
      dfparcsal$eigenc <- median(dfnet$eigenc)
    }
    if (index == "mean")
    {
      dfparcsal$NoOrder <- mean(dfnet$NoOrder)
      dfparcsal$MR <- mean(dfnet$MR)
      dfparcsal$Krisk <- mean(dfnet$Krisk)
      dfparcsal$Kdegree <- mean(dfnet$Kdegree)
      dfparcsal$Degree <- mean(dfnet$Degree)
      dfparcsal$eigenc <- mean(dfnet$eigenc)
    }
    dfsal <- rbind(dfsal, dfparcsal)
  }
  write.csv(dfsal,paste0(directorystr,"/",index,"_all_",dest_method,".csv"))
  return(dfsal)
}


dfsalmedian <- calc_stats(dftot,"median")
dfsalmean <- calc_stats(dftot,"mean")
write.csv(dfsalmean,paste0("extinctions/Destructions_2nd_mean_",dest_method,".csv"))

#boxplot(dfsalmedian[!names(dftot) %in% c("Network")])
boxplot(dfsalmean[!names(dftot) %in% c("Network")])
