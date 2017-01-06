library(dplyr)

method <- "PLdestroymethod"
index <- "mean"

resdest <- read.csv(paste0("python/results/",index,"_all_",method,".csv"))
resdest$NetworkFile <- paste0("M_",resdest$Network,".csv")
# min_giant_component <- 100
# resdest <- resdest[resdest$giant_component > min_giant_component,]

redinfo <- read.csv("extinctions/size_networks.csv")

resdest <- inner_join(resdest,redinfo)

min_giant_component <- 0
resdest <- resdest[resdest$giant_component > min_giant_component,]


resdest$best=0
countMR <- 0
countKrisk <- 0
countKdegree <- 0
countDegree <- 0
counteigenc <- 0
for (i in 1:nrow(resdest))
{
  best = min(resdest[i,3:8])
  resdest$best[i] <- best
  if (resdest$MR[i] == best)
    countMR = countMR +1
  if (resdest$Degree[i] == best)
    countDegree = countDegree+1
  if (resdest$Kdegree[i] == best)
    countKdegree = countKdegree +1
  if (resdest$Krisk[i] == best)
    countKrisk = countKrisk +1
  if (resdest$eigenc[i] == best)
    counteigenc = counteigenc +1
}
print(paste0("Extinction method: ",  method, " index: ",index," min size:",min_giant_component))
print("==============================================================")
print(paste0("MR is the best for ",countMR," networks (",100*countMR/nrow(resdest),"%)"))
print(paste0("Krisk is the best for ",countKrisk," networks (",100*countKrisk/nrow(resdest),"%)"))
print(paste0("Kdegree is the best for ",countKdegree," networks (",100*countKdegree/nrow(resdest),"%)"))
print(paste0("Degree is the best for ",countDegree," networks (",100*countDegree/nrow(resdest),"%)"))
print(paste0("eigenc is the best for ",counteigenc," networks (",100*counteigenc/nrow(resdest),"%)"))
