library(dplyr)

resdest <- read.csv(paste0("extinctions/Destructions_1st_mean.csv"))

min_giant_component <- 0
resdest <- resdest[resdest$giant_component > min_giant_component,]

countKrisk <- 0
countKdegree <- 0
countDegree <- 0
counteigenc <- 0
for (i in 1:nrow(resdest))
{
  best = resdest$best[i]
  if (resdest$degree[i] == best)
    countDegree = countDegree+1
  if (resdest$kdegree[i] == best)
    countKdegree = countKdegree +1
  if (resdest$krisk[i] == best)
    countKrisk = countKrisk +1
  if (resdest$eigenc[i] == best)
    counteigenc = counteigenc +1
}
print(paste0("First algorithm destruction method for networks with size > ",min_giant_component))
print("==============================================================")
print(paste0("Krisk is the best for ",countKrisk," networks (",100*countKrisk/nrow(resdest),"%)"))
print(paste0("Kdegree is the best for ",countKdegree," networks (",100*countKdegree/nrow(resdest),"%)"))
print(paste0("Degree is the best for ",countDegree," networks (",100*countDegree/nrow(resdest),"%)"))
print(paste0("eigenc is the best for ",counteigenc," networks (",100*counteigenc/nrow(resdest),"%)"))
