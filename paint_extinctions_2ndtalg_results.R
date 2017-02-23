# This script compares results of the half giant component destruction and creates plots

library(ggplot2)
library(gridExtra)
library(dplyr)


comparativa <- function(results_by_r,baseIndex = "krisk", bestindexline = FALSE)
{
  plabels <- gsub("M_","",unlist(strsplit(as.character(results_by_r[results_by_r$Index == baseIndex,]$Network),".csv")))
  pl <- ggplot(data = results_by_r) +
  scale_x_continuous(name="",breaks = results_by_r[results_by_r$Index == baseIndex,]$index, labels = plabels)+
  ggtitle("Destruction removing animal species")+
  geom_point(data = results_by_r,
            aes(x = index, y = comp_perf, shape = Index, fill = Index, color = Index), size = 3, alpha= 0.45)+
  scale_color_manual(values  = cols) +
  scale_shape_manual(values = pshapes) +
  scale_y_continuous(name =ytxt)+
    theme_bw() + theme(axis.text.x  = element_text(face="bold", angle=90, hjust= 1,vjust=0.75, size=9),
                       axis.title.x = element_text(face="bold",color="grey30", size=15),
                       axis.title.y = element_text(face="bold",color="grey30", size=15),
                       panel.grid.minor = element_blank(),
                       panel.grid.major = element_line(color="grey30", size=0.5, linetype = 3),
                       legend.position="bottom",
                       legend.title = element_blank(),
                       legend.text = element_text(face="bold", color="grey30", size=15),
                       plot.title = element_text(hjust = 0.5, size=18),
                       axis.text.y = element_text(face="bold", color="grey30", size=15)
                     )
  if (bestindexline){
    results_best <- results_by_r[results_by_r$Index == baseIndex,]
    pl <- pl + geom_line(data = results_best, aes(x = index, y = comp_perf),
                         color = cols[baseIndex], size = 0.5, alpha= 0.45, linetype = 1)
  }
  return(pl)
}

listamethods <- c("GCdestroymethod","PLdestroymethod")

for (method in listamethods)
{
  if (method == "GCdestroymethod")
  {
    ytxt <- "AUC (Giant Component size)"
    xtxt <- "Network"
    xtxt2 <- "Species in Giant Component"
  } else {
    ytxt <- "AUC (suviving vegetal species)"
    xtxt <- "Red"
    xtxt2 <- "Especies en la componente gigante"
  }

  results_ext <- read.csv(paste0("extinctions/Destructions_2nd_mean_",method,".csv"))
  names(results_ext) <- c("X","Network","NoOrder","MusRank","krisk","kdegree","degree","eigenc")
  if (method == "GCdestroymethod")
    results_ext <- results_ext[order(-results_ext$kdegree),]
  results_ext$best <- min(results_ext[!names(results_ext) %in% c("Network","X","NoOrder")])
  results_by_row <- data.frame( Network = c(), Index = c(), performance = c())
  for (i in 1:nrow(results_ext))
  {
    results_by_row <- rbind(results_by_row, data.frame( Network = results_ext$Network[i], Index = "krisk", performance = results_ext$krisk[i]))
    results_by_row <- rbind(results_by_row, data.frame( Network = results_ext$Network[i], Index = "MusRank", performance = results_ext$MusRank[i]))
    results_by_row <- rbind(results_by_row, data.frame( Network = results_ext$Network[i], Index = "degree", performance = results_ext$degree[i]))
    results_by_row <- rbind(results_by_row, data.frame( Network = results_ext$Network[i], Index = "kdegree", performance = results_ext$kdegree[i]))
    results_by_row <- rbind(results_by_row, data.frame( Network = results_ext$Network[i], Index = "eigen", performance = results_ext$eigen[i]))
    results_by_row <- rbind(results_by_row, data.frame( Network = results_ext$Network[i], Index = "best", performance = results_ext$best[i]))
  }

  # alpha_level = 0.5
  # p <- ggplot(data = results_ext) + geom_point(data = results_ext, aes(x = giant_component,y = krisk/giant_component), color = "red", alpha = alpha_level) +
  #      geom_point(data = results_ext, aes(x = giant_component,y = kdegree/giant_component), color = "blue", alpha = alpha_level) +
  #      scale_y_log10()

  if (method == "GCdestroymethod"){
    aux_ord_df <- data.frame(Network = results_ext$Network, value = (results_ext$kdegree))
  } else
    aux_ord_df <- data.frame(Network = results_ext$Network, value = (results_ext$MusRank))

  aux_ord_df <- aux_ord_df[order(aux_ord_df$value),]
  aux_ord_df$index <- seq(1,nrow(aux_ord_df))

  results_by_row$index <- 0
  # maxperf <- max(results_by_row$performance)
  # results_by_row$comp_perf <- maxperf - results_by_row$performance
  results_by_row$comp_perf  <- results_by_row$performance
  for (i in 1:nrow(aux_ord_df))
    results_by_row[results_by_row$Network == aux_ord_df$Network[i],]$index <- aux_ord_df$index[i]

  results_by_row <- results_by_row[order(results_by_row$index),]
  cols <- c("kdegree" = "darkgreen", "eigen" = "black","degree" = "blue", "krisk" = "red", "best" = "forestgreen", "MusRank" = "bisque3")
  pshapes <- c("kdegree" = 17, "eigen" = 18,"degree" = 15, "krisk" = 16, "best" = 24, "MusRank" = 25)

  results_by_r <- results_by_row[is.element(results_by_row$Index, c("krisk","kdegree","degree","eigen","MusRank")),]
  if (method == "GCdestroymethod"){
    todos <- comparativa(results_by_r, baseIndex= "kdegree", bestindexline = TRUE)
  } else
    todos <- comparativa(results_by_r, baseIndex= "MusRank", bestindexline = TRUE)

  # results_by_q <- results_by_row[is.element(results_by_row$Index, c("krisk","degree","kdegree","eigen")),]
  # results_by_q$giant_component <- 0
  # for (i in 1:nrow(results_by_q))
  #   results_by_q$giant_component[i] <- results_ext[results_ext$Network == results_by_q$Network[i],]$giant_component
  #
  # results_o <- inner_join(results_by_q,results_ext,key=Network)
  # is_best <- c()
  # results_by_q$isbest <- FALSE
  #
  # for (m in 1:nrow(results_by_q))
  # {
  #   results_by_q$isbest[m] <- results_o[results_o$Network == results_by_q[m,]$Network,]$performance == results_by_q$performance[m]
  # }
  #
  # results_by_q_best <- results_by_q[results_by_q$isbest,]
  # q <- ggplot(results_by_q, aes(x=giant_component, y = comp_perf, color = Index)) + geom_point(alpha = 0.3, size=2) + scale_x_log10() + xlab(xtxt2)+
  #   scale_color_manual(values  = cols) +
  #     theme_bw()  +
  #    scale_y_continuous(name =ytxt,breaks=c(0,25,50),labels=c("50%","25%","0%"), limits=c(0,50))+
  #   theme(          axis.title.x = element_text(face="bold",color="grey30", size=10),
  #                   axis.title.y = element_text(face="bold",color="grey30", size=10),
  #                   axis.text.x = element_text(face="bold", color="grey30", size=10),
  #                   axis.text.y = element_text(face="bold", color="grey30", size=10))
  # q <- q + facet_grid(Index ~.)
  #
  # qb <- ggplot(results_by_q_best, aes(x=giant_component, y = comp_perf, color = Index, shape = Index)) +
  #   geom_point( size = 3,alpha = 0.4) + scale_x_log10() + xlab(xtxt2)+
  #   scale_color_manual(values  = cols) + ggtitle("Destruction removing species of both classes. Top performer")+
  #   theme_bw()  +
  #   scale_shape_manual(values = pshapes) +
  #   scale_y_continuous(name =ytxt,breaks=c(0,25,50),labels=c("50%","25%","0%"), limits=c(0,50))+
  #   theme(          axis.title.x = element_text(face="bold",color="grey30", size=10),
  #                   axis.title.y = element_text(face="bold",color="grey30", size=10),
  #                   axis.text.x = element_text(face="bold", color="grey30", size=10),
  #                   axis.text.y = element_text(face="bold", color="grey30", size=10),
  #                   plot.title = element_text(hjust = 0.5)
  #                  )
  #
  # mo <- lm(formula = results_by_q$performance ~ log(results_by_q$giant_component))
  # summary(mo)

  ppi <- 300
  png(paste0("graphs/2ndtalg_all_comparison_",method,".png"), width=(15*ppi), height=6*ppi, res=ppi)
  print(todos)
  dev.off()

  resdest <- read.csv(paste0("python/results/mean_all_",method,".csv"))
  resdest$NetworkFile <- paste0("M_",resdest$Network,".csv")
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
  print(paste0("Extinction method: ",  method))
  print("==============================================================")
  print(paste0("MR is the best for ",countMR," networks (",100*countMR/nrow(resdest),"%)"))
  print(paste0("Krisk is the best for ",countKrisk," networks (",100*countKrisk/nrow(resdest),"%)"))
  print(paste0("Kdegree is the best for ",countKdegree," networks (",100*countKdegree/nrow(resdest),"%)"))
  print(paste0("Degree is the best for ",countDegree," networks (",100*countDegree/nrow(resdest),"%)"))
  print(paste0("eigenc is the best for ",counteigenc," networks (",100*counteigenc/nrow(resdest),"%)"))
}
