# This script compares results of the half giant component destruction and creates plots

library(ggplot2)
library(gridExtra)
library(dplyr)


languageEl <- "EN"

if (languageEl == "EN")
{
  ytxt <- "Avg. removal percentage to destroy half GC"
  xtxt <- "Network"
  xtxt2 <- "Species in Giant Component"
} else {
  ytxt <- "Porcentaje de extinciones para destruir media CG"
  xtxt <- "Red"
  xtxt2 <- "Especies en la componente gigante"
}

comparativa <- function(results_by_r,baseIndex = "krisk", bestindexline = FALSE)
{
  plabels <- gsub("M_","",unlist(strsplit(as.character(results_by_r[results_by_r$Index == baseIndex,]$Network),".csv")))
  pl <- ggplot(data = results_by_r) +
    scale_x_continuous(name="",breaks = results_by_r[results_by_r$Index == baseIndex,]$index, labels = plabels)+
    ggtitle("Destruction removing species of both classes")+
    geom_point(data = results_by_r,
               aes(x = index, y = comp_perf, shape = Index, fill = Index, color = Index), size = 2, alpha= 0.45)+
    scale_color_manual(values  = cols) +
    scale_shape_manual(values = pshapes) +
    scale_y_continuous(name =ytxt)+
    theme_bw() + theme(axis.text.x  = element_text(face="bold", margin = unit(c(1, 1, 1, 1), "mm"),angle=90, vjust=0, size=8),
                       axis.title.x = element_text(face="bold",color="grey30", size=14),
                       axis.title.y = element_text(face="bold",color="grey30", size=11),
                       panel.grid.minor = element_blank(),
                       plot.title = element_text(hjust = 0.5),
                       axis.text.y = element_text(face="bold", color="grey30", size=12)
    )
  if (bestindexline){
    results_best <- results_by_r[results_by_r$Index == baseIndex,]
    pl <- pl + geom_line(data = results_best, aes(x = index, y = comp_perf),
                         color = cols[baseIndex], size = 0.5, alpha= 0.45, linetype = 1)
  }
  return(pl)
}
# comparativa <- function(results_by_r,baseIndex = "krisk")
# {
#   plabels <- gsub("M_","",unlist(strsplit(as.character(results_by_r[results_by_r$Index == baseIndex,]$Network),".csv")))
#   pl <- ggplot(data = results_by_r) +
#   #geom_line(aes(x = index, y = comp_perf, fill = Index, color = Index))+
#   scale_x_continuous(name="",breaks = results_by_r[results_by_r$Index == baseIndex,]$index, labels = plabels)+
#   ggtitle("Destruction removing species of both classes")+
#   geom_point(data = results_by_r,
#             aes(x = index, y = comp_perf, shape = Index, fill = Index, color = Index), size = 2, alpha= 0.45)+
#   scale_color_manual(values  = cols) +
#   scale_shape_manual(values = pshapes) +
#   scale_y_continuous(name =ytxt,breaks=c(0,25,50),labels=c("50%","25%","0%"), limits=c(0,55))+
#   theme_bw() + theme(axis.text.x  = element_text(face="bold", margin = unit(c(1, 1, 1, 1), "mm"),angle=90, vjust=0, size=8),
#                      axis.title.x = element_text(face="bold",color="grey30", size=14),
#                      axis.title.y = element_text(face="bold",color="grey30", size=11),
#                      panel.grid.minor = element_blank(),
#                      plot.title = element_text(hjust = 0.5),
#                      axis.text.y = element_text(face="bold", color="grey30", size=12)
#                      )
#   return(pl)
# }
results_ext <- read.csv("extinctions/Destructions_1st_mean.csv")

results_ext <- results_ext[order(-results_ext$krisk),]
results_by_row <- data.frame( Network = c(), Index = c(), performance = c())
for (i in 1:nrow(results_ext))
{
  results_by_row <- rbind(results_by_row, data.frame( Network = results_ext$Network[i], Index = "krisk", performance = 100*results_ext$krisk[i]/results_ext$giant_component[i]))
  results_by_row <- rbind(results_by_row, data.frame( Network = results_ext$Network[i], Index = "degree", performance = 100*results_ext$degree[i]/results_ext$giant_component[i]))
  results_by_row <- rbind(results_by_row, data.frame( Network = results_ext$Network[i], Index = "kdegree", performance = 100*results_ext$kdegree[i]/results_ext$giant_component[i]))
  results_by_row <- rbind(results_by_row, data.frame( Network = results_ext$Network[i], Index = "eigen", performance = 100*results_ext$eigen[i]/results_ext$giant_component[i]))
  results_by_row <- rbind(results_by_row, data.frame( Network = results_ext$Network[i], Index = "best", performance = 100*results_ext$best[i]/results_ext$giant_component[i]))
}

alpha_level = 0.5
p <- ggplot(data = results_ext) + geom_point(data = results_ext, aes(x = giant_component,y = krisk/giant_component), color = "red", alpha = alpha_level) +
     geom_point(data = results_ext, aes(x = giant_component,y = kdegree/giant_component), color = "blue", alpha = alpha_level) +
     scale_y_log10()

aux_ord_df <- data.frame(Network = results_ext$Network,
                         #value = (results_ext$krisk+results_ext$kdegree+results_ext$degree+results_ext$eigen)/results_ext$giant_component)
                         value = (results_ext$krisk)/results_ext$giant_component)

aux_ord_df <- aux_ord_df[order(aux_ord_df$value),]
aux_ord_df$index <- seq(1,nrow(aux_ord_df))

results_by_row$index <- 0
# maxperf <- max(results_by_row$performance)
# results_by_row$comp_perf <- maxperf - results_by_row$performance
results_by_row$comp_perf <- results_by_row$performance
for (i in 1:nrow(aux_ord_df))
  results_by_row[results_by_row$Network == aux_ord_df$Network[i],]$index <- aux_ord_df$index[i]

results_by_row <- results_by_row[order(results_by_row$index),]
cols <- c("kdegree" = "darkgreen", "eigen" = "black","degree" = "blue", "krisk" = "red", "best" = "forestgreen")
pshapes <- c("kdegree" = 17, "eigen" = 8,"degree" = 15, "krisk" = 16, "best" = 24)


results_by_r <- results_by_row[is.element(results_by_row$Index, c("krisk","kdegree","degree","eigen")),]
todos <- comparativa(results_by_r, baseIndex = "krisk", bestindexline = TRUE)

results_by_r <- results_by_row[is.element(results_by_row$Index, c("krisk","kdegree")),]
r <- comparativa(results_by_r)

results_by_r <- results_by_row[is.element(results_by_row$Index, c("krisk","degree")),]
s <- comparativa(results_by_r)

results_by_r <- results_by_row[is.element(results_by_row$Index, c("krisk","eigen")),]
t <- comparativa(results_by_r)

results_by_r <- results_by_row[is.element(results_by_row$Index, c("degree","kdegree")),]
u <- comparativa(results_by_r, baseIndex = "degree")

results_by_r <- results_by_row[is.element(results_by_row$Index, c("krisk","best")),]
v <- comparativa(results_by_r)

results_by_q <- results_by_row[is.element(results_by_row$Index, c("krisk","degree","kdegree","eigen")),]
results_by_q$giant_component <- 0
for (i in 1:nrow(results_by_q))
  results_by_q$giant_component[i] <- results_ext[results_ext$Network == results_by_q$Network[i],]$giant_component

results_o <- inner_join(results_by_q,results_ext,key=Network)
is_best <- c()
results_by_q$isbest <- FALSE

for (m in 1:nrow(results_by_q))
{
  results_by_q$isbest[m] <- results_o[results_o$Network == results_by_q[m,]$Network,]$performance == results_by_q$performance[m]
}

results_by_q_best <- results_by_q[results_by_q$isbest,]
q <- ggplot(results_by_q, aes(x=giant_component, y = comp_perf, color = Index)) + geom_point(alpha = 0.3, size=2) + scale_x_log10() + xlab(xtxt2)+
  scale_color_manual(values  = cols) +
    theme_bw()  +
   scale_y_continuous(name =ytxt,breaks=c(0,25,50),labels=c("50%","25%","0%"), limits=c(0,50))+
  theme(          axis.title.x = element_text(face="bold",color="grey30", size=10),
                  axis.title.y = element_text(face="bold",color="grey30", size=10),
                  axis.text.x = element_text(face="bold", color="grey30", size=10),
                  axis.text.y = element_text(face="bold", color="grey30", size=10))
q <- q + facet_grid(Index ~.)

qb <- ggplot(results_by_q_best, aes(x=giant_component, y = comp_perf, color = Index, shape = Index)) +
  geom_point( size = 3,alpha = 0.4) + scale_x_log10() + xlab(xtxt2)+
  scale_color_manual(values  = cols) + ggtitle("Destruction removing species of both classes. Top performer")+
  theme_bw()  +
  scale_shape_manual(values = pshapes) +
  scale_y_continuous(name =ytxt,breaks=c(0,25,50),labels=c("50%","25%","0%"), limits=c(0,50))+
  theme(          axis.title.x = element_text(face="bold",color="grey30", size=10),
                  axis.title.y = element_text(face="bold",color="grey30", size=10),
                  axis.text.x = element_text(face="bold", color="grey30", size=10),
                  axis.text.y = element_text(face="bold", color="grey30", size=10),
                  plot.title = element_text(hjust = 0.5)
                 )

mo <- lm(formula = results_by_q$performance ~ log(results_by_q$giant_component))
summary(mo)

ppi <- 300
png(paste0("graphs/1stalg_all_comparison.png"), width=(15*ppi), height=5*ppi, res=ppi)
print(todos)
dev.off()

ppi <- 300
png(paste0("graphs/1st_alg_best_sizeGC.png"), width=(15*ppi), height=4.5*ppi, res=ppi)
print(qb)
dev.off()


num_redes <- nrow(results_ext)

bkrisk <- sum(results_ext$krisk == results_ext$best)
print(sprintf("krisk is the best for %d networks (%.2f%%)",bkrisk,100*bkrisk/num_redes))

bkdegree <- sum(results_ext$kdegree == results_ext$best)
print(sprintf("kdegree is the best for %d networks (%.2f%%)",bkdegree,100*bkdegree/num_redes))

bdegree <- sum(results_ext$degree == results_ext$best)
print(sprintf("degree is the best for %d networks (%.2f%%)",bdegree,100*bdegree/num_redes))

beigen <- sum(results_ext$eigen == results_ext$best)
print(sprintf("eigen is the best for %d networks (%.2f%%)",beigen,100*beigen/num_redes))
