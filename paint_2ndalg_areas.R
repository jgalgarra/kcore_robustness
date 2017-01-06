# Normalized area extinction plots for the second destruction algorithm
# Output data at graphs
#
#
# Parameters:
#  red                  Network name
#  languaEl             Language "ES","EN"
#  metodo               "GCdestroymethod" (measures remaining GC), "PLdestroymethod" (measures surviving plant species)
# Requires:
#   Extinctions data at "extinctions/"

library(ggplot2)
library(grid)
library(gridExtra)

pintacurva <- function(datos_f,area, pcolor, texty, titulo = "")
{
  
  ls <- ggplot(data = datos_f, aes(x = Primary*100, y = RemainingGC*100)) + geom_area(color = pcolor,fill=pcolor,alpha=0.3) +
    theme_bw() + ylab(texty)+xlab(paste0(xtext,crit))+ggtitle(titulo)+xlim(c(0,100))+ylim(c(0,100))+
    geom_text(data = datos_f,aes(x = 55, y= 85, label = sprintf("Area = %0.4f",area))
              , color= "black", hjust= 0, size = 5) +
    theme(axis.line.x = element_line(color="black", size = 0.5),
          axis.line.y = element_line(color="black", size = 0.5))+
    theme(panel.border = element_blank(),
          legend.key = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.major.y = element_line(linetype = 2, color="ivory3"),
          panel.grid.major.x = element_blank(), #,element_line(linetype = 2, color="slategray"),
          legend.title = element_text(size=9, face="bold"),
          legend.text = element_text(size=9, face="bold"),
          panel.spacing = unit(22, "mm"),
          plot.title = element_text(hjust = 1.2),
          axis.text = element_text(face="bold", size=12),
          axis.title.x = element_text(face="bold", size=12),
          axis.title.y  = element_text(face="bold", size=12))
  return(ls)
  
}

languageEl <- "EN"
crit <- "MusRank"
metodo <- "PLdestroymethod"
if (metodo == "GCdestroymethod") {
  if (languageEl == "EN"){
    ytext <- "Remaining Giant Component (%)\n"
  } else {
    ytext <- "Componente gigante que queda (%)\n"
  }

} else  {
  if (languageEl == "ES"){
    ytext <- "Plantas sobrevivientes (%)\n"
  } else {
    ytext <- "Surviving plants (%)\n "
  }
}

if (languageEl == "ES"){
  xtext <- "\nAnimales eliminados (%) segun "
} else {
  xtext <- "\nRemoved animal species (%) by "
}

dfGC <- read.csv("python/results/mean_all_GCdestroymethod.csv")
names(dfGC) <- c("X","Network","area_NoOrder","area_MR","area_Krisk","area_Kdegree","area_Degree","area_eigenc")
dfPL <- read.csv("python/results/mean_all_PLdestroymethod.csv")
names(dfPL) <- c("X","Network","area_NoOrder","area_MR","area_Krisk","area_Kdegree","area_Degree","area_eigenc")
dir.create(paste0("graphs/"), showWarnings = FALSE)
dir.create(paste0("graphs/AREAS/"), showWarnings = FALSE)

listaredes <- c("M_PL_002")
listaredes <- paste0("M_",dfGC$Network)  
for (red in listaredes)
{
  if (metodo == "GCdestroymethod") {
    if (languageEl == "EN"){
      ytext <- "Remaining Giant Component (%)\n"
    } else {
      ytext <- "Componente gigante que queda (%)\n"
    }

  } else  {
    if (languageEl == "ES"){
    ytext <- "Plantas sobrevivientes (%)\n"
    } else {
    ytext <- "Surviving plants (%)\n "
    }
  }


  areaPLmr <- dfPL[paste0("M_",dfPL$Network) == red,]$area_MR
  areaPLkd <- dfPL[paste0("M_",dfPL$Network) == red,]$area_Kdegree
  areaGCmr <- dfGC[paste0("M_",dfGC$Network) == red,]$area_MR
  areaGCkd <- dfGC[paste0("M_",dfGC$Network) == red,]$area_Kdegree
  
  metodo <- "GCdestroymethod"
  crit <- "MusRank"
  datos_criterio <- read.table(paste0("python/results/",metodo,"/",red,"_Diam_extin_",crit,".txt"), quote="\"", comment.char="")
  names(datos_criterio) <- c("Primary","RemainingGC")
  p <- pintacurva(datos_criterio,areaGCmr,"pink","Surviving GC (%)")
  crit <- "Kdegree"
  datos_criterio <- read.table(paste0("python/results/",metodo,"/",red,"_Diam_extin_",crit,".txt"), quote="\"", comment.char="")
  names(datos_criterio) <- c("Primary","RemainingGC")
  q <- pintacurva(datos_criterio,areaGCkd,"lightblue","Surviving GC (%)")
  metodo <- "PLdestroymethod"
  crit <- "MusRank"
  datos_criterio <- read.table(paste0("python/results/",metodo,"/",red,"_Diam_extin_",crit,".txt"), quote="\"", comment.char="")
  names(datos_criterio) <- c("Primary","RemainingGC")
  r <- pintacurva(datos_criterio,areaPLmr,"pink","Surviving plant species (%)")
  crit <- "Kdegree"
  datos_criterio <- read.table(paste0("python/results/",metodo,"/",red,"_Diam_extin_",crit,".txt"), quote="\"", comment.char="")
  names(datos_criterio) <- c("Primary","RemainingGC")
  s <- pintacurva(datos_criterio,areaPLkd,"lightblue","Surviving plant species (%)")

  ppi <- 300
  png(paste0("graphs/AREAS/",red,"_extinction_plot.png"), width=(12*ppi), height=10*ppi, res=ppi)
  grid.arrange(p,q,r,s,nrow=2,ncol=2,
               top = textGrob(paste("Network",red), vjust = 1, gp = gpar(fontface = "bold", cex = 1.5)))
  dev.off()
}