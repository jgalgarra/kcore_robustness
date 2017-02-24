# Creates the correlation plot of degree vs kdegree and cumulative distribution plots
# Output data at "graphs/"
#
# Requires: "datos_analisis_condegs.RData"


library(grid)
library(gridExtra)
library(stargazer)
library(igraph)
library(ggplot2)
library(kcorebip)


gen_deg_distribution <- function(red, seq_breaks = seq(1,26, by=5))
{
  result_analysis <- analyze_network(red, directory = "data/", guild_a = "Plant", 
                                     guild_b = "Pollinator", plot_graphs = FALSE)
  ddegree <- igraph::degree(result_analysis$graph,mode = c("out"), loops = TRUE, normalized = FALSE)
  
  occur <- sort(as.numeric(ddegree))
  alpha_level = 0.5
  p = occur/sum(occur)
  dy = rev(cumsum(rev(p)))
  dx = occur
  
  kdegree <- V(result_analysis$graph)$kdegree
  occur <- sort(as.numeric(kdegree))
  p = occur/sum(occur)
  ky = rev(cumsum(rev(p)))
  kx = occur
  
  alfa <- lm(ddegree ~ kdegree)$coefficients[2]
  intercept <- lm(ddegree ~ kdegree)$coefficients[1]
  
  
  red_name <- strsplit(red,".csv")[[1]][1]
  
  auxdkd <- data.frame(degree = ddegree, kdegree = kdegree)
  modelt <- lm(kdegree ~ degree, data = auxdkd)
  fitted_modelt <- data.frame(
    degree = auxdkd$degree, kdegree = auxdkd$kdegree,
    predict(modelt, interval = "confidence")
  )
  
  layer_linet <- geom_line(
    mapping = aes(x =degree, y = fit),
    data = fitted_modelt,
    color = "violetred1",
    alpha = 0.5
  )
  
  textoreg1t <- sprintf("kdegree = %0.3f degree + %0.3f\n",summary(modelt)$coefficients[2, 1] ,summary(modelt)$coefficients[1, 1])
  textoreg2t <- sprintf("R^2 = %0.3f",summary(modelt)$r.squared)
  
  scatterdegkdeg <- ggplot(data = auxdkd, aes(x = degree, y = kdegree)) + geom_point(color = "blue",fill="blue", shape=21)+
    ggtitle(red_name) + ylab(expression(paste(bar(k)["degree"],"\n"))) +
    geom_text(data = auxdkd,aes(x = 0.5, y= 12, 
                                label = paste0(textoreg1t,textoreg2t)
    ), color= "black", hjust= 0, size = 3.5) +  
    layer_linet +
    theme_bw() +
    theme(
      panel.grid.minor = element_blank(),
      axis.title.x = element_text(color="grey30", size=13),
      axis.title.y = element_text(color="grey30", size=13),
      axis.text.x = element_text(face="bold", color="grey30", size=10),
      axis.text.y = element_text(face="bold", color="grey30", size=10)
    )
  
  auxdf <- data.frame(dx,dy,kx,ky)

  auxdf$dxend <- c(auxdf$dx[2:nrow(auxdf)], NA)
  auxdf$dyend <- auxdf$dy
  auxdf$kxend <- c(auxdf$kx[2:nrow(auxdf)], NA)
  auxdf$kyend <- auxdf$ky
  pointdfd <- auxdf[(auxdf$dx != auxdf$dxend) & !is.na(auxdf$dxend),]
  pointdfd <- rbind(pointdfd, data.frame(dx=1,dy=1,kx=1,ky=1,dxend=1,dyend=1,kxend=1,kyend=1))
  pointdfk <- auxdf[(auxdf$kx != auxdf$kxend) & !is.na(auxdf$kxend),]

  yseq_breaks <- c(0.1,0.5,1.0)
  dseq_breaks <- c()
  br <- unique(auxdf$dx)
  for (i in br)
    if ((i<6) | (!(i%%2) | (i == max(br))))
      dseq_breaks <- c(dseq_breaks,i)
  
  dist_deg <- ggplot(auxdf, aes(x=dx, y=dy, xend=dxend, yend=dyend)) +
          geom_point(data = pointdfd, aes(x=dxend, y=dy),  size = 0.7, color = "red", alpha= 0.5) +  # Open points to right
          geom_segment(color = "red", size=0.5, alpha= 0.5)+  # Horizontal line segments
    scale_x_log10(breaks = dseq_breaks) + scale_y_log10(breaks = yseq_breaks) + xlab("degree") + ylab(cumulativetxt) + ggtitle("") +
    theme_bw() + 
    theme(
      panel.grid.minor.x = element_blank(),
      panel.grid.minor.y = element_blank(),
      axis.title.x = element_text(color="grey30", size=13),
      axis.title.y = element_text(color="grey30", size=13),
      axis.text.x = element_text(face="bold", color="grey30", size=10),
      axis.text.y = element_text(face="bold", color="grey30", size=10)
    )
  
  kseq_breaks <- c()
  kseq_breaks <- seq(min(auxdf$kx),max(auxdf$kx), by=2)
  kseq_breaks <- unique(c(1,kseq_breaks))
  
  dist_kdeg <- ggplot(auxdf, aes(x=kx, y=ky, xend=kxend, yend=kyend)) +
    geom_point(data = pointdfk, aes(x=kxend, y=ky),  size = 0.7, color = "palegreen4", alpha = 0.5) +  # Open points to right
    geom_segment(color = "palegreen4", size=0.5, alpha= 0.5)+  # Horizontal line segments
    scale_x_log10(breaks = kseq_breaks) + scale_y_log10(breaks = yseq_breaks) + xlab("kdegree") + 
    ylab(cumulativetxt) + ggtitle("") + ylab("")+
    theme_bw() + 
    theme(
      panel.grid.minor.x = element_blank(),
      panel.grid.minor.y = element_blank(),
      axis.title.x = element_text(color="grey30", size=13),
      axis.title.y = element_text(color="grey30", size=13),
      axis.text.x = element_text(face="bold", color="grey30", size=10),
      axis.text.y = element_text(face="bold", color="grey30", size=10),
      plot.title = element_text(size=16,lineheight=.8, face="bold",hjust = 0.5)
    )
  
  

  
  dist_dkdeg <- ggplot(auxdf, aes(x=kx*alfa, y=ky, xend=kxend*alfa, yend=ky)) +
    geom_point(data = pointdfk, aes(x=kxend*alfa, y=ky),  size = 0.7, color = "palegreen4", alpha = 0.5) +  
    geom_point(data = pointdfd, aes(x=dxend, y=dy),  size = 0.7, color = "red", alpha= 0.5) +
    geom_segment(color = "palegreen4", size=0.5, alpha= 0.5)+  # Horizontal line segments
    geom_segment(aes(x=dx,y=dy,xend=dxend,yend=dyend),color = "red", size=0.5, alpha= 0.5)+
    scale_x_log10(breaks = dseq_breaks) + scale_y_log10(breaks = yseq_breaks) + xlab("degree scale") + ylab(cumulativetxt) + ggtitle("") + ylab("")+
    theme_bw() + 
    theme(
      panel.grid.minor.x = element_blank(),
      panel.grid.minor.y = element_blank(),
      axis.title.x = element_text(color="grey30", size=13),
      axis.title.y = element_text(color="grey30", size=13),
      axis.text.x = element_text(face="bold", color="grey30", size=10),
      axis.text.y = element_text(face="bold", color="grey30", size=10)
    )
  

  
  calc_values <- list("dist_deg" = dist_deg, "dist_kdeg" = dist_kdeg, "dist_dkdeg" = dist_dkdeg, 
                      "scatterdegkdeg" = scatterdegkdeg)
  return(calc_values)
}

red <- "M_PL_001"
languageEl <- "EN"
if (languageEl == "EN"){
  xcorrtxt = "Correlation degree kdegree"
  ycorrtxt = "Number of networks"
  cumulativetxt = "Cumulative distribution"
  xscale = "degree scale"
} else {
  xcorrtxt = "Correlación degree kdegree"
  ycorrtxt = "Número de redes"
  cumulativetxt = "Distribución acumulada"
  xscale = "escala degree"
}


grafs <- gen_deg_distribution(paste0(red,".csv"))

ppi <- 300
png(paste0("graphs/kdegree_over_degree_",red,"_",languageEl,".png"), width=(12*ppi), height=4*ppi, res=ppi)
grid.arrange(grafs$dist_deg,grafs$dist_kdeg,grafs$dist_dkdeg, ncol=3, nrow=1, widths=c(1/3,1/3,1/3) )
dev.off()


load("results/datos_analisis_condegs.RData")
distr2 <- ggplot(data = resultdf, aes(x = sqrt(kdegdegRsq) ))+ 
  geom_histogram(  fill = "palegreen4", 
                   color = "white",  alpha = 0.7) + 
  xlab(xcorrtxt)+ylab(ycorrtxt) +
  theme_bw() +
  theme(
    axis.title.x = element_text(color="grey30", size=13),
    axis.title.y = element_text(color="grey30", size=13),
    axis.text.x = element_text(face="bold", color="grey30", size=10),
    axis.text.y = element_text(face="bold", color="grey30", size=10)
  )

ppi <- 300
png(paste0("graphs/kdegree_vs_degree_",red,"_",languageEl,".png"), width=(12*ppi), height=3*ppi, res=ppi)
grid.arrange(distr2,grafs$scatterdegkdeg, ncol=2, nrow=1, widths=c(1/2,1/2) )
dev.off()

ppi <- 300
png(paste0("graphs/ALL_plots_kdegree_degree_",red,"_",languageEl,".png"), width=(12*ppi), height=8*ppi, res=ppi)
grid.arrange(arrangeGrob(distr2,grafs$scatterdegkdeg, ncol=2, nrow=1, widths=c(1/2,1/2) ), 
             nrow=2, heights=c(0.5,0.5),
             arrangeGrob(grafs$dist_deg,grafs$dist_kdeg,grafs$dist_dkdeg, ncol=3, nrow=1, widths=c(1/3,1/3,1/3) ))
dev.off()