library(ggplot2)
library(dplyr)

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


calc_stats_network <- function(dfinput, network)
{
  #dfsal <- data.frame( Network=NA, giant_component=NA, krisk=NA, degree=NA, kdegree=NA,eigenc=NA, best =NA)[numeric(0), ]
  dfparcsal <- data.frame( Network="", giant_component=0, NoOrder = 0, krisk = 0,
                           degree=0,kdegree=0, eigenc=0, best = 0)
  dfnet <- dfinput[dfinput$Network == network,]
  dfparcsal$Network <- j
  dfparcsal$giant_component <- mean(dfnet$giant_component)
  dfparcsal$NoOrder <- fn(dfnet$NoOrder)
  dfparcsal$krisk <- fn(dfnet$krisk)
  dfparcsal$kdegree <- fn(dfnet$kdegree)
  dfparcsal$degree <- fn(dfnet$degree)
  dfparcsal$eigenc <-fn(dfnet$eigenc)
  dfparcsal$best <- min(c(dfparcsal$NoOrder,dfparcsal$krisk,dfparcsal$kdegree,dfparcsal$degree,dfparcsal$eigenc))
  #dfsal <- rbind(dfsal, dfparcsal)
  #write.csv(dfsal,paste0(directorystr,"/Destructions_1st_",fnname,".csv"))
  return(dfparcsal)
}

networklist <- unique(dftot$Network)
#networklist <- c("M_PL_015.csv")

for (network in networklist)
{
  dfsal <- dftot[dftot$Network == network,]
  gcsize <- dfsal$giant_component[1]
  dfindiv <- data.frame(index = c(), removed =c())
  for (k in 1:nrow(dfsal)){
    dfip <- data.frame( index = "NoOrder", removed = dfsal$NoOrder[k])
    dfindiv <- rbind(dfindiv,dfip)
    dfip <- data.frame( index = "krisk", removed = dfsal$krisk[k])
    dfindiv <- rbind(dfindiv,dfip)
    dfip <- data.frame( index = "degree", removed = dfsal$degree[k])
    dfindiv <- rbind(dfindiv,dfip)
    dfip <- data.frame( index = "kdegree", removed = dfsal$kdegree[k])
    dfindiv <- rbind(dfindiv,dfip)
    dfip <- data.frame( index = "eigenc", removed = dfsal$eigenc[k])
    dfindiv <- rbind(dfindiv,dfip)
  }
    
  dfindiv$index <- as.factor(dfindiv$index)
  dfindiv$removed <- 100*dfindiv$removed/gcsize
  
  dfindiv <- filter(dfindiv, index!="NoOrder")
  o <- dfindiv %>% group_by(index,removed) %>% summarise(total=n())
  p <- dfindiv %>% group_by(index) %>% summarise(media=mean(removed))
  p$pshape = 18
  p[p$media == min(p$media),]$pshape = 19
  p$pcolor = "grey15"
  p[p$media == min(p$media),]$pcolor = "black"
  dfindiv$width <- 0
  for (k in (1:nrow(dfindiv)))
    dfindiv$width[k] <- o[(o$index == dfindiv$index[k]) & (o$removed == dfindiv$removed[k]),]$total/300
  redname <- strsplit(network,".csv")[[1]][1]

  recorrido <- (max(dfindiv$removed)-min(dfindiv$removed))
  gr <- ggplot(dfindiv)+
    geom_hline(yintercept=min(p$media),color = "bisque2", size=0.5, alpha = 0.6) +
    geom_jitter(aes(x=index,y=removed,fill=factor(index),color=factor(index)),width=dfindiv$width,height=recorrido/1000,alpha=0.15)+
    geom_point(data=p,aes(index, y=media),size=3,shape=p$pshape,fill=p$pcolor,color=p$pcolor) +
    
    theme_bw()+
    xlab("")+ylab("Removed species (%)\n")+
    ggtitle(sprintf("%s best mean performance: %0.3f",redname,min(p$media))) +
    theme(legend.position = "none",
          axis.text.x = element_text(face="bold", color="grey30", size=12),
          axis.text.y = element_text(face="bold", color="grey30", size=10),
          panel.grid.minor = element_blank(),
          panel.grid.major.y = element_line(linetype = 3, color="ivory3", size=0.25),
          panel.grid.major.x = element_blank(),
          plot.title = element_text(hjust = 0.5))
  
  
  
  dir.create("graphs/FIRST", showWarnings = FALSE)

  ppi <- 300
  png(paste0("graphs/FIRST/First_performance_",redname,".png"), width=(6*ppi), height=2.5*ppi, res=ppi)
  print(gr)
  dev.off()
}