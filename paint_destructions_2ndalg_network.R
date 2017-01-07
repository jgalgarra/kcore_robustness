library(ggplot2)
library(dplyr)

methodlist <- c("GCdestroymethod","PLdestroymethod")
methodstr <- c("Measuring remaining GC", "Measuring surviving plant species")
names(methodstr) <- c("GCdestroymethod","PLdestroymethod")

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
  return(dfparcsal)
}


for (method in (methodlist))
{
  directorystr <- paste0("python/results/",method)

  if (exists("dftot"))
    rm("dftot")
  ficheros <- Sys.glob(paste0(directorystr,"/DIAM_EXTIN_ALL_",method,"*.txt"))

  for (j in ficheros)
  {
    dfparc <- read.csv(j,sep="\t")
    if (exists("dftot"))
      dftot <- rbind(dftot,dfparc)
    else
      dftot <- dfparc
  }
  dftot$Network = paste0("M_",dftot$X,".csv")
  networklist <- unique(dftot$Network)
  #networklist <- c("M_PL_005.csv")
  size_networks <- read.csv("extinctions/size_networks.csv")

  for (network in networklist)
  {
    dfsal <- dftot[dftot$Network == network,]
    gcsize <- size_networks[size_networks$Network == network,]$giant_component
    dfindiv <- data.frame(index = c(), removed =c())
    for (k in 1:nrow(dfsal)){
      dfip <- data.frame( index = "NoOrder", removed = dfsal$NoOrder[k])
      dfindiv <- rbind(dfindiv,dfip)
      dfip <- data.frame( index = "MusRank", removed = dfsal$MR[k])
      dfindiv <- rbind(dfindiv,dfip)
      dfip <- data.frame( index = "krisk", removed = dfsal$Krisk[k])
      dfindiv <- rbind(dfindiv,dfip)
      dfip <- data.frame( index = "degree", removed = dfsal$Degree[k])
      dfindiv <- rbind(dfindiv,dfip)
      dfip <- data.frame( index = "kdegree", removed = dfsal$Kdegree[k])
      dfindiv <- rbind(dfindiv,dfip)
      dfip <- data.frame( index = "eigenc", removed = dfsal$eigenc[k])
      dfindiv <- rbind(dfindiv,dfip)
    }

    dfindiv$index <- as.factor(dfindiv$index)
    dfindiv$removed <- 100*dfindiv$removed/gcsize

    dfindiv <- filter(dfindiv, index!="NoOrder")
    o <- dfindiv %>% group_by(index,removed) %>% summarise(total=n())
    p <- dfindiv %>% group_by(index) %>% summarise(media=mean(removed))
    p$pshape = 20
    p[p$media == min(p$media),]$pshape = 19
    p$pcolor = "grey15"
    p[p$media == min(p$media),]$pcolor = "black"
    dfindiv$width <- 0
    dfindiv$width <- 0
    recorrido <- (max(dfindiv$removed)-min(dfindiv$removed))
    for (k in (1:nrow(dfindiv))){
      dfindiv$width[k] <- o[(o$index == dfindiv$index[k]) & (o$removed == dfindiv$removed[k]),]$total/400
      dfindiv$height[k] <- 0.00005*recorrido/dfindiv$width[k]
    }
    redname <- strsplit(network,".csv")[[1]][1]

    gr <- ggplot(dfindiv)+
      geom_jitter(aes(x=index,y=removed,fill=index,color=index),width=dfindiv$width,height = dfindiv$height,alpha=0.15)+
      geom_hline(yintercept=min(p$media),color = "bisque2", size=0.5, alpha = 0.6) +
      geom_point(data=p,aes(index, y=media),size=3,shape=p$pshape,fill=p$pcolor,color=p$pcolor) +
      theme_bw()+
      xlab("")+ylab(paste0("AUC ",methodstr[method],"\n"))+
      ggtitle(sprintf("%s best mean performance: %0.5f",redname,min(p$media))) +
      theme(legend.position = "none",
            axis.text.x = element_text(face="bold", color="grey30", size=9),
            axis.text.y = element_text(face="bold", color="grey30", size=9),
            panel.grid.minor = element_blank(),
            panel.grid.major.y = element_line(linetype = 3, color="ivory3", size=0.25),
            panel.grid.major.x = element_blank(),
            plot.title = element_text(hjust = 0.5))


    dir.create(paste0("graphs/python/"), showWarnings = FALSE)
    dir.create(paste0("graphs/python/",method,"/"), showWarnings = FALSE)

    ppi <- 300
    png(paste0("graphs/python/",method,"/",redname,".png"), width=(6*ppi), height=3*ppi, res=ppi)
    print(gr)
    dev.off()
  }
}
