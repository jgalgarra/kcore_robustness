# Auxiliar script to create a bipartite graph in screen

library(kcorebip)
library(igraph)

giant.component <- function(graph) {
  cl <- clusters(graph)
  induced.subgraph(graph, which(cl$membership == which.max(cl$csize)))}

plot_bipartite <- function(bg, aspect_ratio = 9/35, vframecolor = "grey70", vlabelcex = 4,
                           vsize = 4, vcolor = c("cadetblue1","pink2"), labelcolor = c("blue","red"),
                           framedisp = FALSE,  color_link = "grey50", vertical = FALSE, nname = "")
{
  l <- layout.bipartite(bg)
  if (vertical)
    la <- l[, c(2,1)]
  else
    la <- l
  png(paste0("datawip/",nname,".png"), height=600, width=1800)
  #functions to plot your graph

  plot.igraph(bg, layout= la,asp=aspect_ratio,vertex.frame.color=vframecolor,
              vertex.label.cex=vlabelcex,vertex.label.color="black",vertex.label.family="Arial",
              vertex.size=vsize, edge.color= color_link, frame = framedisp,
              vertex.color=vcolor[V(bg)$type+1])
  dev.off()

}

quita_vacios <- function(grafo,lgrados)
{
  for (k in 1:length(lgrados))
  {
    if (lgrados[k] == 0){
      grafo <- delete_vertices(grafo,names(lgrados[k]))
      if ( length(grep("Pollinator",names(lgrados[k]) )) >0 )
        num_species_b <<- num_species_b -1
      else
        num_species_a <<- num_species_a -  1
    }
  }
  return(grafo)
}

redname <- "M_PL_012"
num_species_a <- 0
num_species_b <- 0


result_analysis <- analyze_network(paste0(redname,".csv"), directory = "data/", guild_a = "Plant", guild_b = "Pollinator", plot_graphs = FALSE)
original_GC <- length(V(giant.component(result_analysis$graph)))
print(paste("Original Network",redname))
print(paste("Giant Component",original_GC))
print(paste("Animal species",result_analysis$num_guild_b,"Vegetal species",result_analysis$num_guild_a))
bp <- get_bipartite(result_analysis$graph, plot_graphs = FALSE)
plot_bipartite(bp, aspect_ratio = 1/6,vlabelcex=1.5,vsize = 3.5, vframecolor = "transparent",
               color_link = "grey20", vertical = FALSE, nname = redname)
