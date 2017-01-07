# Auxiliar script to create a bipartite graph in screen

library(kcorebip)

plot_bipartite <- function(bg, aspect_ratio = 9/35, vframecolor = "grey70", vlabelcex = 4,
                           vsize = 4, vcolor = c("lightblue","pink2"), labelcolor = c("blue","red"),
                           framedisp = FALSE,  color_link = "grey50", vertical = FALSE, nname = "")
{
  l <- layout.bipartite(bg)
  if (vertical)
    la <- l[, c(2,1)]
  else
    la <- l
  png(paste0("datawip/",nname,".png"), heigh=600, width=1800)
  #functions to plot your graph

  plot.igraph(bg, layout= la,asp=aspect_ratio,vertex.frame.color=vframecolor,
              vertex.label.cex=vlabelcex,vertex.label.color="black",
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

redname <- "M_PL_007"
num_species_a <- 0
num_species_b <- 0


result_analysis <- analyze_network(paste0(redname,".csv"), directory = "data/", guild_a = "Plant", guild_b = "Pollinator", plot_graphs = FALSE)
original_GC <- length(V(giant.component(result_analysis$graph)))
print(paste("Original Network",redname))
print(paste("Giant Component",original_GC))
print(paste("Animal species",result_analysis$num_guild_b,"Vegetal species",result_analysis$num_guild_a))
bp <- get_bipartite(result_analysis$graph, plot_graphs = FALSE)
plot_bipartite(bp, aspect_ratio = 1/6,vlabelcex=2.0,vsize = 5, vframecolor = "grey20",
               color_link = "grey20", vertical = FALSE, nname = redname)

redname <- "M_PL_007_kd"
result_analysis <- analyze_network(paste0(redname,".csv"), directory = "datawip/", guild_a = "Plant", guild_b = "Pollinator", plot_graphs = FALSE)
kd_GC <- length(V(giant.component(result_analysis$graph)))
result_analysis$degree <- igraph::degree(result_analysis$graph)
num_species_a <- result_analysis$num_guild_a
num_species_b <- result_analysis$num_guild_b
result_analysis$graph <- quita_vacios(result_analysis$graph, result_analysis$degree )
print("Network destroyed following kdegree")
print(paste("Giant Component",kd_GC))
print(paste("Animal species",num_species_b,"Vegetal species",num_species_a))
bp <- get_bipartite(result_analysis$graph, plot_graphs = FALSE)
plot_bipartite(bp, aspect_ratio = 1/6,vlabelcex=2.0,vsize = 5, vframecolor = "grey20",
               color_link = "grey20", vertical = FALSE, nname = redname)

redname <- "M_PL_007_mr"
result_analysis <- analyze_network(paste0(redname,".csv"), directory = "datawip/", guild_a = "Plant", guild_b = "Pollinator", plot_graphs = FALSE)
result_analysis$degree <- igraph::degree(result_analysis$graph)
mr_GC <- length(V(giant.component(result_analysis$graph)))
num_species_a <- result_analysis$num_guild_a
num_species_b <- result_analysis$num_guild_b
result_analysis$graph <- quita_vacios(result_analysis$graph, result_analysis$degree )
print("Network destroyed following MusRank")
print(paste("Giant Component",mr_GC))
print(paste("Animal species",num_species_b,"Vegetal species",num_species_a))
bp <- get_bipartite(result_analysis$graph, plot_graphs = FALSE)
plot_bipartite(bp, aspect_ratio = 1/6,vlabelcex=2.0,vsize = 5, vframecolor = "grey20",
               color_link = "grey20", vertical = FALSE, nname = redname)

# ziggurat_graph("data/","M_PL_007.csv", height_box_y_expand = 0.75,
#                lsize_legend = 7, lsize_core_box = 6,corebox_border_size=1,
#                plotsdir = "datawip/",color_link = "slategray3", alpha_link = 0.5,
#                lsize_kcoremax = 6,lsize_zig = 5,lsize_kcore1 = 5,
#                displace_legend = c(-0.2,0.2),displace_outside_component = c(-0.5,0.6),print_to_file = TRUE)
#
#
# ziggurat_graph("datawip/","M_PL_007_kd.csv", height_box_y_expand = 0.75,
#                lsize_legend = 7, lsize_core_box = 6,corebox_border_size=1,
#                plotsdir = "datawip/",color_link = "slategray3", alpha_link = 0.5,
#                lsize_kcoremax = 6,lsize_zig = 5,lsize_kcore1 = 5,
#                displace_legend = c(-0.2,0.2),displace_outside_component = c(-0.5,0.6),print_to_file = TRUE)
#
#
# ziggurat_graph("datawip/","M_PL_007_mr.csv", height_box_y_expand = 0.75,
#                lsize_legend = 7, lsize_core_box = 6,corebox_border_size=1,
#                plotsdir = "datawip/",color_link = "slategray3", alpha_link = 0.5,
#                lsize_kcoremax = 6,lsize_zig = 5,lsize_kcore1 = 5,
#                displace_legend = c(-0.2,0.2),displace_outside_component = c(-0.5,0.6),print_to_file = TRUE)
