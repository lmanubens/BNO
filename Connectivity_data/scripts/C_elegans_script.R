
library("igraph")
library("plyr")
library(leiden)
library(gmp)
library(Matrix)



#Preparacio dades C.Elegans

# documentacio matriu https://www.wormatlas.org/neuronalwiring.html#NeuronalconnectivityII
# documentacio igraph https://kateto.net/netscix2016.html

archiu <- read.csv("NeuronConnect.csv",sep = ";")

unique(archiu$Type)


ab <- nrow(archiu)
for (i in 1:ab) {
  if(archiu[i,3]=="R" | archiu[i,3]=="Rp"){
    archiu[i,] <- archiu[i, c(2,1,3,4)]
  }
  if(archiu[i,3]=="EJ"){
    
    archiu <- rbind.fill(archiu,archiu[i,])
  }
}

g <- graph_from_edgelist(as.matrix(archiu[,c(1,2)]), directed = TRUE)

E(g)$weight <- archiu$Nbr
E(g)$name <- archiu$Type


gs <- simplify( g , remove.multiple = T, remove.loops = F, edge.attr.comb=c(weight="sum", type="ignore") )


# l <- layout_with_fr(gs)
# plot(gs, rescale=F, layout=l*0.4)
# 
# l <- layout_with_fr(gs)
# 
# V(gs)$color <- colrs[V(gs)$media.type]
# 
# V(gs)$size <- V(gs)$audience.size*0.7
# 
# V(gs)$label.color <- "black"
# 
# V(gs)$label <- NA
# 
# E(gs)$width <- E(gs)$weight/6
# 
# E(gs)$arrow.size <- .2
# 
# E(gs)$width <- 1+E(gs)$weight/12
# 
# plot(gs,rescale=F,layout=l*0.8)
# 
# barplot(triad_census(gs), log="y")
# 
# #Color nodes along the diameter:
# 
# diam <- get_diameter(gs, directed=T)
# 
# as.vector(diam)
# 
# vcol <- rep("gray40", vcount(gs))
# 
# vcol[diam] <- "gold"
# 
# ecol <- rep("gray80", ecount(gs))
# 
# ecol[E(gs, path=diam)] <- "orange" 
# 
# 
# plot(gs, rescale=F,vertex.color=vcol, edge.color=ecol, edge.arrow.mode=0)
# 
# #encontrar el camino m?s corto entre nodos espec?ficos.
# 
# news.path <- shortest_paths(gs, 
#                             
#                             from = V(gs)[as.list(names)=="ADFL"], 
#                             
#                             to  = V(gs)[as.list(names)=="SMBVR"],
#                             
#                             output = "both") # both path nodes and edges
# 
# 
# 
# # Generate edge color variable to plot the path:
# 
# ecol <- rep("gray80", ecount(gs))
# 
# ecol[unlist(news.path$epath)] <- "orange"
# 
# # Generate edge width variable to plot the path:
# 
# ew <- rep(2, ecount(gs))
# 
# ew[unlist(news.path$epath)] <- 4
# 
# # Generate node color variable to plot the path:
# 
# vcol <- rep("gray40", vcount(gs))
# 
# vcol[unlist(news.path$vpath)] <- "gold"
# 
# 
# 
# plot(gs,rescale=F, vertex.color=vcol, edge.color=ecol, 
#      
#      edge.width=ew, edge.arrow.mode=0)
# 
# inc.edges <- incident(gs,  V(gs)[media=="ADFL"], mode="all")
# 
# 
# 
# # Set colors to plot the selected edges.
# 
# ecol <- rep("gray80", ecount(gs))
# 
# ecol[inc.edges] <- "orange"
# 
# vcol <- rep("grey40", vcount(gs))
# 
# vcol[V(gs)$media=="ADFL"] <- "gold"
# 
# plot(gs,rescale=F, vertex.color=vcol, edge.color=ecol)
# 
# cocitation(gs)
# 
# #comunidades
# ceb <- cluster_edge_betweenness(gs) 
# 
# plot(ceb, gs) 
# 
# cfg <- cluster_fast_greedy(as.undirected(gs))
# 
# plot(cfg, as.undirected(gs))
# 
# l <- layout_with_fr(gs)
# 
# plot(cfg, resacle=F, as.undirected(gs), layout=l*1.0)






#----------------------------------------------------------C_elegans--------------------------------------------------------#



deg <- degree(gs, mode="all")

deg.distF <- degree_distribution(gs, cumulative=F, mode="all")

deg.distT <- degree_distribution(gs, cumulative=T, mode="all")



#Triad census
tirad_gs = triad_census(gs)

#Hubs
hs_gs <- hub_score(gs, weights=NA)$vector

#Authorities
as_gs <- authority_score(gs, weights=NA)$vector


#Small-Worldness ratio: 
SM_Ratio_gs = transitivity(gs)/mean_distance(gs)

#modularity y comunitat
adj_mtx_gs  <- igraph::as_adjacency_matrix(gs)
comunitats_gs <- leiden(adj_mtx_gs)
modularity_gs <- modularity(gs, comunitats_gs)


#cliques
cliques_gs <- cliques(gs)

#n nodes
nEdges_gs <- gsize(gs)
cost_gs <- sum(E(gs)$weight)
nNodes_gs <- gorder(gs)


#routing efficiency
am_gs  <- get.adjacency(gs)
distances_gs <- distances(gs , mode = "in", weights=1/E(gs)$weight)


rout_effs_gs <- as.numeric(as.bigz(1/distances_gs)/as.bigz(dim(am_gs)[[1]]*(dim(am_gs)[[1]]-1)))
rout_effs_matrix_gs <- as.matrix(rout_effs_gs)
rout_effs_matrix_gs <- rout_effs_matrix_gs[is.finite(rout_effs_matrix_gs)]
rout_eff_gs <- as.numeric(sum(rout_effs_matrix_gs))

#storage capacity
storage_caps_gs <- numeric()
for(j in 1:dim(am_gs)[[1]]){
  storage_caps_gs <- append(storage_caps_gs, as.integer(2*log2(chooseZ(as.bigz(sum(am_gs[,j]) + nnzero(am_gs[,j]) -1), as.bigz(sum(am_gs[,j]))))))
}
storage_cap_gs <- mean(storage_caps_gs)#,na.rm=T)

celegans_data <- list()
celegans_data[["deg_distF_celegans"]] <- deg.distF
celegans_data[["deg_distT_celegans"]] <- deg.distT
celegans_data[["triad_census_celegans"]] <- tirad_gs
celegans_data[["hub_score_celegans"]] <- hs_gs
celegans_data[["modularity_celegans"]] <- modularity_gs
celegans_data[["authority_score_celegans"]] <- as_gs 
celegans_data[["communities_celegans"]] <- comunitats_gs
celegans_data[["SW_ratio_celegans"]] <- SM_Ratio_gs 
celegans_data[["cliques_celegans"]] <- cliques_gs 
celegans_data[["cost_celegans"]] <- cost_gs
celegans_data[["nEdges_celegans"]] <- nEdges_gs
celegans_data[["nNodes_celegans"]] <- nNodes_gs 
celegans_data[["rout_eff_celegans"]] <- rout_eff_gs 
celegans_data[["storage_cap_celegans"]] <- storage_cap_gs




save(celegans_data, file="celegans_data.RData")

load("celegans_data.RData")



#------------------------------------------------------ Erdos-Renyi random graph model -------------------------------------#

set.seed(2706)

er_xarxa1 <- sample_gnm(n=283, m=3109, directed = TRUE) 
V(er_xarxa1)$name <- V(gs)$name
E(er_xarxa1)$weight <- E(gs)$weight


deg_xarxa1 <- degree(er_xarxa1, mode="all")

deg.dist_xarxaF1 <- degree_distribution(er_xarxa1, cumulative=F, mode="all")

deg.dist_xarxaT1 <- degree_distribution(er_xarxa1, cumulative=T, mode="all")


tirad_er_xarxa1 = triad_census(er_xarxa1)


#Hubs
hs_er_xarxa1 <- hub_score(er_xarxa1, weights=NA)$vector


#Authorities
as_er_xarxa1 <- authority_score(er_xarxa1, weights=NA)$vector


#Small-Worldness ratio: 
SM_Ratio_er_xarxa1 = transitivity(er_xarxa1)/mean_distance(er_xarxa1)


#modularity y comunitats
adj_mtx_er_xarxa1  <- igraph::as_adjacency_matrix(er_xarxa1)
comunitats_er_xarxa1 <- leiden(adj_mtx_er_xarxa1)
modularity_er_xarxa1 <- modularity(er_xarxa1, comunitats_er_xarxa1)


#cliques
cliques_er_xarxa1 <- cliques(er_xarxa1)

#n nodes
nEdges_er_xarxa1 <- gsize(er_xarxa1)
cost_er_xarxa1 <- sum(E(er_xarxa1)$weight)
nNodes_er_xarxa1 <- gorder(er_xarxa1)

#routing efficiency
am_xarxa1  <- igraph::get.adjacency(er_xarxa1)
distances_xarxa1 <- distances(er_xarxa1 , mode = "in", weights=1/E(er_xarxa1)$weight)

rout_effs_xarxa1 <- as.numeric(as.bigz(1/distances_xarxa1)/as.bigz(dim(am_xarxa1)[[1]]*(dim(am_xarxa1)[[1]]-1)))
rout_effs_matrix_xarxa1 <- as.matrix(rout_effs_xarxa1)
rout_effs_matrix_xarxa1 <- rout_effs_matrix_xarxa1[is.finite(rout_effs_matrix_xarxa1)]
rout_eff_xarxa1 <- as.numeric(sum(rout_effs_matrix_xarxa1))



#storage capacity
storage_caps_xarxa1 <- numeric()
for(j in 1:dim(am_xarxa1)[[1]]){
  storage_caps_xarxa1 <- append(storage_caps_xarxa1, as.integer(2*log2(chooseZ(as.bigz(sum(am_xarxa1[,j]) + nnzero(am_xarxa1[,j]) -1), as.bigz(sum(am_xarxa1[,j]))))))
}
storage_cap_xarxa1 <- mean(storage_caps_xarxa1)#,na.rm=T)

celegansER_data <- list()
celegansER_data[["deg_distF_celegansER"]] <- deg.dist_xarxaF1
celegansER_data[["deg_distT_celegansER"]] <- deg.dist_xarxaT1
celegansER_data[["triad_census_celegansER"]] <- tirad_er_xarxa1
celegansER_data[["hub_score_celegansER"]] <- hs_er_xarxa1
celegansER_data[["modularity_celegansER"]] <- modularity_er_xarxa1
celegansER_data[["authority_score_celegansER"]] <- as_er_xarxa1
celegansER_data[["communities_celegansER"]] <- comunitats_er_xarxa1
celegansER_data[["SW_ratio_celegansER"]] <- SM_Ratio_er_xarxa1 
celegansER_data[["cliques_celegansER"]] <- cliques_er_xarxa1 
celegansER_data[["cost_celegansER"]] <- cost_er_xarxa1 
celegansER_data[["nEdges_celegansER"]] <- nEdges_er_xarxa1
celegansER_data[["nNodes_celegansER"]] <- nNodes_er_xarxa1 
celegansER_data[["rout_eff_celegansER"]] <- rout_eff_xarxa1
celegansER_data[["storage_celegansER"]] <- storage_cap_xarxa1


save(celegansER_data, file="celegansER_data.RData")

load("celegansER_data.RData")




#------------------------------------------------ Watts-Strogatz small-world model--------------------------------------------#




sw_xarxa2 <- rewire(graph.lattice(283, nei=11, dir=TRUE, mutual=TRUE), 
                    with = each_edge(p = .1, loops = FALSE))
# sw_xarxa2 <- sample_smallworld(dim=1, size=283, nei=11, p=0.1)
# V(sw_xarxa2)$name <- V(gs)$name
# obtenim pesos aleatoris per les connexions de mes que tenim en sw_xarxa2 en comparacio amb gs
pesos_extra <- sample(E(gs)$weight,length(E(sw_xarxa2))-length(E(gs)))
E(sw_xarxa2)$weight <- c(E(gs)$weight,pesos_extra)


deg_xarxa2 <- degree(sw_xarxa2, mode="all")


deg.dist_xarxaF2 <- degree_distribution(sw_xarxa2, cumulative=F, mode="all")

deg.dist_xarxaT2 <- degree_distribution(sw_xarxa2, cumulative=T, mode="all")


tirad_sw_xarxa2 = triad_census(sw_xarxa2)

#Hubs
hs_sw_xarxa2 <- hub_score(sw_xarxa2, weights=NA)$vector

#Authorities
as_sw_xarxa2 <- authority_score(sw_xarxa2, weights=NA)$vector

#Small-Worldness ratio: 
SM_Ratio_sw_xarxa2 = transitivity(sw_xarxa2)/mean_distance(sw_xarxa2)

#modularity y comunitats
adj_mtx_sw_xarxa2  <- igraph::as_adjacency_matrix(sw_xarxa2)
comunitats_sw_xarxa2 <- leiden(adj_mtx_sw_xarxa2)
modularity_sw_xarxa2 <- modularity(sw_xarxa2, comunitats_sw_xarxa2)

#cliques
cliques_sw_xarxa2 <- cliques(sw_xarxa2)

#n nodes
nEdges_sw_xarxa2 <- gsize(sw_xarxa2)
cost_sw_xarxa2 <- sum(E(sw_xarxa2)$weight)
nNodes_sw_xarxa2 <- gorder(sw_xarxa2)

#routing efficiency
am_sw_xarxa2 <- get.adjacency(sw_xarxa2)
distances_sw_xarxa2 <- distances(sw_xarxa2 , mode = "in", weights=1/E(sw_xarxa2)$weight)

rout_effs_sw_xarxa2 <- as.numeric(as.bigz(1/distances_sw_xarxa2)/as.bigz(dim(am_sw_xarxa2)[[1]]*(dim(am_sw_xarxa2)[[1]]-1)))
rout_effs_matrix_sw_xarxa2 <- as.matrix(rout_effs_sw_xarxa2)
rout_effs_matrix_sw_xarxa2 <- rout_effs_matrix_sw_xarxa2[is.finite(rout_effs_matrix_sw_xarxa2)]
rout_eff_sw_xarxa2 <- as.numeric(sum(rout_effs_matrix_sw_xarxa2))

#storage capacity
storage_caps_sw_xarxa2 <- numeric()
for(j in 1:dim(am_sw_xarxa2)[[1]]){
  storage_caps_sw_xarxa2 <- append(storage_caps_sw_xarxa2, as.integer(2*log2(chooseZ(as.bigz(sum(am_sw_xarxa2[,j]) + nnzero(am_sw_xarxa2[,j]) -1), as.bigz(sum(am_sw_xarxa2[,j]))))))
}
storage_cap_sw_xarxa2 <- mean(storage_caps_sw_xarxa2)#,na.rm=T)

celegansSW_data <- list()
celegansSW_data[["deg_distF_celegansSW"]] <- deg.dist_xarxaF2
celegansSW_data[["deg_distT_celegansSW"]] <- deg.dist_xarxaT2
celegansSW_data[["triad_census_celegansSW"]] <- tirad_sw_xarxa2
celegansSW_data[["hub_score_celegansSW"]] <- hs_sw_xarxa2
celegansSW_data[["modularity_celegansSW"]] <- modularity_sw_xarxa2
celegansSW_data[["authority_score_celegansSW"]] <- as_sw_xarxa2 
celegansSW_data[["communities_celegansSW"]] <- comunitats_sw_xarxa2 
celegansSW_data[["SW_ratio_celegansSW"]] <- SM_Ratio_sw_xarxa2 
celegansSW_data[["cliques_celegansSW"]] <- cliques_sw_xarxa2 
celegansSW_data[["cost_celegansSW"]] <- cost_sw_xarxa2 
celegansSW_data[["nEdges_celegansSW"]] <- nEdges_sw_xarxa2
celegansSW_data[["nNodes_celegansSW"]] <- nNodes_sw_xarxa2 
celegansSW_data[["rout_eff_celegansSW"]] <- rout_eff_sw_xarxa2 
celegansSW_data[["storage_cap_celegansSW"]] <- storage_cap_sw_xarxa2


save(celegansSW_data, file="celegansSW_data.RData")

load("celegansSW_data.RData")



#----------------------------- Barabasi-Albert preferential attachment model for scale-free graphs -------------------------#


ba_xarxa3 <-  sample_pa(n=283, power=1, m=11,  directed=T)
pesos_Ba <- sample(E(gs)$weight,length(E(ba_xarxa3)))
E(ba_xarxa3)$weight <- pesos_Ba



deg_xarxa3 <- degree(ba_xarxa3, mode="all")

deg.dist_xarxaF3 <- degree_distribution(ba_xarxa3, cumulative=F, mode="all")

deg.dist_xarxaT3 <- degree_distribution(ba_xarxa3, cumulative=T, mode="all")

tirad_census_ba_xarxa3 = triad_census(ba_xarxa3)

#Hubs
hs_ba_xarxa3 <- hub_score(ba_xarxa3, weights=NA)$vector

#Authorities
as_ba_xarxa3 <- authority_score(ba_xarxa3, weights=NA)$vector

#Small-Worldness ratio: 
SM_Ratio_ba_xarxa3 = transitivity(ba_xarxa3)/mean_distance(ba_xarxa3)

#modularity y comunitats
adj_mtx_ba_xarxa3  <- igraph::as_adjacency_matrix(ba_xarxa3)
comunitats_ba_xarxa3 <- leiden(adj_mtx_ba_xarxa3)
modularity_ba_xarxa3 <- modularity(ba_xarxa3, comunitats_ba_xarxa3)


#cliques
cliques_ba_xarxa3 <- cliques(ba_xarxa3)

#n nodes
nEdges_ba_xarxa3 <- gsize(ba_xarxa3)
cost_ba_xarxa3 <- sum(E(ba_xarxa3)$weight)
nNodes_ba_xarxa3 <- gorder(ba_xarxa3)

#routing efficiency
am_ba_xarxa3 <- get.adjacency(ba_xarxa3)
distances_ba_xarxa3 <- distances(ba_xarxa3 , mode = "in", weights=1/E(ba_xarxa3)$weight)


rout_effs_ba_xarxa3 <- as.numeric(as.bigz(1/distances_ba_xarxa3)/as.bigz(dim(am_ba_xarxa3)[[1]]*(dim(am_ba_xarxa3)[[1]]-1)))
rout_effs_matrix_ba_xarxa3 <- as.matrix(rout_effs_ba_xarxa3)
rout_effs_matrix_ba_xarxa3 <- rout_effs_matrix_ba_xarxa3[is.finite(rout_effs_matrix_ba_xarxa3)]
rout_eff_ba_xarxa3 <- as.numeric(sum(rout_effs_matrix_ba_xarxa3))

#storage capacity
storage_caps_ba_xarxa3 <- numeric()
for(j in 1:dim(am_ba_xarxa3)[[1]]){
  storage_caps_ba_xarxa3 <- append(storage_caps_ba_xarxa3, as.integer(2*log2(chooseZ(as.bigz(sum(am_ba_xarxa3[,j]) + nnzero(am_ba_xarxa3[,j]) -1), as.bigz(sum(am_ba_xarxa3[,j]))))))
}
storage_cap_sw_xarxa2 <- mean(storage_caps_sw_xarxa2)#,na.rm=T)

celegansBA_data <- list()
celegansBA_data[["deg_distF_celegansBA"]] <- deg.dist_xarxaF3
celegansBA_data[["deg_distT_celegansBA"]] <- deg.dist_xarxaT3
celegansBA_data[["triad_census_celegansBA"]] <- tirad_census_ba_xarxa3
celegansBA_data[["hub_score_celegansBA"]] <- hs_ba_xarxa3
celegansBA_data[["modularity_celegansBA"]] <- modularity_ba_xarxa3
celegansBA_data[["authority_score_celegansBA"]] <- as_ba_xarxa3
celegansBA_data[["communities_celegansBA"]] <- comunitats_ba_xarxa3 
celegansBA_data[["SW_ratio_celegansBA"]] <- SM_Ratio_ba_xarxa3
celegansBA_data[["cliques_celegansBA"]] <- cliques_ba_xarxa3 
celegansBA_data[["cost_celegansBA"]] <- cost_ba_xarxa3
celegansBA_data[["nEdges_celegansBA"]] <- nEdges_ba_xarxa3 
celegansBA_data[["nNodes_celegansBA"]] <- nNodes_ba_xarxa3 
celegansBA_data[["rout_eff_celegansBA"]] <- rout_eff_ba_xarxa3 
celegansBA_data[["storage_cap_celegansBA"]] <- storage_cap_sw_xarxa2 


save(celegansBA_data, file="celegansBA_data.RData")

load("celegansBA_data.RData")




#----------------------------COMPARACIONS------------------------------------------------#


df_plot1 <- data.frame(deg=rep(0:max(deg),4),
                       freq=c(1-deg.distT,
                              c(1-deg.dist_xarxaT1,rep(1,max(deg)-max(deg_xarxa1))),
                              c(1-deg.dist_xarxaT2,rep(1,max(deg)-max(deg_xarxa2))),
                              c(1-deg.dist_xarxaT3[0:(max(deg)+1)])
                       ),
                       type=c(rep("c.elegans",max(deg)+1),
                              rep("Erd?s-Renyi",max(deg)+1),
                              rep("Small-World",max(deg)+1),
                              rep("Barabasi-Alberts",max(deg)+1)
                       )
)

library(ggplot2)
ggplot(df_plot1, aes(x=deg,y=freq,color=type)) + geom_point()

df_plot2 <- data.frame(deg=rep(0:max(deg),4),
                       freq=c(deg.distF,
                              c(deg.dist_xarxaF1,rep(0,max(deg)-max(deg_xarxa1))),
                              c(deg.dist_xarxaF2,rep(0,max(deg)-max(deg_xarxa2))),
                              c(deg.dist_xarxaF3[0:(max(deg)+1)])
                       ),
                       type=c(rep("c.elegans",max(deg)+1),
                              rep("Erd?s-Renyi",max(deg)+1),
                              rep("Small-World",max(deg)+1),
                              rep("Barabasi-Alberts",max(deg)+1)
                       )
)

library(ggplot2)
ggplot(df_plot2, aes(x=deg,y=freq,color=type)) + 
  geom_point() + 
  scale_y_log10() + 
  scale_x_log10(limits=c(5,100)) +
  geom_smooth(aes(fill=type)) +
  theme_classic()
