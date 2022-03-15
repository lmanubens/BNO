library("igraph")
library("plyr")
library(leiden)
library(gmp)
library(Matrix)







#----------------------------FLY EM -------------------------------------------------------------------------



archiuFly <- read.csv("v1.2_exported-traced-adjacencies-v1.2.csv",sep = ",",header = F,skip = 1)

#
archiuFly[,1] <- as.numeric(as.factor(archiuFly[,1]))
archiuFly[,2] <- as.numeric(as.factor(archiuFly[,2]))
archiuFly[,3] <- as.numeric(as.factor(archiuFly[,3]))
gFly <- graph_from_edgelist(as.matrix(archiuFly[,c(1,2)]), directed = T)

gFlyS <- simplify( gFly , remove.multiple = T, remove.loops = F, edge.attr.comb=c(weight="sum", type="ignore") )
 
save(gFlyS, file="gFlyS.RData")
#

load("gFlyS.RData")

#----#

degFly <- degree(gFlyS, mode="all")

deg.degFlyF <- degree_distribution(gFlyS, cumulative=F, mode="all")

deg.degFlyT <- degree_distribution(gFlyS, cumulative=T, mode="all")

#tirad_census
tirad_census_gFlyS = triad_census(gFlyS)


#Hubs
hs_gFlyS <- hub_score(gFlyS, weights=NA)$vector


#Authorities
as_gFlyS<- authority_score(gFlyS, weights=NA)$vector


#Small-Worldness ratio: 
SM_Ratio_gFlyS = transitivity(gFlyS)/mean_distance(gFlyS)


#modularity y comunitats
adj_mtx_gFlyS  <- igraph::as_adjacency_matrix(gFlyS)
camunitats_gFlyS <- leiden(adj_mtx_gFlyS)
modularity_gFlyS <- modularity(gFlyS, camunitats_gFlyS)


#cliques
cliques_gFlyS <- cliques(gFlyS)

#n nodes
nEdges_gFlyS <- gsize(gFlyS)
cost_gFlyS <- sum(E(gFlyS)$weight)
nNodes_gFlyS <- gorder(gFlyS)

#routing efficiency
am_gFlyS  <- get.adjacency(gFlyS)
distances_gFlyS <- distances(gFlyS , mode = "in", weights=1/E(gFlyS)$weight)


rout_effs_gFlyS <- as.numeric(as.bigz(1/distances_gFlyS)/as.bigz(dim(am_gFlyS)[[1]]*(dim(am_gFlyS)[[1]]-1)))
rout_effs_matrix_gFlyS <- as.matrix(rout_effs_gFlyS)
rout_effs_matrix_gFlyS <- rout_effs_matrix_gFlyS[is.finite(rout_effs_matrix_gFlyS)]
rout_eff_gFlyS <- as.numeric(sum(rout_effs_matrix_gFlyS))

#storage capacity
storage_caps_gFlyS <- numeric()
for(j in 1:dim(am_gFlyS)[[1]]){
  storage_caps_gFlyS <- append(storage_caps_gFlyS, as.integer(2*log2(chooseZ(as.bigz(sum(am_gFlyS[,j]) + nnzero(am_gFlyS[,j]) -1), as.bigz(sum(am_gFlyS[,j]))))))
}
storage_cap_gFlyS <- mean(storage_caps_gFlyS)#,na.rm=T)

FlyEM_data <- list()
FlyEM_data[["deg_distF_FlyEM"]] <- deg.degFlyF
FlyEM_data[["deg_distT_FlyEM"]] <- deg.degFlyT
FlyEM_data[["triad_census_FlyEM"]] <- tirad_census_gFlyS
FlyEM_data[["hub_score_FlyEM"]] <- hs_gFlyS
FlyEM_data[["modularity_FlyEM"]] <- modularity_gFlyS
FlyEM_data[["authority_score_FlyEM"]] <- as_gFlyS 
FlyEM_data[["communities_FlyEM"]] <- camunitats_gFlyS 
FlyEM_data[["SW_ratio_FlyEM"]] <- SM_Ratio_gFlyS 
FlyEM_data[["cliques_FlyEM"]] <- cliques_gFlyS
FlyEM_data[["cost_FlyEM"]] <- cost_gFlyS
FlyEM_data[["nEdges_FlyEM"]] <- nEdges_gFlyS
FlyEM_data[["nNodes_FlyEM"]] <- nNodes_gFlyS
FlyEM_data[["rout_eff_FlyEM"]] <- rout_eff_gFlyS
FlyEM_data[["storage_cap_FlyEM"]] <- storage_cap_gFlyS


save(FlyEM_data, file="FlyEM_data.RData")

load("FlyEM_data.RData")










#-------------------------------------------------------------------ER_FLY---------------------------------------------------#
set.seed(2706)

er_Fly <- sample_gnm(n=length(unique(archiuFly[,1])), m=length(E(gFlyS))) #n= NUmero de nodes (unique(archiuFly[,1])) m=numero de links (length(E(gFlyS)))
#plot(er_Fly, vertex.size=6, vertex.label=NA, layout=layout_in_circle)  


deg_erFly <- degree(er_Fly, mode="all")

deg.dist_erFlyF <- degree_distribution(er_Fly, cumulative=F, mode="all")

deg.dist_erFlyT <- degree_distribution(er_Fly, cumulative=T, mode="all")

tirad_census_er_Fly = triad_census(er_Fly)


#Hubs
hs_er_Fly <- hub_score(er_Fly, weights=NA)$vector


#Authorities
as_er_Fly <- authority_score(er_Fly, weights=NA)$vector


#Small-Worldness ratio: 
SM_Ratio_er_Fly = transitivity(er_Fly)/mean_distance(er_Fly)

#modularity y comunitats
adj_mtx_er_Fly  <- igraph::as_adjacency_matrix(er_Fly)
camunitats_er_Fly <- leiden(adj_mtx_er_Fly)
modularity_er_Fly <- modularity(er_Fly, camunitats_er_Fly)


#cliques
cliques_er_Fly <- cliques(er_Fly)

#n nodes
nEdges_er_Fly <- gsize(er_Fly)
cost_er_Fly <- sum(E(er_Fly)$weight)
nNodes_er_Fly <- gorder(er_Fly)

#routing efficiency
am_er_Fly  <- get.adjacency(er_Fly)
distances_er_Fly <- distances(er_Fly , mode = "in", weights=1/E(er_Fly)$weight)


rout_effs_er_Fly <- as.numeric(as.bigz(1/distances_er_Fly)/as.bigz(dim(am_er_Fly)[[1]]*(dim(am_er_Fly)[[1]]-1)))
rout_effs_matrix_er_Fly <- as.matrix(rout_effs_er_Fly)
rout_effs_matrix_er_Fly <- rout_effs_matrix_er_Fly[is.finite(rout_effs_matrix_er_Fly)]
rout_eff_er_Fly <- as.numeric(sum(rout_effs_matrix_er_Fly))

#storage capacity
storage_caps_er_Fly <- numeric()
for(j in 1:dim(am_er_Fly)[[1]]){
  storage_caps_er_Fly <- append(storage_caps_er_Fly, as.integer(2*log2(chooseZ(as.bigz(sum(am_er_Fly[,j]) + nnzero(am_er_Fly[,j]) -1), as.bigz(sum(am_er_Fly[,j]))))))
}
storage_cap_er_Fly <- mean(storage_caps_er_Fly)#,na.rm=T)


FlyEM_ER_data <- list()
FlyEM_ER_data[["deg_distF_FlyEM_ER"]] <- deg.dist_erFlyF
FlyEM_ER_data[["deg_distT_FlyEM_ER"]] <- deg.dist_erFlyT
FlyEM_ER_data[["triad_census_FlyEM_ER"]] <- tirad_census_er_Fly
FlyEM_ER_data[["hub_score_FlyEM_ER"]] <- hs_er_Fly
FlyEM_ER_data[["modularity_FlyEM_ER"]] <- modularity_er_Fly
FlyEM_ER_data[["authority_score_FlyEM_ER"]] <- as_er_Fly 
FlyEM_ER_data[["communities_FlyEM_ER"]] <- camunitats_er_Fly 
FlyEM_ER_data[["SW_ratio_FlyEM_ER"]] <- SM_Ratio_er_Fly
FlyEM_ER_data[["cliques_FlyEM_ER"]] <- cliques_er_Fly 
FlyEM_ER_data[["cost_FlyEM_ER"]] <- cost_er_Fly
FlyEM_ER_data[["nEdges_FlyEM_ER"]] <- nEdges_er_Fly
FlyEM_ER_data[["nNodes_FlyEM_ER"]] <- nNodes_er_Fly
FlyEM_ER_data[["rout_eff_FlyEM_ER"]] <- rout_eff_er_Fly 
FlyEM_ER_data[["storage_cap_FlyEM_ER"]] <- storage_cap_er_Fly


save(FlyEM_ER_data, file="FlyEM_ER_data.RData")

load("FlyEM_ER_data.RData")


#-------------------------------------------------------------------SW_FLY---------------------------------------------------#

sw_Fly <- sample_smallworld(dim=1, size=length(unique(archiuFly[,1])), nei=11, p=0.1)
#plot(sw_Fly, vertex.size=6, vertex.label=NA, layout=layout_in_circle)


deg_swFly <- degree(sw_Fly, mode="all")

deg.dist_swFlyF <- degree_distribution(sw_Fly, cumulative=F, mode="all")

deg.dist_swFlyT <- degree_distribution(sw_Fly, cumulative=T, mode="all")

tirad_census_sw_Fly = triad_census(sw_Fly)

#Hubs
hs_sw_Fly <- hub_score(sw_Fly, weights=NA)$vector

#Authorities
as_sw_Fly <- authority_score(sw_Fly, weights=NA)$vector

#Small-Worldness ratio: 
SM_Ratio_sw_Fly = transitivity(sw_Fly)/mean_distance(sw_Fly)

#modularity y comunitats
adj_mtx_sw_Fly <- igraph::as_adjacency_matrix(sw_Fly)
camunitats_sw_Fly <- leiden(adj_mtx_sw_Fly)
modularity_sw_Fly <- modularity(sw_Fly, camunitats_sw_Fly)


#cliques
cliques_sw_Fly <- cliques(sw_Fly)

#n nodes
nEdges_sw_Fly <- gsize(sw_Fly)
cost_sw_Fly <- sum(E(sw_Fly)$weight)
nNodes_sw_Fly <- gorder(sw_Fly)

#routing efficiency
am_sw_Fly  <- get.adjacency(sw_Fly)
distances_sw_Fly <- distances(sw_Fly , mode = "in", weights=1/E(sw_Fly)$weight)


rout_effs_sw_Fly <- as.numeric(as.bigz(1/distances_sw_Fly)/as.bigz(dim(am_sw_Fly)[[1]]*(dim(am_sw_Fly)[[1]]-1)))
rout_effs_matrix_sw_Fly <- as.matrix(rout_effs_sw_Fly)
rout_effs_matrix_sw_Fly <- rout_effs_matrix_sw_Fly[is.finite(rout_effs_matrix_sw_Fly)]
rout_eff_sw_Fly <- as.numeric(sum(rout_effs_matrix_sw_Fly))

#storage capacity
storage_caps_sw_Fly <- numeric()
for(j in 1:dim(am_sw_Fly)[[1]]){
  storage_caps_sw_Fly <- append(storage_caps_sw_Fly, as.integer(2*log2(chooseZ(as.bigz(sum(am_sw_Fly[,j]) + nnzero(am_sw_Fly[,j]) -1), as.bigz(sum(am_sw_Fly[,j]))))))
}
storage_cap_sw_Fly <- mean(storage_caps_sw_Fly)#,na.rm=T)

FlyEM_SW_data <- list()
FlyEM_SW_data[["deg_distF_FlyEM_SW"]] <- deg.dist_swFlyF
FlyEM_SW_data[["deg_distT_FlyEM_SW"]] <- deg.dist_swFlyT
FlyEM_SW_data[["triad_census_FlyEM_SW"]] <- tirad_census_sw_Fly
FlyEM_SW_data[["hub_score_FlyEM_SW"]] <- hs_sw_Fly
FlyEM_SW_data[["modularity_FlyEM_SW"]] <- modularity_sw_Fly
FlyEM_SW_data[["authority_score_FlyEM_SW"]] <- as_sw_Fly 
FlyEM_SW_data[["communities_FlyEM_SW"]] <- camunitats_sw_Fly 
FlyEM_SW_data[["SW_ratio_FlyEM_SW"]] <- SM_Ratio_sw_Fly 
FlyEM_SW_data[["cliques_FlyEM_SW"]] <- cliques_sw_Fly  
FlyEM_SW_data[["cost_FlyEM_SW"]] <- cost_sw_Fly 
FlyEM_SW_data[["nEdges_FlyEM_SW"]] <- nEdges_sw_Fly 
FlyEM_SW_data[["nNodes_FlyEM_SW"]] <- nNodes_sw_Fly 
FlyEM_SW_data[["rout_eff_FlyEM_SW"]] <- rout_eff_sw_Fly 
FlyEM_SW_data[["storage_cap_FlyEM_SW"]] <- storage_cap_sw_Fly 


save(FlyEM_SW_data, file="FlyEM_SW_data.RData")

load("FlyEM_SW_data.RData")



#-------------------------------------------------------------------BA_FLY---------------------------------------------------#

ba_fly <-  sample_pa(n=length(unique(archiuFly[,1])), power=1, m=11,  directed=F)
#plot(ba_fly, vertex.size=6, vertex.label=NA, layout=layout_in_circle)


deg_baFly <- degree(ba_fly, mode="all")

deg.dist_baFlyF <- degree_distribution(ba_fly, cumulative=F, mode="all")

deg.dist_baFyT <- degree_distribution(ba_fly, cumulative=T, mode="all")

tirad_census_ba_fly = triad_census(ba_fly)

#Hubs
hs_ba_fly <- hub_score(ba_fly, weights=NA)$vector


#Authorities
as_baFly <- authority_score(ba_fly, weights=NA)$vector


#Small-Worldness ratio: 
SM_Ratio_ba_fly= transitivity(ba_fly)/mean_distance(ba_fly)

#modularity y comunitats
adj_mtx_ba_fly <- igraph::as_adjacency_matrix(ba_fly)
camunitats_ba_fly <- leiden(adj_mtx_ba_fly)
modularity_ba_fly <- modularity(ba_fly, camunitats_ba_fly)


#cliques
cliques_ba_fly <- cliques(ba_fly)

#n nodes
nEdges_ba_fly <- gsize(ba_fly)
cost_ba_fly <- sum(E(ba_fly)$weight)
nNodes_ba_fly <- gorder(ba_fly)

#routing efficiency
am_ba_fly  <- get.adjacency(ba_fly)
distances_ba_fly <- distances(ba_fly , mode = "in", weights=1/E(ba_fly)$weight)


rout_effs_ba_fly <- as.numeric(as.bigz(1/distances_ba_fly)/as.bigz(dim(am_ba_fly)[[1]]*(dim(am_ba_fly)[[1]]-1)))
rout_effs_matrix_ba_fly <- as.matrix(rout_effs_ba_fly)
rout_effs_matrix_ba_fly <- rout_effs_matrix_ba_fly[is.finite(rout_effs_matrix_ba_fly)]
rout_eff_ba_fly <- as.numeric(sum(rout_effs_matrix_ba_fly))

#storage capacity
storage_caps_ba_fly <- numeric()
for(j in 1:dim(am_ba_fly)[[1]]){
  storage_caps_ba_fly <- append(storage_caps_ba_fly, as.integer(2*log2(chooseZ(as.bigz(sum(am_ba_fly[,j]) + nnzero(am_ba_fly[,j]) -1), as.bigz(sum(am_ba_fly[,j]))))))
}
storage_cap_ba_fly <- mean(storage_caps_ba_fly)#,na.rm=T)

FlyEM_BA_data <- list()
FlyEM_BA_data[["deg_distF_FlyEM_BA"]] <- deg.dist_baFlyF
FlyEM_BA_data[["deg_distT_FlyEM_BA"]] <- deg.dist_baFyT
FlyEM_BA_data[["triad_census_FlyEM_BA"]] <- tirad_census_ba_fly
FlyEM_BA_data[["hub_score_FlyEM_BA"]] <- hs_ba_fly
FlyEM_BA_data[["modularity_FlyEM_BA"]] <- modularity_ba_fly
FlyEM_BA_data[["authority_score_FlyEM_BA"]] <- as_baFly 
FlyEM_BA_data[["communitySize_FlyEM_BA"]] <- communitySize_ba_fly 
FlyEM_BA_data[["SW_ratio_FlyEM_BA"]] <- SM_Ratio_ba_fly 
FlyEM_BA_data[["cliques_FlyEM_BA"]] <- communitySize_ba_fly
FlyEM_BA_data[["cost_FlyEM_BA"]] <- cost_ba_fly 
FlyEM_BA_data[["nEdges_FlyEM_BA"]] <- nEdges_ba_fly 
FlyEM_BA_data[["nNodes_FlyEM_BA"]] <- nNodes_ba_fly 
FlyEM_BA_data[["rout_eff_FlyEM_BA"]] <- rout_eff_ba_fly 
FlyEM_BA_data[["storage_cap_FlyEM_BA"]] <- storage_cap_ba_fly


save(FlyEM_BA_data, file="FlyEM_BA_data.RData")

load("FlyEM_BA_data.RData")
