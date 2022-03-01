library("igraph")
library("plyr")







#----------------------------FLY EM -------------------------------------------------------------------------



archiuFly <- read.csv("v1.2_exported-traced-adjacencies-v1.2.csv",sep = ",",header = F,skip = 1)

#archiuFly[,1] <- as.numeric(as.factor(archiuFly[,1]))
#archiuFly[,2] <- as.numeric(as.factor(archiuFly[,2]))
#archiuFly[,3] <- as.numeric(as.factor(archiuFly[,3]))
#gFly <- graph_from_edgelist(as.matrix(archiuFly[,c(1,2)]), directed = F)
# 
#gFlyS <- simplify( gFly , remove.multiple = T, remove.loops = F, edge.attr.comb=c(weight="sum", type="ignore") )
# 
#save(gFlyS, file="gFlyS.RData")

load("gFlyS.RData")

#----#

degFly <- degree(gFlyS, mode="all")

deg.degFlyF <- degree_distribution(gFlyS, cumulative=F, mode="all")

deg.degFlyT <- degree_distribution(gFlyS, cumulative=T, mode="all")

tirad_census_gFlyS = triad_census(gFlyS)


#Hubs
hs_gFlyS <- hub_score(gFlyS, weights=NA)$vector


#Authorities
as_gFlyS<- authority_score(gFlyS, weights=NA)$vector


#Small-Worldness ratio: 
SM_Ratio_gFlyS = transitivity(gFlyS)/mean_distance(gFlyS)

#modularity
library(leiden)
adj_mtx_gFlyS  <- igraph::as_adjacency_matrix(gFlyS)
wtc_gFlyS <- leiden(adj_mtx_gFlyS)
modularity_gFlyS <- modularity(gFlyS, wtc_gFlyS)

#comunitats
comunitats_gFlyS <- cluster_louvain(gFlyS)
communitySize_gFlyS <- comunitats_gFlyS

#cliques
cliques_gFlyS <- cliques(gFlyS)

FlyEM_data <- list()
FlyEM_data[["deg_distF_FlyEM"]] <- deg.degFlyF
FlyEM_data[["deg_distT_FlyEM"]] <- deg.degFlyT
FlyEM_data[["triad_census_FlyEM"]] <- tirad_census_gFlyS
FlyEM_data[["hub_score_FlyEM"]] <- hs_gFlyS
FlyEM_data[["modularity_FlyEM"]] <- modularity_gFlyS
FlyEM_data[["authority_score_FlyEM"]] <- as_gFlyS 
FlyEM_data[["communitySize_FlyEM"]] <- communitySize_gFlyS 
FlyEM_data[["SW_ratio_FlyEM"]] <- SM_Ratio_gFlyS 
FlyEM_data[["cliques_FlyEM"]] <- cliques_gFlyS 


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

#modularity
library(leiden)
adj_mtx_er_Fly  <- igraph::as_adjacency_matrix(er_Fly)
wtc_er_Fly <- leiden(adj_mtx_er_Fly)
modularity_er_Fly <- modularity(er_Fly, wtc_er_Fly)

#comunitats
comunitats_er_Fly <- cluster_louvain(er_Fly)
communitySize_er_Fly <- comunitats_er_Fly

#cliques
cliques_er_Fly <- cliques(er_Fly)

FlyEM_ER_data <- list()
FlyEM_ER_data[["deg_distF_FlyEM_ER"]] <- deg.dist_erFlyF
FlyEM_ER_data[["deg_distT_FlyEM_ER"]] <- deg.dist_erFlyT
FlyEM_ER_data[["triad_census_FlyEM_ER"]] <- tirad_census_er_Fly
FlyEM_ER_data[["hub_score_FlyEM_ER"]] <- hs_er_Fly
FlyEM_ER_data[["modularity_FlyEM_ER"]] <- modularity_er_Fly
FlyEM_ER_data[["authority_score_FlyEM_ER"]] <- as_er_Fly 
FlyEM_ER_data[["communitySize_FlyEM_ER"]] <- communitySize_er_Fly 
FlyEM_ER_data[["SW_ratio_FlyEM_ER"]] <- SM_Ratio_er_Fly
FlyEM_ER_data[["cliques_FlyEM_ER"]] <- cliques_er_Fly 


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

#modularity
library(leiden)
adj_mtx_sw_Fly  <- igraph::as_adjacency_matrix(sw_Fly)
wtc_sw_Fly <- leiden(adj_mtx_sw_Fly)
modularity_sw_Fly <- modularity(sw_Fly, wtc_sw_Fly)

#comunitats
comunitats_sw_Fly <- cluster_louvain(sw_Fly)
communitySize_sw_Fly <- comunitats_sw_Fly

#cliques
cliques_sw_Fly <- cliques(sw_Fly)

FlyEM_SW_data <- list()
FlyEM_SW_data[["deg_distF_FlyEM_SW"]] <- deg.dist_swFlyF
FlyEM_SW_data[["deg_distT_FlyEM_SW"]] <- deg.dist_swFlyT
FlyEM_SW_data[["triad_census_FlyEM_SW"]] <- tirad_census_sw_Fly
FlyEM_SW_data[["hub_score_FlyEM_SW"]] <- hs_sw_Fly
FlyEM_SW_data[["modularity_FlyEM_SW"]] <- modularity_sw_Fly
FlyEM_SW_data[["authority_score_FlyEM_SW"]] <- as_sw_Fly 
FlyEM_SW_data[["communitySize_FlyEM_SW"]] <- communitySize_sw_Fly 
FlyEM_SW_data[["SW_ratio_FlyEM_SW"]] <- SM_Ratio_sw_Fly 
FlyEM_SW_data[["cliques_FlyEM_SW"]] <- cliques_sw_Fly 


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

#modularity
library(leiden)
adj_mtx_ba_fly  <- igraph::as_adjacency_matrix(ba_fly)
wtc_ba_fly <- leiden(adj_mtx_ba_fly)
modularity_ba_fly <- modularity(ba_fly, wtc_ba_fly)

#comunitats
comunitats_ba_fly <- cluster_louvain(ba_fly)
communitySize_ba_fly <- comunitats_ba_fly

#cliques
cliques_ba_fly <- cliques(ba_fly)

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


save(FlyEM_BA_data, file="FlyEM_BA_data.RData")

load("FlyEM_BA_data.RData")
