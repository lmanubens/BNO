library("igraph")
library("plyr")
library(leiden)
library(gmp)
library(Matrix)

archiuRata <- read.csv("affinity_map_nonflip_thres30_20211103.csv",sep = ",")

archiuDades <- read.csv("mouse_mesoscale_connectivity_matrix.csv",sep = ",",header=FALSE)
archiuNoms <- read.csv("mouse_mesoscale_regionname.csv")

noms <- as.character(archiuNoms$region_name)

names(archiuDades) <- noms
rownames(archiuDades) <- noms


gRatoli <-graph_from_adjacency_matrix(
  as.matrix(archiuDades),
  mode = "directed",
  weighted = TRUE,
  diag = TRUE,
  add.colnames = NULL,
  add.rownames = NA
)
gsRatoli <- simplify( gRatoli , remove.multiple = T, remove.loops = F, edge.attr.comb=c(weight="sum", type="ignore") )

lRatoli <- layout_with_fr(gsRatoli)

plot(gsRatoli, rescale=F, layout=lRatoli*0.4)


degRatoli <- degree(gsRatoli, mode="all")

deg.distRatoliF <- degree_distribution(gsRatoli, cumulative=F, mode="all")

deg.distRatoliT <- degree_distribution(gsRatoli, cumulative=T, mode="all")



#Triad census
tirad_gsRatoli = triad_census(gsRatoli)

#Hubs
hs_gsRatoli <- hub_score(gsRatoli, weights=NA)$vector

#Authorities
as_gsRatoli <- authority_score(gsRatoli, weights=NA)$vector


#Small-Worldness ratio: 
SM_Ratio_gsRatoli = transitivity(gsRatoli)/mean_distance(gsRatoli)

#modularity y comunitats
adj_mtx_gsRatoli  <- igraph::as_adjacency_matrix(gsRatoli)
camunitats_gsRatoli <- leiden(adj_mtx_gsRatoli)
modularity_gsRatoli <- modularity(gsRatoli, camunitats_gsRatoli)


#cliques
cliques_gsRatoli <- cliques(gsRatoli)

#n nodes
nEdges_gsRatoli <- gsize(gsRatoli)
cost_gsRatoli <- sum(E(gsRatoli)$weight)
nNodes_gsRatoli <- gorder(gsRatoli)

#routing efficiency
am_gsRatoli  <- get.adjacency(gsRatoli)
distances_gsRatoli <- distances(gsRatoli , mode = "in", weights=1/E(gsRatoli)$weight)


rout_effs_gsRatoli <- as.numeric(as.bigz(1/distances_gsRatoli)/as.bigz(dim(am_gsRatoli)[[1]]*(dim(am_gsRatoli)[[1]]-1)))
rout_effs_matrix_gsRatoli <- as.matrix(rout_effs_er_FlyO)
rout_effs_matrix_gsRatoli <- rout_effs_matrix_gsRatoli[is.finite(rout_effs_matrix_gsRatoli)]
rout_eff_gsRatoli <- as.numeric(sum(rout_effs_matrix_gsRatoli))

#storage capacity
storage_caps_gsRatoli <- numeric()
for(j in 1:dim(am_gsRatoli)[[1]]){
  storage_caps_gsRatoli <- append(storage_caps_gsRatoli, as.integer(2*log2(chooseZ(as.bigz(sum(am_gsRatoli[,j]) + nnzero(am_gsRatoli[,j]) -1), as.bigz(sum(am_gsRatoli[,j]))))))
}
storage_cap_gsRatoli <- mean(storage_caps_gsRatoli)#,na.rm=T)

ratoli_data <- list()
ratoli_data[["deg_distF_ratoli"]] <- deg.distRatoliF
ratoli_data[["deg_distT_ratoli"]] <- deg.distRatoliT
ratoli_data[["triad_census_ratoli"]] <- tirad_gsRatoli
ratoli_data[["hub_score_ratoli"]] <- hs_gsRatoli
ratoli_data[["modularity_ratoli"]] <- modularity_gsRatoli
ratoli_data[["communities_ratoli"]] <- camunitats_gsRatoli
ratoli_data[["authority_score_ratoli"]] <- as_gsRatoli 
ratoli_data[["SW_ratio_ratoli"]] <- SM_Ratio_gsRatoli
ratoli_data[["cliques_ratoli"]] <- cliques_gsRatoli 
ratoli_data[["cost_ratoli"]] <- cost_gsRatoli
ratoli_data[["nEdges_ratoli"]] <- nEdges_gsRatoli
ratoli_data[["nNodes_ratoli"]] <- nNodes_gsRatoli 
ratoli_data[["rout_eff_ratoli"]] <- rout_eff_gsRatoli 
ratoli_data[["storage_cap_ratoli"]] <- storage_cap_gsRatoli

save(ratoli_data, file="ratoli_data.RData")

load("ratoli_data.RData")





