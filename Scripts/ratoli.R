
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

#plot(gsRatoli, rescale=F, layout=lRatoli*0.4)

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

#modularity
library(leiden)
adj_mtx_gsRatoli  <- igraph::as_adjacency_matrix(gsRatoli)
wtc_gsRatoli <- leiden(adj_mtx_gsRatoli)
modularity_gsRatoli <- modularity(gsRatoli, wtc_gsRatoli)

#comunitats
comunitats_ba_flyO <- cluster_louvain(ba_flyO)
communitySize_ba_flyO <- comunitats_ba_flyO

ratoli_data <- list()
ratoli_data[["deg_distF_ratoli"]] <- deg.distRatoliF
ratoli_data[["deg_distT_ratoli"]] <- deg.distRatoliT
ratoli_data[["triad_census_ratoli"]] <- tirad_gsRatoli
ratoli_data[["hub_score_ratoli"]] <- hs_gsRatoli
ratoli_data[["modularity_ratoli"]] <- modularity_gsRatoli
ratoli_data[["authority_score_ratoli"]] <- as_gsRatoli 
ratoli_data[["SW_ratio_ratoli"]] <- SM_Ratio_gsRatoli 

save(ratoli_data, file="ratoli_data.RData")

load("ratoli_data.RData")
