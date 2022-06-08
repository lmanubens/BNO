library("igraph")
library("plyr")
library(leiden)
library(gmp)
library(Matrix)

#archiuRata <- read.csv("affinity_map_nonflip_thres30_20211103.csv",sep = ",")

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
  add.rownames = NULL
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
adj_mtx_gsRatoli  <- igraph::as_adjacency_matrix(gsRatoli, attr="weight")
camunitats_gsRatoli <- leiden(adj_mtx_gsRatoli)
modularity_gsRatoli <- modularity(gsRatoli, camunitats_gsRatoli)


#cliques

#n nodes
nEdges_gsRatoli <- gsize(gsRatoli)
#cost_gsRatoli <- sum(E(gsRatoli)$weight)
nNodes_gsRatoli <- gorder(gsRatoli)

#routing efficiency
#am_gsRatoli  <- get.adjacency(gsRatoli, attr="weight")
#distances_gsRatoli <- distances(gsRatoli , mode = "in", weights=1/E(gsRatoli)$weight)


#rout_effs_gsRatoli <- as.numeric(as.bigz(1/distances_gsRatoli)/as.bigz(dim(am_gsRatoli)[[1]]*(dim(am_gsRatoli)[[1]]-1)))
#rout_effs_matrix_gsRatoli <- as.matrix(rout_effs_gsRatoli)
#rout_effs_matrix_gsRatoli <- rout_effs_matrix_gsRatoli[is.finite(rout_effs_matrix_gsRatoli)]
#rout_eff_gsRatoli <- as.numeric(sum(rout_effs_matrix_gsRatoli))

#storage capacity
#storage_caps_gsRatoli <- numeric()
#for(j in 1:dim(am_gsRatoli)[[1]]){
#  storage_caps_gsRatoli <- append(storage_caps_gsRatoli, as.integer(2*log2(chooseZ(as.bigz(sum(am_gsRatoli[,j]) + nnzero(am_gsRatoli[,j]) -1), as.bigz(sum(am_gsRatoli[,j]))))))
#}
#storage_cap_gsRatoli <- mean(storage_caps_gsRatoli)#,na.rm=T)

ratoli_data <- list()
ratoli_data[["deg_distF_ratoli"]] <- deg.distRatoliF
ratoli_data[["deg_distT_ratoli"]] <- deg.distRatoliT
ratoli_data[["triad_census_ratoli"]] <- tirad_gsRatoli
ratoli_data[["hub_score_ratoli"]] <- hs_gsRatoli
ratoli_data[["modularity_ratoli"]] <- modularity_gsRatoli
ratoli_data[["communities_ratoli"]] <- camunitats_gsRatoli
ratoli_data[["authority_score_ratoli"]] <- as_gsRatoli 
ratoli_data[["SW_ratio_ratoli"]] <- SM_Ratio_gsRatoli
#ratoli_data[["cost_ratoli"]] <- cost_gsRatoli
ratoli_data[["nEdges_ratoli"]] <- nEdges_gsRatoli
ratoli_data[["nNodes_ratoli"]] <- nNodes_gsRatoli 
#ratoli_data[["rout_eff_ratoli"]] <- rout_eff_gsRatoli 
#ratoli_data[["storage_cap_ratoli"]] <- storage_cap_gsRatoli

save(ratoli_data, file="ratoli_data.RData")

load("ratoli_data.RData")


#------------------------------------------------------ Erdos-Renyi random graph model -------------------------------------#


set.seed(2706)

er_Ratoli <- sample_gnm(n=length(unique(row.names(archiuDades))), m=length(E(gsRatoli)),directed=T)
V(er_Ratoli)$name <- V(gsRatoli)$name
E(er_Ratoli)$weight <- E(gsRatoli)$weight

deg.dist_er_RatoliF <- degree_distribution(er_Ratoli, cumulative=F, mode="all")

deg.dist_er_RatoliT <- degree_distribution(er_Ratoli, cumulative=T, mode="all")



#Triad census
tirad_er_Ratoli = triad_census(er_Ratoli)

#Hubs
hs_er_Ratoli <- hub_score(er_Ratoli, weights=NA)$vector

#Authorities
as_er_Ratoli <- authority_score(er_Ratoli, weights=NA)$vector


#Small-Worldness ratio: 
SM_Ratio_er_Ratoli = transitivity(er_Ratoli)/mean_distance(er_Ratoli)

#modularity y comunitats
adj_mtx_er_Ratoli  <- igraph::as_adjacency_matrix(er_Ratoli, attr="weight")
camunitats_er_Ratoli <- leiden(adj_mtx_er_Ratoli)
modularity_er_Ratoli <- modularity(er_Ratoli, camunitats_er_Ratoli)


#cliques

#n nodes
nEdges_er_Ratoli <- gsize(er_Ratoli)
#cost_er_Ratoli <- sum(E(er_Ratoli)$weight)
nNodes_er_Ratoli <- gorder(er_Ratoli)

#routing efficiency
#am_er_Ratoli  <- get.adjacency(er_Ratoli, attr="weight")
#distances_er_Ratoli <- distances(er_Ratoli , mode = "in", weights=1/E(er_Ratoli)$weight)


#rout_effs_er_Ratoli <- as.numeric(as.bigz(1/distances_er_Ratoli)/as.bigz(dim(am_er_Ratoli)[[1]]*(dim(am_er_Ratoli)[[1]]-1)))
#rout_effs_matrix_er_Ratoli <- as.matrix(rout_effs_er_Ratoli)
#rout_effs_matrix_er_Ratoli <- rout_effs_matrix_er_Ratoli[is.finite(rout_effs_matrix_er_Ratoli)]
#rout_eff_er_Ratoli <- as.numeric(sum(rout_effs_matrix_er_Ratoli))

#storage capacity
#storage_caps_er_Ratoli <- numeric()
#for(j in 1:dim(am_er_Ratoli)[[1]]){
#  storage_caps_er_Ratoli <- append(storage_caps_er_Ratoli, as.integer(2*log2(chooseZ(as.bigz(sum(am_er_Ratoli[,j]) + nnzero(am_er_Ratoli[,j]) -1), as.bigz(sum(am_er_Ratoli[,j]))))))
#}
#storage_cap_er_Ratoli <- mean(storage_caps_er_Ratoli)#,na.rm=T)

ratoli_ER_data <- list()
ratoli_ER_data[["deg_distF_ratoli_ER"]] <- deg.dist_er_RatoliF
ratoli_ER_data[["deg_distT_ratoli_ER"]] <- deg.dist_er_RatoliT
ratoli_ER_data[["triad_census_ratoli_ER"]] <- tirad_er_gsRatoli
ratoli_ER_data[["hub_score_ratoli_ER"]] <- hs_er_gsRatoli
ratoli_ER_data[["modularity_ratoli_ER"]] <- modularity_er_gsRatoli
ratoli_ER_data[["communities_ratoli_ER"]] <- camunitats_er_gsRatoli
ratoli_ER_data[["authority_score_ratoli_ER"]] <- as_er_Ratoli 
ratoli_ER_data[["SW_ratio_ratoli_ER"]] <- SM_Ratio_er_Ratoli
#ratoli_ER_data[["cost_ratoli_ER"]] <- cost_er_Ratoli
ratoli_ER_data[["nEdges_ratoli_ER"]] <- nEdges_er_Ratoli
ratoli_ER_data[["nNodes_ratoli_ER"]] <- nNodes_er_Ratoli 
#ratoli_ER_data[["rout_eff_ratoli_ER"]] <- rout_eff_er_Ratoli 
#ratoli_ER_data[["storage_cap_ratoli_ER"]] <- storage_cap_er_Ratoli

save(ratoli_ER_data, file="ratoli_ER_data.RData")

load("ratoli_ER_data.RData")


#----------------------------------------------------SW-----------------------------------------------------
sw_ratoli <- rewire(graph.lattice(length(unique(row.names(archiuDades))), nei=200, dir=TRUE, mutual=TRUE), 
                    with = each_edge(p = .02, loops = FALSE))
V(sw_ratoli)$name <- V(gsRatoli)$name

pesos_sw_ratoli <- sample(E(gsRatoli)$weight,length(E(sw_ratoli)))
E(sw_ratoli)$weight <- pesos_sw_ratoli


deg_sw_ratoli <- degree(sw_ratoli, mode="all")


deg.dist_sw_ratoliF <- degree_distribution(sw_ratoli, cumulative=F, mode="all")

deg.dist_sw_ratoliT <- degree_distribution(sw_ratoli, cumulative=T, mode="all")


tirad_sw_ratoli = triad_census(sw_ratoli)

#Hubs
hs_sw_ratoli <- hub_score(sw_ratoli, weights=NA)$vector

#Authorities
as_sw_ratoli <- authority_score(sw_ratoli, weights=NA)$vector

#Small-Worldness ratio: 
SM_Ratio_sw_ratoli = transitivity(sw_ratoli)/mean_distance(sw_ratoli)

#modularity y comunitats
adj_mtx_sw_ratoli  <- igraph::as_adjacency_matrix(sw_ratoli, attr="weight")
comunitats_sw_ratoli <- leiden(adj_mtx_sw_ratoli)
modularity_sw_ratoli <- modularity(sw_ratoli, comunitats_sw_ratoli)

#cliques


#n nodes
nEdges_sw_ratoli <- gsize(sw_ratoli)
#cost_sw_ratoli <- sum(E(sw_ratoli)$weight)
nNodes_sw_ratoli <- gorder(sw_ratoli)

#routing efficiency
#am_sw_ratoli <- get.adjacency(sw_ratoli, attr="weight")
#distances_sw_ratoli <- distances(sw_ratoli , mode = "in", weights=1/E(sw_ratoli)$weight)

#rout_effs_sw_ratoli <- as.numeric(as.bigz(1/distances_sw_ratoli)/as.bigz(dim(am_sw_ratoli)[[1]]*(dim(am_sw_ratoli)[[1]]-1)))
#rout_effs_matrix_sw_ratoli <- as.matrix(rout_effs_sw_ratoli)
#rout_effs_matrix_sw_ratoli <- rout_effs_matrix_sw_ratoli[is.finite(rout_effs_matrix_sw_ratoli)]
#rout_eff_sw_ratoli <- as.numeric(sum(rout_effs_matrix_sw_ratoli))

#storage capacity
#storage_caps_sw_ratoli <- numeric()
#for(j in 1:dim(am_sw_ratoli)[[1]]){
#  storage_caps_sw_ratoli <- append(storage_caps_sw_ratoli, as.integer(2*log2(chooseZ(as.bigz(sum(am_sw_ratoli[,j]) + nnzero(am_sw_ratoli[,j]) -1), as.bigz(sum(am_sw_ratoli[,j]))))))
#}
#storage_cap_sw_ratoli <- mean(storage_caps_sw_ratoli)#,na.rm=T)

ratoli_SW_data <- list()
ratoli_SW_data[["deg_distF_ratoli_SW"]] <- deg.dist_sw_ratoliF
ratoli_SW_data[["deg_distT_ratoli_SW"]] <- deg.dist_sw_ratoliT
ratoli_SW_data[["triad_census_ratoli_SW"]] <- tirad_sw_ratoli
ratoli_SW_data[["hub_score_ratoli_SW"]] <- hs_sw_ratoli
ratoli_SW_data[["modularity_ratoli_SW"]] <- modularity_sw_ratoli
ratoli_SW_data[["authority_score_ratoli_SW"]] <- as_sw_ratoli 
ratoli_SW_data[["communities_ratoli_SW"]] <- comunitats_sw_ratoli 
ratoli_SW_data[["SW_ratio_ratoli_SW"]] <- SM_Ratio_sw_ratoli 
ratoli_SW_data[["cliques_ratoli_SW"]] <- cliques_sw_ratoli 
#ratoli_SW_data[["cost_ratoli_SW"]] <- cost_sw_ratoli
ratoli_SW_data[["nEdges_ratoli_SW"]] <- nEdges_sw_ratoli
ratoli_SW_data[["nNodes_ratoli_SW"]] <- nNodes_sw_ratoli 
#ratoli_SW_data[["rout_eff_ratoli_SW"]] <- rout_eff_sw_ratoli 
#ratoli_SW_data[["storage_cap_ratoli_SW"]] <- storage_cap_sw_ratoli


save(ratoli_SW_data, file="ratoli_SW_data.RData")

load("ratoli_SW_data.RData")
#----------------------------- Barabasi-Albert preferential attachment model for scale-free graphs -------------------------#


ba_ratoli <-  sample_pa(n=length(unique(row.names(archiuDades))), power=1,  m=length(E(gsRatoli)),  directed=T)
pesos_ba_ratoli <- sample(E(gsRatoli)$weight,length(E(ba_ratoli)))
E(ba_ratoli)$weight <- pesos_ba_ratoli



deg_ba_ratoli <- degree(ba_ratoli, mode="all")

deg.dist_ba_ratoliF <- degree_distribution(ba_ratoli, cumulative=F, mode="all")

deg.dist_ba_ratoliT <- degree_distribution(ba_ratoli, cumulative=T, mode="all")

tirad_census_ba_ratoli = triad_census(ba_ratoli)

#Hubs
hs_ba_ratoli <- hub_score(ba_ratoli, weights=NA)$vector

#Authorities
as_ba_ratoli <- authority_score(ba_ratoli, weights=NA)$vector

#Small-Worldness ratio: 
SM_Ratio_ba_ratoli = transitivity(ba_ratoli)/mean_distance(ba_ratoli)

#modularity y comunitats
adj_mtx_ba_ratoli  <- igraph::as_adjacency_matrix(ba_ratoli, attr="weight")
comunitats_ba_ratoli <- leiden(adj_mtx_ba_ratoli)
modularity_ba_ratoli <- modularity(ba_ratoli, comunitats_ba_ratoli)


#cliques


#n nodes
nEdges_ba_ratoli <- gsize(ba_ratoli)
#cost_ba_ratoli <- sum(E(ba_ratoli)$weight)
nNodes_ba_ratoli <- gorder(ba_ratoli)

#routing efficiency
#am_ba_ratoli <- get.adjacency(ba_ratoli, attr="weight")
#distances_ba_ratoli <- distances(ba_ratoli , mode = "in", weights=1/E(ba_ratoli)$weight)


#rout_effs_ba_ratoli <- as.numeric(as.bigz(1/distances_ba_ratoli)/as.bigz(dim(am_ba_ratoli)[[1]]*(dim(am_ba_ratoli)[[1]]-1)))
#rout_effs_matrix_ba_ratoli <- as.matrix(rout_effs_ba_ratoli)
#rout_effs_matrix_ba_ratoli <- rout_effs_matrix_ba_ratoli[is.finite(rout_effs_matrix_ba_ratoli)]
#rout_eff_ba_ratoli <- as.numeric(sum(rout_effs_matrix_ba_ratoli))

#storage capacity
#storage_caps_ba_ratoli <- numeric()
#for(j in 1:dim(am_ba_ratoli)[[1]]){
#  storage_caps_ba_ratoli <- append(storage_caps_ba_ratoli, as.integer(2*log2(chooseZ(as.bigz(sum(am_ba_ratoli[,j]) + nnzero(am_ba_ratoli[,j]) -1), as.bigz(sum(am_ba_ratoli[,j]))))))
#}
#storage_cap_ba_ratoli <- mean(storage_caps_ba_ratoli)#,na.rm=T)

ratoli_BA_data <- list()
ratoli_BA_data[["deg_distF_ratoli_BA"]] <- deg.dist_ba_ratoliF
ratoli_BA_data[["deg_distT_ratoli_BA"]] <- deg.dist_ba_ratoliT
ratoli_BA_data[["triad_census_ratoli_BA"]] <- tirad_census_ba_ratoli
ratoli_BA_data[["hub_score_ratoli_BA"]] <- hs_ba_ratoli
ratoli_BA_data[["modularity_ratoli_BA"]] <- modularity_ba_ratoli
ratoli_BA_data[["authority_score_ratoli_BA"]] <- as_ba_ratoli
ratoli_BA_data[["communities_ratoli_BA"]] <- comunitats_ba_ratoli 
ratoli_BA_data[["SW_ratio_ratoli_BA"]] <- SM_Ratio_ba_ratoli
ratoli_BA_data[["cliques_ratoli_BA"]] <- cliques_ba_ratoli 
#ratoli_BA_data[["cost_ratoli_BA"]] <- cost_ba_ratoli
ratoli_BA_data[["nEdges_ratoli_BA"]] <- nEdges_ba_ratoli 
ratoli_BA_data[["nNodes_ratoli_BA"]] <- nNodes_ba_ratoli
#ratoli_BA_data[["rout_eff_ratoli_BA"]] <- rout_eff_ba_ratoli 
#ratoli_BA_data[["storage_cap_ratoli_BA"]] <- storage_cap_ba_ratoli 


save(ratoli_BA_data, file="ratoli_BA_data.RData")

load("ratoli_BA_data.RData")



#--------------------------------------archiu data---------------#
#OP

load("~/BNO/Connectivity_data/ratoli_data.RData")

write(ratoli_data[["nNodes_ratoli"]], file = "networkRatoli.dat", sep = " ")

dadesDeArchiuRatoli <- get.data.frame(gsRatoli)
dadesDeArchiuRatoli$from <- as.numeric(as.factor(dadesDeArchiuRatoli$from))
dadesDeArchiuRatoli$to <- as.numeric(as.factor(dadesDeArchiuRatoli$to))
dadesDeArchiuRatoli$name <- NULL

write.table(dadesDeArchiuRatoli, file = "networkRatoli.dat", sep = " ",append = TRUE,
            col.names = F,
            row.names = FALSE)

#ER

load("~/BNO/Connectivity_data/ratoli_ER_data.RData")

write(ratoli_ER_data[["nNodes_ratoli_ER"]], file = "networkRatoli_ER.dat", sep = " ")

dadesDeArchiudadesDeArchiuRatoli_ER <- get.data.frame(er_Ratoli)
dadesDeArchiudadesDeArchiuRatoli_ER$from <- as.numeric(as.factor(dadesDeArchiudadesDeArchiuRatoli_ER$from))
dadesDeArchiudadesDeArchiuRatoli_ER$to <- as.numeric(as.factor(dadesDeArchiudadesDeArchiuRatoli_ER$to))
dadesDeArchiudadesDeArchiuRatoli_ER$name <- NULL

write.table(dadesDeArchiudadesDeArchiuRatoli_ER, file = "networkRatoli_ER.dat", sep = " ",append = TRUE,
            col.names = F,
            row.names = FALSE)


#SW


load("~/BNO/Connectivity_data/ratoli_SW_data.RData")

write(ratoli_SW_data[["nNodes_ratoli_SW"]], file = "networkRatoli_SW.dat", sep = " ")

dadesDeArchiudadesDeArchiuRatoli_SW <- get.data.frame(sw_ratoli)
dadesDeArchiudadesDeArchiuRatoli_SW$from <- as.numeric(as.factor(dadesDeArchiudadesDeArchiuRatoli_SW$from))
dadesDeArchiudadesDeArchiuRatoli_SW$to <- as.numeric(as.factor(dadesDeArchiudadesDeArchiuRatoli_SW$to))
dadesDeArchiudadesDeArchiuRatoli_SW$name <- NULL

write.table(dadesDeArchiudadesDeArchiuRatoli_SW, file = "networkRatoli_SW.dat", sep = " ",append = TRUE,
            col.names = F,
            row.names = FALSE)

#BA


load("~/BNO/Connectivity_data/ratoli_BA_data.RData")

write(ratoli_BA_data[["nNodes_ratoli_BA"]], file = "networkRatoli_BA.dat", sep = " ")

dadesDeArchiudadesDeArchiuRatoli_BA <- get.data.frame(ba_ratoli)
dadesDeArchiudadesDeArchiuRatoli_BA$from <- as.numeric(as.factor(dadesDeArchiudadesDeArchiuRatoli_BA$from))
dadesDeArchiudadesDeArchiuRatoli_BA$to <- as.numeric(as.factor(dadesDeArchiudadesDeArchiuRatoli_BA$to))
dadesDeArchiudadesDeArchiuRatoli_BA$name <- NULL

write.table(dadesDeArchiudadesDeArchiuRatoli_BA, file = "networkRatoli_BA.dat", sep = " ",append = TRUE,
            col.names = F,
            row.names = FALSE)