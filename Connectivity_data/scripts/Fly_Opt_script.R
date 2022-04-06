library("igraph")
library("plyr")
library(leiden)
library(gmp)
library(Matrix)



archiuFlyOPtic <- read.table("FlyCircuit.csv",sep = ";",row.names=1,header=T)


gFlyOptic <- graph_from_adjacency_matrix(as.matrix(archiuFlyOPtic),mode = "directed",weighted = TRUE)
gFlyOpticS <- simplify(gFlyOptic , remove.multiple = T, remove.loops = F, edge.attr.comb=c(weight="sum", type="ignore") )
write.graph(gFlyOpticS, file = "gFlyOpticSData.ncol",format = "NCOL")

#Degree distribution

degFlyOptic <- degree(gFlyOpticS, mode="all")

deg.degFlyOF <- degree_distribution(gFlyOpticS, cumulative=F, mode="all")


deg.degFlyOT <- degree_distribution(gFlyOpticS, cumulative=T, mode="all")

#triad_census
tirad_census_FlyOptic = triad_census(gFlyOpticS)


#Hubs
hs_FlyOptic <- hub_score(gFlyOpticS, weights=NA)$vector


#Authorities
as_FlyOptic <- authority_score(gFlyOpticS, weights=NA)$vector

#Small-Worldness ratio: 
SM_Ratio_FlyOptic= transitivity(gFlyOpticS)/mean_distance(gFlyOpticS)

#modularity y comunitats
adj_mtx_FlyOptic  <- igraph::as_adjacency_matrix(gFlyOpticS, attr="weight")
camunitats_FlyOptic <- leiden(adj_mtx_FlyOptic)
modularity_FlyOptic <- modularity(gFlyOpticS, camunitats_FlyOptic)


#cliques



#n nodes
nEdges_gFlyOpticS <- gsize(gFlyOpticS)
#cost_gFlyOpticS <- sum(E(gFlyOpticS)$weight)
nNodes_gFlyOpticS <- gorder(gFlyOpticS)

#routing efficiency
#am_gFlyOpticS  <- get.adjacency(gFlyOpticS, attr="weight")
#distances_gFlyOpticS <- distances(gFlyOpticS , mode = "in", weights=1/E(gFlyOpticS)$weight)


#rout_effs_gFlyOpticS <- as.numeric(as.bigz(1/distances_gFlyOpticS)/as.bigz(dim(am_gFlyOpticS)[[1]]*(dim(am_gFlyOpticS)[[1]]-1)))
#rout_effs_matrix_gFlyOpticS <- as.matrix(rout_effs_gFlyOpticS)
#rout_effs_matrix_gFlyOpticS <- rout_effs_matrix_gFlyOpticS[is.finite(rout_effs_matrix_gFlyOpticS)]
#rout_eff_gFlyOpticS <- as.numeric(sum(rout_effs_matrix_gFlyOpticS))

#storage capacity
#storage_caps_gFlyOpticS <- numeric()
#for(j in 1:dim(am_gFlyOpticS)[[1]]){
#  storage_caps_gFlyOpticS <- append(storage_caps_gFlyOpticS, as.integer(2*log2(chooseZ(as.bigz(sum(am_gFlyOpticS[,j]) + nnzero(am_gFlyOpticS[,j]) -1), as.bigz(sum(am_gFlyOpticS[,j]))))))
#}
#storage_cap_gFlyOpticS <- mean(storage_caps_gFlyOpticS)#,na.rm=T)


FlyOP_data <- list()
FlyOP_data[["deg_distF_FlyOP"]] <- deg.degFlyOF
FlyOP_data[["deg_distT_FlyOP"]] <- deg.degFlyOT
FlyOP_data[["triad_census_FlyOP"]] <- tirad_census_FlyOptic
FlyOP_data[["hub_score_FlyOP"]] <- hs_FlyOptic
FlyOP_data[["modularity_FlyOP"]] <- modularity_FlyOptic
FlyOP_data[["authority_score_FlyOP"]] <- as_FlyOptic 
FlyOP_data[["communities_FlyOP"]] <- camunitats_FlyOptic 
FlyOP_data[["SW_ratio_FlyOP"]] <- SM_Ratio_FlyOptic 
#FlyOP_data[["cost_FlyOP"]] <- cost_gFlyOpticS
FlyOP_data[["nEdges_FlyOP"]] <- nEdges_gFlyOpticS
FlyOP_data[["nNodes_FlyOP"]] <- nNodes_gFlyOpticS
#FlyOP_data[["rout_eff_FlyOP"]] <- rout_eff_gFlyOpticS 
#FlyOP_data[["storage_cap_FlyOP"]] <- storage_cap_gFlyOpticS

save(FlyOP_data, file="FlyOP_data.RData")

load("FlyOP_data.RData")


#------------------------------------------------------ Erdos-Renyi random graph model -------------------------------------#


set.seed(2706)

er_FlyO <- sample_gnm(n=length(unique(row.names(archiuFlyOPtic))), m=length(E(gFlyOpticS)),directed=T)
V(er_FlyO)$name <- V(gFlyOpticS)$name
E(er_FlyO)$weight <- E(gFlyOpticS)$weight

deg_erFlyO <- degree(er_FlyO, mode="all")



deg.dist_erFlyOF <- degree_distribution(er_FlyO, cumulative=F, mode="all")


deg.dist_erFlyOT <- degree_distribution(er_FlyO, cumulative=T, mode="all")


tirad_census_er_FlyO = triad_census(er_FlyO)

#Hubs
hs_er_FlyO <- hub_score(er_FlyO, weights=NA)$vector

#Authorities
as_er_FlyO <- authority_score(er_FlyO, weights=NA)$vector

#Small-Worldness ratio: 
SM_Ratio_er_FlyO= transitivity(er_FlyO)/mean_distance(er_FlyO)


#modularity y comunitats
adj_mtx_er_FlyO  <- igraph::as_adjacency_matrix(er_FlyO, attr="weight")
camunitats_er_FlyO <- leiden(adj_mtx_er_FlyO)
modularity_er_FlyO <- modularity(er_FlyO, camunitats_er_FlyO)


#cliques

#n nodes
nEdges_er_FlyO <- gsize(er_FlyO)
#cost_er_FlyO <- sum(E(er_FlyO)$weight)
nNodes_er_FlyO <- gorder(er_FlyO)

#routing efficiency
#am_er_FlyO  <- get.adjacency(er_FlyO, attr="weight")
#distances_er_FlyO <- distances(er_FlyO , mode = "in", weights=1/E(er_FlyO)$weight)


#rout_effs_er_FlyO <- as.numeric(as.bigz(1/distances_er_FlyO)/as.bigz(dim(am_er_FlyO)[[1]]*(dim(am_er_FlyO)[[1]]-1)))
#rout_effs_matrix_er_FlyO <- as.matrix(rout_effs_er_FlyO)
#rout_effs_matrix_er_FlyO <- rout_effs_matrix_er_FlyO[is.finite(rout_effs_matrix_er_FlyO)]
#rout_eff_er_FlyO <- as.numeric(sum(rout_effs_matrix_er_FlyO))

#storage capacity
#storage_caps_er_FlyO <- numeric()
#for(j in 1:dim(am_er_FlyO)[[1]]){
  #storage_caps_er_FlyO <- append(storage_caps_er_FlyO, as.integer(2*log2(chooseZ(as.bigz(sum(am_er_FlyO[,j]) + nnzero(am_er_FlyO[,j]) -1), as.bigz(sum(am_er_FlyO[,j]))))))
  #}
#storage_cap_er_FlyO <- mean(storage_caps_er_FlyO)#,na.rm=T)


FlyOP_ER_data <- list()
FlyOP_ER_data[["deg_distF_FlyOP_ER"]] <- deg.dist_erFlyOF
FlyOP_ER_data[["deg_distT_FlyOP_ER"]] <- deg.dist_erFlyOT
FlyOP_ER_data[["triad_census_FlyOP_ER"]] <- tirad_census_er_FlyO
FlyOP_ER_data[["hub_score_FlyOP_ER"]] <- hs_er_FlyO
FlyOP_ER_data[["modularity_FlyOP_ER"]] <- modularity_er_FlyO
FlyOP_ER_data[["authority_score_FlyOP_ER"]] <- as_er_FlyO 
FlyOP_ER_data[["communities_FlyOP_ER"]] <- camunitats_er_FlyO 
FlyOP_ER_data[["SW_ratio_FlyOP_ER"]] <- SM_Ratio_er_FlyO 
#FlyOP_ER_data[["cost_FlyOP_ER"]] <- cost_er_FlyO
FlyOP_ER_data[["nEdges_FlyOP_ER"]] <- nEdges_er_FlyO
FlyOP_ER_data[["nNodes_FlyOP_ER"]] <- nNodes_er_FlyO
#FlyOP_ER_data[["rout_eff_FlyOP_ER"]] <- rout_eff_er_FlyO 
#FlyOP_ER_data[["storage_cap_FlyOP_ER"]] <- storage_cap_er_FlyO


save(FlyOP_ER_data, file="FlyOP_ER_data.RData")

load("FlyOP_ER_data.RData")



#------------------------------------------ Watts-Strogatz small-world model--------------------------------------------#


#sw_FlyO <- sample_smallworld(dim=1, size=length(unique(row.names(archiuFlyOPtic))), dir=TRUE, nei=11, p=0.1)
sw_FlyO <- rewire(graph.lattice(length(unique(row.names(archiuFlyOPtic))), nei=20, dir=TRUE, mutual=TRUE), 
                    with = each_edge(p = .02, loops = FALSE))
#plot(sw_Fly, vertex.size=6, vertex.label=NA, layout=layout_in_circle)
V(sw_FlyO)$name <- V(gFlyOpticS)$name

pesos_sw <- sample(E(gFlyOpticS)$weight,length(E(sw_FlyO)))
E(sw_FlyO)$weight <- pesos_sw


deg_swFlyO <- degree(sw_FlyO, mode="all")


deg.dist_swFlyFO <- degree_distribution(sw_FlyO, cumulative=F, mode="all")


deg.dist_swFlyTO <- degree_distribution(sw_FlyO, cumulative=T, mode="all")


tirad_census_sw_FlyO = triad_census(sw_FlyO)

#Hubs
hs_sw_FlyO <- hub_score(sw_FlyO, weights=NA)$vector

#Authorities
as_sw_FlyO <- authority_score(sw_FlyO, weights=NA)$vector

#Small-Worldness ratio: 
SM_Ratio_sw_FlyO= transitivity(sw_FlyO)/mean_distance(sw_FlyO)

#modularity y comunitats
adj_mtx_sw_FlyO  <- igraph::as_adjacency_matrix(sw_FlyO, attr="weight")
camunitats_sw_FlyO <- leiden(adj_mtx_sw_FlyO)
modularity_sw_FlyO <- modularity(sw_FlyO, camunitats_sw_FlyO)


#cliques

#n nodes
nEdges_sw_FlyO <- gsize(sw_FlyO)
#cost_sw_FlyO <- sum(E(sw_FlyO)$weight)
nNodes_sw_FlyO <- gorder(sw_FlyO)

#routing efficiency
#am_sw_FlyO  <- get.adjacency(sw_FlyO, attr="weight")
#distances_sw_FlyO <- distances(sw_FlyO , mode = "in", weights=1/E(sw_FlyO)$weight)


#rout_effs_sw_FlyO <- as.numeric(as.bigz(1/distances_sw_FlyO)/as.bigz(dim(am_sw_FlyO)[[1]]*(dim(am_sw_FlyO)[[1]]-1)))
#rout_effs_matrix_sw_FlyO <- as.matrix(rout_effs_sw_FlyO)
#rout_effs_matrix_sw_FlyO <- rout_effs_matrix_sw_FlyO[is.finite(rout_effs_matrix_sw_FlyO)]
#rout_eff_sw_FlyO <- as.numeric(sum(rout_effs_matrix_sw_FlyO))

#storage capacity
#storage_caps_sw_FlyO <- numeric()
#for(j in 1:dim(am_sw_FlyO)[[1]]){
  #storage_caps_sw_FlyO <- append(storage_caps_sw_FlyO, as.integer(2*log2(chooseZ(as.bigz(sum(am_sw_FlyO[,j]) + nnzero(am_sw_FlyO[,j]) -1), as.bigz(sum(am_sw_FlyO[,j]))))))
  #}
#storage_cap_sw_FlyO <- mean(storage_caps_sw_FlyO)#,na.rm=T)



FlyOP_SW_data <- list()
FlyOP_SW_data[["deg_distF_FlyOP_SW"]] <- deg.dist_erFlyOF
FlyOP_SW_data[["deg_distT_FlyOP_SW"]] <- deg.dist_erFlyOT
FlyOP_SW_data[["triad_census_FlyOP_SW"]] <- tirad_census_sw_FlyO
FlyOP_SW_data[["hub_score_FlyOP_SW"]] <- hs_sw_FlyO
FlyOP_SW_data[["modularity_FlyOP_SW"]] <- modularity_sw_FlyO
FlyOP_SW_data[["authority_score_FlyOP_SW"]] <- as_sw_FlyO 
FlyOP_SW_data[["communities_FlyOP_SW"]] <- camunitats_sw_FlyO
FlyOP_SW_data[["SW_ratio_FlyOP_SW"]] <- SM_Ratio_sw_FlyO 
#FlyOP_SW_data[["cost_FlyOP_SW"]] <- cost_sw_FlyO
FlyOP_SW_data[["nEdges_FlyOP_SW"]] <- nEdges_sw_FlyO
FlyOP_SW_data[["nNodes_FlyOP_SW"]] <- nNodes_sw_FlyO
FlyOP_SW_data[["rout_eff_FlyOP_SW"]] <- 166.249834
#FlyOP_SW_data[["storage_cap_FlyOP_SW"]] <- storage_cap_sw_FlyO 


save(FlyOP_SW_data, file="FlyOP_SW_data.RData")

load("FlyOP_SW_data.RData")




#----------------------------- Barabasi-Albert preferential attachment model for scale-free graphs -------------------------#


ba_flyO <-  sample_pa(n=length(unique(row.names(archiuFlyOPtic))), power=1, m=11,  directed=T)
#plot(ba_fly, vertex.size=6, vertex.label=NA, layout=layout_in_circle)
V(ba_flyO)$name <- V(gFlyOpticS)$name

pesos_Ba <- sample(E(gFlyOpticS)$weight,length(E(ba_flyO)))
E(ba_flyO)$weight <- pesos_Ba

deg_baFlyO <- degree(ba_flyO, mode="all")


deg.dist_baFlyFO <- degree_distribution(ba_flyO, cumulative=F, mode="all")


deg.dist_baFlyTO <- degree_distribution(ba_flyO, cumulative=T, mode="all")


tirad_census_ba_flyO = triad_census(ba_flyO)

#Hubs
hs_ba_flyO <- hub_score(ba_flyO, weights=NA)$vector

#Authorities
as_ba_flyO <- authority_score(ba_flyO, weights=NA)$vector

#Small-Worldness ratio: 
SM_Ratio_ba_flyO= transitivity(ba_flyO)/mean_distance(ba_flyO)

#modularity y comunitats
adj_mtx_ba_flyO  <- igraph::as_adjacency_matrix(ba_flyO, attr="weight")
camunitats_ba_flyO <- leiden(adj_mtx_ba_flyO)
modularity_ba_flyO <- modularity(ba_flyO, camunitats_ba_flyO)


#cliques

#n nodes
nEdges_ba_flyO <- gsize(ba_flyO)
#cost_ba_flyO <- sum(E(ba_flyO)$weight)
nNodes_ba_flyO <- gorder(ba_flyO)

#routing efficiency
#am_ba_flyO  <- get.adjacency(ba_flyO, attr="weight")
#distances_ba_flyO <- distances(ba_flyO , mode = "in", weights=1/E(ba_flyO)$weight)


#rout_effs_ba_flyO <- as.numeric(as.bigz(1/distances_ba_flyO)/as.bigz(dim(am_ba_flyO)[[1]]*(dim(am_ba_flyO)[[1]]-1)))
#rout_effs_matrix_ba_flyO <- as.matrix(rout_effs_ba_flyO)
#rout_effs_matrix_ba_flyO <- rout_effs_matrix_ba_flyO[is.finite(rout_effs_matrix_ba_flyO)]
#rout_eff_ba_flyO <- as.numeric(sum(rout_effs_matrix_ba_flyO))

#storage capacity
#storage_caps_ba_flyO <- numeric()
#for(j in 1:dim(am_ba_flyO)[[1]]){
  #storage_caps_ba_flyO <- append(storage_caps_ba_flyO, as.integer(2*log2(chooseZ(as.bigz(sum(am_ba_flyO[,j]) + nnzero(am_ba_flyO[,j]) -1), as.bigz(sum(am_ba_flyO[,j]))))))
  #}
#storage_cap_ba_flyO <- mean(storage_caps_ba_flyO)#,na.rm=T)


FlyOP_BA_data <- list()
FlyOP_BA_data[["deg_distF_FlyOP_BA"]] <- deg.dist_baFlyFO
FlyOP_BA_data[["deg_distT_FlyOP_BA"]] <- deg.dist_baFlyTO
FlyOP_BA_data[["triad_census_FlyOP_BA"]] <- tirad_census_ba_flyO
FlyOP_BA_data[["hub_score_FlyOP_BA"]] <- hs_ba_flyO
FlyOP_BA_data[["modularity_FlyOP_BA"]] <- modularity_ba_flyO
FlyOP_BA_data[["authority_score_FlyOP_BA"]] <- as_ba_flyO
FlyOP_BA_data[["communities_FlyOP_BA"]] <- camunitats_ba_flyO
FlyOP_BA_data[["SW_ratio_FlyOP_BA"]] <- SM_Ratio_ba_flyO 
#FlyOP_BA_data[["cost_FlyOP_BA"]] <- cost_ba_flyO
FlyOP_BA_data[["nEdges_FlyOP_BA"]] <- nEdges_ba_flyO
FlyOP_BA_data[["nNodes_FlyOP_BA"]] <- nNodes_ba_flyO
FlyOP_BA_data[["rout_eff_FlyOP_BA"]] <- 58.095561
#FlyOP_BA_data[["storage_cap_FlyOP_BA"]] <- storage_cap_ba_flyO 


save(FlyOP_BA_data, file="FlyOP_BA_data.RData")

load("FlyOP_BA_data.RData")


#--------------------------------------archiu data---------------#
# FlyOP

load("~/BNO/Connectivity_data/FlyOP_data.RData")

write(FlyOP_data[["nNodes_FlyOP"]], file = "networkFlyOP.dat", sep = " ")

dadesDeArchiuFlyOp <- get.data.frame(gFlyOpticS)
dadesDeArchiuFlyOp$from <- as.numeric(as.factor(dadesDeArchiuFlyOp$from))
dadesDeArchiuFlyOp$to <- as.numeric(as.factor(dadesDeArchiuFlyOp$to))
dadesDeArchiuFlyOp$name <- NULL

write.table(dadesDeArchiuFlyOp, file = "networkFlyOP.dat", sep = " ",append = TRUE,
            col.names = F,
            row.names = FALSE)

# Fly ER

load("~/BNO/Connectivity_data/FlyOP_ER_data.RData")

write(FlyOP_ER_data[["nNodes_FlyOP_ER"]], file = "networkFlyOP_ER.dat", sep = " ")

dadesDeArchiuFlyOp_ER <- get.data.frame(er_FlyO)
dadesDeArchiuFlyOp_ER$from <- as.numeric(as.factor(dadesDeArchiuFlyOp_ER$from))
dadesDeArchiuFlyOp_ER$to <- as.numeric(as.factor(dadesDeArchiuFlyOp_ER$to))
dadesDeArchiuFlyOp_ER$name <- NULL

write.table(dadesDeArchiuFlyOp_ER, file = "networkFlyOP_ER.dat", sep = " ",append = TRUE,
            col.names = F,
            row.names = FALSE)


# Fly SW


load("~/BNO/Connectivity_data/FlyOP_SW_data.RData")

write(FlyOP_SW_data[["nNodes_FlyOP_SW"]], file = "networkFlyOP_SW.dat", sep = " ")

dadesDeArchiuFlyOp_SW <- get.data.frame(sw_FlyO)
dadesDeArchiuFlyOp_SW$from <- as.numeric(as.factor(dadesDeArchiuFlyOp_SW$from))
dadesDeArchiuFlyOp_SW$to <- as.numeric(as.factor(dadesDeArchiuFlyOp_SW$to))
dadesDeArchiuFlyOp_SW$name <- NULL

write.table(dadesDeArchiuFlyOp_SW, file = "networkFlyOP_SW.dat", sep = " ",append = TRUE,
            col.names = F,
            row.names = FALSE)

# Fly BA


load("~/BNO/Connectivity_data/FlyOP_BA_data.RData")

write(FlyOP_BA_data[["nNodes_FlyOP_BA"]], file = "networkFlyOP_BA.dat", sep = " ")

dadesDeArchiuFlyOp_BA <- get.data.frame(ba_flyO)
dadesDeArchiuFlyOp_BA$from <- as.numeric(as.factor(dadesDeArchiuFlyOp_BA$from))
dadesDeArchiuFlyOp_BA$to <- as.numeric(as.factor(dadesDeArchiuFlyOp_BA$to))
dadesDeArchiuFlyOp_BA$name <- NULL

write.table(dadesDeArchiuFlyOp_BA, file = "networkFlyOP_BA.dat", sep = " ",append = TRUE,
            col.names = F,
            row.names = FALSE)

