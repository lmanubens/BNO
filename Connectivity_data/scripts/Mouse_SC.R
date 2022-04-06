library("igraph")
library("plyr")
library(leiden)
library(gmp)
library(Matrix)


archiuDadesMouse <- read.csv("affinity_map_nonflip_thres30_20211103.csv",sep = ",")
colnames(archiuDadesMouse)[2:13603] <- gsub("X","",colnames(archiuDadesMouse)[2:13603]) 

rownames(archiuDadesMouse) <- archiuDadesMouse[,1]

archiuDadesMouse$X <- NULL

#Afegir files que tinguin nom de columnes que no hi siguin i assignarlos 0s
archiuDadesMouse <- as.matrix(archiuDadesMouse)

matriu0s <- matrix(0,nrow=sum(!colnames(archiuDadesMouse) %in% rownames(archiuDadesMouse)),ncol=ncol(archiuDadesMouse))

archiuDadesMouse <- rbind(archiuDadesMouse,matriu0s)

rownames(archiuDadesMouse) <- c(rownames(archiuDadesMouse)[1:2743],colnames(archiuDadesMouse)[!colnames(archiuDadesMouse) %in% rownames(archiuDadesMouse)])



gMouseSC <-graph_from_adjacency_matrix(
  archiuDadesMouse,
  mode = "directed",
  weighted = TRUE,
  diag = TRUE,
  add.colnames = NULL,
  add.rownames = NULL
)

gsMouseSC <- simplify( gMouseSC , remove.multiple = T, remove.loops = F, edge.attr.comb=c(weight="sum", type="ignore") )


deg.deggMouseSCF <- degree_distribution(gsMouseSC, cumulative=F, mode="all")

deg.deggMouseSCT <- degree_distribution(gsMouseSC, cumulative=T, mode="all")

#triad_census
tirad_census_gMouseSC = triad_census(gsMouseSC)


#Hubs
hs_gMouseSC <- hub_score(gsMouseSC, weights=NA)$vector


#Authorities
as_gMouseSC <- authority_score(gsMouseSC, weights=NA)$vector

#Small-Worldness ratio: 
SM_Ratio_gMouseSC= transitivity(gsMouseSC)/mean_distance(gsMouseSC)

#modularity y comunitats
adj_mtx_gMouseSC  <- igraph::as_adjacency_matrix(gsMouseSC)
camunitats_gMouseSC <- leiden(adj_mtx_gMouseSC)
modularity_gMouseSC <- modularity(gsMouseSC, camunitats_gMouseSC)


#cliques



#n nodes
nEdges_gMouseSC <- gsize(gsMouseSC)
#cost_gMouseSC <- sum(E(gsMouseSC)$weight)
nNodes_gMouseSC <- gorder(gsMouseSC)

#routing efficiency
#am_gMouseSC  <- get.adjacency(gsMouseSC)

#distances_gMouseSC <- distances(gsMouseSC , mode = "in", weights=1/E(gsMouseSC)$weight)


#rout_effs_gMouseSC <- as.numeric(as.bigz(1/distances_gMouseSC)/as.bigz(dim(am_gMouseSC)[[1]]*(dim(am_gMouseSC)[[1]]-1)))
#rout_effs_matrix_gMouseSC <- as.matrix(rout_effs_gMouseSC)
#rout_effs_matrix_gMouseSC <- rout_effs_matrix_gMouseSC[is.finite(rout_effs_matrix_gMouseSC)]
#rout_eff_gMouseSC <- as.numeric(sum(rout_effs_matrix_gMouseSC))

#storage capacity
#storage_caps_gMouseSC <- numeric()
#for(j in 1:dim(am_gMouseSC)[[1]]){
#  storage_caps_gMouseSC <- append(storage_caps_gMouseSC, as.integer(2*log2(chooseZ(as.bigz(sum(am_gMouseSC[,j]) + nnzero(am_gMouseSC[,j]) -1), as.bigz(sum(am_gMouseSC[,j]))))))
#}
#storage_cap_gMouseSC <- mean(storage_caps_gMouseSC)#,na.rm=T)


MouseSC_data <- list()
MouseSC_data[["deg_distF_MouseSC"]] <- deg.deggMouseSCF
MouseSC_data[["deg_distT_MouseSC"]] <- deg.deggMouseSCT
MouseSC_data[["triad_census_MouseSC"]] <- tirad_census_gMouseSC
MouseSC_data[["hub_score_MouseSC"]] <- hs_gMouseSC
MouseSC_data[["modularity_MouseSC"]] <- modularity_gMouseSC
MouseSC_data[["authority_score_MouseSC"]] <- as_gMouseSC 
MouseSC_data[["communities_MouseSC"]] <- camunitats_gMouseSC
MouseSC_data[["SW_ratio_MouseSC"]] <- SM_Ratio_gMouseSC
#MouseSC_data[["cost_MouseSC"]] <- cost_gMouseSC
MouseSC_data[["nEdges_MouseSC"]] <- nEdges_gMouseSC
MouseSC_data[["nNodes_MouseSC"]] <- nNodes_gMouseSC
#MouseSC_data[["rout_eff_MouseSC"]] <- rout_eff_gMouseSC 
#MouseSC_data[["storage_cap_MouseSC"]] <- storage_cap_gMouseSC

save(MouseSC_data, file="MouseSC_data.RData")

load("MouseSC_data.RData")
#------------------------------------------------------ Erdos-Renyi random graph model -------------------------------------#

set.seed(2706)

er_MouseSC <- sample_gnm(n=nNodes_gMouseSC, m=length(E(gsMouseSC)),directed=T)
V(er_MouseSC)$name <- V(gsMouseSC)$name
E(er_MouseSC)$weight <- E(gsMouseSC)$weight

deg_er_MouseSC <- degree(er_MouseSC, mode="all")



deg.dist_er_MouseSCF <- degree_distribution(er_MouseSC, cumulative=F, mode="all")


deg.dist_er_MouseSCT <- degree_distribution(er_MouseSC, cumulative=T, mode="all")


tirad_census_er_MouseSC = triad_census(er_MouseSC)

#Hubs
hs_er_MouseSC <- hub_score(er_MouseSC, weights=NA)$vector

#Authorities
as_er_MouseSC <- authority_score(er_MouseSC, weights=NA)$vector

#Small-Worldness ratio: 
SM_Ratio_er_MouseSC= transitivity(er_MouseSC)/mean_distance(er_MouseSC)


#modularity y comunitats
adj_mtx_er_MouseSC  <- igraph::as_adjacency_matrix(er_MouseSC, attr="weight")
camunitats_er_MouseSC <- leiden(adj_mtx_er_MouseSC)
modularity_er_MouseSC <- modularity(er_MouseSC, camunitats_er_MouseSC)


#cliques

#n nodes
nEdges_er_MouseSC <- gsize(er_MouseSC)
#cost_er_MouseSC <- sum(E(er_MouseSC)$weight)
nNodes_er_MouseSC <- gorder(er_MouseSC)

#routing efficiency
#am_er_MouseSC  <- get.adjacency(er_MouseSC, attr="weight")
#distances_er_MouseSC <- distances(er_MouseSC , mode = "in", weights=1/E(er_MouseSC)$weight)


#rout_effs_er_MouseSC <- as.numeric(as.bigz(1/distances_er_MouseSC)/as.bigz(dim(am_er_MouseSC)[[1]]*(dim(am_er_MouseSC)[[1]]-1)))
#rout_effs_matrix_er_MouseSC <- as.matrix(rout_effs_er_MouseSC)
#rout_effs_matrix_er_MouseSC <- rout_effs_matrix_er_MouseSC[is.finite(rout_effs_matrix_er_MouseSC)]
#rout_eff_er_MouseSC <- as.numeric(sum(rout_effs_matrix_er_MouseSC))

#storage capacity
#storage_caps_er_MouseSC <- numeric()
#for(j in 1:dim(am_er_MouseSC)[[1]]){
#  storage_caps_er_MouseSC <- append(storage_caps_er_MouseSC, as.integer(2*log2(chooseZ(as.bigz(sum(am_er_MouseSC[,j]) + nnzero(am_er_MouseSC[,j]) -1), as.bigz(sum(am_er_MouseSC[,j]))))))
#}
#storage_cap_er_MouseSC <- mean(storage_caps_er_MouseSC)#,na.rm=T)


MouseSC_ER_data <- list()
MouseSC_ER_data[["deg_distF_MouseSC_ER"]] <- deg.dist_er_MouseSCF
MouseSC_ER_data[["deg_distT_MouseSC_ER"]] <- deg.dist_er_MouseSCT
MouseSC_ER_data[["triad_census_MouseSC_ER"]] <- tirad_census_er_MouseSC
MouseSC_ER_data[["hub_score_MouseSC_ER"]] <- hs_er_MouseSC
MouseSC_ER_data[["modularity_MouseSC_ER"]] <- modularity_er_MouseSC
MouseSC_ER_data[["authority_score_MouseSC_ER"]] <- as_er_MouseSC 
MouseSC_ER_data[["communities_MouseSC_ER"]] <- camunitats_er_MouseSC 
MouseSC_ER_data[["SW_ratio_MouseSC_ER"]] <- SM_Ratio_er_MouseSC 
#MouseSC_ER_data[["cost_MouseSC_ER"]] <- cost_er_MouseSC
MouseSC_ER_data[["nEdges_MouseSC_ER"]] <- nEdges_er_MouseSC
MouseSC_ER_data[["nNodes_MouseSC_ER"]] <- nNodes_er_MouseSC
#MouseSC_ER_data[["rout_eff_MouseSC_ER"]] <- rout_eff_er_MouseSC 
#MouseSC_ER_data[["storage_cap_MouseSC_ER"]] <- storage_cap_er_MouseSC


save(MouseSC_ER_data, file="MouseSC_ER_data.RData")

load("MouseSC_ER_data.RData")



#------------------------------------------------ Watts-Strogatz small-world model--------------------------------------------#



sw_MouseSC <- rewire(graph.lattice(nNodes_gMouseSC, nei=4, dir=TRUE, mutual=TRUE), 
                     with = each_edge(p = .02, loops = FALSE))

pesos_ba_sw_MouseSC <- sample(E(gsMouseSC)$weight,length(E(sw_MouseSC)))
E(sw_MouseSC)$weight <- pesos_ba_sw_MouseSC


deg_sw_MouseSC <- degree(sw_MouseSC, mode="all")


deg.dist_sw_MouseSCF <- degree_distribution(sw_MouseSC, cumulative=F, mode="all")


deg.dist_sw_MouseSCT <- degree_distribution(sw_MouseSC, cumulative=T, mode="all")


tirad_census_sw_MouseSC = triad_census(sw_MouseSC)

#Hubs
hs_sw_MouseSC<- hub_score(sw_MouseSC, weights=NA)$vector

#Authorities
as_sw_MouseSC <- authority_score(sw_MouseSC, weights=NA)$vector

#Small-Worldness ratio: 
SM_Ratio_sw_MouseSC = transitivity(sw_MouseSC)/mean_distance(sw_MouseSC)

#modularity y comunitats
adj_mtx_sw_MouseSC  <- igraph::as_adjacency_matrix(sw_MouseSC, attr="weight")
camunitats_sw_MouseSC <- leiden(adj_mtx_sw_MouseSC)
modularity_sw_MouseSC <- modularity(sw_MouseSC, camunitats_sw_MouseSC)


#cliques

#n nodes
nEdges_sw_MouseSC <- gsize(sw_MouseSC)
#cost_sw_MouseSC <- sum(E(sw_MouseSC)$weight)
nNodes_sw_MouseSC <- gorder(sw_MouseSC)

#routing efficiency
#am_sw_MouseSC  <- get.adjacency(sw_MouseSC, attr="weight")
#distances_sw_MouseSC <- distances(sw_MouseSC , mode = "in", weights=1/E(sw_MouseSC)$weight)


#rout_effs_sw_MouseSC <- as.numeric(as.bigz(1/distances_sw_MouseSC)/as.bigz(dim(am_sw_MouseSC)[[1]]*(dim(am_sw_MouseSC)[[1]]-1)))
#rout_effs_matrix_sw_MouseSC <- as.matrix(rout_effs_sw_MouseSC)
#rout_effs_matrix_sw_MouseSC <- rout_effs_matrix_sw_MouseSC[is.finite(rout_effs_matrix_sw_MouseSC)]
#rout_eff_sw_MouseSC <- as.numeric(sum(rout_effs_matrix_sw_MouseSC))

#storage capacity
#storage_caps_sw_MouseSC <- numeric()
#for(j in 1:dim(am_sw_MouseSC)[[1]]){
#  storage_caps_sw_MouseSC <- append(storage_caps_sw_MouseSC, as.integer(2*log2(chooseZ(as.bigz(sum(am_sw_MouseSC[,j]) + nnzero(am_sw_MouseSC[,j]) -1), as.bigz(sum(am_sw_MouseSC[,j]))))))
#}
#storage_cap_sw_MouseSC <- mean(storage_caps_sw_MouseSC)#,na.rm=T)



MouseSC_SW_data <- list()
MouseSC_SW_data[["deg_distF_MouseSC_SW"]] <- deg.dist_sw_MouseSCF
MouseSC_SW_data[["deg_distT_MouseSC_SW"]] <- deg.dist_sw_MouseSCT
MouseSC_SW_data[["triad_census_MouseSC_SW"]] <- tirad_census_sw_MouseSC
MouseSC_SW_data[["hub_score_MouseSC_SW"]] <- hs_sw_MouseSC
MouseSC_SW_data[["modularity_MouseSC_SW"]] <- modularity_sw_MouseSC
MouseSC_SW_data[["authority_score_MouseSC_SW"]] <- as_sw_MouseSC 
MouseSC_SW_data[["communities_MouseSC_SW"]] <- camunitats_sw_MouseSC
MouseSC_SW_data[["SW_ratio_MouseSC_SW"]] <- SM_Ratio_sw_MouseSC 
#MouseSC_SW_data[["cost_MouseSC_SW"]] <- cost_sw_MouseSC
MouseSC_SW_data[["nEdges_MouseSC_SW"]] <- nEdges_sw_MouseSC
MouseSC_SW_data[["nNodes_MouseSC_SW"]] <- nNodes_sw_MouseSC
#MouseSC_SW_data[["rout_eff_MouseSC_SW"]] <- rout_eff_sw_MouseSC
#MouseSC_SW_data[["storage_cap_MouseSC_SW"]] <- storage_cap_sw_MouseSC 


save(MouseSC_SW_data, file="MouseSC_SW_data.RData")

load("MouseSC_SW_data.RData")

#----------------------------- Barabasi-Albert preferential attachment model for scale-free graphs -------------------------#

ba_MouseSC <-  sample_pa(n=nNodes_gMouseSC, power=1,  m=5,  directed=T)
pesos_ba_MouseSC <- sample(E(gsMouseSC)$weight,length(E(ba_MouseSC)))
E(ba_MouseSC)$weight <- pesos_ba_MouseSC



deg_ba_MouseSC <- degree(ba_MouseSC, mode="all")

deg.dist_ba_MouseSCF <- degree_distribution(ba_MouseSC, cumulative=F, mode="all")

deg.dist_ba_MouseSCT <- degree_distribution(ba_MouseSC, cumulative=T, mode="all")

tirad_census_ba_MouseSC = triad_census(ba_MouseSC)

#Hubs
hs_ba_MouseSC <- hub_score(ba_MouseSC, weights=NA)$vector

#Authorities
as_ba_MouseSC <- authority_score(ba_MouseSC, weights=NA)$vector

#Small-Worldness ratio: 
SM_Ratio_ba_MouseSC = transitivity(ba_MouseSC)/mean_distance(ba_MouseSC)

#modularity y comunitats
adj_mtx_ba_MouseSC  <- igraph::as_adjacency_matrix(ba_MouseSC, attr="weight")
comunitats_ba_MouseSC <- leiden(adj_mtx_ba_MouseSC)
modularity_ba_MouseSC <- modularity(ba_MouseSC, comunitats_ba_MouseSC)


#cliques


#n nodes
nEdges_ba_MouseSC <- gsize(ba_MouseSC)
#cost_ba_MouseSC <- sum(E(ba_MouseSC)$weight)
nNodes_ba_MouseSC <- gorder(ba_MouseSC)

#routing efficiency
#am_ba_MouseSC <- get.adjacency(ba_MouseSC, attr="weight")
#distances_ba_MouseSC <- distances(ba_MouseSC , mode = "in", weights=1/E(ba_MouseSC)$weight)


#rout_effs_ba_MouseSC <- as.numeric(as.bigz(1/distances_ba_MouseSC)/as.bigz(dim(am_ba_MouseSC)[[1]]*(dim(am_ba_MouseSC)[[1]]-1)))
#rout_effs_matrix_ba_MouseSC <- as.matrix(rout_effs_ba_MouseSC)
#rout_effs_matrix_ba_MouseSC <- rout_effs_matrix_ba_MouseSC[is.finite(rout_effs_matrix_ba_MouseSC)]
#rout_eff_ba_MouseSC <- as.numeric(sum(rout_effs_matrix_ba_MouseSC))

#storage capacity
#storage_caps_ba_MouseSC <- numeric()
#for(j in 1:dim(am_ba_MouseSC)[[1]]){
#  storage_caps_ba_MouseSC <- append(storage_caps_ba_MouseSC, as.integer(2*log2(chooseZ(as.bigz(sum(am_ba_MouseSC[,j]) + nnzero(am_ba_MouseSC[,j]) -1), as.bigz(sum(am_ba_MouseSC[,j]))))))
#}
#storage_cap_ba_MouseSC <- mean(storage_caps_ba_MouseSC)#,na.rm=T)

MouseSC_BA_data <- list()
MouseSC_BA_data[["deg_distF_MouseSC_BA"]] <- deg.dist_ba_MouseSCF
MouseSC_BA_data[["deg_distT_MouseSC_BA"]] <- deg.dist_ba_MouseSCT
MouseSC_BA_data[["triad_census_MouseSC_BA"]] <- tirad_census_ba_MouseSC
MouseSC_BA_data[["hub_score_MouseSC_BA"]] <- hs_ba_MouseSC
MouseSC_BA_data[["modularity_MouseSC_BA"]] <- modularity_ba_MouseSC
MouseSC_BA_data[["authority_score_MouseSC_BA"]] <- as_ba_MouseSC
MouseSC_BA_data[["communities_MouseSC_BA"]] <- comunitats_ba_MouseSC 
MouseSC_BA_data[["SW_ratio_MouseSC_BA"]] <- SM_Ratio_ba_MouseSC
MouseSC_BA_data[["cliques_MouseSC_BA"]] <- cliques_ba_MouseSC 
#MouseSC_BA_data[["cost_MouseSC_BA"]] <- cost_ba_MouseSC
MouseSC_BA_data[["nEdges_MouseSC_BA"]] <- nEdges_ba_MouseSC 
MouseSC_BA_data[["nNodes_MouseSC_BA"]] <- nNodes_ba_MouseSC
#MouseSC_BA_data[["rout_eff_MouseSC_BA"]] <- rout_eff_ba_MouseSC 
#MouseSC_BA_data[["storage_cap_MouseSC_BA"]] <- storage_cap_ba_MouseSC 


save(MouseSC_BA_data, file="MouseSC_BA_data.RData")

load("MouseSC_BA_data.RData")



#--------------------------------------archiu data---------------#
# Mouse_SC OP

load("~/BNO/Connectivity_data/MouseSC_data.RData")

write(MouseSC_data[["nNodes_MouseSC"]], file = "networkMouseSC.dat", sep = " ")

dadesDeMouseSC <- get.data.frame(gsMouseSC)
dadesDeMouseSC$from <- as.numeric(as.factor(dadesDeMouseSC$from))
dadesDeMouseSC$to <- as.numeric(as.factor(dadesDeMouseSC$to))
dadesDeMouseSC$name <- NULL

write.table(dadesDeMouseSC, file = "networkMouseSC.dat", sep = " ",append = TRUE,
            col.names = F,
            row.names = FALSE)



# Mouse_SC ER

load("~/BNO/Connectivity_data/MouseSC_ER_data.RData")

write(MouseSC_ER_data[["nNodes_MouseSC_ER"]], file = "networkMouseSC_ER.dat", sep = " ")

dadesDeMouseSC_ER <- get.data.frame(er_MouseSC)
dadesDeMouseSC_ER$from <- as.numeric(as.factor(dadesDeMouseSC_ER$from))
dadesDeMouseSC_ER$to <- as.numeric(as.factor(dadesDeMouseSC_ER$to))
dadesDeMouseSC_ER$name <- NULL

write.table(dadesDeMouseSC_ER, file = "networkMouseSC_ER.dat", sep = " ",append = TRUE,
            col.names = F,
            row.names = FALSE)


# Mouse_SC SW


load("~/BNO/Connectivity_data/MouseSC_SW_data.RData")

write(MouseSC_SW_data[["nNodes_MouseSC_SW"]], file = "networkMouseSC_SW.dat", sep = " ")

dadesDeMouseSC_SW <- get.data.frame(sw_MouseSC)
dadesDeMouseSC_SW$from <- as.numeric(as.factor(dadesDeMouseSC_SW$from))
dadesDeMouseSC_SW$to <- as.numeric(as.factor(dadesDeMouseSC_SW$to))
dadesDeMouseSC_SW$name <- NULL

write.table(dadesDeMouseSC_SW, file = "networkMouseSC_SW.dat", sep = " ",append = TRUE,
            col.names = F,
            row.names = FALSE)

# Mouse_SC BA


load("~/BNO/Connectivity_data/MouseSC_BA_data.RData")

write(MouseSC_BA_data[["nNodes_MouseSC_BA"]], file = "networkMouseSC_BA.dat", sep = " ")

dadesDeMouseSC_BA <- get.data.frame(ba_MouseSC)
dadesDeMouseSC_BA$from <- as.numeric(as.factor(dadesDeMouseSC_BA$from))
dadesDeMouseSC_BA$to <- as.numeric(as.factor(dadesDeMouseSC_BA$to))
dadesDeMouseSC_BA$name <- NULL

write.table(dadesDeMouseSC_BA, file = "networkMouseSC_BA.dat", sep = " ",append = TRUE,
            col.names = F,
            row.names = FALSE)
