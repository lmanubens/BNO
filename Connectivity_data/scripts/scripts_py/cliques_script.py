from igraph import *
import networkx as nx

gFly_Optic = Graph.Read_Ncol("gFlyOpticSData.ncol")

print(gFly_Optic.es['weight'])

cliques_gFly = gFly_Optic.cliques()

print(cliques_gFly)
