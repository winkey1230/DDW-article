WeightGR <- function(clusterID,cluster,RR,gal,...){
  # this function generate weight matrix (GR) according to clusters and RR and gal
  # clusterID,cluster,RR seen in function cluster2nb and clusterRR2glist.
  # gal seen in function NeighberGN
  # 2019/3/6
  
  nb0 <- NeighberGN(clusterID = clusterID,cluster = cluster,gal = gal)
  glist0 <- clusterRR2glist(clusterID = clusterID,RR = RR,nb = nb0)
  nr <- nb2listw(neighbours = nb0,glist = glist0,...)
  nr
}
