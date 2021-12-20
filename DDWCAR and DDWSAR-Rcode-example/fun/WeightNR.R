WeightNR <- function(clusterID,cluster,RR,region.id,...){
  # this function generate weight matrix (NR) according to clusters and RR
  # clusterID,cluster,RR,region.id seen in function cluster2nb and clusterRR2glist.
  # 2019/3/6
  
  nb0 <- cluster2nb(clusterID = clusterID,cluster = cluster,region.id = region.id)
  glist0 <- clusterRR2glist(clusterID = clusterID,RR = RR,nb = nb0)
  nr <- nb2listw(neighbours = nb0,glist = glist0,...)
  nr
}
