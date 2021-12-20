WeightGN <- function(clusterID,cluster,gal,...){
  # this function generates GNW
  # clusterID,cluster,gal seen in function NeighberGN
  # 2019/3/7
  
  gnnb <- NeighberGN(clusterID = clusterID,cluster = cluster,gal = gal)
  gnw <- nb2listw(gnnb,glist = NULL,...)
  gnw
}