WeightNG <- function(clusterID,cluster,gal,...){
  # this function generates GNW
  # clusterID,cluster,gal seen in function NeighberGN
  # 2019/3/7
  
  ngnb <- NeighberNG(clusterID = clusterID,cluster = cluster,gal = gal)
  ngw <- nb2listw(ngnb,glist = NULL,...)
  ngw
}