WeightGG <- function(clusterID,cluster,gal,...){
  # this function generates GGW
  # clusterID,cluster,gal seen in function NeighberGG
  # 2019/3/6
  
  nb0 <- NeighberGG(clusterID = clusterID,cluster = cluster,gal = gal)
  ggw <- nb2listw(nb0,glist = NULL,...)
  ggw
}