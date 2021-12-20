WeightNN <- function(clusterID,cluster,region.id,...){
  # this function generates NNW according to clusters
  # clusterID,cluster,region.id seen in function cluster2nb
  # 2019/3/6
  nnnb <- cluster2nb(clusterID = clusterID,cluster = cluster,region.id = region.id)
  nnw <- nb2listw(nnnb,...)
  nnw
}