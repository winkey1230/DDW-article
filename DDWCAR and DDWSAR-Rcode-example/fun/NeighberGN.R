NeighberGN <- function(clusterID,cluster,gal){
  # this function generates neighbor list of GNW according to clusters
  # clusterID and cluster seen in function cluster2nb
  # gal: a list, commonly coming from function read.gal
  # 2019/3/6
  
  region.id <- as.character(attr(gal,"region.id")) 
  ggnb <- NeighberGG(clusterID = clusterID,cluster = cluster,gal = gal)
  nnnb <- cluster2nb(clusterID = clusterID,cluster = cluster,region.id = region.id)
  baseID <- setdiff(region.id,clusterID)
  for (i in 1:length(baseID)) {
    idi <- match(baseID[i],region.id)
    nnnb[[idi]] <- ggnb[[idi]]
  }
  class(nnnb) <- c("nb","nbGN")
  nnnb[which(sapply(nnnb, length)==0)] <- as.integer(0)
  nnnb
}
