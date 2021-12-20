NeighberGG <- function(clusterID,cluster,gal){
  # this function generates neighbor list of GGW according to clusters
  # clusterID and cluster seen in function cluster2nb
  # gal: a list, commonly coming from function read.gal
  # 2019/3/6
  
  # checking
  clusterID <- as.character(clusterID)
  cluster <- as.character(cluster)
  region.id <- attr(gal,"region.id")
  region.id <- as.character(region.id)
  if(length(setdiff(clusterID,region.id))!=0) stop("each element of clusterID must be in region.id")
  if(length(clusterID) != length(cluster)) stop("the length of clusterID must equal to that of cluster")
  if(sum(duplicated(clusterID))) stop("there are duplicated elements in clusterID")
  # merge cluster and base
  clusters <- data.frame(code = clusterID,num = cluster,stringsAsFactors = F)
  All <- merge.data.frame(clusters,data.frame(code = region.id,stringsAsFactors = F),by = "code",all.y = T)
  if(sum(is.na(All$num))!=0){
    if("base" %in% All$num) stop("the elements of cluster can't be 'base'")
    All$num[is.na(All$num)] <- "base"
  }
  # building the neighber of GGW
  for (i in 1:length(region.id)) {
    ocode <- region.id[i]
    nbcode <- region.id[gal[[i]]]
    nbcluster <- All[match(nbcode,All$code),2]
    ocluster <- All[match(ocode,All$code),2]
    newnb <- gal[[i]][nbcluster==ocluster]
    gal[[i]] <- newnb
  }
  class(gal) <- c("nb","nbGG")
  gal[which(sapply(gal, length)==0)] <- as.integer(0)
  gal
}




