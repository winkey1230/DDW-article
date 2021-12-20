cluster2nb <- function(clusterID,cluster,region.id){
  # this function transform the spatial location to neighbor according to clusters,
  # any two locations from the same cluster or base are neighbors, otherwise they are not.
  # clusterID: a vector of character or integer,representing the locations belonging to clusters,
  # cluster: a vector of character or integer, representing which cluster the corresponding belongs to.
  # region.id: all of location ids in the entire area. The length are larger than that of clusterID
  # 2019/3/6
  
  # checking
  clusterID <- as.character(clusterID)
  cluster <- as.character(cluster)
  region.id <- as.character(region.id)
  if(length(setdiff(clusterID,region.id))!=0) stop("each element of clusterID must be in region.id")
  if(length(clusterID) != length(cluster)) stop("the length of clusterID must equal to that of cluster")
  if(sum(duplicated(clusterID),duplicated(region.id))!=0) stop("there are duplicated elements in clusterID or region.id")
  # merge cluster and base
  clusters <- data.frame(code = clusterID,num = cluster,stringsAsFactors = F)
  All <- merge.data.frame(clusters,data.frame(code = region.id,stringsAsFactors = F),by = "code",all.y = T)
  if(sum(is.na(All$num))!=0){
    if("base" %in% All$num) stop("the elements of cluster can't be 'base'")
    All$num[is.na(All$num)] <- "base"
  }
  # finding neignbor
  res <- list()
  for (i in 1:length(region.id)) {
    aa <- All[All$code==region.id[i],2]
    od <- subset(All,num==aa,select=1)[[1]]
    od <- setdiff(od,region.id[i])
    neighber <- match(od,region.id)
    res[[i]] <- sort(neighber)
  }
  class(res) <- c("nb","nbNN")
  res[which(sapply(res, length)==0)] <- as.integer(0)
  attr(res, "region.id") <- region.id
  attr(res, "gis") <- TRUE
  attr(res, "call") <- TRUE
  res <- sym.attr.nb(res)
  res
}
