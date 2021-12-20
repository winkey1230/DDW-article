clusterRR2glist <- function(clusterID,RR,nb){
  # this function generates glist according to the RRs of location in cluster
  # the locations out of clusters are assigned to equal cluster.
  # this glist is commonly used in function WeightNR
  # clusterID,cluster,region.id seen in function cluster2nb
  # RR: a vector of numeric, often coming from the 11th column of gis.txt file from satscan
  # nb is a neighbor commonly from function cluster2nb or NeighberGN
  
  # checking
  clusterID <- as.character(clusterID)
  if(length(clusterID) != length(RR)) 
    stop("the length of clusterID, RR and cluster  must be equal to each other")
  if(sum(duplicated(clusterID))!=0) stop("there are duplicated elements in clusterID")
  if(sum(class(nb) %in% c("nbNN","nbGN"))==0) warning("nb is not the class: 'nbNN' or 'nbGN'")
  region.id <- attr(nb,"region.id")
  region.id <- as.character(region.id)
  if(length(setdiff(clusterID,region.id))!=0) stop("each element of clusterID must be in nb's region.id")
  # construct glist
  glist <- list()
  for (i in 1:length(region.id)) {
    nbcode <- region.id[nb[[i]]]
    ocode <- region.id[i]
    if(!region.id[i] %in% clusterID) {
       wei <- rep(1,length(nbcode))
       glist[[i]] <- wei
    } else{
      oRR <- RR[ocode == clusterID]
      nbRR <- RR[match(nbcode,clusterID)]
      mm <- c((nbRR+1)/(oRR+1),(oRR+1)/(nbRR+1))
      mm <- matrix(mm,ncol = 2)
      wei <- 1/apply(mm, 1, max)
      glist[[i]] <- wei
    }
  }
  glist[which(sapply(glist, length)==0)] <- list(NULL)
  glist
}



