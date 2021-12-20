

####################################################
### AUXILIARY FUNCTIONS
####################################################
 
is.neighbor <- function(i, j, adjacencyMatrix, kNearestNeighbors){
  # Test if there is a connected path between i and j within the kNearestNeighbors group according to
  # the adjacencyMatrix
  # Args:
  #   i:  index of the object to test if i and j are connected
  #   j:  index of the object to test if i and j are connected
  #   adjacencyMatrix: a matrix with 0s and 1s indicating the adjacency of the objectcs
  #   kNearestNeighbors: group of objects to test if there is a path between i and j within this group
  # Returns:
  #   TRUE if there is a connected path between i and j within the kNearestNeighbors group and FALSE otherwise

  # # i is removed from kNearestNeighbors to avoid infinite recursion
  # kNearestNeighbors <- kNearestNeighbors[kNearestNeighbors != i]

  if(adjacencyMatrix[i,j]){
    return(TRUE)
  }
  else{
    # i's neighborhood

    neighborhood <- which(as.logical(adjacencyMatrix[i,]))

    for(v in neighborhood){
      if(v %in% kNearestNeighbors){
        # remove v from kNearestNeighbors
        kNearestNeighbors <- kNearestNeighbors[kNearestNeighbors!=v]
        found <- is.neighbor(v,j,adjacencyMatrix, kNearestNeighbors)
        if(found){
          return(TRUE)
        }
      }
    }
  }
  return(FALSE)
}

###########################################################################################################
# 

k.connected.nearest.neighbors <- function(i, distances, k, adjacencyMatrix){
  # Calculate the group of k nearest neighbors of i which are connectected one each other accordintg to adjacentMatrix
  # Args:
  #   i:                index of the object whose k-connected-nearest-neighbors are wanted to obtain
  #   distnaces:        matrix with the distances between each pair of objects
  #   k:                size of the neighborhood
  #   adjacencyMatrix:  a matrix with 0s and 1s indicating the adjacency of the objectcs
  # Returns:
  #   a vector with the index of the k nearest connected neighbors of i
  
  nObjects <- ncol(distances)
  if(k >= nObjects-1){
    stop("Neighborhood size is greater or equal than the number of objects")
  }
  
  # we include the i-th area in order to remove the distance between i and i from the candidate
  # neighbors
  kGroup <- c(i)
  connectedNeighbors.i <- adjacencyMatrix[i,] 
  distToConnectedNeighbors.i <- connectedNeighbors.i * distances[i,]
  for(j in 1:k){
    nearestNeighbor <- which.min(distToConnectedNeighbors.i)
    kGroup <- c(kGroup, nearestNeighbor)
    
    # Include the neighbors of 'nearestNeighbor' as neighbors of 'i' to find the next nearest
    # connected neighbor
    connectedNeighbors.i <- (connectedNeighbors.i | adjacencyMatrix[nearestNeighbor,])
    # set connectedNeighbors.i[kGroup]=NA since the distances have already been used
    connectedNeighbors.i[kGroup] <- NA
    distToConnectedNeighbors.i <- connectedNeighbors.i * distances[i,]
  }
  # remove the area i from the k-nearest neighbors list
  kGroup <- kGroup[kGroup!=i]
  
  return(kGroup)
}#end k.connected.nearest.neighbors

###########################################################################################################

k.nearest.neighbors <- function(i, distances, k){
  # Calculate the group of k nearest neighbors of i (not necesarily are connectected one each other 
  #  accordintg to a adjacentMatrix, geographical border)
  # Args:
  #   distnaces:        matrix with the distances between each pair of objects
  #   k:                size of the neighborhood  #  
  # Returns:
  #   a vector with the index of the k nearest connected neighbors of i
  
  nObjects <- ncol(distances)
  if(k >= nObjects-1){
    stop("Neighborhood size is greater or equal than the number of objects")
  }
  
  # ordered distances between i and the rest of the objects and the corresponding indices 
  # (distance between i and i is included and shouldn't be taken into consideration)
  di <- sort(distances[i,],index.return=TRUE)
  # remove i from the ordere list. 
  i.index <- which(di$ix == i)
  di$x <- di$x[-i.index]
  di$ix <- di$ix[-i.index]
  
  kGroup <- di$ix[1:k]
  return(kGroup)
}#end k.nearest.neighbors

###########################################################################################################

calculate.distanceMatrix <- function(x, distance="euclidean"){
  
  # Calculate the matrix with the distances between each object
  # Args:
  #   x:      numeric vector or a data.frame with the value(s) for each object. Thus, each object 
  #           can be represented by a single point (univariate) or by several points (multivariate: 
  #           each row of the data.frame)
  # distance: type of distance between points, i.e. "euclidean", "canberra", "minkowski", 
  #           "manhattan", etc. (see dist function for more information)
  #      
  # Returns:
  #   A nxn matrix with the distances between each object in x, being n the number of objects.

  if(is.data.frame(x) | (is.numeric(x)&length(x)>1)){
    x <- as.data.frame(x)
    distances <- as.matrix(dist(x,method=distance,diag=TRUE,upper=TRUE))
    return(distances)
  }
  else{
    stop("No allowed data object in calculate.distanceMatrix")
  }
}#end calculate.distanceMatrix



############################################
vector.2.cluster.matrix <- function(partition){
  # convert a vector with a clustering partition into a clustering matrix.
  # Args:
  #   partition: a vector so that partition[i] contains the cluster number for the i-st object
  # Returns:
  #   a matrix where clusterMat[i,j]=1 if the i-st and the j-st object are in the same cluster and 0
  #   otherwise.
  n <- length(partition)
  clusterMat <- matrix(data=0, ncol=n, nrow=n)
  clusterLabels <- unique(partition)
  for(i in clusterLabels){
    indices <- which(partition == i)
    #mark clusters
    clusterMat[indices,indices] <- 1
  }
  return(clusterMat)
}
############################################

cluster.matrix.to.vector <- function(clusterMatrix){
  # convert a cluster matrix with 0s and 1s where x_{ij}=1 if the i-st and j-st element are clustered
  # toghether and 0 otherwise.
  # Args:
  #   clusterMatrix: a nxn matrix where clusterMatrix[i,j]=1 if the i-st and the j-st object are in the 
  #   same cluster and 0 otherwise.
  # Returns:
  #   partition: a vector with the clustering partition denoting the number of cluster for each object.
  clusters <- apply(clusterMatrix,MARGIN = 1, function(x){which(x==1)})
  clusters <- unique(clusters)
  partition <- numeric(nrow(clusterMatrix))
  for(clusterNo in 1:length(clusters)){
    partition[clusters[[clusterNo]]] <- clusterNo
  }
  return(partition)
}

# ############################################

geometric.mean <- function(x){
  # compute the geometric mean of the values in x
  m <-  prod(x)^(1/length(x))
  return(m)
}
####################################################
## END AUXILIARY FUNCTIONS
####################################################






############################################
#' connectedClusters
#' Split clusters from partition if they are not connected according to adjacency matrix W.
#'
#' @param partition: a vector so that partition[i] contains the cluster number for the i-st object
#' @param W: a matrix with 0s and 1s indicating the adjacency of the objectcs (adjacency matrix)
#' 
#' @return a vector with the new clustering partition where the objects within a cluster are 
#' connected according to W (adjacency matrix)
#' 
#' @examples
connectedClusters <- function(partition, W){
  
  clustering.matrix <- vector.2.cluster.matrix(partition)
  M <- clustering.matrix * W
  
  n <- length(partition)
  cluster.index <- 1
  cluster.partition <- numeric(n)
  visited <- logical(n)
  while(any(visited==FALSE)){
    node <- (which(visited==FALSE))[1]
    cluster.partition[node] <- cluster.index
    visited[node] <- TRUE
    neighbors <- which((M[,node]*(!visited))==1)
    while(length(neighbors) > 0){
      neighbor <- neighbors[1]
      if(length(neighbors)>1){
        neighbors <- neighbors[2:length(neighbors)]  
      }
      else{
        neighbors <- integer(0)
      }
      
      visited[neighbor] = TRUE
      cluster.partition[neighbor] <- cluster.index
      neighbors.aux <- which((M[,neighbor]*(!visited))==1)
      neighbors.aux <- neighbors.aux[!(neighbors.aux %in% neighbors)]
      neighbors <- c(neighbors,neighbors.aux)
    }
    cluster.index <- cluster.index + 1
  }
  return(cluster.partition)
}

###########################################################################################################

assign.objects.to.centers <- function(centers,  distanceMatrix, adjacencyMatrix, forceConnectedNeighbors=TRUE,
                                      distanceToBackground=NULL){
  # Creates a clustering partition by assigning each object (Area) to its closest clustering center given that all 
  # the objects in the same cluster have to be connected
  # Args:
  #   centers:          a vector with the indexes of the clustering centers
  #   distnaces:        matrix with the distances between each pair of objects
  #   adjacencyMatrix:  a matrix with 0s and 1s indicating the adjacency of the objectcs
  #   forceConnectedNeighbors: force that areas within a same cluster share a geographical border
  #   distanceToBackground: a vector with the distances between each object and the background data, if background
  #                         cluster is considered, NULL otherwise.
  # Returns:
  #   a vector indicating the cluster for each object
  
  n <- nrow(distanceMatrix)
  nClusters <- length(centers)
  clusteringPartition <- numeric(n)
  
  if(nClusters == 0){
    clusteringPartition <- 1
  }
  else if(nClusters == 1){
    clusteringPartition[1:n] <- 1 
    ## CUIDADO ESTO HABRÁ QUE CAMBIARLO SI CONSIDERAMOS EL BACKGROUND CLUSTER
  }
  else{
    # initialize the clustering partition with the clustering centers
    clusteringPartition[centers] <- 1:nClusters
    
    if(forceConnectedNeighbors){ #forceConnectedNeighbors=TRUE
      # connected contains a 1 if the jth object shares a border with any object in the ith cluster
      # (i is the row index and j is the column index). At the first iteration, clusters only contain the centers, therefore,
      # only the neighbors of cluster centers contain a 1 in 'connected'. If there is no common border, the cell
      # contains an Inf value
      connected <- adjacencyMatrix[centers,]
      connected <- as.matrix(connected)
      connected[connected==0] <- Inf
      
      clusterDistances <- distanceMatrix[centers,]
        # initialize the distance between a cluster center and itself to an Inf value
      clusterDistances[,centers] <- Inf
      
      used <- logical(n)
      used[centers] <- TRUE
      
      it <- 1
      while(any(!used) & it <= n){
        
        distanceToNeighbors <-  clusterDistances *  connected
        
        index <- which.min(distanceToNeighbors)
        # R aggregates the distanceToNeighbors matrix by columns to obtain the index. Thus indices 1:5
        # corresponds to the distances from the first element to cluster centers 1:5.
        clusterIndex <- index %% nClusters ##resto
        objectIndex <- index %/% nClusters + 1
        # when the reminder is 0, that is because the object is in the last cluster
        if(clusterIndex == 0){
          clusterIndex <- nClusters
          objectIndex <- objectIndex - 1
        }
        
        if(!is.null(distanceToBackground)){
          #background cluster is considered. If the selected object is closer to the background than to the centroid, 
          # the object is assigned to the background cluster and the cluster assignation process finishes since the
          # rest of the objects will be closer to the background than to any other centroid.
          if(distanceToBackground[objectIndex] < distanceMatrix[centers[clusterIndex],objectIndex]){
            # 0 is used as index for backgroundCluster
            clusterIndex <- 0
            #set every object as used since the rest of the objects will be also closer to the backgroundCluster than
            #to any other centroid
            # used[] <- TRUE
          }
        }
        
        # add the object to the corresponding cluster
        clusteringPartition[objectIndex] <- clusterIndex
        # add the neighbors of the object to its corresponding cluster in 'connected'
        connected[clusterIndex,adjacencyMatrix[objectIndex,]==1] <- 1
        # set the object as used
        used[objectIndex] <- TRUE
        # the object has already assigned to a cluster, therefore it cannot be chosen anymore
        # -> we set the distance between the cluster centers and the object to Inf
        clusterDistances[,objectIndex] <- Inf
        
        # print(paste0("it ",it, " : obj",objectIndex, " -> c", clusterIndex, " (Dist.cent: ", 
        #              distanceMatrix[centers[clusterIndex],objectIndex], " -- Dist.back: ",
        #              distanceToBackground[objectIndex], ")"))
        # print(paste0("     UsedObjects: ",sum(used)))
        it <- it+1
      }
      
      if(!is.null(distanceToBackground)){
        # as 0 is used as id for background cluster, we add 1 to clusteringPartition. Thus the first cluster (#1) is
        # the background cluster and then the rest of the clusters acording to the centroids detected by the algorithm. 
        clusteringPartition <- clusteringPartition + 1
      }
      
      
      
    }
    else{
      # forceConnectedNeighbors=FALSE
      if(!is.null(distanceToBackground)){
        # we include the background cluster as cluster #1 for the assignation process
        distanceToCenter <- cbind(distanceToBackground,distanceMatrix[ ,centers])
      }
      else{
        distanceToCenter <- distanceMatrix[ ,centers]  
      }
      clusteringPartition <- apply(distanceToCenter, MARGIN=1, FUN=which.min)
    }
  }
  
  return(clusteringPartition)
}

###########################################################################################################

center.detection.boxplot.outliers <- function(rho, delta, coef=1.5){
  # outlier detection to obtain clustering centroids
  # Args:
  #   rho:   numeric vector calculated using calculate.rho function
  #   delta: numeric vector calculated using calculate.delta function
  #   coef:  coef to determine the 'whiskers' of the boxplot (+- coef*(Q3-Q1)).
  #          All values away form this limits are considered outliers (cluster centers)
  # Returns:
  #   a set with the objects detected as clustering centers
  
  gama <- rho*delta
  outlierValues <- boxplot.stats(gama, coef =coef)$out
  centerIndexes <- which(gama %in% outlierValues)
  return(centerIndexes)
}

###########################################################################################################

calculate.delta <- function(rho, distances){
  # Calculate the delta values used in the clustering algorithm
  # Args:
  #   rho:        rho values for each object
  #   distnaces:  matrix with the distances between each pair of objects
  # Returns:
  #   a vector with the delta values for each object
  
  nObjects <- length(rho)
  delta <- numeric(nObjects)
  rhoOrdered <- sort(rho,index.return=TRUE)
  for(index in 1:(nObjects-1)){
    i <- rhoOrdered$ix[index]
    delta[i] <- min(distances[i,rhoOrdered$ix[(index+1):nObjects]])
  }#end for
  delta[rhoOrdered$ix[nObjects]] <- max(delta)
  return(delta)
}#end calculate.delta

###########################################################################################################

calculateAdjacencyMatrix.neighborhoodLevel <- function(adjacencyMatrix, neighborhoodLevel=2){
  # Given an adjacency matrix, computes a new adjacency matrix taken into account a diferent neighborhood level
  # Args:
  #   adjacencyMatrix:    a matrix with 0s and 1s indicating the adjacency of the objectcs
  #   neighborhoodLevel:  The neighborhood level to be considered. If =1, the adjacencyMatrix is unaltered,
  #                       if =2 2nd order neighborhood is considered to recalculate the adjacencymatrix, ...      
  #      
  # Returns:
  #   A nxn adjacency matrix taking into account the neighborhood level
  if(neighborhoodLevel>1){
    require(Matrix)
    aux <- as(adjacencyMatrix,"dgCMatrix")
    
    for(i in 2:neighborhoodLevel){
      aux <- aux%*%aux
      diag(aux) <- 0
      aux[aux>0] <- 1
    }
    
    W.new <- as.matrix(aux)
    return(W.new)
  } else{
    return(adjacencyMatrix)
  }
}

###########################################################################################################

calculate.rho.geographical.neighbors <- function(distances, adjacencyMatrix, neighborhoodLevel=1){
  # Calculate the rho values (inverse of average distance between each area and the areas in its neighborhood). The 
  # neighborhoodLevel paramater can be used to expand the neighborhood considering the neighbors of
  # the area's neighbors. 
  # Args:
  #   neighborhoodLevel:                neighborhood considered to calculate the neighbors of each area Ai. If neighborhoodLevel=1,
  #                                     the neighbors of each area Ai are those sharing a border with Ai (geographical neighbors). 
  #                                     If neighborhoodLevel=2, the neighbors of Ai include its geographical neighbors and the 
  #                                     neighbors of these neighbors (2nd level neighborhood), and so on for other levels.
  #   distnaces:                        matrix with the distances between each pair of objects
  #   adjacencyMatrix:                  a matrix with 0s and 1s indicating the adjacency of the objectcs
  # Returns:
  #   a vector with the rho values for each object
  
  nObjects <- ncol(distances)
  rho <- numeric(nObjects)
  
  adjacencyMatrix <- calculateAdjacencyMatrix.neighborhoodLevel(adjacencyMatrix, 
                                                                neighborhoodLevel=neighborhoodLevel)
  
  # number of geographical neihgbors for each area
  k <- colSums(adjacencyMatrix)
  # we add 1e-10 to avoid 0 values
  sumDistancesToNeighbors <- colSums(adjacencyMatrix * distances) + 1e-10
  rho <- k / sumDistancesToNeighbors
  return(rho)
}#end calculate.rho

############################################
#' getClusterPartition
#' Obtain a clustering partition of the data using the DBSC clustering algorithm. See Santafe et all (2021) "Dealing with risk
#' discontinuities to estimate cancer mortality risks when the number of small areas is large" 
#'
#' @param x data objects to by clusters
#' @param W adjacency matrix describing geographical adjacency between areas
#' @param l states for the neighborhoodLevel and it is used to obtain the neighborhood of each area Ai, and therefore to 
#' calculate the local density of that area (the local density is a measure used in the clustering algorithm to identify 
#' the cluster centroids). If l=1, only areas that share a border with Ai are considered as neighbors of Ai; if l=2, 
#' areas within a distance of 1 area from Ai are also considered as its neighbors, and so on. 
#' @param adjacencyMatrix.neighborhoodLevel an adjacency matrix as W but considering a specific neighborhood level. 
#' The calculation of the neighborhood of each area Ai for a given neighborhoodLevel (l) may be very time consuming 
#' when the number of areas is large and l is also big. However, these neighborhoods only depend on W and l, therefore 
#' they can be calculated in advance to be used as a parameter in the clustering algorithm. If adjacencyMatrix.neighborhoodLevel 
#' parameter is used, parameter l is ignored.
#' @param distance: type of distance between points, i.e. "euclidean", "canberra", "minkowski", "manhattan", etc. (see 
#' dist function for more information)
#' @param backgroundData: A vector with the same length as the objects in x and also measured on the same scale. This object 
#' represents the centroid of a background cluster. In the clustering process, each object in 'x' can be assigned to the 
#' background cluster if it is closer to the background cluster than to any other centroid. If backgroundData is NULL, the 
#' background cluster is not considered and objects in 'x' are always assigned to the closest centroid. 
#' @param splitBackgroundCluster: when using a background cluster, it results in a unique background that may be geographically 
#' unconnected. If splitBackgroundCluster=FALSE, this unique background cluster is preserved as-is. By contrast, if 
#' splitBackgroundCluster=TRUE, the background cluster is split so that the subsets which are not geographically connected form 
#' new clusters.
#'
#' @return a vector with the clustering partition
#'
#' @examples

getClusterPartition <- function(x, W, l=1, adjacencyMatrix.neighborhoodLevel=NULL, distance="euclidean", backgroundData=NULL,
                                splitBackgroundCluster=TRUE){
  
  
  # no. of elements to be clustered
  n <- nrow(x)
  if(!is.null(backgroundData)){
    # if background cluster is considered, we add the backgroundData to x in order to calculate the distance of each
    # element in x to the background cluster
    auxData <- rbind(x,backgroundData)
    auxDistanceMatrix <- calculate.distanceMatrix(auxData, distance=distance)
    distanceMatrix <- auxDistanceMatrix[1:n,1:n]
    distanceToBackground <- auxDistanceMatrix[n+1,1:n]
  }
  else{
    distanceMatrix <- calculate.distanceMatrix(x, distance=distance)
    distanceToBackground <- NULL
  }
  
  
  if(is.null(adjacencyMatrix.neighborhoodLevel)){
    adjacencyMatrix.neighborhoodLevel <- W
  }
  else{
    #neighborhood are given in adjacencyMatrix.neighborhoodLevel so we use it as W and l=1
    l <- 1
  }
  
  
  rho <- calculate.rho.geographical.neighbors(distances=distanceMatrix, adjacencyMatrix=adjacencyMatrix.neighborhoodLevel, 
                                              neighborhoodLevel=l)
  delta <- calculate.delta(rho,distanceMatrix)
  
  
  centers <- center.detection.boxplot.outliers(rho,delta,coef=2)  
  
  partition <- assign.objects.to.centers(centers, distanceMatrix, adjacencyMatrix=W,
                                         distanceToBackground=distanceToBackground)  
  
  if(splitBackgroundCluster & !is.null(backgroundData)){
    partition <- connectedClusters(partition,W)
  }
  return(partition)
}

  