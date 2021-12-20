DDW <- function(type,clusterID=NULL,cluster=NULL,RR=NULL,gal=NULL,region.id=NULL,...){
  # this function generates the classic Weight matrix (AG) and DDW: GG,GN,GR,NG,NR,NN
  # clusterID,cluster,RR,gal,region.id seen in the six funcitons: WeightGG,WeightGN,...,WeightNN.
  # type: a character, which must be one of the seven characters:AG,GG,GN,GR,NG,NR,NN, representing
  #       which type of DDW you will construct.
  # 2019/3/7
  
  # checking
  if (!type %in% c("AG","GG","GN","GR","NG","NR","NN")) 
    stop("'type' must be one of these characters:AG,GG,GN,GR,NG,NR,NN")
  if (type == "AG"){
    if(is.null(gal)) stop("The parameter 'gal' is necessary")
        else swm <- nb2listw(gal,...)
  }
  if (type == "GG"){
    if(is.null(clusterID)|is.null(cluster)|is.null(gal))
      stop("'clusterID','cluster','gal' are necessary")
    else swm <- WeightGG(clusterID = clusterID,cluster = cluster,gal = gal,...)
  }
  if (type == "GN"){
    if(is.null(clusterID)|is.null(cluster)|is.null(gal))
      stop("'clusterID','cluster','gal' are necessary")
    else swm <- WeightGN(clusterID = clusterID,cluster = cluster,gal = gal,...)
  }
  if (type == "GR"){
    if(is.null(clusterID)|is.null(cluster)|is.null(gal)|is.null(RR))
      stop("'clusterID','cluster','gal','RR' are necessary")
    else swm <- WeightGR(clusterID = clusterID,cluster = cluster,gal = gal,RR = RR,...)
  }
  if (type == "NG"){
    if(is.null(clusterID)|is.null(cluster)|is.null(gal))
      stop("'clusterID','cluster','gal' are necessary")
    else swm <- WeightNG(clusterID = clusterID,cluster = cluster,gal = gal,...)
  }
  if (type == "NN"){
    if(is.null(clusterID)|is.null(cluster)|is.null(region.id))
      stop("'clusterID','cluster','region.id' are necessary")
    else swm <- WeightNN(clusterID = clusterID,cluster = cluster,region.id = region.id,...)
  }
  if (type == "NR"){
    if(is.null(clusterID)|is.null(cluster)|is.null(region.id)|is.null(RR))
      stop("'clusterID','cluster','RR','region.id' are necessary")
    else swm <- WeightNR(clusterID = clusterID,cluster = cluster,region.id = region.id,RR = RR,...)
  }
  return(swm)
}