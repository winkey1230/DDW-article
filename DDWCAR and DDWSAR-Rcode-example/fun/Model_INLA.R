temp.folder <- function(mainDir=NULL){
  ## Create a temporal folder in the main directory. The name of the temporal folder is based on System time
  ## (including milliseconds) to avoid that different models running in parallel use the same temporal folder
  tmp.folder <- format(Sys.time(),"%Y-%m-%d_%H_%M_%OS5")
  
  if(is.null(mainDir)) mainDir <- getwd()
  
  if(!file.exists(tmp.folder)) {
    dir.create(file.path(mainDir, tmp.folder))
  }
  else {
    warning(paste0("Temp folder ", tmp.folder, " already exist. Files inside the folder are removed and it
                   may cause problems if different models are running in pararell"))
    do.call(file.remove, list(list.files(file.path(mainDir, tmp.folder), full.names=TRUE)))
  }
  
  return(file.path(mainDir, tmp.folder))
}


model.selection <- function(cluster.partition, Y.real, E.real, C, Carto,
                            max.cluster=NULL, final.cluster=NULL, start=NULL,
                            plot.dic=TRUE, strategy="simplified.laplace", 
                            tempfolder.path=NULL, tempfolder.remove=TRUE){
                            
     #### Fit a model to each cluster partition candidate (up to "max.cluster")
     #### and choose the best by minimising the Deviance Information Criterion.

tempfolder.dir <- temp.folder(tempfolder.path)


time <- system.time({

    ## Define our uniform prior distributions ##
    sdunif="expression:
      logdens=-log_precision/2;
      return(logdens)"
  
    lunif = "expression:
      a = 1;
      b = 1;
      beta = exp(theta)/(1+exp(theta));
      logdens = lgamma(a+b)-lgamma(a)-lgamma(b)+(a-1)*log(beta)+(b-1)*log(1-beta);
      log_jacobian = log(beta*(1-beta));
      return(logdens+log_jacobian)"
    
    R.Leroux <- diag(dim(C)[1])-C
  
  
    ## N? of areas and partition candidates ##
    n <- length(Y.real)
    n.partition <- nrow(cluster.partition)
    
    if(n.partition==1) plot.dic=FALSE
    
    if(is.null(final.cluster)){
      
      ## Store the DIC values
      if(is.null(max.cluster)){
        max.cluster <- n.partition
      }else{
        if(max.cluster>n.partition) stop("\n WARNING: 'max.cluster' is greater than the number of partition candidates.")  
      }
      
      if(is.null(start)){
        dic.list <- rep(NA, max.cluster)
        waic.list <- rep(NA, max.cluster) 
      }else{
        dic.list <- rep(NA, max.cluster-start+1)
        waic.list <- rep(NA, max.cluster-start+1) 
      }
      
      
      ## Fit a model to each cluster configuration candidate ##
      if(is.null(start)){
        from <- 1
      }else{
        from <- start
      }
      
      
      ## graph device to plot DIC if the plot is required
      if(plot.dic==TRUE){
        if(.Platform$OS.type == "unix"){
          X11()
        }
        else{
          win.graph()
        }
      }

      aux <- 1
      for(i in from:max.cluster){
        
        factor.clust <- cluster.partition[i,]
        
        if(length(unique(factor.clust))==1){
          
          ## Model without cluster structure ##
          data.temp <- data.frame(Y.real=Y.real, E.real=E.real, region=1:n)
        
          R.Leroux <- diag(dim(C)[1])-C
        
          formula <- Y.real ~ f(region, model="generic1", Cmatrix = R.Leroux, constr=TRUE,
                                hyper=list(prec=list(prior=sdunif),beta=list(prior=lunif, initial=0)))
        
          model <-tryCatch(
            inla(formula, family="poisson", data=data.temp, E=E.real,
                        control.fixed=list(mean=0, mean.intercept=0, prec=0.1, prec.intercept=0.001),
                        control.compute=list(dic=TRUE, mlik=TRUE, cpo=TRUE, waic=TRUE, config=TRUE),
                        control.predictor=list(compute=TRUE, cdf=c(log(1))),
                        control.inla=list(strategy=strategy, npoints=21))
            , error=function(e) NULL)
          
          if(is.null(model)){
            dic.list[aux] <- Inf
            waic.list[aux] <- Inf
          }
          else{
            dic.list[aux] <- model$dic$dic
            waic.list[aux] <- model$waic$waic
          }
      
        }else{
      
          ## Model with cluster structure ##
          data.temp <- data.frame(Y.real=Y.real, E.real=E.real, region=1:n,
                                  clust=as.numeric(as.factor(factor.clust)))
      
          ## Cluster configuration map ##
          lista <- factor.clust
          cluster.map <- unionSpatialPolygons(Carto,lista)
          nb2INLA(poly2nb(cluster.map), file=file.path(tempfolder.dir, "cluster_nb.inla"))
          
          g <- inla.read.graph(file.path(tempfolder.dir, "cluster_nb.inla"))
          R.clust = matrix(0, g$n, g$n)
          for (ii in 1:g$n){
            R.clust[ii,ii]=g$nnbs[[ii]]
            R.clust[ii,g$nbs[[ii]]]=-1
          }
          R.clust.Leroux <- diag(dim(R.clust)[1])-R.clust

	        ## INLA formula ##
          all.lc <- inla.make.lincomb(Predictor = rep(1/(n),n))
          names(all.lc) <- "intercept"

	        if(length(unique(data.temp$clust))>5){
	          
	          ## Model with two-level random effects ##
	          formula <- Y.real ~ f(region, model="generic1", Cmatrix = R.Leroux, constr=TRUE,
	                                hyper=list(prec=list(prior=sdunif),beta=list(prior=lunif, initial=0))) +
	                              f(clust, model="generic1", Cmatrix=R.clust.Leroux, constr=TRUE,
	                                hyper=list(prec=list(prior=sdunif),beta=list(prior=lunif, initial=0)))
	          
	        }else{

	          ## Model with random effects (region) + fixed effects (cluster) ##
	          beta <- names(table(data.temp$clust))[-1]
	          for(i in beta){
	            beta.lc <- list(list(list(list(weight=1)),list(list(weight=1))))
	            names(beta.lc) <- paste("beta",i,sep="")
	            names(beta.lc[[1]][[1]]) <- "(Intercept)"
	            names(beta.lc[[1]][[2]]) <- paste("factor(clust)",i,sep="")
	            
	            all.lc <- c(all.lc,beta.lc)
	          }
	          
	          formula <- Y.real ~ f(region, model="generic1", Cmatrix = R.Leroux, constr=TRUE,
	                                hyper=list(prec=list(prior=sdunif),beta=list(prior=lunif, initial=0))) +
	            factor(clust)
	        }
          
          model <-tryCatch(
            inla(formula, family="poisson", data=data.temp, E=E.real,
                 control.fixed=list(mean=0, mean.intercept=0, prec=0.1, prec.intercept=0.001),
                 control.compute=list(dic=TRUE, mlik=TRUE, cpo=TRUE, waic=TRUE, config=TRUE),
                 control.predictor=list(compute=TRUE, cdf=c(log(1))),
                 lincomb=all.lc,
                 control.inla=list(strategy=strategy, npoints=21))
            , error=function(e) NULL)
          
          if(is.null(model)){
            dic.list[aux] <- Inf
            waic.list[aux] <- Inf
          }
          else{
            dic.list[aux] <- model$dic$dic
            waic.list[aux] <- model$waic$waic
          }
          
        }
        
        ## Plot DIC values ##
	      if(plot.dic==TRUE & is.finite(dic.list[aux])) {
	      
	        if(aux==1){
	          y.lim <- dic.list[aux]*c(1/1.01,1.005)
	          plot(seq(1,length(dic.list)), dic.list, type="l", xlab="Number of partition candidates", ylab="DIC",
	               main=paste("Model",aux,"of",max.cluster), ylim=y.lim)
	          points(1,dic.list[aux],pch=19)
	          mtext(paste("Minimun DIC value:",round(min(dic.list,na.rm=T),2)), side=3, line=-2)
	        }
	        
	        if(aux>1){
	          if(dic.list[aux]<y.lim[1]){
	            y.lim[1] <- y.lim[1]-abs(y.lim[1]-dic.list[aux])
	          }
	          if(dic.list[aux]>y.lim[2]){
	            y.lim[2] <- y.lim[2]+abs(y.lim[2]-dic.list[aux])
	          }
	          plot(seq(1,length(dic.list)), dic.list, type="l", xlab="Number of partition candidates", ylab="DIC",
	               main=paste("Model",aux,"of",max.cluster-from+1), ylim=y.lim)
	          points(aux,dic.list[aux],pch=19)
	          mtext(paste("Minimun DIC value:",round(min(dic.list,na.rm=T),2)), side=3, line=-2)
	        }
	      }
        
        aux <- aux+1
      }
    }else{
      if(final.cluster>n.partition) stop("\n WARNING: 'final.cluster' is greater than the number of partition candidates.")  
    }
    
    
    ## FINAL MODEL ##
    if(n.partition==1){
      results <- list(model.final=model, factor.clust=factor.clust, dic.list=NULL, waic.list=NULL)

    }else{
      
      cat("\n\n ************* RUNNING FINAL MODEL ************* \n\n")
      if(is.null(final.cluster)){
        best.model <- from+which.min(dic.list)-1
      }else{
        best.model <- final.cluster
        plot.dic <- FALSE
        dic.list <- NULL
        waic.list <- NULL
      }
      
      factor.clust <- cluster.partition[best.model,]
      
      if(length(unique(factor.clust))==1){
        
        ## Model without cluster structure ##
        data.temp <- data.frame(Y.real=Y.real, E.real=E.real, region=1:n)
        
        R.Leroux <- diag(dim(C)[1])-C
        
        formula <- Y.real ~ f(region, model="generic1", Cmatrix = R.Leroux, constr=TRUE,
                              hyper=list(prec=list(prior=sdunif),beta=list(prior=lunif, initial=0)))
        
        model <- inla(formula, family="poisson", data=data.temp, E=E.real,
                      control.fixed=list(mean=0, mean.intercept=0, prec=0.1, prec.intercept=0.001),
                      control.compute=list(dic=TRUE, mlik=TRUE, cpo=TRUE, waic=TRUE, config=TRUE),
                      control.predictor=list(compute=TRUE, cdf=c(log(1))),
                      control.inla=list(strategy=strategy, npoints=21))
        
      }else{
        
        ## Model with cluster structure ##
        data.temp <- data.frame(Y.real=Y.real, E.real=E.real, region=1:n,
                                clust=as.numeric(as.factor(factor.clust)))
        
        ## Cluster configuration map ##
        lista <- factor.clust
        cluster.map <- unionSpatialPolygons(Carto,lista)
        nb2INLA(poly2nb(cluster.map), file=file.path(tempfolder.dir, "cluster_nb.inla"))
        
        g <- inla.read.graph(file.path(tempfolder.dir, "cluster_nb.inla"))
        R.clust = matrix(0, g$n, g$n)
        for (ii in 1:g$n){
          R.clust[ii,ii]=g$nnbs[[ii]]
          R.clust[ii,g$nbs[[ii]]]=-1
        }
        R.clust.Leroux <- diag(dim(R.clust)[1])-R.clust
        
        ## INLA formula ##
        all.lc <- inla.make.lincomb(Predictor = rep(1/(n),n))
        names(all.lc) <- "intercept"
        
        if(length(unique(data.temp$clust))>5){
          
          ## Model with two-level random effects ##
          formula <- Y.real ~ f(region, model="generic1", Cmatrix = R.Leroux, constr=TRUE,
                                hyper=list(prec=list(prior=sdunif),beta=list(prior=lunif, initial=0))) +
            f(clust, model="generic1", Cmatrix=R.clust.Leroux, constr=TRUE,
              hyper=list(prec=list(prior=sdunif),beta=list(prior=lunif, initial=0)))
          
        }else{
          
          ## Model with random effects (region) + fixed effects (cluster) ##
          beta <- names(table(data.temp$clust))[-1]
          for(i in beta){
            beta.lc <- list(list(list(list(weight=1)),list(list(weight=1))))
            names(beta.lc) <- paste("beta",i,sep="")
            names(beta.lc[[1]][[1]]) <- "(Intercept)"
            names(beta.lc[[1]][[2]]) <- paste("factor(clust)",i,sep="")
            
            all.lc <- c(all.lc,beta.lc)
          }
          
          formula <- Y.real ~ f(region, model="generic1", Cmatrix = R.Leroux, constr=TRUE,
                                hyper=list(prec=list(prior=sdunif),beta=list(prior=lunif, initial=0))) +
            factor(clust)
        }
        
        model <- inla(formula, family="poisson", data=data.temp, E=E.real,
                      control.fixed=list(mean=0, mean.intercept=0, prec=0.1, prec.intercept=0.001),
                      control.compute=list(dic=TRUE, mlik=TRUE, cpo=TRUE, waic=TRUE, config=TRUE),
                      control.predictor=list(compute=TRUE, cdf=c(log(1))),
                      lincomb=all.lc,
                      control.inla=list(strategy=strategy, npoints=21))
      }
    
      if(plot.dic==TRUE){
        points(which.min(dic.list),min(dic.list,na.rm=TRUE), pch=19, col="red")
        lines(rep(which.min(dic.list),2), c(0,min(dic.list,na.rm=TRUE)), lty=2, col="red")
        axis(1, at=which.min(dic.list), labels=which.min(dic.list), col.axis="red")
        mtext(paste("Minimun DIC value:",round(min(dic.list,na.rm=T),2)), side=3, line=-2, col="red")
      }

      results <- list(model.final=model, factor.clust=factor.clust, dic.list=dic.list, waic.list=waic.list)
    }
  })

  ## Remove temp folder if requested
  if(tempfolder.remove){
    unlink(tempfolder.dir, recursive=TRUE)
  }


  results <- append(results, list(cpu.time=time[3]))
  return(results)
}
