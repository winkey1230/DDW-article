samplePosterirorRisks <- function(Model, n.MCMC=1000){
  ## sample the posterior distribution of the risks to obtain `n.MCMC` sampled risks for each area. These values are used to calculate the
  ## DIC score in a set of observed data different to the one used to fit the model.
  # Args:
  #   Model:  inla model already fitted. Note that `config=TRUE` must be set in the `control.compute` parameter when fitting the model.
  #   n.MCMC: number of risks to be sampled for each area from each corresponding posterior distribution
  # Returns:
  #   the same inla model given as paremeter (`Model`) where a new field is included (Model$risk.MCMC). `risk.MCMC` is a n x n.MCMC matrix 
  #   with the sampled risks. Where n is the number of areas and n.MCMC is the number of sampled risks for each area.
  n <- length(Model$marginals.fitted.values)
  model.MCMC <- inla.posterior.sample(n.MCMC, Model)
  risk.MCMC <- array(unlist(lapply(model.MCMC, function(x) exp(x$latent[1:n]))), dim=c(n,n.MCMC))
  ## save the sampled risks in the model
  Model$risk.MCMC <- risk.MCMC
  return(Model)
}

#####################################################################################################################################

approximatedDIC <- function(Model, Obs, Exp){
  ## Compute the DIC of a model given the Observed and Expected data in a time period. The calculation of the DIC score is an approximation
  ## based on sampled values for the posterior risks
  # Args:
  #   Model:  inla model already fitted which includes the field Model$risk.MCMC with a n x n.MCMC matrix with the sampled risks (n is the 
  #     number of areas and n.MCMC is the number of sampled risks for each area). This model can be obtained with the samplePosterirorRisks
  #     function.
  #   Obs: Observed values for each area in a specific time period
  #   Exp: Expected observations for each area in a specific time period
  # Returns:
  #   the value of the DIC score for the model and the given Observed and Expected data
  mu.MCMC <- apply(Model$risk.MCMC, MARGIN=2, FUN=function(x) Exp*x)
  mean.deviance <- mean(apply(mu.MCMC, 2, function(x) -2*sum(log(dpois(Obs,x)))))
  deviance.mean <- -2*sum(log(dpois(Obs,apply(mu.MCMC,1,mean))))
  pD <- mean.deviance-deviance.mean
  DIC <- mean.deviance+pD
  return(DIC)
}


# ALGO FALLA PORQUE NO HACE EL MISMO CÁLCULO QUE POR AÑOS
# # ADEMAS aunque no se usan bucles tarda más tiempo en ejecutarse
# approximatedDIC.byYear2 <- function(Model, Dataset){
#   ## Compute the DIC of a model given the Observed and Expected data in a time period. The calculation of the DIC score is an approximation
#   ## based on sampled values for the posterior risks
#   # Args:
#   #   Model:  inla model already fitted which includes the field Model$risk.MCMC with a n x n.MCMC matrix with the sampled risks (n is the
#   #     number of areas and n.MCMC is the number of sampled risks for each area). This model can be obtained with the samplePosterirorRisks
#   #     function.
#   #   Dataset: data.frame with the data for the areas and years (we need the observed and the expected counts)
# 
#   mu.MCMC <- apply(Model$risk.MCMC, MARGIN=2, FUN=function(x) Dataset$exp*x)
#   # mean.deviance <- mean(apply(mu.MCMC, 2, function(x) -2*sum(log(dpois(Obs,x)))))
#   mean.deviance.byYear <- rowMeans(-2*apply(mu.MCMC,2,function(x) aggregate(log(dpois(Obs,x)),by=list(Category=Dataset$year),FUN=sum)$x))
# 
#   # deviance.mean <- -2*sum(log(dpois(Obs,apply(mu.MCMC,1,mean))))
#   deviance.mean.byYear <- -2*aggregate(log(dpois(Obs,apply(mu.MCMC,1,mean))),by=list(Category=Dataset$year),FUN=sum)$x
# 
#   pD.byYear <- mean.deviance.byYear-deviance.mean.byYear
#   DIC.byYear <- mean.deviance.byYear+pD.byYear
#   names(DIC.byYear) <- unique(Dataset$year)
#   return(DIC.byYear)
# }

approximatedDIC.byYear <- function(Model, Dataset){
  ## Compute the DIC of a model given the Observed and Expected data in a time period. The calculation of the DIC score is an approximation
  ## based on sampled values for the posterior risks
  # Args:
  #   Model:  inla model already fitted which includes the field Model$risk.MCMC with a n x n.MCMC matrix with the sampled risks (n is the
  #     number of areas and n.MCMC is the number of sampled risks for each area). This model can be obtained with the samplePosterirorRisks
  #     function.
  #   Dataset: data.frame with the data for the areas and years (we need the observed and the expected counts)
  years <- unique(Dataset$year)
  n.years <- length(years)
  DIC <- numeric(n.years)
  names(DIC) <- years
  for(i in 1:n.years){
    Obs <- Dataset[Dataset$year==years[i],"obs"]
    Exp <- Dataset[Dataset$year==years[i],"exp"]
    DIC[i] <- approximatedDIC(Model,Obs,Exp)
  }
  return(DIC)
}
  

approximatedDIC.bySims <- function(Model, Dataset, year.index=4){
  ## Compute the DIC of a model given the Observed and Expected data in a time period. The calculation of the DIC score is an approximation
  ## based on sampled values for the posterior risks
  # Args:
  #   Model:  inla model already fitted which includes the field Model$risk.MCMC with a n x n.MCMC matrix with the sampled risks (n is the
  #     number of areas and n.MCMC is the number of sampled risks for each area). This model can be obtained with the samplePosterirorRisks
  #     function.
  #   Dataset: data.frame with the data for the areas, years and simulations (we need the observed and the expected counts)
  years <- unique(Dataset$year)
  sims <- unique(Dataset$sim)
  n.sims <- length(sims)
  DIC <- numeric(n.sims)
  names(DIC) <- sims
  for(i in 1:n.sims){
    Obs <- Dataset[Dataset$year==years[year.index] & Dataset$sim==i,"obs"]
    Exp <- Dataset[Dataset$year==years[year.index] & Dataset$sim==i,"exp"]
    DIC[i] <- approximatedDIC(Model,Obs,Exp)
  }
  return(DIC)
}

# 
# Model.LCAR <- samplePosterirorRisks(Model.LCAR)
# Model.AHC  <- samplePosterirorRisks(Model.AHC)
# Model.mBoxplot.k1$model.final <- samplePosterirorRisks(Model.mBoxplot.k1$model.final)
# Model.mBoxplot.k2$model.final <- samplePosterirorRisks(Model.mBoxplot.k2$model.final)
# Model.mBoxplot.k3$model.final <- samplePosterirorRisks(Model.mBoxplot.k3$model.final)
# 
# years <- unique(Datos$year)
# n.years <- length(years)
# DIC <- array(dim=c(5,n.years))
# row.names(DIC) <- c("LCAR","AHC","DBSC.k1","DBSC.k2","DBSC.k3")
# colnames(DIC) <- years
# for(i in 1:n.years){
#   Obs <- Datos[Datos$year==years[i],"obs"]
#   Exp <- Datos[Datos$year==years[i],"exp"]
#   DIC["LCAR",i] <- approximatedDIC(Model.LCAR,Obs,Exp)
#   DIC["AHC",i] <- approximatedDIC(Model.AHC,Obs,Exp)
#   DIC["DBSC.k1",i] <- approximatedDIC(Model.mBoxplot.k1$model.final,Obs,Exp)
#   DIC["DBSC.k2",i] <- approximatedDIC(Model.mBoxplot.k2$model.final,Obs,Exp)
#   DIC["DBSC.k3",i] <- approximatedDIC(Model.mBoxplot.k3$model.final,Obs,Exp)
# }
# 
# DIC

approximatedDIC_sampled <- function(Model, Exp, Obs){
  ## Compute the DIC of a model given the Observed and Expected data in a time period. The calculation of the DIC score is an approximation
  ## based on sampled values for the posterior risks
  # Args:
  #   Model:  inla model already fitted which includes the field Model$risk.MCMC with a n x n.MCMC matrix with the sampled risks (n is the 
  #     number of areas and n.MCMC is the number of sampled risks for each area). This model can be obtained with the samplePosterirorRisks
  #     function.
  #   Obs: Observed values for each area in a specific time period
  #   Exp: Expected observations for each area in a specific time period
  # Returns:
  #   the value of the DIC score for the model and the given Observed and Expected data
  lambda=Exp*Model$summary.fitted.values$mean
  Obs.sampled <- t(sapply(Obs, function(x) rpois(1000,x)))
  mean.deviance <- mean(apply(Obs.sampled, 2, function(x) -2*sum(log(dpois(x,lambda)))))
  deviance.mean <- -2*sum(log(dpois(round(apply(Obs.sampled,1,mean)),lambda)))
  pD <- mean.deviance-deviance.mean
  DIC <- mean.deviance+pD
  return(DIC)
}

approximatedDIC_sampled2 <- function(Model, Exp, Obs){
  ## Compute the DIC of a model given the Observed and Expected data in a time period. The calculation of the DIC score is an approximation
  ## based on sampled values for the posterior risks
  # Args:
  #   Model:  inla model already fitted which includes the field Model$risk.MCMC with a n x n.MCMC matrix with the sampled risks (n is the 
  #     number of areas and n.MCMC is the number of sampled risks for each area). This model can be obtained with the samplePosterirorRisks
  #     function.
  #   Obs: Observed values for each area in a specific time period
  #   Exp: Expected observations for each area in a specific time period
  # Returns:
  #   the value of the DIC score for the model and the given Observed and Expected data
  DIC <- numeric(100)
  for(i in 1:100){
    print(i)
  # Obs.sampled <- t(sapply(Obs, function(x) rpois(100,x)))
    Obs.sampled <- sapply(Obs, function(x) rpois(1,x))
  mu.MCMC <- apply(Model$risk.MCMC, MARGIN=2, FUN=function(x) Exp*x)
    mean.deviance <- mean(apply(mu.MCMC, 2, function(x) -2*sum(log(dpois(Obs.sampled,x)))))
  deviance.mean <- -2*sum(log(dpois(Obs.sampled,apply(mu.MCMC,1,mean))))
  pD <- mean.deviance-deviance.mean
  DIC[i] <- mean.deviance+pD
  }
  return(DIC)
  
  
  
  mean.deviance <- mean(apply(Obs.sampled, 2, function(x) -2*sum(log(dpois(x,lambda)))))
  deviance.mean <- -2*sum(log(dpois(round(apply(Obs.sampled,1,mean)),lambda)))
  pD <- mean.deviance-deviance.mean
  DIC <- mean.deviance+pD
  return(DIC)
}
