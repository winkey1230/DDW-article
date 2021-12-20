##############################################################################-#
# This file presents the application of DDW in the CAR model and also gives    #
# application of DBSC and LCAR model for comparison                            #
#                                                                              #
# The four function files,i.e.,"AHC_algorithm.R" "ApproximatedDIC.R",          #
# "DBSC_algorithm.R","Model_INLA.R", comes from Santafé, G., et al.'s work,    #
# seen in: "Santafé, G., et al., Dealing with risk discontinuities to estimate #
# cancer mortality risks when the number of small areas is large. Statistical  #
# Methods in Medical Research, 2021. 30(1): p. 6-21."                          #
#                                                                              #
# ############################################################################-#

################ load data and function ########################################
`%+%` <- function(x,y) paste0(x,y)
`%>%` <- magrittr::`%>%` # pipe operator
library(ParSatscan) # The installed package is ./fun/ParSatscan.tar.gz
# detach("package:rsatscan", unload = TRUE)
library(rsatscan)
library(INLA)
library(maptools)
library(spdep)
library(RColorBrewer)
library(tmap)

path <- "E:\\BaiduNetdiskWorkspace\\research\\DDW\\DDWCAR and DDWSAR-Rcode-example\\"
setwd(path)
load(".\\data\\simdat_CAR.Rdat")
Rfun <- dir(".\\fun") %>% .[grepl(pattern = ".R",x = .)]
invisible(lapply(".\\fun\\"%+%Rfun,source))  # loading function


############### get a case data from simulation datasets #######################
scenario <- "simdat_scen2A"
simdat <- simdat_all[[scenario]]
true_risk <- simdat$true_risk
n <- ncol(W) # W represents the spatial adjacence relationship
Qs <- diag(apply(W,2,sum))-W # R
i = 1 # using the 
Datos <- simdat$basedata
Datos$obs <- simdat$simobs[,i]
Datos$SMR <- Datos$obs/Datos$exp


############# get the excess risks using LCAR, DBSC, and DDWCAR models #########
resultsi <- list() 
####### LCAR ------
R.Leroux <- diag(n)-Qs
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
data.INLA <- data.frame(Y.real=Datos[Datos$year==4,"obs"],
                        E.real=Datos[Datos$year==4,"exp"],
                        region=1:n)
formula <- Y.real ~ f(region, model="generic1", Cmatrix = R.Leroux, constr=TRUE,
                      hyper=list(prec=list(prior=sdunif),beta=list(prior=lunif, initial=0)))
model <- inla(formula, family="poisson", data=data.INLA, E=E.real,
              control.compute=list(dic=TRUE, mlik=TRUE, cpo=TRUE, waic=TRUE, config=TRUE),
              control.predictor=list(compute=TRUE, cdf=c(log(1))),
              control.inla=list(strategy="simplified.laplace"))
resultsi$LCAR <- model


####### DBSC -----
# prepare data
multivariate.data <- matrix(log(Datos[,"SMR"]+0.0001),nrow=n)
colnames(multivariate.data) <- unique(Datos$year)
clustering.inputdata <- as.data.frame(multivariate.data)
backgroundData <- rep(log(0.0001),ncol(multivariate.data))

# l=1, without background ----
partition <- getClusterPartition(x=clustering.inputdata, W=W, l=1)
cluster.model <- model.selection(cluster.partition=matrix(partition,nrow=1),
                                 strategy="simplified.laplace",
                                 Y.real=Datos[Datos$year==4,"obs"],
                                 E.real=Datos[Datos$year==4,"exp"],
                                 C=Qs, Carto=Carto.MUN)
resultsi$DBSCl1.nb <- cluster.model$model.final

# l=2, without background ----
partition <- getClusterPartition(x=clustering.inputdata, W=W, l=2)
cluster.model <- model.selection(cluster.partition=matrix(partition,nrow=1),
                                 strategy="simplified.laplace",
                                 Y.real=Datos[Datos$year==4,"obs"],
                                 E.real=Datos[Datos$year==4,"exp"],
                                 C=Qs, Carto=Carto.MUN)
resultsi$DBSCl2.nb <- cluster.model$model.final

# l=3, without background ----
partition <- getClusterPartition(x=clustering.inputdata, W=W, l=3)
cluster.model <- model.selection(cluster.partition=matrix(partition,nrow=1),
                                 strategy="simplified.laplace",
                                 Y.real=Datos[Datos$year==4,"obs"],
                                 E.real=Datos[Datos$year==4,"exp"],
                                 C=Qs, Carto=Carto.MUN)
resultsi$DBSCl3.nb <- cluster.model$model.final

# l=1, with background ----
partition <- getClusterPartition(x=clustering.inputdata, W=W, l=1, backgroundData=backgroundData)
cluster.model <- model.selection(cluster.partition=matrix(partition,nrow=1),
                                 strategy="simplified.laplace",
                                 Y.real=Datos[Datos$year==4,"obs"],
                                 E.real=Datos[Datos$year==4,"exp"],
                                 C=Qs, Carto=Carto.MUN)
resultsi$DBSCl1 <- cluster.model$model.final

# l=2, with background ----
partition <- getClusterPartition(x=clustering.inputdata, W=W, l=2, backgroundData=backgroundData)
cluster.model <- model.selection(cluster.partition=matrix(partition,nrow=1),
                                 strategy="simplified.laplace",
                                 Y.real=Datos[Datos$year==4,"obs"],
                                 E.real=Datos[Datos$year==4,"exp"],
                                 C=Qs, Carto=Carto.MUN)
resultsi$DBSCl2 <- cluster.model$model.final

# l=3, with background ----
partition <- getClusterPartition(x=clustering.inputdata, W=W, l=3, backgroundData=backgroundData)
cluster.model <- model.selection(cluster.partition=matrix(partition,nrow=1),
                                 strategy="simplified.laplace",
                                 Y.real=Datos[Datos$year==4,"obs"],
                                 E.real=Datos[Datos$year==4,"exp"],
                                 C=Qs, Carto=Carto.MUN)
resultsi$DBSCl3 <- cluster.model$model.final

###### DDWCAR -------------
# find cluster using satscan + MCS-P
scandat <- dplyr::group_by(Datos,muni) %>% dplyr::summarise_at(.vars = c("obs","exp"),.funs = sum)
scandat$SMR <- scandat$obs/scandat$exp
scandat[,c("long","lat")] <- simdat$coordinate
class(scandat) <- "data.frame"
scandat$muni <- as.character(scandat$muni)
gallist <- spdep::mat2listw(W)[[2]]
region.id <- scandat$muni
attr(gallist,"region.id") <- region.id
scandat_temp <- data.frame(id = scandat$muni,x = scandat$long,y = scandat$lat,
                           pop = scandat$exp,case = scandat$obs)
scanareas <- 3 # 3 for common diseases, 1 for rare (or much rare) diseases
ssoption <- list(CoordinatesType=0,
                 ReportGiniClusters="n",
                 IterativeScanMaxIterations=11,
                 IterativeScan="n",
                 ScanAreas=scanareas)
windowsizes <- c(1:20)
parmethod <- "MCS-P"
sslocation <- "D:\\Program Files\\SaTScan" # the intalled location of satscan software
res <- satscan_msp(data = scandat_temp,sizes = windowsizes,method = parmethod,
                   gallist = gallist,sslocation = sslocation,ssoption = ssoption)
rescluster <- subset(res$Oresult$gis,P_VALUE <0.05)
rescluster$LOC_ID <- as.character(rescluster$LOC_ID)
clusterIDs <- rescluster[,1]
clusters <- rescluster[,2]
RRs <- rescluster[,12]
if(length(clusters)==0) {
  partition <- rep(1,n)
  AGwlist <- DDW(type = "AG",gal = gallist,zero.policy = T)
  ANwlist <- listw2mat(AGwlist)
  ANwlist[upper.tri(ANwlist)] <- 1/(n-1); ANwlist[lower.tri(ANwlist)] <- 1/(n-1)
  ANwlist <- mat2listw(ANwlist)
} 
if(length(clusters)!=0){
  partition <- rescluster[,1:2]
  names(partition)[1] <- "muni"
  partition <- merge(scandat[,1:2],partition,by = "muni",all.x = T)
  sum(partition$muni == scandat$muni) # check the sequence
  partition <- partition$CLUSTER
  partition[is.na(partition)] <- max(partition,na.rm = T)+1
}

### DDW with fixed cluster effect---------------
partitionf <- as.factor(partition)
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
data.INLA <- data.frame(Y.real=Datos[Datos$year==4,"obs"],
                        E.real=Datos[Datos$year==4,"exp"],
                        region=1:n,partitionf = partitionf)
if(length(clusters)!=0)
formula <- Y.real ~ partitionf + f(region, model="generic1", Cmatrix = R.Leroux, constr=TRUE,
                                   hyper=list(prec=list(prior=sdunif),beta=list(prior=lunif, initial=0)))
if(length(clusters)!=0)
  formula <- Y.real ~ f(region, model="generic1", Cmatrix = R.Leroux, constr=TRUE,
                       hyper=list(prec=list(prior=sdunif),beta=list(prior=lunif, initial=0)))

# NNCAR_fix-------
ddwtype <- "NN"
if(length(clusters)==0) ddwlist <- ANwlist
if(length(clusters)!=0) ddwlist <- DDW(type = ddwtype,clusterID = clusterIDs,cluster = clusters,region.id = region.id,zero.policy = T)
ddwmat <- listw2mat(ddwlist) 
neighbor_n <- apply(ddwmat, 1, function(x) sum(x!=0))
ddwmat <-  diag(neighbor_n) - diag(neighbor_n)%*%ddwmat
R.Leroux <- diag(n) - ddwmat
model <- inla(formula, family="poisson", data=data.INLA, E=E.real,
              control.compute=list(dic=TRUE, mlik=TRUE, cpo=TRUE, waic=TRUE, config=TRUE),
              control.predictor=list(compute=TRUE, cdf=c(log(1))),
              control.inla=list(strategy="simplified.laplace"))
resultsi[[ddwtype%+%"CAR_fix"]] <- model


# NGCAR_fix-------
ddwtype <- "NG"
if(length(clusters)==0) ddwlist <- ANwlist
if(length(clusters)!=0) ddwlist <- DDW(type = ddwtype,clusterID = clusterIDs,cluster = clusters,gal = gallist,zero.policy = T)
ddwmat <- listw2mat(ddwlist) 
neighbor_n <- apply(ddwmat, 1, function(x) sum(x!=0))
ddwmat <-  diag(neighbor_n) - diag(neighbor_n)%*%ddwmat
R.Leroux <- diag(n) - ddwmat
model <- inla(formula, family="poisson", data=data.INLA, E=E.real,
              control.compute=list(dic=TRUE, mlik=TRUE, cpo=TRUE, waic=TRUE, config=TRUE),
              control.predictor=list(compute=TRUE, cdf=c(log(1))),
              control.inla=list(strategy="simplified.laplace"))
resultsi[[ddwtype%+%"CAR_fix"]] <- model

# GNCAR_fix-------
ddwtype <- "GN"
if(length(clusters)==0) ddwlist <- AGwlist
if(length(clusters)!=0) ddwlist <- DDW(type = ddwtype,clusterID = clusterIDs,cluster = clusters,gal = gallist,zero.policy = T)
ddwmat <- listw2mat(ddwlist) 
neighbor_n <- apply(ddwmat, 1, function(x) sum(x!=0))
ddwmat <-  diag(neighbor_n) - diag(neighbor_n)%*%ddwmat
R.Leroux <- diag(n) - ddwmat
model <- inla(formula, family="poisson", data=data.INLA, E=E.real,
              control.compute=list(dic=TRUE, mlik=TRUE, cpo=TRUE, waic=TRUE, config=TRUE),
              control.predictor=list(compute=TRUE, cdf=c(log(1))),
              control.inla=list(strategy="simplified.laplace"))
resultsi[[ddwtype%+%"CAR_fix"]] <- model

# GGCAR_fix-------
ddwtype <- "GG"
if(length(clusters)==0) ddwlist <- AGwlist
if(length(clusters)!=0) ddwlist <- DDW(type = ddwtype,clusterID = clusterIDs,cluster = clusters,gal = gallist,zero.policy = T)
ddwmat <- listw2mat(ddwlist) 
neighbor_n <- apply(ddwmat, 1, function(x) sum(x!=0))
ddwmat <-  diag(neighbor_n) - diag(neighbor_n)%*%ddwmat
R.Leroux <- diag(n) - ddwmat
model <- inla(formula, family="poisson", data=data.INLA, E=E.real,
              control.compute=list(dic=TRUE, mlik=TRUE, cpo=TRUE, waic=TRUE, config=TRUE),
              control.predictor=list(compute=TRUE, cdf=c(log(1))),
              control.inla=list(strategy="simplified.laplace"))
resultsi[[ddwtype%+%"CAR_fix"]] <- model

####### DDWCAR with random cluster effect ------
# NNCAR_random ------
ddwtype <- "NN"
if(length(clusters)==0) ddwlist <- ANwlist
if(length(clusters)!=0) ddwlist <- DDW(type = ddwtype,clusterID = clusterIDs,cluster = clusters,region.id = region.id,zero.policy = T)
ddwmat <- listw2mat(ddwlist) 
neighbor_n <- apply(ddwmat, 1, function(x) sum(x!=0))
ddwmat <-  diag(neighbor_n) - diag(neighbor_n)%*%ddwmat
cluster.model <- model.selection(cluster.partition=matrix(partition,nrow=1),
                                 strategy="simplified.laplace",
                                 Y.real=Datos[Datos$year==4,"obs"],
                                 E.real=Datos[Datos$year==4,"exp"],
                                 C=ddwmat, Carto=Carto.MUN)
resultsi[[ddwtype%+%"CAR_random"]] <- cluster.model$model.final

# NGCAR_random ------
ddwtype <- "NG"
if(length(clusters)==0) ddwlist <- ANwlist
if(length(clusters)!=0) ddwlist <- DDW(type = ddwtype,clusterID = clusterIDs,cluster = clusters,gal = gallist,zero.policy = T)
ddwmat <- listw2mat(ddwlist) 
neighbor_n <- apply(ddwmat, 1, function(x) sum(x!=0))
ddwmat <-  diag(neighbor_n) - diag(neighbor_n)%*%ddwmat
cluster.model <- model.selection(cluster.partition=matrix(partition,nrow=1),
                                 strategy="simplified.laplace",
                                 Y.real=Datos[Datos$year==4,"obs"],
                                 E.real=Datos[Datos$year==4,"exp"],
                                 C=ddwmat, Carto=Carto.MUN)
resultsi[[ddwtype%+%"CAR_random"]] <- cluster.model$model.final

# GNCAR_random ------
ddwtype <- "GN"
if(length(clusters)==0) ddwlist <- AGwlist
if(length(clusters)!=0) ddwlist <- DDW(type = ddwtype,clusterID = clusterIDs,cluster = clusters,gal = gallist,zero.policy = T)
ddwmat <- listw2mat(ddwlist) 
neighbor_n <- apply(ddwmat, 1, function(x) sum(x!=0))
ddwmat <-  diag(neighbor_n) - diag(neighbor_n)%*%ddwmat
cluster.model <- model.selection(cluster.partition=matrix(partition,nrow=1),
                                 strategy="simplified.laplace",
                                 Y.real=Datos[Datos$year==4,"obs"],
                                 E.real=Datos[Datos$year==4,"exp"],
                                 C=ddwmat, Carto=Carto.MUN)
resultsi[[ddwtype%+%"CAR_random"]] <- cluster.model$model.final

# GGCAR_random ------
ddwtype <- "GG"
if(length(clusters)==0) ddwlist <- AGwlist
if(length(clusters)!=0) ddwlist <- DDW(type = ddwtype,clusterID = clusterIDs,cluster = clusters,gal = gallist,zero.policy = T)
ddwmat <- listw2mat(ddwlist) 
neighbor_n <- apply(ddwmat, 1, function(x) sum(x!=0))
ddwmat <-  diag(neighbor_n) - diag(neighbor_n)%*%ddwmat
cluster.model <- model.selection(cluster.partition=matrix(partition,nrow=1),
                                 strategy="simplified.laplace",
                                 Y.real=Datos[Datos$year==4,"obs"],
                                 E.real=Datos[Datos$year==4,"exp"],
                                 C=ddwmat, Carto=Carto.MUN)
resultsi[[ddwtype%+%"CAR_random"]] <- cluster.model$model.final

##### risk maps ------
LS <- unlist(lapply(resultsi, function(x) -sum(log(x$cpo$cpo))))
results <- resultsi[c("LCAR","DBSCl2.nb","GGCAR_random")] # select the model to map risk
mean.risk <- lapply(results, function(x) x$summary.fitted.values$'0.5quant')
mean.risk$true.risk <- true_risk
mean.risk$muni <- Carto.MUN$muni
mean.risk$SMR <- Datos[Datos$year==4,"SMR"]

Carto <- merge(Carto.MUN,mean.risk)
paleta <- brewer.pal(8,"RdYlGn")[8:1]
values <- c(-Inf,0.67,0.77,0.91,1,1.10,1.30,1.50,Inf)
model.names <- names(results)

pdf("results.pdf")
tmap_mode("plot")
Map <- tm_shape(Carto) + 
  tm_polygons(col=c("true.risk","SMR",model.names),
              palette=paleta, title="", legend.show=T, border.alpha=0,
              legend.reverse=T, style="fixed", breaks=values, interval.closure="left",
              labels=c("[0-0.67)","[0.67-0.77)","[0.77-0.91)","[0.91-1)",
                       "[1-1.10)","[1.10-1.30)","[1.30-1.50)","[1.50-Inf)")) +
  tm_grid(n.x=2, n.y=2, alpha=0.2, labels.format=list(scientific=T), labels.inside.frame=F, labels.col="white") +
  tm_layout(main.title="Mean Risks", main.title.position=0.2,
            panel.labels=c("true.risk","SMR",model.names),
            legend.outside=T, legend.outside.position="right", legend.frame=F, legend.outside.size=0.2,
            outer.margins=c(0.02,0.01,0.02,0.01)) + 
  tm_facets(ncol=4, nrow=2)
print(Map)
dev.off()

##### plot the cluster distribution -------------------
cluster.risk <- data.frame(scenario1 = simdat_all$simdat_scen1A$true_risk,
                           scenario2 = simdat_all$simdat_scen2A$true_risk,
                           scenario3 = simdat_all$simdat_scen3A$true_risk)

cluster.risk$muni <- Carto.MUN$muni
Carto <- merge(Carto.MUN,cluster.risk)
paleta <- brewer.pal(8,"RdYlGn")[8:1]
values <- c(-Inf,0.67,0.77,0.91,1,1.10,1.30,1.50,Inf)

pdf("results.pdf")
tmap_mode("plot")
Map <- tm_shape(Carto) + 
  tm_polygons(col=c("scenario1","scenario2","scenario3"),
              palette=paleta, title="", legend.show=T, border.alpha=0,
              legend.reverse=T, style="fixed", breaks=values, interval.closure="left",
              labels=c("[0-0.67)","[0.67-0.77)","[0.77-0.91)","[0.91-1)",
                       "[1-1.10)","[1.10-1.30)","[1.30-1.50)","[1.50-Inf)")) +
  tm_grid(n.x=2, n.y=2, alpha=0.2, labels.format=list(scientific=T), labels.inside.frame=F, labels.col="white") +
  tm_layout(main.title="The clustering distribution", main.title.position=0.2,
            panel.labels=c("scenario1(No clusters)","scenario2","scenario3"),
            legend.outside=T, legend.outside.position="right", legend.frame=F, legend.outside.size=0.2,
            outer.margins=c(0.02,0.01,0.02,0.01),between.margin = -3) + 
  tm_facets(ncol=3, nrow=1)
print(Map)
dev.off()



