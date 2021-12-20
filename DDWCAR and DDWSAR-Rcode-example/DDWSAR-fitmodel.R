##############################################################################-#
# This file presents the application of DDW in the SAR model               
##############################################################################-#

################ load data and function ########################################
`%+%` <- function(x,y) paste0(x,y)
`%>%` <- magrittr::`%>%` # pipe operator
library(ParSatscan)
# detach("package:rsatscan", unload = TRUE)
library(rsatscan)
library(spdep)
library(spatialreg)
path <- "E:\\BaiduNetdiskWorkspace\\research\\DDW\\DDWCAR and DDWSAR-Rcode-example\\"
setwd(path)
load(".\\data\\simdat(A case)_SAR.Rdat")
Rfun <- dir(".\\fun") %>% .[grepl(pattern = ".R",x = .)]
invisible(lapply(".\\fun\\"%+%Rfun,source))  # loading function
region.id <- simdat@data$code # 所有分析以该顺序为准
gallist <- poly2nb(simdat,row.names = region.id)
analysis_dat <- simdat@data
analysis_dat$y <- log(analysis_dat$cases/analysis_dat$pop)

#### detect clusters using satscan------------
analysis_dat[,c("long","lat")] <- analysis_dat[,c("X","Y")]
analysis_dat_temp <- data.frame(id = analysis_dat$code,x = analysis_dat$long,y = analysis_dat$lat,
                           pop = analysis_dat$pop,case = analysis_dat$cases)
# set scan parameters
ssoption <- list(CoordinatesType=0,
                 ReportGiniClusters="n",
                 IterativeScanMaxIterations=11,
                 IterativeScan="n",
                 ScanAreas=1)
windowsizes <- c(1:50)
parmethod <- "MCS-P" # The MCHS-P is also available
sslocation <- "D:\\Program Files\\SaTScan"
res <- satscan_msp(data = analysis_dat_temp,sizes = windowsizes,method = parmethod,
                   gallist = gallist,sslocation = sslocation,ssoption = ssoption)
clustersdat <- subset(res$Oresult$gis,P_VALUE <0.05)
clustersdat$LOC_ID <- as.integer(as.character(clustersdat[,1]))
clusterIDs <- clustersdat[[1]]; clusters <- clustersdat[[2]];RR <- clustersdat[[12]]

#### DDWSAR -------
ddwtype <- "NN"
ddwlist <- DDW(type = ddwtype,clusterID = clusterIDs,cluster = clusters,region.id = region.id,zero.policy = T)
fitmodel <- errorsarlm(y ~ x,analysis_dat,ddwlist)
fitmodel <- lagsarlm(y ~ x,analysis_dat,ddwlist)

#### the AG and other DDW could be constructed as following
ddwtype <- "NG"
ddwlist <- DDW(type = ddwtype,clusterID = clusterIDs,cluster = clusters,gal = gallist,zero.policy = T)

ddwtype <- "GN"
ddwlist <- DDW(type = ddwtype,clusterID = clusterIDs,cluster = clusters,gal = gallist,zero.policy = T)

ddwtype <- "GG"
ddwlist <- DDW(type = ddwtype,clusterID = clusterIDs,cluster = clusters,gal = gallist,zero.policy = T)

ddwtype <- "GR"
ddwlist <- DDW(type = ddwtype,clusterID = clusterIDs,cluster = clusters,RR = RR,gal = gallist,zero.policy = T)

ddwtype <- "NR"
ddwlist <- DDW(type = ddwtype,clusterID = clusterIDs,cluster = clusters,RR = RR,region.id = region.id,zero.policy = T)

ddwtype <- "AG"
ddwlist <- DDW(type = ddwtype,gal = gallist,zero.policy = T)

