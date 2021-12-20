#' @description This function implement satscan using mcs-p or mchs-p based on 
#' ParSatscan package.
#' @param data A dataframe, must include the five variables: id, x, y,pop,case.
#' @param sizes A vector. The candidated maximum window sizes.
#' @param gallist When method = MCS-P, the gallist is not used.When method = MCHS-P
#' gallist must be specified. The details of gallist are seen in function selpar_mchsp.
#' @param ssoption A list for the settings of rsatscan.
#' @param location The location used to save the satscan file.
#' @param others The other parameters are seen in function selpar_mchsp or selpar_mcsp
satscan_msp <- function(data,sizes,gallist = NULL,method = "MCS-P",
                        location = "./temp",sslocation,ssoption,verbose=F,
                        ssbatchfilename = "SaTScanBatch64",...){
  scandat <- as.data.frame(data)
  if(any(!(c("id", "x", "y","pop","case") %in% names(data)))) 
    stop("data must include the five variables: id, x, y,pop,case")
  dir.create(location,showWarnings = F)
  write.geo(scandat[,c("id","x","y")],filename = "Spain",location = location)
  write.cas(scandat[,c("id","case")],filename = "Spain",location = location)
  write.pop(data.frame(scandat[,"id"],"unspecified",scandat[,"pop"]),filename = "Spain",location = location)
  invisible(ss.options(reset = T))
  #select parameters excluding the maximum window sizes, 3.4 has no effect.
  ssoptions <- list(CaseFile= paste0(location,"/Spain.cas"),
                    PrecisionCaseTimes=0,
                    PopulationFile=paste0(location,"/Spain.pop"),
                    CoordinatesFile=paste0(location,"/Spain.geo"))
  ssoptions <- c(ssoptions,ssoption)
  ss.options(ssoptions)
  if(method == "MCS-P") 
    xx <- selpar_mcsp(pop = scandat$pop,case = scandat$case,sizes = sizes,
                       id = scandat$id,sslocation = sslocation,location = location,
                       verbose=F,ssbatchfilename = ssbatchfilename,...)
  else if(method == "MCHS-P") 
    xx <- selpar_mchsp(pop = scandat$pop,case = scandat$case,sizes = sizes,
                       id = scandat$id,sslocation = sslocation,location = location,
                       gallist = gallist,verbose=F,ssbatchfilename = ssbatchfilename,...)
  else stop("Method must be 'MCS-P' or 'MCHS-P'")
  xx
}





