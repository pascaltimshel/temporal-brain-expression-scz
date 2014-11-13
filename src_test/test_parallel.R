
############################## doParallel ###############################
#http://stackoverflow.com/questions/1395309/how-to-make-r-use-all-processors

library(parallel)
library(foreach)
library(doParallel)

#setup parallel backend to use 8 processors
cl<-makeCluster(8)
registerDoParallel(cl)

cl <- makeCluster(detectCores() - 1)
registerDoParallel(cl, cores = detectCores() - 1)

data = foreach(i = 1:length(filenames), .packages = c("ncdf","chron","stats"),
               .combine = rbind) %dopar% {
                 try({
                   # your operations; line 1...
                   # your operations; line 2...
                   # your output
                 })
}


########################### doParallel - simple - using the multicore-like functionality #########
library(doParallel)
registerDoParallel(cores=2)
x<-foreach(i=1:3) %dopar% sqrt(i)
x

########################### doMC #########################
#import packages
library(foreach)
library(doMC)

registerDoParallel(cl)

?registerDoMC
?getDoParWorkers
getDoParWorkers() #---> get number of workers/CPU

foreach(i=1:10) %dopar% {
  #loop contents here
  
}
