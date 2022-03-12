library("spatstat.data",lib.loc="/home/au591455/Rstuff/library")
library("spatstat.geom",lib.loc="/home/au591455/Rstuff/library")
library("spatstat.random",lib.loc="/home/au591455/Rstuff/library")
library("spatstat.core",lib.loc="/home/au591455/Rstuff/library")
library("spatstat.linnet",lib.loc="/home/au591455/Rstuff/library")
library("spatstat",lib.loc="/home/au591455/Rstuff/library")
install.packages(pkgs="tictoc",lib ="/Users/simon/Desktop/TestR",type="source",destdir ="/Users/simon/Desktop/TestR" ) 


fileConn<-file("C:/Users/simon/Desktop/TestR/output.txt")
current_log = readLines(fileConn)
writeLines(c(current_log,"Hello World"), fileConn)
current_log = readLines(fileConn)
Simon = "CAT"
writeLines(c(current_log,paste(5/10,Simon , "Lolk")), fileConn)

close(fileConn)



print(1+1)
setwd("/home/au591455/Rstuff/Results") 
SAVEM = matrix(c(1,2,3,4),2,2)

save(SAVEM,file = "SAVEM.Rdata")

SAVEM = rbind(SAVEM,c(5,6))

save(SAVEM,file = "SAVEM.Rdata")

SAVEM = c(1,2,3)

load("SAVEM.Rdata")
