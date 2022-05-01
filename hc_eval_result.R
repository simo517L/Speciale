liblocation=NULL

library("spatstat.data",lib.loc=liblocation )
library("spatstat.geom",lib.loc=liblocation )
library("spatstat.random",lib.loc=liblocation )
library("spatstat.core",lib.loc=liblocation )
library("spatstat.linnet",lib.loc=liblocation )
library("spatstat",lib.loc=liblocation )
library(parallel)
library("codetools",lib.loc=liblocation )
library("iterators",lib.loc=liblocation )
library("tictoc",lib.loc=liblocation )
library("foreach",lib.loc=liblocation )
library("doParallel",lib.loc=liblocation )
library("cluster",lib.loc=liblocation)
library(utils)

library("ppMeasures",lib.loc=liblocation)

library("foreign",lib.loc=liblocation)
library("maps",lib.loc=liblocation)
library("shapefiles",lib.loc=liblocation)
library("sp",lib.loc=liblocation)
library("fossil",lib.loc=liblocation)
setwd("C:/Users/simon/Desktop/TestR") 
Data =  readRDS(file = "DataPPP.Rdata")
Matern4a = readRDS(file = "Matern_a.Rdata")
Matern4b = readRDS(file = "Matern_b.Rdata")
Clust4a = readRDS(file = "Clust_a.Rdata")
Clust4b = readRDS(file = "Clust_b.Rdata")
Clust4c = readRDS(file = "Clust_c.Rdata")
poistest  = readRDS(file = "poisPPP.Rdata")

hc_hard_nn =rowMeans(readRDS("eval1.Rdata"))
hc_easy_nn =rowMeans(readRDS("eval2.Rdata"))
hc_big_nn =rowMeans(readRDS("eval3.Rdata"))

hc_hard_T4 =rowMeans(readRDS("eval1T4.Rdata"))
hc_easy_T4 =rowMeans(readRDS("eval2T4.Rdata"))
hc_big_T4 =rowMeans(readRDS("eval3T4.Rdata"))

hc_hard_T6 =rowMeans(readRDS("eval1T6.Rdata"))
hc_easy_T6 =rowMeans(readRDS("eval2T6.Rdata"))
hc_big_T6 =rowMeans(readRDS("eval3T6.Rdata"))

hc_hard_T9 =rowMeans(readRDS("eval1T9.Rdata"))
hc_easy_T9 =rowMeans(readRDS("eval2T9.Rdata"))
hc_big_T9 =rowMeans(readRDS("eval3T9.Rdata"))

hc_hard_T12 =rowMeans(readRDS("eval1T12.Rdata"))
hc_easy_T12 =rowMeans(readRDS("eval2T12.Rdata"))
hc_big_T12 =rowMeans(readRDS("eval3T12.Rdata"))

hc_hard_st=rowMeans(readRDS("eval1st.Rdata"))
hc_easy_st =rowMeans(readRDS("eval2st.Rdata"))
hc_big_st =rowMeans(readRDS("eval3st.Rdata"))

M_hard = matrix(c(hc_hard_nn,hc_hard_T4,hc_hard_T6,hc_hard_T9,hc_hard_T12,hc_hard_st),6,9,byrow = T)
signif(M_hard, digits = 5)

signif(c(hc_hard_nn,hc_hard_T4,hc_hard_T6,hc_hard_T9,hc_hard_T12,hc_hard_st), digits = 4)



M_easy = matrix(c(hc_easy_nn,hc_easy_T4,hc_easy_T6,hc_easy_T9,hc_easy_T12,hc_easy_st),6,9,byrow = T)
signif(M_easy, digits = 5)



M_big= matrix(c(hc_big_nn,hc_big_T4,hc_big_T6,hc_big_T9,hc_big_T12,hc_big_st),6,9,byrow = T)
signif(M_big, digits = 5)
