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
library(flextable)
library(magrittr)
library(data.table)
library(scales)
library(officer)

tablefunc = function(MM){
 MM = cbind( c("Neares Point","T (sq=4)","T (sq=6)","T (sq=9)","T (sq=12)","Spike-time" ),MM)
  MM = data.frame(MM)
  ft <- flextable( MM)
  ft <- delete_part(ft, part = "header")
  ft <- add_header_row(ft,
                       values =c("",rep(c("Single","Average","Complete"),3 )),top = F)
  ft <- add_header_row(ft,colwidths = c(1,3, 3,3),
                       values = c(" ","Agglomerative coefficient", "F-measure","Adj. Rand index"),top = T
  )
  ft <- bold(ft, j=1, bold = TRUE, part = "body")
  ft <- align(ft,j=1,align = "center")
  ft <- align(ft,i=1,align = "center",part="header")
  ft <- align(ft,i=2,align = "center",part="header")
  #ft <- autofit(ft, add_w = 0, add_h = 0)
  ft <- theme_vanilla(ft)
  ft <- vline(ft, j = 1, border =  fp_border() , part = "body")
  ft <- border_inner_v(ft, border =  fp_border(), part = "header")
  ft <- vline(ft, j =dim(MM)[2], border =  fp_border() , part = "all")
  ft <- vline_left(ft, border = fp_border() )
  
  return(ft)
}

M_hard = matrix(c(hc_hard_nn,hc_hard_T4,hc_hard_T6,hc_hard_T9,hc_hard_T12,hc_hard_st),6,9,byrow = T)
M_hard = signif(M_hard, digits = 6)
tablefunc(M_hard)




signif(M_hard, digits = 5)

signif(c(hc_hard_nn,hc_hard_T4,hc_hard_T6,hc_hard_T9,hc_hard_T12,hc_hard_st), digits = 4)



M_easy = matrix(c(hc_easy_nn,hc_easy_T4,hc_easy_T6,hc_easy_T9,hc_easy_T12,hc_easy_st),6,9,byrow = T)
M_easy  = signif(M_easy, digits = 6)
tablefunc(M_easy)


M_big= matrix(c(hc_big_nn,hc_big_T4,hc_big_T6,hc_big_T9,hc_big_T12,hc_big_st),6,9,byrow = T)
M_big = signif(M_big, digits = 6)
tablefunc(M_big)


M_membrane = matrix(c(rowMeans(readRDS("evalmembraneNP.Rdata"))
,rowMeans(readRDS("evalmembraneT.Rdata"))
,rowMeans(readRDS("evalmembraneST.Rdata"))
),3,9,byrow = T)
M_membrane = signif(M_membrane, digits = 6)
M_membrane  = cbind( c("Neares Point","T (sq=6)","Spike-time" ),M_membrane )
M_membrane = data.frame(M_membrane )
ft <- flextable( M_membrane )
ft <- delete_part(ft, part = "header")
ft <- add_header_row(ft,
                     values =c("",rep(c("Single","Average","Complete"),3 )),top = F)
ft <- add_header_row(ft,colwidths = c(1,3, 3,3),
                     values = c(" ","Agglomerative coefficient", "F-measure","Adj. Rand index"),top = T
)
ft <- bold(ft, j=1, bold = TRUE, part = "body")
ft <- align(ft,j=1,align = "center")
ft <- align(ft,i=1,align = "center",part="header")
ft <- align(ft,i=2,align = "center",part="header")
#ft <- autofit(ft, add_w = 0, add_h = 0)
ft <- theme_vanilla(ft)
ft <- vline(ft, j = 1, border =  fp_border() , part = "body")
ft <- border_inner_v(ft, border =  fp_border(), part = "header")
ft <- vline(ft, j =dim(M_membrane )[2], border =  fp_border() , part = "all")
ft <- vline_left(ft, border = fp_border() )
ft


