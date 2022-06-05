library(xtable)
#Loads Data 
setwd("C:/Users/simon/Desktop/SpecialeProject/Speciale/Data")# Remember to add your own path for the data 

hc_hard_np =rowMeans(readRDS("eval1.Rdata"))
hc_easy_np =rowMeans(readRDS("eval2.Rdata"))
hc_big_np =rowMeans(readRDS("eval3.Rdata"))

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


# Define a function which uses the xtable function to generate a table for a latex file 
tablefunc = function(MM){
  MM = rbind(  rep(c("Single","Average","Complete"),3 ),MM)
  MM = data.frame(MM,row.names =  c("","Nearest Point","T (sq=4)","T (sq=6)","T (sq=9)","T (sq=12)","Spike-time" ))
  tableMM = xtable(MM,auto=T,align = "llll|lll|lll")
  addtorow <- list()
  addtorow$pos <- list(0)
  addtorow$command <- paste0(paste0('& \\multicolumn{3}{c}{', c("Agglomerative coefficient", "F-measure","Adj. Rand index"), '}', collapse=''), '\\\\')
  print(tableMM,add.to.row=addtorow,include.colnames = FALSE,
        hline.after = c(0,1,nrow(tableMM)),scalebox=0.8)
}


M_easy = matrix(c(hc_easy_np,hc_easy_T4,hc_easy_T6,hc_easy_T9,hc_easy_T12,hc_easy_st),6,9,byrow = T)
M_easy = signif(M_easy, digits = 4)
tablefunc(M_easy)

M_hard = matrix(c(hc_hard_np,hc_hard_T4,hc_hard_T6,hc_hard_T9,hc_hard_T12,hc_hard_st),6,9,byrow = T)
M_hard = signif(M_hard, digits = 4)
tablefunc(M_hard)

M_big= matrix(c(hc_big_np,hc_big_T4,hc_big_T6,hc_big_T9,hc_big_T12,hc_big_st),6,9,byrow = T)
M_big = signif(M_big, digits = 4)
tablefunc(M_big)

# We do not define a new function for the membrane data, as we only need to make one table
M_membrane = matrix(c(rowMeans(readRDS("evalmembraneNP.Rdata"))
,rowMeans(readRDS("evalmembraneT.Rdata"))
,rowMeans(readRDS("evalmembraneST.Rdata"))
),3,9,byrow = T)

M_membrane = signif(M_membrane, digits = 4)
M_membrane = rbind(  rep(c("Single","Average","Complete"),3 ),M_membrane)
M_membrane = data.frame(M_membrane,row.names =  c("","Nearest Point","T (sq=6)","Spike-time" ))
tablemembrane = xtable(M_membrane,auto=T,align = "llll|lll|lll")
addtorow <- list()
addtorow$pos <- list(0)
addtorow$command <- paste0(paste0('& \\multicolumn{3}{c}{', c("Agglomerative coefficient ", "F-measure","Adj. Rand index"), '}', collapse=''), '\\\\')
print(tablemembrane,add.to.row=addtorow,include.colnames = FALSE,
      hline.after = c(0,1,nrow(tablemembrane)),scalebox=0.8)



# We load the data from the Divisive alg. 
hc_hard_np_diana =rowMeans(readRDS("eval1diana.Rdata"))
hc_easy_np_diana =rowMeans(readRDS("eval2diana.Rdata"))
hc_big_np_diana =rowMeans(readRDS("eval3diana.Rdata"))

hc_hard_T4_diana =rowMeans(readRDS("eval1T4diana.Rdata"))
hc_easy_T4_diana =rowMeans(readRDS("eval2T4diana.Rdata"))
hc_big_T4_diana =rowMeans(readRDS("eval3T4diana.Rdata"))

hc_hard_T6_diana =rowMeans(readRDS("eval1T6diana.Rdata"))
hc_easy_T6_diana =rowMeans(readRDS("eval2T6diana.Rdata"))
hc_big_T6_diana =rowMeans(readRDS("eval3T6diana.Rdata"))

hc_hard_T9_diana =rowMeans(readRDS("eval1T9diana.Rdata"))
hc_easy_T9_diana =rowMeans(readRDS("eval2T9diana.Rdata"))
hc_big_T9_diana =rowMeans(readRDS("eval3T9diana.Rdata"))

hc_hard_T12_diana =rowMeans(readRDS("eval1T12diana.Rdata"))
hc_easy_T12_diana =rowMeans(readRDS("eval2T12diana.Rdata"))
hc_big_T12_diana =rowMeans(readRDS("eval3T12diana.Rdata"))

hc_hard_st_diana=rowMeans(readRDS("eval1st_diana.Rdata"))
hc_easy_st_diana =rowMeans(readRDS("eval2st_diana.Rdata"))
hc_big_st_diana =rowMeans(readRDS("eval3st_diana.Rdata"))

M_easy_diana = matrix(c(hc_easy_np_diana,hc_easy_T4_diana,hc_easy_T6_diana,hc_easy_T9_diana,hc_easy_T12_diana,hc_easy_st_diana),6,3,byrow = T)

M_hard_diana = matrix(c(hc_hard_np_diana,hc_hard_T4_diana,hc_hard_T6_diana,hc_hard_T9_diana,hc_hard_T12_diana,hc_hard_st_diana),6,3,byrow = T)

M_big_diana= matrix(c(hc_big_np_diana,hc_big_T4_diana,hc_big_T6_diana,hc_big_T9_diana,hc_big_T12_diana,hc_big_st_diana),6,3,byrow = T)

# Define a function which uses the xtable function to generate a table for a latex fill, but this time for the Divisive alg. 
tablefunc_diana = function(MM){
  MM = signif(MM , digits = 4)
  MM = rbind(  c("Divisive coefficient", "F-measure","Adj. Rand index"),MM)
  MM = data.frame(MM,row.names =  c("","Nearest Point","T (sq=4)","T (sq=6)","T (sq=9)","T (sq=12)", "Spike-time" ))
  tableMM = xtable(MM,auto=T,align = "llll")
  print(tableMM,include.colnames = FALSE,
        hline.after = c(0,1,nrow(tableMM)),scalebox=0.8)
}

tablefunc_diana(M_easy_diana)

tablefunc_diana(M_hard_diana)

tablefunc_diana(M_big_diana)


M_membrane_diana = matrix(c(rowMeans(readRDS("evalmembraneNP_diana.Rdata"))
                      ,rowMeans(readRDS("evalmembraneT_diana.Rdata"))
                      ,rowMeans(readRDS("evalmembraneST_diana.Rdata")) ),3,3,byrow = T)
# We do not define a new function for the membrane data, as we only need to make one table
M_membrane_diana = signif(M_membrane_diana , digits = 4)
M_membrane_diana = rbind(  c("Divisive coefficient", "F-measure","Adj. Rand index"),M_membrane_diana)
M_membrane_diana = data.frame(M_membrane_diana,row.names =  c("","Nearest Point","T (sq=6)", "Spike-time" ))
tableMM = xtable(M_membrane_diana,auto=T,align = "llll")
print(tableMM,include.colnames = FALSE,
      hline.after = c(0,1,nrow(tableMM)),scalebox=0.8)
