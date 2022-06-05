library(spatstat)
library(xtable)
path = "C:/Users/simon/Desktop/SpecialeProject/Speciale/Data" # Remember to set your own path 
setwd(path) 
#We load the data
memb_LOF_control_T =mean( readRDS("memb_LOF_control_T.Rdata")[1:500] <= 0.05)
memb_LOF_acid_T = mean(readRDS("memb_LOF_acid_T.Rdata")[1:500] <= 0.05)
memb_LOF_rotenone_T = mean(readRDS("memb_LOF_rotenone_T.Rdata")[1:500] <= 0.05)

memb_LOF_control_NP =mean( readRDS("memb_LOF_control_NP.Rdata")[1:500] <= 0.05)
memb_LOF_acid_NP = mean(readRDS("memb_LOF_acid_NP.Rdata")[1:500] <= 0.05)
memb_LOF_rotenone_NP = mean(readRDS("memb_LOF_rotenone_NP.Rdata")[1:500] <= 0.05)

memb_LOF_control_ST =mean( readRDS("memb_LOF_control_ST.Rdata") <= 0.05)
memb_LOF_acid_ST = mean(readRDS("memb_LOF_acid_ST.Rdata") <= 0.05)
memb_LOF_rotenone_ST = mean(readRDS("memb_LOF_rotenone_ST.Rdata") <= 0.05)

memb_perm_control =mean( readRDS("memb_perm_control.Rdata")[1:500] <= 0.05)
memb_perm_acid = mean(readRDS("memb_perm_acid.Rdata")[1:500] <= 0.05)
memb_perm_rotenone = mean(readRDS("memb_perm_rotenone.Rdata")[1:500] <= 0.05)

# We place the data into a latex table
MM = matrix(c(memb_perm_control,memb_perm_acid,memb_perm_rotenone,
              memb_LOF_control_T,memb_LOF_acid_T,memb_LOF_rotenone_T,
              memb_LOF_control_ST,memb_LOF_acid_ST,memb_LOF_rotenone_ST,
              memb_LOF_control_NP,memb_LOF_acid_NP,memb_LOF_rotenone_NP),3,4)

MM = rbind( c( "Permutation" , "LOF (T-stat)", "LOF (ST)" , "LOF (NP)"),MM)

MM = data.frame(MM,row.names = c(" ","control", "acid" , "rotenone" ))
tableMM = xtable(MM,auto=T)
print(tableMM,include.colnames = FALSE,
      hline.after = c(0,1,nrow(tableMM)))



