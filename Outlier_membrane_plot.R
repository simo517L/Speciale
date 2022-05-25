library(spatstat)


setwd("C:/Users/simon/Desktop/TestR")
memb_LOF_control_T =mean( readRDS("memb_LOF_control_T.Rdata") <= 0.05)
memb_LOF_acid_T = mean(readRDS("memb_LOF_control_T.Rdata") <= 0.05)
memb_LOF_rotenone_T = mean(readRDS("memb_LOF_control_T.Rdata") <= 0.05)

memb_LOF_control_NP =mean( readRDS("memb_LOF_control_NP.Rdata") <= 0.05)
memb_LOF_acid_NP = mean(readRDS("memb_LOF_acid_NP.Rdata") <= 0.05)
memb_LOF_rotenone_NP = mean(readRDS("memb_LOF_rotenone_NP.Rdata") <= 0.05)

memb_LOF_control_ST =mean( readRDS("memb_LOF_control_ST.Rdata") <= 0.05)
memb_LOF_acid_ST = mean(readRDS("memb_LOF_acid_ST.Rdata") <= 0.05)
memb_LOF_rotenone_ST = mean(readRDS("memb_LOF_rotenone_ST.Rdata") <= 0.05)

memb_perm_control =mean( readRDS("memb_perm_control.Rdata") <= 0.05)
memb_perm_acid = mean(readRDS("memb_perm_acid.Rdata") <= 0.05)
memb_perm_rotenone = mean(readRDS("memb_perm_rotenone.Rdata") <= 0.05)

c(memb_LOF_control_T,memb_LOF_acid_T,memb_LOF_rotenone_T)

MM = cbind( c("LOF (T-stat)", "LOF (NP)" , "LOF (ST)", "Permutation" ),MM)
MM = data.frame(MM)
ft <- flextable( MM)
ft <- delete_part(ft, part = "header")
ft <- add_header_row(ft,
                     values = c("control", "acid" , "rotenone" ),top = T)
ft <- bold(ft, j=1, bold = TRUE, part = "body")
#ft <- autofit(ft, add_w = 0, add_h = 0)
ft <- theme_vanilla(ft)
ft <- vline(ft, j = 1, border =  fp_border() , part = "body")
ft <- border_inner_v(ft, border =  fp_border(), part = "header")
ft <- vline(ft, j =dim(MM)[2], border =  fp_border() , part = "all")
ft <- vline_left(ft, border = fp_border() )
ft <- fontsize(ft, size = 8, part = "all")
ft <- align(ft,align = "center")
ft <- align(ft,i=1,align = "center",part="header")
ft <- align(ft,j=1,align = "center",part="header")
ft <- align(ft,i=2,align = "center",part="header")



