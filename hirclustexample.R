#Simuler grupper af data fra forskellige fordelinger, fx. 3 forskellige normalfordelinger, og gem dem i en samlet dataframe, med gruppenummer. Afprøv clustering algoritmer som du finder i R på dine simulerede eksempler.

r1= rnorm(20)
r2= rnorm(20, mean = 20)
r3= rnorm(20, mean = 30 , sd = 5)
label = c(rep("r1",20),rep("r2",20),rep("r3",20))

df = as.data.frame(c(r1,r2,r3))

d <- dist(df)

hc = hclust(d)

plot(hc)

res.coph <- cophenetic(hc)
cor(res.coph,d)

mean(cutree(hc,k=3) ==  c(rep(1,20),rep(2,20),rep(3,20)))



library(cluster)

df = as.data.frame(c(r1,r2,r3))

d <- dist(df)
hctest= agnes(d)
plot(hctest)


hc2= agnes(d)
cor(cophenetic(hc2),d)
plot(hc2)

hc3= agnes(d,method = "single")
cor(cophenetic(hc3),d)
plot(hc3)

hc4= agnes(d,method = "complete")
cor(cophenetic(hc4),d)
plot(hc4)

mean(cutree(hc2,k=3) ==  c(rep(1,20),rep(2,20),rep(3,20)))

hcd1 = diana(d)
cor(cophenetic(hcd1),d)

plot(hcd1)




library(colorspace)
library(dendextend)

make_dend = function(hc,k,label){
n = length( hc$order)
dend <- as.dendrogram(hc)
dend <- rotate(dend, 1:n)
dend <- color_branches(dend, k=k)

labels_colors(dend) <-
  rainbow_hcl(k)[sort_levels_values(
    as.numeric(as.factor(label))[order.dendrogram(dend)]
  )]

labels(dend) <- paste(as.character(label)[order.dendrogram(dend)],
                      "(",labels(dend),")", 
                      sep = "")
#dend <- hang.dendrogram(dend,hang_height=0.1)
dend <- set(dend, "labels_cex", 0.5)
return(dend)
}

dend = make_dend(hc,k=3,label)

plot(dend,horiz =  TRUE,  nodePar = list(cex = .007))

dend = make_dend(hc2,k=3,label)
plot(dend,horiz =  TRUE,  nodePar = list(cex = .007))

dend = make_dend(hc3,k=3,label)
plot(dend,horiz =  TRUE,  nodePar = list(cex = .007))

dend = make_dend(hc4,k=3,label)
plot(dend,horiz =  TRUE,  nodePar = list(cex = .007))



dend = make_dend(hcd1,k=3,label)
plot(dend,horiz =  TRUE,  nodePar = list(cex = .007))





r1= rnorm(20,mean=10)
r2= rpois(20,lambda = 30)
r3= rexp(20, rate = 1/20)

df = as.data.frame(c(r1,r2,r3))

d <- dist(df)

hc2= agnes(d)
cor(cophenetic(hc2),d)
plot(hc2)

dend = make_dend(hc2,k=3,label)
plot(dend,horiz =  TRUE,  nodePar = list(cex = .007))

#hvor godt virker algoritmerne?
#  Hvad for en kriterium bruger du for at bedømme, om de virker godt?

# vizulaisering, cophenetic correlation coeffcient 

#  hvordan afhænger resultater fra antal af forskellige grupper, og parametrene af de til bunds liggende fordelinge


