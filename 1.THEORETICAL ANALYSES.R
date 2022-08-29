#------------------------------#
#---| THEORETICAL ANALYSES |---#
#------------------------------#

library(readxl)
library(parallel)

# Prepare cluster
ncores<-detectCores(logical=T)
cl<-makeCluster(ncores)

# Import pond coordinates
Pond_xy<-read_excel("Pond_xy.xlsx")

# Transform pond coordinates into matrix form
Pond_xy<-as.matrix(Pond_xy)
rownames(Pond_xy)<-Pond_xy[,1]
Pond_xy<-Pond_xy[,-1]

# Import database
Dataset<-read_excel("Dataset.xlsx")
Dataset<-as.data.frame(Dataset)

# Character vector with those orders included in the analysis
orders<-unique(Dataset[,5])
orders<-orders[c(1,2,4,5,6,9,12,13,22,28,29)]
orders<-sort(orders)

# Estimate mean body size for the orders selected
bs.orders<-bz_x_taxa(taxa = orders, M = Dataset, taxa_en = 5, volumen_en = 4)


#---| (1) SIMULATIONS INCLUDING ALL POSSIBLE COMBINATIONS OF THE THREE SCALING RELATIONSHIPS (Figure 2 and Supporting Information Table SM2) |---#

# DISPERSAL + POOL + DENSITY
simu.all4.N50.2_0.75_0.4<-body.size.diversity.4(bs.taxa=bs.orders, xy=Pond_xy, xy.sampled = T,      
                                               b.disp=0.5, b.pool=-0.4, b.local=-0.75, J.local=50, 
                                               m.pool=0.00001, m.neighbour=0.01, taxa=orders)       
# NONE
simu.none4.N50.2_0.75_0.4<-body.size.diversity.4(bs.taxa=bs.orders, xy=Pond_xy, xy.sampled = T,     
                                                 b.disp=0, b.pool=0, b.local=0, J.local=50,
                                                 m.pool=0.00001, m.neighbour=0.01, taxa=orders)     
# DISPERSAL
simu.disp4.N50.2_0.75_0.4<-body.size.diversity.4(bs.taxa=bs.orders, xy=Pond_xy, xy.sampled = T,      
                                                 b.disp=0.5, b.pool=0, b.local=0, J.local=50,
                                                 m.pool=0.00001, m.neighbour=0.01, taxa=orders)     
# POOL
simu.pool4.N50.2_0.75_0.4<-body.size.diversity.4(bs.taxa=bs.orders, xy=Pond_xy, xy.sampled = T,     
                                                 b.disp=0, b.pool=-0.4, b.local=0, J.local=50,
                                                 m.pool=0.00001, m.neighbour=0.01, taxa=orders)     
# DENSITY
simu.density4.N50.2_0.75_0.4<-body.size.diversity.4(bs.taxa=bs.orders, xy=Pond_xy, xy.sampled = T,  
                                                    b.disp=0, b.pool=0, b.local=-0.75, J.local=50,
                                                    m.pool=0.00001, m.neighbour=0.01, taxa=orders)  
# DISPERSAL + POOL 
simu.disp.pool4.N50.2_0.75_0.4<-body.size.diversity.4(bs.taxa=bs.orders, xy=Pond_xy, xy.sampled = T,  
                                                      b.disp=0.5, b.pool=-0.4, b.local=0, J.local=50, 
                                                      m.pool=0.00001, m.neighbour=0.01, taxa=orders)        
# DISPERSAL + DENSITY
simu.disp.density4.N50.2_0.75_0.4<-body.size.diversity.4(bs.taxa=bs.orders, xy=Pond_xy, xy.sampled = T, 
                                                         b.disp=0.5, b.pool=0, b.local=-0.75, J.local=50, 
                                                         m.pool=0.00001, m.neighbour=0.01, taxa=orders) 
# POOL + DENSITY
simu.pool.density4.N50.2_0.75_0.4<-body.size.diversity.4(bs.taxa=bs.orders, xy=Pond_xy, xy.sampled = T, 
                                                         b.disp=0, b.pool=-0.4, b.local=-0.75, J.local=50, 
                                                         m.pool=0.00001, m.neighbour=0.01, taxa=orders) 


#---| (2) PARAMETER SPACE (Supporting information Figure SM3) |---# 

par.space.1<-espacio.parametros(diversity.metrics = c("alpha.coalescent.std", "beta.add.coalescent.std", "gamma.coalescent"), disp = "k", density = "j", pool = "i", gradient.i = c(0, -0.1, -0.2, -0.3), gradient.j = round(seq(from=0, to=-2, by=-0.05),2), gradient.k = round(seq(from=0, to=2, by=0.1),1), N.replicates = 5)

par.space.2<-espacio.parametros(diversity.metrics = c("alpha.coalescent.std", "beta.add.coalescent.std", "gamma.coalescent"), disp = "k", density = "j", pool = "i", gradient.i = c(-0.4, -0.6), gradient.j = round(seq(from=0, to=-2, by=-0.05),2), gradient.k = round(seq(from=0, to=2, by=0.1),1), N.replicates = 5)


# Close cluster
stopCluster(cl)

#------------------------------#
#---| END OF THE SCRIPT |------#
#------------------------------#
