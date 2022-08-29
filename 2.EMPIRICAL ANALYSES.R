#----------------------------#
#---| EMPIRICAL ANALYSES |---#
#----------------------------#

library(rsq)
library(MASS)

#---| (1) ESTIMATION OF STANDARDIZED DIVERSITY IN THE METACOMMUNITY (Figure 3 a-c and Supporting information Table SM3) |---#

diversity.std<-diversity.std.orders(M = Dataset, comm.col = 1, order.col = 5, species.col = 3, taxa.names = orders)

summary(glm(log(diversity.std[["alpha.mean"]])~log2(bs.orders))) 
rsq(glm(log(diversity.std[["alpha.mean"]])~log2(bs.orders))) 
confint(glm(log(diversity.std[["alpha.mean"]])~log2(bs.orders)))

summary(glm(log(diversity.std[["beta.add"]])~log2(bs.orders)))   
rsq(glm(log(diversity.std[["beta.add"]])~log2(bs.orders)))  
confint(glm(log(diversity.std[["beta.add"]])~log2(bs.orders)))

summary(glm(log(diversity.std[["gamma"]])~log2(bs.orders)))      
rsq(glm(log(diversity.std[["gamma"]])~log2(bs.orders)))
confint(glm(log(diversity.std[["gamma"]])~log2(bs.orders)))


#---| (2) ESTIMATION OF THE RICHNESS OF THE REGIONAL POOL IN THE METACOMMUNITY  (Supporting information Figure SM4) |---#

gamma.pool.orders2<-pool.richness.estimate(M=Dataset, taxa=orders, spp.in = 3, taxa.in = 5)

summary(glm(log(gamma.pool.orders2)~log2(bs.orders)))
rsq(glm(log(gamma.pool.orders2)~log2(bs.orders))) 
confint(glm(log(gamma.pool.orders2)~log2(bs.orders)))


#------------------------------#
#---| END OF THE SCRIPT |------#
#------------------------------#
