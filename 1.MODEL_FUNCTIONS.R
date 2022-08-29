#----------------------------------------------------------#
#---| COALESCNET-BASED METACOMMUNITY MODEL & FUNCTIONS |---#
#----------------------------------------------------------#

#---| (1) COALESCENT-BASED METACOMMUNITY MODEL |---#

# This function runs the coalescent-based metacommunity model

# Arguments: (1) bs.taxa: vector containing the mean body size of each order
#            (2) xy: pond coordinates
#            (3) xy.sampled: logical indicating that only those sampled ponds should be included in the analysis
#            (4) b.disp: scaling exponent from the dispersal ability-body size scaling relationship
#            (5) b.pool: scaling exponent from the pool richness-body size scaling relationship
#            (6) b.local: scaling exponent from the local density-body size scaling relationship
#            (7) J.local: coefficient that multiplies the abundance (J) of all taxa (parameter ´a´ in the manuscript)
#            (8) m.pool: dispersal rate from the regional pool to the metacommunity
#            (9) m.neighbour: dispersal rate among communities
#            (10) taxa: character vector containing the name of the orders included in the analysis

body.size.diversity.4<-function(bs.taxa, xy, xy.sampled, b.disp,  b.pool, b.local, J.local, m.pool, m.neighbour, taxa){
  
  require(vegan)
  require(igraph)
  require(iNEXT)
  
  diversity.bs.coalescent<-NULL
  
  ### (1) CREATE A MATRIX CONTAINING THE EUCLIDEAN DISTANCE BETWEEN EACH PAIR OF PONDS ###
  
  sampled.ponds<-xy[which(xy[,3]==1),]                                                # Select those ponds that were sampled
  dist.pond<-round(as.matrix(vegdist(sampled.ponds, method = "euclidean", diag = F))) # Euclidean distance between ponds
  dist.pond<-1-(dist.pond/max(dist.pond))                                             # Standardization of the Euclidean distance between ponds
  
  
  ### (2) CREATE A VECTOR CONTAINING THE RICHNESS OF THE POOL FOR EACH TAXON ACCORDING TO ITS BODY SIZE ###
  
  gamma.taxon.pool<-(round(as.vector(100*(bs.taxa/max(bs.taxa))^b.pool)))   # Pool richness-body size scaling relationship 
  names(gamma.taxon.pool)<-taxa
  
  
  ### (3) PREPARE COALESCENT DYNAMIC AND APPLY IT TO EACH TAXON ###
  
  for (taxon.t in taxa){
    
    ## (3.0) Select the body size of taxon taxon.t 
    
    BS.taxon.t<-bs.taxa[which(rownames(bs.taxa)==taxon.t),] 
    
    ## (3.1) Create a matrix containing the probability of dispersal between ponds (off-diagonal elements) and the probability of local recruitment (diagonal elements) 
    ##       considering: (i) closeness between ponds; (ii) the effect of body size in the dispersal ability of each taxon (dispersal ability-body size scaling relationship)
    
    M.migra<-m.neighbour*(((BS.taxon.t/max(bs.taxa))^b.disp)*dist.pond)   # Matrix containing the probability of dispersal between ponds 
    M.migra.temp<-M.migra                                                 
    diag(M.migra.temp)<-0                                                 
    diag(M.migra)<-1-(apply(M.migra.temp, 1, sum))                        # Introduce the probability of local recruitment into the diagonal
    
    ## (3.2) Create a vector containing the ID of the species that belong to taxon taxon.t in the regional pool ##
    
    gamma.taxon.t.pool<-gamma.taxon.pool[taxon.t]   
    spp.id<-1:gamma.taxon.t.pool                                   
    spp.pool.id<-rep(taxon.t, gamma.taxon.t.pool)                  
    spp.pool.id<-paste0(spp.pool.id, seq(1:length(spp.pool.id)))
    
    Meta.pool<-rep(1, length(spp.pool.id))
    names(Meta.pool)<-spp.pool.id
    
    ## (3.3)  Create a vector containing the local abundance (J) of each taxon in each community according to its body size (local density-body size scaling relationship) ##  
    
    ifelse(b.local==0, J.taxon.local<-as.matrix(rep(sum(round(J.local*((bs.taxa/max(bs.taxa))^b.local))),length(taxa))),
                       J.taxon.local<-round(J.local*((bs.taxa/max(bs.taxa))^b.local)))  # Vector containing the abundance (J) of each taxon in each community given its body size 
    
    if(b.local==0) rownames(J.taxon.local)<-taxa
    
    J.taxon.t.local<-J.taxon.local[which(rownames(J.taxon.local)==taxon.t),]  # Abundance of taxon taxon.t in each community (J) according to its body size
    
    ## (3.4) Step 1 of the coalescent dynamic: community colonization ## 
    
    Meta.coalescent<-NULL # Metacommunity to be filled (spp-site matrix)
    
    for(i in 1:ncol(M.migra)){
      Meta.coalescent<-cbind(Meta.coalescent, rmultinom(1, 1, Meta.pool))   # Add one different individual randomly selected from the regional pool to each community 
    }
    
    ## (3.5) Step 2 of the coalescent dynamic: addition of individuals to each community until J is reached ##
    
    for (i in 2:J.taxon.t.local){
      cat("coalescent construction in J: ", i, "\n")  
      
      Meta.coalescent.est<-Meta.coalescent/(i-1) # Standardize the relative abundance of each species in each community 
      
      Pool.neighbor<-(Meta.coalescent.est%*%M.migra)+m.pool*Meta.pool # Spp-site matrix containing the probability of selecting an individual from any particular species to be added to a specific community
                                                                      
      new<-NULL
      
      for(k in 1:ncol(Pool.neighbor)){            
        
        new.k<-rmultinom(1,1,Pool.neighbor[,k]) # Select one individual from a particular species to be added to community k
        
        new<-cbind(new, new.k)
        
      }
      
      Meta.coalescent<-Meta.coalescent+new # Add the new recruit to each community 
      
    }
    
    Meta.coalescent<-as.data.frame(Meta.coalescent)
    
    
  ### (4) ESTIMATE DIVERSITY METRICS ###    
    
    ## (4.1) Standardized mean alpha diversity (from Meta.coalescent) ###
    
    alpha_inext<-iNEXT(Meta.coalescent, q=0, datatype="abundance", size=rep(min(J.taxon.local)*2, ncol(Meta.coalescent)))
    alpha_comm.coalescent.std<-lapply((lapply(alpha_inext$iNextEst, `[[`, 4)), unique)  # Standardized alpha diversity per community (columns in Meta.coalescent)
    alpha_comm.coalescent.std<-lapply(alpha_comm.coalescent.std, round)                                  
    alpha_mean.coalescent.std<-mean(unlist(alpha_comm.coalescent.std))                  # Standardized mean alpha diversity
    
    ## (4.1) Non-standardized mean alpha diversity (from Meta.coalescent)
    
    alpha_comm.coalescent<-apply(replace(Meta.coalescent, Meta.coalescent>0, 1), 2, sum)  # Alpha diversity per community (columns in Meta.coalescent)
    alpha_mean.coalescent<-mean(alpha_comm.coalescent)                                    # Non-standardized mean alpha diversity
    
    ## (4.2) Gamma diversity (from Meta.coalescent)
    
    occurrence.coalescent<-apply(replace(Meta.coalescent, Meta.coalescent>0, 1), 1, sum)  
    gamma.coalescent<-sum(ifelse(occurrence.coalescent>=1, 1, 0)) # Gamma diversity
    
    ## (4.3) Non-standardized additive beta diversity and standardized additive beta diversity  
    
    beta_add.coalescent<-gamma.coalescent-alpha_mean.coalescent
    beta_add.coalescent.std<-gamma.coalescent-alpha_mean.coalescent.std
    
    ## (4.4) Bray-Curtis beta diversity 
       
    beta_bray.coalescent<-vegdist(t(Meta.coalescent), method="bray", diag=F) 
    beta_mean_bray.coalescent<-round(mean(beta_bray.coalescent),2)           
    
    ## (4.5) Jaccard beta diversity
     
    beta_jaccard.coalescent<-vegdist(t(Meta.coalescent), method="jaccard", binary=T, diag=F) 
    beta_mean_jaccard.coalescent<-round(mean(beta_jaccard.coalescent),2)                    
    
    #
    
    diversity.coalescent.taxon.t<-c(alpha_mean.coalescent, alpha_mean.coalescent.std, beta_add.coalescent, beta_add.coalescent.std, beta_mean_bray.coalescent, beta_mean_jaccard.coalescent, gamma.coalescent)
    names(diversity.coalescent.taxon.t)<-c("alpha.coalescent", "alpha.coalescent.std", "beta.add.coalescent", "beta.add.coalescent.std", "beta.bray.coalescent", "beta.jaccard.coalescent", "gamma.coalescent")
    
    diversity.bs.coalescent<-rbind(diversity.bs.coalescent, diversity.coalescent.taxon.t)
    
    
  }  # End routine
  
  rownames(diversity.bs.coalescent)<-taxa
  diversity.bs.coalescent<-cbind(diversity.bs.coalescent, log2(bs.taxa))
  colnames(diversity.bs.coalescent)[8]<-"log2(BS)"
  
  return(diversity.bs.coalescent)
  
}


#---| (2) MEAN BODY SIZE |---#

# This function estimates the mean body size of each taxon

# Arguments: (1) taxa: character vector containing the names of the orders included in the analysis
#            (2) M: Database
#            (3) taxa_en: number of the column of Database that contains the names of the orders
#            (4) volumen_en: number of the column of Database that contains the body size of each individual

bz_x_taxa<-function(taxa, M, taxa_en, volumen_en){
  
  out<-NULL
  out.t<-NULL
  
  for (i in taxa) {
    
    taxa_i<-which(M[,taxa_en]==i) 
    bz_i<-M[taxa_i,volumen_en] 
    bz_mean_i<-mean(bz_i,na.rm=T) 
    out<-rbind(out, bz_mean_i)

  }
  
  out.t<-cbind(out.t, out)
  colnames(out.t)<-"Mean body size"
  rownames(out.t)<-taxa
  return(out.t)
}


#---| (3) PARAMETER SPACE |---#

# This function runs the coalescent-based metacommunity model along a gradient in b.pool, b.disp and b.local returning the parameter space  

# Arguments: (1) diversity.metrics: character vector containing the name of the diversity metrics to be included in the analysis
#            (2) disp: any of the letters ´i´, ´j´ or ´k´ indicating which numeric vector (i, j or k) will contain the values of the scaling exponent b.disp
#            (3) density: any of the letters ´i´, ´j´ or ´k´ indicating which numeric vector (i, j or k) will contain the values of the scaling exponent b.local
#            (4) pool: any of the letters ´i´, ´j´ or ´k´ indicating which numeric vector (i, j or k) will contain the values of the scaling exponent b.pool
#            (5) gradient.i: numeric vector containing the values of the scaling exponent named as ´i´
#            (6) gradient.j: numeric vector containing the values of the scaling exponent named as ´j´
#            (7) gradient.k: numeric vector containing the values of the scaling exponent named as ´k´
#            (8) N.replicates: number of times (i.e.replicates) that the model should be run for each combination of i, j and k

espacio.parametros<-function(diversity.metrics, disp, density, pool, gradient.i, gradient.j, gradient.k, N.replicates){
  
  library(lme4)
  library(rlist)
  library(rsq)
  library(parallel)
  
  out.alpha<-list()
  out.beta.add<-list()
  out.gamma<-list()
  
  cl<-makeCluster(N.replicates)
  print(paste("Starting time:", Sys.time()))
  
  # (1) Starts for cycles
  
  for (i in gradient.i) { 
    if(disp=="i")cat("--------------------------------------------> DISPERSION: ", i, "\n")
    if(density=="i")cat("--------------------------------------------> DENSIDAD: ", i, "\n")
    if(pool=="i")cat("--------------------------------------------> POOL: ", i, "\n")
    ii<-which(gradient.i==i)
    
    M.kj.alpha<-matrix(NA, nrow=length(gradient.k), ncol=length(gradient.j)) 
    M.kj.beta.add<-matrix(NA, nrow=length(gradient.k), ncol=length(gradient.j))
    M.kj.gamma<-matrix(NA, nrow=length(gradient.k), ncol=length(gradient.j))
    
    for(j in gradient.j) {
      id.j<-which(gradient.j==j)                                             
      if(disp=="j") cat("-----------------------> Dispersion: ", j, "\n")
      if(density=="j") cat("-----------------------> Densidad: ", j, "\n")
      if(pool=="j") cat("-----------------------> Pool: ", j, "\n")
      
      for (k in gradient.k) { 
        if(disp=="k")cat("---> Dispersion: ", k, "\n")
        if(density=="k")cat("---> Densidad: ", k, "\n")
        if(pool=="k")cat("---> Pool: ", k, "\n")
        id.k<-which(gradient.k==k)                                           
        
        if(disp=="i" & density=="j" & pool=="k"){list.replicates.ijk<-clusterCall(cl, fun=body.size.diversity.4, bs.taxa=bs.orders, xy=Pond_xy, xy.sampled = T, 
                                                                                  b.disp=i, b.local=j, b.pool=k, J.local=50, 
                                                                                  m.pool=0.00001, m.neighbour=0.01, taxa=orders) # Gives a list object with the replicates of the result of the model.
        save.image()
        
        df.replicates.ijk<-list.rbind(list.replicates.ijk)                                    # Transforms list.replicates in a data.frame.
        replicate.id<-rep(1:length(list.replicates.ijk), each=nrow(list.replicates.ijk[[1]])) # Creates a vector that identifies each replicate in the dataframe df.replicates.
        replicates.ijk<-as.data.frame(cbind(df.replicates.ijk, replicate.id)) 
        colnames(replicates.ijk)[8]<-"BS"
        
        }
        
        if(disp=="i" & density=="k" & pool=="j"){list.replicates.ijk<-clusterCall(cl, fun=body.size.diversity.4, bs.taxa=bs.orders, xy=Pond_xy, xy.sampled = T, 
                                                                                  b.disp=i, b.local=k, b.pool=j, J.local=50, 
                                                                                  m.pool=0.00001, m.neighbour=0.01, taxa=orders) # Gives a list object with the replicates of the result of the model.
        save.image()
        
        df.replicates.ijk<-list.rbind(list.replicates.ijk)                                    # Transforms list.replicates in a data.frame.
        replicate.id<-rep(1:length(list.replicates.ijk), each=nrow(list.replicates.ijk[[1]])) # Creates a vector that identifies each replicate in the dataframe df.replicates.
        replicates.ijk<-as.data.frame(cbind(df.replicates.ijk, replicate.id)) 
        colnames(replicates.ijk)[8]<-"BS"
        }
        
        if(disp=="j" & density=="i" & pool=="k"){list.replicates.ijk<-clusterCall(cl, fun=body.size.diversity.4, bs.taxa=bs.orders, xy=Pond_xy, xy.sampled = T, 
                                                                                  b.disp=j, b.local=i, b.pool=k, J.local=50, 
                                                                                  m.pool=0.00001, m.neighbour=0.01, taxa=orders) # Gives a list object with the replicates of the result of the model.
        save.image()
        
        df.replicates.ijk<-list.rbind(list.replicates.ijk)                                    # Transforms list.replicates in a data.frame.
        replicate.id<-rep(1:length(list.replicates.ijk), each=nrow(list.replicates.ijk[[1]])) # Creates a vector that identifies each replicate in the dataframe df.replicates.
        replicates.ijk<-as.data.frame(cbind(df.replicates.ijk, replicate.id)) 
        colnames(replicates.ijk)[8]<-"BS"
        }
        
        if(disp=="j" & density=="k" & pool=="i"){list.replicates.ijk<-clusterCall(cl, fun=body.size.diversity.4, bs.taxa=bs.orders, xy=Pond_xy, xy.sampled = T, 
                                                                                  b.disp=j, b.local=k, b.pool=i, J.local=50, 
                                                                                  m.pool=0.00001, m.neighbour=0.01, taxa=orders) # Gives a list object with the replicates of the result of the model.
        save.image()
        
        df.replicates.ijk<-list.rbind(list.replicates.ijk)                                    # Transforms list.replicates in a data.frame.
        replicate.id<-rep(1:length(list.replicates.ijk), each=nrow(list.replicates.ijk[[1]])) # Creates a vector that identifies each replicate in the dataframe df.replicates.
        replicates.ijk<-as.data.frame(cbind(df.replicates.ijk, replicate.id)) 
        colnames(replicates.ijk)[8]<-"BS"
        }
        
        if(disp=="k" & density=="i" & pool=="j"){list.replicates.ijk<-clusterCall(cl, fun=body.size.diversity.4, bs.taxa=bs.orders, xy=Pond_xy, xy.sampled = T, 
                                                                                  b.disp=k, b.local=i, b.pool=j, J.local=50, 
                                                                                  m.pool=0.00001, m.neighbour=0.01, taxa=orders) # Gives a list object with the replicates of the result of the model.
        save.image()
        
        df.replicates.ijk<-list.rbind(list.replicates.ijk)                                    # Transforms list.replicates in a data.frame.
        replicate.id<-rep(1:length(list.replicates.ijk), each=nrow(list.replicates.ijk[[1]])) # Creates a vector that identifies each replicate in the dataframe df.replicates.
        replicates.ijk<-as.data.frame(cbind(df.replicates.ijk, replicate.id)) 
        colnames(replicates.ijk)[8]<-"BS"
        }
        
        if(disp=="k" & density=="j" & pool=="i"){list.replicates.ijk<-clusterCall(cl, fun=body.size.diversity.4, bs.taxa=bs.orders, xy=Pond_xy, xy.sampled = T, 
                                                                                  b.disp=k, b.local=j, b.pool=i, J.local=50, 
                                                                                  m.pool=0.00001, m.neighbour=0.01, taxa=orders) # Gives a list object with the replicates of the result of the model.
        save.image()
        
        df.replicates.ijk<-list.rbind(list.replicates.ijk)                                    # Transforms list.replicates in a data.frame.
        replicate.id<-rep(1:length(list.replicates.ijk), each=nrow(list.replicates.ijk[[1]])) # Creates a vector that identifies each replicate in the dataframe df.replicates.
        replicates.ijk<-as.data.frame(cbind(df.replicates.ijk, replicate.id)) 
        colnames(replicates.ijk)[8]<-"BS"
        }
        
        # (2) GLM´s between diversity metrics and taxa body size 
        
        for(d in diversity.metrics) {
    
          if (d=="beta.add.coalescent.std") {
            
            if (length(which(replicates.ijk[,d]==0))>0) { M.kj.beta.add[id.k, id.j]<-NA } 
            
            else {  beta.add.glm.replicates<-lmList(log(beta.add.coalescent.std) ~ BS | replicate.id,  replicates.ijk, family = gaussian) # Applies a glm model between beta.add.coalescent.std and body size for each replicate 
            beta.add.coef<-coef(beta.add.glm.replicates)            # Gives coefficients of the linear models.
            beta.add.slopes<-beta.add.coef[2]                       # Gives the slopes of the linear models.
            mean.beta.add.slopes<-mean(as.matrix(beta.add.slopes))
            names(mean.beta.add.slopes)<-"Beta.add.mean.slope"
            M.kj.beta.add[id.k, id.j]<-mean.beta.add.slopes }
          } 
          
          if (d=="alpha.coalescent.std") {
            
            alpha.glm.replicates<-lmList(log(alpha.coalescent.std) ~ BS | replicate.id,  replicates.ijk, family = gaussian) # Applies a glm model between alpha.coalescent.std and body size for each replicate
            alpha.coef<-coef(alpha.glm.replicates) # Gives coefficients of the linear models.
            alpha.slopes<-alpha.coef[2]            # Gives the slopes of the linear models.
            mean.alpha.slopes<-mean(as.matrix(alpha.slopes))
            names(mean.alpha.slopes)<-"Alpha.mean.slopes"
            M.kj.alpha[id.k, id.j]<-mean.alpha.slopes }
          
          if (d=="gamma.coalescent") {
            
            gamma.glm.replicates<-lmList(log(gamma.coalescent) ~ BS | replicate.id,  replicates.ijk, family = gaussian) # Applies a glm model between gamma.coalescent.std and body size for each replicate
            gamma.coef<-coef(gamma.glm.replicates) # Gives coefficients of the linear models.
            gamma.slopes<-gamma.coef[2]            # Gives the slopes of the linear models.
            mean.gamma.slopes<-mean(as.matrix(gamma.slopes))
            names(mean.gamma.slopes)<-"Gamma.mean.slopes"
            M.kj.gamma[id.k, id.j]<-mean.gamma.slopes }
          
        }
        
        save.image()
        
      } # Closes cycle for k
      save.image()
      
    } # Closes cycle for j
    save.image()
    
    
    colnames(M.kj.alpha)<-gradient.j
    rownames(M.kj.alpha)<-gradient.k
    
    colnames(M.kj.beta.add)<-gradient.j
    rownames(M.kj.beta.add)<-gradient.k
    
    colnames(M.kj.gamma)<-gradient.j
    rownames(M.kj.gamma)<-gradient.k
    
    out.alpha[[ii]]<-M.kj.alpha
    out.beta.add[[ii]]<-M.kj.beta.add
    out.gamma[[ii]]<-M.kj.gamma
    
    save.image()
    
  }
  
  names(out.alpha)<-gradient.i
  names(out.beta.add)<-gradient.i
  names(out.gamma)<-gradient.i
  
  
  out<-list("Mean.alpha.slope"=out.alpha, "Mean.beta.add.slope"=out.beta.add, "Mean.gamma.slope"=out.gamma)
  stopCluster(cl)
  print(paste("Ending time:", Sys.time()))
  return(out)
  
}


#---| (4) STANDARDIZED DIVERSITY IN THE METACOMMUNITY OF TEMPORARY PONDS |---#

# This function estimates the standardized diversity in the metacommunity of temporary ponds 

# Arguments: (1) M: Database
#            (2) comm.col: number of the column of Database that contains the ID of the ponds
#            (3) order.col: number of the column of Database that contains the names of the orders
#            (4) species.col: number of the column of Database that contains the names of the species
#            (5) taxa.names: character vector containing the names of the orders included in the analysis

diversity.std.orders<-function(M, comm.col, order.col, species.col, taxa.names){
  
  library(iNEXT)
  
  ### (1) GAMMA DIVERSITY ###
  
  N.obs.orders.spp<-list()
  Gamma.std.orders<-NULL
  Gamma.obs<-NULL 
  
  for(taxon in taxa.names){
    
    taxon.id<-which(taxa.names==taxon)
    M.taxon<-M[which(M[,order.col]==taxa.names[taxon.id]), ] 
    taxon.spp.names<-unique(M.taxon[,species.col])
    S.taxon<-length(taxon.spp.names) 
    
    Gamma.obs<-c(Gamma.obs, S.taxon) 
    
    N.taxon.spp<-NULL
    
    for(taxon.sp.i in taxon.spp.names){
      
      M.taxon.sp.i<-M.taxon[which(M.taxon[,species.col]==taxon.sp.i), ]
      N.taxon.sp.i<-length(M.taxon.sp.i[, species.col])
      N.taxon.spp<-c(N.taxon.spp, N.taxon.sp.i)
      
    }
    
    N.obs.orders.spp[[taxon.id]]<-N.taxon.spp
    
  }
  
  names(Gamma.obs)<-taxa.names 
  names(N.obs.orders.spp)<-taxa.names
  min.abundance<-min(as.numeric(lapply(N.obs.orders.spp,sum)))
  std.N.obs.orders.spp<-iNEXT(N.obs.orders.spp, q=0, datatype = "abundance", size = c(1, min.abundance*2))
  
  for(taxon in taxa.names) {
    
    Gamma.std.taxon<-std.N.obs.orders.spp$iNextEst[[taxon]]["2", "qD"]
    Gamma.std.orders<-c(Gamma.std.orders, Gamma.std.taxon)
    
  }
  
  Gamma.std.orders<-round(Gamma.std.orders)
  names(Gamma.std.orders)<-taxa.names

  ### (2) MEAN ALPHA DIVERSITY ###
  
  communities<-sort(unique(M[, comm.col]))
  
  S.std.mean.orders<-NULL
  
  for(taxon in taxa.names){
    
    S.std.orders<-NULL
    
    taxon.id<-which(taxa.names==taxon)
    M.taxon<-M[which(M[,order.col]==taxa.names[taxon.id]), ] 
    spp.taxon<-unique(M.taxon[,species.col])
    N.taxon.spp.communities<-NULL
    
    for(sp_i in spp.taxon){
      
      M.taxon.sp_i<-M.taxon[which(M.taxon[,species.col]==sp_i), ] 
      N.sp_i.communities<-NULL
      
      for(comm_i in communities){
        
        M.taxon.sp_i.comm_i<-M.taxon.sp_i[which(M.taxon.sp_i[, comm.col]==comm_i),] 
        N.taxon.sp_i.comm_i<-length(M.taxon.sp_i.comm_i[,species.col]) 
        N.sp_i.communities<-c(N.sp_i.communities, N.taxon.sp_i.comm_i) 
      }
      
      N.taxon.spp.communities<-rbind(N.taxon.spp.communities, N.sp_i.communities) 
      
    }
    
    colnames(N.taxon.spp.communities)<-communities
    rownames(N.taxon.spp.communities)<-spp.taxon
    
    empty.comm<-which(apply(N.taxon.spp.communities, 2, sum)==0)

    if(length(empty.comm)>0) N.taxon.spp.communities<-N.taxon.spp.communities[, -empty.comm]

    std.S.taxon<-iNEXT(N.taxon.spp.communities, q = 0, datatype = "abundance", size = c(1, min.abundance*2)) 

    for(comm in colnames(N.taxon.spp.communities)) {
      
      S.std.taxon<-std.S.taxon$iNextEst[[as.character(comm)]][nrow(std.S.taxon$iNextEst[[as.character(comm)]]), "qD"]
      S.std.orders<-c(S.std.orders, S.std.taxon)
      
    }

    alpha.mean.taxon<-(mean(S.std.orders))
    S.std.mean.orders<-c(S.std.mean.orders, alpha.mean.taxon)
  }
  
  ### (3) BETA DIVERSITY ###
  
  #beta.add<-Gamma.std.orders-S.std.mean.orders
  #beta.multi<-Gamma.std.orders/S.std.mean.orders
  
  beta.add<-Gamma.obs-S.std.mean.orders 
  beta.multi<-Gamma.obs/S.std.mean.orders 
  
  return(list("alpha.mean"=S.std.mean.orders, "beta.add"=beta.add, "gamma"=Gamma.obs, "beta.multi"=beta.multi))
  
}


#---| (5) RICHNESS OF THE REGIONAL POOL |---#

# This function estimates the richness of the regional pool for each order in the metacommunity of temporary ponds

# Arguments: (1) M: Database
#            (2) spp.in: number of the column of Database that contains the names of the species
#            (3) taxa.in: number of the column of Database that contains the names of the orders
#            (4) taxa: character vector containing the names of the orders included in the analysis

pool.richness.estimate<-function(M, spp.in, taxa.in, taxa) {
  
  require(iNEXT)
  require(rlist)
  
  S.pool.taxa<-NULL
  
  for(t in taxa){
    
    M.t<-M[which(M[,taxa.in]==t),] #Sub-matrix for taxon t
    spp.t<-unique(M.t[,spp.in])
    
    abund.spp.t<-NULL
    
    for(i in spp.t){
      
      M.t.i<-M.t[which(M.t[,spp.in]==i),] #Sub-matrix for sp i of taxon t
      abund.i.t<-length(M.t.i[,spp.in])
      
      abund.spp.t<-c(abund.spp.t, abund.i.t)
      
    }
    
    Chao.t<-ChaoRichness(x=abund.spp.t, datatype="abundance")
    S.t.pool<-Chao.t[2]
    S.pool.taxa<-c(S.pool.taxa, S.t.pool)
  }
  
  S.pool.taxa<-list.rbind(S.pool.taxa)
  colnames(S.pool.taxa)<-"S.pool"
  rownames(S.pool.taxa)<-taxa
  return(S.pool.taxa)
  
}


#------------------------------#
#---| END OF THE SCRIPT |------#
#------------------------------#


