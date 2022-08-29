#---------------#
#---| PLOTS |---#
#---------------#

library(ggplot2)
library(viridis)
library(patchwork)
library(latticeExtra)
library(sp)
library(lattice)
library(grid)

#---| (1) THEORETICAL DIVERSITY-BODY SIZE RELTIONSHIPS INCLUDING ALL POSSIBLE COMBINATIONS OF THE THREE SCALING RELATIONSHIPS (Figure 2) |---#

just.for.theoretical.plot<-as.data.frame(rbind(simu.all4.N50.2_0.75_0.4, 
                                               simu.none4.N50.2_0.75_0.4, 
                                               simu.disp4.N50.2_0.75_0.4, 
                                               simu.pool4.N50.2_0.75_0.4, 
                                               simu.density4.N50.2_0.75_0.4, 
                                               simu.disp.pool4.N50.2_0.75_0.4, 
                                               simu.disp.density4.N50.2_0.75_0.4, 
                                               simu.pool.density4.N50.2_0.75_0.4))
simu.id<-rep(c("Disp + Pool + Dens", "None", "Dispersal", "Pool richness", "Density", "Disp + Pool", "Disp + Dens", "Pool + Dens"), each=nrow(just.for.theoretical.plot)/length(just.for.theoretical.plot))
just.for.theoretical.plot<-as.data.frame(cbind(just.for.theoretical.plot, simu.id))
just.for.theoretical.plot$simu.id <- factor(just.for.theoretical.plot$simu.id, levels = c("Disp + Pool + Dens", "None", "Dispersal", "Pool richness", "Density", "Disp + Pool", "Disp + Dens", "Pool + Dens"))
colnames(just.for.theoretical.plot)[8]<-"BS"


# (1.1) SCALING RELATIONSHIPS IN ISOLATION #

# Alpha diversity #
plot.theoretical.alpha.isolated<-ggplot(just.for.theoretical.plot[c(which(just.for.theoretical.plot[,"simu.id"]=="Dispersal"),  
                                                                    which(just.for.theoretical.plot[,"simu.id"]=="Pool richness"), 
                                                                    which(just.for.theoretical.plot[,"simu.id"]=="Density")),],
                                        aes(x = BS, 
                                            y = log(alpha.coalescent.std), 
                                            group = as.factor(simu.id), 
                                            col = as.factor(simu.id)))  +
                                        theme_classic() +
                                        xlab("Log(Body Size)") +
                                        ylab("Log(Alpha)") +
                                        coord_cartesian(ylim = c(1.4, 3.5)) +
                                        geom_smooth(method = lm, se = T) +
                                        scale_color_discrete("Scaling relationships") +
                                        theme(legend.position = c(0.75,0.9)) +
                                        theme(axis.text=element_text(size=17),
                                              axis.title=element_text(size=17,face="bold")) +
                                        theme(legend.title = element_text(size = 13, face="bold"),
                                              legend.text = element_text(size = 13)) +
                                        theme(axis.title.x = element_text(vjust=-1), 
                                              axis.title.y = element_text(vjust=4)) +
                                        theme(plot.margin = unit(c(1,1,1,1), "cm"))


# Beta diversity #
plot.theoretical.beta.isolated<-ggplot(just.for.theoretical.plot[c(which(just.for.theoretical.plot[,"simu.id"]=="Dispersal"),  
                                                                   which(just.for.theoretical.plot[,"simu.id"]=="Pool"), 
                                                                   which(just.for.theoretical.plot[,"simu.id"]=="Density")),],
                                       aes(x = BS, 
                                           y = log(beta.add.coalescent.std), 
                                           group = as.factor(simu.id), 
                                           col = as.factor(simu.id)))  +
                                      theme_classic() +
                                      xlab("Log(Body Size)") +
                                      ylab("Log(Beta)") +
                                      geom_smooth(method = lm, se = T) +
                                      theme(legend.position = "none") +
                                      theme(axis.text=element_text(size=17),
                                            axis.title=element_text(size=17,face="bold")) +
                                      theme(axis.title.x = element_text(vjust=-1), 
                                            axis.title.y = element_text(vjust=4)) +
                                      theme(plot.margin = unit(c(1,1,1,1), "cm"))

# Gamma diversity #
plot.theoretical.gamma.isolated<-ggplot(just.for.theoretical.plot[c(which(just.for.theoretical.plot[,"simu.id"]=="Dispersal"),  
                                                                    which(just.for.theoretical.plot[,"simu.id"]=="Pool"), 
                                                                    which(just.for.theoretical.plot[,"simu.id"]=="Density")),],
                                        aes(x = BS, 
                                            y = log(gamma.coalescent), 
                                            group = as.factor(simu.id), 
                                            col = as.factor(simu.id)))  +
                                        theme_classic() +
                                        xlab("Log(Body Size)") +
                                        ylab("Log(Gamma)") +
                                        geom_smooth(method = lm, se = T) +
                                        theme(legend.position = "none") +
                                        theme(axis.text=element_text(size=17),
                                              axis.title=element_text(size=17,face="bold")) +
                                        theme(axis.title.x = element_text(vjust=-1), 
                                              axis.title.y = element_text(vjust=4)) +
                                        theme(plot.margin = unit(c(1,1,1,1), "cm"))



# (1.2) SCALING RELATIONSHIPS IN COMBINATION #

# Alpha diversity #
plot.theoretical.alpha.combinations<-ggplot(just.for.theoretical.plot[c(which(just.for.theoretical.plot[,"simu.id"]=="Disp + Pool + Dens"),  
                                                                        which(just.for.theoretical.plot[,"simu.id"]=="Disp + Pool"), 
                                                                        which(just.for.theoretical.plot[,"simu.id"]=="Disp + Dens"),
                                                                        which(just.for.theoretical.plot[,"simu.id"]=="Pool + Dens")),],
                                            aes(x = BS, 
                                                y = log(alpha.coalescent.std), 
                                                group = as.factor(simu.id), 
                                                col = as.factor(simu.id)))  +
                                            theme_classic() +
                                            xlab("Log(Body Size)") +
                                            ylab("Log(Alpha)") +
                                            coord_cartesian(ylim = c(1.4, 3.5)) +
                                            geom_smooth(method = lm, se = T) +
                                            scale_color_discrete("Scaling relationships") +
                                            theme(legend.position = c(0.8,0.9)) +
                                            theme(axis.text=element_text(size=17),
                                                  axis.title=element_text(size=17,face="bold")) +
                                            theme(legend.title = element_text(size = 13, face="bold"),
                                                  legend.text = element_text(size = 13)) +
                                            theme(axis.title.x = element_text(vjust=-1), 
                                                  axis.title.y = element_text(vjust=4)) +
                                            theme(plot.margin = unit(c(1,1,1,1), "cm"))


# Beta diversity #
plot.theoretical.beta.combinations<-ggplot(just.for.theoretical.plot[c(which(just.for.theoretical.plot[,"simu.id"]=="Disp + Pool + Dens"),  
                                                                       which(just.for.theoretical.plot[,"simu.id"]=="Disp + Pool"), 
                                                                       which(just.for.theoretical.plot[,"simu.id"]=="Disp + Dens"),
                                                                       which(just.for.theoretical.plot[,"simu.id"]=="Pool + Dens")),],
                                           aes(x = BS, 
                                               y = log(beta.add.coalescent.std), 
                                               group = as.factor(simu.id), 
                                               col = as.factor(simu.id)))  +
                                          theme_classic() +
                                          xlab("Log(Body Size)") +
                                          ylab("Log(Beta)") +
                                          geom_smooth(method = lm, se = T) +
                                          scale_color_discrete("Scaling relationship") +
                                          theme(legend.position = "none") +
                                          theme(axis.text=element_text(size=17),
                                                axis.title=element_text(size=17,face="bold")) +
                                          theme(axis.title.x = element_text(vjust=-1), 
                                                axis.title.y = element_text(vjust=4)) +
                                          theme(plot.margin = unit(c(1,1,1,1), "cm"))

# Gamma diversity #
plot.theoretical.gamma.combinations<-ggplot(just.for.theoretical.plot[c(which(just.for.theoretical.plot[,"simu.id"]=="Disp + Pool + Dens"),  
                                                                        which(just.for.theoretical.plot[,"simu.id"]=="Disp + Pool"), 
                                                                        which(just.for.theoretical.plot[,"simu.id"]=="Disp + Dens"),
                                                                        which(just.for.theoretical.plot[,"simu.id"]=="Pool + Dens")),],
                                            aes(x = BS, 
                                                y = log(gamma.coalescent), 
                                                group = as.factor(simu.id), 
                                                col = as.factor(simu.id)))  +
                                            theme_classic() +
                                            xlab("Log(Body Size)") +
                                            ylab("Log(Gamma)") +
                                            geom_smooth(method = lm, se = T) +
                                            scale_color_discrete("Scaling relationship") +
                                            theme(legend.position = "none") +
                                            theme(axis.text=element_text(size=17),
                                                  axis.title=element_text(size=17,face="bold")) +
                                            theme(axis.title.x = element_text(vjust=-1), 
                                                  axis.title.y = element_text(vjust=4)) +
                                            theme(plot.margin = unit(c(1,1,1,1), "cm"))


# (1.3) NO SCALING RELATIONSHIPS #

# Alpha diversity #
plot.theoretical.alpha.no.scaling<-ggplot(just.for.theoretical.plot[which(just.for.theoretical.plot[,"simu.id"]=="None"),],
                                          aes(x = BS, 
                                              y = log(alpha.coalescent.std), 
                                              group = as.factor(simu.id), 
                                              col = as.factor(simu.id)))  +
                                          theme_classic() +
                                          xlab("Log(Body Size)") +
                                          ylab("Log(Alpha)") +
                                          coord_cartesian(ylim = c(2.3, 3.3)) +
                                          geom_smooth(method = lm, se = T) +
                                          scale_color_discrete("Scaling relationships") +
                                          theme(legend.position = c(0.8,0.93)) +
                                          theme(axis.text=element_text(size=17),
                                                axis.title=element_text(size=17,face="bold")) +
                                          theme(legend.title = element_text(size = 13, face="bold"),
                                                legend.text = element_text(size = 13)) +
                                          theme(axis.title.x = element_text(vjust=-1), 
                                                axis.title.y = element_text(vjust=4)) +
                                          theme(plot.margin = unit(c(1,1,1,1), "cm"))

# Beta diversity #
plot.theoretical.beta.no.scaling<-ggplot(just.for.theoretical.plot[which(just.for.theoretical.plot[,"simu.id"]=="None"),],
                                         aes(x = BS, 
                                             y = log(beta.add.coalescent.std), 
                                             group = as.factor(simu.id), 
                                             col = as.factor(simu.id)))  +
                                        theme_classic() +
                                        xlab("Log(Body Size)") +
                                        ylab("Log(Beta)") +
                                        coord_cartesian(ylim = c(1, 3)) +
                                        geom_smooth(method = lm, se = T) +
                                        theme(legend.position = "none") +
                                        theme(axis.text=element_text(size=17),
                                              axis.title=element_text(size=17,face="bold")) +
                                        theme(axis.title.x = element_text(vjust=-1), 
                                              axis.title.y = element_text(vjust=4)) +
                                        theme(plot.margin = unit(c(1,1,1,1), "cm"))

# Gamma diversity #
plot.theoretical.gamma.no.scaling<-ggplot(just.for.theoretical.plot[which(just.for.theoretical.plot[,"simu.id"]=="None"),],
                                          aes(x = BS, 
                                              y = log(gamma.coalescent), 
                                              group = as.factor(simu.id), 
                                              col = as.factor(simu.id)))  +
                                          theme_classic() +
                                          xlab("Log(Body Size)") +
                                          ylab("Log(Gamma)") +
                                          coord_cartesian(ylim = c(2.5, 4)) +
                                          geom_smooth(method = lm, se = T) +
                                          theme(legend.position = "none") +
                                          theme(axis.text=element_text(size=17),
                                                axis.title=element_text(size=17,face="bold")) +
                                          theme(axis.title.x = element_text(vjust=-1), 
                                                axis.title.y = element_text(vjust=4)) +
                                          theme(plot.margin = unit(c(1,1,1,1), "cm"))


# (1.4) ALL PLOTS TOGETHER #

plot.theoretical.alpha.isolated + plot.theoretical.beta.isolated + plot.theoretical.gamma.isolated +
  plot.theoretical.alpha.combinations + plot.theoretical.beta.combinations + plot.theoretical.gamma.combinations +
  plot.theoretical.alpha.no.scaling + plot.theoretical.beta.no.scaling + plot.theoretical.gamma.no.scaling +
  plot_layout(ncol = 3, nrow = 3)



#---| (2) PARAMETER SPACE (Supporting Information Figure SM3) |---#

# Alpha diversity #
par.space.ALPHA_0<-par.space.1$Mean.alpha.slope$`0`
par.space.ALPHA_0.1<-par.space.1$Mean.alpha.slope$`-0.1`
par.space.ALPHA_0.2<-par.space.1$Mean.alpha.slope$`-0.2`
par.space.ALPHA_0.3<-par.space.1$Mean.alpha.slope$`-0.3`
par.space.ALPHA_0.4<-par.space.2$Mean.alpha.slope$`-0.4`
par.space.ALPHA_0.6<-par.space.2$Mean.alpha.slope$`-0.6`

Plot.par.space.ALPHA_0<-levelplot(value ~ Var1 * Var2, reshape2::melt(par.space.ALPHA_0), panel = panel.levelplot.points, cex = 0, contour=T, 
                                  main=list("Scaling in pool richnnes (b.pool = 0)", cex=1.5), xlab = list("Scaling in dispersal (b.disp)", cex=1.5), ylab = list("Scaling in local density (b.local)", cex=1.5),
                                  col.regions = viridis(40), cuts=40, ylim = c(-2.05, 0.05), xlim = c(-0.05, 2.05), 
                                  at=seq(round(min(par.space.ALPHA_0),4)-0.05, max(par.space.ALPHA_0), length.out=40),
                                  scales=list(x=list(cex=1.5), y=list(cex=1.5)),
                                  colorkey = list(labels=list(cex=1.2))) + 
                                  layer_(panel.2dsmoother(..., n = 600))

Plot.par.space.ALPHA_0.1<-levelplot(value ~ Var1 * Var2, reshape2::melt(par.space.ALPHA_0.1), panel = panel.levelplot.points, cex = 0, contour=T, 
                                    main="Scaling in pool richnnes (b.pool = -0.1)", xlab = "Scaling in dispersal (b.disp)", ylab = "Scaling in local density (b.local)",
                                    col.regions = viridis(40), cuts=40, ylim = c(-2.05, 0.05), xlim = c(-0.05, 2.05), 
                                    at=seq(round(min(par.space.ALPHA_0.1),4)-0.05, max(par.space.ALPHA_0.1), length.out=40)) + 
                                    layer_(panel.2dsmoother(..., n = 600))

Plot.par.space.ALPHA_0.2<-levelplot(value ~ Var1 * Var2, reshape2::melt(par.space.ALPHA_0.2), panel = panel.levelplot.points, cex = 0, contour=T, 
                                    main=list("Scaling in pool richnnes (b.pool = -0.2)", cex=1.5), xlab = list("Scaling in dispersal (b.disp)", cex=1.5), ylab = list("Scaling in local density (b.local)", cex=1.5),
                                    col.regions = viridis(40), cuts=40, ylim = c(-2.05, 0.05), xlim = c(-0.05, 2.05), 
                                    at=seq(round(min(par.space.ALPHA_0.2),4)-0.05, max(par.space.ALPHA_0.2), length.out=40),
                                    scales=list(x=list(cex=1.5), y=list(cex=1.5)),
                                    colorkey = list(labels=list(cex=1.2)),
                                    labels=list(col="black", cex=0.8)) + 
                                    layer_(panel.2dsmoother(..., n = 600))

Plot.par.space.ALPHA_0.3<-levelplot(value ~ Var1 * Var2, reshape2::melt(par.space.ALPHA_0.3), panel = panel.levelplot.points, cex = 0, contour=T, 
                                    main=list("Scaling in pool richnnes (b.pool = -0.3)", cex=1.5), xlab = list("Scaling in dispersal (b.disp)", cex=1.5), ylab = list("Scaling in local density (b.local)", cex=1.5),
                                    col.regions = viridis(40), cuts=40, ylim = c(-2.05, 0.05), xlim = c(-0.05, 2.05), 
                                    at=seq(round(min(par.space.ALPHA_0.3),4)-0.05, max(par.space.ALPHA_0.3), length.out=40),
                                    # labels=list("labels"=c(rep("", 7), "-0.19", rep("", 31)), col="white", cex=0.8),
                                    scales=list(x=list(cex=1.5), y=list(cex=1.5)),
                                    colorkey = list(labels=list(cex=1.2))) + 
                                    layer_(panel.2dsmoother(..., n = 200))

Plot.par.space.ALPHA_0.4<-levelplot(value ~ Var1 * Var2, reshape2::melt(par.space.ALPHA_0.4), panel = panel.levelplot.points, cex = 0, contour=T, 
                                    main=list("Scaling in pool richnnes (b.pool = -0.4)", cex=1.5), xlab = list("Scaling in dispersal (b.disp)", cex=1.5), ylab = list("Scaling in local density (b.local)", cex=1.5),
                                    col.regions = viridis(40), cuts=40, ylim = c(-2.05, 0.05), xlim = c(-0.05, 2.05), 
                                    at=seq(round(min(par.space.ALPHA_0.4),4)-0.05, max(par.space.ALPHA_0.4), length.out=40),
                                   # labels=list("labels"=c(rep("", 7), "-0.19", rep("", 31)), col="white", cex=0.8),
                                    scales=list(x=list(cex=1.5), y=list(cex=1.5)),
                                    colorkey = list(labels=list(cex=1.2))) + 
                                    layer_(panel.2dsmoother(..., n = 200))

Plot.par.space.ALPHA_0.6<-levelplot(value ~ Var1 * Var2, reshape2::melt(par.space.ALPHA_0.6), panel = panel.levelplot.points, cex = 0, contour=T, 
                                    main=list("Scaling in pool richnnes (b.pool = -0.6)", cex=1.5), xlab = list("Scaling in dispersal (b.disp)", cex=1.5), ylab = list("Scaling in local density (b.local)", cex=1.5),
                                    col.regions = viridis(40), cuts=40, ylim = c(-2.05, 0.05), xlim = c(-0.05, 2.05), 
                                    at=seq(round(min(par.space.ALPHA_0.6),4)-0.05, max(par.space.ALPHA_0.6), length.out=40),
                                    scales=list(x=list(cex=1.5), y=list(cex=1.5)),
                                    colorkey = list(labels=list(cex=1.2))) + 
                                    layer_(panel.2dsmoother(..., n = 200)) 

# Beta diversity #
par.space.BETA_0<-par.space.1$Mean.beta.add.slope$`0`
par.space.BETA_0.1<-par.space.1$Mean.beta.add.slope$`-0.1`
par.space.BETA_0.2<-par.space.1$Mean.beta.add.slope$`-0.2`
par.space.BETA_0.3<-par.space.1$Mean.beta.add.slope$`-0.3`
par.space.BETA_0.4<-par.space.2$Mean.beta.add.slope$`-0.4`
par.space.BETA_0.6<-par.space.2$Mean.beta.add.slope$`-0.6`

Plot.par.space.BETA_0<-levelplot(value ~ Var1 * Var2, reshape2::melt(par.space.BETA_0), panel = panel.levelplot.points, cex = 0, contour=T, 
                                 main=list("Scaling in pool richnnes (b.pool = 0)", cex=1.5), xlab = list(label="Scaling in dispersal (b.disp)", cex=1.5), ylab = list(label="Scaling in local density (b.local)", cex=1.5),
                                 col.regions = viridis(40), cuts=40, ylim = c(-2.05, 0.05), xlim = c(-0.05, 2.05), 
                                 at=seq(round(min(par.space.BETA_0),4)-0.005, max(par.space.BETA_0), length.out=40),
                                 scales=list(x=list(cex=1.5), y=list(cex=1.5)),
                                 colorkey = list(labels=list(cex=1.2))) + 
                                 layer_(panel.2dsmoother(..., n = 600))

Plot.par.space.BETA_0.1<-levelplot(value ~ Var1 * Var2, reshape2::melt(par.space.BETA_0.1), panel = panel.levelplot.points, cex = 0, contour=T, 
                                   main="Scaling in pool richnnes (b.pool = -0.1)", xlab = "Scaling in dispersal (b.disp)", ylab = "Scaling in local density (b.local)",
                                   col.regions = viridis(40), cuts=40, ylim = c(-2.05, 0.05), xlim = c(-0.05, 2.05), 
                                   at=seq(round(min(par.space.BETA_0.1),4)-0.005, max(par.space.BETA_0.1), length.out=40)) + 
                                   layer_(panel.2dsmoother(..., n = 600))

Plot.par.space.BETA_0.2<-levelplot(value ~ Var1 * Var2, reshape2::melt(par.space.BETA_0.2), panel = panel.levelplot.points, cex = 0, contour=T, 
                                   main=list("Scaling in pool richnnes (b.pool = -0.2)", cex=1.5), xlab = list(label="Scaling in dispersal (b.disp)", cex=1.5), ylab = list(label="Scaling in local density (b.local)", cex=1.5),
                                   col.regions = viridis(40), cuts=40, ylim = c(-2.05, 0.05), xlim = c(-0.05, 2.05), 
                                   at=seq(round(min(par.space.BETA_0.2),4)-0.005, max(par.space.BETA_0.2), length.out=40),
                                   scales=list(x=list(cex=1.5), y=list(cex=1.5)),
                                   colorkey = list(labels=list(cex=1.2)),
                                   labels=list(col="black", cex=0.8)) + 
                                   layer_(panel.2dsmoother(..., n = 600))

Plot.par.space.BETA_0.3<-levelplot(value ~ Var1 * Var2, reshape2::melt(par.space.BETA_0.3), panel = panel.levelplot.points, cex = 0, contour=T, 
                                   main=list("Scaling in pool richnnes (b.pool = -0.3)", cex=1.5), xlab = list("Scaling in dispersal (b.disp)", cex=1.5), ylab = list("Scaling in local density (b.local)", cex=1.5),
                                   col.regions = viridis(40), cuts=40, ylim = c(-2.05, 0.05), xlim = c(-0.05, 2.05), 
                                   at=seq(round(min(par.space.BETA_0.3),4)-0.005, max(par.space.BETA_0.3), length.out=40),
                                    labels=list( col="black", cex=0.8),
                                   scales=list(x=list(cex=1.5), y=list(cex=1.5)),
                                   colorkey = list(labels=list(cex=1.2))) + 
                                   layer_(panel.2dsmoother(..., n = 600))

Plot.par.space.BETA_0.4<-levelplot(value ~ Var1 * Var2, reshape2::melt(par.space.BETA_0.4), panel = panel.levelplot.points, cex = 0, contour=T, 
                                   main=list("Scaling in pool richnnes (b.pool = -0.4)", cex=1.5), xlab = list("Scaling in dispersal (b.disp)", cex=1.5), ylab = list("Scaling in local density (b.local)", cex=1.5),
                                   col.regions = viridis(40), cuts=40, ylim = c(-2.05, 0.05), xlim = c(-0.05, 2.05), 
                                   at=seq(round(min(par.space.BETA_0.4),4)-0.005, max(par.space.BETA_0.4), length.out=40),
                                  # labels=list("labels"=c(rep("", 32),"-0.18", "", "", "", "", "" ), col="black", cex=0.8),
                                   scales=list(x=list(cex=1.5), y=list(cex=1.5)),
                                   colorkey = list(labels=list(cex=1.2))) + 
                                   layer_(panel.2dsmoother(..., n = 600))

Plot.par.space.BETA_0.6<-levelplot(value ~ Var1 * Var2, reshape2::melt(par.space.BETA_0.6), panel = panel.levelplot.points, cex = 0, contour=T, 
                                   main=list("Scaling in pool richnnes (b.pool = -0.6)", cex=1.5), xlab = list(label="Scaling in dispersal (b.disp)", cex=1.5), ylab = list(label="Scaling in local density (b.local)", cex=1.5),
                                   col.regions = viridis(40), cuts=40, ylim = c(-2.05, 0.05), xlim = c(-0.05, 2.05), 
                                   at=seq(round(min(par.space.BETA_0.6),4)-0.005, max(par.space.BETA_0.6), length.out=40),
                                   scales=list(x=list(cex=1.5), y=list(cex=1.5)),
                                   colorkey = list(labels=list(cex=1.2))) + 
                                   layer_(panel.2dsmoother(..., n = 200))


# Gamma diversity #
par.space.GAMMA_0<-par.space.1$Mean.gamma.slope$`0`
par.space.GAMMA_0.1<-par.space.1$Mean.gamma.slope$`-0.1`
par.space.GAMMA_0.2<-par.space.1$Mean.gamma.slope$`-0.2`
par.space.GAMMA_0.3<-par.space.1$Mean.gamma.slope$`-0.3`
par.space.GAMMA_0.4<-par.space.2$Mean.gamma.slope$`-0.4`
par.space.GAMMA_0.6<-par.space.2$Mean.gamma.slope$`-0.6`

Plot.par.space.GAMMA_0<-levelplot(value ~ Var1 * Var2, reshape2::melt(par.space.GAMMA_0), panel = panel.levelplot.points, cex = 0, contour=T, 
                                  main=list("Scaling in pool richnnes (b.pool = 0)", cex=1.5), xlab = list(label="Scaling in dispersal (b.disp)", cex=1.5), ylab = list(label="Scaling in local density (b.local)", cex=1.5),
                                  col.regions = viridis(40), cuts=40, ylim = c(-2.05, 0.05), xlim = c(-0.05, 2.05), 
                                  at=seq(round(min(par.space.GAMMA_0),4)-0.005, max(par.space.GAMMA_0), length.out=40),
                                  scales=list(x=list(cex=1.5), y=list(cex=1.5)),
                                  colorkey = list(labels=list(cex=1.2))) + 
                                  layer_(panel.2dsmoother(..., n = 600))

Plot.par.space.GAMMA_0.1<-levelplot(value ~ Var1 * Var2, reshape2::melt(par.space.GAMMA_0.1), panel = panel.levelplot.points, cex = 0, contour=T, 
                                    main="Scaling in pool richnnes (b.pool = -0.1)", xlab = "Scaling in dispersal (b.disp)", ylab = "Scaling in local density (b.local)",
                                    col.regions = viridis(40), cuts=40, ylim = c(-2.05, 0.05), xlim = c(-0.05, 2.05), 
                                    at=seq(round(min(par.space.GAMMA_0.1),4)-0.005, max(par.space.GAMMA_0.1), length.out=40)) + 
                                    layer_(panel.2dsmoother(..., n = 600))

Plot.par.space.GAMMA_0.2<-levelplot(value ~ Var1 * Var2, reshape2::melt(par.space.GAMMA_0.2), panel = panel.levelplot.points, cex = 0, contour=T, 
                                    main=list("Scaling in pool richnnes (b.pool = -0.2)", cex=1.5), xlab = list(label="Scaling in dispersal (b.disp)", cex=1.5), ylab = list(label="Scaling in local density (b.local)", cex=1.5),
                                    col.regions = viridis(40), cuts=40, ylim = c(-2.05, 0.05), xlim = c(-0.05, 2.05), 
                                    at=seq(round(min(par.space.GAMMA_0.2),4)-0.005, max(par.space.GAMMA_0.2), length.out=40),
                                    scales=list(x=list(cex=1.5), y=list(cex=1.5)),
                                    colorkey = list(labels=list(cex=1.2)),
                                    labels=list(col="black", cex=0.8)) + 
                                    layer_(panel.2dsmoother(..., n = 600))

Plot.par.space.GAMMA_0.3<-levelplot(value ~ Var1 * Var2, reshape2::melt(par.space.GAMMA_0.3), panel = panel.levelplot.points, cex = 0, contour=T, 
                                    main=list("Scaling in pool richnnes (b.pool = -0.3)", cex=1.5), xlab = list(label="Scaling in dispersal (b.disp)", cex=1.5), ylab = list(label="Scaling in local density (b.local)", cex=1.5),
                                    col.regions = viridis(40), cuts=40, ylim = c(-2.05, 0.05), xlim = c(-0.05, 2.05), 
                                    at=seq(round(min(par.space.GAMMA_0.3),4)-0.005, max(par.space.GAMMA_0.3), length.out=40),
                                    # labels=list("labels"=c(rep("", 30),"-0.18", rep("", 8)), col="black", cex=0.8),
                                    scales=list(x=list(cex=1.5), y=list(cex=1.5)),
                                    colorkey = list(labels=list(cex=1.2))) +
                                    layer_(panel.2dsmoother(..., n = 600)) 

Plot.par.space.GAMMA_0.4<-levelplot(value ~ Var1 * Var2, reshape2::melt(par.space.GAMMA_0.4), panel = panel.levelplot.points, cex = 0, contour=T, 
                                    main=list("Scaling in pool richnnes (b.pool = -0.4)", cex=1.5), xlab = list(label="Scaling in dispersal (b.disp)", cex=1.5), ylab = list(label="Scaling in local density (b.local)", cex=1.5),
                                    col.regions = viridis(40), cuts=40, ylim = c(-2.05, 0.05), xlim = c(-0.05, 2.05), 
                                    at=seq(round(min(par.space.GAMMA_0.4),4)-0.005, max(par.space.GAMMA_0.4), length.out=40),
                                   # labels=list("labels"=c(rep("", 30),"-0.18", rep("", 8)), col="black", cex=0.8),
                                    scales=list(x=list(cex=1.5), y=list(cex=1.5)),
                                    colorkey = list(labels=list(cex=1.2))) + 
                                    layer_(panel.2dsmoother(..., n = 600))

Plot.par.space.GAMMA_0.6<-levelplot(value ~ Var1 * Var2, reshape2::melt(par.space.GAMMA_0.6), panel = panel.levelplot.points, cex = 0, contour=T, 
                                    main=list("Scaling in pool richnnes (b.pool = -0.6)", cex=1.5), xlab = list(label="Scaling in dispersal (b.disp)", cex=1.5), ylab = list(label="Scaling in local density (b.local)", cex=1.5),
                                    col.regions = viridis(40), cuts=40, ylim = c(-2.05, 0.05), xlim = c(-0.05, 2.05), 
                                    at=seq(round(min(par.space.GAMMA_0.6),4)-0.005, max(par.space.GAMMA_0.6), length.out=40),
                                    scales=list(x=list(cex=1.5), y=list(cex=1.5)),
                                    colorkey = list(labels=list(cex=1.2))) + 
                                    layer_(panel.2dsmoother(..., n = 200))

# All plots together # 15x25

par(mfrow=c(3,4), oma=c(2,0,2,0))
print(Plot.par.space.ALPHA_0, split=c(1, 1, 4, 3))  
trellis.focus("legend", side="right", clipp.off=TRUE, highlight=FALSE)
grid.text(expression(bold("b"[DBS])), 0.2, 1.14, hjust=0.6, vjust=2, gp=gpar(fontsize=16))
trellis.unfocus()
print(Plot.par.space.ALPHA_0.2, split=c(2, 1, 4, 3), newpage=FALSE)
trellis.focus("legend", side="right", clipp.off=TRUE, highlight=FALSE)
grid.text(expression(bold("b"[DBS])), 0.2, 1.14, hjust=0.6, vjust=2, gp=gpar(fontsize=16))
trellis.unfocus()
print(Plot.par.space.ALPHA_0.4, split=c(3, 1, 4, 3), newpage=FALSE)
trellis.focus("legend", side="right", clipp.off=TRUE, highlight=FALSE)
grid.text(expression(bold("b"[DBS])), 0.2, 1.14, hjust=0.6, vjust=2, gp=gpar(fontsize=16))
trellis.unfocus()
print(Plot.par.space.ALPHA_0.6, split=c(4, 1, 4, 3), newpage=FALSE)
trellis.focus("legend", side="right", clipp.off=TRUE, highlight=FALSE)
grid.text(expression(bold("b"[DBS])), 0.2, 1.14, hjust=0.6, vjust=2, gp=gpar(fontsize=16))
trellis.unfocus()

print(Plot.par.space.BETA_0, split=c(1, 2, 4, 3), newpage=FALSE)
trellis.focus("legend", side="right", clipp.off=TRUE, highlight=FALSE)
grid.text(expression(bold("b"[DBS])), 0.2, 1.14, hjust=0.6, vjust=2, gp=gpar(fontsize=16))
trellis.unfocus()
print(Plot.par.space.BETA_0.2, split=c(2, 2, 4, 3), newpage=FALSE)
trellis.focus("legend", side="right", clipp.off=TRUE, highlight=FALSE)
grid.text(expression(bold("b"[DBS])), 0.2, 1.14, hjust=0.6, vjust=2, gp=gpar(fontsize=16))
trellis.unfocus()
print(Plot.par.space.BETA_0.4, split=c(3, 2, 4, 3), newpage=FALSE)
trellis.focus("legend", side="right", clipp.off=TRUE, highlight=FALSE)
grid.text(expression(bold("b"[DBS])), 0.2, 1.14, hjust=0.6, vjust=2, gp=gpar(fontsize=16))
trellis.unfocus()
print(Plot.par.space.BETA_0.6, split=c(4, 2, 4, 3), newpage=FALSE)
trellis.focus("legend", side="right", clipp.off=TRUE, highlight=FALSE)
grid.text(expression(bold("b"[DBS])), 0.2, 1.14, hjust=0.6, vjust=2, gp=gpar(fontsize=16))
trellis.unfocus()

print(Plot.par.space.GAMMA_0, split=c(1, 3, 4, 3), newpage=FALSE) 
trellis.focus("legend", side="right", clipp.off=TRUE, highlight=FALSE)
grid.text(expression(bold("b"[DBS])), 0.2, 1.14, hjust=0.6, vjust=2, gp=gpar(fontsize=16))
trellis.unfocus()
print(Plot.par.space.GAMMA_0.2, split=c(2, 3, 4, 3), newpage=FALSE)
trellis.focus("legend", side="right", clipp.off=TRUE, highlight=FALSE)
grid.text(expression(bold("b"[DBS])), 0.2, 1.14, hjust=0.6, vjust=2, gp=gpar(fontsize=16))
trellis.unfocus()
print(Plot.par.space.GAMMA_0.4, split=c(3, 3, 4, 3), newpage=FALSE)
trellis.focus("legend", side="right", clipp.off=TRUE, highlight=FALSE)
grid.text(expression(bold("b"[DBS])), 0.2, 1.14, hjust=0.6, vjust=2, gp=gpar(fontsize=16))
trellis.unfocus()
print(Plot.par.space.GAMMA_0.6, split=c(4, 3, 4, 3), newpage=FALSE) 
trellis.focus("legend", side="right", clipp.off=TRUE, highlight=FALSE)
grid.text(expression(bold("b"[DBS])), 0.2, 1.14, hjust=0.6, vjust=2, gp=gpar(fontsize=16))
trellis.unfocus()

#---| (3) EMPIRICAL DIVERSITY-BODY SIZE RELATIONSHIPS IN THE METACOMMUNITY AND ASSOCIATED PARAMETER SPACE (Figure 3) |---#

# (3.1) EMPIRICAL DIVERSITY-BODY SIZE RELATIONSHIPS IN THE METACOMMUNITY #

just.for.empirical.plot<-as.data.frame(cbind(diversity.std$alpha.mean, 
                                             diversity.std$beta.add, 
                                             diversity.std$gamma, 
                                             bs.orders))
colnames(just.for.empirical.plot)<-c("alpha.mean", "beta.add", "gamma", "BS")

# Alpha diversity #
plot.empirical.alpha<-ggplot(just.for.empirical.plot, aes(x=log2(BS), y=log(alpha.mean))) +
                             geom_point(colour="deepskyblue4", size=6 , shape=16, show.legend = F) +
                             geom_point(shape = 1,colour = "deepskyblue4", size=6, show.legend = F) +
                             theme_classic() +
                             scale_size_area(max_size=12) +
                             ylab(label = "Log(Alpha)")+
                             xlab(label = "Log(Body Size)") +
                             theme(axis.text=element_text(size=15),
                                   axis.title=element_text(size=15,face="bold")) +
                             theme(axis.title.x = element_text(vjust=-1), 
                                   axis.title.y = element_text(vjust=4)) +
                             theme(plot.margin = unit(c(1,1,1,1), "cm")) +
                             geom_smooth(method = "glm", colour="black", se = T)  +
                             theme(legend.position = "none") 

# Beta diversity #
plot.empirical.beta<-ggplot(just.for.empirical.plot, aes(x=log2(BS), y=log(beta.add))) +
                            geom_point(colour="deepskyblue4", size=6 , shape=16, show.legend = F) +
                            geom_point(shape = 1,colour = "deepskyblue4", size=6, show.legend = F) +
                            theme_classic() +
                            scale_size_area(max_size=12) +
                            ylab(label = "Log(Beta)")+
                            xlab(label = "Log(Body Size)") +
                            theme(axis.text=element_text(size=15),
                                  axis.title=element_text(size=15,face="bold")) +
                            theme(axis.title.x = element_text(vjust=-1), 
                                  axis.title.y = element_text(vjust=4)) +
                            theme(plot.margin = unit(c(1,1,1,1), "cm")) +
                            geom_smooth(method = "glm", colour="black", se = T)  +
                            theme(legend.position = "none") 

# Gamma diversity #
plot.empirical.gamma<-ggplot(just.for.empirical.plot, aes(x=log2(BS), y=log(gamma))) +
                             geom_point(colour="deepskyblue4", size=6 , shape=16, show.legend = F) +
                             geom_point(shape = 1,colour = "deepskyblue4", size=6, show.legend = F) +
                             theme_classic() +
                             scale_size_area(max_size=12) +
                             ylab(label = "Log(Gamma)")+
                             xlab(label = "Log(Body Size)") +
                             theme(axis.text=element_text(size=15),
                                   axis.title=element_text(size=15,face="bold")) +
                             theme(axis.title.x = element_text(vjust=-1), 
                                   axis.title.y = element_text(vjust=4)) +
                             theme(plot.margin = unit(c(1,1,1,1), "cm")) +
                             geom_smooth(method = "glm", colour="black", se = T)  +
                             theme(legend.position = "none") 

# All plots together #
plot.empirical.alpha + plot.empirical.beta + plot.empirical.gamma 


# (3.2) PARAMETER SPACE ASSOCIATED WITH THE EMPIRICAL RESULTS INCLUDING EMPIRICAL SCALING EXPONENTS GATHERED FROM LITERATURE #

# Alpha diversity #
Plot.par.space.ALPHA_0.3_emp.exp<-Plot.par.space.ALPHA_0.3 +

                                  layer(panel.segments(x0=-0.005, x1=-0.005, y0=-0.04, y1=-1, col="red", lwd=6)) +
                                  layer(panel.segments(x0=2.005, x1=2.005, y0=-0.04, y1=-1, col="red", lwd=6)) +
                                    
                                  layer(panel.segments(x0=-0.005, x1=-0.005, y0=-1.193, y1=-1.193, col="red", lwd=10, lty=1)) +
                                  layer(panel.segments(x0=-0.005, x1=-0.005, y0=-1.396, y1=-1.396, col="red", lwd=10, lty=1)) +
                                  layer(panel.segments(x0=-0.005, x1=-0.005, y0=-1.822, y1=-1.822, col="red", lwd=10, lty=1)) +
                                    
                                  layer(panel.segments(x0=2.005, x1=2.005, y0=-1.193, y1=-1.193, col="red", lwd=10, lty=1)) +
                                  layer(panel.segments(x0=2.005, x1=2.005, y0=-1.396, y1=-1.396, col="red", lwd=10, lty=1)) +
                                  layer(panel.segments(x0=2.005, x1=2.005, y0=-1.822, y1=-1.822, col="red", lwd=10, lty=1)) +
                                    
                                  layer(panel.segments(x0=0.11, x1=0.23, y0=0.005, y1=0.005, col="violet", lwd=10)) +
                                  layer(panel.segments(x0=0.11, x1=0.23, y0=-2.005, y1=-2.005, col="violet", lwd=10)) +
                                    
                                  layer(panel.segments(x0=0.49, x1=0.49, y0=0.005, y1=0.005, col="violet", lwd=10, lty=1)) +
                                  layer(panel.segments(x0=0.49, x1=0.49, y0=-2.005, y1=-2.005, col="violet", lwd=10, lty=1))

# Beta diversity #
Plot.par.space.BETA_0.3_emp.exp<-Plot.par.space.BETA_0.3 +
  
                                  layer(panel.segments(x0=-0.005, x1=-0.005, y0=-0.04, y1=-1, col="red", lwd=6)) +
                                  layer(panel.segments(x0=2.005, x1=2.005, y0=-0.04, y1=-1, col="red", lwd=6)) +
                                  
                                  layer(panel.segments(x0=-0.005, x1=-0.005, y0=-1.193, y1=-1.193, col="red", lwd=10, lty=1)) +
                                  layer(panel.segments(x0=-0.005, x1=-0.005, y0=-1.396, y1=-1.396, col="red", lwd=10, lty=1)) +
                                  layer(panel.segments(x0=-0.005, x1=-0.005, y0=-1.822, y1=-1.822, col="red", lwd=10, lty=1)) +
                                  
                                  layer(panel.segments(x0=2.005, x1=2.005, y0=-1.193, y1=-1.193, col="red", lwd=10, lty=1)) +
                                  layer(panel.segments(x0=2.005, x1=2.005, y0=-1.396, y1=-1.396, col="red", lwd=10, lty=1)) +
                                  layer(panel.segments(x0=2.005, x1=2.005, y0=-1.822, y1=-1.822, col="red", lwd=10, lty=1)) +
                                  
                                  layer(panel.segments(x0=0.11, x1=0.23, y0=0.005, y1=0.005, col="violet", lwd=10)) +
                                  layer(panel.segments(x0=0.11, x1=0.23, y0=-2.005, y1=-2.005, col="violet", lwd=10)) +
                                  
                                  layer(panel.segments(x0=0.49, x1=0.49, y0=0.005, y1=0.005, col="violet", lwd=10, lty=1)) +
                                  layer(panel.segments(x0=0.49, x1=0.49, y0=-2.005, y1=-2.005, col="violet", lwd=10, lty=1))

# Gamma diversity #
Plot.par.space.GAMMA_0.3_emp.exp<-Plot.par.space.GAMMA_0.3 +
  
                                  layer(panel.segments(x0=-0.005, x1=-0.005, y0=-0.04, y1=-1, col="red", lwd=6)) +
                                  layer(panel.segments(x0=2.005, x1=2.005, y0=-0.04, y1=-1, col="red", lwd=6)) +
                                  
                                  layer(panel.segments(x0=-0.005, x1=-0.005, y0=-1.193, y1=-1.193, col="red", lwd=10, lty=1)) +
                                  layer(panel.segments(x0=-0.005, x1=-0.005, y0=-1.396, y1=-1.396, col="red", lwd=10, lty=1)) +
                                  layer(panel.segments(x0=-0.005, x1=-0.005, y0=-1.822, y1=-1.822, col="red", lwd=10, lty=1)) +
                                  
                                  layer(panel.segments(x0=2.005, x1=2.005, y0=-1.193, y1=-1.193, col="red", lwd=10, lty=1)) +
                                  layer(panel.segments(x0=2.005, x1=2.005, y0=-1.396, y1=-1.396, col="red", lwd=10, lty=1)) +
                                  layer(panel.segments(x0=2.005, x1=2.005, y0=-1.822, y1=-1.822, col="red", lwd=10, lty=1)) +
                                  
                                  layer(panel.segments(x0=0.11, x1=0.23, y0=0.005, y1=0.005, col="violet", lwd=10)) +
                                  layer(panel.segments(x0=0.11, x1=0.23, y0=-2.005, y1=-2.005, col="violet", lwd=10)) +
                                  
                                  layer(panel.segments(x0=0.49, x1=0.49, y0=0.005, y1=0.005, col="violet", lwd=10, lty=1)) +
                                  layer(panel.segments(x0=0.49, x1=0.49, y0=-2.005, y1=-2.005, col="violet", lwd=10, lty=1))


# All plots together # 5x18

par(mfrow=c(1,3), oma=c(2,0,2,0)) 

print(Plot.par.space.ALPHA_0.3_emp.exp, split=c(1, 1, 3, 1)) 
trellis.focus("legend", side="right", clipp.off=TRUE, highlight=FALSE)
grid.text(expression(bold("b"[DBS])), 0.2, 1.15, hjust=0.6, vjust=2, gp=gpar(fontsize=16))
trellis.unfocus()
print(Plot.par.space.BETA_0.3_emp.exp, split=c(2, 1, 3, 1), newpage=FALSE)
trellis.focus("legend", side="right", clipp.off=TRUE, highlight=FALSE)
grid.text(expression(bold("b"[DBS])), 0.2, 1.15, hjust=0.6, vjust=2, gp=gpar(fontsize=16))
trellis.unfocus()
print(Plot.par.space.GAMMA_0.3_emp.exp, split=c(3, 1, 3, 1), newpage=FALSE)
trellis.focus("legend", side="right", clipp.off=TRUE, highlight=FALSE)
grid.text(expression(bold("b"[DBS])), 0.2, 1.15, hjust=0.6, vjust=2, gp=gpar(fontsize=16))
trellis.unfocus()


#---| (4) EMPIRICAL POOL RICHNESS-BODY SIZE RELATIONSHIP IN THE METACOMMUNITY (Supporting information Figure SM4) |---#

par(mar=c(5,5,2,2))
plot(log(gamma.pool.orders2)~log2(bs.orders), 
     xlab="Log(Body size)", ylab="Log(Pool richness)", 
     pch=19, cex=2, cex.axis=1.5, cex.lab=1.5)

abline(coef = c(3.8978, -0.3552), col="red", lwd=2)


#------------------------------#
#---| END OF THE SCRIPT |------#
#------------------------------#
