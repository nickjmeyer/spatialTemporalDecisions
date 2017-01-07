rm(list=ls(all=TRUE))

dirName = "../../data/wns/2016-11-03-16-24-26"

spatialPar = read.table(paste(dirName, "/sampStats_gravity2_param_spatial_.txt", sep=""),
                        header = FALSE)

edgePar = read.table(paste(dirName, "/sampStats_edgeToEdge2_param_edge_.txt", sep=""),
                        header = FALSE)


pdf("../../data/plotting/covariate_post_dist_spatial_uninfected.pdf")
boxplot(spatialPar[,2:5], ylim = range(spatialPar[,2:9]),
        names = c("Caves", "Cold days", "Area", "SR"))
dev.off()

pdf("../../data/plotting/covariate_post_dist_spatial_infected.pdf")
boxplot(spatialPar[,6:9], ylim = range(spatialPar[,2:9]),
        names = c("Caves", "Cold days", "Area", "SR"))
dev.off()

pdf("../../data/plotting/covariate_post_dist_edge_uninfected.pdf")
boxplot(edgePar[,2:5], ylim = range(edgePar[,2:9]),
        names = c("Caves", "Cold days", "Area", "SR"))
dev.off()

pdf("../../data/plotting/covariate_post_dist_edge_infected.pdf")
boxplot(edgePar[,6:9], ylim = range(edgePar[,2:9]),
        names = c("Caves", "Cold days", "Area", "SR"))
dev.off()
