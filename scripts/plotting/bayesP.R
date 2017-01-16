rm(list=ls(all=TRUE))

dirName = "../../data/wns/2016-11-03-16-24-26/"

## spatial
ss_sim = read.table(paste(dirName,
                          "sampStats_mean_gravity2_spatial_.txt",
                          sep=""),header=TRUE)

ss_obs = read.table(paste(dirName,
                          "obsStats_spatial_.txt",
                          sep=""),header=TRUE)


for(i in 1:ncol(ss_sim)) {
  allvals = c(ss_sim[,i],ss_obs[,i])
  pdf(paste("../../data/plotting/ppc_spatial_", names(ss_sim)[i], ".pdf", sep=""))
  boxplot(ss_sim[,i],ylim=range(allvals))
  abline(h=ss_obs[,i])
  dev.off()
}


## spatial oos
ss_sim = read.table(paste(dirName,
                          "sampStats_mle_Oos_gravity2_spatial_.txt",
                          sep=""),header=TRUE)

ss_obs = read.table(paste(dirName,
                          "obsStats_Oos_spatial_.txt",
                          sep=""),header=TRUE)


for(i in 1:ncol(ss_sim)) {
  allvals = c(ss_sim[,i],ss_obs[,i])
  pdf(paste("../../data/plotting/ppc_oos_spatial_", names(ss_sim)[i], ".pdf", sep=""))
  boxplot(ss_sim[,i],ylim=range(allvals))
  abline(h=ss_obs[,i])
  dev.off()
}



## edge
ss_sim = read.table(paste(dirName,
                          "sampStats_mean_edgeToEdge2_edge_.txt",
                          sep=""),header=TRUE)

ss_obs = read.table(paste(dirName,
                          "obsStats_spatial_.txt", ## obs stats are same for edge and spatial
                          sep=""),header=TRUE)


for(i in 1:ncol(ss_sim)) {
  allvals = c(ss_sim[,i],ss_obs[,i])
  pdf(paste("../../data/plotting/ppc_edge_", names(ss_sim)[i], ".pdf", sep=""))
  boxplot(ss_sim[,i],ylim=range(allvals))
  abline(h=ss_obs[,i])
  dev.off()
}


## edge oos
ss_sim = read.table(paste(dirName,
                          "sampStats_mle_Oos_edgeToEdge2_edge_.txt",
                          sep=""),header=TRUE)

ss_obs = read.table(paste(dirName,
                          "obsStats_Oos_spatial_.txt", ## obs stats are same for edge and spatial
                          sep=""),header=TRUE)


for(i in 1:ncol(ss_sim)) {
  allvals = c(ss_sim[,i],ss_obs[,i])
  pdf(paste("../../data/plotting/ppc_oos_edge_", names(ss_sim)[i], ".pdf", sep=""))
  boxplot(ss_sim[,i],ylim=range(allvals))
  abline(h=ss_obs[,i])
  dev.off()
}
