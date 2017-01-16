rm(list=ls(all=TRUE))

library(ggplot2)

dirName = "../../data/wns/2016-11-03-16-24-26/"


clean_names = c("Total infections", "Infections in 2007",
                "Infections in 2008", "Infections in 2009", "Infections in 2010",
                "Infections in 2011", "Infections in 2012", "Infections in 2013",
                "Mean infection year", "Mean infection longitude",
                "Mean infection latitude", "Mean distance of infections to origin",
                "Smallest longitude of infections", "Smallest latitude of infections",
                "Largest longitude of infections", "Largest latitude of infections",
                "Max distance of infections to origin")

ss_sim_spat = read.table(paste(dirName,
                          "sampStats_mean_gravity2_spatial_.txt",
                          sep=""),header=TRUE)

ss_sim_edge = read.table(paste(dirName,
                          "sampStats_mean_edgeToEdge2_edge_.txt",
                          sep=""),header=TRUE)

ss_obs = read.table(paste(dirName,
                          "obsStats_spatial_.txt", ## observed stats are same for edge and spatial
                          sep=""),header=TRUE)

stopifnot(length(clean_names) == ncol(ss_obs))


for(i in 1:ncol(ss_sim_spat)) {
  rng = range(c(ss_sim_spat[,i], ss_sim_edge[i], ss_obs[,i]))
  df = as.data.frame(rbind(cbind(ss_sim_spat[,i], 1), cbind(ss_sim_edge[,i], 0)))
  df[,2] = factor(df[,2], levels = c(1, 0), labels = c("Spatial", "Network"))
  names(df) = c("value", "model")

  p = ggplot(df, aes(x = model, y = value))
  p = p + geom_boxplot()
  p = p + coord_cartesian(ylim = rng)
  p = p + geom_hline(yintercept = ss_obs[,i])
  p = p + xlab("")
  p = p + ylab(clean_names[i])
  p = p + theme_classic()
  p = p + theme(text = element_text(size = 15))

  ggsave(filename = paste("../../data/plotting/ppc_", names(ss_obs)[i], ".pdf", sep=""),
         plot = p)
}



## oos
clean_names = c("Total infections", "Infections in 2012", "Infections in 2013",
                "Mean infection year", "Mean infection longitude",
                "Mean infection latitude", "Mean distance of infections to origin",
                "Smallest longitude of infections", "Smallest latitude of infections",
                "Largest longitude of infections", "Largest latitude of infections",
                "Max distance of infections to origin")

ss_sim_spat = read.table(paste(dirName,
                          "sampStats_mle_Oos_gravity2_spatial_.txt",
                          sep=""),header=TRUE)

ss_sim_edge = read.table(paste(dirName,
                          "sampStats_mle_Oos_edgeToEdge2_edge_.txt",
                          sep=""),header=TRUE)

ss_obs = read.table(paste(dirName,
                          "obsStats_Oos_spatial_.txt", ## obs stats are same for edge and spatial
                          sep=""),header=TRUE)


stopifnot(length(clean_names) == ncol(ss_obs))


for(i in 1:ncol(ss_sim_spat)) {
  rng = range(c(ss_sim_spat[,i], ss_sim_edge[i], ss_obs[,i]))
  df = as.data.frame(rbind(cbind(ss_sim_spat[,i], 1), cbind(ss_sim_edge[,i], 0)))
  df[,2] = factor(df[,2], levels = c(1, 0), labels = c("Spatial", "Network"))
  names(df) = c("value", "model")

  p = ggplot(df, aes(x = model, y = value))
  p = p + geom_boxplot()
  p = p + coord_cartesian(ylim = rng)
  p = p + geom_hline(yintercept = ss_obs[,i])
  p = p + xlab("")
  p = p + ylab(clean_names[i])
  p = p + theme_classic()
  p = p + theme(text = element_text(size = 15))

  ggsave(filename = paste("../../data/plotting/ppc_oos_", names(ss_obs)[i], ".pdf", sep=""),
         plot = p)
}
