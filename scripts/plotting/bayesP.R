rm(list=ls(all=TRUE))

dirName = "../../data/wns/2016-11-03-16-24-26"

ss_sim = read.table(paste(dirName,
                          "sampStats_mle_Oos_gravity2_spatial_.txt",
                          sep=""),header=TRUE)

ss_obs = read.table(paste(dirName,
                          "obsStats_Oos_spatial_.txt",
                          sep=""),header=TRUE)


for(i in 1:ncol(ss_sim)) {
  allvals = c(ss_sim[,i],ss_obs[,i])
  boxplot(ss_sim[,i],ylim=range(allvals),main=names(ss_sim)[i])
  abline(h=ss_obs[,i])
  readline("Press [enter]")
}
