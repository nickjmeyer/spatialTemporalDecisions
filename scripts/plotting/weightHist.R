rm(list=ls(all=TRUE))

dirName = "../../data/wns/2017-03-20-11-15-53/"

allFiles = dir(dirName)

weightsFiles = c()
for(f in allFiles) {
  if (grepl("^rank_M1Sp_weights_[0-9]{1,2}_.txt$",f)) {
    weightsFiles = c(weightsFiles,f)
  }
}

weights = c()
for(f in weightsFiles) {
  w = read.table(paste(dirName, f, sep = ""), header = FALSE)
  stopifnot(nrow(w) == 18)
  weights = rbind(weights, w[nrow(w),])
}
weights = as.matrix(weights)

pdf("../../data/plotting/weights.pdf")
boxplot(weights, names = c(expression(eta[1]),expression(eta[2]),expression(eta[3])))
dev.off()
