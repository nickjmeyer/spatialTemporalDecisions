rm(list=ls(all=TRUE))

library(ggplot2)

res = read.csv("../../data/results/simResults.csv", header = TRUE)

res = subset(res, !is.na(res$network))

res$agent = factor(res$agent,
                   levels = c("ps", "proximal", "none", "myopic", "all"),
                   labels = c("Policy Search", "Proximal", "None", "Myopic",
                              "All"))
res$network = factor(res$network,
                     levels = c("grid", "rand", "crp", "scalefree"),
                     labels = c("N1", "N2", "N3", "N4"))
res$model = factor(res$model,
                   levels = c("corr", "miss"),
                   labels = c("Correct", "Misspecified"))

resEdgeMiss = subset(res, res$mode == "edge" & res$model == "Correct")
resSpatialMiss = subset(res, res$mode == "spatial" & res$model == "Misspecified")
resEdgeCorr = subset(res, res$mode == "edge" & res$model == "Correct")
resSpatialCorr = subset(res, res$mode == "spatial" & res$model == "Misspecified")


pEdgeMiss = ggplot(data = resEdgeMiss, aes(y = agent, x = mean, color = model, pch = model,
                                   size = model))
pEdgeMiss = pEdgeMiss + geom_point(stroke = 1)
pEdgeMiss = pEdgeMiss + geom_errorbarh(aes(xmin = mean - sd, xmax = mean + sd),
                                       size = 0.5, height = 0.)
pEdgeMiss = pEdgeMiss + scale_shape_manual("Model Specification", values = c(1, 3))
pEdgeMiss = pEdgeMiss + scale_size_manual("Model Specification", values = c(2, 1.25))
pEdgeMiss = pEdgeMiss + facet_grid(network ~ size)
pEdgeMiss = pEdgeMiss + scale_color_hue("Model Specification")
pEdgeMiss = pEdgeMiss + scale_y_discrete(limits = c("None", "Myopic", "Proximal", "Policy Search", "All"))
pEdgeMiss = pEdgeMiss + scale_x_continuous(breaks = pretty(resEdgeMiss$mean, n = 4))
pEdgeMiss = pEdgeMiss + xlab("Estimated Mean Value")
pEdgeMiss = pEdgeMiss + ylab("Treatment Strategy")
pEdgeMiss = pEdgeMiss + theme(legend.position = "none")

ggsave("../../data/plotting/results_plot_edge_miss.pdf", pEdgeMiss)

pEdgeCorr = ggplot(data = resEdgeCorr, aes(y = agent, x = mean, color = model, pch = model,
                                   size = model))
pEdgeCorr = pEdgeCorr + geom_point(stroke = 1)
pEdgeCorr = pEdgeCorr + geom_errorbarh(aes(xmin = mean - sd, xmax = mean + sd),
                                       size = 0.5, height = 0)
pEdgeCorr = pEdgeCorr + scale_shape_manual("Model Specification", values = c(1, 3))
pEdgeCorr = pEdgeCorr + scale_size_manual("Model Specification", values = c(2, 1.25))
pEdgeCorr = pEdgeCorr + facet_grid(network ~ size)
pEdgeCorr = pEdgeCorr + scale_color_hue("Model Specification")
pEdgeCorr = pEdgeCorr + scale_y_discrete(limits = c("None", "Myopic", "Proximal", "Policy Search", "All"))
pEdgeCorr = pEdgeCorr + scale_x_continuous(breaks = pretty(resEdgeCorr$mean, n = 4))
pEdgeCorr = pEdgeCorr + xlab("Estimated Mean Value")
pEdgeCorr = pEdgeCorr + ylab("Treatment Strategy")
pEdgeCorr = pEdgeCorr + theme(legend.position = "none")

ggsave("../../data/plotting/results_plot_edge_corr.pdf", pEdgeCorr)

pSpatialMiss = ggplot(data = resSpatialMiss, aes(y = agent, x = mean, color = model, pch = model,
                                   size = model))
pSpatialMiss = pSpatialMiss + geom_point(stroke = 1)
pSpatialMiss = pSpatialMiss + geom_errorbarh(aes(xmin = mean - sd, xmax = mean + sd),
                                       size = 0.5, height = 0)
pSpatialMiss = pSpatialMiss + scale_shape_manual("Model Specification", values = c(1, 3))
pSpatialMiss = pSpatialMiss + scale_size_manual("Model Specification", values = c(2, 1.25))
pSpatialMiss = pSpatialMiss + facet_grid(network ~ size)
pSpatialMiss = pSpatialMiss + scale_color_hue("Model Specification")
pSpatialMiss = pSpatialMiss + scale_y_discrete(limits = c("None", "Myopic", "Proximal", "Policy Search", "All"))
pSpatialMiss = pSpatialMiss + scale_x_continuous(breaks = pretty(resSpatialMiss$mean, n = 4))
pSpatialMiss = pSpatialMiss + xlab("Estimated Mean Value")
pSpatialMiss = pSpatialMiss + ylab("Treatment Strategy")
pSpatialMiss = pSpatialMiss + theme(legend.position = "none")

ggsave("../../data/plotting/results_plot_spatial_miss.pdf", pSpatialMiss)

pSpatialCorr = ggplot(data = resSpatialCorr, aes(y = agent, x = mean, color = model, pch = model,
                                   size = model))
pSpatialCorr = pSpatialCorr + geom_point(stroke = 1)
pSpatialCorr = pSpatialCorr + geom_errorbarh(aes(xmin = mean - sd, xmax = mean + sd),
                                       size = 0.5, height = 0)
pSpatialCorr = pSpatialCorr + scale_shape_manual("Model Specification", values = c(1, 3))
pSpatialCorr = pSpatialCorr + scale_size_manual("Model Specification", values = c(2, 1.25))
pSpatialCorr = pSpatialCorr + facet_grid(network ~ size)
pSpatialCorr = pSpatialCorr + scale_color_hue("Model Specification")
pSpatialCorr = pSpatialCorr + scale_y_discrete(limits = c("None", "Myopic", "Proximal", "Policy Search", "All"))
pSpatialCorr = pSpatialCorr + scale_x_continuous(breaks = pretty(resSpatialCorr$mean, n = 4))
pSpatialCorr = pSpatialCorr + xlab("Estimated Mean Value")
pSpatialCorr = pSpatialCorr + ylab("Treatment Strategy")
pSpatialCorr = pSpatialCorr + theme(legend.position = "none")

ggsave("../../data/plotting/results_plot_spatial_corr.pdf", pSpatialCorr)
