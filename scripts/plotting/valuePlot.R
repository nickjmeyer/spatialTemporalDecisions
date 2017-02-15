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

resEdge = subset(res, res$mode == "edge")
resSpatial = subset(res, res$mode == "spatial")


pEdge = ggplot()
pEdge = pEdge + geom_point(data = resEdge, aes(y = agent, x = mean, color = model, pch = model,
                                               size = model), stroke = 1.5)
pEdge = pEdge + scale_shape_manual("Model Specification", values = c(1, 3))
pEdge = pEdge + scale_size_manual("Model Specification", values = c(2, 1.25))
pEdge = pEdge + facet_grid(network ~ size)
pEdge = pEdge + xlab("Estimated Mean Value")
pEdge = pEdge + ylab("Treatment Strategy")
pEdge = pEdge + scale_color_hue("Model Specification")
pEdge = pEdge + theme(legend.position = "bottom")


ggsave("../../data/plotting/results_plot_edge.pdf", pEdge)

pSpatial = ggplot()
pSpatial = pSpatial + geom_point(data = resSpatial, aes(y = agent, x = mean, color = model, pch = model,
                                               size = model), stroke = 1.5)
pSpatial = pSpatial + scale_shape_manual("Model Specification", values = c(1, 3))
pSpatial = pSpatial + scale_size_manual("Model Specification", values = c(2, 1.25))
pSpatial = pSpatial + facet_grid(network ~ size)
pSpatial = pSpatial + xlab("Estimated Mean Value")
pSpatial = pSpatial + ylab("Treatment Strategy")
pSpatial = pSpatial + scale_color_hue("Model Specification")
pSpatial = pSpatial + theme(legend.position = "bottom")

ggsave("../../data/plotting/results_plot_spatial.pdf", pSpatial)
