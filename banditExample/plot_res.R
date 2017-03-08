rm(list = ls(all = TRUE))

library(ggplot2)
library(viridis)

df = read.csv("all_res.csv", header = TRUE)

df$name = factor(df$name, levels = c("2-Step Alternating",
                                     "10-Step Alternating",
                                     "20-Step Alternating",
                                     "Thompson Sampling"))


df = aggregate(df[,c("pull", "value")], by = list("name" = df$name, "time" = df$time), FUN = mean)


p = ggplot()
p = p + geom_line(data = df,
                  aes(x = time, y = pull, color = name, lty = name))
p = p + scale_color_discrete("Strategy")
p = p + scale_linetype_manual("Strategy",
                              values = c("2-Step Alternating" = "solid",
                                         "10-Step Alternating" = "solid",
                                         "20-Step Alternating" = "solid",
                                         "Thompson Sampling" = "dashed"))
p = p + ylab("Estimated probability of mistake")
p = p + xlab("Time")
p = p + theme(legend.position = "bottom")

ggsave("../data/plotting/bandit_example.pdf", p)
