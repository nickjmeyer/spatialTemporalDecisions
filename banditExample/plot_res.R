rm(list = ls(all = TRUE))

library(ggplot2)
library(viridis)
library(plyr)

df = read.csv("all_res.csv", header = TRUE)

df$name = factor(df$name, levels = c("2-Step Alternating",
                                     "10-Step Alternating",
                                     "20-Step Alternating",
                                     "Thompson Sampling"))

df = aggregate(df[,c("pull", "value", "regret")], by = list("name" = df$name, "time" = df$time), FUN = mean)

df = ddply(df, .(name), transform, cum_regret = cumsum(regret))


p = ggplot()
p = p + geom_line(data = df,
                  aes(x = time, y = pull, color = name, lty = name))
p = p + scale_color_discrete("Strategy")
p = p + scale_linetype_manual("Strategy",
                              values = c("2-Step Alternating" = "solid",
                                         "10-Step Alternating" = "solid",
                                         "20-Step Alternating" = "solid",
                                         "Thompson Sampling" = "dashed"))
p = p + ylab("Probability of selecting the wrong arm")
p = p + xlab("Time")
p = p + theme(legend.position = "bottom")

ggsave("../data/plotting/bandit_example_mistakes.pdf", p)


p = ggplot()
p = p + geom_line(data = df,
                  aes(x = time, y = cum_regret, color = name, lty = name))
p = p + scale_color_discrete("Strategy")
p = p + scale_linetype_manual("Strategy",
                              values = c("2-Step Alternating" = "solid",
                                         "10-Step Alternating" = "solid",
                                         "20-Step Alternating" = "solid",
                                         "Thompson Sampling" = "dashed"))
p = p + ylab("Expected cumulative regret")
p = p + xlab("Time")
p = p + theme(legend.position = "bottom")

ggsave("../data/plotting/bandit_example_regret.pdf", p)
