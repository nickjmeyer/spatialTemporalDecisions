rm(list = ls(all = TRUE))

library(ggplot2)
library(viridis)
library(plyr)

df = read.csv("all_res.csv", header = TRUE)


df$name = factor(df$name,
                 levels = c("2-Step Alternating",
                            "10-Step Alternating",
                            "20-Step Alternating",
                            "Thompson Sampling"),
                 labels = c("k = 1",
                            "k = 5",
                            "k = 10",
                            "Thompson Sampling"))

df = aggregate(df[,c("pull", "value", "regret")], by = list("name" = df$name, "time" = df$time,
                                                            "eta" = df$eta), FUN = mean)
df$eta = as.factor(df$eta)

df = ddply(df, .(name, eta), transform, cum_regret = cumsum(regret))

df = ddply(df, .(name, eta), transform, cum_value = cumsum(value))

df$avg_cum_value = df$cum_value / df$time



p = ggplot(data = df, aes(x = time, y = avg_cum_value, color = name))
p = p + geom_line()
p = p + facet_wrap( ~ eta, nrow = 1)
p = p + scale_color_discrete("Strategy")
p = p + xlab("Time")
p = p + ylab("Average cumulative value")
print(p)

ggsave("../data/plotting/bandit_example.pdf", p)
