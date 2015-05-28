## plots of spatial and temporal windows of contrast increments...

library(ggplot2)

x <- seq(-2, 2, l = 100)
t <- seq(0, 600, l = 100)

y_x <- dnorm(x, mean = 0, sd = 0.5)

y_x <- y_x / max(y_x)

y_t <- c(dnorm(t[1:50], mean = 240, sd = 120), dnorm(t[51:100], mean = 360, sd = 120))

y_t[t >= 240 & t <= 360] <- max(y_t)
y_t <- y_t / max(y_t)

# do plots:
fig <- qplot(x, y_x, geom = "line") + 
  xlab("Space (deg)") + ylab("") + 
  scale_y_continuous(breaks = c(0, 1), labels = c(0, "Alpha")) + 
  theme_minimal(base_size=10)

ggsave(file=paste0(getwd(),'/figs/envelope_space.pdf'),width=2,height=1.5)

fig <- qplot(t, y_t, geom = "line") + 
  xlab("Time (ms)") + ylab("") + 
  scale_y_continuous(breaks = c(0, 1), labels = c(0, "Alpha")) + 
  coord_cartesian(ylim = c(-0.05, 1.05)) +
  theme_minimal(base_size=10)

ggsave(file=paste0(getwd(),'/figs/envelope_time.pdf'),width=2,height=1.5)
