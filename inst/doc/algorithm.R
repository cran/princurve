## ----setup, echo=FALSE, message=FALSE------------------------------------
knitr::opts_chunk$set(comment = "#>", collapse = TRUE, echo = FALSE)
set.seed(1)
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(princurve)

## ----iterative, fig.show='animate', ffmpeg.format='gif'------------------
set.seed(1)
z <- sort(runif(100, -1.4 * pi, .4 * pi))
s <- data_frame(
  x = cos(z) * 1.5,
  y = sin(z)
)
x <- s %>% 
  sample_frac(1) %>% 
  mutate(
    x = x + rnorm(length(x), 0, .05),
    y = y + rnorm(length(x), 0, .05)
  )

for (it in seq(0, 10)) {
  fit <- principal_curve(as.matrix(x), maxit = it)
  g <- ggplot() +
    geom_point(aes(x, y), x, colour = "darkgray") +
    geom_path(aes(x, y), as_data_frame(fit$s[fit$ord, ])) +
    theme_bw()+
    coord_cartesian(xlim = c(-1.6, 1.6), ylim = c(-1.1, 1.1)) +
    labs(title = paste0("Iteration ", it)) 
  print(g)
}

## ----initialisation, fig.show='animate', ffmpeg.format='gif'-------------
fit0 <- principal_curve(as.matrix(x), maxit = 0)

ggplot() +
  geom_point(aes(x, y), x) +
  theme_bw() +
  coord_cartesian(xlim = c(-1.6, 1.6), ylim = c(-1.1, 1.1)) +
  labs(title = "Step 0a: Initialise curve with principal component\nStep 0b: Orthogonally project points to curve\nStep 0c: Calculate arc-length of projections w.r.t. the origin of the curve")

ggplot() +
  geom_point(aes(x, y), x, colour = "darkgray") +
  geom_path(aes(x, y), as_data_frame(fit0$s[fit0$ord, ])) +
  theme_bw() +
  coord_cartesian(xlim = c(-1.6, 1.6), ylim = c(-1.1, 1.1)) +
  labs(title = "> Step 0a: Initialise curve with principal component\nStep 0b: Orthogonally project points to curve\nStep 0c: Calculate arc-length of projections w.r.t. the origin of the curve")

ggplot() +
  geom_segment(aes(x = x$x, xend = fit0$s[,1], y = x$y, yend = fit0$s[,2]), linetype = "dashed") +
  geom_path(aes(x, y), as_data_frame(fit0$s[fit0$ord, ]), colour = "darkgray") +
  geom_point(aes(x, y), x, colour = "darkgray") +
  theme_bw() +
  coord_cartesian(xlim = c(-1.6, 1.6), ylim = c(-1.1, 1.1)) +
  labs(x = "x", y = "y", title = "Step 0a: Initialise curve with principal component\n> Step 0b: Orthogonally project points to curve\nStep 0c: Calculate arc-length of projections w.r.t. the origin of the curve")

ggplot() +
  geom_segment(aes(x = x$x, xend = fit0$s[,1], y = x$y, yend = fit0$s[,2]), linetype = "dashed", colour = "lightgray") +
  geom_path(aes(x, y), as_data_frame(fit0$s[fit0$ord, ]), colour = "darkgray") +
  geom_point(aes(x, y), x, colour = "darkgray") +
  geom_point(aes(x, y), as_data_frame(fit0$s)) +
  theme_bw() +
  coord_cartesian(xlim = c(-1.6, 1.6), ylim = c(-1.1, 1.1)) +
  labs(x = "x", y = "y", title = "Step 0a: Initialise curve with principal component\nStep 0b: Orthogonally project points to curve\n> Step 0c: Calculate arc-length of projections w.r.t. the origin of the curve")

## ----smooth, fig.show='animate', ffmpeg.format='gif'---------------------
fit1 <- principal_curve(as.matrix(x), maxit = 1)

xdf <- x %>% 
  mutate(lambda = fit0$lambda) %>%
  arrange(lambda) %>% 
  gather(dimension, value, -lambda) %>% 
  group_by(dimension) %>% 
  mutate(smooth = smooth.spline(lambda + runif(length(lambda), 0, .001), value, df = 5)$y) %>% 
  ungroup()

ggplot(xdf) +
  geom_point(aes(lambda, value)) +
  geom_line(aes(lambda, smooth)) +
  facet_wrap(~dimension, scales = "free") +
  theme_bw()

## ----beforeafter, fig.show='animate', ffmpeg.format='gif'----------------
ggplot() +
  geom_point(aes(x, y), x, colour = "darkgray") +
  geom_path(aes(x, y), as_data_frame(fit0$s[fit0$ord, ])) +
  theme_bw() +
  coord_cartesian(xlim = c(-1.6, 1.6), ylim = c(-1.1, 1.1)) +
  labs(x = "x", y = "y", title = "Before")

ggplot() +
  geom_point(aes(x, y), x, colour = "darkgray") +
  geom_path(aes(x, y), as_data_frame(fit1$s[fit1$ord, ])) +
  theme_bw() +
  coord_cartesian(xlim = c(-1.6, 1.6), ylim = c(-1.1, 1.1)) +
  labs(x = "x", y = "y", title = "After")

## ----approx, fig.show='animate', ffmpeg.format='gif'---------------------
xout <- seq(min(xdf$lambda), max(xdf$lambda), length.out = 100)
xadf <- xdf %>% 
  group_by(dimension) %>% 
  do({
    data_frame(lambda = xout, dimension = .$dimension[[1]], smooth = stats::approx(.$lambda, .$smooth, xout)$y)
  }) %>% 
  ungroup() %>% 
  spread(dimension, smooth)

ggplot() +
  geom_point(aes(x, y), x, colour = "darkgray") +
  geom_path(aes(x, y), as_data_frame(fit1$s[fit1$ord, ])) +
  geom_point(aes(x, y), as_data_frame(fit1$s)) +
  theme_bw() +
  coord_cartesian(xlim = c(-1.6, 1.6), ylim = c(-1.1, 1.1)) +
  labs(x = "x", y = "y", title = "Before")

ggplot() +
  geom_point(aes(x, y), x, colour = "darkgray") +
  geom_path(aes(x, y), xadf) +
  geom_point(aes(x, y), xadf) +
  theme_bw() +
  coord_cartesian(xlim = c(-1.6, 1.6), ylim = c(-1.1, 1.1)) +
  labs(x = "x", y = "y", title = "After")

## ----relambda------------------------------------------------------------
ggplot() +
  geom_segment(aes(x = x$x, xend = fit1$s[,1], y = x$y, yend = fit1$s[,2]), linetype = "dashed", colour = "lightgray") +
  geom_path(aes(x, y), as_data_frame(fit1$s[fit1$ord, ]), colour = "darkgray") +
  geom_point(aes(x, y), x, colour = "darkgray") +
  geom_point(aes(x, y), as_data_frame(fit1$s)) +
  theme_bw() +
  coord_cartesian(xlim = c(-1.6, 1.6), ylim = c(-1.1, 1.1)) +
  labs(x = "x", y = "y")

## ----smooth2, fig.show='animate', ffmpeg.format='gif'--------------------
fit2 <- principal_curve(as.matrix(x), maxit = 2)

xdf2 <- x %>% 
  mutate(lambda = fit1$lambda) %>%
  arrange(lambda) %>% 
  gather(dimension, value, -lambda) %>% 
  group_by(dimension) %>% 
  mutate(smooth = smooth.spline(lambda + runif(length(lambda), 0, .001), value, df = 5)$y) %>% 
  ungroup()

ggplot(xdf2) +
  geom_point(aes(lambda, value)) +
  geom_line(aes(lambda, smooth)) +
  facet_wrap(~dimension, scales = "free") +
  theme_bw()

## ----beforeafter2, fig.show='animate', ffmpeg.format='gif'---------------
ggplot() +
  geom_point(aes(x, y), x, colour = "darkgray") +
  geom_path(aes(x, y), as_data_frame(fit1$s[fit1$ord, ])) +
  theme_bw() +
  coord_cartesian(xlim = c(-1.6, 1.6), ylim = c(-1.1, 1.1)) +
  labs(x = "x", y = "y", title = "Before")

ggplot() +
  geom_point(aes(x, y), x, colour = "darkgray") +
  geom_path(aes(x, y), as_data_frame(fit2$s[fit2$ord, ])) +
  theme_bw() +
  coord_cartesian(xlim = c(-1.6, 1.6), ylim = c(-1.1, 1.1)) +
  labs(x = "x", y = "y", title = "After")

## ----approx2, fig.show='animate', ffmpeg.format='gif'--------------------
xout2 <- seq(min(xdf2$lambda), max(xdf2$lambda), length.out = 100)
xadf2 <- xdf2 %>% 
  group_by(dimension) %>% 
  do({
    data_frame(lambda = xout2, dimension = .$dimension[[1]], smooth = stats::approx(.$lambda, .$smooth, xout2)$y)
  }) %>% 
  ungroup() %>% 
  spread(dimension, smooth)

ggplot() +
  geom_point(aes(x, y), x, colour = "darkgray") +
  geom_path(aes(x, y), as_data_frame(fit2$s[fit2$ord, ])) +
  geom_point(aes(x, y), as_data_frame(fit2$s)) +
  theme_bw() +
  coord_cartesian(xlim = c(-1.6, 1.6), ylim = c(-1.1, 1.1)) +
  labs(x = "x", y = "y", title = "Before")

ggplot() +
  geom_point(aes(x, y), x, colour = "darkgray") +
  geom_path(aes(x, y), xadf2) +
  geom_point(aes(x, y), xadf2) +
  theme_bw() +
  coord_cartesian(xlim = c(-1.6, 1.6), ylim = c(-1.1, 1.1)) +
  labs(x = "x", y = "y", title = "After")

## ----relambda2-----------------------------------------------------------
ggplot() +
  geom_segment(aes(x = x$x, xend = fit2$s[,1], y = x$y, yend = fit2$s[,2]), linetype = "dashed", colour = "lightgray") +
  geom_path(aes(x, y), as_data_frame(fit2$s[fit2$ord, ]), colour = "darkgray") +
  geom_point(aes(x, y), x, colour = "darkgray") +
  geom_point(aes(x, y), as_data_frame(fit2$s)) +
  theme_bw() +
  coord_cartesian(xlim = c(-1.6, 1.6), ylim = c(-1.1, 1.1)) +
  labs(x = "x", y = "y")

