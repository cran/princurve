## ----setup, echo = FALSE, message = FALSE--------------------------------
knitr::opts_chunk$set(comment = "#>", echo = FALSE, fig.width = 8, fig.height = 6)
set.seed(1)
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(princurve)

library(magick)

ggif_list <- function(list, .width = 8, .height = 6, .dpi = 120, .fps = 1, ...) {
  dir <- tempfile("gif_files")
  dir.create(dir)
  on.exit(unlink(dir))
  
  img <- lapply(
    seq_along(list),
    function(i) {
      filename <- paste0(dir, "/image-", i, ".png")
      ggsave(filename, list[[i]], width = .width, height = .height, dpi = .dpi)
      image_read(filename)
    }
  )
  
  image_animate(do.call(c, img), fps = .fps, dispose = "none")
}

ggif_lapply <- function(X, FUN, .width = 8, .height = 6, .dpi = 120, .fps = 1, ...) {
  list <- lapply(X, FUN)
  ggif_list(list, .width = .width, .height = .height, .dpi = .dpi, .fps = .fps, ...)
}

## ----dataset-------------------------------------------------------------
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

## ----iterative-----------------------------------------------------------
ggif_lapply(seq(0, 10), function(it) {
  fit <- principal_curve(as.matrix(x), maxit = it)
  
  curve <- 
    as_data_frame(fit$s) %>% 
    mutate(lambda = fit$lambda, it = it) %>% 
    slice(fit$ord) %>% 
    mutate(pos = seq_len(n()))
  
  ggplot() +
    geom_point(aes(x, y), x, colour = "darkgray") +
    geom_path(aes(x, y), curve) +
    theme_bw() +
    coord_cartesian(xlim = c(-1.6, 1.6), ylim = c(-1.1, 1.1)) +
    labs(title = paste0("Iteration ", it)) 
})

## ----initialisation------------------------------------------------------
fit0 <- principal_curve(as.matrix(x), maxit = 0)

steps <- c(
  "Step 0a: Initialise curve with principal component", 
  "Step 0b: Orthogonally project points to curve",
  "Step 0c: Calculate arc-length of projections w.r.t. the origin of the curve"
)

g0 <- ggplot() +
  geom_point(aes(x, y), x) +
  theme_bw() +
  coord_cartesian(xlim = c(-1.6, 1.6), ylim = c(-1.1, 1.1)) +
  labs(x = "x", y = "y", title = paste0(c("", "", ""), steps, collapse = "\n"))

g1 <- ggplot() +
  geom_point(aes(x, y), x, colour = "darkgray") +
  geom_path(aes(x, y), as_data_frame(fit0$s[fit0$ord, ])) +
  theme_bw() +
  coord_cartesian(xlim = c(-1.6, 1.6), ylim = c(-1.1, 1.1)) +
  labs(x = "x", y = "y", title = paste0(c("> ", "", ""), steps, collapse = "\n"))

g2 <- ggplot() +
  geom_segment(aes(x = x$x, xend = fit0$s[,1], y = x$y, yend = fit0$s[,2]), linetype = "dashed") +
  geom_path(aes(x, y), as_data_frame(fit0$s[fit0$ord, ]), colour = "darkgray") +
  geom_point(aes(x, y), x, colour = "darkgray") +
  theme_bw() +
  coord_cartesian(xlim = c(-1.6, 1.6), ylim = c(-1.1, 1.1)) +
  labs(x = "x", y = "y", title = paste0(c("", "> ", ""), steps, collapse = "\n"))

g3 <- ggplot() +
  geom_segment(aes(x = x$x, xend = fit0$s[,1], y = x$y, yend = fit0$s[,2]), linetype = "dashed", colour = "lightgray") +
  geom_path(aes(x, y), as_data_frame(fit0$s[fit0$ord, ]), colour = "darkgray") +
  geom_point(aes(x, y), x, colour = "darkgray") +
  geom_point(aes(x, y), as_data_frame(fit0$s)) +
  theme_bw() +
  coord_cartesian(xlim = c(-1.6, 1.6), ylim = c(-1.1, 1.1)) +
  labs(x = "x", y = "y", title = paste0(c("", "", "> "), steps, collapse = "\n"))

ggif_list(list(g0, g1, g2, g3))

## ----smooth--------------------------------------------------------------
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

## ----beforeafter---------------------------------------------------------
g0 <- ggplot() +
  geom_point(aes(x, y), x, colour = "darkgray") +
  geom_path(aes(x, y), as_data_frame(fit0$s[fit0$ord, ])) +
  theme_bw() +
  coord_cartesian(xlim = c(-1.6, 1.6), ylim = c(-1.1, 1.1)) +
  labs(x = "x", y = "y", title = "Before")

g1 <- ggplot() +
  geom_point(aes(x, y), x, colour = "darkgray") +
  geom_path(aes(x, y), as_data_frame(fit1$s[fit1$ord, ])) +
  theme_bw() +
  coord_cartesian(xlim = c(-1.6, 1.6), ylim = c(-1.1, 1.1)) +
  labs(x = "x", y = "y", title = "After")

ggif_list(list(g0, g1))

## ----approx--------------------------------------------------------------
xout <- seq(min(xdf$lambda), max(xdf$lambda), length.out = 100)
xadf <- xdf %>% 
  group_by(dimension) %>% 
  do({
    data_frame(lambda = xout, dimension = .$dimension[[1]], smooth = stats::approx(.$lambda, .$smooth, xout)$y)
  }) %>% 
  ungroup() %>% 
  spread(dimension, smooth)

g0 <- ggplot() +
  geom_point(aes(x, y), x, colour = "darkgray") +
  geom_path(aes(x, y), as_data_frame(fit1$s[fit1$ord, ])) +
  geom_point(aes(x, y), as_data_frame(fit1$s)) +
  theme_bw() +
  coord_cartesian(xlim = c(-1.6, 1.6), ylim = c(-1.1, 1.1)) +
  labs(x = "x", y = "y", title = "Before")

g1 <- ggplot() +
  geom_point(aes(x, y), x, colour = "darkgray") +
  geom_path(aes(x, y), xadf) +
  geom_point(aes(x, y), xadf) +
  theme_bw() +
  coord_cartesian(xlim = c(-1.6, 1.6), ylim = c(-1.1, 1.1)) +
  labs(x = "x", y = "y", title = "After")

ggif_list(list(g0, g1))

## ----relambda------------------------------------------------------------
ggplot() +
  geom_segment(aes(x = x$x, xend = fit1$s[,1], y = x$y, yend = fit1$s[,2]), linetype = "dashed", colour = "lightgray") +
  geom_path(aes(x, y), as_data_frame(fit1$s[fit1$ord, ]), colour = "darkgray") +
  geom_point(aes(x, y), x, colour = "darkgray") +
  geom_point(aes(x, y), as_data_frame(fit1$s)) +
  theme_bw() +
  coord_cartesian(xlim = c(-1.6, 1.6), ylim = c(-1.1, 1.1)) +
  labs(x = "x", y = "y")

## ----smooth2-------------------------------------------------------------
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

## ----beforeafter2--------------------------------------------------------
g0 <- ggplot() +
  geom_point(aes(x, y), x, colour = "darkgray") +
  geom_path(aes(x, y), as_data_frame(fit1$s[fit1$ord, ])) +
  theme_bw() +
  coord_cartesian(xlim = c(-1.6, 1.6), ylim = c(-1.1, 1.1)) +
  labs(x = "x", y = "y", title = "Before")

g1 <- ggplot() +
  geom_point(aes(x, y), x, colour = "darkgray") +
  geom_path(aes(x, y), as_data_frame(fit2$s[fit2$ord, ])) +
  theme_bw() +
  coord_cartesian(xlim = c(-1.6, 1.6), ylim = c(-1.1, 1.1)) +
  labs(x = "x", y = "y", title = "After")

ggif_list(list(g0, g1))

## ----approx2-------------------------------------------------------------
xout2 <- seq(min(xdf2$lambda), max(xdf2$lambda), length.out = 100)
xadf2 <- xdf2 %>% 
  group_by(dimension) %>% 
  do({
    data_frame(lambda = xout2, dimension = .$dimension[[1]], smooth = stats::approx(.$lambda, .$smooth, xout2)$y)
  }) %>% 
  ungroup() %>% 
  spread(dimension, smooth)

g0 <- ggplot() +
  geom_point(aes(x, y), x, colour = "darkgray") +
  geom_path(aes(x, y), as_data_frame(fit2$s[fit2$ord, ])) +
  geom_point(aes(x, y), as_data_frame(fit2$s)) +
  theme_bw() +
  coord_cartesian(xlim = c(-1.6, 1.6), ylim = c(-1.1, 1.1)) +
  labs(x = "x", y = "y", title = "Before")

g1 <- ggplot() +
  geom_point(aes(x, y), x, colour = "darkgray") +
  geom_path(aes(x, y), xadf2) +
  geom_point(aes(x, y), xadf2) +
  theme_bw() +
  coord_cartesian(xlim = c(-1.6, 1.6), ylim = c(-1.1, 1.1)) +
  labs(x = "x", y = "y", title = "After")

ggif_list(list(g0, g1))

## ----relambda2-----------------------------------------------------------
ggplot() +
  geom_segment(aes(x = x$x, xend = fit2$s[,1], y = x$y, yend = fit2$s[,2]), linetype = "dashed", colour = "lightgray") +
  geom_path(aes(x, y), as_data_frame(fit2$s[fit2$ord, ]), colour = "darkgray") +
  geom_point(aes(x, y), x, colour = "darkgray") +
  geom_point(aes(x, y), as_data_frame(fit2$s)) +
  theme_bw() +
  coord_cartesian(xlim = c(-1.6, 1.6), ylim = c(-1.1, 1.1)) +
  labs(x = "x", y = "y")

## ----compare, fig.width=8, fig.height=8----------------------------------
data("benchmarks", package = "princurve")

ggplot(benchmarks, aes(num_points, median / 1000)) +
  geom_point() +
  geom_line() +
  facet_wrap(~expr, ncol = 1, scales = "free") +
  theme_bw() +
  labs(x = "Number of rows in dataset", y = "Time (s)") +
  scale_colour_brewer(palette = "Set1")

