"bias.correct.curve" <- function(x, pcurve, ...)
{
# bias correction, as suggested by
#Jeff Banfield
  p <- ncol(x)
  ones <- rep(1, p)
  sbar <- apply(pcurve$s, 2, "mean")
  ray <- drop(sqrt(((x - pcurve$s)^2) %*% ones))
  dist1 <- (scale(x, sbar, FALSE)^2) %*% ones
  dist2 <- (scale(pcurve$s, sbar, FALSE)^2) %*% ones
  sign <- 2 * as.double(dist1 > dist2) - 1
  ray <- sign * ray
  sray <- approx(periodic.lowess(pcurve$lambda, ray, ...)$x,
		 periodic.lowess(pcurve$lambda, ray, ...)$y,
		 pcurve$lambda)$y
  ## AW: changed periodic.lowess() to periodic.lowess()$x and $y
  pcurve$s <- pcurve$s + (abs(sray)/ray) * ((x - pcurve$s))
  get.lam(x, pcurve$s, pcurve$tag, stretch = 0)
}


"get.lam" <- function(x, s, tag, stretch = 2)
{
  storage.mode(x) <- "double"
  storage.mode(s) <- "double"
  storage.mode(stretch) <- "double"
  if(!missing(tag))
    s <- s[tag,  ]
  np <- dim(x)
  if(length(np) != 2)
    stop("get.lam needs a matrix input")
  n <- np[1]
  p <- np[2]
  tt <- .Fortran("getlam",
		 n,
		 p,
		 x,
		 s = x,
		 lambda = double(n),
		 tag = integer(n),
		 dist = double(n),
		 as.integer(nrow(s)),
		 s,
		 stretch,
		 double(p),
		 double(p),
                 PACKAGE = "princurve")[c("s", "tag", "lambda", "dist")]
  tt$dist <- sum(tt$dist)
  class(tt) <- "principal.curve"
  tt
}


"lines.principal.curve" <- function(x, ...)
  lines(x$s[x$tag,  ], ...)


"periodic.lowess"<- function(x, y, f = 0.59999999999999998, ...)
{
  n <- length(x)
  o <- order(x)
  r <- range(x)
  d <- diff(r)/(length(unique(x)) - 1)
  xr <- x[o[1:(n/2)]] - r[1] + d + r[2]
  xl <- x[o[(n/2):n]] - r[2] - d + r[1]
  yr <- y[o[1:(n/2)]]
  yl <- y[o[(n/2):n]]
  xnew <- c(xl, x, xr)
  ynew <- c(yl, y, yr)
  f <- f/2
  fit <- lowess(xnew, ynew, f = f, ...)
  approx(fit$x, fit$y, x[o], rule = 2)
  # AW: changed fit to fit$x, fit$y
}


"plot.principal.curve" <- function(x, ...)
  plot(x$s[x$tag,  ], ..., type = "l")


"points.principal.curve" <- function(x, ...)
  points(x$s, ...)

"principal.curve" <- 
  function(x, start = NULL, thresh = 0.001, plot.true = FALSE, maxit = 10,
	   stretch = 2, smoother = "smooth.spline", trace = FALSE, ...)
{
  smooth.list <- c("smooth.spline", "lowess", "periodic.lowess")
  smoother <- match.arg(smoother, smooth.list)
  stretches <- c(2, 2, 0)
  if(missing(stretch))
    stretch <- stretches[match(smoother, smooth.list)]
  this.call <- match.call()
  dist.old <- sum(diag(var(x)))
  d <- dim(x)
  n <- d[1]
  p <- d[2] 
# You can give starting values for the curve
  if(missing(start)) {
# use largest principal component
    if(smoother == "periodic.lowess") start <- startCircle(x)
    else {
      xbar <- apply(x, 2, "mean")
      xstar <- scale(x, xbar, FALSE)
      svd.xstar <- svd(xstar)
      dd <- svd.xstar$d
      lambda <- svd.xstar$u[, 1] * dd[1]
      tag <- order(lambda)
      s <- scale(outer(lambda, svd.xstar$v[, 1]),  - xbar, FALSE)
      dist <- sum((dd^2)[-1]) * n
      start <- list(s = s, tag = tag, lambda = lambda, dist
		    = dist)
    }
  }
  else if(!inherits(start, "principal.curve")) {
# use given starting curve 
    if(is.matrix(start))
      start <- get.lam(x, start, stretch = stretch)
    else stop("Invalid starting curve: should be a matrix or principal.curve")
  }
  pcurve <- start
  if(plot.true) {
    plot(x[, 1:2], xlim = adjust.range(x[, 1], 1.3999999999999999),
	 ylim = adjust.range(x[, 2], 1.3999999999999999))
    lines(pcurve$s[pcurve$tag, 1:2])
  }
  it <- 0
  if(trace)
    cat("Starting curve---distance^2: ", pcurve$dist, "\n")
  while(abs((dist.old - pcurve$dist)/dist.old) > thresh & it < maxit) {
    it <- it + 1
    s <- NULL
    for(j in 1:p) {
      sj <- switch(smoother,
		   lowess = lowess(pcurve$lambda, x[, j], ...)$y,
		   smooth.spline = smooth.spline(pcurve$lambda,
		       x[, j], ..., df = 5)$y,
		   periodic.lowess = periodic.lowess(pcurve$lambda,
		       x[, j], ...)$y)
      s <- cbind(s, sj)
    }
    dist.old <- pcurve$dist
    pcurve <- get.lam(x, s, stretch = stretch)
    if(smoother == "periodic.lowess")
      pcurve <- bias.correct.curve(x, pcurve, ...)
    if(plot.true) {
      plot(x[, 1:2], xlim = adjust.range(x[, 1], 1.3999999999999999),
	   ylim = adjust.range(x[, 2], 1.3999999999999999))
      lines(pcurve$s[pcurve$tag, 1:2])
    }
    if(trace)
      cat("Iteration ", it, "---distance^2: ", pcurve$dist, "\n")
  }
  structure(list(s = pcurve$s, tag = pcurve$tag, lambda = pcurve$lambda, 
		 dist = pcurve$dist, call = this.call),
	    class = "principal.curve")
}


adjust.range <- function (x, fact)
  {
# AW: written by AW, replaces ylim.scale
    r <- range (x);
    d <- diff(r)*(fact-1)/2;
    c(r[1]-d, r[2]+d)
  }


"startCircle" <- function(x)
{
# the starting circle uses the first two co-ordinates,
# and assumes the others are zero
  d <- dim(x)
  n <- d[1]
  p <- d[2]	# use best circle centered at xbar
  xbar <- apply(x, 2, "mean")
  ray <- sqrt((scale(x, xbar, FALSE)^2) %*% rep(1, p))
  radius <- mean(ray)
  s <- cbind(radius * sin((pi * (1:101))/50),
	     radius * cos((pi * (1:101))/50))
  if(p > 2)
    s <- cbind(s, matrix(0, n, p - 2))
  get.lam(x, s)
}


"whiskers" <- function(from, to)
	segments(from[, 1], from[, 2], to[, 1], to[, 2])
