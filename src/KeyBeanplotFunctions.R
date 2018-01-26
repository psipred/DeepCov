require(beanplot)
#SMK this now accepts an additional parameter, beanline.length.
beanplot <- function(..., bw = "SJ-dpi", kernel = "gaussian", cut = 3, cutmin = -Inf, cutmax = Inf, grownage = 10, what = c(TRUE,
	TRUE, TRUE, TRUE), add = FALSE, col, axes = TRUE, log = "auto", handlelog = NA, ll = 0.16, wd = NA, maxwidth = 0.8,
	maxstripline = 0.96, method = "stack", names, overallline = "mean", beanlines = overallline, horizontal = FALSE, side = "no",
	jitter = NULL, beanlinewd = 2, frame.plot = axes, border = NULL, innerborder = NA, at = NULL, boxwex = 1, ylim = NULL,
	xlim = NULL, show.names = NA, beanline.length = NA) {
	mdensityxy <- function(x) {
		if (length(x) > 0) {
			from <- max(cutmin, (min(mlog(x)) - cut * bw))
			to <- min(cutmax, max(mlog(x)) + cut * bw)
			density(mlog(x), bw = bw, kernel = kernel, from = from, to = to)[c("x", "y")]
		} else list(x = numeric(), y = numeric())
	}
	args <- match.call()
	mcall <- as.list(args)
	method <- pmatch(method, c("overplot", "stack", "jitter"))
	if (is.na(method) || method == 0)
		stop("invalid plotting method")
	beanlines <- pmatch(beanlines, c("mean", "median", "quantiles"))
	if (is.na(beanlines) || beanlines == 0)
		stop("invalid beanlines")
	overallline <- pmatch(overallline, c("mean", "median"))
	if (is.na(overallline) || overallline == 0)
		stop("invalid overallline")
	side <- pmatch(side, c("no", "first", "second", "both"))
	if (is.na(side) || side == 0)
		stop("invalid side")
	groups <- beanplot:::getgroupsfromarguments(args)
	groups <- lapply(groups, na.omit)
	n <- length(groups)
	displayn <- if (side == 4)
		ceiling(n/2)
	else n
	if (n == 0)
		stop("no data found to beanplot")
	if (missing(names)) {
		if (is.null(base::names(groups)))
			attr(groups, "names") = 1:displayn
		names <- base::names(groups)
	} else {
		attr(groups, "names") <- names
		if (is.na(show.names))
			show.names <- TRUE
	}
	if (is.null(at)) {
		at <- 1:displayn
	}
	if ((side == 4) && (length(names) > length(at))) {
		for (i in 1:length(at)) {
			names[i] <- beanplot:::makecombinedname(names[i * 2 - 1], names[i * 2])
		}
		length(names) <- length(at)
	}
	combinedpolygons <- ((side == 4) && (length(border) < 2) && (n%%2 == 0))
	if (missing(col))
		col <- par("fg")
	if (!is.list(col))
		col <- list(col)
	else combinedpolygons <- FALSE
	col <- lapply(col, beanplot:::fixcolorvector)
	col <- rep(col, length.out = n)
	if (!is.null(border))
		border <- rep(border, length.out = n)
	if (!add && log == "auto") {
		if (beanplot:::seemslog(groups)) {
			log <- "y"
			message("log=\"y\" selected")
		} else log <- ""
	}
	if (is.na(handlelog))
		if (add && ((horizontal & par()$xlog) || (!horizontal & par()$ylog)))
			handlelog <- TRUE
		else if (!add && (log != ""))
			handlelog <- TRUE
		else handlelog <- FALSE
	if (handlelog) {
		mlog <- base::log
		mexp <- base::exp
	} else {
		mlog <- function(x) {
			x
		}
		mexp <- mlog
	}
	if (!is.numeric(bw)) {
		bw <- mean(sapply(groups, function(x) {
			ifelse(length(x) > 1, density(mlog(x), kernel = kernel, bw = bw)$bw, NA)
		}), na.rm = TRUE)
		if (is.nan(bw))
			bw <- 0.5
	}
	dens <- sapply(groups, mdensityxy)
	for (i in 1:n) dens[["y", i]] <- dens[["y", i]] * min(1, length(groups[[i]])/grownage)
	if (is.na(wd))
		wd <- maxwidth/max(unlist(dens["y", ]))
	wd2 <- wd * boxwex/2
	axespars <- lapply(mcall[base::names(mcall) %in% c("xaxt", "yaxt", "las", "cex.axis", "col.axis", "format", "tick",
		"xaxp", "yaxp")], eval, parent.frame())
	if (!add) {
		if (!is.numeric(xlim)) {
			if (side == 2)
				xlim <- c(0, displayn)
			else if (side == 3)
				xlim <- c(1, displayn + 1)
			else xlim <- c(0.5, displayn + 0.5)
		}
		if (!is.numeric(ylim))
			ylim <- range(groups, mexp(unlist(dens["x", ])))
		plot.new()
		windowpars <- lapply(mcall[base::names(mcall) %in% c("yaxs", "xaxs")], eval)
		if (horizontal) {
			names(windowpars)[names(windowpars) %in% c("xaxs", "yaxs")] <- rev(names(windowpars)[names(windowpars) %in%
				c("xaxs", "yaxs")])
			if (log == "y")
				log <- "x"
			do.call("plot.window", c(list(xlim = ylim, ylim = xlim, log = log), windowpars))
		} else {
			do.call("plot.window", c(list(xlim = xlim, ylim = ylim, log = log), windowpars))
		}
		if (frame.plot)
			box()
		if (axes)
			do.call("axis", c(list(side = 2 - horizontal), axespars))
	}
	if (axes) {
		if (is.na(show.names))
			show.names <- (n > 1)
		if (show.names)
			do.call("axis", c(list(1 + horizontal, at = at, labels = names), axespars))
	}
	if (what[1]) {
		if (overallline == 2)
			val <- mexp(median(mlog(unlist(groups))))
		else val <- mexp(mean(mlog(unlist(groups))))
		if (horizontal)
			abline(v = val, lty = 3)
		else abline(h = val, lty = 3)
	} else {
		val = "value not calculated, overall line was omitted"
	}
	if (what[2]) {
		# SMK don't use my version
		beanplot:::beanplotpolyshapes(side, dens, at, wd2, combinedpolygons, displayn, n, col, border, horizontal, mlog, mexp)
	}
	if (what[3]) {
		# SMK use my version
		stats = beanplotbeanlines(groups, side, beanlines, beanlinewd, at, boxwex, n, col, horizontal, mlog, mexp, beanline.length)
	} else {
		stats = "not calculated, beanlines were omitted"
	}
	if (what[4]) {
		beanplot:::beanplotscatters(groups, side, method, jitter, dens, at, wd2, boxwex, n, ll, maxstripline, col, horizontal, mlog,
			mexp)
	}
	if (any(!is.na(innerborder))) {
		beanplot:::beanplotinnerborders(innerborder, at, dens, side, displayn, n, horizontal, mexp)
	}
	titlepars <- lapply(mcall[base::names(mcall) %in% c("main", "sub", "xlab", "ylab", "cex.main", "col.main", "cex.lab",
		"col.lab", "cex.sub", "col.sub")], eval, parent.frame())
	do.call("title", titlepars)
	invisible(list(bw = bw, wd = wd, names = names, stats = stats, overall = val, log = log))
}


# beanplot:::beanplotpolyshapes <- function(side, dens, at, wd2, combinedpolygons, displayn, n, col, border, horizontal,
	# mlog, mexp) {
	# if ((side == 4) && (!combinedpolygons))
		# at <- rep(at, each = 2)
	# pborder <- NA
	# for (i in 1:ifelse(combinedpolygons, displayn, n)) {
		# if (combinedpolygons) {
			# x1 <- rev(dens[["y", i * 2]] * wd2) + at[i]
			# x2 <- dens[["y", i * 2 - 1]] * -wd2 + at[i]
			# y1 <- rev(dens[["x", i * 2]])
			# y2 <- dens[["x", i * 2 - 1]]
		# } else {
			# if (side == 2 || ((side == 4) && (i%%2 == 1))) {
				# x1 <- NULL
				# y1 <- NULL
			# } else {
				# x1 <- dens[["y", i]] * wd2 + at[i]
				# y1 <- dens[["x", i]]
			# }
			# if (side == 3 || ((side == 4) && (i%%2 == 0))) {
				# x2 <- NULL
				# y2 <- NULL
			# } else {
				# x2 <- rev(dens[["y", i]]) * -wd2 + at[i]
				# y2 <- rev(dens[["x", i]])
			# }
		# }
		# if (length(x1) > 0) {
			# x1 <- c(at[i], x1, at[i])
			# y1 <- c(y1[1], y1, y1[length(y1)])
		# }
		# if (length(x2) > 0) {
			# x2 <- c(at[i], x2, at[i])
			# y2 <- c(y2[1], y2, y2[length(y2)])
		# }
		# if (combinedpolygons || (side == 1))
			# pborder <- border[[i]]
		# if (horizontal) {
			# polygon(mexp(c(y1, y2)), c(x1, x2), col = col[[i]][1], border = pborder)
			# if ((!combinedpolygons) && (side != 1)) {
				# if (is.null(border[[i]]))
					# lines(mexp(c(y1, y2)), c(x1, x2))
				# else if (!is.na(border[[i]]))
					# lines(mexp(c(y1, y2)), c(x1, x2), col = border[[i]])
			# }
		# } else {
			# polygon(c(x1, x2), mexp(c(y1, y2)), col = col[[i]][1], border = pborder)
			# if ((!combinedpolygons) && (side != 1)) {
				# if (is.null(border[[i]]))
					# lines(c(x1, x2), mexp(c(y1, y2)))
				# else if (!is.na(border[[i]]))
					# lines(c(x1, x2), mexp(c(y1, y2)), col = border[[i]])
			# }
		# }
	# }
# }

# SMK MODIFIED
beanplotbeanlines <-
function (groups, side, beanlines, beanlinewd, at, boxwex, n,
    col, horizontal, mlog, mexp, beanline.length)
{
	if (is.na(beanline.length))
		beanline.length = boxwex
    res = NULL
    if (side == 4)
        at <- rep(at, each = 2)
    for (i in 1:n) {
        if (beanlines == 2) {
            quants <- mexp(median(mlog(groups[[i]])))
        }
        else if (beanlines == 3) {
            quants <- mexp(quantile(mlog(groups[[i]])))
        }
        else {
            quants <- mexp(mean(mlog(groups[[i]])))
        }
        res = c(res, quants)
        x1 <- if (side == 3 || ((side == 4) && (i%%2 == 0)))
            at[i]
        else at[i] - beanline.length/2
        x2 <- if (side == 2 || ((side == 4) && (i%%2 == 1)))
            at[i]
        else at[i] + beanline.length/2
        if (horizontal)
            segments(quants, x1, quants, x2, lwd = beanlinewd,
                col = col[[i]][4])
        else segments(x1, quants, x2, quants, lwd = beanlinewd,
            col = col[[i]][4])
    }
    res
}
