control.presmooth <- 
function(n.boot = 500, q.weight = c(0.2, 0.8), length.grid.bw.plugin = 100, length.grid.ise = 100, pilot.par.ini= NULL, save.data = FALSE, save.mise = FALSE, na.action = na.omit){
	list(n.boot = n.boot, q.weight = q.weight, length.grid.bw.plugin = length.grid.bw.plugin, length.grid.ise = length.grid.ise, pilot.par.ini = pilot.par.ini, save.data = save.data, save.mise = save.mise, na.action = na.action)
}

presmooth <-
function(times, status, dataset = NULL, estimand = c("S", "H", "f", "h"), bw.selec = c("plug-in", "bootstrap", "fixed"), fixed.bw = NULL, grid.bw.pil = NULL, grid.bw = NULL, kernel = c("biweight", "triweight"), bound = c("none", "left", "right", "both"), x.est = NULL, control = control.presmooth()){ 
    na.action <- control$na.action
	# Define the 'internal' data frame
	dfr <-
    	if (is.null(dataset)) na.action(data.frame(times, status))
 		else na.action(dataset[, c(deparse(substitute(times)), deparse(substitute(status)))])
	names(dfr) <- c("t", "d")
	dfr$t <- as.numeric(dfr$t)
	dfr$d <- as.integer(dfr$d)
	dfr <- dfr[order(dfr$t, 1 - dfr$d),]
    if(length(kernel) > 1) kernel <- kernel[1]
  	# Define numeric codes for: kernel, estimated function, boundary correction method and bandwith selection method used
	n.kernel <- pmatch(kernel, c("biweight","triweight"))
    if(is.na(n.kernel)) stop("'kernel' argument must be one of 'biweight', 'triweight'")
    # Constants: c1: second moment of the kernel; c2: integral of the squared kernel, c_K; c3: integral of the squared second derivative of the kernel; c4 
    const <- switch(n.kernel, list(c1 = 1/7, c2 = 5/7, c3 = 45/2, c4 = 15/7, c5 = 25/231), list(c1 = 1/9, c2 = 350/429, c3 = 35, c4 = 35/11, c5 = 245/2574))
    if(length(estimand) > 1) estimand <- estimand[1]
    n.estimand <- match(estimand, c("S", "H", "f", "h"))
    if(is.na(n.estimand)) stop("'estimand' argument must be one of 'S', 'H', 'f', 'h'")
    if(length(bound) > 1) bound <- bound[1]
	n.bound <- pmatch(bound, c("none", "left", "right", "both"))
    if(is.na(n.bound)) stop("'bound' argument must be one of 'none', 'left', 'right', 'both'")
    if(length(bw.selec) > 1) bw.selec <- bw.selec[1]
    n.bw.selec <- pmatch(bw.selec, c("plug-in", "bootstrap", "fixed"))
    if(is.na(n.bw.selec)) stop("'bw.selec' argument must be one of 'plug-in', 'bootstrap', 'fixed'")
	n <- dim(dfr)[[1]]
    if(n.estimand <= 2){
		if(!is.null(x.est)) warning("'x.est' values overriden: only estimates at observed times are given")
		x.est <- dfr$t
	}
	else if(is.null(x.est)) x.est <- seq(min(dfr$t), quantile(dfr$t, 0.9), length.out = 50)
    le.x.est <- length(x.est)
	q.w <- quantile(dfr$t, control$q.weight)
	cond <- (dfr$t > q.w[1]) & (dfr$t < q.w[2])
    range.t <- diff(range(dfr$t))
	forh <- n.estimand >= 3
 	if(((n.bw.selec == 1) & forh) | (n.bw.selec == 2)){ 
	    # Default grid where the presmoothing pilot bandwidth will be selected
    	if(is.null(grid.bw.pil)) grid.bw.pil <- q.w[2]/10 * 10^seq(0, 1, by = 0.01)
     	le.bw.pilot <- length(grid.bw.pil)
        # Presmoothing pilot bandwidth, selected by cross-validation
        ls.cv <- .C("lscv", dfr$t, dfr$d, n, grid.bw.pil, le.bw.pilot, n.kernel, lscv = numeric(le.bw.pilot), PACKAGE = "survPresmooth")$lscv
        pilot.bw1 <- grid.bw.pil[which.min(ls.cv)]
	}
	if(((n.bw.selec == 2) & forh) | (n.bw.selec == 1)){
		# Initial parameter estimates of a simple Weibull model for the observed lifetimes, based on the method of moments
		shape0 <- uniroot(function(x) var(dfr$t) - (mean(dfr$t)/ gamma(1+1/x))^2*(gamma(1+2/x) - (gamma(1+1/x))^2), interval = c(0.02, 20))$root
		scale0 <- mean(dfr$t)/gamma(1+1/shape0)
		pilot.par.ini <- 
			if(is.null(control$pilot.par.ini)) c(1/shape0, scale0, shape0, scale0, shape0 , 1/scale0, 1/3, 1/3)
		   	else control$pilot.par.ini
		pilot.bw2.integrand <- function(x, par, modif = as.integer(0)) .C("pilot2forhintegrand", x, length(x), par, length(par), modif, pilot2hint = numeric(length(x)), PACKAGE = "survPresmooth")$pilot2hint
	}
	# Bandwidth computation
	switch(n.bw.selec,
		{ # Plug-in
			le.plugin <- as.integer(ifelse(control$length.grid.bw.plugin %% 2 == 0, control$length.grid.bw.plugin - 1, control$length.grid.bw.plugin))
			grid.plugin <- seq(q.w[1], q.w[2], length.out = le.plugin)
			step.plugin <- (grid.plugin[le.plugin] - grid.plugin[1])/(le.plugin-1)
			esf <- 1 - ecdf(dfr$t)(grid.plugin) + 1/n
			switch(estimand,
				S =,
				H = {
					# Pilot
					# A logistic model is supposed for the probability of non-censoring
				    logistfit <- glm(d ~ t ,family = binomial, data = dfr)
				  	p <- as.numeric(predict(logistfit, newdata = data.frame(t = grid.plugin), type = "response"))
					# Fit a simple Weibull model for observed lifetimes
					llikplugH <- function(p, x) -sum(dweibull(x, p[1], p[2], log = TRUE))
					estparplug <- optim(c(shape0, scale0), llikplugH, method = "L-BFGS-B", lower = 1e-8, x = dfr$t)$par
					# Computation of C1 and C2 according to formulas (2.154) and (2.155) on page 87 of López de Ullibarri's thesis
					C1fun1 <- function(x, t, n, coef, par) .C("c1integrand1", x, length(x), t, n, coef, par, c1fun1 = numeric(length(x)), PACKAGE = "survPresmooth")$c1fun1
					C1fun2 <- function(x, t, n, coef, par) .C("c1integrand2", x, length(x), t, n, coef, par, c1fun2 = numeric(length(x)), PACKAGE = "survPresmooth")$c1fun2
					C1integrand1 <- C1integrand2 <- numeric(le.plugin)
					for(i in 1:le.plugin){
						C1integrand1[i] <- integrate(C1fun1, lower = if(i==1) 0 else grid.plugin[i-1], upper = grid.plugin[i], t = dfr$t, n = n, coef = logistfit$coefficients, par = estparplug, subdivisions = 10000, stop.on.error = FALSE)$value
						C1integrand2[i] <- integrate(C1fun2, lower = if(i==1) 0 else grid.plugin[i-1], upper = grid.plugin[i], t = dfr$t, n = n, coef = logistfit$coefficients, par = estparplug, subdivisions = 10000, stop.on.error = FALSE)$value
					}
					C1integrand <- cumsum(C1integrand1) * cumsum(C1integrand2)
					C1 <- const$c1 / 2 * .C("simpson", C1integrand, le.plugin, step.plugin, integral = numeric(1), PACKAGE = "survPresmooth")$integral
					C2pre <- p * (1 - p) * dweibull(grid.plugin, estparplug[1], estparplug[2]) / esf^2
	           		C2 <- const$c4 / 4 * .C("simpson", C2pre, le.plugin, step.plugin, integral = numeric(1), PACKAGE = "survPresmooth")$integral
                    # Formula on page 91 of López-de-Ullibarri's thesis, modified by considering an empirical upper bound (trying to avoid the sporadical ocurrence of exceedingly large bandwidths)
                  	pilot.bw1 <- min((C2/C1 * ifelse(C1 < 0, -1 , 3/2) / n)^(1/5), 2 * range.t)
                    # Computation of D1 and D2 according to formulas (2.237) and (2.238) on page 125 of López de Ullibarri's thesis
					Dintegrand <- .C("dintegrand", grid.plugin, le.plugin, t, n, logistfit$coefficients, estparplug, p, d1int = numeric(le.plugin), d2int = numeric(le.plugin), PACKAGE = "survPresmooth")[c("d1int", "d2int")]
					D1 <- const$c1 * .C("simpson", Dintegrand$d1int, le.plugin, step.plugin, integral = numeric(1), PACKAGE = "survPresmooth")$integral
					D2 <- -const$c2 * .C("simpson", Dintegrand$d2int, le.plugin, step.plugin, integral = numeric(1), PACKAGE = "survPresmooth")$integral
                   	# Formula on page 127 of López-de-Ullibarri's thesis, modified by considering an empirical upper bound (as above, trying to avoid the sporadical ocurrence of too large bandwidths)
                  	pilot.bw2 <- min((D2/D1 * ifelse(D1 < 0, 1/2 , -1) / n)^(1/3), 2 * range.t)
					pilot.bw <- c(pilot.bw1, pilot.bw2)
					# Bandwidth (computed according to formulas on page 17 of López-de-Ullibarri's thesis)
   				    p <- .C("nadarayawatson", dfr$t, n, dfr$t, dfr$d, n, pilot.bw[2], n.kernel, phat = numeric(n), PACKAGE = "survPresmooth")$phat
                   	esf2 <- 1 - ecdf(dfr$t)(dfr$t) + 1/n
				    Q <- sum(p[cond] * (1 - p[cond]) / esf2[cond]^2) / n
					funalpha1 <- function(x, t, d, n, bw, kern) .C("alphaintegrand", t, d, n, x, length(x), bw, bw, kern, alphaint = numeric(length(x)), PACKAGE = "survPresmooth")$alphaint
                   	prealpha <- numeric(le.plugin)
                   	for(i in 1:le.plugin)
						prealpha[i] <- integrate(funalpha1, lower = if(i == 1) 0 else grid.plugin[i-1], upper = grid.plugin[i], t = dfr$t, d = dfr$d, n = n, bw = pilot.bw[1], kern = n.kernel, subdivisions = 10000, stop.on.error = FALSE)$value
                 	alpha <- cumsum(prealpha)
					A <- .C("simpson", alpha^2, le.plugin, step.plugin, integral = numeric(1), PACKAGE = "survPresmooth")$integral
					# The same modification than for pilot bandwidths is incorporated
					bw <- min((const$c5 * Q / A / 2 / const$c1^2 / n)^(1/3), 2 * range.t)
				}, 
				f = {
					# Pilot
					# suppressWarnings() is necessary to avoid some annoying but irrelevant warnings...
                    llikplugf <- function(p, x, d) -sum(suppressWarnings(log(d * (p[7] * dweibull(x, p[1], p[2]) + p[8] * dweibull(x, p[3], p[4]) + (1-p[7]-p[8]) * dweibull(x, p[5], p[6])) + (1-d) * (p[7] * pweibull(x, p[1], p[2], lower.tail = FALSE) + p[8] * pweibull(x, p[3], p[4], lower.tail = FALSE) +(1-p[7]-p[8]) * pweibull(x, p[5], p[6], lower.tail = FALSE)))))
                   	# Redefinition of llikplugf for a mixture of only two components
                   	llikplugf2 <- function(p, x, d) -sum(suppressWarnings(log(d * (p[5] * dweibull(x, p[1], p[2]) + (1 - p[5]) * dweibull(x, p[3], p[4])) + (1 - d) * (p[5] * pweibull(x, p[1], p[2], lower.tail = FALSE) + (1-p[5]) * pweibull(x, p[3], p[4], lower.tail = FALSE)))))
                   	# Redefinition of llikplugf according to a simple Weibull model
                   	llikplugf3 <- function(p, x, d) -sum(suppressWarnings(log(d * dweibull(x, p[1], p[2]) + (1 - d) * pweibull(x, p[1], p[2], lower.tail = FALSE))))
					estparplug <- constrOptim(pilot.par.ini, llikplugf, grad = NULL, x = dfr$t, d= rep(1, n), ui = rbind(c(rep(0, 6), -1, -1), c(1, rep(0,  7)), c(0, 1, rep(0, 6)), c(rep(0, 2), 1, rep(0, 5)), c(rep(0, 3), 1, rep(0, 4)), c(rep(0, 4), 1, rep(0, 3)), c(rep(0, 5), 1, 0, 0), c(rep(0, 6), 1, 0), c(rep(0, 7), 1), c(rep(0, 6), -1, 0), c(rep(0, 7), -1)), ci = c(-1, rep(0, 8), -1, -1))$par
			    	# If one component of the mixture is negligible, refit with a mixture of only 2 Weibulls
					if(any(estparplug[7] < 0.1, estparplug[8] < 0.1, 1 - estparplug[7] - estparplug[8] < 0.1))
                    	estparplug <- constrOptim(c(pilot.par.ini[1:4], 0.5), llikplugf2, grad = NULL, x = dfr$t, d = rep(1, n), ui = rbind(c(1, rep(0,  4)), c(0, 1, rep(0, 3)), c(rep(0, 2), 1, rep(0, 2)), c(rep(0, 3), 1, 0), c(rep(0, 4), 1), c(rep(0, 4), -1)), ci = c(rep(0, 5), -1))$par              
			   		# Again, if one of the 2-component-mixture is negligible, refit with only 1 Weibull
					if(estparplug[5] < 0.1)
                    	estparplug <- constrOptim(c(pilot.par.ini[3:4]), llikplugf3, grad = NULL, x = dfr$t, d = rep(1, n), ui = rbind(c(1, 0), c(0, 1)), ci = rep(0, 2))$par
					# Formula (4.76) for b_1 on page 236 of A. Jácome's thesis
					pilot.bw2 <- (const$c3 / const$c1 / n / integrate(pilot.bw2.integrand, q.w[1], q.w[2], par = estparplug)$value)^(1/7)
					estparplug2 <- constrOptim(pilot.par.ini, llikplugf, grad = NULL, x = dfr$t,d = dfr$d, ui = rbind(c(rep(0, 6), -1, -1), c(1, rep(0,  7)), c(0, 1, rep(0, 6)), c(rep(0, 2), 1, rep(0, 5)), c(rep(0, 3), 1, rep(0, 4)), c(rep(0, 4), 1, rep(0, 3)), c(rep(0, 5), 1, 0, 0), c(rep(0, 6), 1, 0), c(rep(0, 7), 1), c(rep(0, 6), -1, 0), c(rep(0, 7), -1)), ci = c(-1, rep(0, 8), -1, -1))$par
					if(any(estparplug2[7] < 0.1, estparplug2[8] < 0.1, 1 - estparplug2[7] - estparplug2[8] < 0.1))
                   	estparplug2 <- constrOptim(c(pilot.par.ini[1:4], 0.5), llikplugf2, grad = NULL, x = dfr$t, d = dfr$d, ui = rbind(c(1, rep(0,  4)), c(0, 1, rep(0, 3)), c(rep(0, 2), 1, rep(0, 2)), c(rep(0, 3), 1, 0), c(rep(0, 4), 1), c(rep(0, 4), -1)), ci = c(rep(0, 5), -1))$par
					if(estparplug2[5] < 0.1)
                    	estparplug2 <- constrOptim(c(pilot.par.ini[3:4]), llikplugf3, grad = NULL, x = dfr$t, d = dfr$d, ui = rbind(c(1, 0), c(0, 1)), ci = rep(0, 2))$par
					intf3sq <- integrate(pilot.bw2.integrand, q.w[1], q.w[2], par = estparplug2)$value
					# If f is replaced by h*p*(1-F)/(1-H), 3rd pilot bandwidth in A. Jácome's thesis is unnecessary
                    km <- .C("presmestim", dfr$t, n, dfr$t, n, NULL, NULL, NULL, as.numeric(dfr$d), as.integer(1), pest = numeric(n), PACKAGE = "survPresmooth")$pest
                    fact <- sum((km[cond] / (1 - ecdf(dfr$t)(sort(dfr$t)[cond])))^2 * (dfr$d[order(dfr$t)][cond] == 1))/n
         		   	pilot.bw3 <- (const$c3 * fact / const$c1 / n / intf3sq)^(1/7)
                    pilot.bw <- c(pilot.bw1, pilot.bw2, pilot.bw3)
					# Bandwidth
					# Computation of the 5 integrals in psi(L)
					funalpha2 <- function(x, t, d, n, bw, kern) .C("alphaintegrand", t, d, n, x, length(x), bw[1], bw[2], kern, alphaint = numeric(length(x)), PACKAGE = "survPresmooth")$alphaint
					prealpha <- numeric(le.plugin)
					for(i in 1:le.plugin)
    		        	prealpha[i] <- integrate(funalpha2, lower = if(i == 1) 0 else grid.plugin[i - 1], upper = grid.plugin[i], t = dfr$t, d = dfr$d, n = n, bw = pilot.bw[1:2], kern = n.kernel, subdivisions = 10000, stop.on.error = FALSE)$value
					alpha <- cumsum(prealpha)
					# In termsmise.c f replaced by h*p*(1-F)/(1-H)
                   	ints1to5 <- .C("termsmise", dfr$t, dfr$d, n, esf, grid.plugin, le.plugin, step.plugin, pilot.bw, n.kernel, n.estimand, alpha, int1 = numeric(1), int2 = numeric(1), int3 = numeric(1), int4 = numeric(1), int5 = numeric(1), PACKAGE = "survPresmooth")[paste("int", 1:5, sep = "")]
                   	# The objective function is formula (3.180) for AMISE in López-de-Ullibarri's thesis
					# When minimizing, an empirical upper bound is considered, trying to avoid the sporadical ocurrence of exceedingly large bandwidths
                    bw <- pmin(constrOptim(pilot.bw[-3], function(pars, n, kern, c1, c2, n.est, ints) .C("funplugin", pars[1], pars[2], n, kern, c1, c2, n.est, ints$int1, ints$int2, ints$int3, ints$int4, ints$int5, fplug = numeric(1), PACKAGE = "survPresmooth")$fplug, grad = NULL, n = n, kern = n.kernel, c1 = const$c1, c2 = const$c2, n.est = n.estimand, ints = ints1to5, ui = rbind(c(1, 0), c(0, 1)), ci = c(0, 0))$par, 2 * range.t)
             	},
				h = {
					# Log-likelihood function (a mixture of three Weibulls is considered as model for the distribution of lifetimes)
					# suppressWarnings() is necessary to avoid some annoying but irrelevant warnings...
					llikplugh <- function(p, x) -sum(suppressWarnings(log(p[7] * dweibull(x, p[1], p[2]) + p[8] * dweibull(x, p[3], p[4]) + (1-p[7]-p[8]) * dweibull(x, p[5], p[6]))))
					# Optimization with constrOptim() allows to take into account the constraint for the mixture weights (p[1] + p[2] <= 1)
					estparplug <- constrOptim(pilot.par.ini, llikplugh, grad = NULL, x = dfr$t, ui = rbind(c(rep(0, 6), -1, -1), c(1, rep(0,  7)), c(0, 1, rep(0, 6)), c(rep(0, 2), 1, rep(0, 5)), c(rep(0, 3), 1, rep(0, 4)), c(rep(0, 4), 1, rep(0, 3)), c(rep(0, 5), 1, 0, 0), c(rep(0, 6), 1, 0), c(rep(0, 7), 1), c(rep(0, 6), -1, 0), c(rep(0, 7), -1)), ci = c(-1, rep(0, 8), -1, -1))$par
			       	# If one component of the mixture is negligible, refit with a mixture of only 2 Weibulls
					if(any(estparplug[7] < 0.1, estparplug[8] < 0.1, 1 - estparplug[7] - estparplug[8] < 0.1)){
    	            	# Redefinition of llikplugh for a mixture of only two components
						llikplugh <- function(p, x) -sum(suppressWarnings(log(p[5] * dweibull(x, p[1], p[2]) + (1 - p[5]) * dweibull(x, p[3], p[4]))))
                       	estparplug <- constrOptim(c(pilot.par.ini[1:4], 0.5), llikplugh, grad = NULL, x = dfr$t, ui = rbind(c(1, rep(0,  4)), c(0, 1, rep(0, 3)), c(rep(0, 2), 1, rep(0, 2)), c(rep(0, 3), 1, 0), c(rep(0, 4), 1), c(rep(0, 4), -1)), ci = c(rep(0, 5), -1))$par   
					}
			        # Again, if one of the 2-component-mixture is negligible, refit with only 1 Weibull
					if(estparplug[5] < 0.1){
					# Redefinition of llikplugh with a simple Weibull model
                    	llikplugh <- function(p, x) -sum(suppressWarnings(log(dweibull(x, p[1], p[2]))))
                        estparplug <- constrOptim(c(pilot.par.ini[3:4]), llikplugh, grad = NULL, x = dfr$t, ui = rbind(c(1, 0), c(0, 1)), ci = rep(0, 2))$par
					}
					# Formula (4.76) for b_1 on page 236 of A. Jácome's thesis
					pilot.bw2 <- (const$c3 / const$c1 / n / integrate(pilot.bw2.integrand, q.w[1], q.w[2], par = estparplug)$value)^(1/7)
                   	pilot.bw <- c(pilot.bw1, pilot.bw2)
					# Bandwidth
					ints1to5 <- .C("termsmise", dfr$t, dfr$d, n, esf, grid.plugin, le.plugin, step.plugin, pilot.bw, n.kernel, n.estimand, NULL, int1 = numeric(1), int2 = numeric(1), int3 = numeric(1), int4 = numeric(1), int5 = numeric(1), PACKAGE = "survPresmooth")[paste("int", 1:5, sep = "")]
                   	# The objective function isformula (3.180) for AMISE in López-de-Ullibarri's thesis
					# When minimizing, an empirical upper bound is considered, trying to avoid the sporadical ocurrence of exceedingly large bandwidths
                    bw <- pmin(constrOptim(pilot.bw[-3], function(pars, n, kern, c1, c2, n.est, ints) .C("funplugin", pars[1], pars[2], n, kern, c1, c2, n.est, ints$int1, ints$int2, ints$int3, ints$int4, ints$int5, fplug = numeric(1), PACKAGE = "survPresmooth")$fplug, grad = NULL, n = n, kern = n.kernel, c1 = const$c1, c2 = const$c2, n.est = n.estimand, ints = ints1to5, ui = rbind(c(1, 0), c(0, 1)), ci = c(0, 0))$par, 2 * range.t)
				}
			) # End of switch(estimand)
		}, # End of plugin
		{# Bootstrap
			# Pilot
     		# Smoothing pilot bandwidth (only if estimand is "f" or "h")
     		pilot.bw2 <- if(forh){ # if 1
       			# A mixture of three Weibulls is considered as model for the distribution of lifetimes
				llikboot <- function(p, x, d) -sum(suppressWarnings(log(d * (p[7] * dweibull(x, p[1], p[2]) + p[8] * dweibull(x, p[3], p[4]) + (1-p[7]-p[8]) * dweibull(x, p[5], p[6])) + (1-d) * (p[7] * pweibull(x, p[1], p[2], lower.tail = FALSE) + p[8] * pweibull(x, p[3], p[4], lower.tail = FALSE) +(1-p[7]-p[8]) * pweibull(x, p[5], p[6], lower.tail = FALSE)))))
				# Optimization with constrOptim() allows to take into account the constraint for the mixture weights (p[1] + p[2] <= 1)
				estparboot <- constrOptim(pilot.par.ini, llikboot, grad = NULL, x = dfr$t,d = dfr$d, ui = rbind(c(rep(0, 6), -1, -1), c(1, rep(0,  7)), c(0, 1, rep(0, 6)), c(rep(0, 2), 1, rep(0, 5)), c(rep(0, 3), 1, rep(0, 4)), c(rep(0, 4), 1, rep(0, 3)), c(rep(0, 5), 1, 0, 0), c(rep(0, 6), 1, 0), c(rep(0, 7), 1), c(rep(0, 6), -1, 0), c(rep(0, 7), -1)), ci = c(-1, rep(0, 8), -1, -1))$par
		        # If one component of the mixture is negligible, refit with a mixture of only 2 Weibulls
				if(any(estparboot[7] < 0.1, estparboot[8] < 0.1, 1 - estparboot[7] - estparboot[8] < 0.1)){
					llikboot <- function(p, x, d) -sum(suppressWarnings(log(d * (p[5] * dweibull(x, p[1], p[2]) + (1 - p[5]) * dweibull(x, p[3], p[4])) + (1 - d) * (p[5] * pweibull(x, p[1], p[2], lower.tail = FALSE) + (1 - p[5]) * pweibull(x, p[3], p[4], lower.tail = FALSE)))))
					estparboot <- constrOptim(c(pilot.par.ini[1:4], 0.5), llikboot, grad = NULL, x = dfr$t, d = dfr$d, ui = rbind(c(1, rep(0,  4)), c(0, 1, rep(0, 3)), c(rep(0, 2), 1, rep(0, 2)), c(rep(0, 3), 1, 0), c(rep(0, 4), 1), c(rep(0, 4), -1)), ci = c(rep(0, 5), -1))$par
				}
			    # Again, if one of the 2-component-mixture is negligible, refit with only 1 Weibull
				if(estparboot[5] < 0.1){
                	llikboot <- function(p, x, d) -sum(suppressWarnings(log(d * dweibull(x, p[1], p[2]) + (1 - d) * pweibull(x, p[1], p[2], lower.tail = FALSE))))
                    estparboot <- constrOptim(c(pilot.par.ini[3:4]), llikboot, grad = NULL, x = dfr$t, d = dfr$d, ui = rbind(c(1, 0), c(0, 1)), ci = rep(0, 2))$par
				}
				#  f: optimal bandwidth for estimating f without presmoothing, given by Sellero et al. (1999);h: formula for g_hat at the bottom of page 235 in A. Jácome's thesis
                km <- .C("presmestim", dfr$t, n, dfr$t, n, NULL, NULL, NULL, as.numeric(dfr$d), as.integer(1), pest = numeric(n), PACKAGE = "survPresmooth")$pest
                fact <- sum((km[cond] / (1 - ecdf(dfr$t)(sort(dfr$t)[cond])))^2 * (dfr$d[order(dfr$t)][cond] == 1))/n
                intden <- integrate(pilot.bw2.integrand, q.w[1], q.w[2], par = estparboot, modif = as.integer(1))$value
				ifelse (n.estimand == 3, (const$c2 * fact / const$c1^2 / n /intden)^0.2, (const$c3 * fact / const$c1 / n / intden)^(1/7))
     		} # End of if 1
         	pilot.bw <- c(pilot.bw1, pilot.bw2)
			# Bandwith
        	n.boot <- as.integer(control$n.boot)
			le.ise <- as.integer(control$length.grid.ise)
			grid.ise <- seq(q.w[1], q.w[2], length.out = le.ise)
			step.ise <- (grid.ise[le.ise] - grid.ise[1])/(le.ise-1)
			if(is.null(grid.bw)){
				grid.bw.1 <- median(dfr$t)/10 * 10^seq(0, 1, by = 0.05)
				grid.bw.2 <- if(forh) grid.bw.1
        	}
        	else {
            	grid.bw.1 <- 
					if(is.list(grid.bw)) grid.bw[[1]]
					else if(is.numeric(grid.bw)) grid.bw
          		grid.bw.2 <- if(forh) grid.bw[[2]]
			}
			le.bw.1 <- length(grid.bw.1)
			le.bw.2 <- length(grid.bw.2)
        	# Estimates (with the pilot bandwidth) of:
			# a) conditional probability of non-censoring
        	p.hat <- .C("nadarayawatson", dfr$t, n, dfr$t, dfr$d, n, pilot.bw[1], n.kernel, phat = numeric(n), PACKAGE = "survPresmooth")$phat
			# b) survival, density, hazard or cumulative hazard
            ps.dfr <- switch(estimand,
				S =,
				H =	.C("presmestim", grid.ise, le.ise, dfr$t, n, NULL, NULL, NULL, p.hat, n.estimand, pest = numeric(le.ise), PACKAGE = "survPresmooth")$pest,
				f = .C("presmdensfast", grid.ise, le.ise, dfr$t, n, pilot.bw[2], n.kernel, p.hat, dhat = numeric(le.ise), PACKAGE = "survPresmooth")$dhat,
				h = .C("presmtwfast", grid.ise, le.ise, dfr$t, n, pilot.bw[2], n.kernel, p.hat, hhat = numeric(le.ise), PACKAGE = "survPresmooth")$hhat
			)
			# Compute bootstrapped presmoothed estimates
			mise <- .C("misevect", dfr$t, dfr$d, n, n.boot, grid.ise, le.ise, grid.bw.1, le.bw.1, grid.bw.2, le.bw.2, n.kernel, n.estimand, p.hat, ps.dfr, mise = numeric(if(forh) le.bw.1*le.bw.2 else le.bw.1), PACKAGE = "survPresmooth")$mise
			# Return the bandwidth minimizing the mise
			pos.min.mise <- which(mise == min(mise))
        	bw <- 
				pmin(if(forh) c(grid.bw.1[(pos.min.mise - 1)%%le.bw.1 + 1], grid.bw.2[(pos.min.mise - 1)%/%le.bw.1 + 1])
		        	else grid.bw.1[pos.min.mise], 2*range.t)
      	}, # End of bootstrap
		bw <- fixed.bw	# fixed 
	) # End of switch(n.bw.selec)
   	# Estimate
	# With fixed presmoothing bandwidth = 0, non-presmoothed estimates are obtained
	p <- 
		if(n.bw.selec == 3 & bw[1] == 0) as.numeric(dfr$d)
        else .C("nadarayawatson", dfr$t, n, dfr$t, dfr$d, n, bw[1], n.kernel, phat = numeric(n), PACKAGE = "survPresmooth")$phat
    if(forh & (n.bound == 4) & (max(x.est) > bw[2]) & all((x.est < bw[2]) | (max(dfr$t) - x.est < bw[2])))
		warning("Possibly anomalous estimate due to a too large smoothing bandwidth: try correcting only one of the two boundary effects")
 	estim <- .C("presmestim", x.est, le.x.est, dfr$t, n, if(forh) bw[2], if(forh) n.kernel, if(forh) n.bound, p, n.estimand, pest = numeric(le.x.est), PACKAGE = "survPresmooth")$pest
	ps <- list(call = match.call(), data = if(control$save.data) dfr, q.weight = q.w, bw.selec = bw.selec, mise = if((n.bw.selec == 2) & control$save.mise) if(forh) matrix(mise*step.ise/n.boot, nrow = le.bw.1, ncol = le.bw.2) else mise*step.ise/n.boot, grid.pil = if(n.bw.selec != 3) grid.bw.pil, pilot.bw = if(n.bw.selec != 3) pilot.bw, bandwidth = bw, grid.bw = if(n.bw.selec == 2) list(grid.bw.1 = grid.bw.1, grid.bw.2 = grid.bw.2), p.hat = p, x.est = x.est, estimand = estimand, estimate = estim)
	class(ps) <- "survPresmooth"
	ps
}

print.survPresmooth <- 
function(x, more = NULL, ...){
    dots <- list(...)
    if(class(x) != "survPresmooth") stop("x object must be of class survPresmooth")
    cat("\nPresmoothed estimation of",
    switch(x$estimand, F = "distribution function, F(t)",
     f = "density function",
     H = "cumulative hazard function, H(t)",
     h = "hazard function, h(t)"),
    "\n\n")
    df <- data.frame(x$x.est, x$estimate)
    colnames(df) <- c("t", paste(x$estimand,"(t)", sep = ""))
    print(df, ...)
    cat("\n\nBandwidth selection method:", x$bw.selec, "\n\n")
    cat("\nBandwidth(s): \n")
    cat("\tpresmoothing:\t", ifelse(is.null(dots$digits), x$bandwidth[1], round(x$bandwidth[1], dots$digits)), "\n")
    if(length(x$bandwidth) == 2) cat("\tsmoothing:\t", ifelse(is.null(dots$digits), x$bandwidth[2], round(x$bandwidth[2], dots$digits)), "\n")
    if(!is.null(more)){
      ext.names <- c(call = "Call:",data = "Data:", q.weight = "Range of the weight function:",bw.selec = "bw.selec:", mise = "Bootstrap MISE:", grid.pil = "Grid used for the pilot bandwidth:", pilot.bw = "Pilot bandwidth(s):", bandwidth = "Bandwidth(s):", grid.bw = "Grid used for bandwidth selection:", p.hat = "Estimate of the conditional probability of uncensoring:", x.est = "x.est:", estimand = "Estimand:")
      inx <- more %in% names(x)
      if(any(inx)){
        cat("\nOther items: \n")
        invisible(sapply(1:length(more[inx]), function(s){
			cat("\t", ext.names[more][s],"\n")
			print(x[[more[inx][s]]], ...)
			}))
      if(!all(inx))
        if(sum(!inx) == 1 )
          warning(paste(more[!inx], "is not a component of x"))
        else
          warning(paste(paste(more[!inx], collapse = " , "), "are not components of x"))
      }
      else
        if(sum(!inx) == 1 )
          warning(paste(more[!inx], "is not a component of x"))
        else
          warning(paste(paste(more[!inx], collapse = " , "), "are not components of x"))
    }
}

