############################################################################
# Spatial Interpolation Lib by Pablo Alfaro (pabloalfaropineyro@gmail.com) #
############################################################################
# A complete library for performing spatial interpolation on space time    #
# series data                                                              #
# Allows incorporation of gridded auxiliary fields to assist on            #
# interpolation modelling via techniques from the Regression Kriging family#
# as well as provides functionalities for plotting and storing the results.#
# This library was developed as the graduation project for the author's    #
# Masters Degree on Mathematical Engineering                               #
############################################################################
# Copyright (C) 2012 Pablo Alfaro                                          #
# This program is free software: you can redistribute it and/or modify it  #
# under the terms of the GNU General Public License as published by the    #
# Free Software Foundation, either version 3 of the License, or (at your   #
# option) any later version.                                               #
# This program is distributed in the hope that it will be useful, but      #
# WITHOUT ANY WARRANTY; without even the implied warranty of               #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General #
# Public License for more details.                                         #
# You should have received a copy of the GNU General Public License along  #
# with this program. If not, see http://www.gnu.org/licenses/.             #
############################################################################

autokrigeVGMPanelMod <- function (x, y, model, subscripts, ...) {
  vgm.panel.xyplot(x, y, model = model, subscripts = TRUE, ...)
  no_digits = function(a) {
    if (a > 10) { return(0)
    } else {
      if (a < 1) { return(2) 
      } else { return(1) }
    }
  }
  
  if (model$model[1] == 'Nug') {
    nugget = sum(model[1, "psill"])
    sill = sum(model[, "psill"])
    range = sum(model[, "range"])
    txt = paste("Model: ", as.character(model[2, "model"]), 
                "\nNugget: ", round(nugget, digits = no_digits(nugget)), 
                "\nSill: ", round(sill, digits = no_digits(sill)), "\nRange: ", 
                round(range, digits = no_digits(range)), sep = "")
    if (model[2, "model"] %in% c("Mat", "Ste")) {
      kappa = model[2, "kappa"]
      txt = paste(txt, "\nKappa: ", round(kappa, digits = no_digits(kappa)), sep = "")
    }
  } else {
    nugget = 0
    sill = sum(model[, "psill"])
    range = sum(model[, "range"])
    txt = paste("Model: ", as.character(model[1, "model"]), 
                "\nNugget: ", round(nugget, digits = no_digits(nugget)), 
                "\nSill: ", round(sill, digits = no_digits(sill)), "\nRange: ", 
                round(range, digits = no_digits(range)), sep = "")
    if (model[1, "model"] %in% c("Mat", "Ste")) {
      kappa = model[1, "kappa"]
      txt = paste(txt, "\nKappa: ", round(kappa, digits = no_digits(kappa)), sep = "")
    }
  }
  ltext(x = max(x) * 0.02, y = max(y) * 0.98, txt, font = 2, adj = c(0, 1), cex = 0.7, col = grey(0.3))    
}

annotatedplot <- function(krigeobj, xlab = "Distance", ylab = "Semi-variance",
                          main = "Experimental variogram and fitted variogram model",
                          col.points='black', col.line='red', cex=1, pch=16, lwd=2,
                          fontsize=16){
  require("lattice")
  if ('variogramCloud' %in% class(krigeobj$exp_var)) {
    vdf <- as.data.frame(krigeobj$exp_var)
    pointLabels <- paste0(vdf$left, ',', vdf$right)
  } else {
    pointLabels <- as.character(krigeobj$exp_var$np)
  }
  
  ylim = c(min(0, 1.04 * min(krigeobj$exp_var$gamma)), 1.04 * max(krigeobj$exp_var$gamma))
  if (ylim[1] == ylim[2]) ylim[2] <- ylim[1] + 1

  if (!is.null(krigeobj$var_model)) {
    annotatedplot <- xyplot(gamma ~ dist, data = krigeobj$exp_var, panel = autokrigeVGMPanelMod,
                            labels = pointLabels, shift = 0.03, model = krigeobj$var_model,
                            direction = c(krigeobj$exp_var$dir.hor[1], krigeobj$exp_var$dir.ver[1]),
                            ylim = ylim,
                            xlim = c(0, 1.04 * max(krigeobj$exp_var$dist)), xlab = xlab, ylab = ylab,
                            main = main, mode = "direct", 
                            par.settings=simpleTheme(col.points=col.points, col.line=col.line, cex = cex, pch=pch, lwd=lwd))
  } else {
    annotatedplot <- xyplot(gamma ~ dist, data = krigeobj$exp_var, labels = pointLabels, shift = 0.03, 
                            direction = c(krigeobj$exp_var$dir.hor[1], krigeobj$exp_var$dir.ver[1]),
                            ylim = ylim, xlim = c(0, 1.04 * max(krigeobj$exp_var$dist)), xlab = xlab, 
                            ylab = ylab, main = main, mode = "direct", 
                            par.settings=simpleTheme(col.points=col.points, col.line=col.line, cex = cex, pch=pch, lwd=lwd))
  }
  
  annotatedplot <- update(annotatedplot, par.settings = list(fontsize = list(text = fontsize, points = 8)))
  
  return(annotatedplot)
}

estimateParameters.idw_mod <- function (object, ...) {
  params = getIntamapParams(object$params, ...)
  idpRange = params$idpRange
  if (is.null(idpRange)) 
    idpRange = seq(0.1, 2.9, 0.1)
  nfold = params$nfold
  if (is.null(nfold)) 
    nfold = 5
  mse = rep(NA_real_, length(idpRange))
  if ("formulaString" %in% names(object)) { formulaString = object$formulaString
  } else { formulaString = as.formula("value ~ 1") }
  dots = list(...)
  if ("nmax" %in% names(dots)) {
    nmax = dots$nmax
  }
  else nmax = object$params$nmax

  if (nfold==nrow(object$observations)) {
    nclus <- params$nclus
    if (is.null(nclus)) { nclus <- 1
    } else if (nclus > length(idpRange)) { nclus <- length(idpRange) }

    # More than one core to use and a minimum amount of work to justify creating the threads
    if ((nclus > 1) & (length(idpRange) * nfold / nclus >= 200)) {
      length(idpRange) * nrow(object$observations) / 8
      cl <- makeCluster(getOption("cl.cores", nclus))
      if (exists(x = 'setMKLthreads')) { clusterEvalQ(cl = cl, expr = setMKLthreads(1)) }
      mse <- parSapplyLB(cl=cl, X=seq(along = idpRange), FUN = function(i, formulaString, object, nfold, nmax, idpRange, verbose) {
        require('gstat')
        return(mean(krige.cv(formulaString, object$observations, 
                             nfold = nfold, nmax = nmax, set = list(idp = idpRange[i]), 
                             verbose = verbose)$residual^2))
      }, formulaString=formulaString, object=object, nfold=nfold, nmax=nmax, idpRange=idpRange, verbose=params$debug.level)
      stopCluster(cl)
      best = which.min(mse)
    } else {
      i <- 1
      minMSE <- .Machine$double.xmax
      oldMSE <- minMSE
      while (i <= length(idpRange)) {
        mse <- mean(krige.cv(formulaString, object$observations, 
                             nfold = nfold, nmax = nmax, set = list(idp = idpRange[i]), 
                             verbose = params$debug.level)$residual^2)
        if (mse < minMSE) {
          minMSE <- mse
          best <- i
        }
        if (mse < oldMSE) {
          i <- i + 1
          oldMSE <- mse
        } else {
          i <- length(idpRange) + 1
        }
      }      
    }    
  } else {
    for (i in seq(along = idpRange)) {
      mse[i] = mean(krige.cv(formulaString, object$observations, 
                             nfold = nfold, nmax = nmax, set = list(idp = idpRange[i]), 
                             verbose = params$debug.level)$residual^2)
    }
    best = which.min(mse)
  }
  
  #plot(x=idpRange, y=mse2)
  object$inverseDistancePower = idpRange[best]
  if (params$debug.level) print(paste("best idp value found is", object$inverseDistancePower, "rmse", sqrt(mse[best])))
  return(object)
}

getExperimentalVariogram <- function(formula, input_data, verbose = FALSE, GLS.model = NA, miscFitOptions = list(), ...) {
  # Check for anisotropy parameters
  if('alpha' %in% names(list(...))) warning('Anisotropic variogram model fitting not supported, see the documentation of autofitVariogram for more details.')
  
  # Take the misc fit options and overwrite the defaults by the user specified ones
  miscFitOptionsDefaults = list(merge.small.bins = TRUE, min.np.bin = 5,
                                orig.behavior = TRUE,
                                num.bins = NA, equal.width.bins = FALSE,
                                equal.np.bins = FALSE)
  miscFitOptions = modifyList(miscFitOptionsDefaults, miscFitOptions)
  
  # Create boundaries
  longlat = !is.projected(input_data)
  if(is.na(longlat)) longlat = FALSE
  diagonal = spDists(t(bbox(input_data)), longlat = longlat)[1,2] * 0.35 # 0.35 times the length of the central axis through the area
  
  ### BEGIN MODIFICATIONS ###
  
  if(miscFitOptions[["orig.behavior"]]){
    if(verbose) cat ("Boundaries as defined by original autofitVariogram...\n\n")
    # compute boundaries the old way
    boundaries = c(2,4,6,9,12,15,25,35,50,65,80,100) * diagonal/100 # Boundaries for the bins in km
    # If you specifiy a variogram model in GLS.model the Generelised least squares sample variogram is constructed
    if(!is(GLS.model, "variogramModel")) {
      experimental_variogram = variogram(formula, input_data, boundaries = boundaries, ...)
    } else {
      if(verbose) cat("Calculating GLS sample variogram\n")
      g = gstat(NULL, "bla", formula, input_data, model = GLS.model, set = list(gls=1))
      experimental_variogram = variogram(g, boundaries = boundaries, ...)
    }
    # request by Jon Skoien
    if(miscFitOptions[["merge.small.bins"]]) {
      # bin the old way
      if(verbose) cat("Checking if any bins have less than 5 points, merging bins when necessary...\n\n")
      while(TRUE) {
        if(length(experimental_variogram$np[experimental_variogram$np < miscFitOptions[["min.np.bin"]]]) == 0 | length(boundaries) == 1) break
        boundaries = boundaries[2:length(boundaries)]
        if(!is(GLS.model, "variogramModel")) {
          experimental_variogram = variogram(formula, input_data, boundaries = boundaries, ...)
        } else {
          experimental_variogram = variogram(g, boundaries = boundaries, ...)
        }
      }
    }
    ### equal-width bins (approximately, variogram does its own binning too) ###
  } else if(miscFitOptions[["equal.width.bins"]]){
    if(verbose) cat("Using equal-width bins...\n")
    if('width' %in% names(list(...))) stop('Cannot pass width when equal.width.bins = TRUE. Supply "init.width" in miscFitOptions instead.')
    # replace diagonal with cutoff, if provided
    if('cutoff' %in% names(list(...))) diagonal <- list(...)[['cutoff']]
    # user must supply either width or num.bins
    if(!'init.width' %in% names(miscFitOptions)){
      if(is.na(miscFitOptions[['num.bins']])) stop('when equal.width.bins = TRUE, user must supply either init.width or num.bins as well.')
      width <- diagonal/miscFitOptions[['num.bins']]
      if(verbose) cat("initial width not provided. Calculating using num.bins.\n")
    } else {
      width <- miscFitOptions[['init.width']]
      if(verbose) cat("initial width provided.\n")
    }
    # get the experimental variogram
    if(!is(GLS.model, "variogramModel")) {
      experimental_variogram = variogram(formula, input_data, width = width, ...)
    } else {
      experimental_variogram = variogram(g, width = width, ...)
    }
    # merge small bins if requested
    if(miscFitOptions[['merge.small.bins']]){
      if(verbose) cat("Checking if any bins have less than ", miscFitOptions[["min.np.bin"]], " points, merging bins when necessary...\n")
      iter <- 0
      maxiter <- 1000
      while(TRUE){
        if(!any(experimental_variogram$np < miscFitOptions[["min.np.bin"]])) break                        
        # increase width by 10% and try again
        width <- width*1.1
        if(!is(GLS.model, "variogramModel")) {
          experimental_variogram = variogram(formula, input_data, width = width, ...)
        } else {
          experimental_variogram = variogram(g, width = width, ...)
        }
        iter <- iter + 1
        if(iter > maxiter){
          cat('maximum number of interations reached. Try decreasing min.np.bin or init.width.\n\n')
          break
        }
      }
    }
    ### equal observation count bins ###
  } else if(miscFitOptions[["equal.np.bins"]]){
    if(verbose) cat("Using bins of equal observation counts...\n")
    if('boundaries' %in% names(list(...))) stop('Cannot pass boundaries when equal.np.bins is TRUE. Pass num.bins or min.np.bin instead.')
    # replace diagonal with cutoff, if provided
    if('cutoff' %in% names(list(...))) diagonal <- list(...)[['cutoff']]
    # get a sorted list of distances
    dists <- sort(spDists(input_data))
    # apply the cutoff
    dists <- dists[dists < diagonal & dists > 0]
    # split the data into bins based on number of observations
    if(is.na(miscFitOptions[['num.bins']])){
      # compute number of bins based on the minimum number of observations per bin
      miscFitOptions[['num.bins']] <- floor(0.5*length(dists)/miscFitOptions[['min.np.bin']])
      if(verbose) cat("num.bins not supplied. Setting num.bins =", miscFitOptions[['num.bins']], 'based on min.np.bin.\n')
    }
    cat("checking bins, decreasing num.bins if necessary... \n")
    while(TRUE){
      # compute interval based on the number of bins
      interval <- length(dists)/miscFitOptions[['num.bins']]
      # define boundaries
      boundaries <- rep(NA_real_, miscFitOptions[['num.bins']])
      for(i in 1:miscFitOptions[['num.bins']]){
        boundaries[i] <- dists[round(i*interval)]
      }
      if(length(boundaries == length(unique(boundaries)))) break
      # reduce number of bins
      miscFitOptions[['num.bins']] <- miscFitOptions[['num.bins']] - 1
    }
    if(!is(GLS.model, "variogramModel")) {
      experimental_variogram = variogram(formula, input_data, boundaries = boundaries, ...)
    } else {
      if(verbose) cat("Calculating GLS sample variogram\n")
      g = gstat(NULL, "bla", formula, input_data, model = GLS.model, set = list(gls=1))
      experimental_variogram = variogram(g, boundaries = boundaries, ...)
    }
  } else {
    # default behavior of variogram
    if(verbose) cat("No binning action specified in miscFitOptions.\n\n")
    if(!is(GLS.model, "variogramModel")) {
      experimental_variogram = variogram(formula, input_data, ...)
      # plot(experimental_variogram$dist, experimental_variogram$gamma)
      # experimental_variogram = variogram(formula, input_data, boundaries = boundaries, cressie = cressie)
    } else {
      if(verbose) cat("Calculating GLS sample variogram\n")
      g = gstat(NULL, "bla", formula, input_data, model = GLS.model, set = list(gls=1))
      experimental_variogram = variogram(g, ...)
    }
  }
  ### END MODIFICATIONS ###
  return (experimental_variogram)
}

isInvalidVariogram <- function(v, minPsill = 1E-3, minRange = 0.05) {
  return((is.null(v) |
   length(v) == 0 |
   any(!v$model %in% c('Nug', "Pow")  & ((v$psill < minPsill) | (v$range <= minRange))) |
   any(v$model == 'Nug' & v$psill < 0) |
   any(v$model == 'Pow' & (v$range > 2 | v$range <= 0))))
}

getModelVariogram <- function(experimental_variogram, formula, input_data = NULL, model = c("Sph", "Exp", "Gau", "Ste"),
                              kappa = c(0.05, seq(0.2, 2, 0.1), 5, 10), fix.values = c(NA,NA,NA), 
                              verbose = FALSE, start_vals = c(NA,NA,NA), fit.method = 7, nPuntosIniciales=3,
                              tryFixNugget=F, data=NULL, quitarNuggetSiEsCero=T, ...) {
  
  # set initial values
  # Nugget
  fit_nugget = is.na(fix.values[1])
  if(!fit_nugget) {
    initial_nugget = fix.values[1]
  } else {
    if (tryFixNugget) {
      if(is.na(start_vals[1])) { 
        initial_nugget = c(0, seq(from = 0, to = min(experimental_variogram$gamma), length.out = nPuntosIniciales))
        fit_nugget <- c(F, rep(T, nPuntosIniciales))
      } else { 
        initial_nugget = c(start_vals[1], seq(from = 0, to = min(experimental_variogram$gamma), length.out = nPuntosIniciales))
        fit_nugget <- c(F, rep(T, nPuntosIniciales))
      }
    } else {
      if(is.na(start_vals[1])) { initial_nugget = seq(from = 0, to = min(experimental_variogram$gamma), length.out = nPuntosIniciales)
      } else { initial_nugget = c(start_vals[1], seq(from = 0, to = min(experimental_variogram$gamma), length.out = nPuntosIniciales)) }
    }
  }

  # Range
  fit_range = is.na(fix.values[2])
  if(!fit_range) {
    initial_range = fix.values[2]
  } else {
    if(is.na(start_vals[2])) {
      initial_range <- c(experimental_variogram$dist[nrow(experimental_variogram)]*0.2, seq(from=experimental_variogram$dist[2], to=experimental_variogram$dist[nrow(experimental_variogram)], length.out = nPuntosIniciales))
    } else {
      initial_range <- c(start_vals[2], experimental_variogram$dist[nrow(experimental_variogram)]*0.1, seq(from=experimental_variogram$dist[2], to=experimental_variogram$dist[nrow(experimental_variogram)], length.out = nPuntosIniciales))
    }
  }
  
  # Partial sill
  fit_sill = is.na(fix.values[3])
  if(!fit_sill) {
    initial_sill = fix.values[3]
  } else {
    if(is.na(start_vals[3])) { # Sill
      initial_sill = c(mean(c(max(experimental_variogram$gamma), median(experimental_variogram$gamma))), seq(from=min(experimental_variogram$gamma), to = max(experimental_variogram$gamma), length.out = nPuntosIniciales))
    } else {
      initial_sill = c(start_vals[3], mean(c(max(experimental_variogram$gamma), median(experimental_variogram$gamma))), seq(from=min(experimental_variogram$gamma), to = max(experimental_variogram$gamma), length.out = nPuntosIniciales))
    }
  }
  
  initVals <- cbind(fit_nugget, expand.grid(initial_nugget=initial_nugget, initial_range=initial_range, initial_sill=initial_sill))
  if ('Pow' %in% model) { initValsPow <- cbind(fit_nugget, expand.grid(initial_nugget=initial_nugget, initial_range=c(0.2, 0.5, 1, 1.5, 2), initial_sill=initial_sill)) }
  
  getModel = function(psill, model, range, kappa, nugget, fit_range, fit_sill, fit_nugget, fit.method, verbose) {
    if(verbose) debug.level = 1 else debug.level = 0
    if(fit.method==5) {
      #fit.variogram.reml(formula, locations, data, model, debug.level = 1, set, degree = 0)
      #        obj = try(fit.variogram.reml(experimental_variogram,
      #                                         model=vgm(psill=psill, model=model, range=range,
      #                                                        nugget=nugget,kappa = kappa),
      #                         fit.ranges = c(fit_range), fit.sills = c(fit_nugget, fit_sill),
      #                         debug.level = 0, fit.method=fit.method), TRUE)                
    } else if(fit.method==8) {
      obj = try(fit.variogram.gls(formula, data, model=vgm(psill = psill, model = model, range = range, nugget = nugget, kappa = kappa), ignoreInitial = F, trace = F))
    } else {
      #if (!fit_nugget & nugget == 0) { model=vgm(psill=psill, model=model, range=range, kappa = kappa)
      #} else { model=vgm(psill=psill, model=model, range=range, nugget=nugget, kappa = kappa) }
      obj = try(fit.variogram(experimental_variogram,
                              model=vgm(psill=psill, model=model, range=range, nugget=nugget, kappa = kappa),
                              fit.ranges = c(fit_range), fit.sills = c(fit_nugget, fit_sill),
                              debug.level = 0, fit.method=fit.method), TRUE)
    }
    if("try-error" %in% class(obj)) {
      #print(traceback())
      warning("An error has occured during variogram fitting. Used:\n",
              "\tnugget:\t", nugget,
              "\n\tmodel:\t", model,
              "\n\tpsill:\t", psill,
              "\n\trange:\t", range,
              "\n\tkappa:\t",ifelse(kappa == 0, NA, kappa),
              "\n as initial guess. This particular variogram fit is not taken into account. \nGstat error:\n", obj)
      return(NULL)
    } else return(obj)
  }
  
  # Automatically testing different models, the one with the smallest sums-of-squares is chosen
  test_models = model
  SSerr_list = c()
  vgm_list = list()
  counter = 1
  
  maxGamma <- max(experimental_variogram$gamma)
  minPsill <- maxGamma * 0.01
  minRange <- 0.05
  
  # verbose <- T
  # plot(experimental_variogram$dist, experimental_variogram$gamma, ylim=c(0, max(experimental_variogram$gamma)))
  # vCloud <- variogram(formula, observaciones, cloud=T, cutoff=max(boundaries))
  # plot(vCloud)
  m <- test_models[1]
  for(m in test_models) {
    bestFit <- NULL
    bestSSerr <- Inf
    
    if (m == 'Pow') { iv <- initValsPow
    } else { iv <- initVals }
    
    #n <- 1
    for (n in 1:nrow(iv)) {
      fit_nugget <- iv[n, 1]
      # test multiple starting values, keep best for each model
      nugget <- iv[n, 2]
      range <- iv[n, 3]
      psill <- iv[n, 4]
      
      if(m != "Mat" && m != "Ste") { # If not Matern and not Stein
        model_fit = getModel(
          psill = psill - nugget, model = m, range = range, kappa = 0, nugget = nugget, 
          fit_range = fit_range, fit_sill = fit_sill, fit_nugget = fit_nugget, verbose = verbose, 
          fit.method=fit.method)
        
        if (!isInvalidVariogram(model_fit, minPsill = minPsill, minRange = minRange) && attr(model_fit, "SSErr") < bestSSerr) { # skip models that failed
          bestFit <- model_fit
          bestSSerr <- attr(model_fit, "SSErr")
        }
      } else { # Else loop also over kappa values
        for(k in kappa) {
          model_fit = getModel(psill = psill - nugget, model = m, range = range, k, nugget = nugget, fit_range, fit_sill, fit_nugget, verbose = verbose, fit.method=fit.method)
          
          if (!isInvalidVariogram(model_fit, minPsill=minPsill, minRange = minRange) && attr(model_fit, "SSErr") < bestSSerr) { # skip models that failed
            bestFit <- model_fit
            bestSSerr <- attr(model_fit, "SSErr")
          }          
        }
      }
    }
    
    
    if (!is.null(bestFit)) {
      vgm_list[[counter]] = bestFit
      SSerr_list = c(SSerr_list, attr(bestFit, "SSErr"))
      counter = counter + 1
    }
  }
  
  # Check for negative values in sill or range coming from fit.variogram
  # and NULL values in vgm_list, and remove those with a warning
  strange_entries = sapply(X = vgm_list, FUN = isInvalidVariogram, minPsill=minPsill, minRange=minRange)
  if(any(strange_entries)) {
    if(verbose) {
      print(vgm_list[strange_entries])
      cat("^^^ ABOVE MODELS WERE REMOVED ^^^\n\n")
    }
    # warning("Some models where removed for being either NULL or having a negative sill/range/nugget, \n\tset verbose == TRUE for more information")
    SSerr_list = SSerr_list[!strange_entries]
    vgm_list = vgm_list[!strange_entries]
  }
  
  if(verbose) {
    if (length(SSerr_list) > 0) {
      cat("Selected:\n")
      print(vgm_list[[which.min(SSerr_list)]])
      cat("\nTested models, best first:\n")
      tested = data.frame("Tested models" = sapply(vgm_list, function(x) as.character(x[2,1])),
                          kappa = sapply(vgm_list, function(x) as.character(x[2,4])),
                          "SSerror" = SSerr_list)
      tested = tested[order(tested$SSerror),]
      print(tested)      
    } else {
      cat("No variogram fit produced a valid model!")
    }
  }
  
  if (length(vgm_list) > 0) { 
    #ordenSSErr <- order(SSerr_list)
    #SSerr_list <- SSerr_list[ordenSSErr]
    #vgm_list <- vgm_list[ordenSSErr]
    #maxErrorTolerable <- SSerr_list[1] + (SSerr_list[length(SSerr_list)] - SSerr_list[1]) * 0.25
    
    #idx <- SSerr_list < maxErrorTolerable
    #vgm_list <- vgm_list[idx]
    #SSerr_list <- SSerr_list[idx]
    
    #rangos <- sapply(vgm_list, getExperimentalRange, maxDist=experimental_variogram$dist[nrow(experimental_variogram)])
    #i <- which.max(rangos)    
    #var_model <- vgm_list[[i]]
    #sserr <- SSerr_list[[i]]
    
    var_model <- vgm_list[[which.min(SSerr_list)]]
    sserr <- min(SSerr_list)
  } else { 
    var_model = NULL 
    sserr <- Inf
  }
  
  if (quitarNuggetSiEsCero) {
    iNugget <- which(var_model == 'Nug')
    if (!is.null(var_model) && (!fit_nugget) && length(iNugget) > 0 && (var_model$psill[iNugget] == 0)) {
      var_model <- var_model[-iNugget,]
      row.names(var_model) <- as.character(1:nrow(var_model))
    }
  }
  
  result = list(exp_var = experimental_variogram, var_model=var_model, sserr = sserr)
  class(result) = c("autofitVariogram","list")
  
  return(result)
}

afvmod = function(formula, input_data, model = c("Sph", "Exp", "Gau", "Ste"),
                  kappa = c(0.05, seq(0.2, 2, 0.1), 5, 10), fix.values = c(NA,NA,NA),
                  verbose = FALSE, GLS.model = NA, start_vals = c(NA,NA,NA),
                  miscFitOptions = list(), fit.method = 7, nPuntosIniciales=3,
                  tryFixNugget=F, quitarNuggetSiEsCero=T, ...) {
  #experimental_variogram <- getExperimentalVariogram(formula = formula, input_data = input_data, verbose = verbose, 
  #                                                   GLS.model = GLS.model, miscFitOptions = miscFitOptions, boundaries=boundaries)
  experimental_variogram <- getExperimentalVariogram(formula = formula, input_data = input_data, verbose = verbose, 
                                                     GLS.model = GLS.model, miscFitOptions = miscFitOptions, ...)
  return (getModelVariogram(experimental_variogram=experimental_variogram, formula = formula, input_data = input_data, 
                            model = model, kappa = kappa, fix.values = fix.values, verbose = verbose, start_vals = start_vals, 
                            fit.method = fit.method, nPuntosIniciales=nPuntosIniciales, tryFixNugget = tryFixNugget, 
                            quitarNuggetSiEsCero=quitarNuggetSiEsCero, ...))
}

getBoundariesPVariogramaEmpiricoV4_MultiTime <- function(formula, input_data, input_data_t, cutoff=NA) {
  nFechas <- nrow(input_data_t[[1]])
  
  spatialData <- do.call("rbind", rep(list(sp::coordinates(input_data)), nFechas))
  spatialData <- data.frame(x=spatialData[,1], y=spatialData[,2])
  spatialData <- do.call("cbind", args = c(list(spatialData), lapply(input_data_t, function(x) { as.vector(t(x))})))
  spatialData <- na.omit(spatialData)
  sp::coordinates(spatialData) <- c('x', 'y')
  proj4string(spatialData) <- proj4string(input_data)
  # anis <- estimateAnisotropy(spatialData, depVar = 'value')
  boundaries <- getBoundariesPVariogramaEmpiricoV4(fml = formula, observaciones = spatialData, proporcionDeLaDiagonalValida = 1, minDivDist = 6, cutoff=cutoff)
}

getExperimentalSpatialVariogramMultiTime <- function(formula, input_data, input_data_t, cutoff=NA) {
  boundaries <- getBoundariesPVariogramaEmpiricoV4_MultiTime(formula, input_data, input_data_t, cutoff)
  
  expVars <- list()
  for (t in 1:nFechas) {
    aux <- input_data
    aux@data <- data.frame(do.call("cbind", args = lapply(input_data_t, function(x, t) { return (x[t,] ) }, t)))
    aux <- aux[complete.cases(aux@data),]
    if (nrow(aux) > ncol(aux) -1) { expVars[[length(expVars) + 1]] <- variogram(object = formula, data = aux, cloud=T, cutoff=boundaries[length(boundaries)]) }
  }

  res <- data.frame(np=0, dist=boundaries, gamma=0, dir.hor=0, dir.ver=0, id=as.factor('var1'))
  t <- 1
  for (t in 1:length(expVars)) {
    ies <- findInterval(expVars[[t]]$dist, res$dist, rightmost.closed = T, all.inside = T)
    res$gamma[ies] <- res$gamma[ies] + expVars[[t]]$gamma
    res$np[ies] <- res$np[ies] + 1
  }
  res$gamma <- res$gamma / nFechas
  res <- res[res$np != 0,]
  #res <- rbind(data.frame(np=nFechas, dist=0, gamma=0, dir.hor=0, dir.ver=0, id=as.factor('var1')), res)
  class(res) <- c("gstatVariogram", class(res))
  
  return(res)
}

getExperimentalSill <- function(variogramModel, maxDist, n=100) {
  if ('Pow' %in% variogramModel$model) { return(quantile(variogramLine(variogramModel, maxdist = maxDist, n=n)$gamma, 0.95))
  } else { return (sum(variogramModel$psill)) }
}

getExperimentalRange <- function(variogramModel, maxDist, n=100) {
  sill <- getExperimentalSill(variogramModel, maxDist, n)
  vgmLine <- variogramLine(variogramModel, maxdist = maxDist, n=n)
  diffs <- abs(vgmLine$gamma - sill * 0.95)
  # me fijo cuales son menores o iguales que el mínimo más 1% para dejar un márgen para variogramas casi chatos
  iMins <- which(diffs <= min(diffs) * 1.01)
  iDist <- iMins[length(iMins)]
  # iDist <- which.min(abs(vgmLine$gamma - sill * 0.95))
  return(vgmLine$dist[iDist])
}

estimateSTAnisotropyEx <- function(empVgm, interval=c(1, max(empVgm$spacelag)), methods = c('metric', 'vgm', 'range', 'linear'), spatialVgm, temporalVgm) {
  anis <- numeric()
  i <- 1
  for (i in 1:length(methods)) {
    obj <- tryCatch(expr = { estiStAni(empVgm = empVgm, interval = interval, method = methods[i], spatialVgm = spatialVgm, temporalVgm = temporalVgm) },
                    error = function(e) { 
                      class(e) <- c(class(e), 'try-error')
                      return (e) },
                    warning= function(w) { 
                      class(w) <- c(class(w), 'try-error')
                      return (w) } )
    
    if (!"try-error" %in% class(obj) && obj >= interval[1] && obj <= interval[2]) {
      anis <- c(anis, obj)
      names(anis)[length(anis)] <- methods[i]
    }
  }
  return (anis)
}

fit.StVariogramEx <- function (object, model, ..., method = "L-BFGS-B", fit.method = 6, stAni = NA, wles, errorLevel=2) {
  if (!inherits(object, "StVariogram")) stop("\"object\" must be of class \"StVariogram\"")
  if (!inherits(model, "StVariogramModel")) stop("\"model\" must be of class \"StVariogramModel\".")
  sunit <- attr(object$spacelag, "units")
  tunit <- attr(object$timelag, "units")
  tu.obj = attr(model, "temporal units")
  if (!is.null(tu.obj)) stopifnot(identical(tunit, tu.obj))
  object <- na.omit(object)
  ret <- model
  
  if (!missing(wles)) {
    if (wles) { fit.method = 1
    } else { fit.method = 6 }
  }
  if (fit.method == 0) {
    attr(ret, "optim.output") <- "no fit"
    attr(ret, "MSE") <- mean((object$gamma - variogramSurface(model, data.frame(spacelag = object$dist, timelag = object$timelag))$model)^2)
    attr(ret, "spatial unit") <- sunit
    attr(ret, "temporal unit") <- tunit
    return(ret)
  }
  if ((fit.method == 7 | fit.method == 11) & is.null(model$stAni) & is.na(stAni)) {
    warning("An uninformed spatio-temporal anisotropy value of '1 (spatial unit)/(temporal unit)' is automatically selected. Consider providing a sensible estimate for stAni or using a different fit.method.")
    stAni <- 1
  }  

  weightingFun <- switch(fit.method, function(obj, ...) obj$np, 
                         function(obj, gamma, ...) obj$np/gamma^2, 
                         function(obj, ...) obj$np, 
                         function(obj, gamma, ...) obj$np/gamma^2, 
                         function(obj, ...) stop("fit.method = 5 (REML), is not yet implemented"), 
                         function(obj, ...) 1, 
                         function(obj, curStAni, ...) { 
                           if (is.na(stAni)) { obj$np/(obj$dist^2 + (curStAni * obj$timelag)^2) 
                           } else { obj$np/(obj$dist^2 + (stAni * obj$timelag)^2) }
                         }, 
                         function(obj, ...) {
                           dist <- obj$dist
                           dist[dist == 0] <- min(dist[dist != 0], na.rm = T)
                           obj$np/dist^2
                         }, 
                         function(obj, ...) {
                           dist <- obj$timelag
                           dist[dist == 0] <- min(dist[dist != 0], na.rm = T)
                           obj$np/dist^2
                         }, 
                         function(obj, gamma, ...) 1/gamma^2, 
                         function(obj, curStAni, ...) {
                           if (is.na(stAni)) { 1/(obj$dist^2 + (curStAni * obj$timelag)^2) 
                           } else { 1/(obj$dist^2 + (stAni * obj$timelag)^2) } 
                         }, 
                         function(obj, ...) {
                           dist <- obj$dist
                           dist[dist == 0] <- min(dist[dist != 0], na.rm = T)
                           1/(obj$dist^2)
                         }, function(obj, ...) {
                           dist <- obj$timelag
                           dist[dist == 0] <- min(dist[dist != 0], na.rm = T)
                           1/(obj$timelag^2)
                         })

  if (is.null(weightingFun)) stop(paste("fit.method =", fit.method, "is not implementend"))
  
  fitFun = function(par, trace = FALSE, ...) {
    curModel <- gstat:::insertPar(par, model)
    gammaMod <- variogramSurface(curModel, data.frame(spacelag = object$dist, timelag = object$timelag))$model
    resSq <- (object$gamma - gammaMod)^2
    resSq <- resSq * weightingFun(object, gamma = gammaMod, curStAni = curModel$stAni)
    if (trace) print(c(par, MSE = mean(resSq)))
    mean(resSq)
  }
  
  result = tryCatch({
    eval(pars.fit <- optim(par = extractPar(model), fn=fitFun, ..., method = method))
  }, error = function(e) {
    if (errorLevel == 2) { stop(e)
    } else if (errorLevel == 1) { warning(e) }
    class(e) <- c(class(e), 'try-error')
    return(e)
  })

  if (!'try-error' %in% class(result)) {
    ret <- gstat:::insertPar(pars.fit$par, model)
    attr(ret, "optim.output") <- pars.fit
  } else {
    ret <- model
    attr(ret, "optim.output") <- "optim error"
    attr(ret, "optim") <- list(par=extractPar(model), value=fitFun(par = extractPar(model), ...))
  }

  attr(ret, "MSE") <- mean((object$gamma - variogramSurface(ret, data.frame(spacelag = object$dist, timelag = object$timelag))$model)^2)
  attr(ret, "spatial unit") <- sunit
  attr(ret, "temporal unit") <- tunit
  return(ret)
}


fit.StVariogramEx2 <- function (object, model, ..., method = "L-BFGS-B", fit.method = 6, stAni = NA, wles, errorLevel=2) {
  if (!inherits(object, "StVariogram")) stop("\"object\" must be of class \"StVariogram\"")
  if (!inherits(model, "StVariogramModel")) stop("\"model\" must be of class \"StVariogramModel\".")
  sunit <- attr(object$spacelag, "units")
  tunit <- attr(object$timelag, "units")
  tu.obj = attr(model, "temporal units")
  if (!is.null(tu.obj)) stopifnot(identical(tunit, tu.obj))
  object <- na.omit(object)
  ret <- model
  if (!missing(wles)) {
    if (wles) { fit.method = 1
    } else { fit.method = 6 }
  }
  if (fit.method == 0) {
    attr(ret, "optim.output") <- "no fit"
    attr(ret, "MSE") <- mean((object$gamma - variogramSurface(model, data.frame(spacelag = object$dist, timelag = object$timelag))$model)^2)
    attr(ret, "spatial unit") <- sunit
    attr(ret, "temporal unit") <- tunit
    return(ret)
  }
  if ((fit.method == 7 | fit.method == 11) & is.null(model$stAni) & is.na(stAni)) {
    warning("An uninformed spatio-temporal anisotropy value of '1 (spatial unit)/(temporal unit)' is automatically selected. Consider providing a sensible estimate for stAni or using a different fit.method.")
    stAni <- 1
  }
  
  weightingFun <- switch(fit.method, 
                         function(obj, ...) obj$np, 
                         function(obj, gamma, ...) obj$np/gamma^2, 
                         function(obj, ...) obj$np, 
                         function(obj, gamma, ...) obj$np/gamma^2, 
                         function(obj, ...) stop("fit.method = 5 (REML), is not yet implemented"), 
                         function(obj, ...) 1, 
                         function(obj, curStAni, ...) { 
                           if (is.na(stAni)) { obj$np/(obj$dist^2 + (curStAni * obj$timelag)^2) 
                             } else { obj$np/(obj$dist^2 + (stAni * obj$timelag)^2) }
                         }, 
                         function(obj, ...) {
                           dist <- obj$dist
                           dist[dist == 0] <- min(dist[dist != 0], na.rm = T)
                           obj$np/dist^2
                         }, 
                         function(obj, ...) {
                          dist <- obj$timelag
                          dist[dist == 0] <- min(dist[dist != 0], na.rm = T)
                          obj$np/dist^2
                         }, 
                         function(obj, gamma, ...) 1/gamma^2, 
                         function(obj, curStAni, ...) {
                           if (is.na(stAni)) { 1/(obj$dist^2 + (curStAni * obj$timelag)^2) 
                           } else { 1/(obj$dist^2 + (stAni * obj$timelag)^2) } 
                         }, 
                         function(obj, ...) {
                           dist <- obj$dist
                           dist[dist == 0] <- min(dist[dist != 0], na.rm = T)
                           1/(obj$dist^2)
                         }, function(obj, ...) {
                           dist <- obj$timelag
                           dist[dist == 0] <- min(dist[dist != 0], na.rm = T)
                           1/(obj$timelag^2)
                         })

  if (is.null(weightingFun)) stop(paste("fit.method =", fit.method, "is not implementend"))
  
  fitFun = function(par, trace = FALSE, ...) {
    curModel <- gstat:::insertPar(par, model)
    gammaMod <- variogramSurface(curModel, data.frame(spacelag = object$dist, timelag = object$timelag))$model
    resSq <- (object$gamma - gammaMod)^2
    resSq <- resSq * weightingFun(object, gamma = gammaMod, curStAni = curModel$stAni)
    if (trace) print(c(par, MSE = mean(resSq)))
    mean(resSq)
  }
  
  simpleVGMGr_nugget = function(simpleVgm, dist) {
    if (nrow(simpleVgm) > 1) { 
      res <- rep(0, length(dist))
      for (i in 1:nrow(simpleVgm)) res <- res + simpleVGMGr_nugget(simpleVgm = simpleVgm[i, ], dist = dist)
      return (res)
    } else if (simpleVgm$model == 'Nug') { return (rep(1, length(dist)))
    } else { return (rep(0, length(dist))) }
  }
  
  simpleVGMGr_psill = function(simpleVgm, dist) {
    if (nrow(simpleVgm) > 1) { 
      res <- rep(0, length(dist))
      for (i in 1:nrow(simpleVgm)) res <- res + simpleVGMGr_psill(simpleVgm = simpleVgm[i, ], dist = dist)
      return (res)
    } else if (simpleVgm$model == 'Nug') { return (rep(0, length(dist)))
    } else if (simpleVgm$model == 'Exp') { return (1 - exp(-dist / simpleVgm$range))
    } else if (simpleVgm$model == 'Sph') {
      gr <- (3 * dist) / (2 * simpleVgm$range) - dist^3 / (2 * simpleVgm$range^3)
      gr[dist >= simpleVgm$range] <- 1
      return(gr)
    } else if (simpleVgm$model == 'Pow') {
      minP <- 0
      maxP <- 2
      pCalc <- (maxP - minP) * (0.5 + atan(simpleVgm$psill) / pi) + minP;
      return(dist^pCalc)
    } else { 
      f <- function(x) {
        model <- simpleVgm
        model$psill <- x
        return (variogramLine(object=model, dist_vector=dist)$gamma)
      }
      return(jacobian(func=f, x=simpleVgm$psill))
    }
  }
  
  simpleVgm <- model$joint[2,]
  simpleVGMGr_range = function(simpleVgm, dist) {
    if (nrow(simpleVgm) > 1) { 
      res <- rep(0, length(dist))
      for (i in 1:nrow(simpleVgm)) res <- res + simpleVGMGr_range(simpleVgm = simpleVgm[i, ], dist = dist)
      return (res)
    } else if (simpleVgm$model == 'Nug') { return (rep(0, length(dist)))
    } else if (simpleVgm$model == 'Exp') { return (-simpleVgm$psill / simpleVgm$range^2 * dist * exp(-dist / simpleVgm$range))
    } else if (simpleVgm$model == 'Sph') {
      gr <- 3 * simpleVgm$psill * (dist^3 - dist * simpleVgm$range^2) / (2 * (simpleVgm$range^4))
      gr[dist >= simpleVgm$range] <- 0
      return (gr)
    } else if (simpleVgm$model == 'Pow') {
      minP <- 0
      maxP <- 2
      pCalc <- (maxP - minP) * (0.5 + atan(simpleVgm$psill) / pi) + minP;
      gr <- simpleVgm$range * (dist^pCalc) * log(dist + 1E-8) * ((maxP - minP) / (pi * (1 + simpleVgm$psill^2)));
    } else { 
      f <- function(x) {
        model <- simpleVgm
        model$range <- x
        return (variogramLine(object=model, dist_vector=dist)$gamma)
      }
      return(jacobian(func=f, x=simpleVgm$range))
    }
  }
  
  dist <- object$dist
  timeDist <- object$timelag

  grFun = function(par, traceGr = FALSE, ...) {
    curModel <- gstat:::insertPar(par, model)
    
    if (is.na(stAni)) { curStAni <- curModel$stAni
    } else { curStAni <- stAni }
    jointDist <- sqrt(dist^2 + (timeDist*curStAni)^2)
    
    if (model$stModel == 'metric') {
      # vst(h, u) = vj(sqrt(h^2+(K*u)^2))
      # vst(h, u) d/dsill = vj(sqrt(h^2+(K*u)^2)) d/dsill
      # vst(h, u) d/drange = vj(sqrt(h^2+(K*u)^2)) d/drange
      # vst(h, u) d/dnugget = vj(sqrt(h^2+(K*u)^2)) d/dnugget
      # vst(h, u) d/danis = 0
      grVST <- cbind(simpleVGMGr_psill(simpleVgm = model$joint, dist = jointDist), 
                     simpleVGMGr_range(simpleVgm = model$joint, dist = jointDist), 
                     simpleVGMGr_nugget(simpleVgm = model$joint, dist = jointDist), 
                     rep(0, length(jointDist)))
    #} else if (model$stModel == 'simpleSumMetric') {
      # vst(h, u) = nug * 1h>0,v>0 + vs(h) + vt(u) + vj(sqrt(h^2+(K*u)^2))
      # vst(h, u) dsill.s = vs(h) dsill
      # vst(h, u) drange.s = vs(h) drange
      # vst(h, u) dsill.t = vt(u) dsill
      # vst(h, u) drange.t = vt(u) drange
      # vst(h, u) dsill.st = vj(sqrt(h^2+(K*u)^2)) dsill
      # vst(h, u) drange.st = vj(sqrt(h^2+(K*u)^2)) drange
      # vst(h, u) dnug = 1h>0,v>0
      # vst(h, u) danis = vj(sqrt(h^2+(K*u)^2))
    } else {
      f <- function(x) {
        curModel <- gstat:::insertPar(x, model)
        return (variogramSurface(curModel, data.frame(spacelag = object$dist, timelag = object$timelag))$model)
      }
      grVST <- jacobian(func=f, x=par)
    }
    
    f <- function(x) {
      curModel <- gstat:::insertPar(x, model)
      return (variogramSurface(curModel, data.frame(spacelag = object$dist, timelag = object$timelag))$model)
    }
    
    #if (max(abs(grVST - jacobian(func=f, x=par))) > 1E-1) print(grVST - jacobian(func=f, x=par))

    # f(par) = 1/(Nu*Nv) * sum_u (sum_v w(par, u, v) * (ev(u,v) - vst(par, u, v))^2)
    # f'(par) = 1/(Nu*Nv) * sum_u (sum_v (w(par, u,v) * 2 * (vst(par, u, v) - ev(u,v)) * vst'(par, u, v)) + (ev(u,v) - vst(par, u, v))^2 * w'(par, u, v))
    # where: f is the objective function to minimize
    # w is the weigthing function
    # ev is te empirical variogram (object$gamma)
    # vst is the space time variogram model surface
    vst <- variogramSurface(curModel, data.frame(spacelag = object$dist, timelag = object$timelag))$model
    w <- weightingFun(object, gamma = vst, curStAni = curStAni)

    # If we go into the if below, the Weighting function is not constant with respect to par (anisotropy) and it's gradient is not zero
    # If we go into the else, the Weighting function is constant with respect to par and it's gradient is zero so the above expression is simplified
    if (is.na(stAni) && fit.method %in% c(2, 4, 7, 10, 11)) {
      grW <- matrix(data = 0, nrow = nrow(object), ncol = length(par))    
      iAnis <- which(names(par) == 'anis')
      if (fit.method %in% c(2, 4)) { grW[, iAnis] <- object$np * grVgmModel[, iAnis] / vst^2
      } else if (fit.method == 7) { grW[, iAnis] <- object$np * - 2 * timeDist^2 * curStAni / (dist^2 + (timeDist*curStAni)^2)^2
      } else if (fit.method == 10) { grW[, iAnis] <- grVgmModel[, iAnis] / vst^2 
      } else if (fit.method == 11) { grW[, iAnis] <- - 2 * timeDist^2 * curStAni / (dist^2 + (timeDist*curStAni)^2)^2 }
      
      gr <- (((w * 2 * (vst - object$gamma)) %*% grVST) + ((object$gamma - vst) ^ 2) %*% grW) / nrow(object)
    } else { gr <- ((w * 2 * (vst - object$gamma)) %*% grVST) / nrow(object) }
    
    if (traceGr) print(gr)
    return (gr)
  }

  result = tryCatch({
    eval(pars.fit <- optim(par = extractPar(model), fn=fitFun, gr=grFun, ..., method = method))
  }, error = function(e) {
    if (errorLevel == 2) { stop(e)
    } else if (errorLevel == 1) { warning(e) }
    class(e) <- c(class(e), 'try-error')
    return(e)
  })
  
  if (!'try-error' %in% class(result)) {
    ret <- gstat:::insertPar(pars.fit$par, model)
    attr(ret, "optim.output") <- pars.fit
  } else {
    ret <- model
    attr(ret, "optim.output") <- "optim error"
    attr(ret, "optim") <- list(par=extractPar(model), value=fitFun(par = extractPar(model), ...))
  }
  
  attr(ret, "MSE") <- mean((object$gamma - variogramSurface(ret, data.frame(spacelag = object$dist, timelag = object$timelag))$model)^2)
  attr(ret, "spatial unit") <- sunit
  attr(ret, "temporal unit") <- tunit
  return(ret)
}

as.STSDF.STFDFEx = function(from) {
  # take out the NA cells and fill the index
  # NOTE: does not (yet) take out empty space/time entities 
  # -- should this be optional?
  n = length(from@sp)
  m = nrow(from@time)
  index = cbind(rep(1:n, m), rep(1:m, each=n))
  # copied from sp:
  if (nrow(from@data) > 1) { sel = apply(sapply(from@data, is.na), 1, function(x) !any(x))
  } else { sel = !any(sapply(from@data, is.na)) }

  index = index[sel,,drop=FALSE]
  STSDF(from@sp, from@time, from@data[sel,,drop=FALSE], index, from@endTime)
}

# STSDF -> STIDF
as.STIDF.STSDFEx = function(from) {
  # replicate the sp and time columns; keeps time always ordered?
  sp = from@sp[from@index[,1],]
  if (is(sp, "SpatialPoints"))
    row.names(sp) = make.unique(as.character(row.names(sp)))
  STIDF(sp, from@time[from@index[,2]], 
        from@data,
        from@endTime[from@index[,2]])
}

fit.variogram.gls_mod <- function(formula, data, model, maxiter = 30, eps = 0.01, trace = TRUE, ignoreInitial = TRUE, cutoff = Inf, plot = FALSE) {
  v = as.data.frame(variogram(formula, data, cloud = TRUE, cutoff = cutoff))
  i = v$left
  j = v$right
  y = v$gamma
  h0 = v$dist
  dists = spDists(data)
  n = length(i)
  iter = 0
  converged = FALSE
  # length(v$dist)^2 * 8 / 1024 / 1024
  
  includeNugget <- model$model[1] == 'Nug'
  if (includeNugget) {
    if (ignoreInitial) { init = c(mean(y)/2, mean(y)/2, median(h0)/4)
    } else { init = c(model$psill, model$range[2]) }
    
    gamfn = function(h, th, m = as.character(model$model[2])) variogramLine(vgm(th[2], m, th[3], th[1]), dist_vector = h)$gamma
    minfuncols = function(theta = rbind(1, 1, 1)) {
      res = y - gamfn(h0, theta)
      sum(res^2)
    }
    
    if (any(model$model == 'Pow')) { 
      upperOptim <- c(max(y), max(y), 2)
      init[3] <- 1
    } else { upperOptim = c(max(y), max(y), max(h0)) }
    
    th = th0 = th.ols = optim(init, minfuncols, gr = NULL, method = "L-BFGS-B", 
                              lower = c(0, 0.000000001, 0.000000001), 
                              upper = upperOptim)$par
    if (trace) 
      print(th)
    while (!converged && iter < maxiter) {
      comb = function(i, j) cbind(rep(i, length(j)), rep(j, each = length(i)))
      cov = matrix(variogramLine(model, dist_vector = dists[comb(i, i)])$gamma + 
                   variogramLine(model, dist_vector = dists[comb(j, j)])$gamma - 
                   variogramLine(model, dist_vector = dists[comb(i, j)])$gamma - 
                   variogramLine(model, dist_vector = dists[comb(j, i)])$gamma, length(j), length(j))
      cov = 0.5 * cov^2
      cov = qr(cov)
      minfuncrange = function(range) {
        res = y - gamfn(h0, c(th[1], th[2], range))
        t(res) %*% solve(cov, res)
      }
      minfuncsill = function(sill) {
        res = y - gamfn(h0, c(th[1], sill, th[3]))
        t(res) %*% solve(cov, res)
      }
      minfuncnugget = function(nugget) {
        res = y - gamfn(h0, c(nugget, th[2], th[3]))
        t(res) %*% solve(cov, res)
      }
      th0 = th
      th[1] = optimize(minfuncnugget, lower = 0, upper = upperOptim[1])$minimum
      th[2] = optimize(minfuncsill, lower = 0.000000001, upper = upperOptim[2])$minimum
      th[3] = optimize(minfuncrange, lower = 0.000000001, upper = upperOptim[3])$minimum
      converged = sum(abs((th - th0)/th0)) < eps
      iter = iter + 1
      if (trace) 
        print(th)
      model$psill = c(th[1], th[2])
      model$range[2] = th[3]
      rm(cov)
    }
    if (th[3]/max(h0) > 0.99) warning("range parameter at search space boundary")
    if (!converged) {
      warning("no convergence, returning OLS solution")
      th = th.ols
      model$psill = c(th[1], th[2])
      model$range[2] = th[3]
    }
    if (plot) 
      plot(variogram(formula, data, cloud = TRUE, cutoff = cutoff), 
           model = model)
    else model
  } else {
    if (ignoreInitial) { init = c(mean(y)/2, median(h0)/4)
    } else { init = c(model$psill, model$range) }
    
    gamfn = function(h, th, m = as.character(model$model)) variogramLine(vgm(psill = th[1], m, range = th[2]), dist_vector = h)$gamma
    minfuncols = function(theta = rbind(1, 1)) {
      res = y - gamfn(h0, theta)
      sum(res^2)
    }
    
    if (any(model$model == 'Pow')) {
      # We use 2 as upper bound for range on power models
      upperOptim <- c(max(y), 2)
      init[2] <- 1
    } else { upperOptim = c(max(y), max(h0)) }
    
    th = th0 = th.ols = optim(init, minfuncols, gr = NULL, method = "L-BFGS-B", 
                              lower = c(0.000000001, 0.000000001), upper = upperOptim)$par
    if (trace) print(th)
    while (!converged && iter < maxiter) {
      comb = function(i, j) cbind(rep(i, length(j)), rep(j, each = length(i)))
      cov = matrix(variogramLine(model, dist_vector = dists[comb(i, i)])$gamma + 
                     variogramLine(model, dist_vector = dists[comb(j, j)])$gamma - 
                     variogramLine(model, dist_vector = dists[comb(i, j)])$gamma - 
                     variogramLine(model, dist_vector = dists[comb(j, i)])$gamma, length(j), length(j))
      cov = 0.5 * cov^2
      cov = qr(cov)
      minfuncrange = function(range) {
        res = y - gamfn(h0, c(th[1], range))
        t(res) %*% solve(cov, res)
      }
      minfuncsill = function(sill) {
        res = y - gamfn(h0, c(sill, th[2]))
        t(res) %*% solve(cov, res)
      }
      th0 = th
      th[1] = optimize(minfuncsill, lower = 0.000000001, upper = upperOptim[1])$minimum
      th[2] = optimize(minfuncrange, lower = 0.000000001, upper = upperOptim[2])$minimum
      converged = sum(abs((th - th0)/th0)) < eps
      iter = iter + 1
      if (trace) print(th)
      model$psill = th[1]
      model$range = th[2]
      rm(cov)
    }
    if (th[2]/max(h0) > 0.99) warning("range parameter at search space boundary")
    if (!converged) {
      warning("no convergence, returning OLS solution")
      th = th.ols
      model$psill = th[1]
      model$range = th[2]
    }
    if (plot) 
      plot(variogram(formula, data, cloud = TRUE, cutoff = cutoff), 
           model = model)
    else model
  }
}

afvGLS <- function(formula, input_data, model, cutoff=Inf, verbose=FALSE, useNugget=TRUE) {
  # verbose <- T
  vc <- variogram(formula, input_data, cloud=T, cutoff=cutoff)

  withWarnings <- function(expr) {
    myWarnings <- NULL
    wHandler <- function(w) {
      myWarnings <<- c(myWarnings, list(w))
      invokeRestart("muffleWarning")
    }
    val <- withCallingHandlers(expr, warning = wHandler)
    list(value = val, warnings = myWarnings)
  }
  
  # i <- 1
  fitModels <- vector(mode = "list", length = length(model))
  mses <- rep(NA_real_, length(model))
  
  limites <- getBoundariesPVariogramaEmpiricoV8(fml=formula, observaciones=input_data, cutoff=cutoff)
  if (useNugget) { fixNugget <- NA
  } else { fixNugget <- 0 }

  for (i in 1:length(model)) {
    # i <- 1
    if (verbose) print(paste0(i, ': ', model[[i]]))
    vgIni <- afvmod(
      formula=formula, input_data=input_data, model=model[[i]], boundaries=limites, 
      miscFitOptions=list(orig.behavior=F), fix.values=c(fixNugget, NA, NA), verbose=verbose, 
      nPuntosIniciales=3, fit.method = 7)$var_model
    
    if (!is.null(vgIni)) {
      try({
        res <- withWarnings(expr = { 
          vg <- fit.variogram.gls_mod(
            formula = formula, data = input_data, trace = verbose, model = vgIni, ignoreInitial = F, 
            maxiter = 100) })
        
        # res <- withWarnings(expr = { vg <- fit.variogram.gls(formula = formula, data = input_data, trace = verbose, model = vgIni, ignoreInitial = T, maxiter = 100) })
        #if (length(res$warnings) == 0) {
        # plot(x=vc$dist, y=(variogramLine(object = vg, dist_vector = vc$dist)$gamma))
        fitModels[[i]] <- res$value
        mses[i] <- mean((variogramLine(object = vg, dist_vector = vc$dist)$gamma - vc$gamma)^2)
        
        if (verbose) {
          print(fitModels[[i]])
          print(mses[i])
        }
      })      
    }
  }
  
  if (any(!is.na(mses))) {
    bestMSE <- min(mses, na.rm = T)
    maxGamma <- max(vc$gamma)
    minPsill <- maxGamma * 0.001
    minRange <- 0.05
    
    # Considero todos los que tengan un MSE de como máximo 110% del mínimo
    iesAConsiderar <- !is.na(mses) & mses / bestMSE < 1.1 &
      !sapply(X = fitModels, FUN = isInvalidVariogram, minPsill=minPsill, minRange=minRange)
    fitModels <- fitModels[iesAConsiderar]
    mses <- mses[iesAConsiderar]
    
    nuggets <- sapply(fitModels, FUN=function(x) { 
      iNugget <- which(x$model == 'Nug')
      if (length(iNugget) > 0) { return(x$range[iNugget])
      } else { return (0) }})
    
    #getExperimentalRange(variogramModel = fitModels[[3]], maxDist = max(vc$dist), n=100)
    rangos <- sapply(fitModels, FUN = getExperimentalRange, maxDist = max(vc$dist))
    # sills <- sapply(fitModels, FUN = getExperimentalSill, maxDist = max(vc$dist))
    # Elijo el que tenga el mínimo nugget, si hay empate, elijo el de máximo rango
    iMejor <- order(nuggets, rangos, decreasing = c(F, T))[1]
    
    result = list(exp_var = vc, var_model = fitModels[[iMejor]], sserr = mses[iMejor])
    class(result) = c("autofitVariogram", "list")
  } else {
    result <- NULL
  }

  return(result)
}

afvGLSV2 <- function(formula, input_data, model, cutoff=NA, verbose=FALSE, useNugget=TRUE) {
  # verbose <- T
  if (!is.na(cutoff)) { vc <- variogram(formula, input_data, cloud=T, cutoff=cutoff)
  } else { vc <- variogram(formula, input_data, cloud=T) }
  
  withWarnings <- function(expr) {
    myWarnings <- NULL
    wHandler <- function(w) {
      myWarnings <<- c(myWarnings, list(w))
      invokeRestart("muffleWarning")
    }
    val <- withCallingHandlers(expr, warning = wHandler)
    list(value = val, warnings = myWarnings)
  }
  
  # i <- 3
  fitModels <- vector(mode = "list", length = length(model))
  mses <- rep(NA, length(model))
  
  for (i in 1:length(model)) {
    if (verbose) print(paste0(i, ': ', model[[i]]))
    if (useNugget) { vgIni <- vgm(psill = 1, model = model[[i]], range = 1, nugget = 0.01)
    } else { vgIni <- vgm(psill = 1, model = model[[i]], range = 1) }
    
    try({
      res <- withWarnings(expr = { vg <- fit.variogram.gls_mod(formula = formula, data = input_data, trace = verbose, model = vgIni, ignoreInitial = T, maxiter = 100) })
      # res <- withWarnings(expr = { vg <- fit.variogram.gls(formula = formula, data = input_data, trace = verbose, model = vgIni, ignoreInitial = T, maxiter = 100) })
      #if (length(res$warnings) == 0) {
      # plot(x=vc$dist, y=(variogramLine(object = vg, dist_vector = vc$dist)$gamma))
      fitModels[[i]] <- res$value
      mses[i] <- mean((variogramLine(object = vg, dist_vector = vc$dist)$gamma - vc$gamma)^2)
      
      if (verbose) {
        print(fitModels[[i]])
        print(mses[i])
      }
      #}
    })
  }
  
  if (any(!is.na(mses))) {
    bestMSE <- min(mses, na.rm = T)
    maxGamma <- max(vc$gamma)
    minPsill <- maxGamma * 0.001
    minRange <- 0.05
    
    # Considero todos los que tengan un MSE de como máximo 110% del mínimo
    iesAConsiderar <- !is.na(mses) & mses / bestMSE < 1.1 &
      !sapply(X = fitModels, FUN = isInvalidVariogram, minPsill=minPsill, minRange=minRange)
    fitModels <- fitModels[iesAConsiderar]
    mses <- mses[iesAConsiderar]
    
    nuggets <- sapply(fitModels, FUN=function(x) { 
      iNugget <- which(x$model == 'Nug')
      if (length(iNugget) > 0) { return(x$range[iNugget])
      } else { return (0) }})
    
    #getExperimentalRange(variogramModel = fitModels[[3]], maxDist = max(vc$dist), n=100)
    rangos <- sapply(fitModels, FUN = getExperimentalRange, maxDist = max(vc$dist))
    # sills <- sapply(fitModels, FUN = getExperimentalSill, maxDist = max(vc$dist))
    # Elijo el que tenga el mínimo nugget, si hay empate, elijo el de máximo rango
    iMejor <- order(nuggets, rangos, decreasing = c(F, T))[1]
    
    result = list(exp_var = vc, var_model = fitModels[[iMejor]], sserr = mses[iMejor])
    class(result) = c("autofitVariogram", "list")
  } else {
    result <- NULL
  }
  
  return(result)
}
