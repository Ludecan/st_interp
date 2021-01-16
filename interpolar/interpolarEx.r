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

# Busca en el stack el path relativo a getwd() del script para poder hacer los source correctamente
# Intenta con el frame actual - 3 primero que es donde estaba siempre cuando se hizo el programa
iFrame <- sys.nframe()
if (iFrame >= 3) { script.dir.interpolarEx <- sys.frame(iFrame - 3)$ofile 
} else { script.dir.interpolarEx <- NULL }
while ((is.null(script.dir.interpolarEx) || is.na(regexpr('interpolarEx.r', script.dir.interpolarEx, fixed=T)[1])) && iFrame >= 0) {
  script.dir.interpolarEx <- sys.frame(iFrame)$ofile
  iFrame <- iFrame - 1
}
if (is.null(script.dir.interpolarEx)) { script.dir.interpolarEx <- ''
} else { script.dir.interpolarEx <- paste0(dirname(script.dir.interpolarEx), '/') }

source(paste0(script.dir.interpolarEx, '../instalarPaquetes/instant_pkgs.r'), encoding = 'WINDOWS-1252')
source(paste0(script.dir.interpolarEx, '../pathUtils/pathUtils.r'), encoding = 'WINDOWS-1252')
source(paste0(script.dir.interpolarEx, '../IO/uIOMatrices.r'), encoding = 'WINDOWS-1252')
source(paste0(script.dir.interpolarEx, '../TryUtils/tryUtils.r'), encoding = 'WINDOWS-1252')
source(paste0(script.dir.interpolarEx, '../cacheFunciones/cacheFunciones.r'), encoding = 'WINDOWS-1252')
source(paste0(script.dir.interpolarEx, './mapearEx.r'), encoding = 'WINDOWS-1252')
source(paste0(script.dir.interpolarEx, 'funcionesAuxiliares.r'), encoding = 'WINDOWS-1252')
source(paste0(script.dir.interpolarEx, 'parsearParamsInterpolarYMapear.r'), encoding = 'WINDOWS-1252')
source(paste0(script.dir.interpolarEx, '../agregacion/agregacion.r'), encoding = 'WINDOWS-1252')
source(paste0(script.dir.interpolarEx, '../sysutils/sysutils.r'), encoding = 'WINDOWS-1252')

instant_pkgs(pkgs = c('unmarked', 'VGAM', 'cli', 'devtools'), silent = TRUE, doCargarPaquetes=FALSE)
instant_pkgs(
  pkgs = c('sp', 'digest', 'rgdal', 'parallel', 'doParallel', 'iterators', 'MASS', 'hash', 'Rcpp', 
           'raster', 'fields', 'xts', 'spacetime', 'lattice', 'numDeriv', 'Rmisc', 'nlme', 'glmnet', 
           'rms', 'leaps', 'AICcmodavg', 'zoo', 'FNN', 'gtools', 'gstat', 'automap', 'evd', 
           'htmltools', 'httr', 'stats', 'float', 'intamap', 'pROC'), 
  silent = TRUE)
library(unmarked)
library(VGAM)
library(cli)
library(devtools)
library(gridExtra)
library(sp)
library(digest)
library(rgdal)
library(parallel)
library(doParallel)
library(iterators)
library(MASS)
library(hash)
library(Rcpp)
library(raster)
library(fields)
library(xts)
library(spacetime)
library(lattice)
library(numDeriv)
library(Rmisc)
library(nlme)
library(glmnet)
library(rms)
library(leaps)
library(AICcmodavg)
library(zoo)
library(FNN)
library(gtools)
library(gstat)
library(automap)
library(evd)
library(htmltools)
library(httr)
library(stats)
library(float)
library(intamap)
library(pROC)

# instant_pkgs_github(reposgithub = 'Ludecan/intamap', minVersions = '1.4-4', silent = TRUE)
# instant_pkgs_github(reposgithub = 'jskoien/intamap', minVersions = '1.4-6', silent = TRUE)

formulaStr <- function(coeficientes, nDecimales=1, quitarCeros=TRUE) {
  if (quitarCeros) coeficientes <- coeficientes[coeficientes != 0 | names(coeficientes) == '(Intercept)']
  return(paste0("y ~ ", 
                paste(sprintf(paste0("%+.", nDecimales,"f*%s "), coeficientes[-1], names(coeficientes[-1])), collapse=""), 
                sprintf(paste0("%+.", nDecimales, "f"), coeficientes[1])))
}

formulaConCoeficientes <- function(modelo, nDecimales=1, quitarCeros=TRUE) {
  return(formulaStr(coefficients(modelo), nDecimales = nDecimales, quitarCeros = quitarCeros))
}

distKmToP4Str <- function(p4str, distKm) {
  if (grepl(pattern = "+proj=longlat", x = p4str, fixed = T)) {
    return(distKm)
  } else {
    unidad <- gsub(pattern = '+units=', fixed=T, replacement = '', 
                   x = regmatches(p4str, regexpr(pattern = '[+]units=[[:alpha:]]+', p4str)))
    
    if (unidad == 'm') { return(distKm * 1000)
    } else if (unidad =='km') { return(distKm)
    } else stop(paste0('interpolarEx.distKmToP4Str: Unidad de distancia no implementada "', unidad, '"'))
  }
}

crearSpatialPointsDataFrame <- function(
    x, y, value, proj4string='+proj=longlat +datum=WGS84', SRS_string, nombresFilas=NULL) {
  # crea un spatialPointsDataFrame a partir de las coordenadas x e y y una columna de valor value
  # por defecto asume que el CRS es lat/long y WGS84 pero se le puede especificar otro si es el caso
  # observaciones <- data.frame(x=x, y=y, value=value)
  if (is.character(x) & is.character(y)) {
    observaciones <- value
    sp::coordinates(observaciones) <- c(x, y)
  } else {
    observaciones <- cbind(data.frame(x=x, y=y), value)
    sp::coordinates(observaciones) <- c('x', 'y')
  }
  
  if (!is.null(nombresFilas)) row.names(observaciones) <- nombresFilas
  proj4string(observaciones) <- sp::CRS(projargs = proj4string, SRS_string=SRS_string)
  return (observaciones)
}

imitarObjetoIntamap <- function(observaciones, formulaString=value~1, predictions, intCRS=proj4string(observaciones), 
                                targetCRS=proj4string(predictions), class='automap', outputWhat=list(mean=T, variance=ncol(predictions@data) > 1)) {
  res <- createIntamapObject(observations=observaciones, formulaString=formulaString, predictionLocations=geometry(predictions),
                             intCRS=intCRS, targetCRS=targetCRS, class=class, outputWhat=outputWhat)
  res$campoMedia <- names(predictions@data)[1]
  if (ncol(predictions@data) > 1) { res$campoVarianza <- names(predictions@data)[2]
  } else { res$campoVarianza <- '' }
  res$predictions <- predictions
  return(res)
}

cargarCoordenadas <- function(
    pathArchivo, proj4string='+proj=longlat +datum=WGS84', SRS_string, 
    incluyeNombresObservaciones=FALSE) {
  if (incluyeNombresObservaciones) {
    observaciones <- read.table(file=pathArchivo, sep=' ', dec='.', header=T,
                                colClasses=c('character', 'numeric', 'numeric') , comment.char='')
    nombresFilas <- rownames(observaciones)
  } else {
    observaciones <- read.table(file=pathArchivo, sep=' ', dec='.', header=T,
                                colClasses=c('numeric', 'numeric'), comment.char='')
    nombresFilas <- NULL
  }
  
  return (crearSpatialPointsDataFrame(
    x=observaciones[,1], y=observaciones[,2], value=rep(x=NA, nrow(observaciones)), 
    proj4string=proj4string, SRS_string=SRS_string, nombresFilas=nombresFilas))
}

cargarObservaciones <- function(
    pathArchivo, proj4string='+proj=longlat +datum=WGS84', SRS_string, 
    incluyeNombresObservaciones=FALSE) {
  if (incluyeNombresObservaciones) {
    observaciones <- read.table(
      file=pathArchivo, sep=' ', dec='.', header=T, 
      colClasses=c('character', 'numeric', 'numeric', 'numeric'), comment.char='')
    nombresFilas <- rownames(observaciones)
  } else {
    observaciones <- read.table(file=pathArchivo, sep=' ', dec='.', header=T,
                                colClasses=c('numeric', 'numeric', 'numeric'), comment.char='')
    nombresFilas <- NULL
  }
  
  return (crearSpatialPointsDataFrame(
    x=observaciones[,2], y=observaciones[,3], value=observaciones[,1], proj4string=proj4string, 
    SRS_string=SRS_string, nombresFilas=nombresFilas))
}

cargarVectorDeBinario <- function(pathArchivo, NAValue=-.Machine$double.xmax, what='numeric', recordSize=8) {
  if (file.exists(pathArchivo)) {
    n <- file.info(pathArchivo)$size / recordSize
    vector <- readBin(pathArchivo, what=what, size=recordSize, n)
    vector[vector == NAValue] <- NA
    return (vector)
  } else {
    return (NULL)
  }
}

cargarMatrizDeBinario <- function(pathArchivo, nFilas, NAValue=-.Machine$double.xmax, what='numeric', recordSize=8) {
  if (file.exists(pathArchivo)) {
    n <- file.info(pathArchivo)$size / recordSize
    matriz <- matrix(readBin(pathArchivo, what=what, size=recordSize, n), nrow=nFilas, byrow=T)
    matriz[matriz == NAValue] <- NA
    return (matriz)
  } else {
    return (NULL)
  }
}

crearObservacionesBinarias <- function(observaciones, zcol=1) {
  return (SpatialPointsDataFrame(
    geometry(observaciones), 
    data=data.frame(value=as.numeric(observaciones@data[,zcol] >= 1E-3))))
}

crearCoordsAInterpolar <- function(
    xs, ys, grid=T, proj4string='+proj=longlat +datum=WGS84', SRS_string) {  
  # crea un spatialPoints o spatialPixels a partir de las coordenadas xs e ys
  # si grid=T es un spatialPixels y xs e ys son los ejes de la grilla
  # si grid=F es un spatialPoints y xs e ys son las coordenadas de los puntos
  # por defecto asume que el CRS es lat/long y WGS84 pero se le puede especificar otro si es el caso
  if (grid) { coordsAInterpolar <- expand.grid(x=xs, y=ys)
  } else { coordsAInterpolar <- data.frame(x=xs, y=ys) }
  sp::coordinates(coordsAInterpolar) <- c('x','y')
  projAInterpolar <- sp::CRS(projargs=proj4string, SRS_string=SRS_string)
  proj4string(coordsAInterpolar) <- projAInterpolar
  gridded(coordsAInterpolar) <- grid
  return (coordsAInterpolar)
}

crearGrillaAuxiliar <- function(
    xs, ys, gridValues, proj4string='+proj=longlat +datum=WGS84', SRS_string) {
  # crea un spatialPointsDataFrame a partir de las coordenadas x e y y una columna de valor value
  # por defecto asume que el CRS es lat/long y WGS84 pero se le puede especificar otro si es el caso
  grilla <- expand.grid(x=xs, y=ys)
  grilla$gridValue <- gridValues
  sp::coordinates(grilla) <- c('x', 'y')
  gridded(grilla) <- T
  proj4string(grilla) <- sp::CRS(projargs=proj4string, SRS_string=SRS_string)
  return (grilla)
}

crearGrillaRectilineaParaArea <- function(
    xMin, xMax, yMin, yMax, deltaX=(xMax-xMin)/250, deltaY=deltaX, 
    proj4stringCoordenadasArea='+proj=longlat +datum=WGS84', SRS_stringCoordenadasArea, 
    proj4stringGrillaSalida, SRS_stringGrillaSalida) {
  
  if (proj4stringCoordenadasArea != proj4stringGrillaSalida || 
      SRS_stringCoordenadasArea != SRS_stringGrillaSalida) {
    # Creamos la caja en la proyección de entrada, proyectamos a la proyección de salida
    # y obtenemos el área en la proyección de salida que engloba las 4 esquinas de la caja en
    # la proyección de entrada
    boundingBox <- matrix(c(xMin, yMin, xMax, yMin, xMin, yMax, xMax, yMax), ncol = 2, byrow = T)
    boundingBox <- SpatialPoints(
      coords=boundingBox, 
      proj4string=CRS(projargs=proj4stringCoordenadasArea, SRS_string=SRS_stringCoordenadasArea))
    boundingBox <- spTransform(
      boundingBox, CRS(projargs=proj4stringGrillaSalida, SRS_stringCoordenadasArea))
    
    xMin <- min(sp::coordinates(boundingBox)[,1])
    xMax <- max(sp::coordinates(boundingBox)[,1])
    yMin <- min(sp::coordinates(boundingBox)[,2])
    yMax <- max(sp::coordinates(boundingBox)[,2])
  }
  
  xs <- seq(from=xMin, to=xMax, by = deltaX)
  ys <- seq(from=yMin, to=yMax, by = deltaY)
  return(crearCoordsAInterpolar(
    xs = xs, ys = ys, grid = T, proj4string = proj4stringGrillaSalida, 
    SRS_string = SRS_stringGrillaSalida))
}

netCDFToSP <- function(
    fname, varName='rfe', p4string="+proj=longlat +datum=WGS84", SRS_string="EPSG:4326", 
    lonDimName='Lon', latDimName='Lat', spObj=NULL, zcol=1) {
  # TO-DO: handle multiple layer/times files
  nc <- nc_open(filename = fname)
  value <- as.vector(ncvar_get(nc, varid = varName))
  
  if (is.null(spObj)) {
    lon <- ncvar_get(nc, lonDimName)
    lat <- ncvar_get(nc, latDimName)
    coords <- as.matrix(expand.grid(lon=lon, lat=lat))
    spObj <- SpatialPixelsDataFrame(
      points = SpatialPoints(
        coords = coords, proj4string = CRS(projargs = p4string, SRS_string = SRS_string)), 
      data = data.frame(value), tolerance = 0.000152602)
    gridded(spObj) <- TRUE
    # Workaround for krigeST bug
    row.names(spObj) <- as.character(row.names(spObj))
  } else { spObj@data[,zcol] <- value }
  
  nc_close(nc)
  return(spObj)
}

extraerValorRegresorSobreSP <- function(i, objSP, pathsRegresor, fn=NULL, zcol=1, silent=T, ...) {
  # i <- 1
  if (!silent) print(i)

  if (!is.na(pathsRegresor[i]) && file.exists(pathsRegresor[i]) && length(objSP) > 0) {
    ext <- getFileExt(pathsRegresor[i])
    if (ext == 'tif') { evaluarConReintentos(regresor <- readGDAL(pathsRegresor[i], silent = silent))
    } else if (ext == 'nc') {
      # TO-DO: receive varName as parameter
      evaluarConReintentos(regresor <- netCDFToSP(fname = pathsRegresor[i]))
    } else { stop(paste0('extraerValorRegresorSobreSP: extensión no soportada "', pathsRegresor[i], '"')) }

    # Obtengo los valores del regresor en las coordenadas de las observaciones
    if (!identicalCRS(x = objSP, y = regresor)) {
      auxObjSP <- spTransform(x = geometry(objSP), CRSobj = regresor@proj4string)
      res <- over(x = auxObjSP, y = regresor, fn = fn, returnList = F, ...=...)[,zcol]
    } else {  res <- over(x = objSP, y = regresor, fn = fn, returnList = F, ...=...)[,zcol]
    }
  } else {
    res <- numeric(length(objSP))
    res[] <- NA
  }
  return(res)
}

extraerValoresRegresorSobreSP <- function(
    objSP, pathsRegresor, iInicial = 1, iFinal = length(pathsRegresor), fn=NULL, zcol=1, silent=T, 
    nCoresAUsar=0, setNames=T, ...) {
  pr <- pathsRegresor[seq.int(from = iInicial, to = iFinal, by = 1)]
  pathsUnicos <- unique(pr)
  
  # netCDF is not reentrant. Force single core
  if (any(getFileExt(pathsUnicos) == 'nc', na.rm=T)) { nCoresAUsar <- 1
  } else if (nCoresAUsar <= 0) { nCoresAUsar <- min(detectCores(T, T), length(pathsUnicos)) }
  
  if (nCoresAUsar > 1) {
    cl <- makeCluster(getOption('cl.cores', nCoresAUsar))
    clusterExport(cl, varlist = c('script.dir.interpolarEx'))
    clusterEvalQ(cl = cl, {
      source(paste0(script.dir.interpolarEx, '../TryUtils/tryUtils.r'), encoding = 'WINDOWS-1252')
      source(paste0(script.dir.interpolarEx, '../pathUtils/pathUtils.r'), encoding = 'WINDOWS-1252')
      require(rgdal)
      require(sp)
    })
    if (exists(x = 'setMKLthreads')) { clusterEvalQ(cl = cl, expr = setMKLthreads(1)) }
    valoresSobreSP <- parSapplyLB(
      cl = cl, X=1:length(pathsUnicos), FUN = extraerValorRegresorSobreSP, objSP=objSP, 
      pathsRegresor=pathsUnicos, fn=fn, zcol=zcol, silent=silent, ...=...)
    stopCluster(cl)
  } else {
    valoresSobreSP <- sapply(
      X=1:length(pathsUnicos), FUN = extraerValorRegresorSobreSP, objSP=objSP, 
      pathsRegresor=pathsUnicos, fn=fn, zcol=zcol, silent=silent, ...=...)
  }
  
  if (length(pathsUnicos) < length(pr)) { 
    iMatch <- match(x = pr, pathsUnicos)
    valoresSobreSP <- t(valoresSobreSP[, iMatch])
  } else {
    valoresSobreSP <- t(valoresSobreSP)
  }
  
  if (setNames) {
    colnames(valoresSobreSP) <- row.names(objSP)
    
    if (!is.null(row.names(pathsRegresor))) {
      rownames(valoresSobreSP) <- rownames(pathsRegresor)
    } else if (!is.null(names(pathsRegresor))) {
      rownames(valoresSobreSP) <- names(pathsRegresor)
    } else {
      warning(paste0(
        'interpolarEx.extraerValoresRegresorSobreSP: using setNames == TRUE but both row.names(',
        'pathsRegresor) and names(pathsRegresor) are NULL. make sure to set them to have appropiate',
        ' result rownames'))
    }
  }
  
  return(valoresSobreSP)
}

extraerValoresRegresoresSobreSP <- function(
    objSP, pathsRegresores, iInicial = 1, iFinal = nrow(pathsRegresores), fn=NULL, zcol=1, silent=T, 
    nCoresAUsar=0, setNames=T, ...) {
  # TO-DO: paralelizar esta función. Hoy se está haciendo paralelo en el tiempo pero si solo se 
  # quiere cargar una fecha para más de un regresor se está serializando innecesariamente. Hay que 
  # hacer que decida si paralelizar por filas o columnas para aprovechar mejor los recursos
  res <- vector(mode = "list", ncol(pathsRegresores))

  # i <- 1
  for (i in 1:ncol(pathsRegresores)) {
    res[[i]] <- extraerValoresRegresorSobreSP(
      objSP=objSP, pathsRegresor=pathsRegresores[, i,drop=F], iInicial=iInicial, iFinal=iFinal, 
      fn=fn, zcol = zcol, silent=silent, nCoresAUsar=nCoresAUsar, setNames = setNames, ...=...)
  }
  
  # Código para verificación
  # i <- 5
  # for (i in 1:nrow(pathsRegresores)) {
  #   pathsRegresores[i,]
  #   valoresRegresoresSobreObservaciones[[1]][i,]
  #   vr <- readGDAL(fname = pathsRegresores[i,], silent = T)
  #   auxSP <- spTransform(coordsObservaciones, vr@proj4string)
  #   print(i)
  #   print(max(abs(valoresRegresoresSobreObservaciones[[1]][i,] - over(auxSP, vr))))
  # }

  names(res) <- colnames(pathsRegresores)
  return(res)
}

getGridIndexes <- function(observaciones, grid) {
  # observaciones must be spatialPoints, grid must be spatialPixelsDataFrame
  # returns the indexes of the grid cells containing the coordinates in observations
  return (over(observaciones, grid))
}

getResiduals <- function(observaciones, grid) {
  gridIndexes <- getGridIndexes(observaciones, grid)
  return (grid[gridIndexes] - observaciones$value)
}

cargarSHP <- function(pathSHP, proj4strSHP=NULL, overrideP4str=FALSE, encoding=NULL) {
  if (file.exists(pathSHP)) {
    dsnSHP <- dirname(pathSHP)
    layerSHP <- nombreArchSinPathNiExtension(pathSHP)

    archPrj <- changeFileExt(pathSHP, '.prj')
    if (!overrideP4str && file.exists(archPrj)) { proj4strSHP <- NULL }
    return (readOGR(dsn=dsnSHP, layer=layerSHP, p4s=proj4strSHP, encoding=encoding))
  } else {
    return (NULL)
  }
}

guardarSHP <- function(shp, path) {
  dsnSHP <- dirname(path)
  layerSHP <- nombreArchSinPathNiExtension(path)
  writeOGR(obj = shp, dsn = dsnSHP, layer = layerSHP, driver = "ESRI Shapefile", overwrite_layer = T)
}

cargarSHPYObtenerMascaraParaGrilla <- function(
    pathSHP, proj4strSHP=NULL, grilla, spSinMascara=NULL, overrideP4str=FALSE, encoding = NULL) {
  if (file.exists(pathSHP)) {
    grilla <- geometry(grilla)
    if (!is.null(spSinMascara)) {
      spSinMascara <- geometry(spSinMascara)
      
      if (!identicalCRS(grilla, spSinMascara)) {
        spSinMascara <- spTransform(spSinMascara, grilla@proj4string)
      }
    }
    
    objCache <- list(pathSHP, proj4strSHP, grilla, spSinMascara, overrideP4str, encoding)
    pathCache <- getPathCache(
      objCache, dirEjecucion = paste0(dirname(pathSHP), '/'), 
      prefijoNombreArchivoCache = paste0(nombreArchSinPathNiExtension(pathSHP), '_'))
    if (!file.exists(pathCache)) {
      shp <- cargarSHP(pathSHP, proj4strSHP, overrideP4str = overrideP4str, encoding = encoding)
      if (!identicalCRS(grilla, shp)) { shp <- spTransform(shp, grilla@proj4string) }
      
      # Recorte al contorno del país
      if (!is.null(grilla)) {
        mask <- !is.na(over(grilla, geometry(shp)))
        if (!is.null(spSinMascara)) { mask[over(spSinMascara, grilla)] <- TRUE }
      } else { mask <- TRUE }

      #aux <- SpatialPixelsDataFrame(points = grilla, data = data.frame(value=as.integer(mask)))
      #spplot(aux, zcol=1)
      shpMask <- list(shp=shp, mask=mask)
      crearDirectoriosSiNoExisten(dirname(pathCache))
      guardarCache(pathCache = pathCache, obj = shpMask)
    } else {
      shpMask <- cargarCache(pathCache = pathCache)
    }
  } else {
    if (!is.null(grilla)) { mask <- rep(TRUE, length(grilla))
    } else { mask <- TRUE }
    shpMask <- list(shp=NULL, mask=mask)
  }
  
  return (shpMask)
}

setMinMaxVal <- function(observacionesValue, params) {
  if (is.null(params$mLimitarValoresInterpolados)) {
    if (is.null(params$minVal)) params$minVal <- NA
    if (is.null(params$maxVal)) params$maxVal <- NA
  } else if (params$mLimitarValoresInterpolados=='NoLimitar') {
    params$minVal <- NA
    params$maxVal <- NA
  } else if (params$mLimitarValoresInterpolados=='LimitarMinimo') {
    params$minVal <- params$minimoLVI
    params$maxVal <- NA
  } else if (params$mLimitarValoresInterpolados=='LimitarMaximo') {
    params$minVal <- NA
    params$maxVal <- params$maximoLVI
  } else if (params$mLimitarValoresInterpolados=='LimitarMinimoyMaximo') {
    params$minVal <- params$minimoLVI
    params$maxVal <- params$maximoLVI
  } else if (params$mLimitarValoresInterpolados=='UsarPromDesvEst') {
    media <- mean(observacionesValue, na.rm=T)
    desvEst <- sd(observacionesValue, na.rm=T)

    params$minVal <- media - params$factorDesvEstLVI * desvEst
    params$maxVal <- media + params$factorDesvEstLVI * desvEst
  } else if (params$mLimitarValoresInterpolados=='UsarPromDesvEstYMinimoYMaximo') {
    media <- mean(observacionesValue, na.rm=T)
    desvEst <- sd(observacionesValue, na.rm=T)

    params$minVal <- max(params$minimoLVI, media - params$factorDesvEstLVI * desvEst, na.rm=T)
    params$maxVal <- min(params$maximoLVI, media + params$factorDesvEstLVI * desvEst, na.rm=T)
  } else if (params$mLimitarValoresInterpolados=='UsarMedianaMAD') {
    mediana <- median(observacionesValue, na.rm=T)
    desvMedAbs <- mad(observacionesValue, na.rm=T)
    
    params$minVal <- mediana - params$factorDesvEstLVI * desvMedAbs
    params$maxVal <- mediana + params$factorDesvEstLVI * desvMedAbs
  } else if (params$mLimitarValoresInterpolados=='UsarMedianaMADYMinimoYMaximo') {
    mediana <- median(observacionesValue, na.rm=T)
    desvMedAbs <- mad(observacionesValue, na.rm=T)
    
    params$minVal <- max(params$minimoLVI, mediana - params$factorDesvEstLVI * desvMedAbs, na.rm=T)
    params$maxVal <- min(params$maximoLVI, mediana + params$factorDesvEstLVI * desvMedAbs, na.rm=T)
  } else {
    stop(paste0('interpolarEx.setMinMaxVal: método para limitar valores interpolados desconocido ',  params$mLimitarValoresInterpolados))
  }
  
  if (!is.na(params$minVal)) { params$minVal <- min(params$minVal, observacionesValue, na.rm=T)}
  if (!is.na(params$maxVal)) { params$maxVal <- max(params$maxVal, observacionesValue, na.rm=T)}
  
  return(params)
}

getDefaultSpatialCutoff <- function(coordsObservaciones, params) {
  dists <- spDists(coordsObservaciones, longlat=!is.projected(coordsObservaciones))  
  dists <- dists[upper.tri(dists)]
  #return(max(dists/2.5))
  if ((params$usarFitVariogramGLS == 'auto' && length(coordsObservaciones) <= 50) || 
      (is.logical(params$usarFitVariogramGLS) && as.logical(params$usarFitVariogramGLS))) { 
    if (length(coordsObservaciones) <= 100) { p <- 2/3
    } else { p <- 0.1 }
  } else { p <- 0.75 }
  return(quantile(dists, probs = p))
}

getVentana <- function(ti, nT, tamanioSemiVentana, tlagsAR=NULL, tlags=NULL) {
  if (length(tlagsAR) > 0) {
    tIniVentana <- max(1, ti - 2 * tamanioSemiVentana)
    if (is.null(tlags)) { tFinVentana <- ti
    } else { tFinVentana <- max(ti, tlags, 2 * tamanioSemiVentana + 1) } 
  } else if (ti + tamanioSemiVentana <= nT) {
    tIniVentana <- max(1, ti - tamanioSemiVentana)
    tFinVentana <- min(nT, ti + tamanioSemiVentana + tIniVentana - (ti - tamanioSemiVentana))
  } else {
    tFinVentana <- nT
    tIniVentana <- max(1, ti - (tamanioSemiVentana + (ti + tamanioSemiVentana - tFinVentana)))
  }
  
  tsVentana <- seq.int(from=tIniVentana, to=tFinVentana, by = 1)
  iTiEnTsVentana <- which(tsVentana == ti)
  return(list(tsVentana=tsVentana, iTiEnTsVentana=iTiEnTsVentana))
}

prepIDWEstimator <- function(interpolacion, params, objParameters) {
  interpolacion$campoMedia <- 'var1.pred'
  interpolacion$campoVarianza <- NULL
  if (is.null(objParameters)) {
    if (!is.na(params$inverseDistancePower)) { interpolacion$inverseDistancePower <- params$inverseDistancePower
    } else { 
      # estimateParameters for idw uses n-fold cross validation which uses random numbers for 
      # sampling. Either set the random seed to achieve consistent results or use 
      # nfolds=length(observaciones) for LOOCV
      # set.seed(31)
      interpolacion <- estimateParameters(
        interpolacion, idpRange=seq(1.5, 4, 0.25), nfolds=length(interpolacion$observations))
      # source(paste0(script.dir.interpolarEx, 'afvmod.r'))
      # interpolacion <- estimateParameters.idw_mod(interpolacion, idpRange=seq(1.5, 4, 0.25), nfolds=length(observaciones))
    }
  } else {
    interpolacion$inverseDistancePower <- objParameters$inverseDistancePower
  }
  return(interpolacion)
}

interpolarEx <- function(
    observaciones, coordsAInterpolar, params, shpMask=NULL, longitudesEnColumnas=T, 
    objParameters=NULL, valoresCampoBaseSobreObservaciones=NULL, valoresCampoBase=NULL) {
  proj4StringAInterpolar <- proj4string(coordsAInterpolar)

  # Saco las observaciones NA y las coordenadas duplicadas
  if (!is.null(valoresCampoBaseSobreObservaciones)) {
    i <- !is.na(observaciones$value) & !is.na(valoresCampoBaseSobreObservaciones)
    observaciones <- observaciones[i, ]
    valoresCampoBaseSobreObservaciones <- valoresCampoBaseSobreObservaciones[i]
    if (!is.null(params$valoresRegresoresSobreObservaciones)) {
      params$valoresRegresoresSobreObservaciones <- lapply(
        params$valoresRegresoresSobreObservaciones, function(x) return(x[, i, drop=F]))
    }
    i <- !duplicated(sp::coordinates(observaciones))
    observaciones <- observaciones[i, ]
    valoresCampoBaseSobreObservaciones <- valoresCampoBaseSobreObservaciones[i]
    if (!is.null(params$valoresRegresoresSobreObservaciones)) {
      params$valoresRegresoresSobreObservaciones <- lapply(
        params$valoresRegresoresSobreObservaciones, function(x) return(x[, i, drop=F]))
    }
  } else {
    i <- !is.na(observaciones$value)
    observaciones <- observaciones[i,]
    if (!is.null(params$valoresRegresoresSobreObservaciones)) {
      params$valoresRegresoresSobreObservaciones <- lapply(
        params$valoresRegresoresSobreObservaciones, function(x) return(x[, i, drop=F]))
    }
    i <- !duplicated(sp::coordinates(observaciones))
    observaciones <- observaciones[i,]
    if (!is.null(params$valoresRegresoresSobreObservaciones)) {
      params$valoresRegresoresSobreObservaciones <- lapply(
        params$valoresRegresoresSobreObservaciones, function(x) return(x[, i, drop=F]))
    }
  }
  
  # Shapefile con el contorno del país / Shapefile for the country boundaries
  if (is.null(shpMask)) {
    shpMask <- cargarSHPYObtenerMascaraParaGrilla(
      pathSHP=params$pathSHPMapaBase, grilla=coordsAInterpolar)
  }
  
  params <- setMinMaxVal(observacionesValue=observaciones$value, params)
  if (!is.null(valoresCampoBase) & !is.null(valoresCampoBaseSobreObservaciones)) {
    observaciones$value <- observaciones$value - valoresCampoBaseSobreObservaciones
    
    if (params$metodoIgualacionDistribuciones == 'regresionLinealRobusta') {
      # Si se hizo una regresión lineal robusta tengo que sacar los residuos extremos que deja detrás porque afectan
      # fuertemente a la interpolación
      i <- filtroMedianaMAD(observaciones$value)
      observaciones <- observaciones[i, ]
      valoresCampoBaseSobreObservaciones <- valoresCampoBaseSobreObservaciones[i]
    }
    
    if (params$simpleKrigingEnRK) {
      if (is.null(params$betaSimpleKriging)) {
        beta <- params$betaSimpleKriging
      } else {
        beta <- mean(observaciones$value)  
      }
      # Si es simple kriging no puedo usar el variograma potencia porque es para campos no estacionarios
      params$modelosVariograma <- params$modelosVariograma[params$modelosVariograma != 'Pow']
    } else {
      beta <- NULL
    }
  } else { 
    beta <- params$betaSimpleKriging
  }
  
  # Si la varianza es muy chica (todos los valores son iguales) el método de interpolación da un error.
  # Manejamos el caso aparte, según los datos observados el campo es constante
  # For small variances (all observations the same) the library gives an error, we handle the case separately
  varianzaObservaciones <- var(observaciones$value)
  if (varianzaObservaciones > 1E-6 && params$interpolationMethod !='none' && length(observaciones) >= 4) {
    if (!identicalCRS(observaciones, coordsAInterpolar)) { 
      observaciones <- spTransform(observaciones, coordsAInterpolar@proj4string) 
    }
  
    mapaConstante <- F
    
    if (params$interpolationMethod == 'automatic') { 
      params$interpolationMethod <- seleccionarMetodoInterpolacion(observaciones$value) 
    }
    
    # exijo una cantidad mínima de puntos pq sino no vale la pena ni preparar los threads
    # Minimum amount of work for multicoreing
    if (params$nCoresAUsar <= 0) { if (length(coordsAInterpolar) > 100) { nCoresAUsar <- detectCores(T, T) } else { nCoresAUsar <- 1 }
    } else { nCoresAUsar <- params$nCoresAUsar }
    
    # crear objeto intamap/create intamap object
    # interpolation has 4 basic steps
    
    # interpolacion <- intamap::preProcess(interpolacion) 
    # interpolacion <- intamap::estimateParameters(interpolacion) and interpolacion <- estimateAnisotropy(interpolacion)
    # interpolacion <- intamap::spatialPredict(interpolacion)
    # interpolacion <- intamap::postProcess(interpolacion)
    
    # You can override any of the steps by adding the required fields to the interpolation object
    # We do this for the automap method to override variogram estimation

    if (params$interpolationMethod != 'idw') { outputWhat <- list(mean=T, variance=T) 
    } else { outputWhat <- list(mean=T) }
    
    interpolacion <- createIntamapObject(
      observations=observaciones, formulaString=value ~ 1, 
      predictionLocations=coordsAInterpolar[shpMask$mask, ],
      intCRS=proj4StringAInterpolar, targetCRS=proj4StringAInterpolar, 
      class=params$interpolationMethod, 
      params=list(nclus=nCoresAUsar, nmin=params$nmin, nmax=params$nmax, 
                  maxdist=params$maxdist, beta=beta, debug.level = 0),
      outputWhat=outputWhat)
    # mapearPuntosGGPlot(observaciones, shpMask$shp, dibujarTexto = length(observaciones) <= 50)
    # checkSetup(interpolacion)
    interpolacion <- intamap::preProcess(interpolacion)
    
    if (params$interpolationMethod == 'automap') {
      interpolacion$campoMedia <- 'var1.pred'
      interpolacion$campoVarianza <- 'var1.var'    
      
      if (is.null(objParameters)) {
        # Variograma
        source(paste0(script.dir.interpolarEx, 'afvmod.r'), encoding = 'WINDOWS-1252')
        require('gstat')
        
        if (params$imitarSurfer) {
          maxDist <- max(spDists(coordsObservaciones, longlat = FALSE))
          cutoff <- maxDist / 2.5
          limites <- seq(from = 0, to = cutoff, length.out = 25)
          variogramas <- afvmod(formula=interpolacion$formulaString, input_data=interpolacion$observations, 
                                model=c('Lin'), boundaries=limites, miscFitOptions=list(orig.behavior=F), 
                                fix.values=c(0, cutoff, NA), verbose=params$verbose, nPuntosIniciales=params$nPuntosIniciales,
                                tryFixNugget=params$tryFixNugget, fit.method = 7)
          variogramas$var_model$range <- maxDist
        } else {
          usarAFVGLS <- (params$usarFitVariogramGLS == 'auto' && length(interpolacion$observations) <= 50) || (is.logical(params$usarFitVariogramGLS) && as.logical(params$usarFitVariogramGLS))
          source(paste0(script.dir.interpolarEx, 'getBoundariesPVariogramaEmpirico.r'), encoding = 'WINDOWS-1252')
          
          if (usarAFVGLS) {
            variogramas <- afvGLS(
              formula=interpolacion$formulaString, input_data=interpolacion$observations, 
              cutoff = params$cutoff, model=params$modelosVariograma, verbose=params$verbose, 
              useNugget=params$usarNugget)
          }
          
          if (!usarAFVGLS || is.null(variogramas$var_model)) {
            # find the "best" amount onf intervals for the empirical variogram
            limites <- getBoundariesPVariogramaEmpiricoV8(fml=interpolacion$formulaString, observaciones=observaciones, cutoff=params$cutoff)
            if (params$usarNugget) { fixNugget <- NA
            } else { fixNugget <- 0 }
            
            #variogramas <- vector(length = length(params$modelosVariograma), mode = 'list')
            #for (i in seq_along(params$modelosVariograma)) {
            #  variogramas[[i]] <- afvmod(formula=interpolacion$formulaString, input_data=interpolacion$observations, 
            #                        model=params$modelosVariograma[i], boundaries=limites, miscFitOptions=list(orig.behavior=F), 
            #                        fix.values=c(fixNugget, NA, NA), verbose=params$verbose, nPuntosIniciales=params$nPuntosIniciales,
            #                        tryFixNugget=params$tryFixNugget, fit.method = 7)  
            #}
            
            variogramas <- afvmod(formula=interpolacion$formulaString, input_data=interpolacion$observations, 
                                  model=params$modelosVariograma, boundaries=limites, miscFitOptions=list(orig.behavior=F), 
                                  fix.values=c(fixNugget, NA, NA), verbose=params$verbose, nPuntosIniciales=params$nPuntosIniciales,
                                  tryFixNugget=params$tryFixNugget, fit.method = 7)
          }          
        }
        
        #png('Resultados/Ejemplos/Variogramas/VariogramaModelo.png', width = 800, height = 600)
        #print(annotatedplot(variogramas, xlab = "Distancia", ylab = "Semivarianza", main = "Variograma Empírico y Variograma Modelo Ajustado"))
        #dev.off()
        #variogramas$var_model <- NULL
        #png('Resultados/Ejemplos/Variogramas/VariogramaEmpirico.png', width = 800, height = 600)
        #print(annotatedplot(variogramas, xlab = "Distancia", ylab = "Semivarianza", main = "Variograma Empírico"))
        #dev.off()
        #variogramas$exp_var <- variogram(interpolacion$formulaString, data=interpolacion$observations, cloud=T, cutoff=params$cutoff)
        #png('Resultados/Ejemplos/Variogramas/NubeVariograma.png', width = 800, height = 600)
        #print(annotatedplot(variogramas, xlab = "Distancia", ylab = "Semivarianza", main = "Nube de variograma"))
        #dev.off()
        
        # annotatedplot(variogramas)
        # mapearPuntosGGPlot(observaciones, shpBase = shpMask$shp, continuo = T, dibujar=F)
        
        #if (!is.na(params$cutoff)) { vc <- variogram(interpolacion$formulaString, data=interpolacion$observations, cloud=T, cutoff=params$cutoff)
        #} else { vc <- variogram(interpolacion$formulaString, data=interpolacion$observations, cloud=T, cutoff=max(limites)) }
        #plot(vc$dist, vc$gamma)
        #sDist <- seq(from=0, to = max(vc$dist), length.out = 200)
        #lines(sDist, variogramLine(object = variogramas3$var_model, dist_vector = sDist)$gamma, col='red')
        #lines(sDist, variogramLine(object = variogramas2$var_model, dist_vector = sDist)$gamma, col='blue')
        #lines(sDist, variogramLine(object = variogramas$var_model, dist_vector = sDist)$gamma, col='black')
        #for (i in 1:length(limites)) {
        #  abline(v=limites[i])
        #}
        #variogramas$var_model <- vgm(psill = 1, range = max(variogramas$exp_var$dist), model = 'Exp')
        # interpolacion$sampleVariogram <- vc
        # interpolacion$variogramModel <- bestVG
        # fechasObservaciones[ti]
        if (params$modoDiagnostico) {
          png(paste0(params$carpetaParaModoDiagnostico, '05.2-Variograma.png'))
          print(annotatedplot(krigeobj = variogramas))
          dev.off()
        }
        interpolacion$sampleVariogram <- variogramas$exp_var
        interpolacion$variogramModel <- variogramas$var_model
        interpolacion$sserr <- variogramas$sserr
        
        if (is.null(interpolacion$variogramModel)) stop('interpolarEx: no se pudo ajustar ningún variograma a los datos')
        interpolacion <- estimateAnisotropy(interpolacion)
        
        # annotatedplot(variogramas)
        # cressie estimator, not giving good results, should investigate further
        # usarCressie <- is.unsorted(interpolacion$sampleVariogram$gamma[1:3]) || (interpolacion$sampleVariogram$gamma[1] > mean(interpolacion$sampleVariogram$gamma) * 0.5)
        # if (usarCressie) {
        # variogramas2 <- afvmod(interpolacion$formulaString,interpolacion$observations,model=modelos,boundaries=limites,miscFitOptions=lMiscFitOptions,cressie=T)
        # usarCressie <- !(is.unsorted(variogramas2$exp_var$gamma[1:3]) || (variogramas2$exp_var$gamma[1] > mean(variogramas2$exp_var$gamma) * 0.5))
        # # interpolacion2 <- estimateParameters(interpolacion,model=modelos,cressie=T)    
        # # usarCressie <- !(is.unsorted(interpolacion2$sampleVariogram$gamma[1:3]) || (interpolacion2$sampleVariogram$gamma[1] > mean(interpolacion2$sampleVariogram$gamma) * 0.5))
        # if (usarCressie) {
        #   interpolacion$sampleVariogram <- variogramas2$exp_var
        #   interpolacion$variogramModel <- variogramas2$var_model
        #   #interpolacion <- interpolacion2
        # }
        # # rm(interpolacion2)
        # rm(variogramas2)
        # } 
        #interpolacion <- estimateAnisotropy(interpolacion)
      } else {
        interpolacion$sampleVariogram <- objParameters$sampleVariogram
        interpolacion$variogramModel <- objParameters$variogramModel
        interpolacion$sserr <- objParameters$sserr
        interpolacion <- estimateAnisotropy(interpolacion)
      }
    } else if (params$interpolationMethod == 'idw') {
      interpolacion <- prepIDWEstimator(
        interpolacion=interpolacion, params=params, objParameters=objParameters)
    } else if (params$interpolationMethod == 'copula') {
      interpolacion$campoMedia <- 'mean'
      interpolacion$campoVarianza <- 'variance'
      
      # anisotropy <- list(lower = c(0, 1), upper = c(pi, Inf), params = c(pi/3, 2))
      # correlation <- list(model = 'Ste', lower = c(0.01, 0.01, 0.01), upper = c(0.99, Inf, 20),
      #                       params = c(0.05, 4, 3))
      # margin <- list(name = 'gev', lower = c(0.01, -Inf), upper = c(Inf, Inf), params = c(30, 0.5))
      # trend <- list(F = as.matrix(rep(1, 196)), lower = -Inf, upper = Inf, params = 40)  
      # copula <-list(method='chisq')
      
      # interpolacion <- estimateParameters(interpolacion, margin=margin, trend=trend, correlation=correlation, anisotropy=anisotropy, copula=copula, tol=0.001)
      # help(copulaEstimation)
      interpolacion <- estimateParameters(interpolacion)
      # warnings()
    }
    
    #delta <- c(-2.5, -1.5, -0.5, 0, 0.5, 1.5, 2.5)
    #block <- expand.grid(x=delta, y=delta)
    if (!is.na(params$block)) { interpolacion <- intamap::blockPredict(interpolacion, block=params$block)
    } else { interpolacion <- intamap::spatialPredict(interpolacion) }
    
    interpolacion <- postProcess(interpolacion)
    #observaciones@data[3,]$value <- NA
    #observaciones <- observaciones[!is.na(observaciones$value),]
    #interpolacion3 <- krige(value~1, observaciones, coordsAInterpolar, model=interpolacion$variogramModel)
    #interpolacion4 <- krige(value~1, observaciones, coordsAInterpolar, model=interpolacion$variogramModel, beta=0)
    #cbind(over(coordsObservaciones, interpolacion3)[1], over(coordsObservaciones, interpolacion4)[1])
    #observaciones$Nombre
    #intamap:::spatialPredict.automap()
    #interpolacion2 <- interpolacion
    #cbind(over(coordsObservaciones, interpolacion$predictions)[1],  over(coordsObservaciones, interpolacion2$predictions)[1])
    # interpolacion2 <- interpolacion
    # if (nCoresAUsar > 1) { closeAllConnections() }
  } else {
    # zero variance, constant map. Highly unlikely except for zero rainfall
    mapaConstante <- T

    nPredictionLocations <- length(coordsAInterpolar[shpMask$mask,])
    if (params$interpolationMethod != 'none') {
      if (length(observaciones) >= 1) {
        var1.pred <- rep(NA_real_, nPredictionLocations)
        if (!is.null(valoresCampoBase)) { var1.pred[!is.na(valoresCampoBase[shpMask$mask])] <- mean(observaciones$value)
        } else { var1.pred[] <- mean(observaciones$value) }
      } else if (!is.null(valoresCampoBase)) {
        var1.pred <- rep(0, times=nPredictionLocations)
      } else{ 
        var1.pred <- rep(NA_real_, nPredictionLocations) 
      }
      
      predictions <- data.frame(var1.pred=var1.pred, var1.var=rep(x=0, times=nPredictionLocations))
    } else if (!is.null(valoresCampoBase)) {
      predictions <- data.frame(var1.pred=rep(0, times=nPredictionLocations), var1.var=rep(x=0, times=nPredictionLocations))
    } else {
      predictions <- data.frame(var1.pred=rep(NA_real_, times=nPredictionLocations), var1.var=rep(x=NA_real_, times=nPredictionLocations))
    }
    
    sp::coordinates(predictions) <- sp::coordinates(coordsAInterpolar[shpMask$mask,])
    proj4string(predictions) <- proj4string(coordsAInterpolar)
    gridded(predictions) <- gridded(coordsAInterpolar[shpMask$mask,])
    class(predictions@data[, 1]) <- 'numeric'
    class(predictions@data[, 2]) <- 'numeric'
    
    interpolacion <- imitarObjetoIntamap(
      observaciones=observaciones, formulaString=value~1, predictions=predictions, 
      intCRS=proj4StringAInterpolar, targetCRS=proj4StringAInterpolar, class='automap', 
      outputWhat=list(mean=T, variance=T))
  }
  
  if (params$modoDiagnostico & (!is.null(valoresCampoBase) & !is.null(valoresCampoBaseSobreObservaciones))) {
    escalaResiduos <- crearEscalaEquiespaciada(
      datos = c(observaciones$value, interpolacion$predictions@data[, interpolacion$campoMedia]), 
      nDigitos = 2, continuo = T)
    # En modo diagnóstico si es RK guardo el mapa con los residuos y su interpolación 
    mapearPuntosGGPlot(puntos = observaciones, shpBase = shpMask$shp, continuo = T, dibujarTexto = T, escala = escalaResiduos,
                       nDigitos = 2, titulo = paste0('Residuos - ', params$strFecha), zcol='value', 
                       nomArchResultados = paste0(params$carpetaParaModoDiagnostico, '05.1-Residuos.png'), 
                       dibujar = F)
    
    if (params$interpolationMethod != 'none') {
      if (gridded(coordsAInterpolar)) {
        mapearGrillaGGPlot(grilla = interpolacion$predictions, shpBase = shpMask$shp, zcol=1, escala = escalaResiduos,
                           titulo = paste0('Interpolación - ', params$strFecha), 
                           nomArchResultados = paste0(params$carpetaParaModoDiagnostico, '05.3-Interpolacion.png'), 
                           dibujar = F, dibujarPuntosObservaciones = T, coordsObservaciones = observaciones)
      } else {
        mapearPuntosGGPlot(puntos = interpolacion$predictions, shpBase = shpMask$shp, zcol=1, escala = escalaResiduos,
                           titulo = paste0('Interpolación - ', params$strFecha), 
                           nomArchResultados = paste0(params$carpetaParaModoDiagnostico, '05.3-Interpolacion.png'), 
                           dibujar = F, dibujarTexto = T)
      }
    }
    rm(escalaResiduos)
  }  
  
  if (!is.null(valoresCampoBase) & !is.null(valoresCampoBaseSobreObservaciones)) {
    interpolacion$predictions@data[, interpolacion$campoMedia] <- interpolacion$predictions@data[, interpolacion$campoMedia] + 
      valoresCampoBase[shpMask$mask]
  }
  
  if (params$modoDiagnostico) {
    if (!is.null(valoresCampoBase) & !is.null(valoresCampoBaseSobreObservaciones)) {
      observaciones$value <- observaciones$value + valoresCampoBaseSobreObservaciones 
    }

    if (!is.null(valoresCampoBase)) {
      # Mapa de campo base
      if (gridded(coordsAInterpolar)) {
        spAux <- SpatialPixelsDataFrame(points = coordsAInterpolar, data = data.frame(value=valoresCampoBase))
        mapearGrillaGGPlot(grilla = spAux, shpBase = shpMask$shp, zcol=1, continuo = params$especEscalaDiagnostico$continuo, titulo = paste0('Ajuste Regresores - ', params$strFecha),  
                           subtitulo = params$formulaRegresionCC, nomArchResultados = paste0(params$carpetaParaModoDiagnostico, '04-AjusteRegresores.png'), 
                           dibujar = F)
      } else {
        spAux <- SpatialPointsDataFrame(coords = coordsAInterpolar, data = data.frame(value=valoresCampoBase))
        mapearPuntosGGPlot(puntos = spAux, shpBase = shpMask$shp, zcol=1, continuo = params$especEscalaDiagnostico$continuo, titulo = paste0('Ajuste Regresores - ', params$strFecha), 
                           subtitulo = params$formulaRegresionCC, nomArchResultados = paste0(params$carpetaParaModoDiagnostico, '04-AjusteRegresores.png'), 
                           dibujar = F, dibujarTexto = T)
      }
    }
    
    if (gridded(coordsAInterpolar)) {
      mapearGrillaGGPlot(
        grilla = interpolacion$predictions, shpBase = shpMask$shp, zcol=interpolacion$campoMedia, 
        continuo = params$especEscalaDiagnostico$continuo, titulo = paste0(params$nombreModelo, ' - ', params$strFecha), 
        nomArchResultados = paste0(params$carpetaParaModoDiagnostico, '06-InterpolacionMasCampoBase.png'), 
        dibujar = F, dibujarPuntosObservaciones = T, coordsObservaciones = observaciones)
    } else {
      mapearPuntosGGPlot(
        puntos = interpolacion$predictions, shpBase = shpMask$shp, zcol=interpolacion$campoMedia, 
        continuo = params$especEscalaDiagnostico$continuo, titulo = paste0(params$nombreModelo, ' - ', params$strFecha), 
        nomArchResultados = paste0(params$carpetaParaModoDiagnostico, '06-InterpolacionMasCampoBase.png'), 
        dibujar = F, dibujarTexto = T)
    }
  }

  #escala <- crearEscalaEnQuantiles(interpolacion$predictions@data[,1], nDigitos = 1, brewerPal = 'Blues', nIntervalos = 8, continuo = F)
  #mapearGrillaGGPlot(grilla = interpolacion$predictions, shpBase = shpMask$shp, escala = escala)
  if (!mapaConstante) {
    if (params$metodoRemocionDeSesgo != 'ninguno' & gridded(coordsAInterpolar) & params$interpolationMethod == 'automap') {
      interpolacion <- simpleBiasAdjustmentEx(
        observaciones = observaciones, interpolacion = interpolacion, interpolationParams = params, 
        errorRelativoParaCorregir = 0.15, shpMask=shpMask)
      
      if (params$modoDiagnostico) {
        if (gridded(coordsAInterpolar)) {
          mapearGrillaGGPlot(
            grilla = interpolacion$predictions, shpBase = shpMask$shp, zcol=interpolacion$campoMedia, 
            continuo = params$especEscalaDiagnostico$continuo, titulo = paste0(params$nombreModelo, ' - ', params$strFecha), 
            nomArchResultados = paste0(params$carpetaParaModoDiagnostico, '07-RemocionDeSesgo.png'), 
            dibujar = F, dibujarPuntosObservaciones = T, coordsObservaciones = observaciones)
        } else {
          mapearPuntosGGPlot(
            puntos = interpolacion$predictions, shpBase = shpMask$shp, zcol=interpolacion$campoMedia, 
            continuo = params$especEscalaDiagnostico$continuo, titulo = paste0(params$nombreModelo, ' - ', params$strFecha), 
            nomArchResultados = paste0(params$carpetaParaModoDiagnostico, '07-RemocionDeSesgo.png'), 
            dibujar = F, dibujarTexto = T)
        }      
      }
    }

    # el segundo chequeo, equivalente a min(observaciones) = 0 para doubles, verifica que al menos 
    # el minimo dato sea efectivamente 0, sino no tiene sentido aplicar la máscara y si es < 0 
    # además aplicarla puede hacer daño. agrego el control como robustez para evitar errores de 
    # parámetros
    if (params$umbralMascaraCeros > 0 & abs(min(observaciones$value)) < 1E-3) {
      interpolacion <- aplicarMascaraRnR(
        observaciones=observaciones, interpolacion=interpolacion, params=params, shpMask=shpMask)
      
      if (params$modoDiagnostico) {
        if (gridded(coordsAInterpolar)) {
          mapearGrillaGGPlot(
            grilla = interpolacion$predictions, shpBase = shpMask$shp, zcol=interpolacion$campoMedia, 
            continuo = params$especEscalaDiagnostico$continuo, titulo = paste0(params$nombreModelo, ' - ', params$strFecha), 
            nomArchResultados = paste0(params$carpetaParaModoDiagnostico, '09-MascaraCerosAplicada.png'), 
            dibujar = F, dibujarPuntosObservaciones = T, coordsObservaciones = observaciones)
        } else {
          mapearPuntosGGPlot(
            puntos = interpolacion$predictions, shpBase = shpMask$shp, zcol=interpolacion$campoMedia, 
            continuo = params$especEscalaDiagnostico$continuo, titulo = paste0(params$nombreModelo, ' - ', params$strFecha), 
            nomArchResultados = paste0(params$carpetaParaModoDiagnostico, '09-MascaraCerosAplicada.png'), 
            dibujar = F, dibujarTexto = T)
        }      
      }
    }
    
    # interpolacion$predictions@data[,interpolacion$campoMedia] es para seleccionar los datos independiente del método, interpolacion$campoMedia es mean para copula y var1.pred para automap
    # force predictions to be in given range
    if (!is.na(params$minVal)) interpolacion$predictions@data[,interpolacion$campoMedia][interpolacion$predictions@data[,interpolacion$campoMedia] < params$minVal] <- params$minVal
    if (!is.na(params$maxVal)) interpolacion$predictions@data[,interpolacion$campoMedia][interpolacion$predictions@data[,interpolacion$campoMedia] > params$maxVal] <- params$maxVal
    
    if (params$modoDiagnostico) {
      escalaVariableObjetivo <- darEscala(
        especificacion = params$especEscalaDiagnostico,
        valores=c(observaciones$value, interpolacion$predictions@data[,interpolacion$campoMedia]))
      
      # En modo diagnóstico guardo el mapa con las observaciones y uno para cada regresor
      mapearPuntosGGPlot(
        # Mapa de observaciones
        puntos = observaciones, shpBase = shpMask$shp, dibujarTexto = T, 
        escala = escalaVariableObjetivo, 
        titulo = paste0('Observaciones - ', params$strFecha), zcol = 'value', 
        nomArchResultados = paste0(params$carpetaParaModoDiagnostico, '01-Observaciones.png'), 
        dibujar = F)
      
      if (gridded(coordsAInterpolar)) {
        mapearGrillaGGPlot(
          grilla = interpolacion$predictions, shpBase = shpMask$shp, zcol=interpolacion$campoMedia, 
          escala = escalaVariableObjetivo, titulo = paste0(params$nombreModelo, ' - ', params$strFecha), 
          nomArchResultados = paste0(params$carpetaParaModoDiagnostico, '10-ProductoFinal.png'), 
          dibujar = F, dibujarPuntosObservaciones = T, coordsObservaciones = observaciones)
      } else {
        mapearPuntosGGPlot(
          puntos = interpolacion$predictions, shpBase = shpMask$shp, zcol=interpolacion$campoMedia, 
          escala = escalaVariableObjetivo, titulo = paste0(params$nombreModelo, ' - ', params$strFecha), 
          nomArchResultados = paste0(params$carpetaParaModoDiagnostico, '10-ProductoFinal.png'), 
          dibujar = F, dibujarTexto = T)
      }      
    }
  }
  
  # mapearPuntosGGPlot(puntos = observaciones, shpBase = shpMask$shp, zcol='value', continuo = T)
  #if (!is.null(valoresCampoBase) & !is.null(valoresCampoBaseSobreObservaciones)) {
  #  aux <- observaciones
  #  aux$value <- observaciones$value + valoresCampoBaseSobreObservaciones
  #  mapearPuntosGGPlot(puntos = aux, shpBase = shpMask$shp)
  #  aux <- interpolacion$predictions
  #  aux$var1.pred <- valoresCampoBase[shpMask$mask]
  #  mapearGrillaGGPlot(grilla = aux, shpBase=shpMask$shp)
  #}
  #mapearGrillaGGPlot(grilla = interpolacion$predictions, shpBase=shpMask$shp, continuo=T)
  
  # Additional data
  interpolacion$shpMask <- shpMask
  if (longitudesEnColumnas) {
    interpolacion$nColCoordsAInterpolar <- length(unique(coordsAInterpolar@coords[,1]))
    interpolacion$nRowCoordsAInterpolar <- length(coordsAInterpolar@coords[,2]) / interpolacion$nColCoordsAInterpolar
  } else {
    interpolacion$nColCoordsAInterpolar <- length(unique(coordsAInterpolar@coords[,2]))
    interpolacion$nRowCoordsAInterpolar <- length(coordsAInterpolar@coords[,1]) / interpolacion$nColCoordsAInterpolar
  }
  interpolacion$longitudesEnColumnas <- longitudesEnColumnas
  interpolacion$lengthCoordsAInterpolar <- length(coordsAInterpolar)
  
  return (interpolacion)
}

interpolarEx2 <- function(
    x, y, value, xs, ys, params, interpolateGrid=T, 
    proj4stringObservaciones='+proj=longlat +datum=WGS84', SRS_stringObservaciones=NULL,
    proj4stringCoordsAInterpolar='+proj=longlat +datum=WGS84', SRS_stringCoordsAInterpolar, 
    shpMask=NULL, longitudesEnColumnas=T) {
  observaciones <- crearSpatialPointsDataFrame(
    x, y, value, proj4stringObservaciones, SRS_stringObservaciones)
  coordsAInterpolar <- crearCoordsAInterpolar(
    xs, ys, interpolateGrid, proj4stringCoordsAInterpolar, SRS_stringCoordsAInterpolar)
  return (interpolarEx(observaciones, coordsAInterpolar, params, shpMask))
}


stUniversalKriging <- function(ti=1, spObservaciones, fechasObservaciones, valoresObservaciones, coordsAInterpolar,
                               valoresRegresoresSobreObservaciones=NULL, pathsRegresoresST=NULL, shpMask=NULL,
                               modelosVariograma=c('Exp', 'Sph', 'Pen'),
                               cutoff=spDists(t(bbox(spObservaciones)), longlat=!is.projected(spObservaciones))[1, 2] * 0.75,
                               tlags=0:6, nTsST=max(tlags), tlagsAR=NULL,
                               fixNugget=NA, tryFixNugget=T, nPuntosIniciales=2, 
                               modelosVariogramaST=c('Separable', 'ProductSum', 'Metric', 'SimpleSumMetric', 'SumMetric'),
                               ventanaIgualacionDistribuciones=8,
                               ndeps=1E-4, fit.method=6, deltaT=difftime(fechasObservaciones[2], fechasObservaciones[1]),
                               pathsSalida=NULL, verbose=F, params) {
  if (!is.null(pathsRegresoresST)) {
    valoresRegresoresSobreCoordsAInterpolar_ti <- extraerValoresRegresoresSobreSP(coordsAInterpolar, pathsRegresores = pathsRegresoresST[ti, , drop=F])
    valoresRegresoresSobreCoordsAInterpolar_ti <- lapply(valoresRegresoresSobreCoordsAInterpolar_ti, FUN = function(x) { as.vector(t(x)) } )
  } else { valoresRegresoresSobreCoordsAInterpolar_ti <- NULL }
  
  return(stUniversalKrigingEx(ti = ti, spObservaciones = spObservaciones, fechasObservaciones = fechasObservaciones, valoresObservaciones = valoresObservaciones,
                              coordsAInterpolar = coordsAInterpolar, valoresRegresoresSobreObservaciones = valoresRegresoresSobreObservaciones, 
                              valoresRegresoresSobreCoordsAInterpolar_ti = valoresRegresoresSobreCoordsAInterpolar_ti, shpMask = shpMask, 
                              modelosVariograma = modelosVariograma, cutoff = cutoff, tlags = tlags, nTsST = nTsST, tlagsAR = tlagsAR, 
                              fixNugget = fixNugget, tryFixNugget = tryFixNugget, nPuntosIniciales = nPuntosIniciales, modelosVariogramaST=modelosVariogramaST,
                              ventanaIgualacionDistribuciones=ventanaIgualacionDistribuciones, ndeps = ndeps, fit.method = fit.method, 
                              pathsSalida = pathsSalida, verbose = verbose, params=params))
}

fitSTVariogramModelI <- function(i, modelosVariogramaST, 
                                 empVgm, spatialEmpVgm, temporalEmpVgm, estSpaceVgm, estTimeVgm, estJointVgm, estTimeAnis, 
                                 fml, tlags, modelosVariograma, fixNugget, tryFixNugget, estJointSill, cutoff, nPuntosIniciales, fit.method, 
                                 ndeps, errorLevel) {
  # i <- 5
  print(paste0(i, ' - ', modelosVariogramaST[i]))
  
  if (modelosVariogramaST[i] == 'Separable') {
    # 1 - Separable
    auxS <- getModelVariogram(experimental_variogram = spatialEmpVgm, formula = fml, model=modelosVariograma, fix.values=c(fixNugget, NA, estJointSill), nPuntosIniciales = nPuntosIniciales, tryFixNugget=tryFixNugget, start_vals = c(0, cutoff, NA), quitarNuggetSiEsCero = F)$var_model
    auxT <- getModelVariogram(experimental_variogram = temporalEmpVgm, formula = fml, model=modelosVariograma, fix.values=c(fixNugget, NA, estJointSill), nPuntosIniciales = nPuntosIniciales, tryFixNugget=tryFixNugget, start_vals = c(0, max(tlags), NA), quitarNuggetSiEsCero = F)$var_model
    #auxS$range[2] <- auxS$range[2] * getExperimentalSill(auxS, cutoff)
    auxS$psill <- auxS$psill / sum(auxS$psill)
    #auxT$range[2] <- auxT$range[2] * getExperimentalSill(auxT, cutoff)
    auxT$psill <- auxT$psill / sum(auxT$psill)
    if ('Pow' %in% auxS$model) auxS$range[2] <- min(auxS$range[2], 2)
    if ('Pow' %in% auxT$model) auxT$range[2] <- min(auxT$range[2], 2)
    
    separableModel <- vgmST('separable', space = auxS, time = auxT, sill = estJointSill)

    minPsill <- 0.01
    maxPsill <- max(empVgm$gamma, na.rm = T)
    maxSNugget <- min(empVgm$gamma[empVgm$timelag == 0], na.rm = T) / maxPsill
    maxTNugget <- min(empVgm$gamma[empVgm$spacelag == 0], na.rm = T) / maxPsill
    if (!is.na(fixNugget) && fixNugget == 0) {
      maxSNugget <- min(maxTNugget, maxSNugget) * 1E-4
      if (maxSNugget <= 0) maxSNugget <- min(empVgm$gamma[empVgm$gamma > 0], na.rm = T) / maxPsill * 1E-4
      maxTNugget <- maxSNugget
    } 
    
    minSRange <- min(empVgm$spacelag[empVgm$spacelag > 0], na.rm = T) * 0.01
    minTRange <- min(empVgm$timelag[empVgm$timelag > 0], na.rm = T) * 0.01
    if ('Pow' %in% auxS$model) { maxSRange <- 2
    } else { maxSRange <- max(spatialEmpVgm$dist)}
    if ('Pow' %in% auxT$model) { maxTRange <- 2
    } else { maxTRange <- max(temporalEmpVgm$dist)}
    
    parscale <- extractPar(separableModel)
    parscale[parscale <= 0.01] <- min(parscale[parscale > 0.1])
    
    separableModel <- fit.StVariogramEx(empVgm, separableModel, fit.method = fit.method,
                                        stAni = estTimeAnis, method = 'L-BFGS-B',
                                        control=list(parscale = parscale, maxit=10000, ndeps=rep(ndeps, 5)),
                                        lower = c(minSRange, 0, minTRange, 0, minPsill),
                                        upper = c(maxSRange, maxSNugget, maxTRange, maxTNugget, maxPsill),
                                        errorLevel = errorLevel)
    # attr(separableModel, 'optim.output')
    # print(plot(empVgm, list(separableModel), all=T, wireframe=T, zlim=c(0, max(empVgm$gamma,na.rm=T) * 1.1), zlab=NULL, xlab=list('distance (km)', rot=30), ylab=list('time lag (days)', rot=-35), scales=list(arrows=F, z = list(distance = 5))))
    rm(auxS, auxT)
    return(separableModel)
  }
  
  if (modelosVariogramaST[i] == 'ProductSum') {
    # 2 - Product Sum
    estSpaceSill <- getExperimentalSill(variogramModel = estSpaceVgm, maxDist = cutoff)
    estTimeSill <- getExperimentalSill(variogramModel = estTimeVgm, maxDist = max(tlags))
    estK <- (estSpaceSill + estTimeSill - estJointSill) / (estSpaceSill * estTimeSill)
    if (estK <= 0) estK <- 1
    
    prodSumModel <- vgmST('productSum', space = estSpaceVgm, time = estTimeVgm, k = estK)
    
    maxSNugget <- min(empVgm$gamma[empVgm$timelag == 0], na.rm = T)
    maxTNugget <- min(empVgm$gamma[empVgm$spacelag == 0], na.rm = T)
    if (!is.na(fixNugget) && fixNugget == 0) {
      maxSNugget <- min(maxTNugget, maxSNugget) * 1E-4
      if (maxSNugget <= 0) maxSNugget <- min(empVgm$gamma[empVgm$gamma > 0], na.rm = T) * 1E-4
      maxTNugget <- maxSNugget
    } 
    
    minSRange <- min(empVgm$spacelag[empVgm$spacelag > 0], na.rm = T) * 0.01
    minTRange <- min(empVgm$timelag[empVgm$timelag > 0], na.rm = T) * 0.01
    if ('Pow' %in% estSpaceVgm$model) { maxSRange <- 2
    } else { maxSRange <- max(spatialEmpVgm$dist)}
    if ('Pow' %in% estTimeVgm$model) { maxTRange <- 2
    } else { maxTRange <- max(temporalEmpVgm$dist)}
    
    r <- range(empVgm$gamma[empVgm$timelag == 0], na.rm = T)
    minSPsill <- r[1] * 0.01
    maxSPsill <- r[2]
    
    r <- range(empVgm$gamma[empVgm$spacelag == 0], na.rm = T)
    minTPsill <- r[1] * 0.01
    maxTPsill <- r[2]
    
    parscale <- extractPar(prodSumModel)
    parscale[parscale <= 0.01] <- min(parscale[parscale > 0.1])
    
    prodSumModel <- fit.StVariogramEx(empVgm, prodSumModel, fit.method = fit.method,
                                      stAni = estTimeAnis, method = 'L-BFGS-B',
                                      control=list(parscale=parscale, maxit=1000, ndeps=rep(ndeps, 7)),
                                      lower = c(minSPsill, minSRange, 0, minTPsill, minTRange, 0, 1E-3), 
                                      upper = c(maxSPsill, maxSRange, maxSNugget, maxTPsill, maxTRange, maxTNugget, estK * 1000),
                                      errorLevel = errorLevel)
    return(prodSumModel)
  }
  
  if (modelosVariogramaST[i] == 'Metric') {
    # 3 - Metric
    metricModel <- vgmST('metric', joint = estJointVgm, stAni = estTimeAnis)
    
    if (!is.na(fixNugget) && fixNugget == 0) { maxNugget <- min(empVgm$gamma[empVgm$gamma > 0], na.rm=T) * 1E-4
    } else { maxNugget <- min(empVgm$gamma[empVgm$gamma > 0], na.rm=T) }
    
    r <- range(empVgm$gamma, na.rm = T)
    minPsill <- r[1] * 0.01
    maxPsill <- r[2]
    
    if ('Pow' %in% estJointVgm$model) { maxRange <- 2
    } else { maxRange <- max(empVgm$dist, na.rm = T)}

    parscale <- extractPar(metricModel)
    parscale[parscale <= 0.01] <- min(parscale[parscale > 0.1])
    
    metricModel <- fit.StVariogramEx(object = empVgm, model = metricModel, fit.method = fit.method,
                                     method = 'L-BFGS-B',
                                     control=list(parscale = parscale, maxit=1000, ndeps=rep(ndeps, 4)),
                                     lower = c(minPsill, 1E-3, 0, 1E-3),
                                     upper = c(maxPsill, maxSRange, maxNugget, estTimeAnis * 1E6), errorLevel = errorLevel)
    #metricModel2 <- fit.StVariogramEx2(object = empVgm, model = metricModel, fit.method = fit.method,
    #                                  method = 'L-BFGS-B', stAni = estTimeAnis, traceGr=T,
    #                                  control=list(parscale = parscale, maxit=1000, ndeps=rep(ndeps, 4)),
    #                                  lower = c(0, 0, 0, 0), errorLevel = errorLevel)  
    #attr(metricModel1, 'optim')
    #attr(metricModel2, 'optim')
    #print(plot(empVgm, list(metricModel1), all=T, wireframe=T, zlim=c(0, max(empVgm$gamma,na.rm=T) * 1.1), zlab=NULL, xlab=list('distance (km)', rot=30), ylab=list('time lag (days)', rot=-35), scales=list(arrows=F, z = list(distance = 5))))
    #print(plot(empVgm, list(metricModel2), all=T, wireframe=T, zlim=c(0, max(empVgm$gamma,na.rm=T) * 1.1), zlab=NULL, xlab=list('distance (km)', rot=30), ylab=list('time lag (days)', rot=-35), scales=list(arrows=F, z = list(distance = 5))))
    return(metricModel)
  }
  
  if (modelosVariogramaST[i] == 'SimpleSumMetric') {
    # 4 - simpleSumMetric
    simpleSumMetricModel <- vgmST('simpleSumMetric', space = estSpaceVgm, time = estTimeVgm, joint = estJointVgm, nugget = estJointVgm$psill[1],
                                  stAni = estTimeAnis)
    
    maxNugget <- min(empVgm$gamma, na.rm = T)
    if (!is.na(fixNugget) && fixNugget == 0) {
      maxNugget <- maxNugget * 1E-4
      if (maxNugget <= 0) maxNugget <- min(empVgm$gamma[empVgm$gamma > 0], na.rm = T) * 1E-4
    } 
    
    minJRange <- min(empVgm$dist[empVgm$dist > 0], na.rm = T) * 0.01
    minSRange <- min(empVgm$spacelag[empVgm$spacelag > 0], na.rm = T) * 0.01
    minTRange <- min(empVgm$timelag[empVgm$timelag > 0], na.rm = T) * 0.01
    if ('Pow' %in% estSpaceVgm$model) { maxSRange <- 2
    } else { maxSRange <- max(spatialEmpVgm$dist)}
    if ('Pow' %in% estTimeVgm$model) { maxTRange <- 2
    } else { maxTRange <- max(temporalEmpVgm$dist)}
    if ('Pow' %in% estJointVgm$model) { maxJRange <- 2
    } else { maxJRange <- max(empVgm$dist, na.rm = T)}
    
    r <- range(empVgm$gamma[empVgm$gamma > 0], na.rm=T)
    minJPsill <- r[1] * 0.01
    maxJPsill <- r[2]
    
    r <- range(empVgm$gamma[empVgm$timelag == 0], na.rm = T)
    minSPsill <- r[1] * 0.01
    maxSPsill <- r[2]
    
    r <- range(empVgm$gamma[empVgm$spacelag == 0], na.rm = T)
    minTPsill <- r[1] * 0.01
    maxTPsill <- r[2]
    
    parscale <- extractPar(simpleSumMetricModel)
    parscale[parscale <= 0.01] <- min(parscale[parscale > 0.1])
    
    simpleSumMetricModel <- fit.StVariogramEx(object = empVgm, model = simpleSumMetricModel,
                                              fit.method = fit.method, method = 'L-BFGS-B',
                                              control=list(parscale=parscale, maxit=1000, ndeps=rep(ndeps, 8)),
                                              lower = c(sill.s = minSPsill, range.s = minSRange, 
                                                        sill.t = minTPsill, range.t = minTRange,
                                                        sill.st= minJPsill, range.st= minJRange, nugget = 0, anis = 1E-3), 
                                              upper = c(sill.s = maxSPsill, range.s = maxSRange, 
                                                        sill.t = maxTPsill, range.t = maxTRange,
                                                        sill.st= maxJPsill, range.st= maxJRange, nugget = maxNugget, 
                                                        anis = estTimeAnis * 1E6), 
                                              errorLevel = errorLevel)
    return(simpleSumMetricModel)
  }
  
  if (modelosVariogramaST[i] == 'SumMetric') {
    # 5 - sumMetric
    sumMetricModel <- vgmST('sumMetric', space = estSpaceVgm, time = estTimeVgm, joint = estJointVgm, stAni = estTimeAnis)
    
    maxSNugget <- min(empVgm$gamma[empVgm$timelag == 0], na.rm = T)
    maxTNugget <- min(empVgm$gamma[empVgm$spacelag == 0], na.rm = T)
    maxJNugget <- min(empVgm$gamma, na.rm = T)
    if (!is.na(fixNugget) && fixNugget == 0) {
      maxSNugget <- min(maxTNugget, maxSNugget, maxJNugget) * 1E-4
      if (maxSNugget <= 0) maxSNugget <- min(empVgm$gamma[empVgm$gamma > 0], na.rm = T) * 1E-4
      maxTNugget <- maxSNugget
      maxJNugget <- maxSNugget
    }
    
    minJRange <- min(empVgm$dist[empVgm$dist > 0], na.rm = T) * 0.01
    minSRange <- min(empVgm$spacelag[empVgm$spacelag > 0], na.rm = T) * 0.01
    minTRange <- min(empVgm$timelag[empVgm$timelag > 0], na.rm = T) * 0.01
    if ('Pow' %in% estSpaceVgm$model) { maxSRange <- 2
    } else { maxSRange <- max(spatialEmpVgm$dist)}
    if ('Pow' %in% estTimeVgm$model) { maxTRange <- 2
    } else { maxTRange <- max(temporalEmpVgm$dist)}
    if ('Pow' %in% estJointVgm$model) { maxJRange <- 2
    } else { maxJRange <- max(empVgm$dist, na.rm = T)}    
    
    r <- range(empVgm$gamma[empVgm$gamma > 0], na.rm=T)
    minJPsill <- r[1] * 0.01
    maxJPsill <- r[2]
    
    r <- range(empVgm$gamma[empVgm$timelag == 0], na.rm = T)
    minSPsill <- r[1] * 0.01
    maxSPsill <- r[2]
    
    r <- range(empVgm$gamma[empVgm$spacelag == 0], na.rm = T)
    minTPsill <- r[1] * 0.01
    maxTPsill <- r[2]
    
    parscale <- extractPar(sumMetricModel)
    parscale[parscale <= 0.01] <- min(parscale[parscale > 0.1])
    
    sumMetricModel <- fit.StVariogram(empVgm, sumMetricModel, fit.method = fit.method,
                                      method = 'L-BFGS-B',
                                      control=list(parscale = parscale, maxit=1000, ndeps=rep(ndeps, 10)),
                                      lower = c(sill.s = minSPsill, range.s = minSRange, nugget.s = 0,
                                                sill.t = minTPsill, range.t = minTRange, nugget.t = 0,
                                                sill.st = minJPsill, range.st = minJRange, nugget.st = 0,
                                                anis = 1E-3),
                                      upper = c(sill.s = maxSPsill, range.s = maxSRange, nugget.s = maxSNugget,
                                                sill.t = maxTPsill, range.t = maxTRange, nugget.t = maxTNugget,
                                                sill.st = maxJPsill, range.st = maxJRange, nugget.st = maxJNugget,
                                                anis = estTimeAnis * 1E6), errorLevel = errorLevel)
    return(sumMetricModel)
  }
}

stUniversalKrigingEx <- function(ti=1, spObservaciones, fechasObservaciones, valoresObservaciones, coordsAInterpolar,
                                 valoresRegresoresSobreObservaciones=NULL, valoresRegresoresSobreCoordsAInterpolar_ti=NULL, 
                                 shpMask=NULL, modelosVariograma=c('Exp', 'Sph', 'Pen'), cutoff=NA, tlags=0:5, nTsST=max(tlags), 
                                 tlagsAR=NULL, fixNugget=NA, tryFixNugget=T, nPuntosIniciales=2, 
                                 modelosVariogramaST=c('Separable', 'ProductSum', 'Metric', 'SimpleSumMetric', 'SumMetric'),
                                 ventanaIgualacionDistribuciones=8, ndeps=1E-4, fit.method=6, 
                                 deltaT=difftime(fechasObservaciones[2], fechasObservaciones[1]), verbose=F, params) {
  # spObservaciones = coordsObservaciones
  # ndeps=1E-4, fit.method=6, deltaT=difftime(fechasObservaciones[2], fechasObservaciones[1])
  modelosVariograma <- modelosVariograma[!(modelosVariograma %in% c('Pow', 'Cir'))]
  if (params$incorporarCoordenadas | params$incorporarTiempo | params$incorporarDistanciaAlAgua | params$incorporarAltitud) {
    regs <- incorporarRegresoresEstaticos(ti = ti, coordsObservaciones = spObservaciones, fechasObservaciones=fechasObservaciones, 
                                          valoresObservaciones=valoresObservaciones, coordsAInterpolar=coordsAInterpolar, params=params, 
                                          valoresRegresoresSobreObservaciones=valoresRegresoresSobreObservaciones, 
                                          valoresRegresoresSobreCoordsAInterpolar_ti=valoresRegresoresSobreCoordsAInterpolar_ti,
                                          incorporarCoordenadas=params$incorporarCoordenadas, formulaCoordenadas=params$formulaCoordenadas, 
                                          incorporarTiempo=params$incorporarTiempo, formulaTiempo=params$formulaTiempo,
                                          incorporarDistanciaAlAgua=params$incorporarDistanciaAlAgua, formulaDistanciaAlAgua=params$formulaDistanciaAlAgua,
                                          incorporarAltitud=params$incorporarAltitud, formulaAltitud=params$formulaAltitud)
    
    valoresRegresoresSobreObservaciones <- regs$valoresRegresoresSobreObservaciones
    valoresRegresoresSobreCoordsAInterpolar_ti <- regs$valoresRegresoresSobreCoordsAInterpolar_ti
    rm(regs)
  }
  
  if (verbose) print(paste0('ti ', ti, ': ', fechasObservaciones[ti]))

  ventana <- getVentana(ti=ti, nT=nrow(valoresObservaciones), tamanioSemiVentana = ventanaIgualacionDistribuciones, tlagsAR=tlagsAR, tlags=tlags)
  tsVentana <- ventana$tsVentana
  iTiEnTsVentana <- ventana$iTiEnTsVentana
  
  fechasVentana <- fechasObservaciones[tsVentana, drop=F]
  valoresVentana <- valoresObservaciones[tsVentana, , drop=F]

  if (is.na(cutoff)) cutoff=quantile(spDists(spObservaciones, longlat=!is.projected(spObservaciones)), probs = 0.75)
  if (nTsST <= 0) nTsST <- length(tsVentana)
  iFechasST <- seq.int(from=max((iTiEnTsVentana - nTsST) * length(spObservaciones) + 1, 1), to = (iTiEnTsVentana) * length(spObservaciones))
  iFechasST <- iFechasST[iFechasST > 0]

  if (sum(!is.na(valoresVentana[iTiEnTsVentana, ])) > 5 && var(valoresVentana[iTiEnTsVentana, ], na.rm=T) > 1E-6) {
    #regresoresConVarEnTi <- lapply(valoresRegresoresSobreObservaciones, FUN = function(x, ti) { var(x[ti,], na.rm = T) }, ti=ti)
    #regresoresConVarEnTi <- !is.na(regresoresConVarEnTi) & regresoresConVarEnTi > 1E-6
    #valoresRegresoresSobreObservaciones <- valoresRegresoresSobreObservaciones[regresoresConVarEnTi]
    #valoresRegresoresSobreCoordsAInterpolar_ti <- valoresRegresoresSobreCoordsAInterpolar_ti[, regresoresConVarEnTi, drop=F]
    
    if (!is.null(valoresRegresoresSobreObservaciones) & length(valoresRegresoresSobreObservaciones) > 0) {
      regs <- c(lapply(valoresRegresoresSobreObservaciones, function(x, tsVentana) { return (x[tsVentana, , drop=F]) }, tsVentana=tsVentana))
      
      df <- data.frame(value=as.vector(t(valoresVentana)), lapply(regs, FUN = function(x) { as.vector(t(x)) } ))
      regresoresConVarianza <- apply(na.omit(df[iFechasST,]), MARGIN = 2, FUN = var) > 1E-6
      df <- df[, regresoresConVarianza, drop=F]
      regresoresConVarianza <- regresoresConVarianza[2:length(regresoresConVarianza)]
      valoresRegresoresSobreObservaciones <- valoresRegresoresSobreObservaciones[regresoresConVarianza]
      valoresRegresoresSobreCoordsAInterpolar_ti <- valoresRegresoresSobreCoordsAInterpolar_ti[, regresoresConVarianza, drop=F]

      rm(regresoresConVarianza)
      #vars <- lapply(regs, function(x, iTiEnTsVentana) { return(var(x[iTiEnTsVentana, ], na.rm = T)) }, iTiEnTsVentana=iTiEnTsVentana)
      #regs <- regs[!is.na(vars) & vars > 1E-6]
      
      hayRegresores <- ncol(df) > 1
      if (hayRegresores) {
        #df <- data.frame(value=as.vector(t(valoresVentana)), lapply(regs, FUN = function(x) { as.vector(t(x)) } ))
        fml <- as.formula(paste(colnames(df)[1], '~', paste('+', colnames(df)[2:ncol(df)], collapse=''), '+1'))  
      } else {
        #df <- data.frame(value=as.vector(t(valoresVentana)))
        fml <- as.formula('value~1')
      }
    } else {
      hayRegresores <- F
      regs <- list()
      df <- data.frame(value=as.vector(t(valoresVentana)))
      fml <- as.formula('value~1')
    }
    
    # gls(model = fml, data = df, na.action = na.omit, correlation = corAR1(form = ~ 1))
    sptData <- STFDF(sp = geometry(spObservaciones), time = fechasVentana, endTime = fechasVentana + deltaT, data = df)
    # acf(na.omit(as(sptData, 'xts')))

    source(paste0(script.dir.interpolarEx, 'afvmod.r'), encoding = 'WINDOWS-1252')
    source(paste0(script.dir.interpolarEx, 'getBoundariesPVariogramaEmpirico.r'), encoding = 'WINDOWS-1252')    
    # Space-Time empirical variogram
    boundaries <- getBoundariesPVariogramaEmpiricoV4_MultiTime(formula=fml, input_data = spObservaciones, input_data_t = c(list(value=valoresVentana), regs), cutoff=cutoff)
    empVgm <- variogramST(formula = fml, locations = sptData, tlags = tlags, boundaries = boundaries, progress=verbose)
    rm(regs)
    
    # Pure spatial variogram
    spatialEmpVgm <- na.omit(empVgm[empVgm$timelag==0, ])
    spatialEmpVgm$dir.hor <- 0
    spatialEmpVgm$dir.ver <- 0
    spatialEmpVgm <- spatialEmpVgm[, c('np', 'dist', 'gamma', 'dir.hor', 'dir.ver', 'id')]
    class(spatialEmpVgm) <- c('gstatVariogram', 'data.frame')
    spatialModelVgm <- getModelVariogram(experimental_variogram = spatialEmpVgm, formula = fml, model=modelosVariograma, fix.values=c(fixNugget, NA, NA), nPuntosIniciales = nPuntosIniciales, tryFixNugget=tryFixNugget, quitarNuggetSiEsCero = F)
    if (verbose) print(annotatedplot(spatialModelVgm), main='Variograma Espacial')
    estSpaceVgm <- spatialModelVgm$var_model  
    
    rm(spatialModelVgm)
    
    # Pure temporal variogram
    temporalEmpVgm <- na.omit(empVgm[empVgm$spacelag==0, ])
    temporalEmpVgm$dist <- temporalEmpVgm$timelag
    temporalEmpVgm$id <- temporalEmpVgm$id[1]
    temporalEmpVgm$dir.hor <- 0
    temporalEmpVgm$dir.ver <- 0
    temporalEmpVgm <- temporalEmpVgm[, c('np', 'dist', 'gamma', 'dir.hor', 'dir.ver', 'id')]
    class(temporalEmpVgm) <- class(spatialEmpVgm)
    temporalModelVgm <- getModelVariogram(experimental_variogram = temporalEmpVgm, formula = fml, model=modelosVariograma, fix.values=c(fixNugget, NA, NA), nPuntosIniciales = nPuntosIniciales, tryFixNugget=tryFixNugget, quitarNuggetSiEsCero = F)
    if (verbose) print(annotatedplot(temporalModelVgm), main='Variograma Temporal')
    estTimeVgm <- temporalModelVgm$var_model
    
    rm(temporalModelVgm)
    
    estTimeAnis <- median(estimateSTAnisotropyEx(empVgm=empVgm, methods=c('metric', 'vgm', 'range', 'linear'), spatialVgm = estSpaceVgm, temporalVgm = estTimeVgm))
    
    # Joint temporal variogram
    # nSLags <- length(unique(empVgm$spacelag))
    jointEmpVgm <- na.omit(empVgm)
    jointEmpVgm$dist <- sqrt(jointEmpVgm$dist^2 + (jointEmpVgm$timelag * estTimeAnis)^2)
    jointEmpVgm <- jointEmpVgm[order(jointEmpVgm$dist),]
    jointEmpVgm$id <- jointEmpVgm$id[1]
    jointEmpVgm$dir.hor <- 0
    jointEmpVgm$dir.ver <- 0
    jointEmpVgm <- jointEmpVgm[, c('np', 'dist', 'gamma', 'dir.hor', 'dir.ver', 'id')]
    class(jointEmpVgm) <- class(spatialEmpVgm)
    jointModelVgm <- getModelVariogram(experimental_variogram = jointEmpVgm, formula = fml, model=modelosVariograma, fix.values=c(fixNugget, NA, NA), nPuntosIniciales = nPuntosIniciales, tryFixNugget=tryFixNugget, quitarNuggetSiEsCero = F)
    if (verbose) print(annotatedplot(jointModelVgm), main='Variograma Conjunto')
    estJointVgm <- jointModelVgm$var_model
    estJointSill <- getExperimentalSill(estJointVgm, cutoff)
    
    rm(jointEmpVgm, jointModelVgm)
    
    # Joint marginal variogram  
    #jointEmpVgm2 <- na.omit(empVgm[seq(from=1, length.out = nSLags, by=nSLags) + tlags[1:nSLags], ])
    #jointEmpVgm2$dist <- sqrt(jointEmpVgm2$dist^2 + (jointEmpVgm2$timelag * estTimeAnis)^2)
    #jointEmpVgm2$id <- jointEmpVgm2$id[1]
    #jointEmpVgm2$dir.hor <- 0
    #jointEmpVgm2$dir.ver <- 0
    #jointEmpVgm2 <- jointEmpVgm2[, c('np', 'dist', 'gamma', 'dir.hor', 'dir.ver', 'id')]
    #class(jointEmpVgm2) <- class(spatialEmpVgm)
    
    #require(matrixStats)
    #jointEmpVgm2$gamma <- pmax(0, jointEmpVgm2$gamma - rowMaxs(cbind(empVgm$gamma[2:(nrow(jointEmpVgm2)+1)], empVgm$gamma[seq.int(from=1, length.out = nrow(jointEmpVgm2)+1, by = nSLags)[-1]])))
    #jointModelVgm2 <- getModelVariogram(experimental_variogram = jointEmpVgm2, formula = fml, model=modelosVariograma, fix.values=c(fixNugget, NA, NA), nPuntosIniciales = nPuntosIniciales, tryFixNugget=tryFixNugget, quitarNuggetSiEsCero = F)
    #if (verbose) print(annotatedplot(jointModelVgm2), main='Variograma Conjunto 2')
    #estJointVgm2 <- jointModelVgm2$var_model
    
    #rm(jointEmpVgm2, jointModelVgm2)
    errorLevel <- 0
    
    if (params$nCoresAUsar <= 0) {nCoresAUsar <- min(length(modelosVariogramaST), detectCores(T, T)) 
    } else { nCoresAUsar <- min(length(modelosVariogramaST), params$nCoresAUsar) }
    
    if (nCoresAUsar > 1) {
      cl <- makeCluster(getOption('cl.cores', nCoresAUsar))
      clusterExport(cl, varlist = c('script.dir.interpolarEx'))
      clusterEvalQ(cl = cl, {
        require('gstat')
        source(paste0(script.dir.interpolarEx, 'afvmod.r'), encoding = 'WINDOWS-1252')
      })
      if (exists(x = 'setMKLthreads')) { clusterEvalQ(cl = cl, expr = setMKLthreads(1)) }
      modelos  <- parSapplyLB(cl = cl, X = 1:length(modelosVariogramaST), FUN = fitSTVariogramModelI, modelosVariogramaST=modelosVariogramaST, 
                              empVgm=empVgm, spatialEmpVgm=spatialEmpVgm, temporalEmpVgm=temporalEmpVgm, estSpaceVgm=estSpaceVgm, 
                              estTimeVgm=estTimeVgm, estJointVgm=estJointVgm, estTimeAnis=estTimeAnis, fml=fml, tlags=tlags, modelosVariograma=modelosVariograma, 
                              fixNugget=fixNugget, tryFixNugget=tryFixNugget, estJointSill=estJointSill, cutoff=cutoff, nPuntosIniciales=nPuntosIniciales, fit.method=fit.method, ndeps=ndeps, 
                              errorLevel=errorLevel)
      stopCluster(cl)
    } else {
      modelos <- sapply(X = 1:length(modelosVariogramaST), FUN = fitSTVariogramModelI, modelosVariogramaST=modelosVariogramaST, 
                        empVgm=empVgm, spatialEmpVgm=spatialEmpVgm, temporalEmpVgm=temporalEmpVgm, estSpaceVgm=estSpaceVgm, 
                        estTimeVgm=estTimeVgm, estJointVgm=estJointVgm, estTimeAnis=estTimeAnis, fml=fml, tlags=tlags, modelosVariograma=modelosVariograma, 
                        fixNugget=fixNugget, tryFixNugget=tryFixNugget, estJointSill=estJointSill, cutoff=cutoff, nPuntosIniciales=nPuntosIniciales, fit.method=fit.method, ndeps=ndeps, 
                        errorLevel=errorLevel)
    }
    
    iModelo <- which.min(lapply(X = modelos, FUN = function(x) { attr(x, 'optim')$value }))
    
    if (verbose && interactive()) {
      names(modelos) <- sapply(X=modelos, FUN = function(x) { x$stModel })
      print(plot(empVgm, modelos,
                 all=T, wireframe=T, zlim=c(0, max(empVgm$gamma, na.rm=T) * 1.1), zlab=NULL, 
                 xlab=list('distance (km)', rot=30),
                 ylab=list('time lag (days)', rot=-35),
                 scales=list(arrows=F, z = list(distance = 5))))
      
      aux <- t(sapply(X = modelos, FUN = function(x) { c(x$stModel, attr(x, 'optim')$value) }))
      print(aux[order(aux[,2]), ])
      nomModelo <- sapply(X=modelos, FUN = function(x) { x$stModel })[iModelo]
      MSEModelo <- round(sapply(X=modelos, FUN = function(x) { attr(x, 'MSE') })[iModelo], 2)
      OptimModelo <- sapply(X=modelos, FUN = function(x) { attr(x, 'optim')$value })[iModelo]
      
      #nVecesModelos[nomModelo] <- nVecesModelos[nomModelo] + 1
      #print(paste0('MejorModelo= ', nomModelo, '. MSE= ', MSEModelo, '. Optim= ', OptimModelo, '. nVecesModelos=', paste(nVecesModelos, collapse = ' ')))
    }
    
    mejorModelo <- modelos[[iModelo]]
    rm(modelos)
    
    if (!is.null(mejorModelo$stAni)) estTimeAnis <- mejorModelo$stAni
    if (is.null(shpMask)) { msk <- T 
    } else { msk <- shpMask$mask }
    
    if (hayRegresores) {
      dfPredLocs <- data.frame(valoresRegresoresSobreCoordsAInterpolar_ti)[msk, , drop=F]
      #dfPredLocs[!msk,] <- NA
      puedoInterpolar <- !all(apply(X = dfPredLocs, MARGIN = 1, FUN = function(x) { any(is.na(x))} ))
      if (puedoInterpolar) {
        predLocs <- STFDF(sp = coordsAInterpolar[msk,], time = fechasObservaciones[ti], endTime = fechasObservaciones[ti] + deltaT, data = dfPredLocs)    
        predLocs <- as.STSDF.STFDFEx(from=predLocs)
      }
      rm(dfPredLocs)
    } else {
      predLocs <- STF(sp = coordsAInterpolar[msk,], time = fechasObservaciones[ti], endTime = fechasObservaciones[ti] + deltaT)
      puedoInterpolar <- TRUE
    }

    if (puedoInterpolar) {
      aux <- as.STSDF.STFDFEx(STFDF(sp = geometry(spObservaciones), time = fechasVentana[max(iTiEnTsVentana - nTsST + 1,1):iTiEnTsVentana], 
                                    endTime = fechasVentana[max(iTiEnTsVentana - nTsST + 1,1):iTiEnTsVentana] + deltaT, data = df[iFechasST, ,drop=F]))
      aux <- spTransform(aux, predLocs@proj4string)
      
      puedoInterpolar <- nrow(aux@data) > 5
    }

    rm(empVgm, spatialEmpVgm, temporalEmpVgm, estSpaceVgm, estTimeVgm, estJointVgm, estTimeAnis)
    # rm(estJointVgm2)
  } else { 
    fml <- as.formula('value~1')
    puedoInterpolar <- FALSE 
  }

  if (puedoInterpolar) {
    interpolacionST <- krigeST(formula = fml, data = aux, newdata = predLocs, modelList = mejorModelo, stAni = estTimeAnis)
    if (inherits(interpolacionST, 'STSDF')) { interpolacionST <- spacetime:::as.STFDF.STSDF(interpolacionST) }
    
    if (gridded(coordsAInterpolar)) { interpolacionST <- SpatialPixelsDataFrame(interpolacionST@sp, data=interpolacionST@data)
    } else { interpolacionST <- SpatialPointsDataFrame(interpolacionST@sp, data=interpolacionST@data) }
    
    # mapearGrillaGGPlot(grilla = interpolacionST, shpBase = shpMask$shp)
    rm(aux, predLocs)
  } else {
    meanVal <- mean(valoresVentana[iTiEnTsVentana, ], na.rm=T)
    if (gridded(coordsAInterpolar)) { interpolacionST <- SpatialPixelsDataFrame(coordsAInterpolar[shpMask$mask,], data=data.frame(value=rep(meanVal, sum(shpMask$mask))))
    } else { interpolacionST <- SpatialPointsDataFrame(coordsAInterpolar, data=data.frame(value=rep(meanVal, length(coordsAInterpolar)))) }
    class(interpolacionST@data[,1]) <- 'numeric'
  }

  spObservaciones$value <- valoresVentana[iTiEnTsVentana,]
  interpolacionST <- imitarObjetoIntamap(observaciones = spObservaciones, formulaString = fml, predictions = interpolacionST)
  
  if (puedoInterpolar) {
    if (params$metodoRemocionDeSesgo != 'ninguno' & gridded(coordsAInterpolar)) {
      obs <- SpatialPointsDataFrame(coords = spObservaciones, data = data.frame(value=valoresVentana[iTiEnTsVentana, ]), proj4string = proj4string(spObservaciones))
      interpolacionST <- simpleBiasAdjustmentEx(
        observaciones=obs, interpolacion=interpolacionST, interpolationParams=params,
        errorRelativoParaCorregir=0.15, inverseDistancePower=3, shpMask=shpMask)
      rm(obs)
    }
    
    # TO-DO: hay que ver como incluir el campo base en la mascara
    if (params$umbralMascaraCeros > 0)
      interpolacionST <- aplicarMascaraRnR(
        observaciones=spObservaciones, interpolacion=interpolacionST, params=params, 
        shpMask=shpMask)
  }
  
  # Additional data
  interpolacionST$shpMask <- shpMask
  interpolacionST$longitudesEnColumnas <- T
  if (interpolacionST$longitudesEnColumnas) {
    interpolacionST$nColCoordsAInterpolar <- length(unique(spObservaciones@coords[,1]))
    interpolacionST$nRowCoordsAInterpolar <- length(spObservaciones@coords[,2]) / interpolacionST$nColCoordsAInterpolar
  } else {
    interpolacionST$nColCoordsAInterpolar <- length(unique(spObservaciones@coords[,2]))
    interpolacionST$nRowCoordsAInterpolar <- length(spObservaciones@coords[,1]) / interpolacionST$nColCoordsAInterpolar
  }
  interpolacionST$lengthCoordsAInterpolar <- length(spObservaciones)
  # TO-DO: buscar la forma de ponerle los coeficientes a la fórmula
  interpolacionST$formulaRegresionCC <- fml
  rm(valoresVentana)
  return(interpolacionST)
}

universalGridding <- function(
    ti, coordsObservaciones, fechasObservaciones, valoresObservaciones, pathsRegresores=NULL, 
    valoresRegresoresSobreObservaciones=NULL, coordsAInterpolar, params, 
    iObservacionesEnCoordsAInterpolar=NULL, shpMask=NULL, longitudesEnColumnas=T, 
    paramsParaRellenoRegresores=NULL, pathsRegresoresParaRellenoRegresores=NULL) {
  # - pathsRegresores debe ser una matriz de nTxnU elementos, con una fila por fecha a interpolar y una columna por 
  # regresor. pathsRegresores[i, j] debe contener el path al archivo que contiene el raster con los datos del regresor j 
  # para la fecha i de las observaciones (pueden haber elementos repetidos en regresores[, j])
  # para ver la descripción de los demás parámetros ver la función universalGriddingEx de esta misma unidad
  # universalGriddingEx recibe los valores de los regresores sobre las coordenadas a interpolar para la fecha ti ya cargados
  # esta función los carga para facilitar su uso

  # Cargo los valores de los regresores sobre toda la grilla a interpolar para la fecha ti en una matriz nCIxnU
  # con las las filas variando según la coordenada a interpolar y las columnas variando según el regresor
  if (!is.null(pathsRegresores) && !is.null(valoresRegresoresSobreObservaciones) && ncol(pathsRegresores) > 0) {
    #valoresRegresoresSobreCoordsAInterpolar_ti <- extraerValoresRegresoresSobreSP(iInicial = ti, iFinal = ti, objSP = coordsObservaciones, pathsRegresores = pathsRegresores)
    valoresRegresoresSobreCoordsAInterpolar_ti <- matrix(NA, nrow=length(coordsAInterpolar), ncol=ncol(pathsRegresores))
    
    j <- 1
    for (j in 1:ncol(pathsRegresores)) {
      if (!is.na(pathsRegresores[ti, j])) {
        ext <- getFileExt(pathsRegresores[ti, j])
        if (ext == 'tif') { evaluarConReintentos(regresor <- readGDAL(pathsRegresores[ti, j], silent = T))
        } else if (ext == 'nc') {
          # TO-DO: receive varName as parameter
          evaluarConReintentos(regresor <- netCDFToSP(fname = pathsRegresores[ti, j]))
        } else { stop(paste0('extraerValorRegresorSobreSP: extensión no soportada "', pathsRegresores[ti, j], '"')) }
        
        if (!identicalCRS(coordsAInterpolar, regresor)) {
          aux <- spTransform(coordsAInterpolar, regresor@proj4string)
          valoresRegresoresSobreCoordsAInterpolar_ti[, j] <- over(x = aux, regresor, returnList = F)[,1]
          rm(aux)
        } else {
          valoresRegresoresSobreCoordsAInterpolar_ti[, j] <- over(x = coordsAInterpolar, regresor, returnList = F)[,1]
        }
      } else 
        valoresRegresoresSobreCoordsAInterpolar_ti[, j] <- NA
    }
    colnames(valoresRegresoresSobreCoordsAInterpolar_ti) <- colnames(pathsRegresores)
  } else { valoresRegresoresSobreCoordsAInterpolar_ti <- NULL }
  
  
  if (params$modoDiagnostico) {
    if (!is.null(valoresRegresoresSobreCoordsAInterpolar_ti)) {
      i <- 1
      if (gridded(coordsAInterpolar)) {
        for (i in (1:ncol(valoresRegresoresSobreCoordsAInterpolar_ti))) {
          escala <- darEscala(
            especificacion = params$especEscalaDiagnostico, 
            valores=valoresRegresoresSobreCoordsAInterpolar_ti[, i])
          spAux <- SpatialPixelsDataFrame(points = coordsAInterpolar, data = data.frame(value=valoresRegresoresSobreCoordsAInterpolar_ti[, i]))
          mapearGrillaGGPlot(grilla = spAux, shpBase = shpMask$shp, zcol=1, escala = escala, titulo = paste0(colnames(valoresRegresoresSobreCoordsAInterpolar_ti)[i], ' - ', params$strFecha), 
                             nomArchResultados = paste0(params$carpetaParaModoDiagnostico, '02.', i, '-Regresor-', colnames(valoresRegresoresSobreCoordsAInterpolar_ti)[i], '.png'), 
                             dibujar = F)
        }
      } else {
        for (i in (1:ncol(valoresRegresoresSobreCoordsAInterpolar_ti))) {
          escala <- darEscala(
            especificacion = params$especEscalaDiagnostico, 
            valores=valoresRegresoresSobreCoordsAInterpolar_ti[, i])
          spAux <- SpatialPointsDataFrame(coords = coordsAInterpolar, data = data.frame(value=valoresRegresoresSobreCoordsAInterpolar_ti[, i]))
          mapearPuntosGGPlot(puntos = spAux, shpBase = shpMask$shp, zcol=1, escala = escala, titulo = paste0(colnames(valoresRegresoresSobreCoordsAInterpolar_ti)[i], ' - ', params$strFecha), 
                             nomArchResultados = paste0(params$carpetaParaModoDiagnostico, '02.', i, '-Regresor-', colnames(valoresRegresoresSobreCoordsAInterpolar_ti)[i], '.png'), 
                             dibujar = F, dibujarTexto = T)
        }
      }
    }
  }

  if (params$rellenarRegresores && !is.null(valoresRegresoresSobreObservaciones)) {
    iNA <- is.na(valoresRegresoresSobreCoordsAInterpolar_ti[shpMask$mask])
    if (any(iNA)) {
      iObservaciones <- over(geometry(coordsObservaciones), coordsAInterpolar)
      if (params$invertirAjusteRegresores) {
        invRegs <- ajusteRegresores(ti = ti, coordsObservaciones = coordsObservaciones, fechasObservaciones = fechasObservaciones, 
                                    valoresObservaciones = valoresObservaciones, coordsAInterpolar=coordsAInterpolar, params=params, 
                                    valoresRegresoresSobreObservaciones=valoresRegresoresSobreObservaciones, 
                                    valoresRegresoresSobreCoordsAInterpolar_ti=valoresRegresoresSobreCoordsAInterpolar_ti,
                                    incorporarCoordenadas = params$incorporarCoordenadas, 
                                    formulaCoordenadas = params$formulaCoordenadas,
                                    incorporarTiempo=params$incorporarTiempo,
                                    formulaTiempo=params$formulaTiempo,
                                    incorporarDistanciaAlAgua = params$incorporarDistanciaAlAgua,
                                    formulaDistanciaAlAgua = params$formulaDistanciaAlAgua, 
                                    incorporarAltitud = params$incorporarAltitud,
                                    formulaAltitud = params$formulaAltitud, invertir=T)
        
        ventana <- getVentana(ti=ti, nT=nrow(valoresObservaciones), tamanioSemiVentana = params$ventanaIgualacionDistribuciones, tlagsAR = params$tlagsAR)
        
        regsInvertidos <- matrix(invRegs$valoresCampoBaseSobreObservaciones, nrow=length(ventana$tsVentana), ncol = ncol(valoresObservaciones))
        iNA2 <- is.na(valoresRegresoresSobreObservaciones[[1]][ventana$tsVentana,])
        valoresRegresoresSobreObservaciones[[1]][ventana$tsVentana,][iNA2] <- regsInvertidos[iNA2]
        
        iNA2 <- is.na(valoresRegresoresSobreCoordsAInterpolar_ti[iObservaciones, 1])
        
        #sp2 <- SpatialPixelsDataFrame(points = coordsAInterpolar, data = data.frame(valoresRegresoresSobreCoordsAInterpolar_ti))
        #mapearGrillaGGPlot(grilla = sp2, shpBase = shpMask$shp, zcol=1, continuo=T, titulo = paste0(fechasObservaciones[ti], ': LST'), 
        #                   nomArchResultados = 'Resultados/Ejemplos/RellenoRegresores/LST_Con_Huecos.png')
        
        valoresRegresoresSobreCoordsAInterpolar_ti[iObservaciones, 1][iNA2] <- regsInvertidos[ventana$iTiEnTsVentana,][iNA2]
        #valoresRegresoresSobreCoordsAInterpolar_ti[iObservaciones, 1] <- regsInvertidos[ventana$iTiEnTsVentana,]
        
        #sp3 <- SpatialPixelsDataFrame(points = coordsAInterpolar, data = data.frame(valoresRegresoresSobreCoordsAInterpolar_ti))
        #mapearGrillaGGPlot(grilla = sp3, shpBase = shpMask$shp, zcol=1, continuo=T, titulo = paste0(fechasObservaciones[ti], ': LST con regresión inversa'),
        #                   nomArchResultados = 'Resultados/Ejemplos/RellenoRegresores/LST_Con_Huecos_RegresionInversa.png')
      }
      
      sp <- SpatialPointsDataFrame(coords = coordsAInterpolar, data = data.frame(value=valoresRegresoresSobreCoordsAInterpolar_ti[,1]))
      # mapearGrillaGGPlot(grilla=SpatialPixelsDataFrame(points = coordsAInterpolar, data = data.frame(value=valoresRegresoresSobreCoordsAInterpolar_ti[,1])), shpBase = shpMask$shp, continuo = T, dibujar = F)
      
      #system.time(sp <- rellenarSP(sp = sp, mascara = shpMask$mask, metodo = paramsParaRellenoRegresores$interpolationMethod, nMuestras = 2000, 
      #                    nRepeticiones = 3, zcol = 1, params = paramsParaRellenoRegresores, 
      #                    pathsRegresores = pathsRegresoresParaRellenoRegresores[ti, ,drop=F],
      #                    minNCuadrantesDisponibles = 3, maxNMuestrasPorCuadrante=125))
      
      #mascara = shpMask$mask
      #metodo = paramsParaRellenoRegresores$interpolationMethod
      #nRepeticiones = 10
      #zcol = 1
      #distanciaMaximaMuestras=0
      #nCuadrantesX=4
      #nCuadrantesY=round(nCuadrantesX * diff(bbox(sp)[2,]) / diff(bbox(sp)[1,]))
      #nCuadrantesZ=5
      #nMuestras = 3 * nCuadrantesX * nCuadrantesY * nCuadrantesZ
      #iEsObligatoriosEnLaMuestra=NULL
      #params = paramsParaRellenoRegresores
      #pathsRegresores = pathsRegresoresParaRellenoRegresores[ti, ,drop=F]
      #minNCuadrantesDisponibles = 3
      #minNDatosDisponibles=30 
      #sustituirConPrimerRegresorSiNoHayDatosDisponibles=TRUE 
      #maxNMuestrasPorCuadrante = round(1.1 * nMuestras / (nCuadrantesX * nCuadrantesY))
      sp <- rellenarSP(sp = sp, mascara = shpMask$mask, metodo = paramsParaRellenoRegresores$interpolationMethod, 
                       nMuestras = 240, nRepeticiones = 10, zcol = 1, params = paramsParaRellenoRegresores, 
                       pathsRegresores = pathsRegresoresParaRellenoRegresores[ti, ,drop=F], minNCuadrantesDisponibles = 3)
      valoresRegresoresSobreCoordsAInterpolar_ti[, 1] <- sp@data[,1]

      # mapearGrillaGGPlot(grilla=SpatialPixelsDataFrame(points = coordsAInterpolar, data = sp@data), shpBase = shpMask$shp, continuo = T, dibujar = F)
      rm(sp)
      
      iObsNA <- is.na(valoresRegresoresSobreObservaciones[[1]][ti,])
      valoresRegresoresSobreObservaciones[[1]][ti, iObsNA] <- valoresRegresoresSobreCoordsAInterpolar_ti[iObservaciones][iObsNA]
    }
  }

  return(universalGriddingEx(
    ti = ti, coordsObservaciones = coordsObservaciones, fechasObservaciones=fechasObservaciones, 
    valoresObservaciones = valoresObservaciones, coordsAInterpolar = coordsAInterpolar, 
    params=params, valoresRegresoresSobreObservaciones=valoresRegresoresSobreObservaciones,
    valoresRegresoresSobreCoordsAInterpolar_ti=valoresRegresoresSobreCoordsAInterpolar_ti, 
    iObservacionesEnCoordsAInterpolar=iObservacionesEnCoordsAInterpolar, shpMask=shpMask, 
    longitudesEnColumnas = longitudesEnColumnas))
}

getPesosVentana <- function(tamanioVentana, iTiEnTsVentana, tsVentana, ventanaIgualacionDistribuciones, 
                            iFilasCompletas, nO, tipo=c('inversoDistanciaATi', 'potencia', 'lineales', 'gaussianos')) {
  # ventanaIgualacionDistribuciones <- 8
  # tamanioVentana <- 2 * ventanaIgualacionDistribuciones + 1
  # iTiEnTsVentana <- trunc(tamanioVentana / 2) + 1
  # tsVentana <- 1:tamanioVentana
  if (tamanioVentana > 1) {
    if (tipo[1] == 'inversoDistanciaATi') {
      aux <- 1 / (abs(iTiEnTsVentana - (tsVentana - tsVentana[1] + 1))+1)
    } else if (tipo[1] == 'gaussianos') {
      # Pesos Gaussianos
      delta <- 8 / (tamanioVentana - 1)
      aux <- dnorm(seq(from=-(iTiEnTsVentana - 1) * delta, to=(tamanioVentana - iTiEnTsVentana) * delta, length.out = tamanioVentana), sd=1)
    } else if (tipo[1] == 'potencia') {
      aux <- 0.7 ^ abs(iTiEnTsVentana - (tsVentana - tsVentana[1] + 1))
    } else {
      # Pesos lineales hasta 0.25
      aux <- (abs(iTiEnTsVentana - (tsVentana - tsVentana[1] + 1)) / ventanaIgualacionDistribuciones)
      aux <- 1 - 0.75 * aux / (max(aux))
    }
    
    aux <- (aux / sum(aux))
    pesos <- rep(aux, nO)
    #round(aux, 2)
    #plot(x=seq_along(aux), y = aux, col ='red')
    #points(x=seq_along(aux), y = aux, col= 'green')
    #points(x=seq_along(aux), y = aux, col= 'blue')
    #points(x=seq_along(aux), y = aux, col= 'violet')
    #points(x=seq_along(aux), y = aux, col= 'black')
  } else { pesos <- rep(1, nO) }
  pesos <- pesos[iFilasCompletas]
}

cachearRegresoresEstaticos <- function(coordsObservaciones, coordsAInterpolar, nCoresAUsar=0) {
  # Guardo en disco los valores de los regresores estaticos para las coordenadas de las observaciones
  # y las coordenadas a interpolar
  require(rgeos)
  
  cuerposDeAgua <- NULL
  altitudTerreno <- NULL
  
  parGDistance <- function(spgeom1, spgeom2=NULL, byid=FALSE, haussdorf=FALSE, densifyFrac=NULL, nCoresAUsar=0) {
    #spgeom1 <- as(coordsAInterpolar, 'SpatialPoints')
    #spgeom2 <- cuerposDeAgua
    #byid=T
    #haussdorf=FALSE
    #densifyFrac=NULL
    #i<-1
    
    auxFunc <- function(x, spgeom1, spgeom2, byid, haussdorf, densifyFrac) { 
      return(gDistance(spgeom1[x,], spgeom2, byid, haussdorf, densifyFrac)) 
    }
    if (nCoresAUsar <= 0) nCoresAUsar <- min(length(spgeom1), detectCores(T, T))
    
    if (nCoresAUsar > 1) {
      cl <- makeCluster(getOption("cl.cores", nCoresAUsar))
      clusterEvalQ(cl, {
        require('rgeos')}
      )
      iSplit <- clusterSplit(cl, 1:length(spgeom1))
      #rm(resultado)
      #resultado <- clusterApply(cl = cl, x = iSplit, fun = auxFunc, spgeom1=spgeom1, spgeom2=spgeom2, 
      #                          byid=byid, haussdorf=haussdorf, densifyFrac=densifyFrac)
      resultado <- parSapply(cl = cl, X = iSplit, FUN = auxFunc, spgeom1, spgeom2, byid, haussdorf, densifyFrac)
      stopCluster(cl)
      
      return(do.call(cbind, resultado))      
    } else { return(gDistance(spgeom1, spgeom2, byid, haussdorf, densifyFrac)) }
  }
  
  # Regresores estáticos sobre coordenadas a interpolar
  source(paste0(script.dir.interpolarEx, '../cacheFunciones/cacheFunciones.r'), encoding = 'WINDOWS-1252')
  pathCacheDatosCoordsAInterpolar <- getPathCache(objParametros = list(coordsAInterpolar, version=2))
  if (!file.exists(pathCacheDatosCoordsAInterpolar) | (file.info(pathCacheDatosCoordsAInterpolar)$size <= 0)) {
    dfDatosCoordsAInterpolar <- list()
    
    cInterp <- sp::coordinates(coordsAInterpolar)
    dfDatosCoordsAInterpolar$x <- cInterp[, 1]
    dfDatosCoordsAInterpolar$y <- cInterp[, 2]
    
    archiCuerposDeAgua <- paste0(script.dir.interpolarEx, 'datasets/cuerposDeAguaUy_25MKm2.shp')
    if (!file.exists(archiCuerposDeAgua)) stop(paste0('interpolarEx.r.cachearRegresoresEstaticos: No se encuentra el archvo ', archiCuerposDeAgua, sep =''))
    cuerposDeAgua <- cargarSHP(pathSHP = archiCuerposDeAgua)
    if (!identicalCRS(cuerposDeAgua, coordsAInterpolar)) {
      cuerposDeAgua <- spTransform(cuerposDeAgua, coordsAInterpolar@proj4string)
    }
    cuerposDeAgua <- gUnaryUnion(cuerposDeAgua)
    dfDatosCoordsAInterpolar$dist <- t(gDistance(spgeom1 = as(coordsAInterpolar, 'SpatialPoints'), spgeom2 = cuerposDeAgua, byid = T))[,1]

    archiAltitud <- paste0(script.dir.interpolarEx, 'datasets/GMTED2010_Mean_Uy.tif')
    if (!file.exists(archiAltitud)) stop(paste0('interpolarEx.r.cachearRegresoresEstaticos: No se encuentra el archvo ', archiAltitud, sep =''))
    altitudTerreno <- readGDAL(archiAltitud, silent=T)
    puntosCoordsAInterpolar <- spTransform(as(coordsAInterpolar, 'SpatialPoints'), altitudTerreno@proj4string)
    dfDatosCoordsAInterpolar$alt <- over(puntosCoordsAInterpolar, altitudTerreno)[,1]

    dfDatosCoordsAInterpolar <- as.data.frame(dfDatosCoordsAInterpolar)
    
    #lala <- SpatialPixelsDataFrame(points = coordsAInterpolar, data = dfDatosCoordsAInterpolar)
    #lala$dist2 <- (lala$dist ^ (1/2)) ^(1/3)
    #mapearGrillaGGPlot(lala, shpMask$shp, continuo = T, zcol=4)
    
    guardarCache(pathCache = pathCacheDatosCoordsAInterpolar, obj = dfDatosCoordsAInterpolar)
  }
  
  # Regresores estáticos sobre coordenadas de observaciones
  dfDatosCoordsObservaciones <- list()
  pathCacheDatosCoordsObservaciones <- getPathCache(objParametros = list(coordsObservaciones, version=2))
  if (!file.exists(pathCacheDatosCoordsObservaciones) | (file.info(pathCacheDatosCoordsObservaciones)$size <= 0)) {
    dfDatosCoordsObservaciones <- list()
    
    cObs <- sp::coordinates(coordsObservaciones)
    dfDatosCoordsObservaciones$x <- cObs[, 1]
    dfDatosCoordsObservaciones$y <- cObs[, 2]
    
    if (is.null(cuerposDeAgua)) {
      archiCuerposDeAgua <- paste0(script.dir.interpolarEx, 'datasets/cuerposDeAguaUy_25MKm2.shp')
      if (!file.exists(archiCuerposDeAgua))
        stop(paste0('interpolarEx.r.cachearRegresoresEstaticos: No se encuentra el archvo ', archiCuerposDeAgua, sep =''))
      cuerposDeAgua <- cargarSHP(archiCuerposDeAgua)
    } 
    if (!identicalCRS(cuerposDeAgua, coordsObservaciones)) {
      cuerposDeAgua <- spTransform(cuerposDeAgua, coordsAInterpolar@proj4string)
    }
    cuerposDeAgua <- gUnaryUnion(cuerposDeAgua)
    dfDatosCoordsObservaciones$dist <- t(gDistance(spgeom1 = as(coordsObservaciones, 'SpatialPoints'), spgeom2 = cuerposDeAgua, byid = T))[,1]
    
    if (is.null(altitudTerreno)) {
      archiAltitud <- paste0(script.dir.interpolarEx, 'datasets/GMTED2010_Mean_Uy.tif')
      if (!file.exists(archiAltitud))
        stop(paste0('interpolarEx.r.cachearRegresoresEstaticos: No se encuentra el archvo ', archiAltitud, sep =''))
      altitudTerreno <- readGDAL(archiAltitud, silent=T)
    }
      
    puntosObsAux <- spTransform(as(coordsObservaciones, 'SpatialPoints'), altitudTerreno@proj4string)
    dfDatosCoordsObservaciones$alt <- over(puntosObsAux, altitudTerreno)[,1]
    
    dfDatosCoordsObservaciones <- as.data.frame(dfDatosCoordsObservaciones)
    
    # lala <- SpatialPointsDataFrame(coords = coordsObservaciones, data = dfDatosCoordsObservaciones)
    # mapearPuntosConEtiquetasGGPlot(lala, shpMask$shp, zcol=1)
    
    guardarCache(pathCache = pathCacheDatosCoordsObservaciones, obj = dfDatosCoordsObservaciones)
  }
  return(NULL)
}

incorporarRegresoresEstaticos <- function(
    ti, coordsObservaciones, fechasObservaciones, valoresObservaciones, coordsAInterpolar, params, 
    valoresRegresoresSobreObservaciones=NULL, valoresRegresoresSobreCoordsAInterpolar_ti=NULL,
    incorporarCoordenadas=FALSE, formulaCoordenadas='x + y',
    incorporarTiempo=FALSE, formulaTiempo='t',
    incorporarDistanciaAlAgua=FALSE, formulaDistanciaAlAgua='I(dist^0.125)',
    incorporarAltitud=FALSE, formulaAltitud='alt') {
  #ventana <- getVentana(ti=ti, nT=nrow(valoresObservaciones), tamanioSemiVentana = params$ventanaIgualacionDistribuciones, tlagsAR = params$tlagsAR)

  source(paste0(script.dir.interpolarEx, '../cacheFunciones/cacheFunciones.r'), encoding = 'WINDOWS-1252')
  # El valor de version es para poder cambiar el código y forzar un recálculo de los caches
  
  geomCoordsInterp <- geometry(coordsAInterpolar)
  geomCoordsObs <- geometry(coordsObservaciones)
  pathCacheDatosCoordsAInterpolar <- getPathCache(objParametros = list(geomCoordsInterp, version=2))
  pathCacheDatosCoordsObservaciones <- getPathCache(objParametros = list(geomCoordsObs, version=2))
  if (!(file.exists(pathCacheDatosCoordsAInterpolar) & (file.info(pathCacheDatosCoordsAInterpolar)$size > 0)) |
      !(file.exists(pathCacheDatosCoordsObservaciones) & (file.info(pathCacheDatosCoordsObservaciones)$size > 0))) {
    cachearRegresoresEstaticos(geomCoordsObs, geomCoordsInterp, params$nCoresAUsar)
  }
  
  dfDatosCoordsAInterpolar <- cargarCache(pathCacheDatosCoordsAInterpolar)
  dfDatosCoordsObservaciones <- cargarCache(pathCacheDatosCoordsObservaciones)
  
  # Escalado y centrado de x e y para evitar problemas de precisión numérica
  if (incorporarCoordenadas && nrow(dfDatosCoordsAInterpolar) > 1) {
    minX <- min(dfDatosCoordsAInterpolar$x)
    rangeX <- diff(range(dfDatosCoordsAInterpolar$x))
    minY <- min(dfDatosCoordsAInterpolar$y)
    rangeY <- diff(range(dfDatosCoordsAInterpolar$y))
    maxRange <- max(rangeX, rangeY)
  
    dfDatosCoordsAInterpolar$x <- (dfDatosCoordsAInterpolar$x - minX) / maxRange
    dfDatosCoordsObservaciones$x <- (dfDatosCoordsObservaciones$x - minX) / maxRange
    dfDatosCoordsAInterpolar$y <- (dfDatosCoordsAInterpolar$y - minY) / maxRange
    dfDatosCoordsObservaciones$y <- (dfDatosCoordsObservaciones$y - minY) / maxRange
  }
  
  dfDatosCoordsAInterpolar <- as.data.frame(dfDatosCoordsAInterpolar)
  dfDatosCoordsObservaciones <- as.data.frame(dfDatosCoordsObservaciones)
  
  if (incorporarCoordenadas) { formulaCompleta <- formulaCoordenadas
  } else { formulaCompleta <- '' }
  if (incorporarDistanciaAlAgua) formulaCompleta <- paste(formulaCompleta, formulaDistanciaAlAgua, sep=' + ')
  if (incorporarAltitud) formulaCompleta <- paste(formulaCompleta, formulaAltitud, sep=' + ')
  
  terminos <- trim(strsplit(formulaCompleta, split = '+', fixed = T)[[1]])
  
  # - valoresRegresoresSobreObservaciones debe ser una lista de largo nU, cuyos elementos son matrices nTxnO, como valoresObservaciones, 
  # pero con los valores de los nU regresores en las ubicaciones de las observaciones. Su valor se puede obtener usando 
  # la funcion extraerValoresRegresoresSobreSP de esta misma unidad    

  if (F) {
    ihs <- function(y, theta=1/2) {
      return(asinh(theta * y) / theta)
    }
    
    par(mfrow=c(2,2),oma=c(0,0,2,0))
    for (i in 1:length(fechasObservaciones)) {
      plot(x=dfDatosCoordsObservaciones$x, y=valoresObservaciones[i,], xlab='x')
      plot(x=dfDatosCoordsObservaciones$y, y=valoresObservaciones[i,], xlab='y')
      plot(x=dfDatosCoordsObservaciones$dist, y=valoresObservaciones[i,], xlab='dist')
      plot(x=dfDatosCoordsObservaciones$alt, y=valoresObservaciones[i,], xlab='alt')
      title(fechasObservaciones[i], outer=TRUE)
    }
    par(mfrow=c(1,1))
    
    dfAux <- dfDatosCoordsObservaciones
    dfAux$temp <- valoresObservaciones[1,]
    instant_pkgs('car')
    #help(boxCox)
    dfAux <- dfAux[complete.cases(dfAux),]
    modelo <- lm(formula = 'temp~x+y+dist+alt', data = dfAux)
    bc <- boxCox(modelo, family="yjPower", lambda = seq(-10, 10, 0.1), plotit = T)
    lambda <- round(with(bc, x[which.max(y)]))
    depvar.transformed <- yjPower(dfAux, lambda)
  }
  
  terminos <- terminos[terminos!='']
  # Aplico la formula a los data frames
  regresoresEstaticosCO <- vapply(terminos, function(x) { eval(parse(text = x), envir = dfDatosCoordsObservaciones) }, FUN.VALUE = numeric(length(coordsObservaciones)))
  # Separo el dataframe en una lista con arrays 1xnO
  regresoresEstaticosCO <- split(regresoresEstaticosCO, rep(1:ncol(regresoresEstaticosCO), each = nrow(regresoresEstaticosCO)))
  
  # Agrego los regresores estáticos a valoresRegresoresSobreObservaciones  
  if (is.null(valoresRegresoresSobreObservaciones)) { valoresRegresoresSobreObservaciones <- list() }
  oldNu <- length(valoresRegresoresSobreObservaciones)
  length(valoresRegresoresSobreObservaciones) <- oldNu + length(regresoresEstaticosCO)
  
  # Para cada elemento de la lista repito el array anterior para volverlo nTxnO
  for (i in 1:length(regresoresEstaticosCO)) {
    valoresRegresoresSobreObservaciones[[oldNu + i]] <- matrix(regresoresEstaticosCO[[i]], ncol = ncol(valoresObservaciones), nrow = nrow(valoresObservaciones), byrow = T)
  }
  names(valoresRegresoresSobreObservaciones)[(oldNu+1):length(valoresRegresoresSobreObservaciones)] <- terminos
  
  if (incorporarTiempo) {
    terminosTiempo <- trim(strsplit(formulaTiempo, split = '+', fixed = T)[[1]])
    oldNU2 <- length(valoresRegresoresSobreObservaciones)
    length(valoresRegresoresSobreObservaciones) <- oldNU2 + length(terminosTiempo)
    tiempo <- data.frame(t=seq.int(from = 1, to = nrow(valoresObservaciones), by = 1))
    tiempo <- vapply(terminosTiempo, function(x) { eval(parse(text = x), envir = tiempo) }, FUN.VALUE = numeric(nrow(valoresObservaciones)))
    tiempo <- split(tiempo, rep(1:ncol(tiempo), each = nrow(tiempo)))
    for (i in 1:length(tiempo)) {
      valoresRegresoresSobreObservaciones[[oldNU2 + i]] <- matrix(tiempo[[i]], ncol = ncol(valoresObservaciones), nrow = nrow(valoresObservaciones), byrow = F)
    }
    names(valoresRegresoresSobreObservaciones)[(oldNU2+1):length(valoresRegresoresSobreObservaciones)] <- terminosTiempo
    
    terminos <- c(terminos, terminosTiempo)    
    dfDatosCoordsAInterpolar$t <- ti
  }
  
  # - valoresRegresoresSobreCoordsAInterpolar_ti debe ser una matrix nCxnU con los valores de todos los regresores (en columnas) sobre
  # todas las coordenadas a interpolar.
  regresoresEstaticosCI <- vapply(terminos, function(x) { eval(parse(text = x), envir = dfDatosCoordsAInterpolar) }, 
                                  FUN.VALUE = numeric(length(coordsAInterpolar)))
  if (!is.matrix(regresoresEstaticosCI)) regresoresEstaticosCI <- as.matrix(t(regresoresEstaticosCI))
  if (is.null(valoresRegresoresSobreCoordsAInterpolar_ti)) { 
    valoresRegresoresSobreCoordsAInterpolar_ti <- regresoresEstaticosCI
  } else { valoresRegresoresSobreCoordsAInterpolar_ti <- cbind(valoresRegresoresSobreCoordsAInterpolar_ti, regresoresEstaticosCI) }
  # Creo que no es necesario
  # colnames(valoresRegresoresSobreCoordsAInterpolar_ti)[(oldNu+1):length(valoresRegresoresSobreObservaciones)] <- terminos
  
  return(list(valoresRegresoresSobreObservaciones=valoresRegresoresSobreObservaciones,
              valoresRegresoresSobreCoordsAInterpolar_ti=valoresRegresoresSobreCoordsAInterpolar_ti))
}

ajusteRegresores <- function(
    ti, coordsObservaciones, fechasObservaciones, valoresObservaciones, coordsAInterpolar, params, 
    valoresRegresoresSobreObservaciones=NULL, valoresRegresoresSobreCoordsAInterpolar_ti=NULL,
    incorporarCoordenadas=FALSE, formulaCoordenadas='x + y', #formulaCoordenadas='I(x^2) + I(y^2) + I(x*y) + x + y',
    incorporarTiempo=FALSE, formulaTiempo='t',
    incorporarDistanciaAlAgua=FALSE, formulaDistanciaAlAgua='I(dist^0.125)',
    incorporarAltitud=FALSE, formulaAltitud='alt',
    descartarCoordenadasNoSignificativas=FALSE,
    invertir=FALSE, shpMask=NULL) {
  if (invertir && length(valoresRegresoresSobreObservaciones) == 1) {
    aux <- valoresRegresoresSobreObservaciones[[1]]
    valoresRegresoresSobreObservaciones[[1]] <- valoresObservaciones
    valoresObservaciones <- aux
  }
  
  # Paso el método a integer para que la comparación sea más rápida abajo
  metodoIgualacionDistribuciones <- which(c('ninguna', 'regresionLineal', 'regresionLinealRobusta', 
                                            'regresionLinealConEliminacionDeOutliers', 'CDFMatching',
                                            'Lasso', 'GLS') == params$metodoIgualacionDistribuciones)
  if (length(metodoIgualacionDistribuciones) == 0) {
    stop(paste0('interpolarEx.r.universalGridding: metodoIgualacionDistribuciones desconocido "', params$metodoIgualacionDistribuciones, '"')) }
  
  oldNu <- length(valoresRegresoresSobreObservaciones)
  origIncorporarCoordenadas <- incorporarCoordenadas
  if (metodoIgualacionDistribuciones == 7) {
    incorporarCoordenadas <- TRUE
  }
  
  if (incorporarCoordenadas | incorporarTiempo | incorporarDistanciaAlAgua | incorporarAltitud) {
    regs <- incorporarRegresoresEstaticos(
      ti = ti, coordsObservaciones = coordsObservaciones, fechasObservaciones = fechasObservaciones, 
      valoresObservaciones = valoresObservaciones, coordsAInterpolar = coordsAInterpolar, 
      params = params, valoresRegresoresSobreObservaciones = valoresRegresoresSobreObservaciones, 
      valoresRegresoresSobreCoordsAInterpolar_ti = valoresRegresoresSobreCoordsAInterpolar_ti,
      incorporarCoordenadas = incorporarCoordenadas, formulaCoordenadas = formulaCoordenadas, 
      incorporarTiempo = incorporarTiempo, formulaTiempo = formulaTiempo,
      incorporarDistanciaAlAgua = incorporarDistanciaAlAgua, formulaDistanciaAlAgua = formulaDistanciaAlAgua,
      incorporarAltitud = incorporarAltitud, formulaAltitud = formulaAltitud)
    valoresRegresoresSobreObservaciones <- regs$valoresRegresoresSobreObservaciones
    valoresRegresoresSobreCoordsAInterpolar_ti <- regs$valoresRegresoresSobreCoordsAInterpolar_ti
    
    seleccionRegresores <- descartarCoordenadasNoSignificativas
    if (seleccionRegresores) { 
      forceIn <- c(rep(1, oldNu), rep(0, ncol(valoresRegresoresSobreCoordsAInterpolar_ti) - oldNu))
    } else { 
      forceIn <- rep(1, ncol(valoresRegresoresSobreCoordsAInterpolar_ti))
    }
    
    if (!origIncorporarCoordenadas) {
      forceIn[names(valoresRegresoresSobreObservaciones) %in% c('x', 'y')] <- -1
    }
  } else { 
    forceIn <- rep(1, oldNu)
    seleccionRegresores <- FALSE
  }
  
  if (params$modoDiagnostico && !is.null(valoresRegresoresSobreCoordsAInterpolar_ti) && oldNu < ncol(valoresRegresoresSobreCoordsAInterpolar_ti)) {
    # i <- oldNu + 1
    # i <- i + 1
    for (i in (oldNu+1):ncol(valoresRegresoresSobreCoordsAInterpolar_ti)) {
      escala <- darEscala(
        especificacion = params$especEscalaDiagnostico, 
        valores=valoresRegresoresSobreCoordsAInterpolar_ti[, i])
      nomArchResultados <- paste0(params$carpetaParaModoDiagnostico, '02.', i, '-Regresor-', 
                                  gsub(pattern = '*', replacement = '', x = colnames(valoresRegresoresSobreCoordsAInterpolar_ti)[i], fixed = T), '.png')
      if (!file.exists(nomArchResultados)) {
        if (gridded(coordsAInterpolar)) {
          spAux <- SpatialPixelsDataFrame(points = coordsAInterpolar, data = data.frame(value=valoresRegresoresSobreCoordsAInterpolar_ti[, i]))
          mapearGrillaGGPlot(grilla = spAux, shpBase = shpMask$shp, zcol=1, escala = escala, titulo = paste0(colnames(valoresRegresoresSobreCoordsAInterpolar_ti)[i], ' - ', params$strFecha), 
                             nomArchResultados = nomArchResultados, dibujar = F)
        } else {
          spAux <- SpatialPointsDataFrame(coords = coordsAInterpolar, data = data.frame(value=valoresRegresoresSobreCoordsAInterpolar_ti[, i]))
          mapearPuntosGGPlot(puntos = spAux, shpBase = shpMask$shp, zcol=1, escala = escala, titulo = paste0(colnames(valoresRegresoresSobreCoordsAInterpolar_ti)[i], ' - ', params$strFecha), 
                             nomArchResultados = nomArchResultados, dibujar = F, dibujarTexto = T)
        }
      }
    }
  }

  nT <- nrow(valoresObservaciones)
  nO <- ncol(valoresObservaciones)
  nU <- length(valoresRegresoresSobreObservaciones)
  if (is.null(nU) || is.na(nU)) nU <- 0
  formulaRegresionCC <- ''
  if (nU > 0) {
    # Obtengo la ventana de igualacion de distribuciones
    ventana <- getVentana(ti=ti, nT=nT, tamanioSemiVentana = params$ventanaIgualacionDistribuciones, 
                          tlagsAR = params$tlagsAR)
    iTiEnTsVentana <- ventana$iTiEnTsVentana
    tsVentana <- ventana$tsVentana
    tamanioVentana <- length(tsVentana)

    # Obtengo los valores de los regresores sobre las observaciones para las fechas de la ventana
    valsRegresoresObservaciones <- matrix(data=NA, nrow = nO * tamanioVentana, ncol = nU)
    for (j in seq.int(1, nU, by=1)) { 
      valsRegresoresObservaciones[, j] <- as.numeric(valoresRegresoresSobreObservaciones[[j]][tsVentana, ])
    }
    colnames(valsRegresoresObservaciones) <- names(valoresRegresoresSobreObservaciones)
    
    varianzaRegresores <- apply(X = valsRegresoresObservaciones, MARGIN = 2, FUN = function(x) { return(var(x, na.rm=T)) })
    bRegresoresConVarianza <- !is.na(varianzaRegresores) & varianzaRegresores > 1E-6
    if (T) {
      n <- length(valsRegresoresObservaciones[, !bRegresoresConVarianza])
      valsRegresoresObservaciones[, !bRegresoresConVarianza] <- valsRegresoresObservaciones[, !bRegresoresConVarianza] + rnorm(n, sd=1e-3)
      iRegresoresConVarianza <- 1:ncol(valsRegresoresObservaciones)
    } else {
      iRegresoresConVarianza <- which(bRegresoresConVarianza)
    }

    # x <- valoresRegresoresSobreCoordsAInterpolar_ti[,1]
    iRegresoresNoNA_ti <- which(apply(X = valoresRegresoresSobreCoordsAInterpolar_ti, MARGIN = 2, FUN = function(x) { return(!all(is.na(x))) }))

    iRegresoresAConservar <- intersect(iRegresoresConVarianza, iRegresoresNoNA_ti)
    if (length(iRegresoresAConservar) < ncol(valsRegresoresObservaciones)) {
      forceIn <- forceIn[iRegresoresAConservar]
      valsRegresoresObservaciones <- valsRegresoresObservaciones[, iRegresoresAConservar, drop=F]
      if (!invertir) valoresRegresoresSobreCoordsAInterpolar_ti <- valoresRegresoresSobreCoordsAInterpolar_ti[, iRegresoresAConservar, drop=F]
    }
    
    if (ncol(valsRegresoresObservaciones) > 0) {
      # Si no es la regresión inversa, los indices de la observacion ti en la ventana
      # Si si es inversa toda la ventana porque se va a precisar
      if (!invertir) { iesVals <- iTiEnTsVentana + seq.int(from=0, to=nO-1, by=1) * tamanioVentana
      } else { iesVals <- 1:length(valoresObservaciones[tsVentana, ]) }
      
      # Igualacion de Distribuciones
      if (metodoIgualacionDistribuciones != 1 && any(!is.na(valoresObservaciones[tsVentana, ]))) {
        # Obtengo los valores de los regresores sobre las observaciones
        # Obtengo un vector concatenando los valores de las observaciones de cada estacion para cada fecha en la ventana 
        # algo como c(obsEst1, obsEst2, ..., obsEstN) donde obsEsti son los valores de todas las fechas en una estacion
        valoresObservacionesTsVentana <- as.numeric(valoresObservaciones[tsVentana, ])
        
        if (params$preECDFMatching && 
            max(valoresRegresoresSobreCoordsAInterpolar_ti, na.rm=T) / max(valoresObservacionesTsVentana, na.rm=T) >= 2) {
          regECDF <- ecdf(valoresRegresoresSobreCoordsAInterpolar_ti[, 1])
          valsRegresoresObservaciones[, 1] <- quantile(
            valoresObservacionesTsVentana, probs=regECDF(valsRegresoresObservaciones[, 1]), na.rm=T)
          iesVals <- iTiEnTsVentana + seq.int(from=0, to=nO-1, by=1) * tamanioVentana
          if (!invertir) {
            valoresRegresoresSobreCoordsAInterpolar_ti[, 1] <- quantile(
              valoresObservacionesTsVentana, 
              probs=regECDF(valoresRegresoresSobreCoordsAInterpolar_ti[, 1]), 
              na.rm = T)
          }  
        }
        
        # plot(x=sp::coordinates(coordsObservaciones)[,1], y=as.numeric(valoresObservaciones[ti, ]))
        # plot(x=sp::coordinates(coordsObservaciones)[,2], y=as.numeric(valoresObservaciones[ti, ]))
        # plot(x=valsRegresoresObservaciones[,1], y=valoresObservacionesTsVentana)
        # mapearPuntosGGPlot(puntos = coordsObservaciones, shpBase = shpMask$shp, continuo = T)
        
        if (metodoIgualacionDistribuciones != 5) {
          iFilasCompletas <- !is.na(valoresObservacionesTsVentana) & 
            !apply(valsRegresoresObservaciones, MARGIN = 1, FUN = function(x) { any(is.na(x)) })

          # Chequeo que tenga suficientes observaciones no nulas para hacer la regresión
          if (sum(iFilasCompletas) > nU + 1 && 
              nrow(unique(valsRegresoresObservaciones[iFilasCompletas, forceIn >= 0, drop=F])) >= nU + 1) {
            pesos <- getPesosVentana(tamanioVentana = tamanioVentana, iTiEnTsVentana = iTiEnTsVentana, tsVentana = tsVentana, 
                                     ventanaIgualacionDistribuciones = params$ventanaIgualacionDistribuciones, iFilasCompletas = iFilasCompletas, 
                                     nO = nO)
            
            if (metodoIgualacionDistribuciones %in% c(2, 3, 4, 7)) {
              # Regresion lineal
              df <- data.frame(value=valoresObservacionesTsVentana[iFilasCompletas], 
                               valsRegresoresObservaciones[iFilasCompletas, , drop=F])

              if (params$verbose) {
                nombresFilas <- rownames(valoresObservaciones[tsVentana, , drop=F])
                if (is.null(nombresFilas)) nombresFilas <- paste0('F', tsVentana)
                nombresColumnas <- colnames(valoresObservaciones)
                if (is.null(nombresColumnas)) nombresColumnas <- paste0('C', 1:ncol(valoresObservaciones))
                
                rownames(df) <- apply(expand.grid(nombresFilas, nombresColumnas),1,function(x) paste(x,collapse=" "))[iFilasCompletas]
                rm(nombresFilas, nombresColumnas)
              }
              
              if (oldNu == 1) {
                # Avoid extrapolation when adjusting regression.
                # If we've sampled less than params$minRatioRangosParaExtrapolacion of the (only)
                # regressors range, we'll add artificial samples with 0 error, taken straight from
                # the regressors unsampled values
                rSobreObs <- range(df[, 2])
                rSobreCoordsAInterpolar <- range(range(valoresRegresoresSobreCoordsAInterpolar_ti[, 1]))
                
                # Compare how much of the range of the regressor we have covered with the samples
                # at the observation locations. If we haven't covered at least
                # params$minRatioRangosParaExtrapolacion, complete the range trusting the regressor
                # ie, we'll add half as many samples as we already have aligned with the x = y line
                
                ratioMuestreado <- diff(rSobreObs) / diff(rSobreCoordsAInterpolar)
                if (ratioMuestreado < params$minRatioRangosParaExtrapolacion) {
                  iNoMuestreados <- which(
                      rSobreObs[1] > valoresRegresoresSobreCoordsAInterpolar_ti[, 1] |
                      rSobreObs[2] < valoresRegresoresSobreCoordsAInterpolar_ti[, 1])
                
                  ordenNoMuestreados <- order(valoresRegresoresSobreCoordsAInterpolar_ti[iNoMuestreados])
                  
                  nNuevasMuestras <- min(round(params$proporcionNuevasMuestras * (1 - ratioMuestreado) * nrow(df)), 
                                         length(iNoMuestreados))
                  iNuevasMuestras <- iNoMuestreados[ordenNoMuestreados[
                    as.integer(round(seq(0, 1, length.out = nNuevasMuestras) * (length(iNoMuestreados)-1) + 1))]]
                  
                  valsNuevasMuestras <- valoresRegresoresSobreCoordsAInterpolar_ti[
                    iNuevasMuestras, 
                    c(1, seq.int(1, ncol(valoresRegresoresSobreCoordsAInterpolar_ti))), 
                    drop=F]
                  colnames(valsNuevasMuestras) <- colnames(df)
                  df <- rbind(df, valsNuevasMuestras)
                  
                  #mapearPuntosGGPlot(puntos = SpatialPointsDataFrame(coords = coordsAInterpolar[iNoMuestreados,], data=data.frame(value=valoresRegresoresSobreCoordsAInterpolar_ti[iNoMuestreados, 1])), shpBase = shpBase)
                  #mapearPuntosGGPlot(puntos = SpatialPointsDataFrame(coords = coordsAInterpolar[iNuevasMuestras,], data=data.frame(value=valoresRegresoresSobreCoordsAInterpolar_ti[iNuevasMuestras, 1])), shpBase = shpBase)
                  #plot(x=df[,2], y=df[,1])
                  # TODO: calculate pesos for the actual time window, this assumes 
                  # params$ventanaIgualacionDistribuciones == 0
                  pesos <- c(pesos, rep(1, nNuevasMuestras))
                }
              }
              
              iUnicos <- !duplicated(df[, forceIn >= 0])
              df <- df[iUnicos, ]
              pesos <- pesos[iUnicos]

              terminoIndependiente='+1'
              if (params$formulaRegresion == '') {
                iColsRegresores <- (2:ncol(df))[forceIn >= 0]
                forceIn <- as.logical(forceIn[forceIn >= 0])
                formulas <- list(as.formula(
                  paste(colnames(df)[1], '~', paste('+', colnames(df)[iColsRegresores], collapse=''), 
                        terminoIndependiente)))
                
                if (seleccionRegresores) {
                  require(gtools)
                  terminos <- labels(terms(formulas[[1]]))
                  if (is.null(params$combinacionesDeRegresoresValidas)) {
                    terminosFijos <- paste(terminos[forceIn], collapse = '+')
                    if (terminosFijos != '') { terminosFijos <- paste(paste0(colnames(df)[1], '~', terminosFijos, '+', sep = ''))
                    } else { terminosFijos <- paste0(colnames(df)[1], '~', sep = '') }
                    terminosLibres <- terminos[!forceIn]
                    
                    #write(x = c('1', terminosLibres, '\n'), file = paste0('F:/Tesis/debugThreads/', Sys.getpid(), '.txt'), append = T)
                    
                    formulas <- list()
                    formulas[[1]] <- paste0(terminosFijos, '1')
                    for (i in seq_along(terminosLibres)) {
                      formulas <- c(formulas, 
                                    paste0(terminosFijos, apply(X = combinations(n = length(terminosLibres), r = i, v=terminosLibres, repeats.allowed = FALSE), MARGIN = 1, FUN = paste, collapse='+'), 
                                           terminoIndependiente))
                    }
                  } else {
                    formulas <- list()
                    terminosLibres <- terminos[(oldNu+1):nU]
                    
                    #write(x = c('2', terminosLibres, '\n'), file = paste0('F:/Tesis/debugThreads/', Sys.getpid(), '.txt'), append = T)
                    
                    for (j in seq_along(params$combinacionesDeRegresoresValidas)) {
                      terminosFijos <- paste(terminos[params$combinacionesDeRegresoresValidas[[j]]], collapse = '+')
                      terminosFijos <- paste(paste0(colnames(df)[1], '~', terminosFijos, '+', sep = ''))
                      
                      formulas[[length(formulas) + 1]] <- paste0(terminosFijos, '1')
                      for (i in seq_along(terminosLibres)) {
                        formulas <- c(formulas, 
                                      paste0(terminosFijos, apply(X = combinations(n = length(terminosLibres), r = i, v=terminosLibres, repeats.allowed = FALSE), MARGIN = 1, FUN = paste, collapse='+'), 
                                             terminoIndependiente))
                      }
                    }
                  }
                  rm(terminos, terminosFijos, terminosLibres)
                }
                #if (seleccionRegresores) {
                #  require('leaps')
                #  sumSel <- summary(regsubsets(x = formulaModelo, data = df, weights = pesos, intercept = terminoIndependiente == '+1'))
                #  # sumSel <- summary(regsubsets(x = as.matrix(df[,2:ncol(df)]), y=as.matrix(df[,1,drop=F]), weights = pesos, force.in = forceIn, intercept = terminoIndependiente == '+1'))
                #  nombresFormula <- colnames(df)[sumSel$which[which.min(sumSel$bic),]]
                #}
              } else {
                if (!grepl(pattern='~', x=params$formulaRegresion)) params$formulaRegresion <- paste0('y~', params$formulaRegresion)
                formulas <- list(as.formula(params$formulaRegresion))
              }
              
              if (metodoIgualacionDistribuciones == 2) {
                modelos <- lapply(X = formulas, FUN = function(x) { return(lm(formula = as.formula(x), data = df, weights = pesos)) })
              } else if (metodoIgualacionDistribuciones == 3) { 
                modelos <- lapply(X = formulas, FUN = function(x) { 
                  res <- tryExpr(modelo <- rlm(formula = as.formula(x), data = df, weights = pesos, maxit=50), silent = T)
                  if (res) { return(modelo) 
                  } else { return(NULL) } 
                })
              } else if (metodoIgualacionDistribuciones == 7) {
                # test <- leaps(y = lm(formula=formulaModelo, data = df, weights = pesos)$residuals, x = df[,2:ncol(df)], int = TRUE, nbest = 1)
                # formulaCorStruct <- paste0('~',paste(labels(terms(formulaModelo)), collapse = '+'),sep='')
                #formulaCorStruct <- as.formula(paste('~', paste('+', colnames(df)[2:ncol(df)], collapse=''), terminoIndependiente))

                modelos <- lapply(X = formulas, FUN = function(x) {
                  #print(x)
                  formulaCorStruct <- as.formula('~x+y')
                  x <- as.formula(x)
                  res1 <- tryExpr(modelo1 <- rms::Gls(model = x, data = df, correlation = corExp(form = formulaCorStruct, nugget = params$usarNugget), control = glsControl(maxIter = 1000)), silent = T)
                  res2 <- tryExpr(modelo2 <- rms::Gls(model = x, data = df, correlation = corSpher(form = formulaCorStruct, nugget = params$usarNugget), control = glsControl(maxIter = 1000)), silent = T)
                  if (res1) {
                    if (res2) {
                      if (AICcmodavg::AICc(modelo1) < AICcmodavg::AICc(modelo2)) { return(modelo1)
                      } else { return(modelo2) }
                    } else { return(modelo1) }
                  } else if (res2) { return(modelo2)
                  } else {
                    res <- tryExpr(modelo <- rlm(formula=x, data = df, weights = pesos), silent = T)
                    if (res) { return(modelo) 
                    } else { res <- tryExpr(modelo <- lm(formula=x, data = df, weights = pesos), silent = T) 
                      if (res) { return(modelo)
                      } else { return(NULL) }
                    }
                  }
                })
              } else {
                modelos <- lapply(X = formulas, FUN = function(x, data, weights) { 
                  modelo <- lm(formula=as.formula(x), data = df, weights = pesos)
                  w <- abs(rstudent(modelo)) < 3 & abs(cooks.distance(modelo)) < 4/nrow(modelo$model)
                  modelo <- update(modelo, weights=as.numeric(w))
                  return(modelo)
                })
              }
              
              modelosValidos <- !sapply(modelos, FUN= is.null)
              modelos <- modelos[modelosValidos]
              formulas <- formulas[modelosValidos]
              
              if (!is.null(params$signosValidosRegresores)) {
                # Filtro los modelos que hayan resultado con signos inválidos en alguno de sus regresores
                
                coeficientesValidosModelo <- function(modelo, signosValidosRegresores) {
                  coeficientes <- coefficients(modelo)
                  coeficientes <- coeficientes[!is.na(coeficientes)]
                  iCoeficientes <- match(names(coeficientes), names(signosValidosRegresores))
                  iCoeficientesNoNA <- !is.na(iCoeficientes)
                  return(all(sign(coeficientes[iCoeficientesNoNA]) == signosValidosRegresores[iCoeficientes[iCoeficientesNoNA]]))
                }
                
                modelosValidos <- sapply(modelos, coeficientesValidosModelo, signosValidosRegresores=params$signosValidosRegresores)
                modelos <- modelos[modelosValidos]
                formulas <- formulas[modelosValidos]
              }
                
              if (length(modelos) > 0) {
                if (length(modelos) == 1) { iModelo <- 1
                } else { iModelo <- which.min(lapply(modelos, FUN = AICcmodavg::AICc)) }
                formulaModelo <- as.formula(formulas[[iModelo]])
                modelo <- modelos[[iModelo]]
                coeficientes <- coefficients(modelo)
                
                if (length(labels(terms(formulaModelo))) > 0) { modeloAjustado <- TRUE
                } else { modeloAjustado <- FALSE }
                rm(modelos)
              } else {
                modeloAjustado <- FALSE
              }
                
              if (modeloAjustado) {
                coeficientes <- coefficients(modelo)
                
                if (metodoIgualacionDistribuciones == 7 && !incorporarCoordenadas) {
                  newdata <- data.frame(valsRegresoresObservaciones[iesVals, , drop=F],
                                   x=sp::coordinates(coordsObservaciones)[iesVals, 1],
                                   y=sp::coordinates(coordsObservaciones)[iesVals, 2])
                } else {
                  newdata <- data.frame(valsRegresoresObservaciones[iesVals, , drop=F])
                }
                
                valoresCampoBaseSobreObservaciones <- predict(
                  object = modelo, newdata = newdata, na.action=na.pass)
                if (!invertir) {
                  if (metodoIgualacionDistribuciones == 7 && !incorporarCoordenadas) {
                    newdata <- data.frame(valoresRegresoresSobreCoordsAInterpolar_ti,
                                          x=sp::coordinates(coordsAInterpolar)[, 1],
                                          y=sp::coordinates(coordsAInterpolar)[, 2])
                  } else {
                    newdata <- data.frame(valoresRegresoresSobreCoordsAInterpolar_ti)
                  } 
                  valoresCampoBase <- predict(
                   object = modelo, newdata = newdata, na.action=na.pass)
                }
              } else {
                # si la formula no tiene terminos (regresores), el mejor predictor es la media de las observaciones
                coeficientes <- mean(valoresObservacionesTsVentana, na.rm=T)
                names(coeficientes) <- 'Intercept'
                  
                if (params$interpolationMethod == 'automap') {
                  # Si la interpolación se hace con Kriging uso Kriging Ordinario, sin regresores
                  valoresCampoBaseSobreObservaciones <- NULL
                  if (!invertir) valoresCampoBase <- NULL
                } else {
                  valoresCampoBaseSobreObservaciones <- rep(coeficientes, length(coordsObservaciones))
                  if (!invertir) valoresCampoBase <- rep(coeficientes, length(coordsAInterpolar))
                }
              }
                
              # plot(x=valsRegresoresObservaciones[,1], y=valoresObservacionesTsVentana)
              # plot(x=valoresRegresoresSobreObservaciones[[1]][ti,], y=valoresObservaciones[ti,])
              # plot(x=valoresCampoBaseSobreObservaciones, y=valoresObservaciones[ti, ])
              # abline(reg = modelo)
              
              # mean(valoresCampoBaseSobreObservaciones)
              # mean(valoresObservaciones[ti, ], na.rm=T)
              
              # cbind(as.numeric(valoresObservaciones[ti, ]), valoresCampoBaseSobreObservaciones, valsRegresoresObservaciones[seq.int(from=1, to=nO, by=1) + nO * (iTiEnTsVentana - 1),iRegresoresNoNA_ti])
              # plot(valoresObservaciones[ti, ], valoresCampoBaseSobreObservaciones)
              # dfbetas(modelo)
              # cooks.distance(modelo) > 4 / (length(modelo$residuals) -  modelo$df.residual - 1)
              
              if (params$modoDiagnostico) {
                png(paste0(params$carpetaParaModoDiagnostico, '03.1-Scatterplots.png'), width = 1080, height = 1080)
                pairs(df)
                dev.off()
                
                if (metodoIgualacionDistribuciones != 7 && modeloAjustado) {
                  tryCatch({
                    png(paste0(params$carpetaParaModoDiagnostico, '03.2-GraficosDiagnosticoLM.png'), width = 1080, height = 1080)
                    par(mfrow=c(2,2)) # Change the panel layout to 2 x 2
                    plot(modelo)
                    par(mfrow=c(1,1)) # Change back to 1 x 1
                    dev.off()
                  }, error = function(e) {
                    dev.off()
                  })
                  
                  capture.output(summary(modelo), file = paste0(params$carpetaParaModoDiagnostico, '03.3-SummaryLM.txt'))
                }
              }
              
              rm(iesVals, df, pesos, iUnicos)
              if (metodoIgualacionDistribuciones != 7 && modeloAjustado) rm(modelo, formulaModelo)
            } else { 
              # metodoIgualacionDistribuciones == 6
              # Lasso
              require('glmnet')
              cvfit <- cv.glmnet(x = valsRegresoresObservaciones[iFilasCompletas, , drop=F], y = valoresObservacionesTsVentana[iFilasCompletas], weights = pesos)
              coeficientes = as.matrix(coef(cvfit, s = "lambda.min"))[,1]
              valoresCampoBaseSobreObservaciones <- as.numeric(predict(cvfit, newx = valsRegresoresObservaciones[iesVals, , drop=F], s = "lambda.min"))
              if (!invertir) valoresCampoBase <- as.numeric(predict(cvfit, newx = valoresRegresoresSobreCoordsAInterpolar_ti))
            }
            
            if (tamanioVentana > 1) {
              # Corrección (cruda) de sesgo introducido por la ventana
              dif <- na.omit(valoresObservaciones[ti, ] - valoresCampoBaseSobreObservaciones)
              if (length(dif) >= 8) {
                # Si tengo al menos 8 observaciones con sus predicciones en base a los regresores en la fecha ti, 
                # como usamos la ventana para calibrar, resto el ME al termino independiente del modelo como forma cruda de eliminar el sesgo
                sesgo <- mean(dif)
                if (abs(sesgo) > 0.1) {
                  coeficientes["(Intercept)"] <- coeficientes["(Intercept)"] + sesgo
                  #mean(valoresObservaciones[ti, ], na.rm=T)
                  #mean(valoresCampoBaseSobreObservaciones)
                  #mean(valoresCampoBase, na.rm=T)
                  valoresCampoBaseSobreObservaciones <- valoresCampoBaseSobreObservaciones + sesgo
                  if (!invertir) valoresCampoBase <- valoresCampoBase + sesgo
                }
              }
            }
            formulaRegresionCC <- formulaStr(coeficientes = coeficientes[!is.na(coeficientes)], nDecimales = 2)
          } else {
            # Si no hay suficientes observaciones no nulas para estimar los coeficientes de los 
            # regresores hago la interpolacion sin campo base
            valoresCampoBase <- NULL
            valoresCampoBaseSobreObservaciones <- NULL
          }
        } else { # if (metodoIgualacionDistribuciones == 5) {
          # CDF Matching. El 1 es porque para CDF matching solo se usa el primer regresor
          iPrimerRegresor <- iRegresoresNoNA_ti[1]
          # iNoNA <- !is.na(valoresObservacionesTsVentana)
        
          #cdfMatching <- function(target_obs, source_obs) {
          #  sourceECDF <- ecdf(source_obs)
          #  return(quantile(target_obs, sourceECDF(source_obs), na.rm=T))
          #}
          
          # regECDF <- ecdf(valsRegresoresObservaciones[, iPrimerRegresor])
          regECDF <- ecdf(valoresRegresoresSobreCoordsAInterpolar_ti[, iPrimerRegresor])
          valoresCampoBaseSobreObservaciones <- quantile(
            valoresObservacionesTsVentana, probs=regECDF(valsRegresoresObservaciones[, iPrimerRegresor]), na.rm=T)
          iesVals <- iTiEnTsVentana + seq.int(from=0, to=nO-1, by=1) * tamanioVentana
          if (!invertir) {
            valoresCampoBase <- quantile(
              valoresObservacionesTsVentana, 
              probs=regECDF(valoresRegresoresSobreCoordsAInterpolar_ti[, iPrimerRegresor]), 
              na.rm = T)
          }
          
          #plot(x=valoresRegresoresSobreCoordsAInterpolar_ti[, iPrimerRegresor], y=valoresCampoBase)
          #plot(x=valoresCampoBaseSobreObservaciones, y=valoresObservacionesTsVentana)
            
          rm(regECDF, iesVals)
        }
        rm(valoresObservacionesTsVentana)
      } else {
        # Si no hay igualacion de distribuciones uso el valor del primer predictor que tenga algún valor no NA como campo base
        iPrimerRegresor <- iRegresoresNoNA_ti[1]
        valoresCampoBaseSobreObservaciones <- valsRegresoresObservaciones[iesVals, iPrimerRegresor]
        class(valoresCampoBaseSobreObservaciones) <- 'numeric'
        if (!invertir) { 
          valoresCampoBase <- valoresRegresoresSobreCoordsAInterpolar_ti[, iPrimerRegresor]
          class(valoresCampoBase) <- 'numeric'
        }
      }
      
      # plot(x=valoresCampoBaseSobreObservaciones, y=valoresObservacionesTi)
      # abline(modelo)
    } else {
      # Si todos los regresores son nulos hago la interpolacion sin campo base
      valoresCampoBase <- NULL
      valoresCampoBaseSobreObservaciones <- NULL
    }
    
    if (invertir) valoresCampoBase <- NULL
    rm(valsRegresoresObservaciones)
  } else {
    # Si no hay regresores hago la interpolacion sin campo base
    valoresCampoBase <- NULL
    valoresCampoBaseSobreObservaciones <- NULL
  }
  
  return(list(valoresCampoBase=valoresCampoBase, valoresCampoBaseSobreObservaciones=valoresCampoBaseSobreObservaciones,
              formulaRegresionCC=formulaRegresionCC))
}

universalGriddingEx <- function(ti, coordsObservaciones, fechasObservaciones, valoresObservaciones, coordsAInterpolar, params, 
                                valoresRegresoresSobreObservaciones=NULL, valoresRegresoresSobreCoordsAInterpolar_ti=NULL, 
                                iObservacionesEnCoordsAInterpolar=NULL, shpMask=NULL, longitudesEnColumnas=T) {
  params <- setMinMaxVal(observacionesValue=as.numeric(valoresObservaciones[ti,]), params)
  if (params$interpolationMethod != 'stUniversalKriging') {
    # La idea es poder pedir el mapeo del campo de una variable usando una serie de campos auxiliares
    # grillados como predictores de su media
    
    # El grillado se aplica en 3 pasos:
    # - Igualación de Distribuciones: puede ser ninguna, regresionLineal o CDFMatching
    # - Interpolación de Residuos: Puede ser ninguna o alguno de los métodos de interpolación (IDW, Kriging, etc)
    # - Postprocesamiento: Por ahora puede ser ninguna o la aplicación de la máscara de lluvia/no lluvia
    # Si todos los pasos son nada se retorna la grilla del primer regresor disponible si hay o un campo de NAs en su defecto.
    # Si no hay regresores el método se reduce a la interpolación habitual
    
    # Se quiere obtener Z(x, y, ti) = F(Z(xi, yi, ti), U1(x, y, tu1i), U2(x, y, tu2i), ..., Un(x, y, tuni))
    # Donde Z es el campo a interpolar, ti son los tiempos conocidos de Z, (xi, yi) son las posiciones conocidas de Z,
    # tji es el tiempo disponible del regresor j correspondiente a ti, (x, y) pertenciente a (X,Y) es una posición arbitraria, 
    # Uj son una serie de rasters auxiliares sobre todo el dominio de (X,Y) y F es la composición de los 3 pasos anteriores
    
    # Si la igualación de distribuciones es igualación de CDF, n tiene que ser 1 (solo un raster auxiliar)
    # La máscara de lluvia/no lluvia solo aplica a precipitación en principio pero
    # podría servir para detección de eventos extremos de otras variables
    
    # Sea:
    # - nO el número de puntos de observacion (xi, yi)
    # - nT el numero de puntos de tiempo ti
    # - nU el número de regresores disponibles
    # - nC el número de coordenadas a interpolar
    # Entonces:
    # - coordsObservaciones debe ser un objeto spatialPoints* con nO puntos, conteniendo la proyeccion y las coordenadas de los puntos de observacion
    # - valoresObservaciones debe ser una matriz nTxnO, con una fila por fecha a interpolar y una columna por punto de observacion
    # - valoresRegresoresSobreObservaciones debe ser una lista de largo nU, cuyos elementos son matrices nTxnO, como valoresObservaciones, 
    # pero con los valores de los nU regresores en las ubicaciones de las observaciones. Su valor se puede obtener usando 
    # la funcion extraerValoresRegresoresSobreSP de esta misma unidad
    # - valoresRegresoresSobreCoordsAInterpolar_ti debe ser una matrix nCxnU con los valores de todos los regresores (en columnas) sobre
    # todas las coordenadas a interpolar.
    # - pathsRegresores debe ser una matriz de nTxnU elementos, con una fila por fecha a interpolar y una columna por 
    # regresor. pathsRegresores[i, j] debe contener el path al archivo que contiene el raster con los datos del regresor j 
    # para la fecha i de las observaciones (pueden haber elementos repetidos en regresores[, j])
    # - coordsAInterpolar debe ser un objeto spatialPoints* o spatialPixels*, según si se quiere interpolar puntos o grillas,
    # con la proyeccion y las coordenadas de los puntos a interpolar
    # params debe ser un objeto creado con createParamsUniversalGridding
    # params$descartarCoordenadasNoSignificativas <- FALSE

    regs <- ajusteRegresores(
      ti = ti, coordsObservaciones = coordsObservaciones, fechasObservaciones = fechasObservaciones, 
      valoresObservaciones = valoresObservaciones, coordsAInterpolar=coordsAInterpolar, 
      params=params, valoresRegresoresSobreObservaciones=valoresRegresoresSobreObservaciones, 
      valoresRegresoresSobreCoordsAInterpolar_ti=valoresRegresoresSobreCoordsAInterpolar_ti,
      incorporarCoordenadas = params$incorporarCoordenadas, 
      formulaCoordenadas = params$formulaCoordenadas,
      incorporarTiempo = params$incorporarTiempo,
      formulaTiempo = params$formulaTiempo,
      incorporarDistanciaAlAgua = params$incorporarDistanciaAlAgua,
      formulaDistanciaAlAgua = params$formulaDistanciaAlAgua, 
      incorporarAltitud = params$incorporarAltitud,
      formulaAltitud = params$formulaAltitud, 
      descartarCoordenadasNoSignificativas = params$descartarCoordenadasNoSignificativas,
      invertir = FALSE,
      shpMask = shpMask)
    # mean(valoresObservaciones, na.rm = T) + c(-3, 3) * sd(valoresObservaciones, na.rm = T)
    # range(valoresObservaciones, na.rm = T)
    
    # Interpolacion
    coordsObservaciones$value <- as.numeric(valoresObservaciones[ti, ])
    params$formulaRegresionCC <- regs$formulaRegresionCC
    if (params$umbralMascaraCeros > 0) {
      # We store these values to be able to use the regressors in the Rain/NoRain mask interpolation
      params$valoresRegresoresSobreObservaciones <- lapply(
        valoresRegresoresSobreObservaciones, FUN = function(x) return(x[ti, , drop=F]))
      params$valoresRegresoresSobreCoordsAInterpolar_ti <- valoresRegresoresSobreCoordsAInterpolar_ti
    }
    
    # observaciones <- coordsObservaciones
    # longitudesEnColumnas <- T
    # valoresCampoBaseSobreObservaciones <- regs$valoresCampoBaseSobreObservaciones
    # valoresCampoBase <- regs$valoresCampoBase
    interpolacion <- interpolarEx(
      observaciones = coordsObservaciones, coordsAInterpolar = coordsAInterpolar, params = params, 
      shpMask = shpMask, valoresCampoBaseSobreObservaciones=regs$valoresCampoBaseSobreObservaciones, 
      valoresCampoBase=regs$valoresCampoBase, longitudesEnColumnas=longitudesEnColumnas)
    interpolacion$formulaRegresionCC <- regs$formulaRegresionCC
    
    # mapearPuntosGGPlot(puntos = observaciones, shpBase = shpMask$shp, continuo=T, zcol='value')
    # aux <- interpolacion$predictions
    # aux$var1.pred<- regs$valoresCampoBase[shpMask$mask]
    # mapearGrillaGGPlot(grilla = aux, shpBase = shpMask$shp, continuo = T)
    # mapearGrillaGGPlot(grilla = interpolacion$predictions, shpBase = shpMask$shp, continuo=T)
    rm(regs)
  } else {
    # params <- paramsIyM
    if (params$usarNugget) { fixNugget <- NA
    } else { fixNugget <- 0 }
    
    interpolacion <- stUniversalKrigingEx(ti=ti, spObservaciones=coordsObservaciones, fechasObservaciones=fechasObservaciones, valoresObservaciones=valoresObservaciones,
                                          coordsAInterpolar=coordsAInterpolar, valoresRegresoresSobreObservaciones=valoresRegresoresSobreObservaciones, 
                                          valoresRegresoresSobreCoordsAInterpolar_ti = valoresRegresoresSobreCoordsAInterpolar_ti, shpMask=shpMask, 
                                          modelosVariograma=params$modelosVariograma, cutoff=params$cutoff, tlags=params$tlags, nTsST=params$nTsST, tlagsAR=params$tlagsAR, 
                                          fixNugget=fixNugget, tryFixNugget=params$tryFixNugget, nPuntosIniciales=params$nPuntosIniciales, 
                                          modelosVariogramaST = params$modelosVariogramaST, ventanaIgualacionDistribuciones=params$ventanaIgualacionDistribuciones,
                                          verbose=params$verbose, params=params)
  }

  if (nrow(interpolacion$predictions) > 1 && !is.null(params$tlagsAR)) {
    interpolacion$predictions <- rellenarSP(sp = interpolacion$predictions, metodo = 'multiTps')
  }
  
  # interpolacion$predictions@data[,interpolacion$campoMedia] es para seleccionar los datos independiente del método, interpolacion$campoMedia es mean para copula y var1.pred para automap
  # force predictions to be in given range
  campoMedia <- interpolacion$campoMedia
  if (!is.na(params$minVal)) interpolacion$predictions@data[,campoMedia][interpolacion$predictions@data[,campoMedia] < params$minVal] <- params$minVal
  if (!is.na(params$maxVal)) interpolacion$predictions@data[,campoMedia][interpolacion$predictions@data[,campoMedia] > params$maxVal] <- params$maxVal
  
  return (interpolacion)
}

universalGriddingCV_i <- function(iObservacion=1, tIni=1, tFin=nrow(valoresObservaciones), coordsObservaciones, fechasObservaciones, 
                                  valoresObservaciones, params, valoresRegresoresSobreObservaciones=NULL, longitudesEnColumnas=T, 
                                  eliminarSerieTemporalCompleta=T, estimarNAs=F) {
  # valoresRegresoresSobreObservaciones lista de matrices con valoresRegresoresSobreObservaciones[[i]][j, k] los valores del regresor i en la fecha j estacion k
  # valoresRegresoresSobreCoordsAInterpolar_ti matriz con los valores de los regresores para la fecha ti. valoresRegresoresSobreCoordsAInterpolar_ti[i, j] tiene el valor de la estacion i, regresor j
  
  # iObservacion=3
  # tIni=1
  # tFin=nrow(valoresObservaciones)
  # estimarNAs=F
  iNAs <- !estimarNAs & is.na(valoresObservaciones[,iObservacion])
  
  ts <- seq.int(from = tIni, to=tFin, by = 1)
  res <- valoresObservaciones[ts, iObservacion, drop=F]
  res[] <- NA

  coordsAInterpolar_i <- coordsObservaciones[iObservacion, ]
  if (eliminarSerieTemporalCompleta) { valoresObservaciones[, iObservacion] <- NA 
  } else { oldValores <- valoresObservaciones[, iObservacion] }
  
  if (!is.null(valoresRegresoresSobreObservaciones)) {
    valoresRegresoresSobreCoordsAInterpolar_ti <- matrix(data=NA, nrow = 1, ncol = length(valoresRegresoresSobreObservaciones))
    colnames(valoresRegresoresSobreCoordsAInterpolar_ti) <- names(valoresRegresoresSobreObservaciones)
  } else { valoresRegresoresSobreCoordsAInterpolar_ti <- NULL }
  
  # valoresRegresoresSobreObservaciones[[1]][1,, drop=F]
  # Creo una máscara dummy para que no la cree (vacía) todos los pasos de tiempo
  shpMask <- cargarSHPYObtenerMascaraParaGrilla(pathSHP='', proj4strSHP='', grilla=coordsAInterpolar_i)
  
  
  # ti <- which(fechasObservaciones==as.POSIXct('2019-07-14', tz=tz(fechasObservaciones[1])))
  # ti <- ts[302]
  # ti <- 365
  for (ti in ts) {
    print(c(iObservacion, ti))
    if (!iNAs[ti]) {
      if (!is.null(valoresRegresoresSobreCoordsAInterpolar_ti)) {
        for (j in seq_along(along.with = valoresRegresoresSobreObservaciones))
          valoresRegresoresSobreCoordsAInterpolar_ti[1, j] <- valoresRegresoresSobreObservaciones[[j]][ti, iObservacion]
      }
      
      if (!eliminarSerieTemporalCompleta) valoresObservaciones[ti, iObservacion] <- NA
      
      # coordsAInterpolar = coordsAInterpolar_i
      # iObservacionesEnCoordsAInterpolar=NULL
      interp <- universalGriddingEx(
        ti = ti, coordsObservaciones = coordsObservaciones, fechasObservaciones=fechasObservaciones, 
        valoresObservaciones = valoresObservaciones, coordsAInterpolar = coordsAInterpolar_i, 
        params = params, valoresRegresoresSobreObservaciones = valoresRegresoresSobreObservaciones,
        valoresRegresoresSobreCoordsAInterpolar_ti = valoresRegresoresSobreCoordsAInterpolar_ti,
        iObservacionesEnCoordsAInterpolar=NULL, shpMask = shpMask, 
        longitudesEnColumnas=longitudesEnColumnas)
      
      if (!eliminarSerieTemporalCompleta) valoresObservaciones[ti, iObservacion] <- oldValores[ti]
      
      res[ti - tIni + 1] <- interp$predictions@data[1, interp$campoMedia]
    } else {
      res[ti - tIni + 1] <- NA
    }
  }
  
  return (res)
}

universalGriddingCV <- function(
    coordsObservaciones, fechasObservaciones, valoresObservaciones, params, pathsRegresores=NULL, 
    longitudesEnColumnas=T, iesAEstimar=1:ncol(valoresObservaciones), 
    eliminarSerieTemporalCompleta=TRUE, estimarNAs=FALSE) {
  # longitudesEnColumnas=T
  # iesAEstimar=1:ncol(valoresObservaciones)
  # eliminarSerieTemporalCompleta=TRUE
  # eliminarSerieTemporalCompleta=FALSE
  # estimarNAs=FALSE
  
  coordsObservaciones <- geometry(coordsObservaciones)
  # Hago un override de los parámetros que son para interpolación de píxeles para que queden para interpolación de puntos. 
  params$block <- NA
  params$coordsAInterpolarSonGrilla <- FALSE
  params$proj4StringAInterpolar <- proj4string(coordsObservaciones)
  params$pathSHPMapaBase <- ''
  params$modoDiagnostico <- FALSE
  if (!is.null(pathsRegresores)) { valoresRegresoresSobreObservaciones <- extraerValoresRegresoresSobreSP(objSP = coordsObservaciones, pathsRegresores = pathsRegresores)
  } else { valoresRegresoresSobreObservaciones <- NULL }

  if (params$nCoresAUsar <= 0) { 
    nCoresAUsar <- min(getAvailableCores(maxCoresPerGB = 1), ncol(valoresObservaciones))
  } else { nCoresAUsar <- params$nCoresAUsar }
  
  if (length(params$tlagsAR) <= 0) {
    if (nCoresAUsar > 1) {
      if (F) {
        # Código para testear el bug de library usando procesos en R 3.4.0
        # Si se usa R 3.4.0, las llamadas a library en instant_pkgs a veces fallan y hacen que no ejecute
        # bien el source de interpolarEx.r en algunos de los procesos, lo que hace que esos procesos
        # no tengan la definición de todas las funciones necesarias y falla el universalGriddingCV_i
        setwd('D:/Workspace/MCH2/MCH/MCHLib/R/Scripts')
        system('D:/Workspace/MCH2/MCH/MCHLib/R/Scripts/deployScriptsSoloLocal.bat')
        setwd('C:/mch/scripts/R/interpolar')
        
        while (T) {
          print(1)
          unlink(dir('D:/testsMCH/instantPkgs', full.names = T))
          cl <- makeCluster(getOption('cl.cores', nCoresAUsar))
          clusterExport(cl, varlist = c('script.dir.interpolarEx'))
          clusterEvalQ(cl = cl, expr = { 
            source(paste0(script.dir.interpolarEx, 'interpolarEx.r'), encoding = 'WINDOWS-1252') 
            if (exists(x = 'setMKLthreads')) { setMKLthreads(1) }
          })
          #clusterExport(cl, varlist = c('coordsObservaciones'))
          #clusterEvalQ(cl = cl, expr = { cargarSHPYObtenerMascaraParaGrilla(pathSHP = '', proj4strSHP = '', grilla = coordsObservaciones) })
          stopCluster(cl)
        }
      }

      cl <- makeCluster(getOption('cl.cores', nCoresAUsar))
      clusterExport(cl, varlist = c('script.dir.interpolarEx'))
      clusterEvalQ(cl = cl, expr = { 
        source(paste0(script.dir.interpolarEx, 'interpolarEx.r'), encoding = 'WINDOWS-1252') 
        if (exists(x = 'setMKLthreads')) { setMKLthreads(1) }
      })
      res <- valoresObservaciones
      res[,] <- parSapplyLB(
          cl=cl, X=iesAEstimar, FUN=universalGriddingCV_i, coordsObservaciones=coordsObservaciones, 
          fechasObservaciones=fechasObservaciones, valoresObservaciones=valoresObservaciones, params=params, 
          valoresRegresoresSobreObservaciones=valoresRegresoresSobreObservaciones, 
          longitudesEnColumnas=longitudesEnColumnas, eliminarSerieTemporalCompleta=eliminarSerieTemporalCompleta,
          estimarNAs=estimarNAs)
      stopCluster(cl)
    } else {
      res <- valoresObservaciones
      res[,] <- NA
      
      # iObservacion <- which(colnames(valoresObservaciones) == 'PRESA.BAYGORRIA.RHT')
      # iObservacion <- iesAEstimar[1]
      for (iObservacion in iesAEstimar) {
        res[,iObservacion] <- universalGriddingCV_i(
          iObservacion, coordsObservaciones=coordsObservaciones, 
          fechasObservaciones=fechasObservaciones, valoresObservaciones=valoresObservaciones, 
          params=params, valoresRegresoresSobreObservaciones=valoresRegresoresSobreObservaciones,
          longitudesEnColumnas=longitudesEnColumnas, 
          eliminarSerieTemporalCompleta=eliminarSerieTemporalCompleta, estimarNAs=estimarNAs)
      }
      #res <- sapply(X=iesAEstimar, FUN=universalGriddingCV_i, 
      #              coordsObservaciones=coordsObservaciones, fechasObservaciones=fechasObservaciones, valoresObservaciones=valoresObservaciones, 
      #              pathsRegresores=pathsRegresores, params=params, longitudesEnColumnas=longitudesEnColumnas)
    }
  } else {
    nT <- nrow(valoresObservaciones)
    tUltimoValorSinAR <- min(2 * params$ventanaIgualacionDistribuciones + max(params$tlagsAR) - 1, nT)
    res <- matrix(data=NA, nrow=nrow(valoresObservaciones), ncol=length(iesAEstimar))
    
    # Si tengo un componente auto-regresivo tengo que ir calculando los pasos de tiempo en orden. Se puede obtener un cierto nivel
    # de paralelismo, paralelizando en las observaciones.
    if (nCoresAUsar > 1) {
      cl <- makeCluster(getOption('cl.cores', nCoresAUsar))
      clusterExport(cl, varlist = c('script.dir.interpolarEx'))
      if (exists(x = 'setMKLthreads')) { clusterEvalQ(cl = cl, expr = setMKLthreads(1)) }
      
      res[1:tUltimoValorSinAR, ] <- parSapplyLB(cl=cl, X=iesAEstimar, FUN=universalGriddingCV_i, tIni = 1, tFin=tUltimoValorSinAR, 
                                                coordsObservaciones=coordsObservaciones, fechasObservaciones=fechasObservaciones, 
                                                valoresObservaciones=valoresObservaciones, params=params, 
                                                valoresRegresoresSobreObservaciones=valoresRegresoresSobreObservaciones, 
                                                longitudesEnColumnas=longitudesEnColumnas, eliminarSerieTemporalCompleta=eliminarSerieTemporalCompleta,
                                                estimarNAs=TRUE)
      
      nOrig <- length(valoresRegresoresSobreObservaciones)
      length(valoresRegresoresSobreObservaciones) <- length(valoresRegresoresSobreObservaciones) + length(params$tlagsAR)
      i <- 1
      for (i in 1:length(params$tlagsAR)) {
        valoresRegresoresSobreObservaciones[[nOrig + i]] <- matrix(data=NA, nrow=nrow(valoresObservaciones), ncol=ncol(valoresObservaciones))
        names(valoresRegresoresSobreObservaciones)[nOrig + i] <- paste0('T_', params$tlagsAR[i])
        for (ti in 1:tUltimoValorSinAR + 1) {
          if (ti + params$tlagsAR[i] < nrow(res)) {
            valoresRegresoresSobreObservaciones[[nOrig + i]][ti + params$tlagsAR[i], iesAEstimar] <- res[ti, ]
            valoresRegresoresSobreObservaciones[[nOrig + i]][ti + params$tlagsAR[i], -iesAEstimar] <- valoresObservaciones[ti, -iesAEstimar]
          }
        }        
      }
      
      ti<-tUltimoValorSinAR+1
      for (ti in (tUltimoValorSinAR+1):nT) {
        #print(ti)
        res[ti, ] <- parSapplyLB(cl=cl, X=iesAEstimar, FUN=universalGriddingCV_i, tIni = ti, tFin=ti, coordsObservaciones=coordsObservaciones, 
                                 fechasObservaciones=fechasObservaciones, valoresObservaciones=valoresObservaciones, params=params, 
                                 valoresRegresoresSobreObservaciones=valoresRegresoresSobreObservaciones, longitudesEnColumnas=longitudesEnColumnas, 
                                 eliminarSerieTemporalCompleta=eliminarSerieTemporalCompleta, estimarNAs=TRUE)

        for (i in 1:length(params$tlagsAR)) {
          if (ti + params$tlagsAR[i] < nrow(res)) {
            valoresRegresoresSobreObservaciones[[nOrig + i]][ti + params$tlagsAR[i], iesAEstimar] <- res[ti, ]
            valoresRegresoresSobreObservaciones[[nOrig + i]][ti + params$tlagsAR[i], -iesAEstimar] <- valoresObservaciones[ti, -iesAEstimar]
          }
        }
      }
      
      stopCluster(cl)
    } else {
      i <- 1
      for (i in 1:length(iesAEstimar)) {
        res[1:tUltimoValorSinAR, i] <- universalGriddingCV_i(iObservacion=iesAEstimar[i], tIni = 1, tFin=tUltimoValorSinAR, coordsObservaciones=coordsObservaciones, 
                                                             fechasObservaciones=fechasObservaciones, valoresObservaciones=valoresObservaciones, 
                                                             params=params, valoresRegresoresSobreObservaciones=valoresRegresoresSobreObservaciones,
                                                             longitudesEnColumnas=longitudesEnColumnas, eliminarSerieTemporalCompleta=eliminarSerieTemporalCompleta,
                                                             estimarNAs=TRUE)
      }
      #params$tlagsAR <- c(1, 2)
      #length(valoresRegresoresSobreObservaciones) <- 2
      
      #save(res, file = 'auxiliar/res.txt')
      #res <- attach(what = 'auxiliar/res.txt')$res
      
      # Extiendo la lista de regresores para incluir los lags temporales 
      nOrig <- length(valoresRegresoresSobreObservaciones)
      length(valoresRegresoresSobreObservaciones) <- length(valoresRegresoresSobreObservaciones) + length(params$tlagsAR)
      i <- 1
      for (i in 1:length(params$tlagsAR)) {
        # Cargo los valores de los pasos de inicialización en la matriz de cada lag temporal
        valoresRegresoresSobreObservaciones[[nOrig + i]] <- matrix(data=NA, nrow=nrow(valoresObservaciones), ncol=ncol(valoresObservaciones))
        names(valoresRegresoresSobreObservaciones)[nOrig + i] <- paste0('T_', params$tlagsAR[i])
        for (ti in 1:tUltimoValorSinAR + 1) {
          if (ti + params$tlagsAR[i] < nrow(res)) {
            valoresRegresoresSobreObservaciones[[nOrig + i]][ti + params$tlagsAR[i], iesAEstimar] <- res[ti, ]
            valoresRegresoresSobreObservaciones[[nOrig + i]][ti + params$tlagsAR[i], -iesAEstimar] <- valoresObservaciones[ti, -iesAEstimar]
          }
        }        
      }

      ti<-tUltimoValorSinAR+1
      i <- 1
      for (ti in (tUltimoValorSinAR+1):nT) {
        for (i in 1:length(iesAEstimar)) {
          res[ti, i] <- universalGriddingCV_i(iObservacion=iesAEstimar[i], tIni = ti, tFin=ti, coordsObservaciones=coordsObservaciones, 
                                              fechasObservaciones=fechasObservaciones, valoresObservaciones=valoresObservaciones, 
                                              params=params, valoresRegresoresSobreObservaciones=valoresRegresoresSobreObservaciones,
                                              longitudesEnColumnas=longitudesEnColumnas, eliminarSerieTemporalCompleta=eliminarSerieTemporalCompleta,
                                              estimarNAs=TRUE)
        }
          
        for (i in 1:length(params$tlagsAR)) {
          if (ti + params$tlagsAR[i] < nrow(res)) {
            valoresRegresoresSobreObservaciones[[nOrig + i]][ti + params$tlagsAR[i], iesAEstimar] <- res[ti, ]
            valoresRegresoresSobreObservaciones[[nOrig + i]][ti + params$tlagsAR[i], -iesAEstimar] <- valoresObservaciones[ti, -iesAEstimar]
          }
        }
      }
    }
  }
  
  colnames(res) <- colnames(valoresObservaciones)[iesAEstimar]
  rownames(res) <- rownames(valoresObservaciones)
  return(res)
}

aplicarMascaraRnR <- function(observaciones, interpolacion, params, shpMask) {
  rango <- range(observaciones$value, na.rm = T)
  if (rango[1] < 1E-3 & rango[2] > 1E-3) {
    obsBinarias <- crearObservacionesBinarias(observaciones=observaciones, zcol = 'value')
    valoresObservacionesBinarias <- t(obsBinarias@data)
    binInterpParams <- params
    binInterpParams$modoDiagnostico <- F
    binInterpParams$mLimitarValoresInterpolados <- 'LimitarMinimoyMaximo'
    binInterpParams$minimoLVI <- 0
    binInterpParams$maximoLVI <- 1
    binInterpParams$umbralMascaraCeros <- 0
    binInterpParams$metodoRemocionDeSesgo='ninguno'
    binInterpParams$metodoIgualacionDistribuciones <- 'GLS'
    binInterpParams$ventanaIgualacionDistribuciones <- 1
    binInterpParams$incorporarCoordenadas <- TRUE
    binInterpParams$descartarCoordenadasNoSignificativas <- FALSE
    binInterpParams$interpolationMethod <- 'automap'
    binInterpParams$minRatioRangosParaExtrapolacion <- 0
    if (!is.null(params$valoresRegresoresSobreObservaciones)) {
      binInterpParams$signosValidosRegresores <- rep(
        x = 1, length = length(params$valoresRegresoresSobreObservaciones))
      names(binInterpParams$signosValidosRegresores) <- names(
        params$valoresRegresoresSobreObservaciones)
    }
    #binInterpParams$betaSimpleKriging <- 0
    binInterpShpMask <- interpolacion$shpMask

    if (!is.null(binInterpShpMask)) binInterpShpMask$mask <- rep(
      TRUE, length(interpolacion$predictionLocations))
    interpBinaria <- universalGriddingEx(
      ti=1, fechasObservaciones='', coordsObservaciones=obsBinarias, 
      valoresObservaciones=valoresObservacionesBinarias, 
      valoresRegresoresSobreObservaciones=params$valoresRegresoresSobreObservaciones,
      valoresRegresoresSobreCoordsAInterpolar_ti=params$valoresRegresoresSobreCoordsAInterpolar_ti,
      coordsAInterpolar=interpolacion$predictionLocations, params=binInterpParams, 
      shpMask=binInterpShpMask)
    
    if (params$umbralMascaraCeros == 'auto') {
      binCV <- universalGriddingCV(obsBinarias, 1, matrix(obsBinarias$value, ncol=nrow(obsBinarias)), 
                                   params = binInterpParams)
      
      umbrales <- sort(unique(binCV))
      # umbrales <- sort(unique(over(coordsObservaciones, interpBinaria$predictions)[, 1]))
      umbrales <- (umbrales[1:(length(umbrales)-1)] + umbrales[2:(length(umbrales))]) * 0.5
      # umbralI <- 0.3
      # auc(roc(response = obsBinarias$value, predictor = as.integer(binCV >= umbralI)))
      aucs <- sapply(umbrales, FUN = function(umbralI) {
        roc_obj <- roc(response = obsBinarias$value, predictor = as.integer(binCV >= umbralI))
        #print(plot(roc_obj, print.auc=TRUE, main=umbralI))
        return(auc(roc_obj))
      })
      
      umbral <- umbrales[which.max(aucs)] + 1e-3
    } else {
      umbral <- params$umbralMascaraCeros
    }
    
    if (params$modoDiagnostico) {
      if (gridded(interpolacion$predictionLocations)) {
        mapearGrillaGGPlot(
          grilla = interpBinaria$predictions, shpBase=shpMask$shp, zcol=1, 
          titulo = paste0('Máscara de Ceros. Interpolación Binaria - ', params$strFecha),  
          nomArchResultados = paste0(params$carpetaParaModoDiagnostico, '08.1-MascaraCerosInterpBinaria.png'), 
          dibujar=F, dibujarPuntosObservaciones=T, coordsObservaciones=obsBinarias)
      } else {
        mapearPuntosGGPlot(
          puntos=interpBinaria$predictions, shpBase=shpMask$shp, zcol=1,
          titulo = paste0('Máscara de Ceros. Interpolación Binaria - ', params$strFecha),  
          nomArchResultados = paste0(params$carpetaParaModoDiagnostico, '08.1-MascaraCerosInterpBinaria.png'), 
          dibujar = F, dibujarTexto = T)
      }
    }
    
    # mapearGrillaGGPlot(grilla = interpBinaria$predictions, shpBase, dibujar=F, isolineas=T)
    # interpBinaria$predictions$aux <- as.integer(interpBinaria$predictions@data[,interpBinaria$campoMedia]  >= umbral)
    # mapearGrillaGGPlot(interpBinaria$predictions, shpBase, continuo = F, zcol='aux', dibujar=F)
    #spplot(interpBinaria$predictions, zcol='aux')
    iValidos <- !is.na(interpolacion$predictions@data[,interpolacion$campoMedia]) & !is.na(interpBinaria$predictions@data[,interpBinaria$campoMedia]) 
    #iValidos <- iValidos & interpolacion$predictions@data[iValidos, interpolacion$campoMedia] <= 0.3
    interpolacion$predictions@data[iValidos,interpolacion$campoMedia] <- interpolacion$predictions@data[iValidos,interpolacion$campoMedia] * (interpBinaria$predictions@data[iValidos,interpBinaria$campoMedia] >= umbral) 
  }
  return (interpolacion)
}

simpleBiasAdjustmentEx <- function(observaciones, interpolacion, interpolationParams, zcol=1, gridIndexes=NULL, 
                                   errorRelativoParaCorregir=0.15, inverseDistancePower=NA, shpMask) {
  if (interpolationParams$metodoRemocionDeSesgo != 'ninguno') {
    origModoDiagnostico <- interpolationParams$modoDiagnostico
    interpolationParams$modoDiagnostico <- F
    iValidObs <- !is.na(observaciones@data[,zcol])
    observaciones <- observaciones[iValidObs,]
    
    if (is.null(gridIndexes)) {  
      if (!identicalCRS(observaciones, interpolacion$predictionLocations)) {
        observaciones <- spTransform(observaciones, interpolacion$predictionLocations@proj4string)
      }
      gridIndexes <- getGridIndexes(observaciones=observaciones, grid=interpolacion$predictionLocations)
    } else { gridIndexes <- gridIndexes[iValidObs] }
    
    # Creo el objeto espacial con los valores de los residuos en los puntos de las estaciones
    # iOverObs dice cuales de las observaciones están dentro del área interpolada
    iOverObs <- !is.na(gridIndexes)
  
    # Si el error relativo de alguno de los valores mayores o iguales a 75% del máximo valor observado es 
    # mayor a errorRelativoParaCorregir se aplica la corrección
    maxVal <- max(observaciones$value[iOverObs])
    if (maxVal > 0) {
      spResiduals <- observaciones[iOverObs,]
      spResiduals$value <- observaciones$value[iOverObs] - interpolacion$predictions@data[gridIndexes[iOverObs], zcol]
      if (interpolationParams$metodoRemocionDeSesgo == 'IDW_ResiduosNegativos') { 
        interpolationParams$interpolationMethod <- 'idw'
        interpolationParams$inverseDistancePower <- inverseDistancePower
        spResiduals$value[spResiduals$value > 0] <- 0
      } else if (interpolationParams$metodoRemocionDeSesgo == 'IDW_ResiduosPositivos') { 
        interpolationParams$interpolationMethod <- 'idw'
        interpolationParams$inverseDistancePower <- inverseDistancePower
        spResiduals$value[spResiduals$value < 0] <- 0
      } else if (interpolationParams$metodoRemocionDeSesgo == 'IDW_ResiduosNegativosYPositivos') {
        interpolationParams$interpolationMethod <- 'idw'
        interpolationParams$inverseDistancePower <- inverseDistancePower
      } else {
        stop(paste0('simpleBiasAdjustmentEx: Método de Remoción de Sesgo Desconocido "', interpolationParams$interpolationMethod, '"'))
      }
      #mapearPuntosGGPlot(puntos = spResiduals, shpBase = shpMask$shp, zcol='value', continuo = T)
      
      iAComparar <- which(observaciones$value[iOverObs] >= maxVal * 0.8)
      errorRelativo <- spResiduals$value[iAComparar] / observaciones$value[iOverObs][iAComparar]
      # Solo corrijo si el máximo error relativo es mayor a errorRelativoParaCorregir
      if (max(errorRelativo[!is.infinite(errorRelativo) & !is.nan(errorRelativo)], na.rm=T) >= errorRelativoParaCorregir) {
        interpolationParams$minVal <- min(spResiduals$value, na.rm=T)
        interpolationParams$maxVal <- max(spResiduals$value, na.rm=T)
        interpolationParams$umbralMascaraCeros <- 0
        interpolationParams$metodoRemocionDeSesgo='ninguno'
        if (!is.null(shpMask)) shpMask$mask <- rep(TRUE, length(interpolacion$predictionLocations))
        residualInterpolation <- interpolarEx(observaciones=spResiduals, coordsAInterpolar=interpolacion$predictionLocations, params=interpolationParams, shpMask)
        # mapearGrillaGGPlot(grilla=residualInterpolation$predictions, shpBase = shpBase, continuo=F)
        
        interpolacion$predictions@data[,1] <- interpolacion$predictions@data[,zcol] + residualInterpolation$predictions@data[,zcol]
        if (origModoDiagnostico) {
          if (gridded(interpolacion$predictionLocations)) {
            mapearGrillaGGPlot(
              grilla = residualInterpolation$predictions, shpBase=shpMask$shp, zcol=1, 
              titulo = paste0('Remoción de Sesgo - ', interpolationParams$strFecha),  
              nomArchResultados = paste0(interpolationParams$carpetaParaModoDiagnostico, '07.1-RemocionSesgo.png'), 
              dibujar=F, dibujarPuntosObservaciones=T, coordsObservaciones=spResiduals)
          } else {
            mapearPuntosGGPlot(
              puntos=residualInterpolation$predictions, shpBase=shpMask$shp, zcol=1,
              titulo = paste0('Remoción de Sesgo - ', interpolationParams$strFecha),  
              nomArchResultados = paste0(interpolationParams$carpetaParaModoDiagnostico, '07.1-RemocionSesgo.png'), 
              dibujar = F, dibujarTexto = T)
          }
        }
        #mapearGrillaGGPlot(grilla=interpolacion$predictions, shpBase = shpMask$shp, continuo=F)
        #mapearGrillaGGPlot(grilla=residualInterpolation$predictions, shpBase = shpMask$shp, continuo=T)
      }    
    }
  }
  
  return (interpolacion)
}

getPredictionVector <- function(interpolacion, soloSobreSHP=FALSE) {
  if (soloSobreSHP) {
    return (interpolacion$predictions@data[,interpolacion$campoMedia])
  } else {
    datos <- numeric(length=interpolacion$lengthCoordsAInterpolar)
    datos[!interpolacion$shpMask$mask] <- NA
    datos[interpolacion$shpMask$mask] <- interpolacion$predictions@data[,interpolacion$campoMedia]
    return (datos)
  }
}

getPredictionMatrix <- function(interpolacion) {
  datos <- getPredictionVector(interpolacion, FALSE)
  return (matrix(datos, nrow=interpolacion$nRowCoordsAInterpolar, ncol=interpolacion$nColCoordsAInterpolar, byrow=interpolacion$longitudesEnColumnas))
}

getVarianceVector <- function(interpolacion, soloSobreSHP=FALSE) {
  if (soloSobreSHP) {
    return (interpolacion$predictions@data[,interpolacion$campoVarianza])
  } else {
    datos <- numeric(length=interpolacion$lengthCoordsAInterpolar)
    datos[!interpolacion$shpMask$mask] <- NA
    datos[interpolacion$shpMask$mask] <- interpolacion$predictions@data[,interpolacion$campoVarianza]
    return (datos)
  }
}

getVarianceMatrix <- function(interpolacion) {
  datos <- getVarianceVector(interpolacion, FALSE)
  return (matrix(datos, nrow=interpolacion$nRowCoordsAInterpolar, ncol=interpolacion$nColCoordsAInterpolar, byrow=T))
}

seleccionarMetodoInterpolacion <- function(valoresObservados) {
  # adaptado de la descripción de interpolate de la documentación de intamap
  dataObs <- valoresObservados
  minObs <- min(dataObs)
  if (minObs <= 0)
    dataObs <- dataObs + minObs + sd(dataObs)

  g = boxcox(dataObs ~ 1,lambda=seq(-2.5,2.5,len=101),plotit=FALSE)$y
  tkfn <- fivenum(dataObs)
  rangoIntercuartil <- IQR(dataObs)
  test <- length(boxplot.stats(dataObs)$out)/length(dataObs) > 0.1 ||
    tkfn[3] - tkfn[2] < rangoIntercuartil/3 ||
    tkfn[4] - tkfn[3] < rangoIntercuartil/3 ||
    g[71] < sort(g)[91]
  
  if (test) { return('copula')
  } else { return('automap') }
}

salvarInterpolacion <- function(baseNomArchResultados, interpolacion, formatoSalida = 'binary', salvarPrediccion=T, salvarVarianza=F, 
                                NAValue=-.Machine$double.xmax) {
  if (formatoSalida == 'binary') {
    if (salvarPrediccion) {
      nomArchResultados <- changeFileExt(appendToFileName(baseNomArchResultados, '_Pred'), '.bin')
      if (gridded(interpolacion$predictionLocations)) {
        datos <- getPredictionMatrix(interpolacion=interpolacion)
        salvarMatrizABinarioSinDimensiones(nomArchResultados, matriz=datos, NAValue=NAValue)
      } else { 
        datos <- getPredictionVector(interpolacion=interpolacion) 
        salvarVectorABinarioSinDimensiones(nomArchResultados, vector=datos, NAValue=NAValue)
      }
    }
    
    if (salvarVarianza && !is.null(interpolacion$campoVarianza)) {
      nomArchResultados <- changeFileExt(appendToFileName(baseNomArchResultados, '_Var'), '.bin')
      if (gridded(interpolacion$predictionLocations)) { 
        datos <- getVarianceMatrix(interpolacion=interpolacion)
        salvarMatrizABinarioSinDimensiones(nomArchResultados, matriz=datos, NAValue=NAValue)
      } else { 
        datos <- getVarianceVector(interpolacion=interpolacion)
        salvarVectorABinarioSinDimensiones(nomArchResultados, vector=datos, NAValue=NAValue)
      }
    }
  } else if (formatoSalida == 'netCDF') {
    # TO-DO
    stop('salvarInterpolacion: formatoSalida "netCDF" no implementado')
  } else if (formatoSalida == 'GeoTiff') {
    #if (i==9) { #para guardar el IBH en raster
    source('../grillas/uIOGrillas.r', encoding = 'WINDOWS-1252')
    if (salvarPrediccion) {
      guardarGrillaGDAL(changeFileExt(baseNomArchResultados, '.tif'), interpolacion$predictions)
    }    
  }
}

expandirPathsRegresor <- function(
    fechasRegresor=funcFechasRegresor(pathsRegresor, tz(fechasRegresando[1])), pathsRegresor, 
    fechasRegresando, expandir=T, funcFechasRegresor=findDateInText) {
  # retorna los paths de un regresor, repitiendo sus valores si tiene mayor periodicidad que el regresando.
  # Por ejemplo si el regresor es cada 8 días y el regresando cada 1, se repite el valor del mismo regresor
  # para los 8 días que abarca su período
  # TO-DO: actualmente retorna en cada ti el último valor previo a ti, hacer una opción para hacerlo centrado
  if (expandir) {
    nuevoPathsRegresor <- character(length(fechasRegresando))
    for (ti in 1:length(fechasRegresando)) {
      fechasPrevias <- which(fechasRegresor <= fechasRegresando[ti])
      
      if (length(fechasPrevias) > 0) {
        iFechaTiEnRegresores <- max(fechasPrevias)
        nuevoPathsRegresor[ti] <- pathsRegresor[iFechaTiEnRegresores]
      } else {
        nuevoPathsRegresor[ti] <- NA
      }
    }
  } else {
    iFechas <- match(fechasRegresando, fechasRegresor)
    nuevoPathsRegresor <- pathsRegresor[iFechas]
  }

  return(nuevoPathsRegresor)
}

nDatosNoNulosEnRegresores_ti <- function (ti=1, pathsRegresores, shpMask=NULL) {
  require('rgdal')
  
  if (is.null(shpMask)) { iPixeles <- T
  } else { iPixeles <- shpMask$mask }
  
  i <- 1
  while (i <= ncol(pathsRegresores) && length(iPixeles) > 0) {
    if (!is.na(pathsRegresores[ti, i])) {
      regresor <- try(readGDAL(pathsRegresores[ti, i], silent=T))
      if ('try-error' %in% class(regresor)) { regresor <- NULL }
    } else { regresor <- NULL }
    
    if (is.null(regresor) || is.na(regresor)) { iPixeles <- integer(0)
    } else { iPixeles <- iPixeles & !is.na(regresor@data[, 1]) }
    i <- i + 1
  }
  
  return(sum(iPixeles))
}

nDatosNoNulosEnRegresores <- function (pathsRegresores, shpMask=NULL) {
  nCoresAUsar <- detectCores(T, T)
  if (nCoresAUsar > 1) {
    cl <- makeCluster(getOption('cl.cores', nCoresAUsar))
    if (exists(x = 'setMKLthreads')) { clusterEvalQ(cl = cl, expr = setMKLthreads(1)) }
    nDatosNoNulos <- parSapplyLB(cl=cl, X=1:nrow(pathsRegresores), FUN=nDatosNoNulosEnRegresores_ti, pathsRegresores=pathsRegresores, shpMask=shpMask)
    stopCluster(cl)
  } else {
    nDatosNoNulos <- sapply(X=1:nrow(pathsRegresores), FUN=nDatosNoNulosEnRegresores_ti, pathsRegresores=pathsRegresores, shpMask=shpMask)
  }
  
  return(nDatosNoNulos)
}

rellenarSP <- function(sp, mascara=rep(TRUE, length(sp)), metodo='automap', nMuestras=3 * nCuadrantesX * nCuadrantesY * nCuadrantesZ, 
                       nRepeticiones = 10, zcol=1, distanciaMaximaMuestras=0, 
                       nCuadrantesX=4, nCuadrantesY=round(nCuadrantesX * diff(bbox(sp)[2,]) / diff(bbox(sp)[1,])), nCuadrantesZ=5, 
                       iEsObligatoriosEnLaMuestra=NULL, pathsRegresores=NULL, params=NULL, minNDatosDisponibles=30, 
                       minNCuadrantesDisponibles=round(nCuadrantesX*nCuadrantesY*0.25),
                       sustituirConPrimerRegresorSiNoHayDatosDisponibles=TRUE, 
                       maxNMuestrasPorCuadrante = round(1.1 * nMuestras / (nCuadrantesX * nCuadrantesY))) {
  require('fields')
  require('stats')

  iAInterpolar <- is.na(sp@data[, zcol]) & mascara
  if (any(iAInterpolar)) {
    nDatosDisponiblesPorCuadrante <- contarNoNulosPorCuadrantesEnObjSP(objSP = sp, nCuadrantesX = nCuadrantesX, nCuadrantesY = nCuadrantesY, zcol = zcol,  shpMask = NULL)

    #mapearGrillaGGPlot(as(sp, 'SpatialPixelsDataFrame'), shpBase = shpMask$shp, continuo = T, dibujar = F)
    if (!is.null(params) & params$difMaxFiltradoDeOutliersRLM > 0 &
        sum(nDatosDisponiblesPorCuadrante >= minNDatosDisponibles) >= minNCuadrantesDisponibles) {
      listaMapasAux <- createDefaultListaMapas(paramsIyM = params, fechasObservaciones = NULL, nObservacionesTemporales = 1)
      outliersRLM <- deteccionOutliersRLM(coordsObservaciones = sp, fechasObservaciones = NULL, valoresObservaciones = matrix(sp@data[,zcol], nrow=1),
                                          params=params, pathsRegresores=pathsRegresores, listaMapas=listaMapasAux, 
                                          factorMADHaciaAbajo=params$difMaxFiltradoDeOutliersRLM, factorMADHaciaArriba = 8)
      
      if (nrow(outliersRLM) > 0) {
        sp@data[outliersRLM$iOutlier, zcol] <- NA
        iAInterpolar <- is.na(sp@data[, zcol]) & mascara
        nDatosDisponiblesPorCuadrante <- contarNoNulosPorCuadrantesEnObjSP(objSP = sp, nCuadrantesX = nCuadrantesX, nCuadrantesY = nCuadrantesY, zcol = zcol,  shpMask = NULL)
      }
    }
    
    if (sum(nDatosDisponiblesPorCuadrante >= minNDatosDisponibles) >= minNCuadrantesDisponibles) {
      # Si hay algún valor nulo en el área a interpolar
      iMuestras <- !is.na(sp@data[, zcol])
      #mapearGrillaGGPlot(as(sp, 'SpatialPixelsDataFrame'), shpBase = shpMask$shp, continuo = T, dibujar = F)
      nDatosDisponibles <- sum(nDatosDisponiblesPorCuadrante)  
            
      # Si hay al menos minNDatosDisponibles valores no nulos en toda el área disponible
      if (distanciaMaximaMuestras > 0) {
        aux1 <- sp::coordinates(sp)[iMuestras,]
        aux2 <- sp::coordinates(sp)[iAInterpolar,]
        
        distancia <- numeric(nrow(aux1))  
        nBloque <- trunc(2^31 / (nrow(aux2) * 8))
        nMax <- trunc(nrow(aux1) / nBloque)
        
        for (i in 1:nMax) {
          ies <- ((i-1) * nBloque + 1):((i) * nBloque)
          distancia[ies] <- rowMins(rdist(aux1[ies, ], aux2))
        }
        ies <- (nMax * nBloque + 1):length(distancia)
        distancia[ies] <- rowMins(rdist(aux1[ies, ], aux2))
        
        iMuestras[iMuestras] <- distancia <= distanciaMaximaMuestras
        rm(aux1, aux2, nBloque, nMax, ies, distancia)
      }
      
      if (!is.null(iEsObligatoriosEnLaMuestra)) {
        iEsObligatoriosEnLaMuestra <- iEsObligatoriosEnLaMuestra[!is.na(sp@data[iEsObligatoriosEnLaMuestra, zcol])]
        nMuestras <- nMuestras - length(iEsObligatoriosEnLaMuestra)
      }
        
      iMetodo <- which(metodo == c('Tps', 'fastTps', 'loess', 'idw', 'automap'))

      # Si tengo menos valores disponibles que la cantidad de muestras por la cantidad de repeticiones que quiero usar. 
      # Hago todas las repeticiones que me de para hacer y reparto las muestras en partes iguales entre ellas
      if (nDatosDisponibles < nMuestras * nRepeticiones) {
        nMuestras <- max(minNDatosDisponibles, trunc(nDatosDisponibles / nRepeticiones))
        nRepeticiones <- ceiling(nDatosDisponibles / nMuestras)
      } 
      
      invNRepeticiones <- 1 / nRepeticiones
      
      aux <- rep(0, sum(iAInterpolar))
      
      iEsMuestreados <- integer(0)
      #i <- 1
      #i <- i+1
      for (i in 1:nRepeticiones) {
        # Hago el muestreo
        ies <- muestrearEnCuadrantesYECDF(sp = sp, size=nMuestras, nCuadrantesX = nCuadrantesX, nCuadrantesY = nCuadrantesY, 
                                          nCuadrantesZ = nCuadrantesZ, zcol=zcol, iEsAExcluir = iEsMuestreados, 
                                          maxNMuestrasPorCuadrante=maxNMuestrasPorCuadrante)
        iEsMuestreados <- c(iEsMuestreados, ies)
        if (!is.null(iEsObligatoriosEnLaMuestra)) ies <- unique(c(ies, iEsObligatoriosEnLaMuestra))
        
        # Ejecuto el método y acumulo el promedio en aux
        if (iMetodo == 1) {
          # Tps
          suppressWarnings(tps <- Tps(x = sp::coordinates(sp)[ies, ], Y = sp@data[ies, zcol], lon.lat = !is.projected(sp), miles=FALSE))
          aux <- aux + predict(tps, x=sp::coordinates(sp)[iAInterpolar,]) * invNRepeticiones
        } else if (iMetodo == 2) {
          # fastTps
          suppressWarnings(ftps <- fastTps(x=sp::coordinates(sp)[ies,], Y = sp@data[iMuestras, zcol], lon.lat = !is.projected(sp), theta = 10000, lambda = 0.2))
          aux <- aux + predict.fastTps(object = ftps, xnew = sp::coordinates(sp)[iAInterpolar,]) * invNRepeticiones
        } else if (iMetodo == 3) {
          # loess
          data <- data.frame(x=sp::coordinates(sp)[iMuestras, 1], y=sp::coordinates(sp)[ies, 2], z=sp@data[iMuestras, zcol])
          data.loess = loess(z ~ x * y, data=data)
          data.fit = data.frame(x=sp::coordinates(sp)[iAInterpolar, 1], y=sp::coordinates(sp)[iAInterpolar, 2])
          aux <- aux + predict(data.loess, newdata =data.fit) * invNRepeticiones
        } else if (iMetodo %in% c(4, 5)) {
          # idw, automap
          # Preparo los objetos espaciales de esta forma así su geometría no cambia y no es necesario
          # volver a cachear los regresores estáticos
          sp2 <- sp
          sp2@data[-ies, zcol] <- NA
          
          shpMask2 <- list(shp=NULL, mask=mascara & iAInterpolar)

          if (is.null(params)) { params <- createParamsInterpolarYMapear(coordsAInterpolarSonGrilla = gridded(interpLocations), interpolationMethod = metodo, inverseDistancePower = 2)
          } else { params$interpolationMethod <- metodo }
          params$rellenarRegresores <- F
          params$nCoresAUsar <- 1
          params$pathSHPMapaBase <- ''
          
          if (!is.null(pathsRegresores)) {
            valoresRegresoresSobreObservaciones <- extraerValoresRegresoresSobreSP(objSP = sp2, pathsRegresores = pathsRegresores, nCoresAUsar = params$nCoresAUsar)
          } else { valoresRegresoresSobreObservaciones <- NULL }
          valsObs <- matrix(sp2@data[,zcol], nrow=1)

          #ti = 1
          #coordsObservaciones = sp2
          #fechasObservaciones = NULL
          #valoresObservaciones = valsObs
          #coordsAInterpolar = geometry(sp)
          #shpMask = shpMask2
          #paramsParaRellenoRegresores = NULL
          #pathsRegresoresParaRellenoRegresores = NULL
          #mapearPuntosGGPlot(puntos = SpatialPointsDataFrame(coords = sp2, data = data.frame(value=valsObs[1,])), shpBase = shpMask$shp, dibujar = F, continuo = T)

          #if (params$difMaxFiltradoDeOutliersRLM > 0) {
          #  listaMapasAux <- createDefaultListaMapas(paramsIyM = params, fechasObservaciones = NULL, nObservacionesTemporales = 1)
          #  outliersRLM <- deteccionOutliersRLM(coordsObservaciones = sp2, fechasObservaciones = NULL, valoresObservaciones = valsObs,
          #                                      params=params, pathsRegresores=pathsRegresores, listaMapas=listaMapasAux, 
          #                                      factorMADHaciaAbajo=params$difMaxFiltradoDeOutliersRLM, factorMADHaciaArriba = 8)
          #  if (nrow(outliersRLM) > 0) valsObs[outliersRLM$iOutlier] <- NA
          #}
          
          regresorRellenado <- universalGridding(ti = 1, coordsObservaciones = sp2, fechasObservaciones = NULL, pathsRegresores = pathsRegresores, 
                                                 valoresObservaciones = valsObs, coordsAInterpolar = geometry(sp), shpMask = shpMask2,
                                                 params = params, valoresRegresoresSobreObservaciones = valoresRegresoresSobreObservaciones,
                                                 paramsParaRellenoRegresores = NULL, pathsRegresoresParaRellenoRegresores = NULL)
          
          if (F) {
            # Método viejo, recacheaba los regresores estáticos todos los pasos pero un poquito más rápido si no hay regresores estáticos
            obs <- sp[ies,]
            names(obs)[zcol] <- 'value'
            interpLocations <- sp[iAInterpolar, ]
            # mapearPuntosGGPlot(obs, shpBase = shpMask$shp, continuo = T, tamaniosPuntos = 1, dibujar=F)
            
            if (is.null(params)) { params <- createParamsInterpolarYMapear(coordsAInterpolarSonGrilla = gridded(interpLocations), interpolationMethod = metodo, inverseDistancePower = 2)
            } else { params$interpolationMethod <- metodo }
            params$rellenarRegresores <- F
            params$nCoresAUsar <- 1
            params$pathSHPMapaBase <- ''
            
            if (!is.null(pathsRegresores)) {
              valoresRegresoresSobreObservaciones <- extraerValoresRegresoresSobreSP(objSP = obs, pathsRegresores = pathsRegresores, nCoresAUsar = params$nCoresAUsar)
            } else { valoresRegresoresSobreObservaciones <- NULL }
            #ti = 1
            #coordsObservaciones = obs
            #fechasObservaciones = NULL
            #valoresObservaciones = matrix(obs@data[,zcol], nrow=1)
            #coordsAInterpolar = interpLocations 
            #paramsParaRellenoRegresores = NULL
            #pathsRegresoresParaRellenoRegresores = NULL
            regresorRellenado <- universalGridding(ti = 1, coordsObservaciones = obs, fechasObservaciones = NULL, pathsRegresores = pathsRegresores, 
                                                   valoresObservaciones = matrix(obs@data[,zcol], nrow=1), coordsAInterpolar = interpLocations, 
                                                   params = params, valoresRegresoresSobreObservaciones = valoresRegresoresSobreObservaciones,
                                                   paramsParaRellenoRegresores = NULL, pathsRegresoresParaRellenoRegresores = NULL)
          }
          
          aux <- aux + regresorRellenado$predictions@data[, regresorRellenado$campoMedia] * invNRepeticiones
        }
      }

      # mapearGrillaGGPlot(grilla = sp, shpBase = shpMask$shp, continuo=T, titulo = 'Original', dibujar=F)
      # sp3 <- SpatialPixelsDataFrame(points = sp, data = data.frame(value=sp@data[,zcol]))
      # sp3@data[iAInterpolar, zcol] <- aux / invNRepeticiones
      # mapearGrillaGGPlot(sp3, shpBase = shpMask$shp, continuo=T, titulo = 'Rellenado', dibujar=F)
      
      sp@data[iAInterpolar, zcol] <- aux
    } else if (sustituirConPrimerRegresorSiNoHayDatosDisponibles && ncol(pathsRegresores) >= 1) {
      sp@data[, zcol] <- extraerValorRegresorSobreSP(i = 1, objSP = geometry(sp), pathsRegresor = pathsRegresores, zcol = zcol, silent = T)
      # sp2 <- SpatialPixelsDataFrame(points = sp, data = sp@data)
      # mapearGrillaGGPlot(sp2, shpBase = shpMask$shp, continuo=T, titulo = 'Rellenado')
      # mapearGrillaGGPlot(readGDAL(changeFileExt(listaMapas$nombreArchivo[iTi], nuevaExtensionConPunto = '.tif')), shpBase = shpMask$shp, continuo=T, titulo = 'Rellenado')
    }
  }
  
  return(sp)
}

getIEsACombinarReduccionSeries <- function(coordsObservaciones, radioReduccionSeriesKm=1) {
  radioReduccionSeries <- distKmToP4Str(p4str = proj4string(coordsObservaciones), distKm = radioReduccionSeriesKm)
  
  iesACombinar <- list()
  
  if (radioReduccionSeries > 0) {
    dists <- spDists(x = coordsObservaciones)
    i <- 1
    for (i in 1:(nrow(dists)-1)) {
      iesI <- numeric()
      for (j in (i+1):ncol(dists)) {
        if (dists[i, j] < radioReduccionSeries) {
          iesI <- c(iesI, j)
        }
      }
      if (length(iesI) > 0) {
        iesI <- c(iesI, i)
        iesACombinar[[length(iesACombinar) + 1]] <- iesI
      }
    }
  }
  return(iesACombinar)
}

reducirSpatialPointsDataFrame <- function(coordsObservaciones, radioReduccionSeriesKm=1, funcionReduccionSeries='mean', iesACombinar=NULL, 
                                          zcolsValores = sapply(X = 1:ncol(coordsObservaciones@data), FUN = function(x, data) is.numeric(data[,x]), data=coordsObservaciones@data)) {
  if (is.null(iesACombinar)) { iesACombinar <- getIEsACombinarReduccionSeries(coordsObservaciones = coordsObservaciones, radioReduccionSeriesKm = radioReduccionSeriesKm) }
  if (length(iesACombinar) > 0) {
    rowNames <- row.names(coordsObservaciones)
    coords <- sp::coordinates(coordsObservaciones)
    
    rowNamesReducidos <- sapply(iesACombinar, FUN = function(x) {return(paste0(row.names(coordsObservaciones)[x], sep=',', collapse = ','))})
    coordsReducidas <- t(sapply(iesACombinar, FUN = function(x) { apply(coords[x, , drop=F], MARGIN = 2, FUN = naSiTodosNAFuncSiNo,
                                                                        func=mean, na.rm=T)} ))
    vals <- as.data.frame(coordsObservaciones@data)
    valsReducidos <- vals[sapply(iesACombinar, FUN = function(x) x[1]), , drop=F]
                     
    if (!is.null(zcolsValores) & any(zcolsValores)) {
      valsZCols <- as.data.frame(coordsObservaciones@data[,zcolsValores, drop=F])
      valsReducidosZCols <- as.data.frame(sapply(iesACombinar, FUN = function(x) { apply(valsZCols[x, , drop=F], MARGIN = 2, FUN = naSiTodosNAFuncSiNo,
                                                                                    func=get(funcionReduccionSeries), na.rm=T)} ))
      valsReducidosZCols <- as.numeric(apply(valsReducidosZCols, MARGIN = 2, FUN = function(x) { 
        x[is.nan(x)] <- NA
        return(x) }
      ))
      
      valsReducidos[, zcolsValores] <- valsReducidosZCols
    } 
    
    iesAEliminar <- integer()
    i <- 1
    for (i in seq_along(iesACombinar)) {
      iCombinado <- iesACombinar[[i]][length(iesACombinar[[i]])]
      
      rowNames[iCombinado] <- rowNamesReducidos[i]
      coords[iCombinado, ] <- coordsReducidas[i, ]
      vals[iCombinado,] <- valsReducidos[i, ]
      
      iesAEliminar <- c(iesAEliminar, iesACombinar[[i]][1:(length(iesACombinar[[i]])-1)])
    }
    rowNames <- rowNames[-iesAEliminar]
    coords <- coords[-iesAEliminar, ]
    vals <- vals[-iesAEliminar, ,drop=F]
    
    res <- SpatialPointsDataFrame(coords = coords, data = vals, proj4string = coordsObservaciones@proj4string)
  } else {
    res <- coordsObservaciones
  }
  return(res)
}

reducirSpatialPointsDataFrameYMatrizObservaciones <- function(coordsObservaciones, valoresObservaciones, radioReduccionSeriesKm=1, funcionReduccionSeries='mean', iesACombinar=NULL) {
  if (is.null(iesACombinar)) { iesACombinar <- getIEsACombinarReduccionSeries(coordsObservaciones = coordsObservaciones, radioReduccionSeriesKm = radioReduccionSeriesKm) }
  coordsObservaciones <- reducirSpatialPointsDataFrame(coordsObservaciones = coordsObservaciones, radioReduccionSeriesKm = radioReduccionSeriesKm, 
                                                       funcionReduccionSeries = funcionReduccionSeries, iesACombinar = iesACombinar)  
  if (length(iesACombinar) > 0) {
    #i <- 2
    #for (i in (1:length(iesACombinar))) {
    #  print(i)
    #  mean(valoresObservaciones[, iesACombinar[[i]]], na.rm=T)
    #  print(valoresObservaciones[, iesACombinar[[i]]])
    #  apply(valoresObservaciones[, x], MARGIN = 1, FUN = naSiTodosNAFuncSiNo, func=get(funcionReduccionSeries), na.rm=T)
    #}
    # x <- iesACombinar[[1]]
    datosReducidos <- matrix(data=sapply(iesACombinar, 
                                         FUN = function(x) { apply(valoresObservaciones[, x, drop=F], MARGIN = 1, FUN = naSiTodosNAFuncSiNo,
                                                             func=get(funcionReduccionSeries), na.rm=T) }), 
                             nrow = nrow(valoresObservaciones), ncol = length(iesACombinar))
    datosReducidos[is.nan(datosReducidos)] <- NA
    
    iesAEliminar <- integer()
    for (i in 1:length(iesACombinar)) {
      valoresObservaciones[,iesACombinar[[i]][length(iesACombinar[[i]])]] <- datosReducidos[, i]
      iesAEliminar <- c(iesAEliminar, iesACombinar[[i]][1:(length(iesACombinar[[i]])-1)])
    }

    valoresObservaciones <- valoresObservaciones[, -iesAEliminar, drop=F]
  }
  return(list(coordsObservaciones=coordsObservaciones, valoresObservaciones=valoresObservaciones))
}

reducirSpatialPointsDataFrameYSeries <- function(coordsObservaciones, series, radioReduccionSeriesKm=1, funcionReduccionSeries='mean', iesACombinar=NULL) {
  if (is.null(iesACombinar)) {iesACombinar <- getIEsACombinarReduccionSeries(coordsObservaciones = coordsObservaciones, radioReduccionSeriesKm = radioReduccionSeriesKm)}
  coordsObservaciones <- reducirSpatialPointsDataFrame(coordsObservaciones = coordsObservaciones, radioReduccionSeriesKm = radioReduccionSeriesKm, 
                                                       funcionReduccionSeries = funcionReduccionSeries, iesACombinar = iesACombinar)  
  if (length(iesACombinar) > 0) {
    # lapply(iesACombinar, function(x) {series$estaciones$Observacion[x]})
    
    #i <- 2
    #for (i in (1:length(iesACombinar))) {
    #  print(i)
    #  mean(series$datos[, iesACombinar[[i]]], na.rm=T)
    #  print(series$datos[, iesACombinar[[i]]])
    #  apply(series$datos[, x], MARGIN = 1, FUN = naSiTodosNAFuncSiNo, func=get(funcionReduccionSeries), na.rm=T)
    #}
    datosReducidos <- sapply(iesACombinar, FUN = function(x) { apply(series$datos[, x, drop=F], MARGIN = 1, FUN = naSiTodosNAFuncSiNo,
                                                                     func=get(funcionReduccionSeries), na.rm=T)} )
    datosReducidos[is.nan(datosReducidos)] <- NA
    
    iesAEliminar <- integer()
    for (i in 1:length(iesACombinar)) {
      series$datos[,iesACombinar[[i]][length(iesACombinar[[i]])]] <- datosReducidos[, i]
      iesAEliminar <- c(iesAEliminar, iesACombinar[[i]][1:(length(iesACombinar[[i]])-1)])
    }
    
    series$estaciones <- series$estaciones[-iesAEliminar]
    series$datos <- series$datos[, -iesAEliminar]
  }
  return(list(coordsObservaciones=coordsObservaciones, series=series))
}

calcOutlyingnessMediaSD <- function(x) {
  medias <- rowMeans(x = x, na.rm = T)
  desvs <- rowSds(x = x, na.rm = T)
  return((x - medias) / desvs)
}

calcOutlyingnessMedianaMAD <- function(x, xEstimarMedianaYMAD=x, porFilas=T, overlap=0) {
  if (porFilas) {
    if (overlap > 0) {
      medianas <- numeric(nrow(x))
      desvs <- numeric(nrow(x))
      for (i in 1:nrow(x)) {
        ies <- (i-overlap):(i+overlap)
        ies <- ies[ies >= 1 & ies <= nrow(x)]
        medianas[i] <- median(xEstimarMedianaYMAD[ies,], na.rm = T)
        desvs[i] <- mad(xEstimarMedianaYMAD[ies,], center = medianas[i], na.rm = T)
      }
    } else {
      medianas <- rowMedians(xEstimarMedianaYMAD, na.rm = T)
      desvs <- rowMads(xEstimarMedianaYMAD, center = medianas, na.rm = T)
    }
    return((x - medianas) / desvs)
  } else {
    mediana <- median(xEstimarMedianaYMAD, na.rm = T)
    desv <- mad(xEstimarMedianaYMAD, center = mediana, na.rm = T)
    return((x - mediana) / desv)
  }
}

deteccionOutliersMediaSD <- function(
    x, factorSDHaciaAbajo=3.5, factorSDHaciaArriba=factorSDHaciaAbajo, sdMin=NA) {
  estimados <- as.numeric(rowMeans(x, na.rm = T))
  stdDifs <- as.numeric(t(apply(x, MARGIN = 1, FUN = outlyingnessMediaSD, sdMin=sdMin)))
  
  return(createDFTestsConEstimadosYStdDifs(
    x, estimados=estimados, stdDifs=stdDifs, factorHaciaAbajo=factorSDHaciaAbajo, 
    factorHaciaArriba=factorSDHaciaArriba))
}

deteccionOutliersMedianaMAD <- function(
    x, factorMADHaciaAbajo=3.5, factorMADHaciaArriba=factorMADHaciaAbajo, 
    desvMedAbsMin=NA) {
  
  estimados <- as.numeric(rowMedians(x, na.rm = T, keep.names = F))
  stdDifs <- as.numeric(t(apply(
    x, MARGIN = 1, FUN = outlyingnessMedianaMAD, desvMedAbsMin=desvMedAbsMin)))
  
  return(createDFTestsConEstimadosYStdDifs(
    x, estimados=estimados, stdDifs=stdDifs, factorHaciaAbajo=factorMADHaciaAbajo, 
    factorHaciaArriba=factorMADHaciaArriba))
}

deteccionOutliersRLM <- function(
    coordsObservaciones, fechasObservaciones, valoresObservaciones, params, pathsRegresores, 
    listaMapas, 
    factorMADHaciaAbajo=3.5, factorMADHaciaArriba=factorMADHaciaAbajo, desvMedAbsMin=NA,
    factorSDHaciaAbajo=NA, factorSDHaciaArriba=factorSDHaciaAbajo, sdMin=NA, 
    returnTestDF=FALSE) {
  if (!is.na(factorMADHaciaAbajo) && !is.na(factorSDHaciaAbajo)) {
    stop('interpolarEx.deteccionOutliersRLM: only one of factorMADHaciaAbajo and factorSDHaciaAbajo can be specified')
  }
  
  paramsAux <- params
  paramsAux$interpolationMethod = 'none'
  paramsAux$metodoIgualacionDistribuciones <- 'regresionLinealRobusta'
  paramsAux$difMaxFiltradoDeOutliersRLM <- 0
  paramsAux$difMaxFiltradoDeOutliersCV <- 0
  paramsAux$radioReduccionSeriesKm <- 0
  paramsAux$pathSHPMapaBase <- ''
  paramsAux$tlagsAR <- NULL
  paramsAux$modoDiagnostico <- FALSE
  # paramsAux$nCoresAUsar <- 1
  listaMapasAux <- listaMapas
  listaMapasAux$nombreArchivo <- agregarCarpetaAlFinal(listaMapasAux$nombreArchivo, 'RLM')
  listaMapasAux$dibujarObservacionesEscalaFija <- F
  listaMapasAux$dibujarEscalaFija <- F
  listaMapasAux$dibujarObservacionesEscalaAdaptada <- F
  listaMapasAux$dibujarEscalaAdaptada <- F
  listaMapasAux$generarThumbnailFija <- F
  listaMapasAux$generarThumbnailAdaptada <- F
  listaMapasAux$incluirIsolineaFija <- F
  listaMapasAux$incluirIsolineaAdaptada <- F
  listaMapasAux$dibujarPuntosObservacionesAdaptada <- F
  listaMapasAux$salvarGeoTiff <- F
  listaMapasAux$salvarBin <- F
  listaMapasAux$titulo <- ''
  listaMapasAux$incluirSubtitulo <- F
  listaMapasAux$recalcularSiYaExiste <- F
  
  if (length(coordsObservaciones) <= 100) {
    # Si tengo pocas observaciones hago CV, sino hago una única regresión
    # params = paramsAux
    pred <- universalGriddingCV(
      coordsObservaciones = coordsObservaciones, fechasObservaciones = fechasObservaciones, 
      valoresObservaciones = valoresObservaciones, params = paramsAux, 
      pathsRegresores=pathsRegresores)
  } else {
    #coordsAInterpolar = coordsObservaciones
    #paramsIyM = paramsAux
    #shpMask = NULL
    #listaMapas = listaMapasAux
    #paramsParaRellenoRegresores = NULL 
    #pathsRegresoresParaRellenoRegresores = NULL
    #espEscalaFija = NULL
    #espEscalaAdaptada = NULL
    #returnInterpolacion = T
    #mapearGrillaGGPlot(SpatialPixelsDataFrame(points = coordsObservaciones, data = data.frame(value=valoresObservaciones[1,])), shpBase = shpMask$shp, continuo = T, dibujar = F)
    interp <- interpolarYMapear(coordsObservaciones = coordsObservaciones, fechasObservaciones = fechasObservaciones,
                                valoresObservaciones =  valoresObservaciones, pathsRegresores = pathsRegresores, 
                                coordsAInterpolar = coordsObservaciones, xyLims = xyLims, paramsIyM = paramsAux, 
                                shpMask = NULL, listaMapas = listaMapasAux, paramsParaRellenoRegresores = NULL, 
                                pathsRegresoresParaRellenoRegresores = NULL, espEscalaFija = NULL, espEscalaAdaptada = NULL,
                                returnInterpolacion = T)
    #mapearGrillaGGPlot(grilla = SpatialPixelsDataFrame(points = coordsObservaciones, data = interp[[1]]$predictions@data), shpBase = shpMask$shp, continuo = T, dibujar = F)
    pred <- do.call(rbind, lapply(interp, FUN = function(x) { x$predictions@data[,x$campoMedia] }))
    rm(interp)
  }
  residuos <- valoresObservaciones - pred
  
  # spResiduos <- SpatialPixelsDataFrame(points = coordsObservaciones, data=data.frame(value=residuos[1,]))
  # spResiduos@data[!is.na(spResiduos@data[,1]) & abs(spResiduos@data[,1]) < 3.5, 1] <- NA
  # mapearGrillaGGPlot(spResiduos, shpBase = shpMask$shp, continuo = T, dibujar = F)
  # plot(density(residuos, na.rm=T))
  # sp@data[,'pred'] <- as.vector(pred)
  # sp@data[,'res'] <- as.vector(residuos)
  # sp@data[,'out'] <- as.vector(outlyingness)
  # sp@data[,'lala'] <- as.numeric((is.na(outlyingness) | abs(outlyingness) < 3))
  
  #mapearGrillaGGPlot(lala, shpBase = shpMask$shp, continuo = T, dibujar = F)
  #mapearGrillaGGPlot(lala, shpBase = shpMask$shp, continuo = T, dibujar = F, zcol=2)
  #mapearGrillaGGPlot(lala, shpBase = shpMask$shp, continuo = T, dibujar = F, zcol=3)
  #mapearGrillaGGPlot(lala, shpBase = shpMask$shp, continuo = T, dibujar = F, zcol=4)
  #lala[!is.na(outlyingness) & abs(outlyingness) > 2.8, ] <- NA
  
  #lala <- SpatialPixelsDataFrame(points = sp, data = sp@data)
  #plot(density(outlyingness, na.rm=T))
  #sum(outlyingness < -2.5, na.rm=T)
  
  #factorMADHaciaAbajo <- 2.5
  #factorMADHaciaArriba <- factorMADHaciaAbajo
  
  #plot(apply(valoresObservaciones, 1, FUN = function(x) { sum(!is.na(x)) }))
  #which.min(apply(valoresObservaciones, 1, FUN = function(x) { sum(!is.na(x)) }))
  #lctype <- raster(x = 'D:/LCType_PuntaDelEsteCorregido/LCType_PuntaDelEsteCorregido.tif')
  #lctypeAux <- extract(lctype, coordsObservaciones)
  #cbind(estaciones$Nombre, lctypeAux)
  
  
  if (!is.na(factorMADHaciaAbajo)) {
    # outlyingnessMedianaMAD(x = residuos['2019-01-08',])
    outlyingness <- t(apply(residuos, MARGIN = 1, FUN = outlyingnessMedianaMAD, 
                            desvMedAbsMin=desvMedAbsMin))
    factorHaciaAbajo <- factorMADHaciaAbajo
    factorHaciaArriba <- factorMADHaciaArriba
  } else {
    outlyingness <- t(apply(residuos, MARGIN = 1, FUN = outlyingnessMediaSD, sdMin=sdMin))
    factorHaciaAbajo <- factorSDHaciaAbajo
    factorHaciaArriba <- factorSDHaciaArriba
  }
  
  if (!returnTestDF) {
    iEsAFiltrar <- outlyingness < -factorHaciaAbajo | outlyingness > factorHaciaArriba
    iEsAFiltrar[is.na(outlyingness)] <- FALSE
    iOutliers <- which(iEsAFiltrar)
    
    if (F) {
      i <- 1
      i <- which(fechasObservaciones==as.POSIXct('2014-01-01', tz=tz(fechasObservaciones[1])))
      i <- i + 1
      x <- residuos[i,]
      iEsAConservar[i,!iEsAConservar[i,],drop=F]
      valoresObservaciones[i, !iEsAConservar[i,],drop=F]
      residuos[i, ]
      valoresObservaciones[i, ]
      
      puntosObservaciones <- SpatialPointsDataFrame(coords = coordsObservaciones, data=data.frame(value=valoresObservaciones[i,]))
      mapearPuntosGGPlot(puntos = puntosObservaciones, shpBase = shpMask$shp, dibujar = F, dibujarTexto = T)
    }
    idx <- arrayInd(iOutliers, dim(valoresObservaciones))
    iFecha <- idx[,1]
    iEstacion <- idx[,2]
    
    fechas <- rownames(valoresObservaciones)[iFecha]
    if (is.null(fechas)) fechas <- as.character(iFecha)
    estaciones <- colnames(valoresObservaciones)[iEstacion]
    if (is.null(estaciones)) estaciones <- as.character(iEstacion)
    
    dfOutliers <- data.frame(iOutlier=iOutliers, iFecha=iFecha, iEstacion=iEstacion, fecha=fechas, 
                             estacion=estaciones, observacion=valoresObservaciones[iOutliers], 
                             estimacion=pred[iOutliers], outlyingness=outlyingness[iOutliers])
    dfOutliers <- dfOutliers[order(dfOutliers[,'iFecha'], abs(dfOutliers[,'outlyingness'])), ]
    
    # valoresObservaciones[!iEsAConservar] <- NA
    # grabarDatos(pathArchivoDatos = paste0(pathDatos, 'Estaciones/TempAireMinTodas_SRT_x3_RLM.txt'), fechas = fechasObservaciones, datos = valoresObservaciones, na = '-99')
    return(dfOutliers)
  } else {
    estaciones <- as.character(sapply(colnames(valoresObservaciones), FUN = function(x) { 
      rep(x, nrow(valoresObservaciones)) }))
    fechas <- rep(rownames(valoresObservaciones), ncol(valoresObservaciones))
    valores <- as.numeric(valoresObservaciones)
    estimados <- as.numeric(pred)
    stdDifs <- as.numeric(outlyingness)
    tiposOutliers <- rep(TTO_SinProblemasDetectados, length(valores))
    tiposOutliers[stdDifs < -factorHaciaAbajo] <- TTO_OutlierPorLoBajo
    tiposOutliers[stdDifs > factorHaciaArriba] <- TTO_OutlierPorLoAlto
    
    return(createDFTests(
      estacion = estaciones, fecha = fechas, valor = valores, estimado = estimados, 
      tipoOutlier = tiposOutliers, stdDif = stdDifs, reemplazar = FALSE))
  }
}

deteccionOutliersUniversalGriddingCV <- function(
    coordsObservaciones, fechasObservaciones, valoresObservaciones, params, pathsRegresores, 
    maxOutlyingness=3.5, maxNIters=5) {
  paramsAux <- params
  paramsAux$difMaxFiltradoDeOutliersRLM <- 0
  paramsAux$difMaxFiltradoDeOutliersCV <- 0
  # Hago dos pasadas porque al hacer validación cruzada la presencia del outlier puede hacer que un dato
  # bueno se identifique incorrectamente como outlier
  # En la primera pasada identifico todos los puntos que se apartan de su estimación por CV
  # La intención de la iteración es determinar si una estación es outlier por culpa de otra detectada como outlier o no.
  # Por ejemplo si hay dos estaciones muy cercanas en una región que ronde los 15 grados, una de ellas tiene 17 grados y la otra 5.
  # En este caso la cercanía de la estación de 5 hará que la estimacionCV de la de 17 sea muy baja, detectándola como outlier.
  # Análogamente la estación de 5 tendrá una estimación alta y también será detectada como outlier.
  # Para definir el problema, elimino ambos datos y los estimo con las demás estaciones. Al no estár la medida 5 la estimación de la 
  # de 17 será cercana a 15 con lo cual no saldrá detectada como outlier, sin embargo la estimación de la de 5 será también cercana a 15
  # resultando en una detección de outlier.
  # Si entre un paso y otro una estación sigue siendo outlier, quiere decir que no eran los demás posibles outliers que la condicionaban
  # Para evitar problemas de convergencia, una vez que un outlier entra al conjunto y permanece un paso, no lo saco más
  iPasada <- 1
  iOutliers <- integer(0)
  iOutliersPermanecidos <- integer(0)
  
  # Busco el mínimo overlap que me permita tener al menos 25 observaciones en todas las fechas para 
  # calcular las medianas y desviaciones
  overlap <- 0
  numDatos <- numeric(nrow(valoresObservaciones))
  repeat {
    ies <- vapply(1:nrow(valoresObservaciones), FUN=function(x, overlap) { return((x-overlap):(x+overlap))}, FUN.VALUE = integer(2 * overlap + 1), overlap)
    if (!is.matrix(ies)) ies <- matrix(ies, nrow=1, byrow = T)
    # Cuento la cantidad de datos disponibles con el overlap
    numDatos <- apply(ies, MARGIN = 2, FUN=function(x) {
      x <- x[x >= 1 & x <= nrow(valoresObservaciones)]
      return(sum(!is.na(valoresObservaciones[x,])))
    })
    
    if (min(numDatos) > min(25, sum(!is.na(valoresObservaciones)))) break
    overlap <- overlap + 1
  }
  
  diffs <- matrix(data = NA, nrow = nrow(valoresObservaciones), ncol = 0)
  params$verbose<- T
  repeat {
    iOutliersAnterior <- iOutliers
    aux <- valoresObservaciones
    aux[iOutliers] <- NA
    
    # valoresObservaciones = aux
    # params = paramsAux
    # eliminarSerieTemporalCompleta=F
    # estimarNAs = iPasada > 1
    # paramsAux$nCoresAUsar <- 1
    
    estimacionCV <- universalGriddingCV(coordsObservaciones=coordsObservaciones, fechasObservaciones=fechasObservaciones, valoresObservaciones = aux, 
                                        params = paramsAux, pathsRegresores = pathsRegresores, eliminarSerieTemporalCompleta = F, 
                                        estimarNAs = iPasada > 1)

    x <- valoresObservaciones - estimacionCV
    diffs <- cbind(diffs, x)
    outlyingness <- calcOutlyingnessMedianaMAD(x, xEstimarMedianaYMAD = diffs, porFilas = T, overlap = overlap)
    
    #estimacionCV[,'Rocha']
    #outlyingness[,'Rocha']
    
    # plot(1:length(outlyingness), outlyingness)
    # hist(valoresObservaciones - estimacionCV)
    iOutliers <- unique(c(iOutliersPermanecidos, which(abs(outlyingness) > maxOutlyingness)))
    iOutliersPermanecidos <- unique(c(iOutliersPermanecidos, base::intersect(iOutliersAnterior, iOutliers)))

    idx <- arrayInd(iOutliers, dim(valoresObservaciones))
    iFecha <- idx[,1]
    iEstacion <- idx[,2]
    iFechasConMasDeUnOutlier <- unique(iFecha[duplicated(iFecha)])
    
    if (params$verbose) {
      fechas <- rownames(valoresObservaciones)[iFecha]
      if (is.null(fechas)) fechas <- as.character(iFecha)
      estaciones <- colnames(valoresObservaciones)[iEstacion]
      if (is.null(estaciones)) estaciones <- as.character(iEstacion)
      
      dfOutliers <- data.frame(iOutlier=iOutliers, iFecha=iFecha, iEstacion=iEstacion, fecha=fechas, 
                               estacion=estaciones, observacion=valoresObservaciones[iOutliers], 
                               estimacion=estimacionCV[iOutliers], outlyingness=outlyingness[iOutliers])
      dfOutliers <- dfOutliers[order(dfOutliers[,'iFecha'], -abs(dfOutliers[,'outlyingness'])), ]
      
      print(iPasada)
      print(dfOutliers)
    }
    iPasada <- iPasada + 1

    if (identical(iOutliersAnterior, iOutliers) | length(iFechasConMasDeUnOutlier) == 0 | iPasada > maxNIters) break
  }
  
  if (length(iOutliers) > 0) {
    idx <- arrayInd(iOutliers, dim(valoresObservaciones))
    iFecha <- idx[,1]
    iEstacion <- idx[,2]
    
    fechas <- rownames(valoresObservaciones)[iFecha]
    if (is.null(fechas)) fechas <- as.character(iFecha)
    estaciones <- colnames(valoresObservaciones)[iEstacion]
    if (is.null(estaciones)) estaciones <- as.character(iEstacion)
    
    dfOutliers <- data.frame(iOutlier=iOutliers, iFecha=iFecha, iEstacion=iEstacion, fecha=fechas, 
                             estacion=estaciones, observacion=valoresObservaciones[iOutliers], 
                             estimacion=estimacionCV[iOutliers], outlyingness=outlyingness[iOutliers])
    dfOutliers <- dfOutliers[order(dfOutliers[,'iFecha'], abs(dfOutliers[,'outlyingness'])), ]
  } else { dfOutliers <- NULL }

  return(dfOutliers)
}

getPoligonoBoundingBox <- function(
    objSP, caja=bbox(objSP), outputCRS=NULL, factorExtensionX=1, 
    factorExtensionY=factorExtensionX) {
  if (factorExtensionX != 1) {
    # Agrando la caja (factorExtensionX - 1) * 0.5 hacia cada lado
    extX <- diff(caja[1,]) * (factorExtensionX - 1) * 0.5
    caja[1,] <- caja[1,] + c(-extX, extX)
  }
  if (factorExtensionY != 1) {
    # Agrando la caja (factorExtensionY - 1) * 0.5 hacia cada lado
    extY <- diff(caja[2,]) * (factorExtensionY - 1) * 0.5
    caja[2,] <- caja[2,] + c(-extY, extY)
  }
  
  boundingBox <- matrix(c(caja[1,1], caja[2,1], #xMin, yMin
                          caja[1,2], caja[2,1], #xMax, yMin
                          caja[1,2], caja[2,2], #xMax, yMax
                          caja[1,1], caja[2,2], #xMin, yMax
                          caja[1,1], caja[2,1]), ncol = 2, byrow = T)
  boundingBox <- Polygon(coords = boundingBox)
  boundingBoxPolygon <- SpatialPolygons(
    list(Polygons(list(boundingBox), ID = "BoundingBox")), proj4string = objSP@proj4string)
  if (is.null(outputCRS)) {
    return(boundingBoxPolygon)
  } else {
    return(spTransform(boundingBoxPolygon, outputCRS))
  }
}

grillaSobreBoundingBox <- function(
    objSP, caja=bbox(getPoligonoBoundingBox(objSP = objSP, outputCRS = objSP@proj4string)), 
    largoDimensiones=diff(t(caja)), nCeldasX=100, 
    nCeldasY = round(nCeldasX * largoDimensiones[2] / largoDimensiones[1])) {
  cellsDim <- c(nCeldasX, nCeldasY)
  cellSize <- as.numeric(largoDimensiones / cellsDim)
  cellcentreOffset <- as.numeric(caja[, 1] + cellSize * 0.5)
  
  return(SpatialGrid(GridTopology(
    cellcentre.offset = cellcentreOffset, cellsize = cellSize, cells.dim = cellsDim),
    proj4string = objSP@proj4string))
}

grillaPixelesSobreBoundingBox <- function(
    objSP, caja=bbox(getPoligonoBoundingBox(objSP = objSP, outputCRS = objSP@proj4string)), 
    largoDimensiones=diff(t(caja)), nCeldasX=100, 
    nCeldasY = round(nCeldasX * largoDimensiones[2] / largoDimensiones[1])) {
  return(as(grillaSobreBoundingBox(
      objSP=objSP, caja=caja, largoDimensiones=largoDimensiones, nCeldasX=nCeldasX, 
      nCeldasY = nCeldasY), 
    'SpatialPixels'))
}

filtradoRegresoresDiscontinuosPorVGM <- function(pathsRegresoresATestear, shpBase, nMuestras=5000, q1Dists=1/3, q2Dists=1-q1Dists, umbral=3) {
  # nMuestras=2500
  # q1Dists=1/3
  # q2Dists=1-q1Dists
  # umbral=3
  # pathsRegresoresATestear <- listaRegresores[[7]]
  
  iReg <- 3
  for (iReg in 1:ncol(pathsRegresoresATestear)) {
    iPrimerNoNA <- which(!is.na(pathsRegresoresATestear[, iReg]))[1]
    rasterI <- readGDAL(pathsRegresoresATestear[iPrimerNoNA, iReg], silent=T)
    
    # mapearGrillaGGPlot(rasterI, shpBase = shpBase, continuo = T, dibujar = F)
    puntosRasterI <- as(object = rasterI, Class = 'SpatialPoints')
    
    # iEsSobreSHPBase <- which(!is.na(over(geometry(rasterI), geometry(shpBase))))
    
    muestras <- muestrearEnCuadrantes(rasterI, size = nMuestras, nCuadrantesX = 5)
    distancias <- spDists(puntosRasterI[muestras, ])
    
    dists <- numeric(length(muestras) * (length(muestras) - 1) / 2)
    k <- 1
    i <- 1
    for (i in 1:(length(muestras)-1)) {
      js <- (i+1):length(muestras)
      ks <- k:(k+length(js)-1)
      
      dists[ks] <- distancias[i, js]
      
      k <- ks[length(ks)] + 1
    }
    
    qs <- quantile(dists, probs=c(q1Dists, q2Dists))
    q1 <- qs[1]
    q2 <- qs[2]
    iDistsQ1 <- which(dists < q1)
    iDistsQ2 <- which(dists > q2)
    
    diffs <- numeric(length(dists))
    
    
    k <- 1
    i <- 1
    for (i in 1:(length(muestras)-1)) {
      js <- (i+1):length(muestras)
      ks <- k:(k+length(js)-1)
      
      diffs[ks] <- ((rasterI@data[muestras[i],1] - rasterI@data[muestras[js],1])^2)/2
      
      k <- ks[length(ks)] + 1
    }
    #plot(dists, diffs)
    rm(distancias, dists)
    
    #lala <- lm(formula = y ~ x + 1, data = data.frame(x=dists, y=diffs))
    #lala
    
    median(diffs[iDistsQ2], na.rm=T) / median(diffs[iDistsQ1], na.rm=T)
  }
}
