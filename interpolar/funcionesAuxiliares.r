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
if (iFrame >= 3) { script.dir.funcionesAuxiliares <- sys.frame(iFrame - 3)$ofile 
} else { script.dir.funcionesAuxiliares <- NULL }
while ((is.null(script.dir.funcionesAuxiliares) || is.na(regexpr('funcionesAuxiliares.r', script.dir.funcionesAuxiliares, fixed=T)[1])) && iFrame >= 0) {
  script.dir.funcionesAuxiliares <- sys.frame(iFrame)$ofile
  iFrame <- iFrame - 1
}
if (is.null(script.dir.funcionesAuxiliares)) { script.dir.funcionesAuxiliares <- ''
} else { script.dir.funcionesAuxiliares <- paste(dirname(script.dir.funcionesAuxiliares), '/', sep='') }

source(paste0(script.dir.funcionesAuxiliares, '../instalarPaquetes/instant_pkgs.r'))
source(paste0(script.dir.funcionesAuxiliares, '../cacheFunciones/cacheFunciones.r'))
source(paste0(script.dir.funcionesAuxiliares, '../Graficas/graficas.r'))
source(paste0(script.dir.funcionesAuxiliares, '../sysutils/sysutils.r'))
instant_pkgs(
  pkgs = c('matrixStats', 'stringi', 'stringr', 'lubridate', 'gridExtra', 'sp', 'rgdal', 'Rmisc', 
           'gstat', 'grid'))


graficarEstaciones <- function(estaciones, serieComparacion, tempAireMin, tempAireMax, fechaMin=NA, fechaMax=NA, 
                               lWidth=0.5, carpetaSalida='Resultados/Graficos/', prefijoArchivos='', nombresSeries=c('MODIS', 'TMin', 'TMax'),
                               width=NA, height=NA, tipoDeArchivoSalida='.png') {
  if (is.na(width)) {
    if (tipoDeArchivoSalida == '.png') { width <- 1200 
    } else { width <- 16 }
  }
  
  if (is.na(height)) {
    if (tipoDeArchivoSalida == '.png') { height <- 660 
    } else { height <- 9 }
  }
  
  if (!is.na(fechaMin)) { tIni <- which(tempAireMin$fechas >= as.POSIXct(fechaMin, tz=tz(tempAireMin$fechas[1])))[1]
  } else { tIni <- 1 }
  if (!is.na(fechaMax)) { 
    tFin <- which(tempAireMin$fechas <= as.POSIXct(fechaMax, tz=tz(tempAireMin$fechas[1])))
    tFin <- tFin[length(tFin)]
  } else { tFin <- length(tempAireMin$fechas) }
  ts <- seq.int(tIni, tFin)
  
  fechas <- tempAireMin$fechas[ts]
  serieComp <- serieComparacion[ts, ]
  TMin <- tempAireMin$datos[ts, ]
  TMax <- tempAireMax$datos[ts, ]
  
  
  crearDirectoriosSiNoExisten(carpetaSalida)
  for (iEstacion in 1:nrow(estaciones)) {
    if (!all(is.na(TMin[, iEstacion])) || !all(is.na(TMax[, iEstacion]))) {
      yMin <- min(serieComp[, iEstacion], TMin[, iEstacion], TMax[, iEstacion], na.rm = T)
      yMax <- max(serieComp[, iEstacion], TMin[, iEstacion], TMax[, iEstacion], na.rm = T) + 9
      
      if (prefijoArchivos != '') { archSalida <- paste(carpetaSalida, prefijoArchivos, '_', estaciones$Nombre[iEstacion], tipoDeArchivoSalida, sep='')
      } else { archSalida <- paste(carpetaSalida, estaciones$Nombre[iEstacion], tipoDeArchivoSalida, sep='') }
      
      iNoNA <- !is.na(serieComp[, iEstacion]) & !is.na(TMin[, iEstacion])
      if (any(iNoNA)) {
        CORR_TMin <- round(cor(serieComp[, iEstacion], TMin[, iEstacion], use = "complete.obs"), 2)
        m <- lm(formula = 'y~x+1', data = data.frame(x=serieComp[, iEstacion], y=TMin[, iEstacion]))
        R2_TMin <- round(summary(m)$adj.r.squared, 2)
      } else {
        CORR_TMin <- NA
        R2_TMin <- NA
      }
      
      iNoNA <- !is.na(serieComp[, iEstacion]) & !is.na(TMax[, iEstacion])
      if (any(iNoNA)) {
        CORR_TMax <- round(cor(serieComp[, iEstacion], TMax[, iEstacion], use = "complete.obs"), 2)
        m <- lm(formula = 'y~x+1', data = data.frame(x=serieComp[, iEstacion], y=TMax[, iEstacion]))
        R2_TMax <- round(summary(m)$adj.r.squared, 2)
      } else {
        CORR_TMax <- NA
        R2_TMax <- NA
      }
      
      if (tipoDeArchivoSalida == '.emf') { emf(archSalida, width = width, height = height)
      } else if (tipoDeArchivoSalida == '.svg') { svg(archSalida, width = width, height = height)
      } else { png(archSalida, width = width, height = height) }
      
      if (prefijoArchivos != '') { titulo <- paste(estacionesINUMET$Nombre[iEstacion], ' ', prefijoArchivos, sep='')
      } else { titulo <- estacionesINUMET$Nombre[iEstacion] }
      
      plot(x = fechas, y = serieComp[, iEstacion], type = 'l', col = "green", lwd=lWidth,
           ylim = c(yMin, yMax), main = titulo, xlab = 'Fecha', ylab = 'Temperatura')
      lines(x = fechas, y = TMin[, iEstacion], col = "blue", lwd=lWidth)
      lines(x = fechas, y = TMax[, iEstacion], col = "red", lwd=lWidth)
      
      legend(tempAireMin$fechas[ts[1]], yMax, # places a legend at the appropriate place
             c(nombresSeries[1], 
               paste(nombresSeries[2], '. Corr=', CORR_TMin, '. R^2=', R2_TMin, sep=''),
               paste(nombresSeries[3], '. Corr=', CORR_TMax, '. R^2=', R2_TMax, sep='')), # puts text in the legend
             lty=c(1, 1, 1), # gives the legend appropriate symbols (lines)
             pch=c(-1, -1, -1),
             lwd=c(lWidth, lWidth, lWidth), col=c('green', 'blue', 'red')) # gives the legend lines the correct color and width
      
      dev.off()
    } 
  }
}

graficarEstacionesSepararPorAnios <- function(estaciones, serieComparacion, tempAireMin, tempAireMax, lWidth=0.5, 
                                              carpetaSalida='Resultados/Graficos/', nombresSeries=c('MODIS', 'TMin', 'TMax'),
                                              width=NA, height=NA, tipoDeArchivoSalida='.png') {
  anios <- unique(year(tempAireMin$fechas))
  
  iAnio <- anios[1]
  for (iAnio in anios) {
    graficarEstaciones(estaciones = estaciones, serieComparacion = serieComparacion, tempAireMin = tempAireMin, tempAireMax = tempAireMax, 
                       lWidth = lWidth,  fechaMin = paste(iAnio, '-01-01', sep=''), fechaMax = paste(iAnio, '-12-31', sep=''), 
                       carpetaSalida=carpetaSalida, prefijoArchivos = as.character(iAnio), nombresSeries = nombresSeries, 
                       width = width, height = height, tipoDeArchivoSalida = tipoDeArchivoSalida)
  }
  
  graficarEstaciones(estaciones = estaciones, serieComparacion=serieComparacion, tempAireMin=tempAireMin, lWidth = lWidth, width = width, height = height, 
                     tipoDeArchivoSalida = tipoDeArchivoSalida, tempAireMax=tempAireMax, carpetaSalida=carpetaSalida, nombresSeries = nombresSeries)
  
}

nombreModelo <- function(params, pathsRegresores=NULL) {
  if (!is.null(pathsRegresores) && !is.null(colnames(pathsRegresores))) { strRegresores <- paste(colnames(pathsRegresores), collapse = '+')
  } else { strRegresores <- '' }
  
  if (!is.null(params$tlagsAR)) {
    if (strRegresores == '') { strRegresores <- paste('AR(', paste(params$tlagsAR, collapse=','),')',sep='')
    } else { strRegresores <- paste(strRegresores, paste('AR(', paste(params$tlagsAR, collapse=','),')',sep=''), sep='+') }
  }
  
  if (params$incorporarAltitud) {
    if (strRegresores != '') { strRegresores <- paste(strRegresores, gsub(pattern = ' ', replacement = '', x = params$formulaAltitud, fixed = T), sep='+')
    } else { strRegresores <- gsub(pattern = ' ', replacement = '', x = params$formulaAltitud, fixed = T) }
  }
    
  if (params$incorporarCoordenadas)  {
    if (strRegresores != '') { strRegresores <- paste(strRegresores, gsub(pattern = ' ', replacement = '', x = params$formulaCoordenadas, fixed = T), sep='+')
    } else { strRegresores <- gsub(pattern = ' ', replacement = '', x = params$formulaCoordenadas, fixed = T) }
  }

  if (params$incorporarDistanciaAlAgua)  {
    if (strRegresores != '') { strRegresores <- paste(strRegresores, gsub(pattern = ' ', replacement = '', x = params$formulaDistanciaAlAgua, fixed = T), sep='+')
    } else { strRegresores <- gsub(pattern = ' ', replacement = '', x = params$formulaDistanciaAlAgua, fixed = T) }
  }
  
  if (params$metodoIgualacionDistribuciones == 'ninguna' || (is.null(pathsRegresores) &&  
      !(params$incorporarCoordenadas || params$incorporarTiempo || params$incorporarDistanciaAlAgua || 
        params$incorporarAltitud))) { strMetodoIgualacion <- ''
  } else if (params$metodoIgualacionDistribuciones == 'regresionLineal') { strMetodoIgualacion <- 'R'
  } else if (params$metodoIgualacionDistribuciones == 'regresionLinealRobusta') { strMetodoIgualacion <- 'RR'
  } else if (params$metodoIgualacionDistribuciones == 'regresionLinealConEliminacionDeOutliers') { strMetodoIgualacion <- 'REO'
  } else if (params$metodoIgualacionDistribuciones == 'CDFMatching') { strMetodoIgualacion <- 'CDF'
  } else if (params$metodoIgualacionDistribuciones == 'GLS') { strMetodoIgualacion <- 'GR'
  } else { strMetodoIgualacion <- params$metodoIgualacionDistribuciones }
  
  if (params$interpolationMethod == 'none') { strInterpolacion <- ''
  } else if (params$interpolationMethod == 'idw') { strInterpolacion <- 'IDW'
  } else if (params$interpolationMethod == 'automap') { 
    if (is.logical(params$usarFitVariogramGLS) && !params$usarFitVariogramGLS) { strInterpolacion <- 'K_NVGLS'
    } else { strInterpolacion <- 'K' }
  } else if (params$interpolationMethod == 'copula') { strInterpolacion <- 'Copula'
  } else if (params$interpolationMethod == 'stUniversalKriging') {
    strInterpolacion <- 'STK'
    if (strRegresores =='') { strMetodoIgualacion <- ''
    } else { strMetodoIgualacion <- 'R' }
  } else { strInterpolacion <- params$interpolationMethod
  }
  
  if (!is.null(params$preECDFMatching) && params$preECDFMatching) { strECDF <- '_ECDF'
  } else { strECDF <- '' }
  
  strMetodo <- paste0(strMetodoIgualacion, strInterpolacion)
  
  if (strMetodo != '') {
    if (any(strRegresores != '')) { return(paste(strMetodo, paste(strRegresores, collapse = '+'), sep='-'))
    } else {return(strMetodo)}
  } else { return(strRegresores) }
}

findDateInText <- function(x, tz) {
  matches <- gregexpr(
    pattern='(19|20)[0-9][0-9][- /._]*(0[1-9]|1[012])[- /._]*(0[1-9]|[12][0-9]|3[01])', text=x)
  fechasReg <- gsub(
    pattern = '[ /._]', replacement = '-', x = unlist(regmatches(x = x, m = matches)))
  idx <- !grepl(pattern = '-', x=fechasReg, fixed = TRUE)
  fechasReg[idx] <- paste(
    substr(fechasReg[idx], 1, 4), '-', 
    substr(fechasReg[idx], 5, 6), '-', 
    substr(fechasReg[idx], 7, 8), sep='')
  return(as.POSIXct(fechasReg, tz=tz)) 
}

cargarRegresor <- function(carpetaRegresor, fechasRegresando, funcFechasRegresor=findDateInText, 
                           expandir=F) {
  pathsRegresor <- sort(dir(carpetaRegresor, pattern = '*.tif$', full.names = T, include.dirs = F))
  if (length(pathsRegresor) > 0) {
    fechasRegresor <- funcFechasRegresor(x=pathsRegresor, tz = tz(fechasRegresando[1])) 
    res <- matrix(expandirPathsRegresor(
      fechasRegresor = fechasRegresor, pathsRegresor = pathsRegresor, 
      fechasRegresando = fechasRegresando, expandir=expandir), ncol = 1)
  } else { res <- matrix(NA, nrow = length(fechasRegresando), ncol = 1) }
  colnames(res) <- ultimaCarpeta(carpetaRegresor, removerTrailingSlash = F)
  return(res)
}

cargarRegresores <- function(carpetaRegresores, fechasRegresando, funcFechasRegresor=findDateInText, 
                             expandir=F) {
  carpetasRegresores <- setdiff(list.dirs(carpetaRegresores, recursive = F, full.names = F), 'RCache')
  
  return(cargarRegresoresEx(
    carpetaRegresores=carpetaRegresores, carpetasRegresores=carpetasRegresores, 
    fechasRegresando=fechasRegresando, funcFechasRegresor=funcFechasRegresor, expandir=expandir))
}

cargarRegresoresEx <- function(carpetaRegresores, carpetasRegresores, fechasRegresando, 
                               funcFechasRegresor=findDateInText, expandir=F) {
  stopifnot(length(carpetasRegresores)==length(expandir) || length(expandir)==1)
  
  if (length(expandir) != length(carpetasRegresores)) expandir <- rep(expandir[1], length(carpetasRegresores))
  
  res <- matrix(ncol=length(carpetasRegresores), nrow = length(fechasRegresando))
  i <- 2
  for (i in 1:length(carpetasRegresores)) {
    pathReg <- paste(carpetaRegresores, '/', carpetasRegresores[i], sep='')
    
    pathsRegresor <- sort(dir(pathReg, pattern = '*.tif$', full.names = T, include.dirs = F))
    if (length(pathsRegresor) > 0) {
      fechasRegresor <- funcFechasRegresor(x=pathsRegresor, tz = tz(fechasRegresando[1])) 
      res[, i] <- expandirPathsRegresor(
        fechasRegresor=fechasRegresor, pathsRegresor=pathsRegresor, 
        fechasRegresando=fechasRegresando, expandir=expandir[i])
    } else {
      res[, i] <- rep(NA_character_, length(fechasRegresando))
    }
  }
  row.names(res) <- as.character(fechasRegresando)
  colnames(res) <- carpetasRegresores
  return(res)
}

errorEnRegresor <- function(i, pathsRegresores) {
  if (!is.na(pathsRegresores[i])) {
    require('rgdal')
    res <- try(readGDAL(pathsRegresores[i], silent=T))
    return("try-error" %in% class(res))
  } else { return(FALSE) }
}

erroresEnRegresores <- function(pathsRegresores) {
  nCoresAUsar <- detectCores(T, T)
  if (nCoresAUsar > 1) {
    cl <- makeCluster(getOption("cl.cores", nCoresAUsar))
    if (exists(x = 'setMKLthreads')) { clusterEvalQ(cl = cl, expr = setMKLthreads(1)) }
    errores <- parSapplyLB(cl=cl, X=1:length(pathsRegresores), FUN=errorEnRegresor, pathsRegresores=pathsRegresores)
    stopCluster(cl)
  } else {
    errores <- sapply(X=1:length(pathsRegresores), FUN=errorEnRegresor, pathsRegresores=pathsRegresores)
  }
  return(pathsRegresores[errores])
}

chequearYActualizarRegresores <- function(pathsRegresores, carpetasOriginalesRegresores) {
  errores <- character(0)
  if (ncol(pathsRegresores) != length(carpetasOriginalesRegresores)) stop('carpetasOriginalesRegresores erroneo')
  for (i in 1:length(pathsRegresores)) {
    if (!is.na(pathsRegresores[i])) {
      res <- try(readGDAL(pathsRegresores[i], silent=T))
      if ("try-error" %in% class(res)) {
        errores <- c(errores, pathsRegresores[i])
        print(pathsRegresores[i])
        iCol <- trunc((i - 1) / nrow(pathsRegresores)) + 1
        pathOrig <- paste(carpetasOriginalesRegresores[iCol], basename(pathsRegresores[i]), sep='') 
        file.copy(from=pathOrig, to=pathsRegresores[i], overwrite = T)
      }
    }
  }
  return(errores)
}

conteoPixelesRegresor <- function(pathsRegresor, shpMask=NULL, zcol=1) {
  # pathsRegresor <- pathsRegresores[, 1]
  # shpMask <- shpMaskNoNulos
  
  i <- 1
  evaluarConReintentos(regresor <- readGDAL(pathsRegresor[1], silent=T))
  while (is.null(regresor)) {
    i <- i + 1
    evaluarConReintentos(regresor <- readGDAL(pathsRegresor[i], silent=T))
  }
  
  if (is.null(shpMask)) { mask <- rep(TRUE, nrow(regresor))
  } else { mask <- shpMask$mask }
  
  
  esNulo <- is.na(regresor@data[mask, zcol])
  numDiasNoNulos <- as.integer(!esNulo)
  maxCantDiasNulosConsecutivos <- numDiasNoNulos
  cantDiasNulosConsecutivos <- maxCantDiasNulosConsecutivos
  
  for (i in (i+1):length(pathsRegresor)) {
    print(i)
    evaluarConReintentos(regresor <- readGDAL(pathsRegresor[i], silent=T))
    
    if (!is.null(regresor)) { esNulo <- is.na(regresor@data[mask, zcol])
    } else { esNulo <- rep(TRUE, length(esNulo)) }
    numDiasNoNulos <- numDiasNoNulos + as.integer(!esNulo)
    
    cantDiasNulosConsecutivos <- (cantDiasNulosConsecutivos + as.integer(esNulo)) * as.integer(esNulo)
    maxCantDiasNulosConsecutivos <- apply(X = cbind(maxCantDiasNulosConsecutivos, cantDiasNulosConsecutivos), FUN = max, MARGIN = 1)
  }
  
  regresor@data[mask, zcol] <- numDiasNoNulos
  regresor@data[!mask, zcol] <- NA
  i <- which(regresor@data[, zcol] < 2200)
  regresor@data[i, zcol] <- NA
  
  regresor@data[, zcol] <- (regresor@data[, zcol] / length(pathsRegresor)) * 100
  
  range(regresor@data[, zcol], na.rm = T)
  escala <- crearEscalaEquiespaciada(datos = regresor@data[,zcol], nDigitos = 0, nIntervalos = 11)
  mapearGrillaGGPlot(grilla = regresor, shpBase = shpMask$shp, zcol = zcol, escala = escala, titulo = 'Porcentaje de Días con Datos', 
                     nomArchResultados = 'Resultados/PorcDiasConDatos.png')
  
  regresor@data[mask, zcol] <- maxCantDiasNulosConsecutivos
  regresor@data[!mask, zcol] <- NA
  i <- which(regresor@data[, zcol] > 1000)
  regresor@data[i, zcol] <- NA
  
  range(regresor@data[, zcol], na.rm = T)
  escala <- crearEscalaEquiespaciada(datos = regresor@data[,zcol], nDigitos = 0, nIntervalos = 11)
  mapearGrillaGGPlot(grilla = regresor, shpBase = shpMask$shp, zcol = zcol, escala = escala, titulo = 'Máxima Cantidad de Días Nulos Consecutivos', 
                     nomArchResultados = 'Resultados/MaxNDiasNulosConsecutivos.png')
  
  regresor@data[mask, zcol] <- maxCantDiasNulosConsecutivos
  regresor@data[!mask, zcol] <- NA
  i <- which(regresor@data[, zcol] > 9)
  regresor@data[i, zcol] <- NA
  
  range(regresor@data[, zcol], na.rm = T)
  nIntervalos <- max(regresor@data[, zcol], na.rm = T) - min(regresor@data[, zcol], na.rm = T) + 1
  
  escala <- crearEscalaEquiespaciada(datos = c(5, 9.1), nDigitos = 0, nIntervalos = nIntervalos)
  mapearGrillaGGPlot(grilla = regresor, shpBase = shpMask$shp, zcol = zcol, escala = escala, titulo = 'Máxima Cantidad de Días Nulos Consecutivos (< 10 Días)',
                     nomArchResultados = 'Resultados/MaxNDiasNulosConsecutivosMenorA10.png')
}


lala <- function() {
  #install.packages("spacetime")
  
  #install.packages("GSIF")
  #library(GSIF)
  #URI = "http://wps.worldgrids.org/pywps.cgi"
  #server <- list(URI=URI, request="execute", version="version=1.0.0", service.name="service=wps",
  #               identifier="identifier=sampler_local1pt_nogml")
  #TWISRE3 <- new("WPS", server=server, inRastername="twisre3")
  #str(TWISRE3)
  #rm(TWISRE3)
  
  twisre3 <- readGDAL(paste(pathDatos, 'TWISRE3/TWISRE3a.tif', sep=''))
  twisre3_uy <- spTransform(coordsAInterpolar, twisre3@proj4string)
  res <- over(x = twisre3_uy, twisre3, returnList = F)[,1]
  twisre3_uy <- SpatialPixelsDataFrame(coordsAInterpolar, data = data.frame(TWI=(res/10+10)))
  shpMask <- cargarSHPYObtenerMascaraParaGrilla(pathSHP = params$pathSHPMapaBase, grilla = twisre3_uy, spSinMascara = estacionesINUMET)
  mapearGrillaGGPlot(twisre3_uy, shpMask$shp, nomArchResultados = paste(pathDatos, 'TWISRE3/twisre3_uy.png', sep=''), continuo = T)
  
  estacionesINUMET <- spTransform(estacionesINUMET, twisre3_uy@proj4string)
  iNA <- which(is.na(twisre3_uy@data[shpMask$mask, ]))
  coordsNA <- SpatialPoints(sp::coordinates(aux2)[iNA,], proj4string = twisre3_uy@proj4string)
  
  twisre3_uy <- rellenarSP(sp = twisre3_uy, mascara = shpMask$mask, metodo='idw', zcol=1)
  mapearGrillaGGPlot(twisre3_uy, shpMask$shp, nomArchResultados = paste(pathDatos, 'TWISRE3/twisre3_uy.png', sep=''), continuo = T)
  
  if (FALSE) {
    iEstacionesNA <- which(is.na(over(estacionesINUMET, twisre3_uy)))
    coordsEstacionesNA <- sp::coordinates(estacionesINUMET)[iEstacionesNA, , drop=F]
    coordsEstacionesNA <- SpatialPoints(
      coords = coordsEstacionesNA, proj4string = estacionesINUMET@proj4string)
    aux3 <- SpatialPoints(coords = sp::coordinates(twisre3_uy), proj4string = aux2@proj4string)
    dist <- gDistance(coordsEstacionesNA, aux3, byid = T)
    i <- order(dist)
    twisre3_uy@data[over(coordsEstacionesNA, geometry(twisre3_uy)), 1] <- twisre3_uy@data[which.min(i), 1]
    over(estacionesINUMET, twisre3_uy)
    twisre3_uy <- SpatialPixelsDataFrame(coordsAInterpolar, data = data.frame(TWI=twisre3_uy@data$TWI))
    mapearGrillaGGPlot(twisre3_uy, shpMask$shp, nomArchResultados = paste(pathDatos, 'TWISRE3/twisre3_uy.png', sep=''), continuo = T)
  }
  
  writeGDAL(twisre3_uy, paste(pathDatos, 'TWISRE3/twisre3_uy.tif', sep=''), options = c('COMPRESS=DEFLATE', 'PREDICTOR=2', 'ZLEVEL=9'))
  
  lalala <- readGDAL(paste(pathDatos, 'TWISRE3/twisre3_uy.tif', sep=''))
  mapearGrillaGGPlot(lalala, shpMask$shp, nomArchResultados = paste(pathDatos, 'TWISRE3/twisre3_uy.png', sep=''), continuo = T)
  
  rm(twisre3, twisre3_uy, res, coordsPDelEste, aux3, lalala)
  
  proj4StringLatLong <- '+proj=longlat +datum=WGS84'
  wktLatLong <- "EPSG:4326"
  coordsAInterpolarLatLong <- spTransform(
    coordsAInterpolar, CRSobj = CRS(projargs = proj4StringLatLong, SRS_string = wktLatLong))
  
  lat <- sp::coordinates(coordsAInterpolarLatLong)[, 2] * pi / 180
  range(cos(lat * pi / 180))
  
  minT <- 1000
  maxT <- -1000
  
  t <- 180
  for (t in 1:365) {
    fi <- (t - 18) * 2 * pi / 365 + 2 ^ (1 - sign(lat)) * pi
    tMinGeom <- 24.2 * cos(lat) - 15.7 * (1 - cos(fi)) * sin(abs(lat))
    
    minT <- min(minT, min(tMinGeom))
    maxT <- max(maxT, max(tMinGeom))
    
    #aux2 <- SpatialPixelsDataFrame(coordsAInterpolar, data = data.frame(tMinGeom=tMinGeom))
    #writeGDAL(aux2, paste('Datos/TGeomTrend/', t, '_uy.tif'))
  }
  
  escalaTMinGeom <- crearEscalaEquiespaciada(c(minT, maxT), nDigitos = 1, nIntervalos = 11)
  for (t in 1:365) {
    fi <- (t - 18) * 2 * pi / 365 + 2 ^ (1 - sign(lat)) * pi
    tMinGeom <- 24.2 * cos(lat) - 15.7 * (1 - cos(fi)) * sin(abs(lat))
    
    aux2 <- SpatialPixelsDataFrame(coordsAInterpolar, data = data.frame(tMinGeom=tMinGeom))
    writeGDAL(aux2, paste(pathDatos, 'TGeomTrend/', t, '_uy.tif', sep=''), options = c('COMPRESS=DEFLATE', 'PREDICTOR=2', 'ZLEVEL=9'))
    mapearGrillaGGPlot(aux2, shpMask$shp, escala = escalaTMinGeom, nomArchResultados =  paste(pathDatos, 'TGeomTrend/plots/', t, '_uy.png', sep=''), continuo = T)
  }
}

testsMultiRegresores <- function() {
  aux <- pathsRegresores[, c(5, 18, 19)]
  colnames(aux)
  
  res <- matrix(nrow=ncol(pathsRegresores), ncol = 3)
  colnames(res) <- c('Pearson', 'Spearman', 'Adj. R^2')
  rownames(res) <- colnames(pathsRegresores)
  
  serie <- tempAireMin$datos
  valoresRegresores <- extraerValoresRegresoresSobreSP(objSP = estacionesINUMET, pathsRegresores = aux, silent = F)
  
  iNoNaRegresores <- lapply(valoresRegresores, FUN = function(x) { !is.na(x) })
  
  i<-1
  iNoNa <- !is.na(serie)
  for (i in 1:length(iNoNaRegresores))
    iNoNa <- iNoNa & iNoNaRegresores[[i]]
  
  if (any(iNoNa)) {
    s <- serie[iNoNa]
    r <- matrix(nrow = length(s), ncol = length(valoresRegresores))
    for (i in 1:length(valoresRegresores))
      r[, i] <- valoresRegresores[[i]][iNoNa] 
    
    colnames(r) <- colnames(aux)
    
    df <- data.frame(y=s, x=r)
    colnames(df)[2:ncol(df)] <- colnames(aux)
    formulaModelo <- as.formula(paste(colnames(df)[1], '~', paste('+', colnames(df)[2:ncol(df)], collapse=''), '+1'))
    
    m <- lm(formula = formulaModelo, data = df)
    smry <- summary(m)
    pred <- predict(m, newdata = data.frame(r))
    
    res[i, 1] <- cor(s, pred)
    res[i, 2] <- cor(s, pred, method = 'spearman')
    res[i, 3] <- smry$adj.r.squared
    
    
    plot(x=s, y=pred)
    
    #arch <- paste(carpeta, nombreSerie, '_vs_', colnames(pathsRegresores)[i], '.png', sep='')
    
    #png(arch)
    #plot(x=r, y=s, xlab=colnames(pathsRegresores)[i], ylab=nombreSerie)
    #dev.off()
  }
  smry
  res[i, ]
}

cargarValores <- function(pathsRaster, zcol=1) {
  # pathsRaster <- dir(paste(pathDatos, "LST_Night_Combinada_Clim", sep=""), full.names = T)
  require(rgdal)
  vals <- numeric()
  for (i in 1:length(pathsRaster)) {
    if (!is.na(pathsRaster[i])) {
      rasterI <- readGDAL(pathsRaster[i], silent=T)
      vals <- c(vals, rasterI@data[, zcol])
    }
  }
  return(vals)
}

muestrearValoresI <- function(i, pathsRaster, zcol=1, nMuestras=5000) {
  print(i)
  require(rgdal)
  if (!is.na(pathsRaster[i]) & file.exists(pathsRaster[i])) {
    rasterI <- readGDAL(pathsRaster[i], silent=T)
    iMuestra <- sample.int(nrow(rasterI), size=min(nrow(rasterI), nMuestras), replace = FALSE)
    return(rasterI@data[iMuestra, zcol])
  } else {
    return(rep(NA_real_, nMuestras))
  }
}

muestrearValores <- function(pathsRaster, zcol=1, nMuestras=5000, nCoresAUsar=0) {
  if (nCoresAUsar <= 0) { nCoresAUsar <- min(detectCores(T, T), length(pathsRaster)) } 
  
  if (nCoresAUsar > 1) {
    cl <- makeCluster(getOption("cl.cores", nCoresAUsar))
    if (exists(x = 'setMKLthreads')) { clusterEvalQ(cl = cl, expr = setMKLthreads(1)) }
    clusterEvalQ(cl = cl, expr = { require(rgdal) })
    valores <- parSapplyLB(cl=cl, X=1:length(pathsRaster), FUN=muestrearValoresI, pathsRaster=pathsRaster, zcol=zcol, nMuestras=nMuestras)
    stopCluster(cl)
  } else {
    valores <- sapply(X=1:length(pathsRaster), FUN=muestrearValoresI, pathsRaster=pathsRaster, zcol=zcol, nMuestras=nMuestras)
  }
  return(valores)
}

crearEscalaEquiespaciadaMultiRasters <- function(pathsRasters, nIntervalos=10, nDigitos=1, colores=NULL, brewerPal='Spectral', continuo=T, 
                                                 usarMuestreo=length(pathsRasters) > 1000, nMuestrasPorRaster=25000, lowTrimPU=0, highTrimPU=lowTrimPU) {
  if (usarMuestreo) { muestras <- muestrearValores(pathsRaster = as.vector(pathsRasters), nMuestras = nMuestrasPorRaster)
  } else { muestras <- cargarValores(as.vector(pathsRasters)) }
  muestras <- sort(muestras[!is.na(muestras)])
  
  iLowTrim <- round(length(muestras) * lowTrimPU + 1)
  iHighTrim <- round(length(muestras) * (1 - highTrimPU))
  
  return(crearEscalaEquiespaciada(datos = muestras[c(iLowTrim, iHighTrim)], nDigitos = nDigitos, nIntervalos = nIntervalos, colores = colores, 
                                  brewerPal = brewerPal, continuo = continuo))
}

plotRasterI <- function(i, pathsRaster, shpBase=NULL, escala=NULL, carpetaSalida='plots', titulos=rep('', length(pathsRaster)), replot=F, 
                        widthPx=1920, heightPx=1080, alturaEscalaContinua = unit(x=0.7, units = 'in')) {
  print(i)
  # pathsRaster <- dir(paste(pathDatos, 'LST_Night_Combinada_Clim', sep=''), full.names = T)
  # carpetaSalida <- paste(pathDatos, 'LST_Night_Combinada_Clim/plots/', sep='')
  # shpMask <- cargarSHPYObtenerMascaraParaGrilla(params$pathSHPMapaBase, grilla = coords)
  # carpetaSalida <- paste(pathDatos, 'LST_Night_Combinada_Clim/plots/', sep='')
  
  if (!is.na(pathsRaster[i])) {
    archivoSalida <- paste(carpetaSalida, changeFileExt(filename = basename(pathsRaster[i]), nuevaExtensionConPunto = '.png'), sep='')
    
    if (replot || !file.exists(archivoSalida) || file.info(archivoSalida)$size == 0) {
      rasterI <- readGDAL(pathsRaster[i], silent = T)
      if (!is.null(shpBase) && !identicalCRS(rasterI, shpBase)) { 
        shpBase <- spTransform(shpBase, rasterI@proj4string)
      }
      if (!is.null(escala) && !all(is.na(rasterI@data[, 1]))) {
        escala <- ajustarExtremosEscala(escala = escala, datos = rasterI@data[, 1], redondear = F)
      }
      
      mapearGrillaGGPlot(
        grilla = rasterI, shpBase = shpBase, escala = escala, nomArchResultados = archivoSalida, 
        continuo = T, dibujar=F, titulo = titulos[i], widthPx = widthPx, heightPx = heightPx, 
        alturaEscalaContinua = alturaEscalaContinua)
    }
  }
  return(NULL)
}

plotRasters <- function(pathsRaster, shpBase=NULL, carpetaSalida=paste(dirname(pathsRaster[1]), '/plots/', sep=''), escala=NULL, replot=F, 
                        widthPx=1080, heightPx=1080, alturaEscalaContinua = unit(x=0.7, units = 'in'),
                        titulos = paste(rownames(pathsRaster), '-', colnames(pathsRaster)[1])) {
  # pathsRaster <- pathsRegresores[, 'MOD11A1_LST_Night', drop=F]
  # pathsRaster <- pathsRegresores[, 'MOD11A1_LST_Night_9', drop=F]
  # pathsRaster <- pathsRegresores[, 'MOD11A1_LST_Night_R', drop=F]
  # pathsRaster <- pathsRegresores[, 'MOD11A1_LST_Night_FR', drop=F]
  # pathsRaster <- pathsRegresores[, 'MOD11A1_LST_Night_FRv2', drop=F]
  # pathsRaster <- pathsRegresores[, 'MYD11A1_LST_Night_FRv2', drop=F]
  # pathsRaster <- pathsRegresores[, 'MYD11A1_LST_Night', drop=F]
  # pathsRaster <- dir(paste(pathDatos, 'LST_Night_Combinada_Filtrada', sep=''), pattern = '*.tif', full.names = T)
  
  # titulos <- paste(rownames(pathsRaster), '-', colnames(pathsRaster)[1])
  
  # pathsRaster <- dir(paste(pathDatos, 'LST_Night_Combinada_Clim_mean', sep=''), pattern = '*.tif$', full.names = T)
  # titulos <- paste('Climatología Media de LST - Día ', 1:length(pathsRaster), ' del año', sep='')
  # pathsRaster <- dir(paste(pathDatos, 'LST_Night_Combinada_Clim_median', sep=''), pattern = '*.tif$', full.names = T)
  # titulos <- paste('Climatología Mediana de LST - Día ', 1:length(pathsRaster), ' del año', sep='')
  # pathsRaster <- dir(paste(pathDatos, 'LST_Night_Combinada_Clim_sd', sep=''), pattern = '*.tif$', full.names = T)
  # titulos <- paste('Climatología Desviación Estándar de LST - Día ', 1:length(pathsRaster), ' del año', sep='')

  # carpetaSalida <- paste(dirname(pathsRaster[1]), '/plots/', sep='')
  # shpBase = shpMask$shp
  # escala <- crearEscalaEquiespaciadaMultiRasters(pathsRaster, nIntervalos=10, nDigitos=1, colores=NULL, brewerPal='Spectral', continuo=T, usarMuestreo=length(pathsRaster) > 1000, nMuestrasPorRaster=25000, lowTrimPU=0.5 / 100, highTrimPU=0.05 / 100)
  # widthPx=1080
  # heightPx=1080
  # alturaEscalaContinua = unit(x=0.7, units = 'in')
  # replot=F
  
  if (is.null(escala)) { escala <- crearEscalaEquiespaciadaMultiRasters(pathsRasters = pathsRaster) }

  dir.create(carpetaSalida, showWarnings = F, recursive = T)
  
  nCoresAUsar <- detectCores(T, T)
  if (nCoresAUsar > 1) {
    cl <- makeCluster(getOption("cl.cores", nCoresAUsar))
    if (exists(x = 'setMKLthreads')) { clusterEvalQ(cl = cl, expr = setMKLthreads(1)) }
    clusterEvalQ(cl = cl, expr = {
      require(rgdal)
      source(paste0(script.dir.funcionesAuxiliares, '/mapearEx.r'))
      source(paste0(script.dir.funcionesAuxiliares, '../pathUtils/pathUtils.r'))
    })
    parSapplyLB(cl=cl, X=1:length(pathsRaster), FUN=plotRasterI, pathsRaster=pathsRaster, shpBase=shpBase, escala=escala, 
                carpetaSalida=carpetaSalida, titulos=titulos, replot=replot, widthPx=widthPx, heightPx=heightPx, 
                alturaEscalaContinua = alturaEscalaContinua)
    stopCluster(cl)
  } else {
    sapply(X=1:length(pathsRaster), FUN=plotRasterI, pathsRaster=pathsRaster, shpBase=shpBase, escala=escala, 
           carpetaSalida=carpetaSalida, titulos=titulos,replot=replot, widthPx=widthPx, heightPx=heightPx, 
           alturaEscalaContinua = alturaEscalaContinua)
  }
}

plotMultiRasters <- function(pathsRasters, carpetaSalida='/plots/', shpBase=NULL, escalas=NULL, replot=F) {
  # Hace gráficos independientes para cada columna en pathsRasters por paso temporal 
  
  # pathsRasters <- pathsRegresores
  # pathsRasters <- unique(pathsRegresores[, 'LST_Night_Combinada_Clim', drop=F])
  # pathsRasters <- pathsRegresores[,1:5]
  # shpBase <- shpMask$shp
  # carpetaSalida='/plots/'
  # escalas=rep(NULL, ncol(pathsRasters))
  # replot=F
  
  if (is.null(escalas)) {
    escala <- crearEscalaEquiespaciadaMultiRasters(pathsRasters)
    escalas <- vector(mode = "list", length = ncol(pathsRasters))
    for (i in 1:length(escalas)) escalas[[i]] <- escala
  }
  
  #i<-1
  #pathsRaster = unique(pathsRasters[, i])
  #carpetaSalida <- paste(dirname(pathsRasters[1, i]), carpetaSalida, sep='')
  #escala=escalas[i]
  
  for (i in 1:ncol(pathsRasters)) {
    print(i)
    cSalida <- paste(dirname(pathsRasters[1, i]), carpetaSalida, sep='')
    plotRasters(pathsRaster = unique(pathsRasters[, i]), shpBase = shpBase, carpetaSalida = cSalida, escala=escalas[[i]], replot=replot)
  }
}

plotMultiRastersEnPanelesI <- function(i, pathsRasters, fechasRasters, shpBase=NULL, escalas=NULL, carpetaSalida='plots', 
                                       postFijoNomArchSalida='', replot=F, nCols, widthPx=1920, heightPx=1080 * widthPx / 1920, 
                                       alturaEscalaContinua = unit(x=2 * heightPx / 1080, units = 'in')) {
  print(i)
  
  # i <- 1
  archivoSalida <- paste(carpetaSalida, fechasRasters[i], postFijoNomArchSalida, '.png', sep='')
  iNoNA <- which(!is.na(pathsRasters[i,]))
  if (length(iNoNA) > 0 && (replot || !file.exists(archivoSalida) || file.info(archivoSalida)$size == 0)) {
    if (is.null(escalas)) {
      escala <- crearEscalaEquiespaciadaMultiRasters(pathsRasters[i, , drop=F])
      escalas <- vector("list", length = ncol(pathsRasters))
      for (j in 1:ncol(pathsRasters)) escalas[[j]] <- escala
    }
    
    gs <- vector(mode = "list", length = ncol(pathsRasters))
    j <- 1
    #j <- j + 1
    for (j in 1:ncol(pathsRasters)) {
      print(j)
      if (!is.na(pathsRasters[i, j])) { grilla <- readGDAL(pathsRasters[i, j], silent = T)
      } else { 
        grilla <- readGDAL(pathsRasters[i, iNoNA[1]], silent = T) 
        grilla@data[,1] <- NA
      }
      
      if (!is.null(shpBase) && !identicalCRS(grilla, shpBase))  {
        shpBase <- spTransform(shpBase, grilla@proj4string)
      }
      if (!is.null(escalas[[j]]) && !all(is.na(grilla@data[, 1]))) {
        escalas[[j]] <- ajustarExtremosEscala(escala = escalas[[j]], datos = grilla@data[, 1], redondear = F)
      } 
      
      gs[[j]] <- mapearGrillaGGPlot(
        grilla = grilla, shpBase = shpBase, escala = escalas[[j]], dibujar=F, 
        titulo = paste(as.character(fechasRasters[i]), ': ', colnames(pathsRasters)[j], sep=''), 
        alturaEscalaContinua = alturaEscalaContinua)
        #, dibujarPuntosObservaciones = T, coordsObservaciones = coordsObservaciones, tamanioFuentePuntos = 4)
    }
    oldSciPen <- getOption("scipen")
    options(scipen=15)
    png(archivoSalida, width = widthPx, height = heightPx)
    tryCatch(expr = print(multiplot(plotlist=gs, cols=nCols)), finally = dev.off())
    options(scipen = oldSciPen)
  }

  return(NULL)
}

permutacionColumnasParaGraficarPorFilas <- function(nGraficos, nColsPlot=2) {
  nFilas <- ceiling(nGraficos / nColsPlot)
  # 1, 2, 3, 4, 5, 6 en una matriz 2x3
  # 1 4 2 5 3 6
  # 1, 2, 3, 4, 5, 6 en una matriz 3x2
  # 1 3 5 2 4 6
  # 1, 2, 3, 4 en una matriz 2x2
  # 1 3 2 4
  
  iCols <- seq.int(from = 1, to = nColsPlot, by = 1)
  permutacionPorFilas <- as.integer(t(sapply(X = seq.int(from = 0, to = nFilas * nColsPlot - 1, by = nColsPlot), FUN = function(x) x + iCols)))
  permutacionPorFilas <- permutacionPorFilas[permutacionPorFilas <= nGraficos]
  return(permutacionPorFilas)
}

reordenarColumnasParaGraficarPorFilas <- function(paths, nColsPlot=2) {
  permutacionPorFilas <- permutacionColumnasParaGraficarPorFilas(nGraficos = ncol(paths), nColsPlot = nColsPlot)
  return(paths[, permutacionPorFilas, drop=F])
}

plotMultiRastersEnPaneles <- function(pathsRasters, fechasRasters, carpetaSalida='plots/', postFijoNomArchSalida='', shpBase=NULL, 
                                      escalas=vector("list", ncol(pathsRasters)), nCols=2, replot=T, byRow=T, nCoresAUsar=0, widthPx=1920, 
                                      heightPx=1080 * widthPx / 1920, alturaEscalaContinua=unit(x=2 * heightPx / 1080, units = 'in')) {
  # Hace un gráfico por paso temporal con un panel para cada columna en pathsRasters
  
  #pathsRasters <- matrix(c('Resultados/2-GrilladoYCV/K/2014_01_12_K.tif', 'Resultados/2-GrilladoYCV/GRK-MOD11A1_LST_Night_FR+MYD11A1_LST_Night_FR+LST_Night_Combinada_Clim_mean+x+y/2014_01_12_GRK-MOD11A1_LST_Night_FR+MYD11A1_LST_Night_FR+LST_Night_Combinada_Clim_mean+x+y.tif'), ncol = 2)
  #fechasRasters <- fechasObservaciones[which(fechasObservaciones == as.POSIXct('2014-01-12', tz=tz(fechasObservaciones[1])))]
  #nCols <- 2
  #colnames(pathsRasters) <- c('K', 'RK')
  #carpetaSalida='F:/Tesis/'
  #widthPx=630 * nCols
  #heightPx=630
  
  # pathsRasters <- pathsRegresores[, c("MOD11A1_LST_Night_filtrado", "MOD11A1_LST_Night_FR", "MOD11A1_LST_Night_FRv2","MYD11A1_LST_Night_filtrado", "MYD11A1_LST_Night_FR", "MYD11A1_LST_Night_FRv2")]
  # pathsRasters <- pathsRegresores[, c("MOD11A1_LST_Night_filtrado", "MOD11A1_LST_Night_FR", "MOD11A1_LST_Night_FRv2","MYD11A1_LST_Night_filtrado", "MYD11A1_LST_Night_FR", "MYD11A1_LST_Night_FRv2")]
  # fechasRasters <- fechasObservaciones
  #carpetaSalida <- 'Resultados/RellenosEspaciales/'
  #shpBase <- shpMask$shp
  #escalas <- NULL
  #nCols <- 3
  #replot=F
  #byRow=T
  #nCoresAUsar=0
  #widthPx=1920
  #heightPx=1080 * widthPx / 1920
  #alturaEscalaContinua=unit(x=2 * heightPx / 1080, units = 'in')
  #postFijoNomArchSalida=''
  
  #pathsRasters <- pathsRegresores[,c(7, 12, 1)]
  #carpetaSalida <- 'Resultados/LST_MODIS_Combinada/'
  # pathsRasters <- pathsRegresores[,c(1:5, 20)]
  #carpetaSalida <- 'Resultados/1-Exploracion/LST_MODIS/'
  #nCols <- 3
  
  #LST_Night_Combinada_Rellenado <- expandirPathsRegresor(pathsRegresor = dir(path = paste(pathDatos, 'LST_Night_Combinada_Rellenado_RK-LST_Night_Combinada_Clim/', sep=''), pattern='*.tif$', full.names = T), 
  #                                                       fechasRegresando = fechasObservaciones, expandir=F)
  #pathsRasters <- cbind(pathsRegresores[,1, drop=F], LST_Night_Combinada_Rellenado, pathsRegresores[,c(5, 20)])
  #nCols <- 2
  
  if (byRow) {
    permutacionPorFilas <- permutacionColumnasParaGraficarPorFilas(nGraficos = ncol(pathsRasters), nColsPlot = nCols)
    pathsRasters <- pathsRasters[, permutacionPorFilas, drop=F]
    escalas <- escalas[permutacionPorFilas]
    rm(permutacionPorFilas)
  }
  
  #carpetaSalida <- 'resultados/1-Exploracion/LST_MODIS_Rellenada/'
  #fechasRasters <- fechasObservaciones
  #shpBase <- shpMask$shp
  #escalas=rep(NULL, ncol(pathsRasters))
  #replot=F
  
  if (length(escalas) > 0) {
    for (i in 1:length(escalas)) {
      if (is.null(escalas[[i]])) {
        muestras <- muestrearValores(pathsRaster = pathsRasters[, i], nMuestras = 25000)
        #muestras <- muestrearValores(pathsRaster = pathsRaster, nMuestras = 25000)
        muestras <- muestras[!is.na(muestras)]
        sMuestras <- sort(muestras)
        lowTrim <- 1
        highTrim <- 0
        lowTrim <- length(sMuestras) * 0.5 / 100
        highTrim <- length(sMuestras) * 0.05 / 100
        #plot(ecdf(sMuestras[lowTrim:(length(sMuestras)-highTrim)]))
        #f <- ecdf(sMuestras[lowTrim:(length(sMuestras)-highTrim)])
        #f(30)
        #f(-3)
        #quantile(sMuestras[lowTrim:(length(sMuestras)-highTrim)], probs=seq(from=0, to=1, length.out = 11))
        #seq(from=sMuestras[iTrim], to = sMuestras[lowTrim:(length(sMuestras)-highTrim)], length.out=11)
        escalas[[i]] <- crearEscalaEquiespaciada(sMuestras[c(lowTrim,(length(sMuestras)-highTrim))], nIntervalos = 10, continuo = T)
        #escala <- crearEscalaEquiespaciadaEnQuantiles(sMuestras[lowTrim:(length(sMuestras)-highTrim)], nIntervalos = 10, continuo = T)
        rm(muestras, sMuestras)
      }
    }
  }
  
  if (F) {
    muestras <- muestrearValores(pathsRaster = pathsRasters, nMuestras = 25000)
    muestras <- muestras[!is.na(muestras)]
    sMuestras <- sort(muestras)
    lowTrim <- 1
    highTrim <- 0
    lowTrim <- length(sMuestras) * 0.5 / 100
    highTrim <- length(sMuestras) * 0.05 / 100
    escalas[[1]] <- crearEscalaEquiespaciada(sMuestras[c(lowTrim,(length(sMuestras)-highTrim))], nIntervalos = 10, continuo = T)
    escalas[[2]] <- escalas[[1]]
    escalas[[3]] <- escalas[[1]]
    #escala <- crearEscalaEquiespaciadaEnQuantiles(sMuestras[lowTrim:(length(sMuestras)-highTrim)], nIntervalos = 11, continuo = T)
    rm(muestras, sMuestras)    
  }
  
  if (!dir.exists(carpetaSalida)) dir.create(carpetaSalida, showWarnings = F)
  
  if (nCoresAUsar <= 0) nCoresAUsar <- min(detectCores(T, T), nrow(pathsRasters))
  if (nCoresAUsar > 1) {
    cl <- makeCluster(getOption("cl.cores", nCoresAUsar))
    if (exists(x = 'setMKLthreads')) { clusterEvalQ(cl = cl, expr = setMKLthreads(1)) }
    clusterExport(cl, varlist = 'script.dir.funcionesAuxiliares')
    clusterEvalQ(cl = cl, expr = {
      require(rgdal)
      require(Rmisc)
      source(paste0(script.dir.funcionesAuxiliares, '/funcionesAuxiliares.r'))
      source(paste0(script.dir.funcionesAuxiliares, '/mapearEx.r'))
    })
    parSapplyLB(cl=cl, X=1:nrow(pathsRasters), FUN=plotMultiRastersEnPanelesI, pathsRasters=pathsRasters, fechasRasters=fechasRasters, shpBase=shpBase, escalas=escalas, 
                carpetaSalida=carpetaSalida, postFijoNomArchSalida=postFijoNomArchSalida, replot=replot, nCols=nCols, widthPx=widthPx, heightPx=heightPx,
                alturaEscalaContinua=alturaEscalaContinua)
    stopCluster(cl)
  } else {
    sapply(X=1:nrow(pathsRasters), FUN=plotMultiRastersEnPanelesI, pathsRasters=pathsRasters, fechasRasters=fechasRasters, shpBase=shpBase, escalas=escalas, 
           carpetaSalida=carpetaSalida, postFijoNomArchSalida=postFijoNomArchSalida, replot=replot, nCols=nCols, widthPx=widthPx, heightPx=heightPx,
           alturaEscalaContinua=alturaEscalaContinua)
  }
}

filtrarRasterI <- function(i, pathsRasters, minVal=NA, maxVal=NA, carpetaSalida='filtrados', zcol=1, gdalOptions = c('COMPRESS=DEFLATE', 'PREDICTOR=2', 'ZLEVEL=9')) {
  grilla <- readGDAL(pathsRasters[i], silent = T)
  if (!is.na(minVal)) {
    ies <- !is.na(grilla@data[,zcol]) & grilla@data[,zcol] < minVal
    grilla@data[ies, zcol] <- NA
  }
  
  if (!is.na(maxVal)) {
    ies <- !is.na(grilla@data[,zcol]) & grilla@data[,zcol] > maxVal
    grilla@data[ies, zcol] <- NA
  }
  
  nomArchSalida <- paste(carpetaSalida, basename(pathsRasters[i]), sep='')
  writeGDAL(dataset = grilla, fname = nomArchSalida, options = gdalOptions)
}

filtrarRasters <- function(pathsRasters, minVal=NA, maxVal=NA, carpetaSalida='filtrados', zcol=1, gdalOptions = c('COMPRESS=DEFLATE', 'PREDICTOR=2', 'ZLEVEL=9')) {
  i <- 1
  pathsRasters <- pathsRegresores[,'LST_Night_Combinada', drop=F]
  minVal <- -8
  maxVal <- 32
  carpetaSalida <- paste(pathDatos, ultimaCarpeta(pathsRasters[i]), '_Filtrada/', sep='') 
  zcol=1
  gdalOptions = c('COMPRESS=DEFLATE', 'PREDICTOR=2', 'ZLEVEL=9')
  nCoresAUsar <- detectCores(T, T)
  
  
  dir.create(carpetaSalida, showWarnings = F, recursive = T)
  if (nCoresAUsar > 1) {
    cl <- makeCluster(getOption("cl.cores", nCoresAUsar))
    if (exists(x = 'setMKLthreads')) { clusterEvalQ(cl = cl, expr = setMKLthreads(1)) }
    clusterEvalQ(cl = cl, expr = {
      require(rgdal)
    })
    parSapplyLB(cl=cl, X=1:length(pathsRasters), FUN=filtrarRasterI, pathsRasters=pathsRasters, minVal=minVal, maxVal=maxVal, carpetaSalida=carpetaSalida, 
                zcol=zcol, gdalOptions=gdalOptions)
    stopCluster(cl)
  } else {
    sapply(X=1:nrow(pathsRasters), FUN=filtrarRasterI, pathsRasters=pathsRasters, minVal=minVal, maxVal=maxVal, carpetaSalida=carpetaSalida, 
           zcol=zcol, gdalOptions=gdalOptions)
  }  
}

crearDFLeonardo <- function() {
  estacs <- spTransform(
    estaciones, CRS(projargs = "+proj=longlat +datum=WGS84", SRS_string = "EPSG:4326"))
  coords <- sp::coordinates(estacs)
  
  # nro estacion | lon | lat | altura | año | mes | día | modis 1 | modis 2 |  temperaturas registradas ese día|
  iEstac <- integer()
  lon <- numeric()
  lat <- numeric()
  altitud <- numeric()
  anio <- integer()
  mes <- integer()
  dia <- integer()
  
  i<-1
  for (i in 1:length(estaciones)) {
    iEstac <- c(iEstac, rep(i, length(fechasObservaciones)))
    lon <- c(lon, rep(coords[i, 'Longitud'], length(fechasObservaciones)))
    lat <- c(lat, rep(coords[i, 'Latitud'], length(fechasObservaciones)))
    altitud <- c(altitud, rep(estacs$Altitud[i], length(fechasObservaciones)))
    anio <- c(anio, year(fechasObservaciones))
    mes <- c(mes, month(fechasObservaciones))
    dia <- c(dia, day(fechasObservaciones))
  }
  
  modis1 <- extraerValoresRegresorSobreSP(objSP = coordsObservaciones, pathsRegresor = pathsRegresores[, 'MOD11A1_LST_Night', drop=F])
  modis2 <- extraerValoresRegresorSobreSP(objSP = coordsObservaciones, pathsRegresor = pathsRegresores[, 'MYD11A1_LST_Night', drop=F])
  
  dataframe <- data.frame(nroEstacion=iEstac, lon=lon, lat=lat, altitud=altitud, anio=anio, mes=mes, dia=dia, 
                          modis1=as.vector(modis1), modis2=as.vector(modis2), tmin=as.vector(valoresObservaciones), stringsAsFactors = F)
  
  write.table(dataframe, file = "D:/Tesis/Datos/dfLeonardo.tsv", append = FALSE, quote = TRUE, sep = "\t",
              eol = "\n", na = "NA", dec = ".", row.names = F, col.names = TRUE, qmethod = c("escape", "double"),
              fileEncoding = "")
}


rellenarRegresor_ti <- function(ti=1, pathsRegresor, carpetaSalida, shpMask, metodo='fastTps', zcol=1) {
  #ti <- 6
  #print(ti)

  sp <- readGDAL(pathsRegresor[ti])
  sp <- rellenarSP(sp=sp, shpMask = shpMask, metodo = metodo, zcol=zcol)
  writeGDAL(dataset = sp, fname = paste(carpetaSalida, basename(pathsRegresor[ti]), sep=''), options = c('COMPRESS=DEFLATE', 'PREDICTOR=2', 'ZLEVEL=9'))
}

rellenarRegresores <- function(pathsRegresores, carpetasSalida=paste(apply(X = pathsRegresores, MARGIN = 2, FUN = function(x) { unique(dirname(x[!is.na(x)])) }), '_rellenado/', sep=''),
                               pathShpMascara, proj4stringSHPMascara, metodo='fastTps', zcols=rep(1, ncol(pathsRegresores))) {
  pathsRegresores <- pathsRegresores[, c('LST_Night_Combinada5', 'LST_Night_Combinada9')]
  carpetasSalida <- paste(apply(X = pathsRegresores, MARGIN = 2, FUN = function(x) { unique(dirname(x[!is.na(x)])) }), '_rellenado/', sep='')
  pathShpMascara <- params$pathSHPMapaBase
  metodo='loess'
  zcols=rep(1, ncol(pathsRegresores))
  nCoresAUsar <- 1
  nCoresAUsar <- detectCores(T, T)
  i<-1
  for (i in 1:ncol(pathsRegresores)) {
    dir.create(carpetasSalida[i], showWarnings = FALSE)
    sp <- readGDAL(pathsRegresores[1, i])
    shpMask <- cargarSHPYObtenerMascaraParaGrilla(pathSHP=pathShpMascara, grilla = sp)
    
    if (nCoresAUsar > 1) {
      cl <- makeCluster(getOption("cl.cores", nCoresAUsar))
      clusterExport(cl, varlist = c('script.dir.funcionesAuxiliares'))
      clusterEvalQ(cl = cl, expr = {
        require('rgdal')
        source(paste0(script.dir.funcionesAuxiliares, 'funcionesAuxiliares.r'))
      })
      if (exists(x = 'setMKLthreads')) { clusterEvalQ(cl = cl, expr = setMKLthreads(1)) }
      parSapplyLB(cl=cl, X=1:length(pathsRegresores[, i]), FUN=rellenarRegresor_ti, 
                  pathsRegresor=pathsRegresores[, i], carpetaSalida=carpetasSalida[i], shpMask=shpMask, metodo=metodo, zcol=zcols[i])
      stopCluster(cl)
    } else {
      sapply(X=1:length(pathsRegresores[, i]), FUN=rellenarRegresor_ti, pathsRegresor=pathsRegresores[, i],
             carpetaSalida=carpetasSalida[i], shpMask=shpMask, metodo=metodo, zcol=zcols[i])
    }
  }
}

coarsenGridEx <- function(object, coarse = 2, offset = sample(c(0:(coarse - 1)), 2, replace = TRUE)) {
  objCor = as(object, "SpatialPoints")
  xss = unique(sp::coordinates(objCor)[, 1])
  yss = unique(sp::coordinates(objCor)[, 2])
  xi = c(1:(length(xss)/coarse))
  xs = sort(xss)[(xi - 1) * coarse + offset[1] + 1]
  yi = c(1:(length(yss)/coarse))
  ys = sort(yss)[(yi - 1) * coarse + offset[2] + 1]
  sel = sp::coordinates(objCor)[, 1] %in% xs & sp::coordinates(objCor)[, 2] %in% ys
  newPoints = objCor[sel, ]
  if ("data" %in% names(getSlots(class(object)))) 
    newPoints = SpatialPointsDataFrame(newPoints, data = object@data[sel, , drop=F])
  gridded(newPoints) = TRUE
  newPoints
}

dividirEnCuadrantes <- function(
    object, nCuadrantesX = 2, nCuadrantesY=2, claseObjetoRetornado=c('SpatialPolygons', 'SpatialGrid')) {
  boundingBox <- bbox(object)
  
  cellsizeX <- (boundingBox[1, 2] - boundingBox[1, 1]) / nCuadrantesX
  cellsizeY <- (boundingBox[2, 2] - boundingBox[2, 1]) / nCuadrantesY
  cellsize <- c(cellsizeX, cellsizeY)
  cellcentre.offset <- c(boundingBox[1, 1], boundingBox[2, 1]) + cellsize / 2
  
  cells.dim <- c(nCuadrantesX, nCuadrantesY)
  newTopology <- GridTopology(cellcentre.offset, cellsize, cells.dim)
  
  if (claseObjetoRetornado[1] != 'SpatialGrid') { 
    return(as.SpatialPolygons.GridTopology(grd=newTopology, proj4string=object@proj4string))
  } else { 
    return(SpatialGrid(newTopology, proj4string = object@proj4string))
  }
}

muestrearEnCuadrantes <- function(sp, size, nCuadrantesX = 2, nCuadrantesY=2, zcol=-1) {
  cuadrantes <- dividirEnCuadrantes(sp, nCuadrantesX=nCuadrantesX, nCuadrantesY=nCuadrantesY)
  
  # Obtengo los índices de sp sobre cada cuadrante en una lista con un array de índices por cada cuadrante
  iES <- over(cuadrantes, geometry(sp), returnList = T)
  nMuestrasPorCuadrante <- trunc((size / length(iES)))
  
  if (zcol > 0) {
    # Saco los píxeles que tienen valores NA
    iNAs <- which(is.na(sp@data[,zcol]))
    iES <- lapply(iES, FUN= function(x) { base::setdiff(x, iNAs) })
  }
  
  # over puede retornar pixeles en dos polígonos distintos (si el pixel cae en el borde)
  # con esto fuerzo a que cada pixel pertenezca solo al primer polígono en el que aparece
  for (i in seq.int(from = 1, to = length(iES) -1, by = 1)) {
    for (j in seq.int(from = i+1, to = length(iES), by = 1)) {
      iES[[j]] <- base::setdiff(x = iES[[j]], y = iES[[i]])
    }
  }
  
  # iTodos son índices de los cuadrantes que tienen menos píxeles que nMuestrasPorCuadrante.
  # En esos cuadrantes elijo todos los píxeles, y reparto los píxeles que les sobraron (nMuestrasPorCuadrante - length(cuadrante))
  # entre los demás cuadrantes.
  # Una vez que el cuadrante está completo lo saco de iTodos y sigo la ronda
  muestras <- integer(0)
  iTodos <- which(sapply(iES, function(x) { length(x) <= nMuestrasPorCuadrante }))
  while (length(iTodos) > 0) {
    for (i in length(iTodos):1) {
      muestras <- c(muestras, iES[[iTodos[i]]])
      iES <- iES[-iTodos[i]]
    }
    
    if (length(iES) > 0) {
      nMuestrasPorCuadrante <- trunc(((size - length(muestras)) / length(iES)))  
      iTodos <- which(sapply(iES, function(x) { length(x) < nMuestrasPorCuadrante }))
    } else {
      iTodos <- integer(0)
    }
  }
  
  for (i in 1:length(iES)) {
    muestras <- c(muestras, sample(iES[[i]], size = nMuestrasPorCuadrante, replace = F))
  }
  return(muestras)
}

muestrearEnCuadrantesYECDF <- function(sp, size, nCuadrantesX = 4, nCuadrantesY=round(nCuadrantesX * diff(bbox(sp)[2,]) / diff(bbox(sp)[1,])), 
                                       nCuadrantesZ=trunc(size/(nCuadrantesX*nCuadrantesY)), zcol=1, iEsAExcluir=integer(0),
                                       maxNMuestrasPorCuadrante = round(1.5 * size / (nCuadrantesX * nCuadrantesY))) {
  #sp <- readGDAL(pathsRegresores[17, 1])
  #size = 2000
  # Esta versión reparte las muestras en nCuadrantesX * nCuadrantesY en el espacio y además entre una agrupación de nCuadrantesZ conjuntos
  # de los cuantiles de la variable. El objetivo es obtener un muestreo de toda la distribución de los valores en cada cuadrante
  cuadrantes <- dividirEnCuadrantes(object = sp, nCuadrantesX=nCuadrantesX, nCuadrantesY=nCuadrantesY)
  
  # Obtengo los índices de sp sobre cada cuadrante en una lista con un array de índices por cada cuadrante
  iES <- over(x = cuadrantes, y = geometry(sp), returnList = T)
  
  # Asigno la misma cantidad de muestras a cada cuadrante, si me quedan muestras sin asignar las reparto
  # aleatoriamente, acotado por maxNMuestrasPorCuadrante
  nMuestrasPorCuadrante <- min(trunc((size / length(iES))), maxNMuestrasPorCuadrante)
  nMuestrasCuadrantes <- rep(nMuestrasPorCuadrante, length(iES))
  nMuestrasTotales <- nMuestrasPorCuadrante * length(iES)
  if (nMuestrasTotales < size & nMuestrasPorCuadrante < maxNMuestrasPorCuadrante) {
    #iESConMasMuestras <- c(3, 16)
    iESConMasMuestras <- sample(seq_along(iES), size = size - nMuestrasTotales, replace = FALSE)
    nMuestrasCuadrantes[iESConMasMuestras] <- nMuestrasCuadrantes[iESConMasMuestras] + 1
  }
  
  # Saco los píxeles que tienen valores NA
  iNAsOAExcluir <- c(which(is.na(sp@data[, zcol])), iEsAExcluir)
  iES <- lapply(iES, FUN= function(x) { base::setdiff(x, iNAsOAExcluir) })
  
  # over puede retornar pixeles en dos polígonos distintos (si el pixel cae en el borde)
  # con esto fuerzo a que cada pixel pertenezca solo al primer polígono en el que aparece
  for (i in seq.int(from = 1, to = length(iES) - 1, by = 1)) {
    for (j in seq.int(from = i + 1, to = length(iES), by = 1)) {
      iES[[j]] <- base::setdiff(x = iES[[j]], y = iES[[i]])
    }
  }
  # La cantidad de pixeles en todos los grupos de iES debe ser igual a la cantidad de pixeles totales
  # sum(sapply(iES, FUN=length)) == length(sp)
  
  # iTodos son índices de los cuadrantes que tienen menos píxeles que nMuestrasPorCuadrante.
  # En esos cuadrantes elijo todos los píxeles, y reparto los píxeles que les sobraron (nMuestrasPorCuadrante - length(cuadrante))
  # entre los demás cuadrantes.
  # Una vez que el cuadrante está completo lo saco de iTodos y sigo la ronda
  muestras <- integer(0)
  iTodos <- which(sapply(X = seq_along(iES), FUN = 
                         function(x, iES, nMuestrasCuadrantes) { 
                           return(length(iES[[x]]) <= nMuestrasCuadrantes[[x]]) }, 
                         iES=iES, nMuestrasCuadrantes=nMuestrasCuadrantes))
  while (length(iTodos) > 0) {
    i <- 1
    muestras <- c(muestras, unlist(iES[iTodos], use.names = FALSE))
    iES <- iES[-iTodos]
    
    if (length(iES) > 0) {
      nMuestrasPorCuadrante <- min(trunc((size - length(muestras)) / length(iES)), maxNMuestrasPorCuadrante)
      nMuestrasCuadrantes <- rep(nMuestrasPorCuadrante, length(iES))
      nMuestrasTotales <- nMuestrasPorCuadrante * length(iES) + length(muestras)
      if (nMuestrasTotales < size & nMuestrasPorCuadrante < maxNMuestrasPorCuadrante) {
        iESConMasMuestras <- sample(seq_along(iES), size = size - nMuestrasTotales, replace = FALSE)
        nMuestrasCuadrantes[iESConMasMuestras] <- nMuestrasCuadrantes[iESConMasMuestras] + 1
      }      
      
      iTodos <- which(sapply(X = seq_along(iES), FUN = 
                               function(x, iES, nMuestrasCuadrantes) { 
                                 return(length(iES[[x]]) <= nMuestrasCuadrantes[[x]]) }, 
                             iES=iES, nMuestrasCuadrantes=nMuestrasCuadrantes))
    } else {
      iTodos <- integer(0)
    }
  }
  
  if (length(iES) > 0) {
    # Ordeno los indices de cada cuadrante de acuerdo al valor en el objeto sp para obtener las ECDFs
    iESCuadrantesOrdenados <- lapply(X = iES, FUN = function(x, sp) { x[order(sp@data[x, zcol])] }, sp=sp)
    i <- 1
    i <- i + 1
    for (i in seq_along(iESCuadrantesOrdenados)) {
      # Ordeno los índices del cuadrante según la distribución de z (sp@data[, zcol])
      # Reparto las muestras proporcionalmente entre los cuantiles de z agrupándolos en nCuadrantesZ
      # P.ej: si nCuadrantesZ es 3 reparto las muestras del intervalo en 1 / 3 por cada tercil,
      #       si nCuadrantesZ es 4 reparto las muestras del intervalo en 1 / 4 por cada cuartil, etc
      iESCuadranteOrdenados <- iESCuadrantesOrdenados[[i]]

      # Al menos tengo que tener una muestra por cuadranteZ si no tengo suficiente, disminuyo los 
      # cuadrantesZ para el i-ésimo cuadranteXY
      if (length(iESCuadranteOrdenados) > nCuadrantesZ) { nCuadrantesZAux <- nCuadrantesZ 
      } else { nCuadrantesZAux <- length(iESCuadranteOrdenados) }
            
      lengthIntervalosZ <- length(iESCuadranteOrdenados) / nCuadrantesZAux
      nMuestrasPorIntervaloZ <- trunc(nMuestrasCuadrantes[[i]] / nCuadrantesZAux)
      nMuestrasIntervalosZ <- rep(nMuestrasPorIntervaloZ, nCuadrantesZAux)
      nMuestrasTotales <- nMuestrasPorIntervaloZ * nCuadrantesZAux
      
      iIntervalosZ <- lapply(X = seq_along(nMuestrasIntervalosZ), FUN = function(x) {
        iIniIntervaloZ <- round((x - 1) * lengthIntervalosZ + 1)
        iFinIntervaloZ <- round(x * lengthIntervalosZ)
        return(seq.int(from = iIniIntervaloZ, to = iFinIntervaloZ, by = 1))
      })
      
      if (nMuestrasTotales < nMuestrasCuadrantes[[i]]) {
        iESQuePuedenTenerMasMuestras <- which(sapply(iIntervalosZ, FUN = length) > nMuestrasIntervalosZ)
        if (nMuestrasCuadrantes[[i]] - nMuestrasTotales < length(iESQuePuedenTenerMasMuestras)) {
          iESConMasMuestras <- sample(iESQuePuedenTenerMasMuestras, size = nMuestrasCuadrantes[[i]] - nMuestrasTotales, replace = FALSE)
        } else { iESConMasMuestras <- iESQuePuedenTenerMasMuestras }
        nMuestrasIntervalosZ[iESConMasMuestras] <- nMuestrasIntervalosZ[iESConMasMuestras] + 1
      }
      
      j <- 1
      j <- j + 1
      for (j in seq_along(nMuestrasIntervalosZ)) {
        if (nMuestrasIntervalosZ[j] > 0) {
          #iIniIntervaloZ <- round((j - 1) * lengthIntervalosZ + 1)
          #iFinIntervaloZ <- round(j * lengthIntervalosZ)
          #iIntervaloZ <- seq.int(from = iIniIntervaloZ, to = iFinIntervaloZ, by = 1)
          #print(paste(paste(i, '.', j, '= ', sep=''),  paste(iIntervaloZ, collapse = ','), sep = ''))
          
          muestras <- c(muestras, sample(iESCuadranteOrdenados[iIntervalosZ[[j]]], size = nMuestrasIntervalosZ[j], replace = F))
        }
      }
    }
  }

  if (F) {
    # Gráficos para mostrar los resultados
    if (gridded(sp)) {
      g1 <- mapearGrillaGGPlot(grilla = sp, continuo = T, dibujar = F, titulo = 'LST MODIS: 2002-07-20')  
    } else {
      g1 <- mapearPuntosGGPlot(sp, tamaniosPuntos = 2, continuo = T, dibujar = F, titulo = )
    }
    
    g2 <- mapearPuntosGGPlot(sp[muestras,], tamaniosPuntos = 2, continuo = T, dibujar = F, titulo = 'Muestreo LST MODIS: 2002-07-20. n=2000')
    png('Resultados/Ejemplos/MuestreoPorCuadrantesYECDF_Mapas.png', width = 1920, height = 1017)
    tryCatch(expr = print(multiplot(plotlist = list(g1, g2), cols=2)), finally = dev.off())

    auxSP <- SpatialPoints(coords = sp::coordinates(sp)[muestras,], proj4string = sp@proj4string) 
    auxSP <- SpatialPointsDataFrame(coords = auxSP, data = data.frame(value=over(auxSP, sp))) 
    
    i1 <- over(cuadrantes, sp, returnList = T)
    i2 <- over(cuadrantes, auxSP, returnList = T)
    
    nombresCuadrantes <- apply(expand.grid(0:(nCuadrantesX-1), 0:(nCuadrantesY-1)), 1, function(x) paste(x,collapse="-"))
    
    png('Resultados/Ejemplos/MuestreoPorCuadrantesYECDF_ECDF_Cuadrantes.png', width = 1920, height = 1017)
    tryCatch(expr = {  
      par(mfrow=c(nCuadrantesX,nCuadrantesY))
      for (i in 1:length(i1)) {
        plot(ecdf(i1[[i]][,1]), main = paste('ECDF LST MODIS: 2002-07-20. Cuadrante ', nombresCuadrantes[i], sep=''), xlab='LST [°C]', xlim=range(sp@data[,zcol], na.rm=T))
        lines(ecdf(i2[[i]][,1]), col='blue', lwd=2, pch='')
        legend("topleft", inset=.05, cex = 2, c("Original","Muestra n=2000"), horiz=F, lty=c(1,1,1), lwd=c(2,2), col=c("black","blue"), bg="grey96")
      }
    }, finally = { 
      par(mfrow=c(1,1))
      dev.off() })
    
    png('Resultados/Ejemplos/MuestreoPorCuadrantesYECDF_ECDF.png', width = 1920, height = 1017)
    tryCatch(expr = {  
      plot(ecdf(sp@data[,zcol]), main = 'ECDF LST MODIS: 2002-07-20', xlab='LST [°C]', xlim=range(sp@data[,zcol], na.rm=T))
      lines(ecdf(sp@data[muestras,zcol]), col='blue', lwd=2, pch='')
      legend("topleft", inset=.05, cex = 2, c("Original","Muestra n=2000"), horiz=F, lty=c(1,1), lwd=c(2,2), col=c("black","blue"), bg="grey96")
    }, finally = dev.off())
  }
  
  return(muestras)
}

muestrearEnCuadrantesYECDF_V2 <- function(sp, size, nCuadrantesX = 2, nCuadrantesY=2, nCuadrantesZ=10, zcol=1) {
  # sp <- readGDAL(pathsRegresores[17, 1])
  # size = 2000
  # NO-USAR.
  # Otro intento, la idea era asignar una clase a cada dato según el cuadrante y el intervalo de la ECDF a la que correspondiera y repartir las 
  # muestras entre esas clases pero no tuvo exito. La ECDF general sigue sesgada y además se sesgan las ECDFs de cada cuadrante
  cuadrantes <- dividirEnCuadrantes(sp, nCuadrantesX=nCuadrantesX, nCuadrantesY=nCuadrantesY)
  
  # Obtengo los índices de sp sobre cada cuadrante en una lista con un array de índices por cada cuadrante
  iES <- over(cuadrantes, geometry(sp), returnList = T)
  
  # Saco los píxeles que tienen valores NA
  iNAs <- which(is.na(sp@data[, zcol]))
  iES <- lapply(iES, FUN= function(x) { base::setdiff(x, iNAs) })
  
  # over puede retornar pixeles en dos polígonos distintos (si el pixel cae en el borde)
  # con esto fuerzo a que cada pixel pertenezca solo al primer polígono en el que aparece
  for (i in seq.int(from = 1, to = length(iES) -1, by = 1)) {
    for (j in seq.int(from = i+1, to = length(iES), by = 1)) {
      iES[[j]] <- base::setdiff(x = iES[[j]], y = iES[[i]])
    }
  }
  
  todosLosIES <- unlist(iES, use.names = F)
  
  iCuadrante <- integer(0)
  for (i in 1:length(iES)) {
    iCuadrante <- c(iCuadrante, rep(i, length(iES[[i]])))
  }
  
  nMuestrasPorIntervaloZ <- length(todosLosIES) / nCuadrantesZ
  lims <- matrix(NA, nrow=nCuadrantesZ, ncol=2)
  lims[, 1] <- 1 + (0:(nCuadrantesZ-1)) * nMuestrasPorIntervaloZ
  lims[, 2] <- 1:nCuadrantesZ * nMuestrasPorIntervaloZ
  lims <- round(lims)
  
  iOrden <- order(sp@data[todosLosIES, zcol])
  plot(x=1:length(iOrden), y=sp@data[todosLosIES[iOrden], zcol])
  iCuadranteZ <- apply(X = lims, MARGIN =1, FUN = function(x) { x[1] <= iOrden & iOrden <= x[2] })
  iCuadranteZ <- apply(iCuadranteZ, MARGIN = 1, FUN = which)
  
  clase <- sprintf('%.2d_%.2d', iCuadrante, iCuadranteZ)
  clases <- sort(unique(clase))
  
  iESPorClase <- vector(mode = "list", length = length(clases))
  i <- 5
  for (i in 1:length(clases)) {
    iESPorClase[[i]] <- todosLosIES[clase == clases[i]]
    plot(x=1:length(iESPorClase[[i]]), y=sp@data[iESPorClase[[i]], zcol])
  }
  
  #all(todosLosIES %in% unlist(iESPorClase))
  #length(todosLosIES) == length(unlist(iESPorClase))
  
  iES <- iESPorClase
  nMuestrasPorClase <- size / length(clases)
  
  cbind(clases, sapply(iES, length))
  
  # iTodos son índices de los cuadrantes que tienen menos píxeles que nMuestrasPorCuadrante.
  # En esos cuadrantes elijo todos los píxeles, y reparto los píxeles que les sobraron (nMuestrasPorCuadrante - length(cuadrante))
  # entre los demás cuadrantes.
  # Una vez que el cuadrante está completo lo saco de iTodos y sigo la ronda
  muestras <- integer(0)
  iTodos <- which(sapply(iESPorClase, function(x) { length(x) < nMuestrasPorClase }))
  while (length(iTodos) > 0) {
    for (i in length(iTodos):1) {
      muestras <- c(muestras, iES[[iTodos[i]]])
      iES <- iES[-iTodos[i]]
    }
    
    if (length(iES) > 0) {
      nMuestrasPorClase <- trunc(((size - length(muestras)) / length(iES)))  
      iTodos <- which(sapply(iES, function(x) { length(x) < nMuestrasPorClase }))
    } else {
      iTodos <- integer(0)
    }
  }
  
  for (i in 1:length(iES)) {
    muestras <- c(muestras, sample(iES[[i]], size = nMuestrasPorClase, replace = F))
  }
  
  
  return(muestras)
}

splitP4str <- function (p4str) {
  # split proj4string into tokens and remove white spaces
  sp4str <- strsplit(p4str, split = ' *[+]', fixed=F)[[1]]
  sp4str <- sp4str[sp4str != '']
  # split tokens into name=value pairs
  sp4str <- strsplit(sp4str, split = ' *= *', fixed = F)
  # format as matrix
  mp4str <- matrix(data='', nrow = length(sp4str), ncol=2)
  for (i in 1:length(sp4str)) {
    for (j in 1:length(sp4str[[i]])) {
      mp4str[i, j] <- sp4str[[i]][j]
    }
  }
  return(mp4str)
}

sameProjection <- function(p4str1, p4str2) {
  sp4str1 <- splitP4str(p4str1)
  sp4str2 <- splitP4str(p4str2)
  
  if (nrow(sp4str1) == nrow(sp4str2)) {
    sp4str1 <- sp4str1[sp4str1[,1] != 'units', ]
    sp4str2 <- sp4str2[sp4str2[,1] != 'units', ]
    iMatch <- match(sp4str1[, 1], sp4str2[, 1])
    return (!any(is.na(iMatch)) && sp4str1[iMatch, 2] == sp4str2[, 2])
  } else return (FALSE)
}

comprimirGeoTiffs <- function(carpeta, nivelCompresion=9, predictor=2) {
  pathsGeoTiffs <- dir(path = carpeta, pattern = '*.tif$', all.files = T, full.names = T, recursive = T, ignore.case = T)
  pathsGeoTiffs <- pathsGeoTiffs[substr(x =pathsGeoTiffs, start = nchar(pathsGeoTiffs) - 2, stop = nchar(pathsGeoTiffs)) == 'tif']
  pathsGeoTiffs <- pathsGeoTiffs[!grepl(pattern = 'compressed', x = pathsGeoTiffs, fixed = T)]
  
  for (i in 1:length(pathsGeoTiffs)) {
    fsal <- agregarCarpetaAlFinal(pathsGeoTiffs[i], 'compressed')
    if (!file.exists(fsal)) {
      gt <- readGDAL(fname = pathsGeoTiffs[i])
      dir.create(dirname(fsal), showWarnings = F)
      writeGDAL(gt, fsal, options = c('COMPRESS=DEFLATE', 'PREDICTOR=2', 'ZLEVEL=9'))
    }
  }
}

contarNoNulosPorCuadrantesEnObjSP <- function(
    objSP, objCuadrantes = NULL, nCuadrantesX=2, nCuadrantesY=nCuadrantesX, zcol=1, shpMask=NULL) {
  if (is.null(objCuadrantes) || !identicalCRS(objSP, objCuadrantes)) {
    objCuadrantes <- dividirEnCuadrantes(
      object = objSP, nCuadrantesX = nCuadrantesX, nCuadrantesY = nCuadrantesY)
  }
  ies <- over(geometry(objSP), objCuadrantes, returnList = F)
  
  if (!is.null(shpMask)) { iNoNAs <- !is.na(objSP@data[, zcol]) & shpMask$mask
  } else { iNoNAs <- !is.na(objSP@data[, zcol]) }
  ies <- ies[iNoNAs]
  
  # x <- 1
  ns <- sapply(1:length(objCuadrantes), FUN = function(x) { sum(ies == x) })
  ns[is.na(ns)] <- 0
  
  return(ns)
}

contarNoNulosPorCuadrantesTi <- function(pathGeoTiff, objCuadrantes = NULL, nCuadrantesX=2, nCuadrantesY=nCuadrantesX, zcol=1, shpMask=NULL) {
  #print(pathGeoTiff)
  #pathGeoTiff <- "datos/satelites/GPM/20200531.tif"
  #pathGeoTiff <- pathsGeoTiffs[1,]
  if (!is.na(pathGeoTiff)) {
    return(contarNoNulosPorCuadrantesEnObjSP(
      objSP = readGDAL(fname = pathGeoTiff, silent = T), objCuadrantes = objCuadrantes, 
      nCuadrantesX=nCuadrantesX, nCuadrantesY=nCuadrantesY, zcol=zcol, shpMask=shpMask))
  } else if (!is.null(objCuadrantes)) { return(rep(NA_integer_, length(objCuadrantes))) 
  } else { return(NA_integer_) }
}

contarNoNulosPorCuadrantes <- function(
    pathsGeoTiffs, nCuadrantesX=2, nCuadrantesY=nCuadrantesX, zcol=1, nCoresAUsar=0, shpMask=NULL, 
    carpetaSalida='Resultados/1-Exploracion/') {
  # pathsGeoTiffs <- pathsRegresores[, 1, drop=F]
  reg <- colnames(pathsGeoTiffs)[1]
  
  if (nCoresAUsar <= 0) { 
    nCoresAUsar <- min(detectCores(all.tests = T, logical = T), length(pathsGeoTiffs))
  } else { 
    nCoresAUsar <- min(nCoresAUsar, length(pathsGeoTiffs)) 
  }
  
  iPrimerNoNA <- which(!is.na(pathsGeoTiffs))[1]
  evaluarConReintentos(rasterI <- readGDAL(fname = pathsGeoTiffs[iPrimerNoNA], silent = T))
  objCuadrantes <- dividirEnCuadrantes(
    object = rasterI, nCuadrantesX = nCuadrantesX, nCuadrantesY = nCuadrantesY)
  
  if (nCoresAUsar > 1) {
    cl <- makeCluster(getOption('cl.cores', nCoresAUsar))
    clusterExport(cl, varlist = c('script.dir.funcionesAuxiliares'))
    if (exists(x = 'setMKLthreads')) { clusterEvalQ(cl = cl, expr = setMKLthreads(1)) }
    clusterEvalQ(cl, {
      require('rgdal')
      source(paste0(script.dir.funcionesAuxiliares, 'funcionesAuxiliares.r'))
    })
    disponiblesPorCuadrantes <- parSapplyLB(
      cl = cl, X=as.vector(pathsGeoTiffs), FUN = contarNoNulosPorCuadrantesTi, 
      objCuadrantes=objCuadrantes, shpMask=shpMask)
    stopCluster(cl)
  } else {
    disponiblesPorCuadrantes <- sapply(
      X=as.vector(pathsGeoTiffs), FUN = contarNoNulosPorCuadrantesTi, objCuadrantes=objCuadrantes, 
      shpMask=shpMask, simplify = T)
  }
  disponiblesPorCuadrantes[is.na(disponiblesPorCuadrantes)] <- 0
  cuadrantes <- apply(expand.grid(0:(nCuadrantesX-1), 0:(nCuadrantesY-1)), 1, function(x) paste(x,collapse="-"))
  
  # head(t(disponiblesPorCuadrantes))
  df <- as.data.frame(t(disponiblesPorCuadrantes), row.names = as.character(1:ncol(disponiblesPorCuadrantes)))
  df$fecha <- rownames(pathsGeoTiffs)
  #head(df)
  
  #colMaxs(as.matrix(df[,1:4]))
  
  if (!is.null(carpetaSalida)) {
    longDF <- reshape(df, idvar = "fecha",  varying = list(1:4), direction = "long", sep='', v.names = 'V')
    colnames(longDF)[2] <- 'Cuadrante'
    longDF$Cuadrante <- as.factor(cuadrantes[longDF$Cuadrante])
    colnames(df)[1:length(cuadrantes)] <- cuadrantes
    
    # yMax <- max(rowSums(t(disponiblesPorCuadrantes)))
    
    seqAnios <- range(year(longDF$fecha))
    seqAnios <- (seqAnios[1]+1):seqAnios[2]
    gs <- vector(mode = "list", length = length(seqAnios))
    
    breaks <- longDF$Cuadrante
    labels <- as.character(breaks)
    colores <- brewer.pal(4, 'Set1')
    
    for (i in seq(along=seqAnios)) {
      dfI <- longDF[year(longDF$fecha) == seqAnios[i],]
      gs[[i]] <- ggplot(data = dfI, aes(x=fecha, y=V, group=Cuadrante)) + 
        geom_area(aes(fill=Cuadrante), position = 'stack') + 
        scale_fill_manual(breaks=breaks, drop=F, labels=labels, values=colores, na.value="gray95") + 
        labs(x = "", y = "")
      ggtitle(paste('Cantidad de Datos Disponible por Cuadrante -', seqAnios[i])) + 
        theme(plot.title = element_text(size = 18, face = "bold", colour = "black", vjust = 1))
    }
    require('gridExtra')
    
    png(filename = paste(carpetaSalida, 'CantDatosDisponibles_', reg, '.png', sep=''), width = 1920, height = 1050)
    tryCatch(expr = do.call("grid.arrange", c(gs, list(top = paste("Cantidad de Datos Disponibles en ", reg, sep = ''),
                                                       layout_matrix = matrix(1:length(gs), ncol=3, byrow=TRUE)))), 
             finally = dev.off())
  }
  return(df)
}

extraerRangoSPDataFrameI <- function(i, pathsSPDataFrames, zcol=1, na.rm=T) {
  objSP <- readGDAL(fname = pathsSPDataFrames[i])
  return(range(objSP@data, na.rm=na.rm))
}

extraerRangoSPDataFrames <- function(pathsSPDataFrames, zcol=1, na.rm=T) {
  nCoresAUsar <- min(detectCores(T, T), length(pathsSPDataFrames))
  
  if (nCoresAUsar > 1) {
    cl <- makeCluster(getOption('cl.cores', nCoresAUsar))
    if (exists(x = 'setMKLthreads')) { clusterEvalQ(cl = cl, expr = setMKLthreads(1)) }
    clusterEvalQ(cl = cl, expr = require(rgdal))
    res <- parSapplyLB(cl=cl, X=1:length(pathsSPDataFrames), FUN=extraerRangoSPDataFrameI, pathsSPDataFrames=pathsSPDataFrames, zcol=zcol, na.rm=na.rm)
    stopCluster(cl)
  } else {
    res <- sapply(X=1:length(pathsSPDataFrames), FUN=extraerRangoSPDataFrameI, pathsSPDataFrames=pathsSPDataFrames, zcol=zcol, na.rm=na.rm)
  }
  return(res)
}

resumenSPDataFrame <- function(x, na.rm=T) {
  cuantiles <- quantile(x, probs=c(0, 0.005, 0.5, 0.999, 1), na.rm = na.rm)
  res <- c(cuantiles, mean(x, na.rm=na.rm), sd(x, na.rm=na.rm), diff(range(x, na.rm=na.rm)), sum(!is.na(x))/length(x))
  res[is.nan(res) | is.infinite(res)] <- NA
  names(res) <- c('Min', '0.5%', 'Mediana', '99.9%', 'Max', 'media', 'sd', 'amplitud', 'cobertura')
  return(res)
}

sapplySPDataFrameI <- function(i, pathsSPDataFrames, funcion=base::mean, zcol=1, na.rm=T) {
  objSP <- readGDAL(pathsSPDataFrames[i])
  #x <- objSP@data[,zcol]
  return(eval(funcion(objSP@data[,zcol], na.rm = na.rm)))
}

sapplySPDataFrames <- function(pathsSPDataFrames, funcion=resumenSPDataFrame, zcol=1, na.rm=T) {
  # pathsSPDataFrames <- pathsRegresor
  nCoresAUsar <- min(detectCores(T, T), length(pathsSPDataFrames))
  
  if (nCoresAUsar > 1) {
    cl <- makeCluster(getOption('cl.cores', nCoresAUsar))
    if (exists(x = 'setMKLthreads')) { clusterEvalQ(cl = cl, expr = setMKLthreads(1)) }
    clusterEvalQ(cl = cl, expr = require(rgdal))
    res <- parSapplyLB(cl=cl, X=1:length(pathsSPDataFrames), FUN=sapplySPDataFrameI, pathsSPDataFrames=pathsSPDataFrames, funcion=funcion, zcol=zcol, na.rm=na.rm)
    stopCluster(cl)
  } else {
    res <- sapply(X=1:length(pathsSPDataFrames), FUN=sapplySPDataFrameI, pathsSPDataFrames=pathsSPDataFrames, funcion=funcion, zcol=zcol, na.rm=na.rm)
  }
  
  return(res)
}

filtrarSPDataFrameIConClimatologia <- function(i, pathsRasters, clasesClimatologia, pathsClimatologiaMin, pathsClimatologiaMax, minValAbs=NA, maxValAbs=NA, 
                                               carpetaSalida='filtrados', zcol=1, gdalOptions = c('COMPRESS=DEFLATE', 'PREDICTOR=2', 'ZLEVEL=9'), reRun=F) {
  #i <- 11
  print(i)
  nomArchSalida <- paste(carpetaSalida, basename(pathsRasters[i]), sep='')
  if (reRun || !file.exists(nomArchSalida) || file.info(nomArchSalida)$size <= 0) {
    spI <- readGDAL(pathsRasters[i])
    # escala <- crearEscalaEquiespaciada(spI@data[, zcol], continuo = T)
    # shpBase <- spTransform(shpMask$shp, spI@proj4string)
    # mapearGrillaGGPlot(spI, shpBase = shpBase, escala = escala, dibujar = F)
    
    iNoNA <- !is.na(spI@data[,zcol])
    if (any(iNoNA)) {
      climatologiaMin <- readGDAL(pathsClimatologiaMin[clasesClimatologia[i]])
      climatologiaMax <- readGDAL(pathsClimatologiaMax[clasesClimatologia[i]])
      
      iMenoresAMin <- iNoNA & ((spI@data[, zcol] < minValAbs) | (!is.na(climatologiaMin@data[,zcol]) & (spI@data[, zcol] < climatologiaMin@data[, zcol])))
      iMayoresAMax <- iNoNA & ((spI@data[, zcol] > maxValAbs) | (!is.na(climatologiaMax@data[,zcol]) & (spI@data[, zcol] > climatologiaMax@data[, zcol])))
      spI@data[iMenoresAMin | iMayoresAMax, zcol] <- NA
      
      #mapearGrillaGGPlot(climatologiaMin, shpBase = shpBase, escala = escala, dibujar = F)
      #mapearGrillaGGPlot(climatologiaMax, shpBase = shpBase, escala = escala, dibujar = F)
      
      # mapearGrillaGGPlot(spI, shpBase = shpBase, escala = escala, dibujar = F)
    }
    writeGDAL(dataset = spI, fname = nomArchSalida, options = gdalOptions)
  }
}

filtrarSPDataFramesConClimatologia <- function(pathsRasters, clasesClimatologia, pathsClimatologiaMin, pathsClimatologiaMax, minValAbs=NA, maxValAbs=NA, 
                                               carpetaSalida='filtrados', zcol=1, gdalOptions = c('COMPRESS=DEFLATE', 'PREDICTOR=2', 'ZLEVEL=9'), nCoresAUsar=0,
                                               reRun=F) {
  #zcol=1
  #gdalOptions = c('COMPRESS=DEFLATE', 'PREDICTOR=2', 'ZLEVEL=9')
  #nCoresAUsar=0
  #reRun=F
  
  #pathsRasters <- dir(paste(pathDatos, 'MODIS/LST_Night_Combinada', sep=''), pattern = '*.tif$', full.names = T)
  #clasesClimatologia = diasDelAnio
  #minValAbs <- -8
  #maxValAbs <- 31
  
  #pathsClimatologiaMin <- paste(pathDatos, sprintf('LST_Night_Combinada_Clim_q0.005/%03d.tif', 1:365), sep='')
  #pathsClimatologiaMax <- paste(pathDatos, sprintf('LST_Night_Combinada_Clim_q0.9999/%03d.tif', 1:365), sep='')
  #carpetaSalida <- paste(pathDatos, ultimaCarpeta(pathsRasters[1]), '_0.0050_0.9999/', sep='')
  
  #pathsClimatologiaMin <- paste(pathDatos, sprintf('LST_Night_Combinada_Clim_q0.010/%03d.tif', 1:365), sep='')
  #pathsClimatologiaMax <- paste(pathDatos, sprintf('LST_Night_Combinada_Clim_q0.9995/%03d.tif', 1:365), sep='')
  #carpetaSalida <- paste(pathDatos, ultimaCarpeta(pathsRasters[1]), '_0.0100_0.9995/', sep='')
  
  #pathsClimatologiaMin <- paste(pathDatos, sprintf('LST_Night_Combinada_Clim_q0.015/%03d.tif', 1:365), sep='')
  #pathsClimatologiaMax <- paste(pathDatos, sprintf('LST_Night_Combinada_Clim_q0.9990/%03d.tif', 1:365), sep='')
  #carpetaSalida <- paste(pathDatos, ultimaCarpeta(pathsRasters[1]), '_0.0150_0.9990/', sep='')
  
  if (nCoresAUsar <= 0) nCoresAUsar <- min(length(pathsRasters), detectCores(T, T))
  
  dir.create(carpetaSalida, showWarnings = F, recursive = T)
  if (nCoresAUsar > 1) {
    cl <- makeCluster(getOption("cl.cores", nCoresAUsar))
    if (exists(x = 'setMKLthreads')) { clusterEvalQ(cl = cl, expr = setMKLthreads(1)) }
    clusterEvalQ(cl = cl, expr = {    
      require(rgdal)
    })
    parSapplyLB(cl=cl, X=1:length(pathsRasters), FUN=filtrarSPDataFrameIConClimatologia, pathsRasters=pathsRasters, clasesClimatologia=clasesClimatologia, 
                pathsClimatologiaMin=pathsClimatologiaMin, pathsClimatologiaMax=pathsClimatologiaMax, minValAbs=minValAbs, maxValAbs=maxValAbs, 
                carpetaSalida=carpetaSalida, zcol=zcol, gdalOptions = gdalOptions, reRun=reRun)
    stopCluster(cl)
  } else {
    sapply(X=1:length(pathsRasters), FUN=filtrarSPDataFrameIConClimatologia, pathsRasters=pathsRasters, clasesClimatologia=clasesClimatologia, 
           pathsClimatologiaMin=pathsClimatologiaMin, pathsClimatologiaMax=pathsClimatologiaMax, minValAbs=minValAbs, maxValAbs=maxValAbs, 
           carpetaSalida=carpetaSalida, zcol=zcol, gdalOptions = gdalOptions, reRun=reRun)
  }
}


filtrarRasterIGradienteAbrupto <- function(i, pathsRasters, pathsRastersCentrado=NULL, pathsRastersEscalado=NULL, minValAbs=NA, maxValAbs=NA, 
                                           carpetaSalida='filtrados', zcol=1, gdalOptions = c('COMPRESS=DEFLATE', 'PREDICTOR=2', 'ZLEVEL=9'), 
                                           nCoresAUsar=0, reRun=F, iPixelesSiempreNulos, iPixelesEnBolaChica, iPixelesEnBolaGrande, doPlots=F, 
                                           shpBase=NULL) {
  
  #doPlots=T
  # i <- 1
  #i <- which(pathsRasters == 'D:/Tesis/Datos/MODIS/LST_Night_Combinada/Combinada_2002-07-08.LST_Night_1km.tif')
  #i <- which(endsWith(pathsRasters, '2002-07-11.LST_Night_1km.tif'))
  #i <- which(endsWith(pathsRasters, '2002-07-13.LST_Night_1km.tif'))
  #i <- which(pathsRasters == 'D:/Tesis/Datos/MODIS/LST_Night_Combinada/Combinada_2002-07-19.LST_Night_1km.tif')
  #i <- which(pathsRasters == 'D:/Tesis/Datos/MODIS/LST_Night_Combinada/Combinada_2002-07-11.LST_Night_1km.tif')
  #i <- which(pathsRasters == 'D:/Tesis/Datos/MODIS/LST_Night_Combinada/Combinada_2002-10-13.LST_Night_1km.tif')
  #i <- which(pathsRasters == 'D:/Tesis/Datos/MODIS/LST_Night_Combinada/Combinada_2002-10-21.LST_Night_1km.tif')
  # i <- which(pathsRasters == 'D:/Tesis/Datos/MODIS/LST_Night_Combinada/Combinada_2004-01-30.LST_Night_1km.tif')
  #i <- which(pathsRasters == 'D:/Tesis/Datos/MODIS/LST_Night_Combinada/Combinada_2008-10-05.LST_Night_1km.tif')
  #i <- which(pathsRasters == 'D:/Tesis/Datos/MODIS/MOD11A1_LST_Night/MOD11A1_2002-07-19.LST_Night_1km.tif')
  print(i)
  if (!is.na(pathsRasters[i])) {
    nomArchSalida <- paste(carpetaSalida, basename(pathsRasters[i]), sep='')
    source(paste0(script.dir.funcionesAuxiliares, '../TryUtils/tryUtils.r'))
    
    run <- reRun || !file.exists(nomArchSalida) || file.info(nomArchSalida)$size <= 0 || !evaluarConReintentos(readGDAL(fname = nomArchSalida, silent = T), maxNIntentos = 1)
  
    if (run) {
      rasterI <- readGDAL(pathsRasters[i], silent = T)
      rasterI@data[!is.na(rasterI@data[, zcol]) & (rasterI@data[, zcol] < minValAbs | rasterI@data[, zcol] > maxValAbs), zcol] <- NA
      #plot(ecdf(x = rasterI@data[, zcol]))
      iNoNA <- !is.na(rasterI@data[,zcol])
      if (any(iNoNA) & length(unique(rasterI@data[iNoNA, zcol])) > 100) {
        if (doPlots) {
          escala <- crearEscalaEquiespaciada(rasterI@data[, zcol], continuo = T)
          if (!is.null(shpBase)) shpBase <- spTransform(shpBase, rasterI@proj4string)
          plots <- list()
          alturaEscalaContinua <- unit(0.075, 'npc')
          rango <- rangoExtendidoANDigitos(rasterI@data[, zcol], na.rm = T)
          
          plots[[length(plots) + 1]] <- mapearGrillaGGPlot(rasterI, shpBase = shpBase, escala = escala, dibujar = F, titulo = paste('Dato Original. [', rango[1], ', ', rango[2], ']', sep=''), alturaEscalaContinua=alturaEscalaContinua)
        }
  
        # Los píxeles que tienen al menos 10% de valores nulos en la bola chica son sospechosos
        # Los píxeles de arriba que tienen al menos 66.6% de nulos en la bola chica se descartan
        # Los píxeles de arriba que tienen una anomalía respecto al promedio en la bola grande sin la bola interior de más de 1 sd se descartan
        # Los píxeles en donde no tengo al menos 66.6% de los píxeles de la la bola grande sin la bola interior para calcular el promedio se descartan
        pctNulosEnBolaChica <- sapply(iPixelesEnBolaChica, 
                                      FUN = function(x, rasterI, zcol) {
                                        if (length(x) > 0) { return(sum(is.na(rasterI@data[x, zcol])) / length(x))
                                        } else { return(0) }
                                      } , rasterI, zcol)
        rasterI@data[pctNulosEnBolaChica > 2/3, zcol] <- NA
        
        anoms <- rasterI
        varAnoms <- var(anoms@data[,zcol], na.rm=T)
        if (!is.na(varAnoms) && varAnoms > 1E-3) {
          if (!is.null(pathsRastersCentrado) ) {
            rasterCentrado <- readGDAL(pathsRastersCentrado[i], silent = T)
            anoms@data[,zcol] <- anoms@data[, zcol] - rasterCentrado@data[, zcol]
            if (doPlots & F) plots[[length(plots) + 1]] <- mapearGrillaGGPlot(rasterCentrado, shpBase = shpBase, dibujar = F, escala = escala, titulo = 'Media Climatológica', alturaEscalaContinua=alturaEscalaContinua)
            rm(rasterCentrado)
            
            rasterICentrado <- anoms
          } else { rasterICentrado <- rasterI }
          limitesControlAnomalias <- quantile(rasterICentrado@data[, zcol], probs=c(1/3, 2/3), na.rm=T)
          
          if (!is.null(pathsRastersEscalado)) {
            rasterEscalado <- readGDAL(pathsRastersEscalado[i], silent = T)
            anoms@data[,zcol] <- anoms@data[, zcol] / rasterEscalado@data[, zcol]
            if (doPlots & F) {
              escalaSD <- crearEscalaEquiespaciada(rasterEscalado@data[,zcol], brewerPal = 'Reds', continuo = T)
              plots[[length(plots) + 1]] <- mapearGrillaGGPlot(rasterEscalado, shpBase = shpBase, dibujar = F, escala = escalaSD, titulo = 'Desv. Estándar Climatológica', alturaEscalaContinua=alturaEscalaContinua)
            }
            rm(rasterEscalado)
          }
          
          if (doPlots) { 
            if (!is.null(pathsRastersCentrado) && !is.null(pathsRastersEscalado)) { titulo <- 'Eliminación de Píxeles Aislados, Centrado y Escalado'
            } else if (!(is.null(pathsRastersCentrado))) { titulo <- 'Eliminación de Píxeles Aislados y Centrado'
            } else if (!(is.null(pathsRastersEscalado))) { titulo <- 'Eliminación de Píxeles Aislados Y Escalado'
            } else { titulo <- 'Eliminación de Píxeles Aislados' }
            escalaCentradoYEscalado <- crearEscalaEquiespaciada(anoms@data[, zcol], continuo = T)
            plots[[length(plots) + 1]] <- mapearGrillaGGPlot(anoms, shpBase = shpBase, dibujar = F, escala=escalaCentradoYEscalado, titulo = titulo, alturaEscalaContinua=alturaEscalaContinua)
          }
          
          pctNulosEnBolaGrande <- sapply(iPixelesEnBolaGrande, 
                                         FUN = function(x, anoms, zcol) {
                                           if (length(x) > 0) { return(sum(is.na(anoms@data[x, zcol])) / length(x)) 
                                           } else { return(0) }
                                         }, anoms, zcol)
          
          promedioEnBolaGrande <- sapply(seq_along(iPixelesEnBolaGrande), 
                                         FUN = function(x, anoms, zcol) {
                                           #return(mean(x = anoms@data[iPixelesEnBolaGrande[[x]], zcol], na.rm=T, trim=0.05))
                                           return(median(x = anoms@data[iPixelesEnBolaGrande[[x]], zcol], na.rm=T)) 
                                         }, anoms, zcol)
          promedioEnBolaGrande[pctNulosEnBolaGrande > 3/4 | is.nan(promedioEnBolaGrande)] <- NA
          
          if (doPlots) {
            mediaLocal <- rasterI
            mediaLocal@data[, zcol] <- promedioEnBolaGrande
            plots[[length(plots) + 1]] <- mapearGrillaGGPlot(mediaLocal, shpBase = shpBase, dibujar = F, escala = escalaCentradoYEscalado, titulo = 'Media Local', alturaEscalaContinua=alturaEscalaContinua)
            rm(mediaLocal)
          }
          
          anomsLocales <- anoms
          anomsLocales@data[, zcol] <- anoms@data[, zcol] - promedioEnBolaGrande
          # plot(density(anoms@data[,zcol], na.rm=T))
          # mediaAnoms <- mean(anoms@data[, zcol], na.rm=T)
          # sdAnoms <- sd(anoms@data[, zcol], na.rm=T)
          if (doPlots) {
            # escalaAnom <- crearEscalaEquiespaciada(anoms@data[, zcol], continuo = T)
            varAnoms <- var(anomsLocales@data[, zcol], na.rm = T)
            
            if (!is.na(varAnoms) & varAnoms > 1E-3) { escalaAnom <- crearEscalaEnQuantiles(datos = anomsLocales@data[, zcol], probs = c(0, 0.05, 0.1, 0.15, 0.2, 0.995, 0.999, 0.9995, 0.9999,1), continuo = F, brewerPal = 'RdYlBu', nDigitos = 1)
            } else { escalaAnom <- crearEscalaEquiespaciada(c(0, 1), continuo = F, brewerPal = 'RdYlBu') }
            
            plots[[length(plots) + 1]] <- mapearGrillaGGPlot(anomsLocales, shpBase = shpBase, escala = escalaAnom, dibujar = F, titulo = 'Anomalías', alturaEscalaContinua=alturaEscalaContinua)
          }
          
          lims <- quantile(anomsLocales@data[, zcol], probs=c(0.1, 0.9995), na.rm=T)
          lowerLim <- lims[1]
          upperLim <- lims[2]
          # medianaAnoms <- median(anomsLocales@data[, zcol], na.rm = T)
          # madAnoms <- mad(anomsLocales@data[, zcol], center = medianaAnoms, na.rm = T)
          # lowerLim <- medianaAnoms - 1.5 * madAnoms
          # upperLim <- medianaAnoms + 4 * madAnoms
          
          # Los píxeles que tengan menos de 5% de los píxeles nulos en la bola grande  se excluyen del control
          # anoms@data[pctNulosEnBolaChica <= 0.05 | (!is.na(rasterI@data[, zcol]) & rasterI@data[, zcol] > limitesControlAnomalias[1] & rasterI@data[, zcol] < limitesControlAnomalias[2]), zcol] <- NA
          # Los píxeles tengan anomalías negativas y su valor original sea mayor que limitesControlAnomalias[1] o 
          # los píxeles tengan anomalías positivas y su valor original sea menor que limitesControlAnomalias[2] se excluyen del control 
          anomsLocales@data[!is.na(anomsLocales@data[, zcol]) & 
                            (pctNulosEnBolaGrande <= 0.05 |
                            ((anomsLocales@data[, zcol] < 0 & rasterICentrado@data[, zcol] > limitesControlAnomalias[1]) | 
                             (anomsLocales@data[, zcol] > 0 & rasterICentrado@data[, zcol] < limitesControlAnomalias[2]))), zcol] <- NA
          iFiltroAnomaliasLocales <- !is.na(anomsLocales@data[,zcol]) & (anomsLocales@data[, zcol] < lowerLim | anomsLocales@data[,zcol] > upperLim) | is.na(promedioEnBolaGrande)
          if (doPlots) plots[[length(plots) + 1]] <- mapearGrillaGGPlot(anomsLocales, shpBase = shpBase, escala = escalaAnom, dibujar = F, titulo = 'Anomalías. Píxeles Sospechosos')
          anomsLocales@data[iFiltroAnomaliasLocales, zcol] <- NA
          
          rm(pctNulosEnBolaChica, anomsLocales)
          
          anoms[iFiltroAnomaliasLocales, zcol] <- NA
          mediana <- median(x = anoms@data[,zcol], na.rm = T)
          invDesvMedAbs <- 1/mad(anoms@data[,zcol], center = mediana, na.rm = T)
          zs <- (anoms@data[,zcol] - mediana)*invDesvMedAbs
          iFiltroAnomaliaGlobal <- !is.na(anoms@data[,zcol]) & ((pctNulosEnBolaGrande >= 0.05 & (zs < -2.5)) | (zs < -4 | zs > 5))
          # aux <- anoms
          # aux@data[,zcol] <- (anoms@data[,zcol] - mediana)*invDesvMedAbs
          # mapearGrillaGGPlot(aux, shpBase = shpBase, dibujar = F, continuo = F)
          rm(anoms, mediana, invDesvMedAbs, zs, pctNulosEnBolaGrande)
          
          filtrado <- rasterI
          filtrado@data[iFiltroAnomaliasLocales, zcol] <- NA
          filtrado@data[iFiltroAnomaliaGlobal, zcol] <- NA
          rm(iFiltroAnomaliaGlobal, iFiltroAnomaliasLocales)
          
          if (doPlots) {
            rango <- rangoExtendidoANDigitos(filtrado@data[, zcol], na.rm = T)
            
            plots[[length(plots) + 1]] <- mapearGrillaGGPlot(filtrado, shpBase = shpBase, dibujar = F, escala = escala, 
                                                             titulo = paste('Dato Filtrado. [', rango[1], ', ', rango[2], ']', sep=''),
                                                             alturaEscalaContinua=alturaEscalaContinua)
            if (length(plots) == 8) { 
              plots <- plots[c(1, 5, 2, 6, 3, 7, 4, 8)]
              cols <- 4
            } else { 
              plots <- plots[c(1, 4, 2, 5, 3, 6)] 
              cols <- 3
            }
            
            # print(multiplot(plotlist = plots, cols = cols))
            
            png(agregarCarpetaAlFinal(changeFileExt(nomArchSalida, '.png'), 'plots'), width = 1920, height = 1048)
            tryCatch(expr = print(multiplot(plotlist = plots, cols = cols)), finally = dev.off())
            rm(plots)
          }
        } else {
          filtrado <- anoms
        }
      } else { filtrado <- rasterI }
      
      writeGDAL(dataset = filtrado, fname = nomArchSalida, options = gdalOptions)
    }
  }
  return(NULL)
}

filtrarRastersGradienteAbrupto <- function(pathsRasters, pathsRastersCentrado=NULL, pathsRastersEscalado=NULL, minValAbs=NA, maxValAbs=NA, 
                                           carpetaSalida='filtrados', zcol=1, gdalOptions = c('COMPRESS=DEFLATE', 'PREDICTOR=2', 'ZLEVEL=9'), 
                                           nCoresAUsar=0, reRun=F, shpBase=NULL, r1=7000, r2=20000) {
  #r1=7000
  #r2=20000
  if (F) {
    zcol=1
    gdalOptions = c('COMPRESS=DEFLATE', 'PREDICTOR=2', 'ZLEVEL=9')
    nCoresAUsar=0
    reRun=F
    shpBase=shpMask$shp
    doPlots <-T
    
    #pathsRasters <- pathsRegresores[, 'LST_Night_Combinada']
    pathsRasters <- pathsRegresores[, 'MOD11A1_LST_Night']
    pathsRasters <- pathsRegresores[, 'MYD11A1_LST_Night']
    pathsRastersCentrado <- paste(pathDatos, sprintf('LST_Night_Combinada_Clim_mean/%03d.tif', diasDelAnio), sep='')
    pathsRastersEscalado <- paste(pathDatos, sprintf('LST_Night_Combinada_Clim_sd/%03d.tif', diasDelAnio), sep='')
    #pathsRastersEscalado <- NULL
    minValAbs <- -7
    maxValAbs <- 31
    carpetaSalida <- paste('Resultados/', ultimaCarpeta(pathsRasters[1]), '_filtrado/', sep='')
  }
  if (nCoresAUsar <= 0) nCoresAUsar <- min(length(pathsRasters), detectCores(T, T))
  
  # raster_Uno <- readGDAL(pathsRasters[which(pathsRasters == "D:/Tesis/Datos/MODIS/LST_Night_Combinada/Combinada_2002-10-22.LST_Night_1km.tif")])
  # raster_Dos <- readGDAL(pathsRasters[which(pathsRasters == "D:/Tesis/Datos/MODIS/LST_Night_Combinada/Combinada_2002-10-23.LST_Night_1km.tif")])
  # raster_Tres <- readGDAL(pathsRasters[which(pathsRasters == "D:/Tesis/Datos/MODIS/LST_Night_Combinada/Combinada_2002-11-07.LST_Night_1km.tif")])
  # raster_Cuatro <- readGDAL(pathsRasters[which(pathsRasters == "D:/Tesis/Datos/MODIS/LST_Night_Combinada/Combinada_2002-11-12.LST_Night_1km.tif")])
  # which(pathsRasters == "D:/Tesis/Datos/MODIS/LST_Night_Combinada/Combinada_2003-09-03.LST_Night_1km.tif")
  # which(pathsRasters == 'D:/Tesis/Datos/MODIS/MOD11A1_LST_Night/MOD11A1_2002-10-22.LST_Night_1km.tif')
  raster_Uno <- readGDAL(pathsRasters[107])
  raster_Dos <- readGDAL(pathsRasters[108])
  raster_Tres <- readGDAL(pathsRasters[123])
  raster_Cuatro <- readGDAL(pathsRasters[128])
  raster_Cinco <- readGDAL(pathsRasters[288])
  raster_Seis <- readGDAL(pathsRasters[389])
  raster_Siete <- readGDAL(pathsRasters[423])
  iPixelesSiempreNulos <- is.na(raster_Uno@data[,zcol]) & is.na(raster_Dos@data[,zcol]) & is.na(raster_Tres@data[,zcol]) & is.na(raster_Cuatro@data[,zcol]) & 
                          is.na(raster_Cinco@data[,zcol]) & is.na(raster_Seis@data[,zcol]) & is.na(raster_Siete@data[,zcol])
  rm(raster_Uno, raster_Dos, raster_Tres, raster_Cuatro, raster_Cinco, raster_Seis, raster_Siete)
  
  # Calculo de bolas
  rasterI <- readGDAL(pathsRasters[1], silent=T)
  puntosRaster <- as(geometry(rasterI), 'SpatialPoints')
  
  pathCache <- paste(dirname(agregarCarpetaAlFinal(pathsRasters[1], carpeta = 'RCache')), '/', digest(object = list(puntosRaster, r1, r2, version=1)), '.RData', sep='')
  if (!file.exists(pathCache) | (file.info(pathCache)$size <= 0)) {
    parOver <- function(nCoresAUsar, spgeom1, spgeom2, returnList=FALSE, fn=NULL, ...) {
      auxFunc <- function(x, spgeom1, spgeom2, returnList = FALSE, fn = NULL, ...) { 
        return(over(spgeom1[x,], spgeom2, returnList, fn, ...)) 
      }
      if (nCoresAUsar <= 0) nCoresAUsar <- min(length(spgeom1), detectCores(T, T))
      
      if (nCoresAUsar > 1) {
        cl <- makeCluster(getOption("cl.cores", nCoresAUsar))
        clusterEvalQ(cl, {require('sp')})
        iSplit <- clusterSplit(cl, 1:length(spgeom1))
        
        resultado <- parSapply(cl = cl, X = iSplit, FUN = auxFunc, spgeom1, spgeom2, returnList, fn, ...)
        #resultado <- parSapply(cl = cl, X = iSplit, FUN = auxFunc, spgeom1, spgeom2, returnList, fn)
        stopCluster(cl)
        
        # x <- iSplit[[1]]
        # over(spgeom1[x,], spgeom2, returnList, fn, ...)
        
        return(do.call(c, resultado))
      } else { return(over(spgeom1, spgeom2, returnList, fn, ...)) }
    }
    
    #bolaInterior <- gBuffer(puntosRaster, byid = T, width = 1000 * 1)
    #iPixelesEnBolaInterior <- parOver(nCoresAUsar = nCoresAUsar, spgeom1 = bolaInterior, spgeom2 = puntosRaster, returnList = T)
    #rm(bolaInterior)
    
    nCoresAUsar <- 1
    
    bolaChica <- gBuffer(puntosRaster, byid = T, width = r1)
    system.time(iPixelesEnBolaChica <- over(bolaChica, puntosRaster, returnList = T))
    #iPixelesEnBolaChica <- parOver(nCoresAUsar = nCoresAUsar, spgeom1 = bolaChica, spgeom2 = puntosRaster, returnList = T)
    rm(bolaChica)
    
    # Saco los siempre nulos de la bola chica y el píxel en cuestión para que no se consideren al buscar sospechosos
    cl <- makeCluster(getOption("cl.cores", nCoresAUsar))
    iPixelesEnBolaChica <- parLapplyLB(cl, seq_along(iPixelesEnBolaChica),
                                       fun = function(x, iPixelesEnBolaChica, iPixelesSiempreNulos) {
                                         iX <- setdiff(iPixelesEnBolaChica[[x]], x)
                                         return(iX[!iPixelesSiempreNulos[iX]])
                                       }, iPixelesEnBolaChica=iPixelesEnBolaChica, iPixelesSiempreNulos=iPixelesSiempreNulos)
    stopCluster(cl)
    gc(T)
    print('Bola chica pronta')
    
    bolaGrande <- gBuffer(puntosRaster, byid = T, width = r2)
    system.time(iPixelesEnBolaGrande <- over(bolaGrande, puntosRaster, returnList = T))
    #iPixelesEnBolaGrande <- parOver(nCoresAUsar = nCoresAUsar, spgeom1 = bolaGrande, spgeom2 = puntosRaster, returnList = T)
    rm(bolaGrande)
    gc(T)
    
    print('Bola grande creada')
    # Saco los siempre nulos de la bola grande para que trate de usarlos sin que estén cada vez que hace los promedios (performance)
    cl <- makeCluster(getOption("cl.cores", nCoresAUsar))
    iPixelesEnBolaGrande <- parLapplyLB(cl, seq_along(iPixelesEnBolaGrande),
                                        fun = function(x, iPixelesEnBolaGrande, iPixelesSiempreNulos) {
                                          iX <- setdiff(iPixelesEnBolaGrande[[x]], x)
                                          return(iX[!iPixelesSiempreNulos[iX]])
                                        }, iPixelesEnBolaGrande=iPixelesEnBolaGrande, iPixelesSiempreNulos=iPixelesSiempreNulos)
    stopCluster(cl)
    print('Bola grande pronta')
    
    # Obtengo el conjunto de píxeles de la bola grande sin la bola interior para calcular los promedios
    #iPixelesEnBolaGrande <- parLapplyLB(cl, 1:length(iPixelesEnBolaGrande), 
    #                                    fun = function(x, iPixelesEnBolaGrande, iPixelesEnBolaInterior) { 
    #                                                   setdiff(x = iPixelesEnBolaGrande[[x]], y = c(iPixelesEnBolaInterior[[x]], x))
    #                                    }, iPixelesEnBolaGrande=iPixelesEnBolaGrande, iPixelesEnBolaInterior=iPixelesEnBolaInterior)
  
    #invDistanciasPixelesEnBolaGrande <- parLapplyLB(cl, 1:length(iPixelesEnBolaGrandeSinBolaInterior), 
    #                                                function(x, puntosRaster, iPixelesEnBolaGrandeSinBolaInterior) { 1/gDistance(puntosRaster[x,], puntosRaster[iPixelesEnBolaGrandeSinBolaInterior[[x]], ], byid = T) }, 
    #                                                puntosRaster=puntosRaster, iPixelesEnBolaGrandeSinBolaInterior=iPixelesEnBolaGrandeSinBolaInterior)

    bolas <- list(iPixelesEnBolaChica, iPixelesEnBolaGrande)
    guardarCache(pathCache = pathCache, obj = bolas)
  } else {
    bolas <- cargarCache(pathCache)
    iPixelesEnBolaChica <- bolas[[1]]
    iPixelesEnBolaGrande <- bolas[[2]]
  }
  rm(bolas)
  gc(T)
  
  dir.create(carpetaSalida, showWarnings = F, recursive = T)
  dir.create(paste(carpetaSalida, 'plots/', sep=''), showWarnings = F, recursive = T)
  
  if (nCoresAUsar > 1) {
    cl <- makeCluster(getOption("cl.cores", nCoresAUsar))
    if (exists(x = 'setMKLthreads')) { clusterEvalQ(cl = cl, expr = setMKLthreads(1)) }
    clusterExport(cl, varlist = c('script.dir.funcionesAuxiliares'))
    clusterEvalQ(cl = cl, expr = {
      require(rgdal)
      require(sp)
      require(Rmisc)
      source(paste0(script.dir.funcionesAuxiliares, 'mapearEx.r'))
    })
    parSapplyLB(cl=cl, X=1:length(pathsRasters), FUN=filtrarRasterIGradienteAbrupto, pathsRasters=pathsRasters, pathsRastersCentrado=pathsRastersCentrado, 
                pathsRastersEscalado=pathsRastersEscalado, minValAbs=minValAbs, maxValAbs=maxValAbs, carpetaSalida=carpetaSalida, zcol=zcol, 
                gdalOptions = gdalOptions, reRun=reRun, iPixelesSiempreNulos=iPixelesSiempreNulos, iPixelesEnBolaChica=iPixelesEnBolaChica, 
                iPixelesEnBolaGrande=iPixelesEnBolaGrande, doPlots=doPlots, shpBase=shpBase)
    stopCluster(cl)
  } else {
    sapply(X=1:length(pathsRasters), FUN=filtrarRasterIGradienteAbrupto, pathsRasters=pathsRasters, pathsRastersCentrado=pathsRastersCentrado, 
           pathsRastersEscalado=pathsRastersEscalado, minValAbs=minValAbs, maxValAbs=maxValAbs, carpetaSalida=carpetaSalida, zcol=zcol, 
           gdalOptions = gdalOptions, reRun=reRun, iPixelesSiempreNulos=iPixelesSiempreNulos, iPixelesEnBolaChica=iPixelesEnBolaChica, 
           iPixelesEnBolaGrande=iPixelesEnBolaGrande, doPlots=doPlots, shpBase=shpBase)
  }
}

detectarDiscontinuidadesEnRasters <- function(pathsRasters, pathsRastersCentrado=NULL, umbralNuggetSobrePsill= 0.66) {
  #pathsRasters <- pathsRegresores[, 'MOD11A1_LST_Night_FRv2']
  #pathsRasters <- pathsRegresores[, 'MYD11A1_LST_Night_FRv2']
  #pathsRastersCentrado <- paste(pathDatos, sprintf('LST_Night_Combinada_Clim_mean/%03d.tif', diasDelAnio), sep='')
  #umbralNuggetSobrePsill= 0.5
    
  #i <- which(fechasObservaciones == as.POSIXct('2014-01-01', tz=tz(fechasObservaciones[1])))
  #i <- which(fechasObservaciones == as.POSIXct('2014-01-02', tz=tz(fechasObservaciones[1])))
  #i <- which(fechasObservaciones == as.POSIXct('2014-12-25', tz=tz(fechasObservaciones[1])))
  
  elaborarBDRastersI <- function(i, pathsRasters, pathsRastersCentrado) {
    rasterI <- readGDAL(pathsRasters[i], silent = T)
    rasterCentradoI <- readGDAL(pathsRastersCentrado[i], silent = T)
    rasterCentradoI <- SpatialGridDataFrameEnMtoSpatialGridDataFrameEnKm(rasterCentradoI)
    
    # mapearGrillaGGPlot(rasterI, shpBase = shpBase, dibujar = F, continuo = T)
    
    rasterI <- raster(rasterI)
    rasterCentradoI <- raster(rasterCentradoI)
    rasterICentrado <- rasterI - rasterCentradoI
    
    var3x3 <- focal(rasterICentrado, w=matrix(1,3,3), fun=var)
    var9x9 <- focal(rasterICentrado, w=matrix(1,9,9), fun=var)
    var51x51 <- focal(rasterICentrado, w=matrix(1,51,51), fun=var)
    
    plot(var3x3)
    plot(var9x9)
    plot(var51x51)
    
    hist(getValues(var51x51))
    
    mean(var3x3[var3x3>1.5])
    sum(lele < 0.1 | lele > 1.5) / length(lele)
    
  }
  
  detectarDiscontinuidadesEnRastersI <- function(i, pathsRasters, pathsRastersCentrado, umbralNuggetSobrePsill) {
    rasterI <- readGDAL(pathsRasters[i], silent = T)
    rasterCentradoI <- readGDAL(pathsRastersCentrado[i], silent = T)
    
    lala <- coarsenGrid(rasterI, coarse = 6)
    length(lala)
    
    centroides <- SpatialPoints(rasterI, proj4string = rasterI@proj4string)
    if (!identicalCRS(centroides, rasterCentradoI)) {
      centroides <- spTransform(centroides, rasterCentradoI@proj4string)
    }
    
    iCentroides <- over(centroides, geometry(rasterCentradoI))
    
    centroides <- SpatialPointsDataFrame(centroides, data=data.frame(value=rasterI@data[,1] - rasterCentradoI@data[iCentroides,1]))
    # centroides <- SpatialPointsDataFrame(centroides, data=data.frame(value=rasterI@data[,1]))
    iNoNA <- !is.na(centroides@data$value)
    centroides <- centroides[iNoNA,]
    
    #mapearGrillaGGPlot(grilla = rasterI, dibujar=F)
    #mapearGrillaGGPlot(grilla = rasterCentradoI, dibujar=F)
    #rasterICentrado <- rasterI
    #rasterICentrado@data[,1] <- rasterI@data[,1] - rasterCentradoI@data[iCentroides,1]
    #mapearGrillaGGPlot(grilla = rasterICentrado, dibujar=F)
    #hist(rasterI@data[,1])
    #hist(rasterICentrado@data[,1])
    
    iMuestras <- muestrearEnCuadrantesYECDF(centroides, size = 3000)
  
    # mapearPuntosGGPlot(centroides[iMuestras,], dibujar = F)
    
    system.time(vc <- variogram(object=value~1, data = centroides[iMuestras,], cloud=T))
    
    qs <- quantile(vc$dist, probs=c(0.15, 0.9))
    mean(vc$gamma[vc$dist <= qs[1]])
    mean(vc$gamma[vc$dist >= qs[2]])
    
    mean(vc$gamma[vc$dist <= qs[1]]) / mean(vc$gamma[vc$dist >= qs[2]])
  }
  
  detectarDiscontinuidadesEnRastersI_v2 <- function(i, pathsRasters, pathsRastersCentrado, umbralNuggetSobrePsill) {
    rasterI <- readGDAL(pathsRasters[i], silent = T)
    rasterCentradoI <- readGDAL(pathsRastersCentrado[i], silent = T)
    rasterCentradoI <- SpatialGridDataFrameEnMtoSpatialGridDataFrameEnKm(rasterCentradoI)
    
    mapearGrillaGGPlot(rasterI, shpBase = shpBase, dibujar = F, continuo = T)
    
    rasterI <- raster(rasterI)
    rasterCentradoI <- raster(rasterCentradoI)
    rasterICentrado <- rasterI - rasterCentradoI
    plot(rasterICentrado)
    
    lala <- focal(rasterICentrado, w=matrix(1,3,3), fun=sd)
    lele <- na.omit(getValues(lala))
    
    sort(unique(lele))
    
    plot(lala)
    hist(lele)
    
      mean(lele[lele>1.5])
    sum(lele < 0.1 | lele > 1.5) / length(lele)
  }
  
  
  
  
  

  hist(focal(rasterICentrado, w=matrix(1,3,3), fun=sd), freq=FALSE)
  plot(focal(rasterICentrado, w=matrix(c(0,1,0,1,-4,1,0,1,0), nrow=3)))
  plot(focal(rasterICentrado, w=matrix(c(1,2,1,0,0,0,-1,-2,-1) / 4, nrow=3)))
  
  lala <- focal(rasterICentrado, w=matrix(c(2,4,5,4,2,4,9,12,9,4,5,12,15,12,5,4,9,12,9,4,2,4,5,4,2), nrow = 5))
  plot(focal(lala, w=matrix(c(0,1,0,1,-4,1,0,1,0), nrow=3)))
  
  
  
  
  hist(focal(rasterICentrado, w=matrix(c(1,2,1,0,0,0,-1,-2,-1) / 4, nrow=3)))
  
  plot( boundaries(rasterICentrado, type='inner') )
  plot( boundaries(rasterICentrado, type='outer') )
  plot( boundaries(rasterICentrado, classes=TRUE) )
  
  
}

outlyingnessMediaSD <- function(x, na.rm=T, nMinParaEstimar = 10, sdMin=NA) {
  if (na.rm) { xSinNA <- x[!is.na(x)]
  } else { xSinNA <- x }
  
  if (length(xSinNA) >= nMinParaEstimar) {
    media <- mean(xSinNA, na.rm = FALSE)
    stDev <- sd(xSinNA, na.rm = FALSE)
    
    if (is.na(sdMin) || stDev > sdMin) {
      return((x - media) / stDev)
    } else {
      res <- rep(NA_real_, length(x))
      names(res) <- names(x)
      return(res) 
    }
  } else { 
    res <- rep(NA_real_, length(x))
    names(res) <- names(x)
    return(res) 
  }
}

filtroMediaSD <- function(x, factorSDHaciaAbajo=3, factorSDHaciaArriba=factorSDHaciaAbajo, na.rm=T, nMinParaAplicarFiltro = 10) {
  if (na.rm) { xSinNA <- x[!is.na(x)]
  } else { xSinNA <- x }
  
  if (length(xSinNA) >= nMinParaAplicarFiltro) {
    media <- mean(x, na.rm = na.rm)
    desvEst <- sd(x, na.rm = na.rm)
    
    # (x - media) / desvEst
    
    res <- x > media - factorSDHaciaAbajo * desvEst & x < media + factorSDHaciaArriba * desvEst
    res[is.na(x)] <- TRUE
    return(res)
  } else {     
    res <- rep(TRUE, length(x))
    names(res) <- names(x)
    return(res) 
  }
}

outlyingnessMedianaMAD <- function(x, na.rm=T, constanteMAD=1.4826, nMinParaEstimar = 10, 
                                   desvMedAbsMin=NA) {
  if (na.rm) { xSinNA <- x[!is.na(x)]
  } else { xSinNA <- x }
  
  if (length(xSinNA) >= nMinParaEstimar) {
    mediana <- median(xSinNA, na.rm = FALSE)
    desvMedAbs <- mad(xSinNA, center = mediana, na.rm = FALSE, constant = constanteMAD)
    
    if (is.na(desvMedAbsMin) || desvMedAbs > desvMedAbsMin) {
      return((x - mediana) / desvMedAbs)
    } else {
      res <- rep(NA_real_, length(x))
      names(res) <- names(x)
      return(res) 
    }
  } else { 
    res <- rep(NA_real_, length(x))
    names(res) <- names(x)
    return(res) 
  }
}

filtroMedianaMAD <- function(x, factorMADHaciaAbajo=3, factorMADHaciaArriba=factorMADHaciaAbajo, na.rm=T, constanteMAD=1.4826, nMinParaAplicarFiltro = 10) {
  if (na.rm) { xSinNA <- x[!is.na(x)]
  } else { xSinNA <- x }
  
  if (length(xSinNA) >= nMinParaAplicarFiltro) {
    mediana <- median(xSinNA, na.rm = FALSE)
    desvMedAbs <- mad(xSinNA, center = mediana, na.rm = FALSE, constant = constanteMAD)
    
    #factorMADHaciaAbajo <- 2.5
    #factorMADHaciaArriba <- factorMADHaciaAbajo
    
    # (x - mediana) / desvMedAbs
    
    res <- x > mediana - factorMADHaciaAbajo * desvMedAbs & x < mediana + factorMADHaciaArriba * desvMedAbs
    res[is.na(x)] <- TRUE
    return(res)
  } else { 
    res <- rep(TRUE, length(x))
    names(res) <- names(x)
    return(res) 
  }
}

seleccionarMejorRegresorPorFechasI <- function(i, fechasObservaciones, valoresObservaciones, coordsObservaciones, regresoresCandidatos, minNoNAs=10) {
  print(i)
  # i <- 4415
  regresoresCandidatosI <- regresoresCandidatos[i,]
  coordsObservaciones$value <- valoresObservaciones[i, ]
  mses <- numeric(length = length(regresoresCandidatosI))
  
  j <- 1
  for (j in seq_along(regresoresCandidatosI)) {
    if (!is.na(regresoresCandidatosI[j])) {
      reg <- readGDAL(regresoresCandidatosI[j], silent = T)
      if (!identicalCRS(reg, coordsObservaciones)) { 
        auxSP <- spTransform(coordsObservaciones, reg@proj4string)
      } else { 
        auxSP <- coordsObservaciones 
      }
      
      regVals <- over(geometry(auxSP), reg)[1]
      xMat <- cbind(regVals, sp::coordinates(auxSP))
      
      iNoNAsReg <- complete.cases(xMat)
      df <- na.omit(data.frame(y = valoresObservaciones[i, ], xMat))
      if (sum(!is.na(valoresObservaciones[i,])) <= 3) {
        # Si no hay al menos 4 observaciones no puedo hacer regresión. Hago el mse sin ajustar
        mses[j] <- mean((df$y - df$band1)^2)
      } else if (sum(iNoNAsReg) >= minNoNAs) {
        modelo <- rlm(formula = y ~ ., data = df)
        mses[j] <- mean(modelo$residuals^2)
      } else { 
        mses[j] <- NA
      }
    } else { mses[j] <- NA }
  }
  if (any(!is.na(mses))) { return(regresoresCandidatosI[which.min(mses)])
  } else { return(regresoresCandidatosI[1]) }
}

seleccionarMejorRegresorPorFechas <- function(fechasObservaciones, valoresObservaciones, coordsObservaciones, regresoresCandidatos, minNoNAs=10, nCoresAUsar=0, 
                                              nombreRegresorResultado = 'Combinado') {
  #regresoresCandidatos <- pathsRegresores[, c( 'LST_Night_Combinada_Clim_mean', 'MOD11A1_LST_Night_FR', 'MYD11A1_LST_Night_FR', 'LST_Night_Combinada_Clim_median')]
  pathCache <- getPathCache(c(fechasObservaciones, valoresObservaciones, coordsObservaciones, regresoresCandidatos, minNoNAs, version=3))
  nCoresAUsar=0
  if (!file.exists(pathCache)) {
    if (nCoresAUsar <= 0) { nCoresAUsar <- min(detectCores(T, T), nrow(regresoresCandidatos)) }
    
    if (nCoresAUsar > 1) {
      cl <- makeCluster(getOption('cl.cores', nCoresAUsar))
      clusterEvalQ(cl, expr = { 
        require(rgdal)
        require(sp)
        require(MASS)
      })
      res <- parSapplyLB(cl = cl, X=1:nrow(regresoresCandidatos), FUN = seleccionarMejorRegresorPorFechasI, 
                         fechasObservaciones=fechasObservaciones, coordsObservaciones=coordsObservaciones, valoresObservaciones=valoresObservaciones, 
                         regresoresCandidatos=regresoresCandidatos, minNoNAs=minNoNAs)
      stopCluster(cl)
    } else {
      res <- sapply(X=1:nrow(regresoresCandidatos), FUN = seleccionarMejorRegresorPorFechasI, 
                    fechasObservaciones=fechasObservaciones, coordsObservaciones=coordsObservaciones, valoresObservaciones=valoresObservaciones, 
                    regresoresCandidatos=regresoresCandidatos, minNoNAs=minNoNAs)
    }
    
    resultado <- matrix(data = res, ncol = 1, dimnames = list(rownames(regresoresCandidatos), nombreRegresorResultado))
    guardarCache(pathCache, obj = resultado)
  } else { resultado <- cargarCache(pathCache = pathCache) }
  
  return(resultado)  
}

calcularIndicatrizCobertura <- function(pathsGeoTiffs, nCuadrantesX=2, nCuadrantesY=nCuadrantesX, 
                                        dispMinPorCuadrante=0.25) {
  # pathsGeoTiffs <- pathsRegresores[, 'MOD11A1_LST_Night_filtrado', drop=F]
  nNoNulosPorCuadrantes <- contarNoNulosPorCuadrantes(pathsGeoTiffs = pathsGeoTiffs, nCuadrantesX = nCuadrantesX, 
                                                      nCuadrantesY = nCuadrantesY, carpetaSalida = NULL)
  
  evaluarConReintentos(rasterI <- readGDAL(fname = pathsGeoTiffs[1], silent = T))
  objCuadrantes <- dividirEnCuadrantes(object = rasterI, nCuadrantesX = nCuadrantesX, nCuadrantesY = nCuadrantesY)
  ies <- over(geometry(rasterI), objCuadrantes, returnList = F)
  
  nPixelesPorCuadrante <- sapply(seq.int(from=1, to = length(objCuadrantes), by = 1), FUN = function(x, ies) { return(sum(ies==x))}, ies=ies)
  rm(rasterI, objCuadrantes, ies)
  
  # res <- rowMaxs(x = as.matrix(nNoNulosPorCuadrantes[,1:4] / nPixelesPorCuadrante)) >= dispMinPorCuadrante
  return(rowMaxs(x = as.matrix(nNoNulosPorCuadrantes[,1:4] / nPixelesPorCuadrante)) >= dispMinPorCuadrante)
}

calcularIndicatrizRangoMaximo <- function(pathsGeoTiffs, pathsGeoTiffsCentrado=NULL, pathsGeoTiffsEscalado=NULL, 
                                          rangoMax) {
  pathsGeoTiffs <- pathsRegresores[, 'MOD11A1_LST_Night_filtrado', drop=F]
  nNoNulosPorCuadrantes <- contarNoNulosPorCuadrantes(pathsGeoTiffs = pathsGeoTiffs, nCuadrantesX = nCuadrantesX, 
                                                      nCuadrantesY = nCuadrantesY, carpetaSalida = NULL)
  
  evaluarConReintentos(rasterI <- readGDAL(fname = pathsGeoTiffs[1], silent = T))
  objCuadrantes <- dividirEnCuadrantes(object = rasterI, nCuadrantesX = nCuadrantesX, nCuadrantesY = nCuadrantesY)
  ies <- over(geometry(rasterI), objCuadrantes, returnList = F)
  
  nPixelesPorCuadrante <- sapply(seq.int(from=1, to = length(objCuadrantes), by = 1), FUN = function(x, ies) { return(sum(ies==x))}, ies=ies)
  rm(rasterI, objCuadrantes, ies)
  
  # res <- rowMaxs(x = as.matrix(nNoNulosPorCuadrantes[,1:4] / nPixelesPorCuadrante)) >= dispMinPorCuadrante
  return(rowMaxs(x = as.matrix(nNoNulosPorCuadrantes[,1:4] / nPixelesPorCuadrante)) >= dispMinPorCuadrante)
}


convertirUnidadCoordenadasSpatialGrid <- function(spGrid, objSPProjDestino) {
  spGrid <- rasterI
  objSPProjDestino <- shpMask$shp
  
  tokensP4strSPGrid <- getTokens(proj4string(spGrid), separadorTokens=' ')
  tokensP4strObjDest <- getTokens(proj4string(objSPProjDestino), separadorTokens=' ')
  
  u1 <- tokensP4strSPGrid[tokensP4strSPGrid[,1] == '+units', 2]
  u2 <- tokensP4strObjDest[tokensP4strObjDest[,1] == '+units', 2]
  
  if (u1 == 'm') {
    if (u2 == 'km') {
      
    }
  }
  
  sp::coordinates(spGrid) <- sp::coordinates(spGrid) / 1000
}

SpatialGridDataFrameEnMtoSpatialGridDataFrameEnKm <- function(grillaM) {
  grillaKm <- grillaM@grid
  grillaKm@cellcentre.offset <- grillaKm@cellcentre.offset * 0.001
  grillaKm@cellsize <- grillaKm@cellsize * 0.001
  # TODO: Agregar el wkt de la proj4string
  return(SpatialGridDataFrame(
    grid=grillaKm, 
    proj4string=gsub(pattern='+units=m', replacement='+units=km', x=proj4string(grillaM), fixed=T), 
    data=grillaM@data))
}

mapearRastersCentradoYEscaladoI <- function(i, pathsRasters, pathsRastersCentrado, pathsRastersEscalado, pathsRastersRellenado,
                                            shpBase, carpetaSalida, convertirRastersDeMaKm=T, replot=T) {
  #i <- 1
  print(i)
  
  if (!is.na(pathsRasters[i])) {
    nomArchSalida <- paste(carpetaSalida, rownames(pathsRasters)[i], '.png', sep='')
    if (replot | !file.exists(nomArchSalida)) {
      evaluarConReintentos(rasterI <- readGDAL(fname = pathsRasters[i], silent = T))
      evaluarConReintentos(rasterCentradoI <- readGDAL(fname = pathsRastersCentrado[i], silent = T))
      evaluarConReintentos(rasterEscaladoI <- readGDAL(fname = pathsRastersEscalado[i], silent = T))
      evaluarConReintentos(rasterRellenadoI <- readGDAL(fname = pathsRastersRellenado[i], silent = T))
      
      if (convertirRastersDeMaKm) rasterI <- SpatialGridDataFrameEnMtoSpatialGridDataFrameEnKm(grillaM = rasterI)
      
      escala <- crearEscalaEquiespaciada(c(rasterI@data[,1], rasterRellenadoI@data[,1]), nIntervalos = 10, continuo = T)

      rango <- rangoExtendidoANDigitos(x = rasterI@data[,1], nDigitos = 1, na.rm = T)
      cobertura <- round((sum(!is.na(rasterI@data[,1])) / length(rasterI)) * 100, 1)
      varianza <- round(var(rasterI@data[,1], na.rm = T), 1)
      g1 <- mapearGrillaGGPlot(grilla = rasterI, shpBase = shpBase, escala = escala,
                               titulo = paste(rownames(pathsRasters)[i], '-', colnames(pathsRasters)[1]), 
                               subtitulo = paste('Rango = [', rango[1], ', ', rango[2], ']. Diff = ', rango[2] - rango[1], '. Cobertura = ', cobertura, '%. Var = ', varianza, '.' , sep=''),
                               continuo = T, dibujar = F)
      
      rasterI@data[, 1] <- (rasterI@data[, 1] - rasterCentradoI@data[, 1]) / rasterEscaladoI@data[, 1]
      rango <- rangoExtendidoANDigitos(x = rasterI@data[,1], nDigitos = 1, na.rm = T)
      g2 <- mapearGrillaGGPlot(grilla = rasterI, shpBase = shpBase, 
                               titulo = paste(rownames(pathsRasters)[i], '-', colnames(pathsRasters)[1], '. Centrado y Escalado.'), 
                               subtitulo = paste('Rango = [', rango[1], ', ', rango[2], ']. Diff = ', rango[2] - rango[1], '. Cobertura = ', cobertura, '%.', sep=''),
                               continuo = T, dibujar = F)
      
      rango <- rangoExtendidoANDigitos(x = rasterRellenadoI@data[,1], nDigitos = 1, na.rm = T)
      iRellenados <- !is.na(rasterRellenadoI@data[,1]) & is.na(over(geometry(rasterRellenadoI), rasterI))
      varianzaRellenados <- round(var(rasterRellenadoI@data[iRellenados, 1], na.rm = T), 1)
      
      # varianzaRellenados <- round(var(rasterRellenadoI@data[, 1], na.rm = T), 1)
      g3 <- mapearGrillaGGPlot(grilla = rasterRellenadoI, shpBase = shpBase, escala = escala,
                               titulo = paste(rownames(pathsRasters)[i], '-', colnames(pathsRastersRellenado)[1]), 
                               subtitulo = paste('Rango = [', rango[1], ', ', rango[2], ']. Diff = ', rango[2] - rango[1], '. Var = ', varianzaRellenados, '. PctVar = ', round((varianzaRellenados / varianza) * 100, 1), '%.', sep=''),
                               continuo = T, dibujar = F)
      
      dir.create(carpetaSalida, showWarnings = F, recursive = T)
      oldSciPen <- getOption("scipen")
      options(scipen=15)
      png(nomArchSalida, width = 1920, height = 660)
      tryCatch(expr = print(multiplot(plotlist=list(g1,g2,g3), cols=3)), finally = dev.off())
      options(scipen = oldSciPen)
    }
  }
}

mapearRastersCentradoYEscalado <- function(pathsRasters, pathsRastersCentrado, pathsRastersEscalado, pathsRastersRellenado, 
                                           shpBase, carpetaSalida, convertirRastersDeMaKm=T, nCoresAUsar=0, replot=T) {
  # i <- 1
  # pathsRasters <- pathsRegresores[, 'MOD11A1_LST_Night_filtrado', drop=F]
  # pathsRastersCentrado <- paste(pathDatos, sprintf('LST_Night_Combinada_Clim_mean/%03d.tif', diasDelAnio), sep='')
  # pathsRastersEscalado <- paste(pathDatos, sprintf('LST_Night_Combinada_Clim_sd/%03d.tif', diasDelAnio), sep='')
  # pathsRastersRellenado <- pathsRegresores[, 'MOD11A1_LST_Night_FR', drop=F]
  # carpetaSalida <- paste(dirname(pathsRasters[1]), '/plotsCentradoYEscalado/', sep='')
  # shpBase <- shpMask$shp
  # convertirRastersDeMaKm <- T
  # replot <- F
  # nCoresAUsar <- 0
  
  # tail(cbind(pathsRasters, pathsRastersCentrado, pathsRastersEscalado, pathsRastersRellenado))
  
  if (nCoresAUsar <= 0) { nCoresAUsar <- min(detectCores(T, T), nrow(pathsRasters)) }
  
  if (nCoresAUsar > 1) {
    cl <- makeCluster(getOption('cl.cores', nCoresAUsar))
    clusterExport(cl, varlist = c('script.dir.funcionesAuxiliares'))
    clusterEvalQ(cl = cl, expr = {
      source(paste0(script.dir.funcionesAuxiliares, '../TryUtils/tryUtils.r'))
      source(paste0(script.dir.funcionesAuxiliares, 'mapearEx.r'))
      source(paste0(script.dir.funcionesAuxiliares, 'funcionesAuxiliares.r'))
    })
    
    res <- parSapplyLB(cl = cl, X=1:nrow(pathsRasters), FUN = mapearRastersCentradoYEscaladoI, pathsRasters=pathsRasters, 
                       pathsRastersCentrado=pathsRastersCentrado, pathsRastersEscalado=pathsRastersEscalado, pathsRastersRellenado=pathsRastersRellenado,
                       shpBase=shpBase, carpetaSalida=carpetaSalida, convertirRastersDeMaKm=convertirRastersDeMaKm, replot=replot)
    stopCluster(cl)
  } else {
    res <- sapply(X=1:nrow(pathsRasters), FUN = mapearRastersCentradoYEscaladoI, pathsRasters=pathsRasters, 
                  pathsRastersCentrado=pathsRastersCentrado, pathsRastersEscalado=pathsRastersEscalado,  pathsRastersRellenado=pathsRastersRellenado,
                  shpBase=shpBase, carpetaSalida=carpetaSalida, convertirRastersDeMaKm=convertirRastersDeMaKm, replot=replot)
  }
}
