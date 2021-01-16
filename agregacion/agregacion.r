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

iFrame <- sys.nframe()
if (iFrame >= 3) { script.dir.agregacion <- sys.frame(iFrame - 3)$ofile 
} else { script.dir.agregacion <- NULL }
while ((is.null(script.dir.agregacion) || is.na(regexpr('agregacion.r', script.dir.agregacion, fixed=T)[1])) && iFrame >= 0) {
  script.dir.agregacion <- sys.frame(iFrame)$ofile
  iFrame <- iFrame - 1
}
if (is.null(script.dir.agregacion)) { script.dir.agregacion <- ''
} else { script.dir.agregacion <- paste0(dirname(script.dir.agregacion), '/') }

source(paste0(script.dir.agregacion, '../instalarPaquetes/instant_pkgs.r'), encoding = 'WINDOWS-1252')
instant_pkgs(c('stats', 'sp', 'Rcpp', 'raster', 'rgdal'))
library(stats)
library(sp)
library(Rcpp)
library(raster)
library(rgdal)

source(paste0(script.dir.agregacion, '../GrADS/ReadGrADS.r'), encoding = 'WINDOWS-1252')
source(paste0(script.dir.agregacion, '../sysutils/sysutils.r'), encoding = 'WINDOWS-1252')
source(paste0(script.dir.agregacion, '../pathUtils/pathUtils.r'), encoding = 'WINDOWS-1252')

naSiTodosNAFuncSiNo <- function(x, func, ...) {
  x <- x[!is.na(x)]
  if (length(x)==0) { return (NA)
  } else { return (func(x, ...)) }
}

naSiPorcNAsFuncSiNo <- function(x, func, porcDatosNoNulos=0.75, ...) {
  iNoNA <- !is.na(x)
  if (sum(iNoNA) > length(x) * porcDatosNoNulos) { return(func(x[iNoNA], ...))
  } else { return(NA) }
}

ocurrenciaUmbral <- function(x, umbral=0, comparador='=') {
  # Retorna 1 si existe un valor que cumpla con umbral la relacion definida por comparador 
  # Retorna 0 sino
  i <- !is.na(x)
  if (any(i)) { 
    if (comparador == '=') { return (any(x[i] == umbral))
    } else if (comparador == '>') { return (any(x[i] > umbral))
    } else if (comparador == '<') { return (any(x[i] < umbral))
    } else if (comparador == '>=') { return (any(x[i] >= umbral))
    } else if (comparador == '<=') { return (any(x[i] <= umbral))
    } else { stop(paste0('Comparador desconocido ', comparador))
    }
  } else { return(FALSE) }
}

contarUmbral <- function(x, umbral=0, comparador='=') {
  # Retorna la cantidad de valores que cumplan con umbral la relacion definida por comparador 
  i <- !is.na(x)
  if (any(i)) { 
    if (comparador == '=') { return (sum(x[i] == umbral))
    } else if (comparador == '>') { return (sum(x[i] > umbral))
    } else if (comparador == '<') { return (sum(x[i] < umbral))
    } else if (comparador == '>=') { return (sum(x[i] >= umbral))
    } else if (comparador == '<=') { return (sum(x[i] <= umbral))
    } else { stop(paste0('Comparador desconocido ', comparador))
    }
  } else { return(0) }
}

porcentajeUmbral <- function(x, umbral=0, comparador='=') {
  # Retorna el porcentaje de valores que cumplan con umbral la relacion definida por comparador 
  i <- !is.na(x)
  if (any(i)) { 
    if (comparador == '=') { return (sum(x[i] == umbral) / sum(i) * 100)
    } else if (comparador == '>') { return (sum(x[i] > umbral) / sum(i) * 100)
    } else if (comparador == '<') { return (sum(x[i] < umbral) / sum(i) * 100)
    } else if (comparador == '>=') { return (sum(x[i] >= umbral) / sum(i) * 100)
    } else if (comparador == '<=') { return (sum(x[i] <= umbral) / sum(i) * 100)
    } else { stop(paste0('Comparador desconocido ', comparador))
    }
  } else { return(FALSE) }
}

parseFuncionAgregacion <- function(funcionAgregacion) {
  #x <- seq(1, 10, length.out=10)
  #funcionAgregacion <- 'max'
  #funcionAgregacion <- 'ocurrenciaUmbral;umbral=1;comparador=='
  #lala <- parseFuncionAgregacion(funcionAgregacion)
  #formals(lala)
  #lala(x)  
  
  paramsSplit <- unlist(strsplit(funcionAgregacion, ";"))
  funcion <- get(paramsSplit[1])
  paramsFuncion <- formals(funcion)
  if (!is.null(paramsFuncion)) {
    paramsSplit <- paramsSplit[-1]
    for (i in 1:length(paramsSplit)) {
      posIgual <- regexpr('=', paramsSplit[i], fixed=T)[1]
      nomParam <- substr(x=paramsSplit[i], start=1, stop=posIgual - 1)
      valParam <- type.convert(substr(x=paramsSplit[i], start=posIgual+1, stop=nchar(paramsSplit[i])), as.is=T)
      
      formals(funcion)[nomParam] <- valParam
    }
  }
  return (funcion)
}

agregar <- function(x, funcionAgregacion, claseIndiceI=rep(1, nrow(x)), ordenarPorClases=F, reducirDimensionesUnitarias=T, na.rm=T, ...) {
  # reduce la primer dimension de x de acuerdo a las clases en claseIndiceI usando funcionAgregacion
  # retorna una matriz length(unique(claseIndiceI)) x dim(x)[2:length(dim(x))] o length(unique(claseIndiceI)) x 1 si x es un vector
  # por ejemplo si x = (1, 2, 3, 4, 5, 6), claseIndiceI=(1,1, 2,2, 3,3) y funcionAgregacion=max
  # retorna una matriz de 3x1 = (2, 4, 6)
  # si x = (1, 2, 
  #         3, 4, 
  #         5, 6), claseIndiceI=(1, 1, 2) y funcionAgregacion=max 
  # retorna una matriz de 2x2 = (3, 4, 
  #                              5, 6)
  # el valor por defecto de claseIndiceI reduce la primer dimension a un único valor

  #funcionAgregacion <- max
  #claseIndiceI <- c(2,2, 1, 1, 1)
  if (is.character(funcionAgregacion))
    funcionAgregacion <- parseFuncionAgregacion(funcionAgregacion)
  
  if (na.rm) { res <- aggregate(x, by=list(claseIndiceI), FUN=naSiTodosNAFuncSiNo, func=funcionAgregacion, ...)
  } else { res <- aggregate(x, by=list(claseIndiceI), FUN=funcionAgregacion, ...) }
  res <- res[,-1]
  if (ordenarPorClases) {
    clases <- unique(claseIndiceI)
    res[order(clases),]
  }

  if (!is.null(dim(x))) {
    dimRes <- c(length(unique(claseIndiceI)), dim(x)[2:length(dim(x))])
    if (reducirDimensionesUnitarias) dimRes <- dimRes[dimRes > 1]

    if (length(dimRes) == 2) { res <- matrix(unlist(res), nrow=dimRes[1], ncol=dimRes[2])
    } else { res <- array(unlist(res), dim=dimRes) }
    
    if (!is.null(colnames(x)) & ncol(x) == ncol(res)) colnames(res) <- colnames(x)
  }
  
  return(res)
}

agregacionEspacialAPoligonos <- function(
    spObj, shpPoligonos, funcionAgregacion, zcol=1, na.rm=T, useRaster=FALSE) {
  if (is.character(funcionAgregacion)) {
    funcionAgregacion <- parseFuncionAgregacion(funcionAgregacion)
  }
  if (!identicalCRS(spObj, shpPoligonos)) {
    shpPoligonos <- spTransform(shpPoligonos, spObj@proj4string)
  }
  
  if (useRaster) {
    return (
      raster::extract(x = raster(spObj), y = shpPoligonos, fun=funcionAgregacion, na.rm=na.rm))
  } else {
    if (na.rm) {
      return (as.numeric(sp::over(
        x=shpPoligonos, y=spObj, fn=naSiTodosNAFuncSiNo, func=funcionAgregacion)[, zcol]))
    } else {
      return(sapply(vals, FUN = function(x) { funcionAgregacion(x[, zcol])}))
    }
  }
  
  
  #if (na.rm) {
  #  return (as.numeric(
  #    over(x=shpPoligonos, y=spObj[,1], fn=naSiTodosNAFuncSiNo, func=funcionAgregacion)))
  #} else {
  #  return(as.numeric(over(x=shpPoligonos, y=spObj[,1], fn=funcionAgregacion)))
  #}
}

agregacionEspacialAPoligonosDesdeArchivo <- function(
    pathSpObj, shpPoligonos, funcionAgregacion, zcol=1, na.rm=T, guardarCSV=FALSE,
    retornarResultados=TRUE, useRaster=FALSE) {
  # pathSpObj <- pathsSpObjs[1]
  spObj <- readGDAL(pathSpObj, silent = T)
  res <- agregacionEspacialAPoligonos(
    spObj=spObj, shpPoligonos=shpPoligonos, funcionAgregacion=funcionAgregacion, zcol=zcol, 
    na.rm=na.rm, useRaster=useRaster)
  
  if (guardarCSV) {
    write.table(x=matrix(res, ncol=1), 
                file=changeFileExt(pathSpObj, nuevaExtensionConPunto = '.csv'), sep = ',',
                row.names = F, col.names = F)
  }
  
  if (retornarResultados) {
    return(res)
  } else {
    return(NULL)
  }
}

agregacionEspacialAPoligonosDesdeArchivos <- function(
    pathsSpObjs, shpPoligonos, funcionAgregacion, zcol=1, na.rm=T, nCoresAUsar=0, 
    guardarCSV=FALSE, retornarResultados=TRUE, useRaster=FALSE) {
  if (nCoresAUsar <= 0) {
    nCoresAUsar <- min(getAvailableCores(maxCoresPerGB = 1), length(pathsSpObjs))
  }
  
  pathsSpObjs <- changeFileExt(listaMapas$nombreArchivo, '.tif')
  
  if (nCoresAUsar > 1) {
    cl <- makeCluster(getOption('cl.cores', nCoresAUsar))
    clusterExport(cl, varlist = c('script.dir.agregacion'))
    clusterEvalQ(cl = cl, expr = {
      source(paste0(script.dir.agregacion, 'agregacion.r'), encoding = 'WINDOWS-1252')
      if (exists(x = 'setMKLthreads')) { setMKLthreads(1) }
    })
    
    resultados <- parSapplyLB(
      cl=cl, X=pathsSpObjs, FUN=agregacionEspacialAPoligonosDesdeArchivo, shpPoligonos=shpPoligonos,
      funcionAgregacion=funcionAgregacion, zcol=zcol, na.rm=na.rm, guardarCSV=guardarCSV, 
      retornarResultados=retornarResultados, useRaster=useRaster)
    stopCluster(cl)
  } else {
    resultados <- sapply(
      X=pathsSpObjs, FUN=agregacionEspacialAPoligonosDesdeArchivo, shpPoligonos=shpPoligonos, 
      funcionAgregacion=funcionAgregacion, zcol=zcol, na.rm=na.rm, guardarCSV=guardarCSV, 
      retornarResultados=retornarResultados, useRaster=useRaster)
  }
  
  return(resultados)
}

testAgregacion <- function() {
  funcionAgregacion <- max
  
  x1 <- c(1, 2, NA, NA, 5, 6)
  claseIndiceI1 <- c(1,1, 2,2, 3,3)
  x1
  claseIndiceI1
  agregar(x1, funcionAgregacion, claseIndiceI1, na.rm=T)
  
  x2 <- matrix(data=c(1, 2, NA, 4, NA, 6), nrow=3, ncol=2, byrow=T)
  claseIndiceI2 <- c(1, 1, 2)
  x2
  claseIndiceI2
  agregar(x2, funcionAgregacion, claseIndiceI2, na.rm=T)
  
  x3 <- array(c(rnorm(30)), dim=c(3, 5, 2), dimnames=c('T', 'Y', 'X'))
  x3[1,1,1] <- NA
  x3[2,1,1] <- NA
  x3[3,1,1] <- NA
  claseIndiceI3 <- c(1, 1, 2)
  x3[1,,]
  x3[2,,]
  x3[3,,]
  claseIndiceI3
  res <- agregar(x3, funcionAgregacion, claseIndiceI3, na.rm=F)
  res[1,,]
  res[2,,] == x3[3,,]
  
  xdf <- data.frame(x=x1, y=x1+1, z=x1*2)
  claseIndiceIdf <- c(1, 2, 1, 2, 3, 3)
  xdf
  claseIndiceIdf
  agregar(xdf, funcionAgregacion, claseIndiceIdf)
}


agregacionTemporalGrillada_ti <- function(
  ti=1, fechas, pathsRegresor, nFechasAAgregar, minNfechasParaAgregar, funcionAgregacion, 
  formatoNomArchivoSalida, paramsCTL, shpBase, iOver, borrarOriginales, overlap, funcEscalado) {
  # fechas <- tempAireMin$fechas
  # pathsRegresor <- pathsRegresores[, 1]
  # formatoNomArchivoSalida <- paste0('Datos/MODIS/MOD11A1_LST_Day_3/MOD11A1_%Y-%m-%d.LST_Day_1km_', nFechasAAgregar, '.tif')
  # ti <- which(fechas==as.POSIXct('2002-09-24', tz=tz(fechas[1])))
  # ti <- tSeq[1]

  if (overlap) {
    tiMin <- max(1, ti - trunc(nFechasAAgregar / 2))
    tiMax <- min(length(pathsRegresor), ti + trunc(nFechasAAgregar / 2))
  } else {
    tiMin <- ti
    tiMax <- ti + nFechasAAgregar - 1
  }

  regresorTs <- vector(mode = "list", tiMax - tiMin + 1)
  n <- 1

  #i <- 1
  nPixeles <- 0
  for (i in tiMin:tiMax) {
    # Solo cargo los no nulos
    if (!is.na(pathsRegresor[i])) {
      if (!is.null(paramsCTL)) { 
        regresorTs[[n]] <- try(
          readXYGridSP(ctl = paramsCTL$ctl, dsetOverride = pathsRegresor[i], 
                       grillaXY = paramsCTL$grilla))
      } else { 
        regresorTs[[n]] <- try(readGDAL(pathsRegresor[i], silent=T))
      }
      if (!('try-error' %in% class(regresorTs[[n]]))) {
        if (!is.null(iOver)) {
          if (is(regresorTs[[n]], 'SpatialGridDataFrame')) {
            arrInd <- arrayInd(ind = iOver, .dim = regresorTs[[n]]@grid@cells.dim)
            xs <- unique(arrInd[,2])
            ys <- unique(arrInd[,1])
            # shpBaseAux <- spTransform(shpBase, regresorTs[[n]]@proj4string)
            # mapearGrillaGGPlot(regresorTs[[n]][xs, ys, 1], shpBase = shpBaseAux)
            
            regresorTs[[n]] <- regresorTs[[n]][xs, ys, 1]
          } else {
            regresorTs[[n]] <- regresorTs[[n]][iOver,]
          }
        }
        
        nPixeles <- nrow(regresorTs[[n]])
        n <- n + 1  
      }
    }
  }
  
  if (length(regresorTs) != n - 1) length(regresorTs) <- n - 1

  if (length(regresorTs) >= minNfechasParaAgregar && nPixeles > 0) {
    res <- regresorTs[[1]]
  
    # spplot(regresorTs[[1]])
    
    noEsNulo <- matrix(nrow = nPixeles, ncol = length(regresorTs))
    for (j in 1:length(regresorTs)) noEsNulo[, j] <- !is.na(regresorTs[[j]]@data[, 1])
    nNoNulos <- rowSums(noEsNulo)
    
    valsPixeles <- matrix(nrow = nPixeles, ncol = length(regresorTs))
    for (j in 1:length(regresorTs)) valsPixeles[, j] <- regresorTs[[j]]@data[, 1]
    
    if (!is.null(funcEscalado)) { valsPixeles <- funcEscalado(valsPixeles) }
    
    # Hardcodeo por performance
    if (identical(x = funcionAgregacion, y = base::mean)) { res@data[, 1] <- rowMeans(x = valsPixeles, na.rm = T) 
    } else if (identical(x = funcionAgregacion, y = base::sum)) { res@data[, 1] <- rowSums(x = valsPixeles, na.rm = T) 
    } else { res@data[, 1] <- apply(X = valsPixeles, MARGIN = 1, FUN = funcionAgregacion, na.rm=T) }
    res@data[nNoNulos < minNfechasParaAgregar, 1] <- NA
    
    # spplot(res)
    
    nomArch <- format(x = fechas[ti], formatoNomArchivoSalida)
    writeGDAL(dataset = res, fname = nomArch, options = c('COMPRESS=DEFLATE', 'PREDICTOR=2', 'ZLEVEL=9'))
    
    return(nomArch)
  } else { return(NA)}
}

agregacionTemporalGrillada <- function(
    fechas, pathsRegresor, formatoNomArchivoSalida=paste0('%.4d-%.2d-%.2d_', nFechasAAgregar, '.tif'), 
    nFechasAAgregar=3, minNfechasParaAgregar=max(trunc(nFechasAAgregar/2), 1), tIni=1, 
    tFin=length(pathsRegresor), funcionAgregacion=base::mean, ctl=NULL, shpBase=NULL,
    borrarOriginales=FALSE, overlap=TRUE, funcEscalado=NULL, nCoresAUsar=0) {
  # Para calcular agregaciones temporales de una serie temporal de un mismo regresor
  # pathsRegresor es una vector de rasters
  # Para cada fecha fi, se toman los píxeles de las fechas entre fi-trunc(nFechasAAgregar/2) y fi+trunc(nFechasAAgregar/2) y se
  # calcula funcionAgregacion con ellos
  # Si no hay al menos minNfechasParaAgregar píxeles disponibles el píxel se devuelve nulo
  if (nCoresAUsar <= 0) {
    nCoresAUsar <- min(getAvailableCores(maxCoresPerGB = 1), tFin - tIni + 1)
  }
  
  dir.create(dirname(formatoNomArchivoSalida), showWarnings = F, recursive = T)
  
  if (overlap) { tSeq <- tIni:tFin
  } else { tSeq <- seq.int(tIni, tFin, by=nFechasAAgregar) }
  
  if (!is.null(ctl)) { paramsCTL <- list(ctl=ctl, grilla=getGrillaNativa(ctl))
  } else { paramsCTL <- NULL }
  
  if (!is.null(shpBase)) {
    iAux <- 2
    iAux <- which.min(!is.na(pathsRegresor))
    if (!is.null(paramsCTL)) {
      regresorAux <- try(readXYGridSP(ctl = paramsCTL$ctl, dsetOverride = pathsRegresor[iAux],
                                      grillaXY = paramsCTL$grilla))
    } else {
      regresorAux <- try(readGDAL(pathsRegresor[iAux], silent=T))
    }
    shpBaseAux <- spTransform(shpBase, regresorAux@proj4string)
    bbaux <- getPoligonoBoundingBox(objSP = shpBaseAux, factorExtensionX = 1.1)
    iOver <- which(!is.na(over(regresorAux, bbaux)))
    rm(iAux, regresorAux, shpBaseAux, bbaux)
  }

  if (nCoresAUsar > 1) {
    cl <- makeCluster(getOption('cl.cores', nCoresAUsar))
    clusterExport(cl, varlist = c('script.dir.agregacion'))
    clusterEvalQ(cl = cl, expr = {
      require('rgdal')
      source(paste0(script.dir.agregacion, '../GrADS/ReadGrADS.r'), encoding = 'WINDOWS-1252')
      source(paste0(script.dir.agregacion, '../interpolar/interpolarEx.r'), encoding = 'WINDOWS-1252')
      if (exists(x = 'setMKLthreads')) { setMKLthreads(1) }
    })

    parSapplyLB(cl=cl, X=tSeq, FUN=agregacionTemporalGrillada_ti,
                fechas=fechas, pathsRegresor=pathsRegresor, nFechasAAgregar=nFechasAAgregar,
                minNfechasParaAgregar=minNfechasParaAgregar, funcionAgregacion=funcionAgregacion,
                formatoNomArchivoSalida=formatoNomArchivoSalida, paramsCTL=paramsCTL,
                shpBase=shpBase, iOver=iOver, borrarOriginales=borrarOriginales, overlap=overlap, 
                funcEscalado=funcEscalado)
    stopCluster(cl)
  } else {
    sapply(X=tSeq, FUN=agregacionTemporalGrillada_ti,
           fechas=fechas, pathsRegresor=pathsRegresor, nFechasAAgregar=nFechasAAgregar,
           minNfechasParaAgregar=minNfechasParaAgregar, funcionAgregacion=funcionAgregacion,
           formatoNomArchivoSalida=formatoNomArchivoSalida, paramsCTL=paramsCTL, shpBase=shpBase,
           iOver=iOver, borrarOriginales=borrarOriginales, overlap=overlap, funcEscalado=funcEscalado)
  }
}

agregacionTemporalGrillada2_ti <- function(ti=1, fechas, pathsRegresores, nFechasAAgregar=1, minNfechasParaAgregar=max(minNfechasParaAgregar=trunc(nFechasAAgregar/2),1), 
                                           funcionAgregacion=mean, formatoNomArchivoSalida=paste0('%Y-%m-%d_', nFechasAAgregar, '.tif')) {
  # fechas <- tempAireMin$fechas
  # pathsRegresor <- pathsRegresores[, 1]
  # formatoNomArchivoSalida <- paste0('Datos/MODIS/MOD11A1_LST_Day_3/MOD11A1_%Y-%m-%d.LST_Day_1km_', nFechasAAgregar, '.tif')
  
  tiMin <- max(1, ti - trunc(nFechasAAgregar / 2))
  tiMax <- min(nrow(pathsRegresores), ti + trunc(nFechasAAgregar / 2))
  
  regresores <- vector(mode = "list", (tiMax - tiMin + 1) * ncol(pathsRegresores))
  
  n <- 1
  
  nPixeles <- 0
  for (i in tiMin:tiMax) {
    for (j in 1:ncol(pathsRegresores)) {	
      # Solo cargo los no nulos
      if (!is.na(pathsRegresores[i, j])) {
        regresores[[n]] <- try(readGDAL(pathsRegresores[i, j], silent=T))
        if (!('try-error' %in% class(regresores[[n]]))) { 
          nPixeles <- nrow(regresores[[n]])
          n <- n + 1  
        }
      }
    }
  }
  if (length(regresores) != n - 1) length(regresores) <- n - 1
  
  if (length(regresores) >= minNfechasParaAgregar && nPixeles > 0) {
    res <- regresores[[1]]
    
    noEsNulo <- matrix(nrow = nPixeles, ncol = length(regresores))
    for (j in 1:length(regresores)) noEsNulo[, j] <- !is.na(regresores[[j]]@data[, 1])
    nNoNulos <- rowSums(noEsNulo)
    
    valsPixeles <- matrix(nrow = nPixeles, ncol = length(regresores))
    for (j in 1:length(regresores)) valsPixeles[, j] <- regresores[[j]]@data[, 1]
    
    # Hardcodeo por performance
    if (identical(x = funcionAgregacion, y = base::mean)) { res@data[, 1] <- rowMeans(x = valsPixeles, na.rm = T) 
    } else if (identical(x = funcionAgregacion, y = base::sum)) { res@data[, 1] <- rowSums(x = valsPixeles, na.rm = T) 
    } else { res@data[, 1] <- apply(X = valsPixeles, MARGIN = 1, FUN = funcionAgregacion, na.rm=T) }
    res@data[nNoNulos < minNfechasParaAgregar, 1] <- NA
    
    source(paste0(script.dir.agregacion, '../pathUtils/pathUtils.r'), encoding = 'WINDOWS-1252')
    
    nomArch <- format(x = fechas[ti], formatoNomArchivoSalida)
    writeGDAL(dataset = res, fname = nomArch, options = c('COMPRESS=DEFLATE', 'PREDICTOR=2', 'ZLEVEL=9'))
    
    return(nomArch)
  } else { return(NA) }
}

agregacionTemporalGrillada2 <- function(fechas, pathsRegresores, formatoNomArchivoSalida=paste0('%.4d-%.2d-%.2d_', nFechasAAgregar, '.tif'), 
                                        nFechasAAgregar=1, minNfechasParaAgregar=max(trunc(nFechasAAgregar/2), 1),
                                        tIni=1, tFin=nrow(pathsRegresores), funcionAgregacion=base::mean) {
  # Para calcular agregaciones temporales combinando varios rasters de origen
  # pathsRegresores es una matriz de rasters
  # Para cada fecha fi, se toman los píxeles de las fechas entre fi-trunc(nFechasAAgregar/2) y fi+trunc(nFechasAAgregar/2) de 
  # todos los regresores y se calcula funcionAgregacion con ellos
  # Si no hay al menos minNfechasParaAgregar píxeles disponibles el píxel se devuelve nulo
  nCoresAUsar <- min(getAvailableCores(maxCoresPerGB = 1), tFin - tIni + 1)
  dir.create(dirname(formatoNomArchivoSalida), showWarnings = F, recursive = T) 
  
  if (nCoresAUsar > 1) {
    cl <- makeCluster(getOption('cl.cores', nCoresAUsar))
    clusterExport(cl, varlist = c('script.dir.agregacion'))
    clusterEvalQ(cl = cl, expr = require('rgdal'))
    if (exists(x = 'setMKLthreads')) { clusterEvalQ(cl = cl, expr = setMKLthreads(1)) }
    parSapplyLB(cl=cl, X=tIni:tFin, FUN=agregacionTemporalGrillada2_ti,
                fechas=fechas, pathsRegresores=pathsRegresores, nFechasAAgregar=nFechasAAgregar, 
                minNfechasParaAgregar=minNfechasParaAgregar, funcionAgregacion=funcionAgregacion, 
                formatoNomArchivoSalida=formatoNomArchivoSalida)
    stopCluster(cl)
  } else {
    sapply(X=tIni:tFin, FUN=agregacionTemporalGrillada2_ti,
           fechas=fechas, pathsRegresores=pathsRegresores, nFechasAAgregar=nFechasAAgregar, 
           minNfechasParaAgregar=minNfechasParaAgregar, funcionAgregacion=funcionAgregacion, 
           formatoNomArchivoSalida=formatoNomArchivoSalida)
    
    #for (i in tIni:tFin) {
    #  print(i)
    #  agregacionTemporalGrillada_ti(ti=i, fechas=fechas, pathsRegresor=pathsRegresor, nFechasAAgregar=nFechasAAgregar, 
    #                                minNfechasParaAgregar=minNfechasParaAgregar, funcionAgregacion=funcionAgregacion, 
    #                                formatoNomArchivoSalida=formatoNomArchivoSalida)
    #}
  }
}

agregacionTemporalGrillada3_claseI <- function(iClase=1, fechas, pathsRegresor, claseFechaI=1:length(pathsRegresor), clases=sort(unique(claseFechaI)),
                                               nomArchivosSalidaClaseI=paste0(clases, '.tif'), minNfechasParaAgregar=0, 
                                               funcionAgregacion=base::mean, interpolarFaltantes='No', overlap=0, 
                                               pathShpMask=NULL, proj4stringShpMask=NULL, spSinMascara=NULL, recalcularSiYaExiste=T, ...) {
  #iClase <- 6
  if (recalcularSiYaExiste || !file.exists(nomArchivosSalidaClaseI[iClase]) || file.info(nomArchivosSalidaClaseI[iClase])$size <= 0) {
    print(iClase)
    
    iClaseIni <- ((iClase-overlap-1) %% length(clases) + 1)
    iClaseFin <- ((iClase+overlap-1) %% length(clases) + 1)
    if (iClaseIni < iClaseFin) { iClases <- iClaseIni:iClaseFin
    } else { iClases <- c(iClaseIni:length(clases),1:iClaseFin) }
    
    pathsClaseI <- pathsRegresor[claseFechaI %in% iClases]
    pathsClaseI <- pathsClaseI[!is.na(pathsClaseI)]
    
    regresores <- raster::stack(pathsClaseI)
    # Esto es para hacer los cálculos en bloques y no exceder la RAM pero no está funcionando. Revisar después.
    #funcionAgregacion<-mean
    #blockFunction <- function(x, minNfechasParaAgregar, funcionAgregacion) {
    #  res <- funcionAgregacion(x, na.rm = T)
    #  if (minNfechasParaAgregar > 0) {
    #    res[sum(!is.na(x)) < minNfechasParaAgregar] <- NA
    #  }
    #  res
    #}
    #calc(x = regresores, fun = blockFunction, minNfechasParaAgregar=minNfechasParaAgregar, funcionAgregacion=funcionAgregacion)
    res <- calc(x=regresores, fun=funcionAgregacion, na.rm=T, ...=...)
    #res <- calc(x=regresores, fun=funcionAgregacion, na.rm=T)
    res[is.nan(res)] <- NA
    
    if (minNfechasParaAgregar > 0) { res[sum(!is.na(regresores)) < minNfechasParaAgregar] <- NA }
    # plot(res)
    
    # Si hay pathShpMask se interpolan todos los píxeles y solo se conservan los píxeles dentro de el shapefile en el
    # Si no, si hay spSinMascara se interpolan todos los píxeles que estén dentro de spSinMascara pero se conservan los que estén afuera también
    # Sino se interpola toda el área
    iNA <- is.na(getValues(res))
    if (!is.null(pathShpMask) && file.exists(pathShpMask)) {
      grilla <- as(res, 'SpatialGrid')
      source(paste0(script.dir.agregacion, 'interpolarEx.r'), encoding = 'WINDOWS-1252')
      shpMask <- cargarSHPYObtenerMascaraParaGrilla(pathSHP = pathShpMask, proj4strSHP = proj4stringShpMask, grilla = grilla, spSinMascara = spSinMascara)
      mascara <- shpMask$mask
    } else if (!is.null(spSinMascara)) {
      grilla <- as(res, 'SpatialGrid')
      mascara <- !iNA
      if (proj4string(grilla) != proj4string(spSinMascara)) {
        spSinMascara <- spTransform(spSinMascara, grilla@proj4string)
      }
      if (any(grepl(pattern = 'SpatialPoints', x = class(spSinMascara)))) { 
        mascara[over(spSinMascara, grilla)] <- TRUE
      } else { 
        mascara[!is.na(over(grilla, geometry(spSinMascara)))] <- TRUE 
      }
    } else {
      mascara <- T
    }
    
    if (interpolarFaltantes != 'No' && any(iNA & mascara)) {
      spCoords <- SpatialPoints(coordinates(res), proj4string = res@proj4string)
      spRes <- SpatialPixelsDataFrame(points = spCoords, data = data.frame(value=getValues(res)))
      paramsRellenarSP <- createParamsInterpolarYMapear(coordsAInterpolarSonGrilla = gridded(spRes), interpolationMethod = interpolarFaltantes,
                                                        nmax=20, inverseDistancePower = 2)
      spRes <- as(rellenarSP(sp = spRes, mascara = mascara, metodo = 'idw', nMuestras = sum(!iNA), nRepeticiones = 1, 
                             zcol = 1, nCuadrantesX = 1, nCuadrantesY = 1, nCuadrantesZ = 1, params=paramsRellenarSP), 'SpatialGridDataFrame')
    } else {
      spCoords <- SpatialPoints(coordinates(res), proj4string = res@proj4string)
      spRes <- as(SpatialPixelsDataFrame(points = spCoords, data = data.frame(value=getValues(res))), 'SpatialGridDataFrame')
    }
    spRes@data[!mascara,] <- NA
    # vals <- getValues(res)
    # vals[!mascara] <- NA
    # res <- setValues(res, vals)
    # plot(res)
    
    writeGDAL(dataset = spRes, fname = nomArchivosSalidaClaseI[iClase], options = c('COMPRESS=DEFLATE', 'PREDICTOR=2', 'ZLEVEL=9'))
    #writeRaster(x = res, filename = nomArchivosSalidaClaseI[iClase], overwrite=T, options = c('COMPRESS=DEFLATE', 'PREDICTOR=2', 'ZLEVEL=9'))
  }
}

agregacionTemporalGrillada3 <- function(
    fechas, pathsRegresor, claseFechaI=1:length(pathsRegresor), clases=sort(unique(claseFechaI)),
    nomArchivosSalidaClaseI=paste0(clases, '.tif'), minNfechasParaAgregar=0, 
    funcionAgregacion=base::mean, interpolarFaltantes='No', overlap=0, pathShpMask=NULL, 
    proj4stringShpMask=NULL, spSinMascara=NULL, recalcularSiYaExiste=T, nCoresAUsar=0, ...) {
  # Para calcular climatologías.
  # Obtiene el conjunto de fechas F que tengan la misma claseFechaI que la fecha fi a calcular
  # Calcula los píxeles de fi y aplicando funcionAgregacion a los mismos píxelespero para las fechas en F
  # Si no hay al menos minNfechasParaAgregar píxeles en F con datos no nulos el píxel será nulo
  # Si se indica overlap > 0 se toman las clases aledañas a las de fi
  # Se puede pasar en interpolarFaltantes el nombre de alguno de los métodos de interpolación para rellenar los píxeles faltantes a la salida
  
  #pathsRegresor <- dir('D:/Tesis/Datos/LST_Night_Combinada_Filtrada', pattern = '*.tif$', full.names = T)
  
  # pathsRegresor <- pathsRegresores[, colnames(pathsRegresores) == 'LST_Night_Combinada']
  #fechas <- fechasObservaciones
  #claseFechaI <- yday(fechas)
  #claseFechaI[claseFechaI == 366] <- 365
  #clases=sort(unique(claseFechaI))
  #nomArchivosSalidaClaseI=paste0(pathDatos, 'LST_Night_Combinada_Clim_mean/', clases, '.tif')
  #minNfechasParaAgregar <- 5
  #funcionAgregacion=base::mean
  #interpolarFaltantes='idw'
  #overlap <- 1
  #pathShpMask=NULL
  #proj4stringShpMask=NULL
  #spSinMascara <- shpMask$shp
  
  if (nCoresAUsar <= 0) {
    nCoresAUsar <- min(getAvailableCores(maxCoresPerGB = 1), length(clases))
  }
  dir.create(unique(dirname(nomArchivosSalidaClaseI)), showWarnings = F, recursive = T) 
  
  if (nCoresAUsar > 1) {
    cl <- makeCluster(getOption('cl.cores', nCoresAUsar))
    clusterExport(cl, varlist = c('script.dir.agregacion'))
    if (exists(x = 'setMKLthreads')) { clusterEvalQ(cl = cl, expr = setMKLthreads(1)) }
    clusterEvalQ(cl, expr = {
      require('raster')
      require('sp')
      require('rgdal')
      source(paste0(script.dir.agregacion, '../interpolar/interpolarEx.r', encoding = 'WINDOWS-1252'))
    })
    parSapplyLB(
      cl=cl, X=1:length(clases), FUN=agregacionTemporalGrillada3_claseI, fechas=fechas, 
      pathsRegresor=pathsRegresor, claseFechaI=claseFechaI, clases=clases, 
      nomArchivosSalidaClaseI=nomArchivosSalidaClaseI, minNfechasParaAgregar=minNfechasParaAgregar, 
      funcionAgregacion=funcionAgregacion, interpolarFaltantes=interpolarFaltantes, overlap=overlap, 
      pathShpMask=pathShpMask, proj4stringShpMask=proj4stringShpMask, spSinMascara=spSinMascara, 
      recalcularSiYaExiste=recalcularSiYaExiste, ...=...)
    stopCluster(cl)
  } else {
    sapply(
      X=1:length(clases), FUN=agregacionTemporalGrillada3_claseI, fechas=fechas, 
      pathsRegresor=pathsRegresor, claseFechaI=claseFechaI, clases=clases, 
      nomArchivosSalidaClaseI=nomArchivosSalidaClaseI, minNfechasParaAgregar=minNfechasParaAgregar, 
      funcionAgregacion=funcionAgregacion, interpolarFaltantes=interpolarFaltantes, overlap=overlap,
      pathShpMask=pathShpMask, proj4stringShpMask=proj4stringShpMask, spSinMascara=spSinMascara, 
      recalcularSiYaExiste=recalcularSiYaExiste, ...=...)
    
    #for (i in 1:length(clases)) {
    #  agregacionTemporalGrillada3_claseI(iClase = i, fechas=fechas, pathsRegresor=pathsRegresor, claseFechaI=claseFechaI, clases=clases,
    #                                nomArchivosSalidaClaseI=nomArchivosSalidaClaseI, minNfechasParaAgregar=minNfechasParaAgregar, 
    #                                funcionAgregacion=funcionAgregacion, interpolarFaltantes=interpolarFaltantes, overlap=overlap,
    #                                pathShpMask=pathShpMask, proj4stringShpMask=proj4stringShpMask, spSinMascara=spSinMascara)
    #}
  }
  return(nomArchivosSalidaClaseI)
}