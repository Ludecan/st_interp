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
if (iFrame >= 3) { script.dir.GetMODISEx <- sys.frame(iFrame - 3)$ofile 
} else { script.dir.GetMODISEx <- NULL }
while ((is.null(script.dir.GetMODISEx) || is.na(regexpr('mapearEx.r', script.dir.GetMODISEx, fixed=T)[1])) && iFrame >= 0) {
  script.dir.GetMODISEx <- sys.frame(iFrame)$ofile
  iFrame <- iFrame - 1
}
if (is.null(script.dir.GetMODISEx)) { script.dir.GetMODISEx <- ''
} else { script.dir.GetMODISEx <- paste0(dirname(script.dir.GetMODISEx), '/') }

# check for uninstalled dependencies, install them and load them
source(paste0(script.dir.GetMODISEx, '../parsearParams/parsearParamsUtils.r'), encoding = 'WINDOWS-1252')
source(paste0(script.dir.GetMODISEx, '../instalarPaquetes/instant_pkgs.r'), encoding = 'WINDOWS-1252')
source(paste0(script.dir.GetMODISEx, '../pathUtils/pathUtils.r'), encoding = 'WINDOWS-1252')
source(paste0(script.dir.GetMODISEx, '../cacheFunciones/cacheFunciones.r'), encoding = 'WINDOWS-1252')
source(paste0(script.dir.GetMODISEx, '../tryUtils/tryUtils.r'), encoding = 'WINDOWS-1252')
instant_pkgs(c('digest', 'sp', 'rgdal', 'Rcpp', 'raster', 'RCurl', 'lubridate', 'doParallel', 'rts'))
# source(paste(script.dir.GetMODISEx, 'ModisDownload.r', sep=''))
#shpGrillaMODISSinusoidalV5 <- readOGR('modis_sinusoidal', 'modis_sinusoidal_grid_world')
#plot(shpGrillaMODISSinusoidalV5)  
#shpGrillaMODISSinusoidalV5

createParamsGetMODIS <- function(pathEjecucion, pathProceso='./', producto='MOD09Q1',
                                 # coordenadas límite en la proyección de salida del area a recortar
                                 xMin, xMax, yMin, yMax,
                                 # baja todos los datos entre estas dos fechas en la periodicidad del producto. El separador debe ser '.'
                                 fechaIniUTC='2011.05.01', fechaFinUTC='2011.05.01',
                                 version='006', 
                                 # las bandas a obtener del dataset, se deben indicar con 0 las bandas a excluir y con 1 las bandas a incluir. 
                                 # Si se especifican menos bandas que las que hay en el dataset las restantes se asumen 0
                                 bandas='0 0 1',
                                 # ruta a la carpeta bin dentro del directorio de instalación de MRT
                                 MRTBinPath='C:/MRT/bin',
                                 # string en formato proj4 de la proyección de salida de los datos
                                 proj4stringResultados,
                                 escala=1, offset=0, minRangoValido=NA, maxRangoValido=NA,
                                 pathArchivosResultado='./') {
  # funcion auxiliar para saber todos los parámetros que se tienen que pasar en params y valores
  # por defecto de algunos parámetros
  return(list(pathEjecucion=pathEjecucion,
              pathProceso=pathProceso,
              producto=producto,
              xMin=xMin, xMax=xMax, yMin=yMin, yMax=yMax,
              fechaIniUTC=fechaIniUTC, fechaFinUTC=fechaFinUTC,
              version=version, bandas=bandas,
              MRTBinPath=MRTBinPath,
              proj4stringResultados=proj4stringResultados,
              escala=escala, offset=offset, minRangoValido=minRangoValido, maxRangoValido=maxRangoValido,
              pathArchivosResultado=pathArchivosResultado))
}

parsearParamsGetMODIS <- function(params) {
  return(getParamValuesFromConstructorParams(params, funcCrearParams=createParamsGetMODIS))
}

xy2MODISTile <- function(x, y, CellSize=926.62543305, TileSize=1200 * 926.62543305, 
                         ULx=-20015109.354, ULy=10007554.677) {
  # (x, y) debe estar en la proyección sinusoidal de la grilla de modis
  # Los valores por defecto de los parámetros son para la grilla sinusoidal v5 de MODIS
  v <- trunc(-(y - ULy) / TileSize);
  h <- trunc((x - ULx) / TileSize);
  return (list(h=h, v=v))
}

getHsYVsBoundingBox <- function(
    proj4stringResultados, SRS_stringResultados, xMin, xMax, yMin, yMax) {
  #xMin <- 0
  #xMax <- 1
  #yMin <- 2
  #yMax <- 3
  #proj4stringResultados <- '+proj=longlat +datum=WGS84'
  
  # Upper Left, Upper Right, Bottom Right, Bottom Left
  coords <- matrix(data = c(xMin, yMax, xMax, yMax, xMax, yMin, xMin, yMin), 
                   nrow = 4, ncol = 2, byrow = T)
  esquinas <- SpatialPoints(
    coords = coords, 
    proj4string = CRS(projargs = proj4stringResultados, SRS_string = SRS_stringResultados))
  
  # spTransform(esquinas, CRS('+proj=longlat +datum=WGS84', SRS_string="EPSG:4326"))
  p4strMODISSinusoidalv5 <- '+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs'
  SRS_MODISSinusoidalv5 <- "SR-ORG:6842"
  esquinas <- spTransform(
    esquinas, CRS(projargs = p4strMODISSinusoidalv5, SRS_string = SRS_MODISSinusoidalv5))
  
  #esquinas <- spTransform(esquinas, shpGrillaMODISSinusoidalV5@proj4string)
  #over(esquinas, shpGrillaMODISSinusoidalV5)
  
  HyV <- xy2MODISTile(coordinates(esquinas)[1,1], coordinates(esquinas)[1,2])
  minH <- HyV$h
  maxH <- minH
  minV <- HyV$v
  maxV <- HyV$v
  for (i in 2:4) {
    HyV <- xy2MODISTile(coordinates(esquinas)[i,1], coordinates(esquinas)[i,2])
    if (HyV$h < minH) { minH <- HyV$h
    } else if (HyV$h > maxH) { maxH <- HyV$h }
    if (HyV$v < minV) { minV <- HyV$v
    } else if (HyV$v > maxV) { maxV <- HyV$v }
  }
  
  hs <- seq.int(from=minH, to=maxH, by=1)
  vs <- seq.int(from=minV, to=maxV, by=1)
  return(list(hs=hs, vs=vs))
}

addDate_temporalResolution <- function(date, temporalResolution) {
  #date <- ymd(date)
  if (temporalResolution$units=='min') {
    minute(date) <- minute(date) + temporalResolution$increment
  } else if (temporalResolution$units=='day') {
    day(date) <- day(date) + temporalResolution$increment
  } else if (temporalResolution$units=='month') {
    month(date) <- month(date) + temporalResolution$increment
  } else if (temporalResolution$units=='year') {
    year(date) <- year(date) + temporalResolution$increment
  } else {
    stop(paste('ModisDownload.subDate_temporalResolution: unknown unit: ', temporalResolution$units, sep=''))
  }
  
  return(date)
}

getMODISEx <- function(
    producto='MOD09Q1',
    # coordenadas límite en la proyección de salida del area a recortar
    xMin, xMax, yMin, yMax,
    # baja todos los datos entre estas dos fechas en la periodicidad del producto. El separador debe ser '.'
    fechaIniUTC='2014.09.06', fechaFinUTC='2014.09.06',
    # las bandas a obtener del dataset, se deben indicar con 0 las bandas a excluir y con 1 las bandas a incluir. 
    # Si se especifican menos bandas que las que hay en el dataset las restantes se asumen 0
    bandas='0 0 1',
    # versión del producto MODIS a descargar
    version='005',                       
    # ruta a la carpeta bin dentro del directorio de instalación de MRT
    MRTBinPath='C:/MRT/bin',
    # string en formato proj4 de la proyección de salida de los datos
    proj4stringResultados,
    SRS_stringResultados,
    pathArchivosResultado='./',
    escala=0.02, offset=-273.15, graficar=F, forceReDownload=F, comprimirGeoTiffs=T) {
  HsYVs <- getHsYVsBoundingBox(proj4stringResultados, SRS_stringResultados, xMin, xMax, yMin, yMax)
  
  tokensP4str <- getTokens(proj4stringResultados, separadorTokens=' ')
  proyeccion <- toupper(getParamValue(tokensP4str, paramName='+proj', classOfValue='character', obligatorio=T))
  if (proyeccion == 'UTM') {
    utm_zone <- getParamValue(tokensP4str, paramName='+zone', classOfValue='character', obligatorio=T)
    proj_params <- '0 0 0 0 0 0 0 0 0 0 0 0 0 0 0' 
    if ('+south' %in% tokensP4str[,2]) {
      # La opción +south no es soportada por MRT, pero lo que hace es sumar un false_northing de 10000000 
      # a la coordenada sin la opción south, por lo tanto le restamos el false northing y obtenemos
      # la correspondiente coordenada sin la opción
      southing <- -10000000 
      utmSur <- TRUE
    } else {
      southing <- 0
      utmSur <- FALSE
    }
  } else {
    utmSur <- FALSE
    southing <- 0
    stop(paste('GetMODISEx.getMODISEx: proyección ', proyeccion, ' no implementada.', sep=''))
  }
  datum <- getParamValue(tokensP4str, paramName='+datum', classOfValue='character', obligatorio=F, paramDefaultValue='WGS84')
  doMosaic <- length(HsYVs$hs) > 1 || length(HsYVs$vs) > 1

  pathProductoMODIS <- subirDirectorio(pathArchivosResultado)
  dir.create(pathArchivosResultado, showWarnings = F, recursive = T)
  oldWD <- getwd()
  setwd(pathProductoMODIS)
  unlink('resample.log')

  if (interactive()) {
    x=producto
    h=HsYVs$hs
    v=HsYVs$vs
    dates=c(fechaIniUTC,fechaFinUTC)
    MRTpath=MRTBinPath
    mosaic=doMosaic
    proj=T
    UL=c(xMin, yMax+southing)
    LR=c(xMax, yMin+southing)
    bands_subset=bandas 
    proj_type=proyeccion
    ncore='auto'
  }
  
  # setwd('E:/mch/datosEnGrilla/fdg/MOD11A1/006')
  system.time(
  evaluarConReintentos(
  ModisDownload(x=producto, h=HsYVs$hs, v=HsYVs$vs, dates=c(fechaIniUTC,fechaFinUTC), version=version, MRTpath=MRTBinPath,
                mosaic=doMosaic, proj=T, UL=c(xMin, yMax+southing), LR=c(xMax, yMin+southing), bands_subset=bandas, 
                proj_type=proyeccion, proj_params=proj_params, utm_zone=utm_zone, datum=datum, forceReDownload=forceReDownload, 
                ncore='auto'))
  )
  
  # Busco en el log las lineas con los nombres de archivos de salida.
  lineasLog <- readLines('resample.log')
  archivos <- lineasLog[grep(lineasLog, pattern = 'output_filename:', fixed=T)]
  archivos <- sub(x = archivos, pattern = 'output_filename:[[:blank:]]*', replacement='')
  archivos <- unique(archivos)

  iBandas <- as.integer(unlist(strsplit(bandas, split = ' ', fixed = T)))
  nBandas <- sum(iBandas)                     
  iBandas <- which(iBandas == 1)
  
  iLineaDescBandas <- grep(lineasLog, pattern = 'band[[:blank:]]*select[[:blank:]]*type[[:blank:]]*lines[[:blank:]]*smpls[[:blank:]]*pixsiz[[:blank:]]*min[[:blank:]]*max[[:blank:]]*fill')[1]
  
  if (doMosaic) {
    # El mosaic tool ya hace la extracción de las bandas así que si se hace mosaic solo quedan las bandas que interesan
    lineasBandas <- strsplit(lineasLog[iLineaDescBandas + 1:nBandas], split = '[[:blank:]]+')
    unlink('Mosaic_*.hdf')
  } else {
    lineasBandas <- strsplit(lineasLog[iLineaDescBandas + iBandas], split = '[[:blank:]]+')
  }

  nombresBandas <- array('', dim = length(lineasBandas))
  minBandas <- array(NA, dim = length(lineasBandas))
  maxBandas <- array(NA, dim = length(lineasBandas))
  for (i in 1:length(lineasBandas)) {
    nombresBandas[i] <- lineasBandas[[i]][3]
    minBandas[i] <- as.numeric(lineasBandas[[i]][9])
    maxBandas[i] <- as.numeric(lineasBandas[[i]][10])
  }
  rm(lineasLog, iBandas, iLineaDescBandas, lineasBandas)
  # Le agrego a los archivos de salida el nombre de la banda para obtener el nombre de archivo final
  archivos <- paste(sub(archivos, pattern = '.tif', replacement = ''), '.', nombresBandas, '.tif', sep = '')
  file.rename(archivos, paste(pathArchivosResultado, archivos, sep=''))
  archivos <- paste(pathArchivosResultado, sort(archivos), sep='')
  
  # No se usa porque hay algún bug en raster que a veces guarda corruptos los archivos comprimidos
  #postProcesarArchivoMODIS <- function(i, archivos, nombresBandas, minBandas, maxBandas, escala, offset, utmSur, proj4stringResultados) {
  #  require('raster')
  #  campo <- raster(archivos[[i]])
  #  
  #  iBanda <- 1
  #  while (iBanda <= length(nombresBandas) & !(grepl(pattern = nombresBandas[iBanda], x = archivos[[i]]))) iBanda <- iBanda + 1
  #  
  #  valores <- getValues(campo)
  #  if (all(is.integer(valores))) {
  #    valores[valores < minBandas[iBanda] | valores > maxBandas[iBanda]] <- NA
  #    iNoNA <- !is.na(valores)
  #    valores[iNoNA] <- valores[iNoNA] * escala + offset
  #    campo <- setValues(campo, valores)
  #    if (utmSur) campo <- projectRaster(campo, crs = CRS(proj4stringResultados))
  #    writeRaster(x=campo, filename=archivos[[i]], overwrite=T, options = c('COMPRESS=DEFLATE', 'PREDICTOR=2', 'ZLEVEL=9'))
  #    #writeRaster(x=campo, filename='D:/lele.tif', overwrite=T, options = c('COMPRESS=DEFLATE', 'PREDICTOR=2', 'ZLEVEL=9'))
  #    #campo <- raster('//192.168.1.223/mch/datosEnGrilla/fdg/MYD11A1/MYD11A1_LST_Night/MYD11A1_2002-07-08.LST_Night_1km.tif')
  #  }
  #}
  
  i <- 1
  postProcesarArchivoMODIS_V2 <- function(i, archivos, nombresBandas, minBandas, maxBandas, escala, offset, utmSur, proj4stringResultados, comprimirGeoTiffs) {
    print(i)
    #lockFile <- paste(tempdir(), '/flock_GDAL.lock', sep='')
    #lock(lockFile, exclusive = FALSE)
    campo <- readGDAL(archivos[[i]], silent = T)
    #unlock(lockFile)
    if (is.integer(campo@data[,])) {
      iBanda <- 1
      while (iBanda <= length(nombresBandas) & !(grepl(pattern = nombresBandas[iBanda], x = archivos[[i]]))) iBanda <- iBanda + 1
      
      iNoNA <- !is.na(campo@data[,1])
      campo@data[iNoNA & (campo@data[,1] < minBandas[iBanda] | campo@data[,1] > maxBandas[iBanda]), 1] <- NA
      campo@data[iNoNA, 1] <- campo@data[iNoNA, 1] * escala + offset
      
      if (utmSur) {
        coords <- coordinates(campo)
        coords[,2] <- coords[,2] + 1E7
        centroides <- SpatialPoints(
          coords, proj4string=CRS(projargs=proj4stringResultados, SRS_string=SRS_stringResultados))
        campo <- SpatialGridDataFrame(
          grid=points2grid(centroides), data=campo@data, 
          proj4string=CRS(projargs=proj4stringResultados, SRS_string=SRS_stringResultados))
      }
      
      nombreArchivo <- appendToFileName(archivos[[i]], postFijo = '_tempCompressed') 
      #lock(lockFile, exclusive = TRUE)
      if (comprimirGeoTiffs) { options = c('COMPRESS=DEFLATE', 'PREDICTOR=2', 'ZLEVEL=9')
      } else { options = NULL }
      writeGDAL(dataset = campo, fname = nombreArchivo, options = options)
      aux <- readGDAL(nombreArchivo, silent = T)
      #unlock(lockFile)
  
      if (all((is.na(campo@data) & is.na(aux@data)) | (!is.na(campo@data) & !is.na(aux@data) & campo@data - aux@data < 1E-4))) {
        file.rename(from = nombreArchivo, to = archivos[[i]])
      } else {
        if (!interactive()) unlink(nombreArchivo)
        stop(paste('Fallo al guardar el archivo ', i, ' ', archivos[[i]], sep =''))
      }
    }
  }
  
  nCoresAUsar <- min(detectCores(T, T), length(archivos))
  if (nCoresAUsar > 1) {
    cl <- makeCluster(getOption('cl.cores', nCoresAUsar))
    clusterExport(cl, varlist = c('appendToFileName','getFileExt'))
    clusterEvalQ(cl, expr = {
      require('rgdal')
      #require('flock')
    })
    if (exists(x = 'setMKLthreads')) { clusterEvalQ(cl = cl, expr = setMKLthreads(1)) }
    parSapplyLB(cl=cl, X=1:length(archivos), FUN=postProcesarArchivoMODIS_V2, archivos=archivos, nombresBandas=nombresBandas, escala=escala, 
                offset= offset, minBandas=minBandas, maxBandas=maxBandas, utmSur=utmSur, proj4stringResultados=proj4stringResultados, 
                comprimirGeoTiffs=comprimirGeoTiffs)
    stopCluster(cl)
  } else {
    sapply(X=1:length(archivos), FUN=postProcesarArchivoMODIS_V2, archivos=archivos, nombresBandas=nombresBandas, escala=escala, 
           offset= offset, minBandas=minBandas, maxBandas=maxBandas, utmSur=utmSur, proj4stringResultados=proj4stringResultados, 
           comprimirGeoTiffs=comprimirGeoTiffs)
  }
  
  matchesFechas <- regexpr(archivos, pattern = '(19|20)[0-9][0-9][- /.](0[1-9]|1[012])[- /.](0[1-9]|[12][0-9]|3[01])')
  fechas <- substr(archivos, start = matchesFechas, stop = matchesFechas + attr(matchesFechas, 'match.length') - 1)
  fechas <- gsub(pattern = '-', replacement = '/', x = fechas, fixed = T)  
  
  setwd(oldWD)
    
  return(data.frame(fechas=fechas, archivos=archivos, stringsAsFactors = F))
}

getModisListEx <- function(producto, version, HsYVs, fechaIniUTC, fechaFinUTC, pathArchivosResultado, forceReDownload=F) {
  oldWD <- getwd()
  
  pathProductoMODIS <- subirDirectorio(pathArchivosResultado)
  setwd(pathProductoMODIS)
  
  modisList <- .getModisList(x=.modisHTTP(x=producto, v=version), h=HsYVs$hs, v=HsYVs$vs, dates=c(fechaIniUTC, fechaFinUTC), forceReDownload=forceReDownload)
  
  setwd(oldWD)
  
  return(modisList)
}

getFechasArchivos <- function(archivos) {
  matchesFechas <- regexpr(archivos, pattern = '(19|20)[0-9][0-9][- /.](0[1-9]|1[012])[- /.](0[1-9]|[12][0-9]|3[01])')
  fechas <- substr(archivos, start = matchesFechas, stop = matchesFechas + attr(matchesFechas, 'match.length') - 1)
  fechas <- gsub(pattern = '-', replacement = '/', x = fechas, fixed = T)  
  return(data.frame(fechas=fechas, archivos=archivos, stringsAsFactors = F))
}

getFechasYArchivosMODIS <- function(
    producto='MOD09Q1',
    # coordenadas límite en la proyección de salida del area a recortar
    xMin, xMax, yMin, yMax,
    # baja todos los datos entre estas dos fechas en la periodicidad del producto. El separador debe ser '.'
    fechaIniUTC='2014.09.06', fechaFinUTC='2014.09.06',
    # las bandas a obtener del dataset, se deben indicar con 0 las bandas a excluir y con 1 las bandas a incluir. 
    # Si se especifican menos bandas que las que hay en el dataset las restantes se asumen 0
    bandas='0 0 1',
    # versión del producto MODIS a descargar
    version='005',
    # ruta a la carpeta bin dentro del directorio de instalación de MRT
    MRTBinPath='C:/MRT/bin',
    # string en formato proj4 de la proyección de salida de los datos
    proj4stringResultados,
    SRS_stringResultados,
    pathArchivosResultado='./',
    escala=0.02, offset=-273.15, graficar=F, comprimirGeoTiffs=T) {
  objParams <- c(producto, xMin, xMax, yMin, yMax, fechaIniUTC, fechaFinUTC, bandas, version, proj4stringResultados, escala, offset, version=2)
  pathCache <- getPathCache(objParams, dirEjecucion=pathArchivosResultado)
  if (!file.exists(pathCache)) {
    fechasArchivos <- getMODISEx(
      producto=producto, xMin=xMin, xMax=xMax, yMin=yMin, yMax=yMax, fechaIniUTC=fechaIniUTC, 
      fechaFinUTC=fechaFinUTC, bandas=bandas, version=version, MRTBinPath=MRTBinPath, 
      proj4stringResultados=proj4stringResultados, SRS_stringResultados=SRS_stringResultados,
      pathArchivosResultado=pathArchivosResultado, graficar=graficar, 
      comprimirGeoTiffs=comprimirGeoTiffs)
    
    temporalRes <- rts:::getNativeTemporalResolution(producto)
    fechaAux <- addDate_temporalResolution(as.Date(fechasArchivos$fechas[nrow(fechasArchivos)], format='%Y/%m/%d'), temporalRes)
    fechaAux <- format(fechaAux, format='%Y.%m.%d')
    # Si la condicion se cumple quiere decir que pude bajar el período al que corresponde fechaFinUTC
    # por lo tanto puedo cachear los resultados
    #if (fechaFinUTC < fechaAux & nrow(fechasArchivos) == length(modisList)) guardarCache(pathCache, fechasArchivos)
    if (fechaFinUTC < fechaAux) guardarCache(pathCache, fechasArchivos)
  } else {
    fechasArchivos <- cargarCache(pathCache)
  }

  # Verifico que los archivos hayan descargado correctamente (si los puedo leer con readGDAL)
  # Si no, espero 15 segundos, vuelvo a descargarlos y reproyectarlos hasta 10 veces
  i <- 1
  for (i in seq_along(fechasArchivos$archivos)) {
    hecho <- FALSE
    nIntentos <- 0
    while (!hecho & nIntentos < 10) {
      result = tryCatch({
        readGDAL(fechasArchivos$archivos[[i]], silent=T)
      }, error = function(e) {
        print(paste('MY_ERROR: ', e, fechasArchivos$archivos[[i]]))
        result <- list()
        class(result) <- 'errorRetry'
        return(result)
      })
      
      if (class(result) == 'errorRetry') {
        nIntentos <- nIntentos + 1
        Sys.sleep(15)
        
        fecha <- gsub(pattern = '/', replacement = '.', x = fechasArchivos$fechas[i], fixed = T)
        
        getMODISEx(
          producto=producto, xMin=xMin, xMax=xMax, yMin=yMin, yMax=yMax, fechaIniUTC=fecha, 
          fechaFinUTC=fecha, bandas=bandas, MRTBinPath=MRTBinPath, 
          proj4stringResultados=proj4stringResultados, SRS_stringResultados=SRS_stringResultados,
          pathArchivosResultado=pathArchivosResultado, graficar=graficar, 
          forceReDownload=T, comprimirGeoTiffs=comprimirGeoTiffs)
      } else {
        hecho <- TRUE
      }
    }
  }
  
  return(fechasArchivos)
}
