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
if (iFrame >= 3) { script.dir.descargadorEx <- sys.frame(iFrame - 3)$ofile 
} else { script.dir.interpolarEx <- NULL }
while ((is.null(script.dir.descargadorEx) || is.na(regexpr('descargadorEx.r', script.dir.descargadorEx, fixed=T)[1])) && iFrame >= 0) {
  script.dir.descargadorEx <- sys.frame(iFrame)$ofile
  iFrame <- iFrame - 1
}
if (is.null(script.dir.descargadorEx)) { script.dir.descargadorEx <- ''
} else { script.dir.descargadorEx <- paste(dirname(script.dir.descargadorEx), '/', sep='') }

source(paste(script.dir.descargadorEx, '../instalarPaquetes/instant_pkgs.r', sep=''))
source(paste(script.dir.descargadorEx, '../PathUtils/pathUtils.r', sep=''))
instant_pkgs(c('RCurl', 'parallel', 'digest', 'data.table', 'lubridate'))

threadHandle <- getCurlHandle()

getURLEx <- function(url) { return(getURL(url, curl = curlHandle, async = FALSE)) }

extraerEnlaces <- function(urls, patronEnlaces='href="[[:print:]]"', concatenarUrlBase=FALSE,
                           maxNConexiones=4) {
  maxNConexiones <- min(maxNConexiones, length(urls))
  if (maxNConexiones > 1) {
    cl <- makeCluster(getOption("cl.cores", maxNConexiones))
    clusterEvalQ(cl, expr = { 
      require('RCurl')
      curlHandle <- getCurlHandle()
    })
    clusterExport(cl = cl, varlist = c('getURLEx'))
    textoshtmls <- parSapplyLB(cl=cl, X=urls, FUN = getURLEx)
    stopCluster(cl)    
  } else {
    curlHandle <- getCurlHandle()
    textoshtmls <- getURL(url = urls, async = FALSE, curl = curlHandle)
    rm(curlHandle)
  }
  textoshtmls <- strsplit(x = textoshtmls, split = '\r*\n')
  textohtml <- unlist(textoshtmls)

  res <- regexpr(pattern = patronEnlaces, text = textohtml)
  urls <- regmatches(textohtml, res)
  if (concatenarUrlBase) {
    urlsBase <- unlist(lapply(seq_along(textoshtmls), FUN = function(x) { rep(names(textoshtmls)[x], length(textoshtmls[[x]]))}))
    urlsBase <- urlsBase[res > 0]
    urls <- paste(urlsBase, urls, sep='')
  }
  return(urls)
}

generarNombresArchivo_CHIRP <- function(fechaMin, fechaMax=fechaMin,
                                        urlBaseCHIRP='ftp://ftp.chg.ucsb.edu/pub/org/chg/products/CHIRP/daily/',
                                        carpetaSalida='D:/testsMCH/descargador/CHIRP/') {
  #fechaMin <- parse_date_time('1981-01-01', orders = 'YmdHMS', tz='UTC', truncated = 5)
  #fechaMax <- Sys.time()
  
  #fechaMin <- parse_date_time('2018-03-01', orders = 'YmdHMS', tz='UTC', truncated = 5)
  #fechaMax <- Sys.time()
  #carpetaSalida='I:/INIA/Datos/CHIRP/'
  
  anioMin <- year(fechaMin)
  anioMax <- year(fechaMax)
  anios <- anioMin:anioMax
  
  urls <- paste(urlBaseCHIRP, anios, '/', sep = '')
  
  urls <- extraerEnlaces(urls = urls, 
                         patronEnlaces = 'chirp\\.[[:digit:]]{4}\\.[[:digit:]]{2}\\.[[:digit:]]{2}\\.tif[\\.gz]?',
                         concatenarUrlBase = TRUE)
  
  nombresArchivosDestino <- paste(carpetaSalida, basename(urls), sep='')

  pathArchivo <- nombresArchivosDestino[1]
  procesarArchivoTifGz <- function(pathArchivo, path7Zip='\"C:/Program Files/7-Zip/7z.exe\"',
                                   pathSHPRecorte=NULL) {
     #D:\testsMCH\descargador\CHIRP\1981\chirp.1981.01.01.tif.gz
    oldwd <- getwd()
    setwd(dirname(pathArchivo))

    system(paste(path7Zip, ' e -aoa ', pathArchivo, sep=''))
    pathArchivoDescomprimido <- substr(pathArchivo, 1, nchar(pathArchivo) - 3)
    
    require(raster)
    
    capaRaster <- raster(pathArchivoDescomprimido)
    #plot(capaRaster)
    
    if (!is.null(pathSHPRecorte)) {
      shpBase <- cargarSHP('C:/mch/ArchivosProcesosLocales/CartografiaBase/uruguay_departamentos.shp')
      shpBase <- spTransform(shpBase, CRS(proj4string(capaRaster)))
      areaRecorte <- getPoligonoBoundingBox(objSP = shpBase, factorExtensionX = 1.2)
      capaRaster <- crop(capaRaster, extent(areaRecorte))
      
      #plot(capaRasterRecortada)
    }
    
    writeGDAL(dataset = as(object = capaRasterRecortada, Class = 'SpatialGridDataFrame'), 
              fname = appendToFileName(pathArchivoDescomprimido, '_recortadoUy'), 
              options = c('COMPRESS=DEFLATE', 'PREDICTOR=2', 'ZLEVEL=9'))
    
    unlink(pathArchivo)
    unlink(pathArchivoDescomprimido)
    
    setwd(oldwd)
  }
    
  descargarArchivos(urls = urls, nombresArchivosDestino = nombresArchivosDestino)
  
  return(list(urls=urls, nombresArchivosDestino=nombresArchivosDestino))
}


generarNombresArchivo_RAMMB <- function(fechaEjecucion=Sys.time(), 
                                        urlBaseRAMMB='http://rammb.cira.colostate.edu/ramsdis/online/',
                                        carpetaSalida='D:/testsMCH/descargador/RAMMB/') {
  if (minute(fechaEjecucion) %% 2 == 0) {
    # Para pruebas
    #textohtml <- unlist(strsplit(x = getURLAsynchronous(url = 'http://rammb.cira.colostate.edu/ramsdis/online/archive.asp?data_folder=rmtc/rmtcsasec1vis04&width=640&height=480'), 
    #                                                    split = '\r*\n'))
    textohtml <- unlist(strsplit(x = getURLAsynchronous(url = paste(urlBaseRAMMB, 
                                                                    c('archive.asp?data_folder=rmtc/rmtcsasec1vis04&width=640&height=480',
                                                                      'archive.asp?data_folder=rmtc/rmtcsasec1ir304&width=640&height=480',
                                                                      'archive.asp?data_folder=rmtc/rmtcsasec1ir404&width=640&height=480',
                                                                      'archive.asp?data_folder=rmtc/rmtcsasec1ir204&width=640&height=480'), sep = '')),
                                 split = '\r*\n'))
    
    urls <- regmatches(textohtml, regexpr('images/rmtc/[[:print:]]+.gif', textohtml))
    if (length(urls) > 0) {
      urls <- paste(urlBaseRAMMB, urls, sep='')
      nombresArchivosDestino <- paste(carpetaSalida, 'RAMMB/', pathDesdeUltimaCarpeta(urls), sep='')
    } else { nombresArchivosDestino <- character(0) }
    return(list(urls=urls, nombresArchivosDestino=nombresArchivosDestino))
  } else return(NULL)
}

generarNombresArchivo_SMN_Satelite <- function(fechaEjecucion=Sys.time(), 
                                               urlBaseSMN='http://www.smn.gov.ar/vmsr/',
                                               carpetaSalida='D:/testsMCH/descargador/') {
  if (minute(fechaEjecucion) %% 2 == 0) {
    textohtml <- unlist(strsplit(x = getURLAsynchronous(url = paste(urlBaseSMN, 
                                                                    c('general.php?dir=YVcxaFoyVnVaWE12WVhKblpXNTBhVzVoTDJsdVpnPT0=',
                                                                      'general.php?dir=YVcxaFoyVnVaWE12WVhKblpXNTBhVzVoTDNSdQ==',
                                                                      'general.php?dir=YVcxaFoyVnVaWE12WVhKblpXNTBhVzVoTDNaaA==',
                                                                      'general.php?dir=YVcxaFoyVnVaWE12WVhKblpXNTBhVzVoTDNacGN3PT0='), sep = '')),
                                 split = '\r*\n'))
    
    prefix = '<script language="javascript" type="text/javascript">'
    textohtml <- textohtml[substring(textohtml, 1, nchar(prefix)) == prefix]
    
    #textohtml <- textohtml[startsWith(textohtml, prefix = '<script language="javascript" type="text/javascript">')]
    
    varDirs <- paste(gsub(x=gsub(x = regmatches(textohtml, regexpr('var dir=\"[[:print:]]+\";</script> ', textohtml)), pattern = 'var dir=\"', replacement = '', fixed = T),
                     pattern = '\";</script> ', replacement = '', fixed = T), '/', sep='')
    
    textohtml <- gsub(x=gsub(x = textohtml, pattern = '<script language=\"javascript\" type=\"text/javascript\"> var img=\"', replacement = '', fixed = T),
                      pattern = '|\"; var dir=\"imagenes/argentina/[[:alpha:]]{2,3}\";</script> ', replacement = '')
    urls <- strsplit(x = textohtml, split = '|', fixed=T)
    if (length(urls) > 0) {
      urls <- unlist(lapply(seq_along(urls), FUN = function(x, urls, varDirs) {return(paste(varDirs[x], urls[[x]], sep=''))}, urls=urls, varDirs=varDirs))
      urls <- paste(urlBaseSMN, urls, sep='')
      nombresArchivosDestino <- paste(carpetaSalida, 'SMN/satelite/', pathDesdeUltimaCarpeta(urls), sep='')
    } else { nombresArchivosDestino <- character(0) }
    
    return(list(urls=urls, nombresArchivosDestino=nombresArchivosDestino))
  } else return(NULL)
}

generarNombresArchivo_SMN_Radar <- function(fechaEjecucion=Sys.time(), 
                                            urlBaseSMNRadar='http://www.smn.gov.ar/radar/',
                                            carpetaSalida='D:/testsMCH/descargador/') {
  if (minute(fechaEjecucion) %% 5 == 0) {
    urls <- c(paste(urlBaseSMNRadar, 'CMAX_240_ZE_1_01.png', sep=''),
              paste(urlBaseSMNRadar, 'ETOP_240_ZE_1_01.png', sep=''),
              paste(urlBaseSMNRadar, 'CMAX_480_ZE_1_01.png', sep=''),
              paste(urlBaseSMNRadar, 'PPI1_120_VE_1_01.png', sep=''))
    
    strFechaDescarga <- format(x = Sys.time(), "_%Y%m%d_%H%M%S")
    nombresArchivosDestino <- paste(carpetaSalida, 'SMN/', appendToFileName(pathDesdeUltimaCarpeta(urls), strFechaDescarga), sep='')
    return(list(urls=urls, nombresArchivosDestino=nombresArchivosDestino))
  } else return(NULL)
}

generarNombresArchivo_DatosMCH2_Horarios <- function(fechaEjecucion=Sys.time(), 
                                            urlBaseDatosMCH2='https://www.inumet.gub.uy/reportes/estadoActual/',
                                            carpetaSalida='D:/testsMCH/descargador/') {
  if (minute(fechaEjecucion) %% 5 == 0) {
    urls <- c(paste(urlBaseDatosMCH2, 'estadoActualDatosHorarios.mch', sep=''))
    
    strFechaDescarga <- paste(format(x = fechaEjecucion, "%Y%m%d_%H"), '0000_', sep='')
    nombresArchivosDestino <- paste(carpetaSalida, 'MCH2/', prependToFileName(filename = pathDesdeUltimaCarpeta(urls), strFechaDescarga), sep='')
    return(list(urls=urls, nombresArchivosDestino=nombresArchivosDestino))
  } else return(NULL)
}

generarNombresArchivo_DatosMCH2_12Horarios <- function(fechaEjecucion=Sys.time(), 
                                                       urlBaseDatosMCH2='https://www.inumet.gub.uy/reportes/estadoActual/',
                                                       carpetaSalida='D:/testsMCH/descargador/') {
  if (hour(fechaEjecucion) %in% c(0, 12) && minute(fechaEjecucion) %% 5 == 0) {
    urls <- c(paste(urlBaseDatosMCH2, 'estadoActualDatosDiariosUltimaSemana.mch', sep=''))
    
    strFechaDescarga <- paste(format(x = fechaEjecucion, "%Y%m%d_%H"), '0000_', sep='')
    nombresArchivosDestino <- paste(carpetaSalida, 'MCH2/', prependToFileName(filename = pathDesdeUltimaCarpeta(urls), strFechaDescarga), sep='')
    return(list(urls=urls, nombresArchivosDestino=nombresArchivosDestino))
  } else return(NULL)
}
  
generarNombresArchivo_DatosMCH2_Diarios <- function(fechaEjecucion=Sys.time(), 
                                                    urlBaseDatosMCH2='https://www.inumet.gub.uy/reportes/estadoActual/',
                                                    carpetaSalida='D:/testsMCH/descargador/') {
  # Por ahora solo se suben diarios los datos de R3, cuando se suban más datos ajustar las fechas de descarga
  if (hour(fechaEjecucion) %in% c(13, 14) && minute(fechaEjecucion) > 30 && minute(fechaEjecucion) %% 5 == 0) {
    urls <- c(paste(urlBaseDatosMCH2, 'estadoActualDatosDiariosUltimaSemanaR3.mch', sep=''))
    
    strFechaDescarga <- paste(format(x = fechaEjecucion, "%Y%m%d"), '_000000_', sep='')
    nombresArchivosDestino <- paste(carpetaSalida, 'MCH2/', prependToFileName(filename = pathDesdeUltimaCarpeta(urls), strFechaDescarga), sep='')
    return(list(urls=urls, nombresArchivosDestino=nombresArchivosDestino))
  } else return(NULL)
}

generarNombresArchivo_WWLLN <- function(fechaEjecucion=Sys.time(), 
                                        urlBaseWWLLN='http://webflash.ess.washington.edu/',
                                        carpetaSalida='D:/testsMCH/descargador/') {
  if (minute(fechaEjecucion) %% 5 == 0) {
    urls <- paste(urlBaseWWLLN, 'AmericaL_plot_weather_map.jpg', sep='')
    nombresArchivosDestino <- paste(carpetaSalida, 'WWLLN/', appendToFileName(filename = basename(urls), postFijo = format(x = Sys.time(), "_%Y%m%d_%H%M%S")), sep='')
    return(list(urls=urls, nombresArchivosDestino=nombresArchivosDestino))
  } else return(NULL)
}

generarNombresArchivo_PrecipGFS <- function(fechaEjecucion=Sys.time(), 
                                            urlBaseGFS='http://mag.ncep.noaa.gov/data/gfs/',
                                            zonaGFS='samer',
                                            carpetaSalida='D:/testsMCH/descargador/',
                                            horasPrecipADescargar=c(1, 3, 6, 24)) {
  fechaEjecucionUTC <- with_tz(fechaEjecucion, tzone = 'UTC')
  horaEjec <- hour(fechaEjecucionUTC)
  if (horaEjec %% 12 == 0 && minute(fechaEjecucionUTC) == 30) {
    horasPrecip <- horasPrecipADescargar[3]
    
    urls <- character(0)
    nombresArchivosDestino <- character(0)
    for (horasPrecip in horasPrecipADescargar) {
      strHorasPrecip <- sprintf('%02d', horasPrecip)
      if (horasPrecip == 1) { horasPron <- sprintf('%03d', seq.int(from=1, to=72, by=1))
      } else { horasPron <- sprintf('%03d', seq.int(from=horasPrecip, to=72, by=3)) }
      
      nuevosUrls <- paste(urlBaseGFS, horaEjec, '/', zonaGFS, '/precip_p', strHorasPrecip, '/gfs_', zonaGFS, '_', horasPron, '_precip_p', strHorasPrecip, '.gif', sep='')
      urls <- c(urls, nuevosUrls)
      nombresArchivosDestino <- c(nombresArchivosDestino, paste(carpetaSalida, 'GFS/precip/', format(fechaEjecucionUTC, '%Y%m%d_%H'), '/precip_p', strHorasPrecip, '/', basename(nuevosUrls), sep=''))
    }
    return(list(urls=urls, nombresArchivosDestino=nombresArchivosDestino))
  } else return(NULL)
}

postprocesarArchivosSinFecha <- function(nombresArchivosDestino) {
  carpetas <- unique(dirname(nombresArchivosDestino))
  i<-1
  for (i in seq_along(carpetas)) {
    archivosEnCarpeta <- setdiff(dir(carpetas[[i]], full.names = F, recursive = F, include.dirs = F, 
                                     no.. = T), 'hashes.csv')
    
    archiHashes <- paste(carpetas[[i]], '/hashes.csv', sep='')
    if (file.exists(archiHashes)) {
      hashes <- fread(archiHashes, sep = ',', header = T, stringsAsFactors = F, key = 'hash', data.table = T)
      iNuevos <- which(!archivosEnCarpeta %in% hashes$filename)
      #archivosEnCarpeta[iNuevos]
      if (length(iNuevos) > 0) {
        hashesNuevos <- sapply(X = paste(carpetas[[i]], '/', archivosEnCarpeta[iNuevos], sep=''), FUN = digest, algo = 'md5', file = TRUE)
        hashesNoRepetidos <- !hashesNuevos %in% hashes$hash
        iHashesRepetidos <- which(!hashesNoRepetidos)
        
        if (length(iHashesRepetidos) > 0) {
          file.remove(paste(carpetas[[i]], '/', archivosEnCarpeta[iNuevos[iHashesRepetidos]], sep=''))
          iNuevos <- iNuevos[hashesNoRepetidos]
          hashesNuevos <- hashesNuevos[hashesNoRepetidos]
        }
        
        if (length(hashesNuevos) > 0) {
          # Agrego los archivos nuevos al hash  
          hashes <- rbindlist(list(hashes, 
                                   list(filename=archivosEnCarpeta[iNuevos], 
                                        hash=hashesNuevos)))
          hayNuevos <- T
        } else {
          hayNuevos <- F
        }
      } else {
        hayNuevos <- F
      }
    } else {
      hashes <- data.table(filename=archivosEnCarpeta, 
                           hash=sapply(X = paste(carpetas[[i]], '/', archivosEnCarpeta, sep=''), FUN = digest, algo = 'md5', file = TRUE),
                           key='hash')
      
      iDuplicados <- duplicated(hashes, by = 'hash')
      if (any(iDuplicados)) {
        file.remove(paste(carpetas[[i]], '/', archivosEnCarpeta[iDuplicados], sep=''))
        hashes <- hashes[!iDuplicados]
      }
      hayNuevos <- T
    }
    if (hayNuevos) {
      if (is.unsorted(hashes$filename)) hashes <- hashes[order(filename)]
      fwrite(x = hashes, file = archiHashes)
    }
  }
}

isCompressed <- function(paths) {
  return( sapply(paths, FUN = function(x) {return(R.utils::isGzipped(x) | isBzipped(x))}) )
}

descargarArchivo <- function(
    i, urls, nombresArchivosDestino, forzarReDescarga=FALSE, maxRetries=5, segundosEntreIntentos=15, 
    curlOpts=NULL, do_unzip=isCompressed(nombresArchivosDestino)) {
  # i <- 1
  nRetries <- 0

  # Check if unzipped file exists
  if (do_unzip[i]) {
    unzippedFilename <- nombreArchSinExtension(nombresArchivosDestino[i])
    if (forzarReDescarga) {
      downloadExists <- FALSE
      unzipExists <- FALSE
    } else {
      downloadExists <- file.exists(nombresArchivosDestino[i]) && 
                        length(readBin(nombresArchivosDestino[i], what='raw')) > 0
      unzipExists <- file.exists(unzippedFilename) && 
                     length(readBin(unzippedFilename, what='raw')) > 0
    }
  } else {
    downloadExists <- !forzarReDescarga && file.exists(nombresArchivosDestino[i]) && 
                      length(readBin(nombresArchivosDestino[i], what='raw')) > 0
    unzipExists <- FALSE
  }
  
  do_download <- forzarReDescarga || (!downloadExists && !unzipExists)
  if (do_download) {
    is_ftp = startsWith(urls[i], 'ftp://') || startsWith(urls[i], 'ftps://')
    is_http = !is_ftp && (startsWith(urls[i], 'http://') || startsWith(urls[i], 'https://'))
    results <- 1
    
    while (do_download & nRetries < maxRetries) {
      f = CFILE(nombresArchivosDestino[i], mode="wb")
      er2 <- try(er <- curlPerform(url = urls[i], curl=threadHandle, writedata = f@ref))
      close(f)
      if (class(er2) == "try-error" || 
          er != 0 || 
          !file.exists(nombresArchivosDestino[i]) || length(readBin(nombresArchivosDestino[i], what='raw')) <= 0) {
        
        if (is_ftp && grepl(pattern = '5[[:digit:]]{2}', er2)) {
          # Handling FTP permanent error cases. Needs improvement
          nRetries <- maxRetries
        } else if (is_http && grepl(pattern = '4[[:digit:]]{2}', er2)) {
          # Handling HTTP client error cases. Needs improvement
          nRetries <- maxRetries
        } else {
          nRetries <- nRetries + 1
          Sys.sleep(segundosEntreIntentos)
        }
        
        # Reinicio el handle si hay algún problema
        threadHandle <- getCurlHandle(.opts = curlOpts)
      } else { do_download <- FALSE }
    }
  } else {
    results <- 2
  }
  do_unzip_i <- !do_download && do_unzip[i] && (forzarReDescarga || !unzipExists)
  if (do_unzip_i) {
    ext <-  getFileExt(nombresArchivosDestino[i])
    if (ext == "gz") { decompFunc <- gzfile 
    } else if (ext == "bz2") { decompFunc <- bzfile 
    } else { decompFunc <- NULL }
    R.utils::decompressFile(filename=nombresArchivosDestino[i], ext=ext, overwrite=TRUE, 
                            FUN=decompFunc)
  }
  
  if (do_download) { 
    results <- 0
    unlink(x = nombresArchivosDestino[i])
  }
  return(results)
}

descargarArchivo_SinKeepalive <- function(i, urls, nombresArchivosDestino, maxRetries=5, segundosEntreIntentos=15) {
  nRetries <- 0
  success <- FALSE
  while (!success & nRetries < maxRetries) {
    f = CFILE(nombresArchivosDestino[i], mode="wb")
    er2 <- try(er <- curlPerform(url = urls[i], writedata = f@ref))
    close(f)
    
    if (class(er2) == "try-error" || (class(er) == "integer" && er != 0)) {
      nRetries <- nRetries + 1
      Sys.sleep(segundosEntreIntentos)
    } else { success <- TRUE }
  }  
  return(success)
}

descargarArchivo_curl <- function(i, urls, nombresArchivosDestino, maxRetries=5, segundosEntreIntentos=15) {
  nRetries <- 0
  success <- FALSE
  while (!success & nRetries < maxRetries) {
    req <- curl_fetch_disk(url = urls[i], path = nombresArchivosDestino[i])
    if (req$status_code >= 400) {
      nRetries <- nRetries + 1
      Sys.sleep(segundosEntreIntentos)
    } else { success <- TRUE }
  }  
  return(success)
}

descargarArchivos <- function(urls, nombresArchivosDestino=paste(pathSalida, basename(urls), sep=''), 
                              nConexionesSimultaneas=4, forzarReDescarga=FALSE, curlOpts=NULL, 
                              pathSalida='', do_unzip=isCompressed(nombresArchivosDestino)) {
  # Retorna:
  # 0 si no se pudo bajar el archivo
  # 1 si se bajo un archivo nuevo
  # 2 si el archivo ya existía
  #nConexionesSimultaneas=4
  #forzarReDescarga=FALSE
  #archivos <- generarNombresArchivoSMN()
  #urls <- archivos$urls
  #nombresArchivosDestino <- archivos$nombresArchivosDestino
  
  if (length(urls) > 0) {
    sapply(unique(dirname(nombresArchivosDestino)), dir.create, showWarnings = F, recursive=T)
    
    # if (length(urls) * 4 < nConexionesSimultaneas) nConexionesSimultaneas <- max(trunc(length(urls) / 4), 1)
    nConexionesAUsar <- min(nConexionesSimultaneas, length(urls))
    
    # nConexionesAUsar descargas paralelas con Keepalive usando RCurl
    if (nConexionesAUsar > 1) {
      cl <- makeCluster(getOption('cl.cores', nConexionesAUsar))
      clusterExport(cl, varlist = c('curlOpts', 'script.dir.descargadorEx'), envir = environment())
      clusterEvalQ(cl, expr = {
        source(paste(script.dir.descargadorEx, '../PathUtils/pathUtils.r', sep=''))
        require('RCurl')
        threadHandle <- getCurlHandle(.opts = curlOpts)
      })
      results <- parSapplyLB(cl = cl, X=seq_along(urls), FUN = descargarArchivo, urls=urls, 
                         nombresArchivosDestino=nombresArchivosDestino, 
                         forzarReDescarga=forzarReDescarga, curlOpts = curlOpts, do_unzip=do_unzip)
      stopCluster(cl)
    } else {
      assign("threadHandle", getCurlHandle(.opts = curlOpts), envir = .GlobalEnv)
      # threadHandle <- getCurlHandle(.opts = curlOpts)
      results <- sapply(X=seq_along(urls), FUN = descargarArchivo, urls=urls,
                        nombresArchivosDestino=nombresArchivosDestino, 
                        forzarReDescarga=forzarReDescarga, curlOpts = curlOpts, do_unzip=do_unzip)
    }
    # descargas secuenciales con Keepalive usando RCurl
    #cl <- makeCluster(getOption('cl.cores', 1))
    #system.time({
    #  clusterEvalQ(cl, expr = { 
    #    require('RCurl')
    #    threadHandle <- getCurlHandle() })
    #  parLapplyLB(cl = cl, X=seq_along(urls), fun = descargarArchivo, urls=urls, nombresArchivosDestino=nombresArchivosDestino, 
    #              nConexionesSimultaneas=nConexionesAUsar)
    #})
    #stopCluster(cl)
    
    # nConexionesAUsar descargas paralelas con Keepalive usando curl
    #cl <- makeCluster(getOption('cl.cores', nConexionesAUsar))
    #system.time({
    #  clusterEvalQ(cl, expr = { require('curl') })
    #  parLapplyLB(cl = cl, X=seq_along(urls), fun = descargarArchivo_curl, urls=urls, nombresArchivosDestino=nombresArchivosDestino)
    #})
    #stopCluster(cl)
  } else {
    results <- integer(0)
  }
  return(results)
}
