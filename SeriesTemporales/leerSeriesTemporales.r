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
if (iFrame >= 3) { script.dir.lecturaDatos <- sys.frame(iFrame - 3)$ofile 
} else { script.dir.lecturaDatos <- NULL }
while ((is.null(script.dir.lecturaDatos) || is.na(regexpr('lecturaDatos.r', script.dir.lecturaDatos, fixed=T)[1])) && iFrame >= 0) {
  script.dir.lecturaDatos <- sys.frame(iFrame)$ofile
  iFrame <- iFrame - 1
}
if (is.null(script.dir.lecturaDatos)) { script.dir.lecturaDatos <- ''
} else { script.dir.lecturaDatos <- paste(dirname(script.dir.lecturaDatos), '/', sep='') }

source(paste(script.dir.lecturaDatos, '../instalarPaquetes/instant_pkgs.r', sep=''))
instant_pkgs(pkgs = c('Rcpp', 'stringi', 'lubridate', 'jsonlite', 'openxlsx'))

setIdsEstaciones <- function(dfEstaciones, colId=1) {
  dfEstaciones[, colId] <- make.names(trimws(dfEstaciones[, colId]))
  dfEstaciones[, colId] <- gsub(pattern='\\.{2,}', replacement = '\\.', x=dfEstaciones[, colId])
  dfEstaciones[, colId] <- gsub(pattern='\\.$', replacement = '', x=dfEstaciones[, colId])
  rownames(dfEstaciones) <- dfEstaciones[, colId]
  return(dfEstaciones)
}

limpiarDatos <- function(dfDatos, dfEstaciones, colIdEstaciones=-1L, header=T, overrideHeader=F,
                         formatoFechas='YmdHMS', truncated=5, tzFechas='UTC') {
  fechas <- parse_date_time(as.character(dfDatos[, 1]), orders = formatoFechas, tz=tzFechas, 
                            truncated = truncated)
  iFechasValidas <- !is.na(fechas)
  fechas <- fechas[iFechasValidas]
  dfDatos <- data.matrix(dfDatos[iFechasValidas, -1])
  
  if (is.numeric(colIdEstaciones) && colIdEstaciones <= 0) { idsEstaciones <- row.names(dfEstaciones)
  } else { idsEstaciones <- make.names(dfEstaciones[,colIdEstaciones]) }
  if (overrideHeader) colnames(dfDatos) <- idsEstaciones 
  
  if (header && !is.null(dfEstaciones)) {
    ies <- na.omit(match(idsEstaciones, colnames(dfDatos)))
    datos <- matrix(data = NA_real_, nrow = nrow(dfDatos), ncol = nrow(dfEstaciones))
    datos[, ies] <- dfDatos[,ies]
  } else datos <- dfDatos
  
  rownames(datos) <- as.character(fechas)
  if (!is.null(dfEstaciones)) colnames(datos) <- idsEstaciones
  
  return (list(fechas=fechas, datos=datos))
}

leerEstaciones <- function(pathArchivoEstaciones, columnaId=1, fileEncoding = '', sep=',') {
  estaciones <- read.table(pathArchivoEstaciones, header=T, sep=sep, dec='.', stringsAsFactors=F, fileEncoding = fileEncoding)
  # Saco los espacios a los nombres de estaciones y agrego una X a las que empiezan con números porque dificultan el trabajo con los dataframes
  if (columnaId > 0) {
    estaciones[,columnaId] <- limpiarIdsEstaciones(estaciones[,columnaId])
    rownames(estaciones) <- estaciones[,columnaId]
  }
  
  return (estaciones)
}

leerDatos <- function(pathArchivoDatos, dfEstaciones, skip=0L, formatoFechas='YmdHMS', truncated=5, tzFechas='UTC', header=T, sep=',', dec='.', na.strings='-9999', fileEncoding = '', colIdEstaciones=-1L) {
  dfDatos <- read.table(pathArchivoDatos, header=header, sep=sep, skip=skip, dec=dec, stringsAsFactors=F, na.strings=na.strings, fileEncoding = fileEncoding)
  return (limpiarDatos(dfDatos = dfDatos, dfEstaciones = dfEstaciones, 
                       colIdEstaciones = colIdEstaciones, header = header, 
                       formatoFechas = formatoFechas, truncated = truncated, tzFechas = tzFechas))
}

grabarDatos <- function(pathArchivoDatos, fechas, datos, sep='\t', dec='.', na='-9999', append=F, col.names=T, formatoFechas='%Y-%m-%d') {
  write.table(x=data.frame(Fechas=format(fechas, format = formatoFechas), datos), file=pathArchivoDatos, append=append, sep=sep, dec=dec, na=na, row.names=F, col.names=col.names)
}

leerSeriesArchivoUnico <- function(pathArchivoDatos, nFilasEstaciones=10, filaId=1, skip=0, formatoFechas='YmdHMS', truncated=5, tzFechas='UTC', headerDatos=F, sep='\t', dec='.', na.strings=c('-9999','NA'), fileEncoding = '') {
  estaciones <- read.table(pathArchivoDatos, header=F, sep=sep, dec=dec, na.strings=na.strings, nrows=nFilasEstaciones, skip=skip, stringsAsFactors=F, fileEncoding = fileEncoding, row.names = 1)
  estaciones <- as.data.frame(t(estaciones), stringsAsFactors=F, sep=sep)
  for (i in 1:ncol(estaciones)) {
    trimEstacionesI <- trimws(estaciones[,i])
    upperTrimEstacionesI <- toupper(trimEstacionesI)
    
    aux <- type.convert(upperTrimEstacionesI, na.strings=na.strings, as.is=T, dec=dec)
    if (is.logical(aux)) { estaciones[,i] <- aux
    } else { estaciones[,i] <- type.convert(trimEstacionesI, na.strings=na.strings, as.is=T, dec=dec) }
  }
  
  estaciones <- setIdsEstaciones(estaciones, filaId)

  datos <- leerDatos(pathArchivoDatos, dfEstaciones=estaciones, skip=skip+nFilasEstaciones, formatoFechas=formatoFechas, truncated = truncated, tzFechas=tzFechas, header=headerDatos, sep=sep, dec=dec, na.strings=na.strings, fileEncoding = fileEncoding)
  colnames(datos$datos) <- rownames(estaciones)
  return(list(estaciones=estaciones, fechas=datos$fechas, datos=datos$datos))
}

leerSeriesXLSX <- function(pathArchivoDatos, hojaEstaciones='InfoPluvios', headerEstaciones=T,
                           colsEstaciones=1:4, colId=2, hojaDatos='Medidas', formatoFechas='YmdHMS',
                           truncated=5, tzFechas='UTC', headerDatos=T, fileEncoding = '',
                           na.strings=-1111) {
  
  
  dfEstaciones <- read.xlsx(xlsxFile = pathArchivoDatos, sheet = hojaEstaciones, 
                            cols = colsEstaciones, colNames = headerEstaciones)
  dfEstaciones <- setIdsEstaciones(dfEstaciones, colId = colId)

  dfDatos <- read.xlsx(xlsxFile = pathArchivoDatos, sheet = hojaDatos, colNames = headerDatos)
  dfDatos <- limpiarDatos(dfDatos = dfDatos, dfEstaciones = dfEstaciones, colIdEstaciones = colId, 
                          header = headerDatos, overrideHeader = headerDatos, 
                          formatoFechas = formatoFechas, truncated = truncated, tzFechas = tzFechas)
  dfDatos$datos[apply(dfDatos$datos, MARGIN = 2, FUN = function(x) {x %in% na.strings})] <- NA
  return(list(estaciones=dfEstaciones, fechas=dfDatos$fechas, datos=dfDatos$datos))
}

grabarSeriesArchivoUnico <- function(pathArchivoDatos, estaciones, fechas, datos, sep='\t', dec='.', na='-9999',
                                     colsEstacionesAEscribir=1:ncol(estaciones), formatoFechas='%Y-%m-%d') {
  unlink(pathArchivoDatos)
  
  # i <- colsEstacionesAEscribir[1]
  for (i in colsEstacionesAEscribir) {
    x <- t(estaciones[, i])
    rownames(x) <- colnames(estaciones)[i]
    write.table(x=x, file=pathArchivoDatos, sep=sep, dec=dec,na=na, row.names=T, col.names=F, append = (i != colsEstacionesAEscribir[1]))
  }
  
  rownames(datos) <- format(fechas, format=formatoFechas)
  write.table(x=datos, file=pathArchivoDatos, append=T, sep=sep, dec=dec, na=na, row.names=T, col.names=F)
}

recolectarSalidasCDT <- function(carpetaSalidaCDT='F:/Tesis/CDT/QcTemp_DatosEstacionesCorreccionCoords_ParaCDT/CorrectedData',
                                 estaciones=NULL, colIdEstaciones=1) {
  pathsEstaciones <- dir(carpetaSalidaCDT, include.dirs = T, recursive = F)
  if (is.null(estaciones)) { idsEstaciones <- pathsEstaciones
  } else { 
    idsEstaciones <- estaciones[,colIdEstaciones]
    pathsEstaciones <- pathsEstaciones[match(idsEstaciones, pathsEstaciones)]
  }
  datos1 <- read.table(paste(carpetaSalidaCDT, '/', pathsEstaciones[1], '/', pathsEstaciones[1], '.txt', sep=''))
  
  fechas <- strptime(datos1[,1], format = '%Y%m%d', tz = 'UTC')
  
  res <- matrix(NA, nrow = nrow(datos1), ncol=length(pathsEstaciones))
  res[, 1] <- datos1[,2]
  row.names(res) <- format(fechas, format = '%Y-%m-%d')
  colnames(res) <- pathsEstaciones
  
  i <- 2
  for (i in seq.int(from = 2, to = length(pathsEstaciones), by = 1)) {
    datosI <- read.table(paste(carpetaSalidaCDT, '/', pathsEstaciones[i], '/', pathsEstaciones[i], '.txt', sep=''))
    res[, i] <- datosI[, 2]
  }
  
  return(list(fechas=fechas, datos=res))
}

leerDatosVarsEstacionesFechasDeJSON <- function(pathArchivoDatos, formatoFechas='YmdHMSz!*', truncated=6) {
  # retorna una lista con 4 componentes:
  # - estaciones es un data.frame con la información de las nE estaciones. Lat, Long, automática, pluviométrica, etc.
  # str(datos$estaciones)
  # - fechas es un vector de POSIXct con las nT fechas de los datos
  # str(datos$fechas)
  # - variables es un data.frame con la información de las nV variables. Periodicidad, unidad, tipo de medida, etc.
  # str(datos$variables)
  # - datos es una lista de matrices de largo nV, con los datos de cada variable en variables
  # - Si todas las variables tienen la misma periodicidad, cada matriz es de tamaño nT x nE. 
  # - Las filas de las matrices corresponden a distintas fechas y las columnas a distintas estaciones.
  # - datos[[iV]][iT, iE] contiene la observación de la variable iV en variables, la fecha iT en fechas y la 
  # - estación iE en estaciones
  # - Si hay diferentes periodicidades, las variables con frecuencias de medición menores tendrán su nT más 
  # - pequeño. Usando rownames(datos[[iV]]) se obtienen las fechas correspondientes
  # str(datos$observaciones$datos)
  
  datos <- fromJSON(pathArchivoDatos)
  datos$fechas <- parse_date_time(datos$fechas, orders = formatoFechas, truncated = truncated)

  names(datos$observaciones$datos) <- datos$variables$variable
  
  datos$estaciones$NombreEstacionR <- make.names(datos$estaciones$NombreEstacion, unique = TRUE)
  datos$estaciones <- datos$estaciones[, c(ncol(datos$estaciones),1:(ncol(datos$estaciones)-1))]
  strsFechas <- as.character(datos$fechas)
  
  i <- 1
  for (i in seq_along(datos$observaciones$datos)) {
    datos$observaciones$datos[[i]] = t(datos$observaciones$datos[[i]])
    colnames(datos$observaciones$datos[[i]]) <- datos$estaciones$nombreEstacionR
    rownames(datos$observaciones$datos[[i]]) <- strsFechas[datos$observaciones$iFechas[[i]]+1]
  }
  return(datos)
}

extraerVariablesFechasDeEstacion <- function(datosVarsEstacionesFechas, estacion) {
  # Una estacion, fechas en filas, variables en columnas
  return(do.call(cbind, lapply(datosVarsEstacionesFechas, function(x, estacion) { return(x[, estacion, drop=T]) }, estacion)))
}

extraerVariablesEstacionesDeFecha <- function(datosVarsEstacionesFechas, fecha) {
  # Una fecha, variables en filas, estaciones en columnas
  if (!(is.integer(fecha) | is.numeric(fecha))) { fecha <- as.character(fecha) }
  
  return(do.call(rbind, lapply(seq_along(datosVarsEstacionesFechas), function(x, datosVarsEstacionesFechas, fecha) {
    res = datosVarsEstacionesFechas[[x]][fecha, , drop=F]
    rownames(res) <- names(datosVarsEstacionesFechas)[x]
    return(res)
  }, datosVarsEstacionesFechas, fecha)))
}
