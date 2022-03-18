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
if (iFrame >= 3) { script.dir.ReadGrADS <- sys.frame(iFrame - 3)$ofile 
} else { script.dir.ReadGrADS <- NULL }
while ((is.null(script.dir.ReadGrADS) || is.na(regexpr('ReadGrADS.r', script.dir.ReadGrADS, fixed=T)[1])) && iFrame >= 0) {
  script.dir.ReadGrADS <- sys.frame(iFrame)$ofile
  iFrame <- iFrame - 1
}
if (is.null(script.dir.ReadGrADS)) { script.dir.ReadGrADS <- ''
} else { script.dir.ReadGrADS <- paste0(dirname(script.dir.ReadGrADS), '/') }

source(paste0(script.dir.ReadGrADS, '../instalarPaquetes/instant_pkgs.r'), encoding = 'WINDOWS-1252')
instant_pkgs(c('stringi', 'rgdal', 'sp', 'akima', 'RANN'))

formaWGS84 <- '+ellps=WGS84 +datum=WGS84 +units=m'
formaEsferica <- '+a=6370000 +towgs84=0,0,0,0,0,0,0 +units=m'

formaDeLaTierra <- formaWGS84

parsePDef <- function(PDef, formaTierra=formaDeLaTierra) {
  iSize <- as.integer(PDef[2])
  jSize <- as.integer(PDef[3])
  tipoProyeccion <- PDef[4]
  if (tipoProyeccion %in% c('lcc', 'lccr')) {
    latRef <- as.numeric(PDef[5]);
    lonRef <- as.numeric(PDef[6]);
    iRef <- as.numeric(PDef[7]);
    jRef <- as.numeric(PDef[8]);
    STrueLat <- as.numeric(PDef[9]);
    NTrueLat <- as.numeric(PDef[10]);
    sLon <- as.numeric(PDef[11]);
    dx <- as.numeric(PDef[12]);
    dy <- as.numeric(PDef[13]);
    # WRF no usa el datum  WGS84 sino que usa un datum esférico
    # Más información acá: http://www.pkrc.net/wrf-lambert.html
    proj4string <- paste0('+proj=lcc +lat_1=', STrueLat, ' +lat_2=', NTrueLat, 
                         ' +lat_0=', latRef, ' ', ' +lon_0=', sLon, ' ', 
                         formaTierra, ' +no_defs')
    
    # TODO: agregar parámetro SRS_string al CRS y cambiar proj4string por wkt
    puntoRef <- SpatialPoints(coords = matrix(c(lonRef, latRef), nrow = 1), 
                              proj4string = CRS(paste0('+proj=longlat ', formaTierra, ' +no_defs')))
    puntoRef <- spTransform(puntoRef, proj4string)
    proj4string <- paste0('+proj=lcc +lat_1=', STrueLat, ' +lat_2=', NTrueLat, 
                         ' +lat_0=', latRef, ' ', ' +lon_0=', sLon, ' ', formaTierra, 
                         ' +x_0=', -sp::coordinates(puntoRef)[1,1], 
                         ' +y_0=', -sp::coordinates(puntoRef)[1,2], ' +no_defs')
    
    return (list(tipoProyeccion=tipoProyeccion, iSize=iSize, jSize=jSize, latRef=latRef, lonRef=lonRef, 
                 iRef=iRef, jRef=jRef, STrueLat=STrueLat, NTrueLat=NTrueLat, sLon=sLon, dx=dx, dy=dy, 
                 proj4string=proj4string))
  } else {
    stop(paste0('Error leyendo el campo PDEF. Tipo de proyección no implementada "', tipoProyeccion, '"'))
  }
}

parseDimDef <- function(dimDef, allowGaussian = FALSE) {
  # xdef 1440 linear -179.875 0.25
  # ydef 480  linear -59.875 0.25
  # zdef 1 levels 1
  dimDef[3] <- tolower(dimDef[3])
  if (dimDef[3] == 'linear') {
    res <- list(type='linear', from=as.numeric(dimDef[4]), by=as.numeric(dimDef[5]), n=as.integer(dimDef[2]))
    res$vals <- seq(from=as.numeric(dimDef[4]), by=as.numeric(dimDef[5]), length.out=as.integer(dimDef[2]))
    return (res)
  } else if (dimDef[3] == 'levels') {
    res <- list(type='levels', nLevels=as.integer(dimDef[2]))
    res$vals <- as.numeric(dimDef[4:length(dimDef)])
    return (res)    
  } else {
    # implement other types
    stop(paste0('parseDimDef: Error unknown dimdef type "', dimDef[3], '"'))
  }
}

parseGradsAbsoluteTime <- function(s) {
  # s <- TDef[4]
  # hh:mmZddmmmyyyy 
  posZ <- regexpr("[Zz]", s)[[1]]
  if (posZ > 0) {
    strDiaMesAnio <- substr(s, posZ + 1, nchar(s))
    if (posZ == 6) {
      # hh:mmZ...
      # 123456
      hora <- as.numeric(substring(s, 1, 2))
      minuto <- as.numeric(substring(s, 4, 5))
    } else if (posZ == 3)  {
      # hhZ...
      # 123
      hora <- as.numeric(substring(s, 1, 2))
      minuto <- 0
    } else { stop(paste0('Fecha inválida: "', s, '"')) }
  } else  {
    hora <- 0
    minuto <- 0
    strDiaMesAnio <- s
  }
  
  strs3DigitosMeses <- c('jan', 'feb', 'mar', 'apr', 'may', 'jun', 'jul', 'aug', 'sep', 'oct', 'nov', 'dec')
  
  if (nchar(strDiaMesAnio) == 9) {
    # ddmmmyyyy
    # 123456789
    dia <- as.numeric(substring(strDiaMesAnio, 1, 2))
    mes <- which(tolower(substring(strDiaMesAnio, 3, 5)) == strs3DigitosMeses)
    anio <- as.numeric(substring(strDiaMesAnio, nchar(strDiaMesAnio)-3, nchar(strDiaMesAnio)))
  } else if (nchar(strDiaMesAnio) == 8) {
    # dmmmyyyy
    # 12345678
    dia <- as.numeric(substring(strDiaMesAnio, 1, 1))
    mes <- which(substring(strDiaMesAnio, 2, 4) == strs3DigitosMeses)
    anio <- as.numeric(substring(strDiaMesAnio, nchar(strDiaMesAnio) - 3, nchar(strDiaMesAnio)))
  }  else if (nchar(strDiaMesAnio) == 7) {
    # ddmmmyy o mmmyyyy
    # 1234567   1234567
    firstChar <- substring(strDiaMesAnio[1], 1, 1)
    if (grepl("^[[:digit:]]", firstChar)) {
      dia <- as.numeric(substring(strDiaMesAnio, 1, 2))
      mes <- which(substring(strDiaMesAnio, 3, 5) == strs3DigitosMeses)
      anio <- as.numeric(substring(strDiaMesAnio, nchar(strDiaMesAnio) - 1, nchar(strDiaMesAnio))) + 1950
    } else {
      dia <- 1
      mes <- which(substring(strDiaMesAnio, 3, 5) == strs3DigitosMeses)
      anio <- as.numeric(substring(strDiaMesAnio, nchar(strDiaMesAnio) - 3, nchar(strDiaMesAnio)))
    }
  }  else if (nchar(strDiaMesAnio) == 6) {
    # dmmmyy
    # 123456
    dia <- as.numeric(substring(strDiaMesAnio, 1, 1))
    mes <- which(substring(strDiaMesAnio, 3, 5) == strs3DigitosMeses)
    anio <- as.numeric(substring(strDiaMesAnio, nchar(strDiaMesAnio) - 1, nchar(strDiaMesAnio))) + 1950
  } else if (nchar(strDiaMesAnio) == 5) {
    # mmmyy
    # 12345
    dia <- 1
    mes <- which(substring(strDiaMesAnio, 3, 5) == strs3DigitosMeses)
    anio <- as.numeric(substring(strDiaMesAnio, nchar(strDiaMesAnio) - 1, nchar(strDiaMesAnio))) + 1950
  }
  
  return (ISOdatetime(year=anio, month=mes, day=dia, hour=hora, min=minuto, sec=0, tz='GMT'))
}

parseTDef <- function(TDef) {
  # 1    2    3      4            5
  # tdef 9999 linear 21Z23sep2006 3hr
  res <- list(type='linear')
  n <- as.integer(TDef[2])
  from <- parseGradsAbsoluteTime(s = TDef[4])
  increment <- as.numeric(substr(TDef[5], 1, nchar(TDef[5]) - 2))
  unit <- tolower(substr(TDef[5], start=nchar(TDef[5])-1, nchar(TDef[5])))
  
  if (unit == 'mn') { deltaDT <- increment * 60
  } else if (unit == 'hr') { deltaDT <- increment * 60 * 60
  } else if (unit == 'dy') { deltaDT <- increment * 60 * 60 * 24
  } else {
    deltaDT <- 0
    stop(paste0('parseTDef: Error unknown unit "', unit, '"'))
  }
  
  res$from <- from
  res$by <- deltaDT
  res$n <- n
  res$vals <- seq(from=from, by=deltaDT, length.out=n)
  res$increment <- increment
  res$unit <- unit
  return(res)
}

parseVarDef <- function(varDef) {
  varStrSplit <- unlist(strsplit(varDef[3], ","))
  if (varStrSplit[1] == '-1') {
    varUnits <- varStrSplit[1]
    varStructure <- as.numeric(varStrSplit[2])
    
    if (length(varStrSplit) > 2) {
      varArg1 <- as.numeric(varStrSplit[3])
      if (length(varStrSplit) > 3) {
        varArg2 <- as.numeric(varStrSplit[4])
      }  else {
        varArg2 <- ''
      }
    } else {
      varArg1 <- ''
      vararg2 <- ''
    } 
  } else {
    varUnits <- varDef[3]
    varStructure <- ''
    varArg1 <- ''
    varArg2 <- ''
  }
  return (list(varname=varDef[1], levs=as.numeric(varDef[2]), units=varUnits, structure=varStructure, 
               description=paste(varDef[4:length(varDef)], collapse=' '), arg1=varArg1, arg2=varArg2))
}

getIFecha <- function(ctl, fecha) {
  return (as.integer((difftime(fecha, ctl$tdef$from, 'days') / ctl$tdef$by)))
}

readXYGrid <- function(ctl, dsetOverride=NA, idxFecha=1, idxVar=1, idxNivelZ=1) {
  if (is.na(dsetOverride)) {
    if (startsWith(ctl$dset, '^')) { 
      if (.Platform$OS.type == "windows") {
        binFile <- paste0(dirname(ctl$ctlFile), '/', gsub(pattern = ':', x = substr(ctl$dset, 2, nchar(ctl$dset)), replacement = '', fixed = T), collapse='')
      } else {
        binFile <- paste0(dirname(ctl$ctlFile), '/', substr(ctl$dset, 2, nchar(ctl$dset)), collapse='')  
      }
    } else { 
      if(.Platform$OS.type == "windows") {
        binFile <- gsub(pattern = ':', x = ctl$dset, replacement = '', fixed = T)
      } else {
        binFile <- ctl$dset
      }
    }
  } else { binFile <- dsetOverride }

  binFile <- file(binFile, open="rb")
  iPos <- ctl$bytesHastaHiperCuboT[idxFecha] + ctl$bytesHastaVariableINivelZ[[idxVar]][idxNivelZ]
  if (iPos > 0) { seek(binFile, iPos) }
    
  if (ctl$vars[[idxVar]]$structure == 40) {
    signed <- ctl$vars[[idxVar]]$arg2 != ''
    bytesPorCelda <- ctl$vars[[idxVar]]$arg1
    what="integer"
  } else {
    signed <- T
    bytesPorCelda <- 4
    what="numeric"
  }

  gridData <- matrix(readBin(con = binFile, what = what, n = ctl$numCells, size = bytesPorCelda, 
                             signed = signed, endian = ctl$endian), 
                     nrow = ctl$numRows, ncol = ctl$numCols)
  close(binFile)
  
  gridData[gridData == ctl$undef] <- NA
  
  if ('yrev' %in% ctl$options) gridData <- gridData[,ncol(gridData):1]
  
  return (gridData)
}

getGrillaNativa <- function(ctl) {
  if (!is.null(ctl$pDef)) {
    if (ctl$pDef$tipoProyeccion %in% c('lcc', 'lccr')) {
      #coords <- expand.grid(x=seq(from = 0, by = ctl$pDef$dx, length.out = ctl$pDef$iSize), 
      #                      y=seq(from = 0, by = ctl$pDef$dy, length.out = ctl$pDef$jSize))
      #sp::coordinates(coords) <- c('x', 'y')
      #grilla <- points2grid(coords)
      
      grilla <- GridTopology(cellcentre.offset = c(ctl$pDef$dx * (-ctl$pDef$iRef + 1), 
                                                   -ctl$pDef$dy * (-ctl$pDef$jRef + 1)),
                             cellsize = c(ctl$pDef$dx, -ctl$pDef$dy), 
                             cells.dim = c(ctl$pDef$iSize, ctl$pDef$jSize))
      
      names(grilla@cellsize) <- c('x', 'y')
      names(grilla@cellcentre.offset) <- names(grilla@cellsize)
      names(grilla@cells.dim) <- names(grilla@cellsize)
      
      return(SpatialGrid(grid = grilla, proj4string = CRS(ctl$pDef$proj4string)))
    } else
      stop(paste0('getGrillaNativa: Error unimplemented GrADS PDef type "', ctl$pDef$tipoProyeccion, '"'))
  } else return(getGrillaRectilineaLatLong(ctl))
}

getGrillaRectilineaLatLong <- function(ctl, formaTierra=formaDeLaTierra) {
  coords <- expand.grid(lon=ctl$xdef$vals, lat=ctl$ydef$vals)
  sp::coordinates(coords) <- c('lon', 'lat')
  grilla <- points2grid(coords)
  # proj4string <- CRS('+proj=longlat +ellps=WGS84 +no_defs')
  # WRF no usa el datum  WGS84 sino que usa un datum esférico
  # Más información acá: http://www.pkrc.net/wrf-lambert.html
  return(SpatialGrid(grid = grilla, proj4string = CRS(paste0('+proj=longlat ', formaTierra, ' +no_defs'))))
}

sampleSP <- function(origSP_DF, newSP_DF, samplingMethod = c('NearestNeighbor', 'Bilinear', 'Bicubic'), 
                     zcolOrigSP=1, zcolNewSP=1) {
  # TO-DO
  stop('TO-DO')
  samplingMethod <- 'NearestNeighbor'
  samplingMethod <- 'Bilinear'
  samplingMethod <- 'Bicubic'
  # origSP_DF <- readXYGridSP(ctl, idxFecha = 1, idxVar = 18, idxNivelZ = 1)
  # newSP_DF <- getGrillaRectilineaLatLong(ctl)
  # newSP_DF <- SpatialGridDataFrame(grid = newSP_DF, data = data.frame(value=rep(NA_real_, length(newSP_DF))))
  if (!is(object = origSP_DF, class2 = 'SpatialPointsDataFrame')) {
    puntosOrigSP <- SpatialPointsDataFrame(
      sp::coordinates(origSP_DF), data = origSP_DF@data, proj4string = origSP_DF@proj4string)
  } else { puntosOrigSP <- origSP_DF }
  
  if (!is(object = newSP_DF, class2 = 'SpatialPointsDataFrame')) {
    samplePoints <- SpatialPointsDataFrame(
      sp::coordinates(newSP_DF), data = data.frame(value=rep(NA_real_, length(newSP_DF))), 
      proj4string = newSP_DF@proj4string)
  } else { samplePoints <- newSP_DF }
  
  
  if (!identicalCRS(puntosOrigSP, samplePoints)) {
    samplePoints <- spTransform(x = samplePoints, CRSobj = puntosOrigSP@proj4string)
  }
  
  #samplingMethod <- 'Bicubic'
  if (samplingMethod[1] == 'NearestNeighbor') {
    idxs <- nn2(data = sp::coordinates(puntosOrigSP), query = sp::coordinates(samplePoints), k = 1)
    newSP_DF@data[,zcolNewSP] <- puntosOrigSP@data[idxs$nn.idx, zcolOrigSP]
  } else if (samplingMethod[1] == 'Bilinear') {
    vals <- interpp(x= puntosOrigSP, z=colnames(puntosOrigSP@data)[zcolOrigSP], xo = samplePoints, 
                    linear = TRUE)
    newSP_DF@data[,zcolNewSP] <- as.numeric(vals$z)
  } else if (samplingMethod[1] == 'Bilinear') {
    coordsOrig <- sp::coordinates(puntosOrigSP)
    coordsSamplePoints <- sp::coordinates(samplePoints)
    
    help(bilinear)
    
    gridded(origSP_DF)
    
    x <- unique(coordsOrig[,1])
    y <- unique(coordsOrig[,2])
    z <- matrix(data = puntosOrigSP@data[,zcolOrigSP], nrow = length(x), ncol = length(y))
    
    vals <- bicubic(x = x, y = y, z = z, 
                    x0 = coordsSamplePoints[,1], y0 = coordsSamplePoints[,2])
    
    vals <- interpp(x= puntosOrigSP, z=colnames(puntosOrigSP@data)[zcolOrigSP], 
                    xo = samplePoints, linear = FALSE)
    vals <- interp(x= puntosOrigSP, z=colnames(puntosOrigSP@data)[zcolOrigSP], 
                   xo = samplePoints, linear = FALSE)
    length(vals)
    newSP_DF@data[,zcolNewSP] <- vals@data[,1]
  }
  
  class(newSP_DF)
  
  mapearGrillaGGPlot(grilla = newSP_DF, continuo=TRUE, dibujar=F, titulo = samplingMethod[1], zcol = zcolNewSP)
  
  return(newSP_DF)
}

extraerSeriesTemporales <- function(ctl, dsetOverride=NA, tIni=1, tFin=ctl$tdef$n, idxVar=1, idxNivelZ=1, 
                                    spUbicaciones, samplingMethod = c('NearestNeighbor', 'Bilinear', 'Bicubic'), 
                                    zColNames=1) {
  if (!is.null(ctl$pDef)) {
    spUbicaciones <- spTransform(spUbicaciones, CRS(ctl$pDef$proj4string))
    grillaXY <- getGrillaNativa(ctl)
    
    x <- unique(sp::coordinates(grillaXY)[,1])
    y <- unique(sp::coordinates(grillaXY)[,2])
  } else {
    x <- ctl$xdef$vals
    y <- ctl$ydef$vals
  }
  
  coords <- sp::coordinates(spUbicaciones)
  
  tis <- seq.int(from = tIni, to = tFin, by = 1)
  res <- matrix(data = NA_real_, nrow = length(tis), ncol = length(spUbicaciones))
  
  for (i in seq_along(tis)) {
    ti <- tis[i]
    z <- readXYGrid(ctl, dsetOverride = dsetOverride, ti, idxVar = idxVar, idxNivelZ = idxNivelZ)
    
    if (samplingMethod[1] == 'NearestNeighbor') {
      idxs <- nn2(data = expand.grid(x, y), query = coords, k = 1)
      res[i, ] <- z[idxs$nn.idx]
    } else if (samplingMethod[1] == 'Bilinear') {
      res[i, ] <- bilinear(x = x, y = y, z = z, x0 = coords[,1], y0 = coords[,2])$z
    } else if (samplingMethod[1] == 'Bicubic') {
      res[i, ] <- bicubic(x = x, y = y, z = z, x0 = coords[,1], y0 = coords[,2])$z
    }
  }

  colnames(res) <- spUbicaciones@data[,zColNames]
  rownames(res) <- as.character(ctl$tdef$vals)
  return(res)
}

readXYGridSP <- function(ctl, dsetOverride=NA, idxFecha=1, idxVar=1, idxNivelZ=1, grillaXY=getGrillaNativa(ctl)) {
  # ctl$pDef$proj4string
  # ctl$pDef$proj4string <- "+proj=lcc +lat_1=-39.765 +lat_2=-39.765 +lat_0=-33.171  +lon_0=-65.478 +a=6370000 +units=m +x_0=-706587.442842593 +y_0=29740.3609605322 +no_defs"
  # ctl$pDef$proj4string <- "+proj=lcc +lat_1=-39.765 +lat_2=-39.765 +lat_0=-33.171  +lon_0=-65.478 +ellps=WGS84 +datum=WGS84 +units=m +x_0=-706587.442842593 +y_0=29740.3609605322 +no_defs"
  # dsetOverride <- 'F:/ADME/precip_rionegro/datos/satelites/GSMaP/originales/gsmap_gauge.20191214.1000.dat'
  datos <- readXYGrid(ctl = ctl, dsetOverride = dsetOverride, idxFecha = idxFecha, idxVar = idxVar, idxNivelZ = idxNivelZ)

  if (ctl$convert360to180) {
    mid <- nrow(datos) %/% 2
    idxLeft <- 1:mid
    idxRight <- (mid + 1):nrow(datos)
    aux <- datos[idxLeft, ]
    datos[idxLeft, ] <- datos[idxRight, ]
    datos[idxRight, ] <- aux
  }
  
  if (is(grillaXY, 'SpatialGrid')) {
    res <- SpatialGridDataFrame(grid = grillaXY, data = data.frame(value=as.numeric(datos)))
  } else if (is(grillaXY, 'SpatialPixels')) {
    res <- SpatialPixelsDataFrame(points = grillaXY, data = data.frame(value=as.numeric(datos)))
  } else {
    stop(paste0('ReadGrADS.readXYGridSP: unknown grid class ', class(grillaXY)))
  }

  return(res)

  if (F) {
    # Obtengo los centroides de la grilla nativa
    puntosGrillaNativa <- SpatialPointsDataFrame(
      sp::coordinates(grillaXY), data = data.frame(value=as.numeric(datos)), 
      proj4string = grillaXY@proj4string)
    
    if (F) {
      lats <- readXYGrid(ctl = ctl, dsetOverride = dsetOverride, idxFecha = idxFecha, idxVar = 1, idxNivelZ = idxNivelZ)
      longs <- readXYGrid(ctl = ctl, dsetOverride = dsetOverride, idxFecha = idxFecha, idxVar = 2, idxNivelZ = idxNivelZ)
      
      prjLongLatWRF <- paste0('+proj=longlat ', formaDeLaTierra, ' +no_defs')
      # TODO: el SRS_string está hardcodeado para formaDeLaTierra WGS84
      SRS_string <- "EPSG:4326"
      
      puntosGrillaNativa <- spTransform(
        puntosGrillaNativa, CRS(projargs = prjLongLatWRF, SRS_string = SRS_string))
      sp::coordinates(puntosGrillaNativa)[c(1, 174), ]
      matrix(c(range(longs)[1], range(lats)[1], range(longs)[2], range(lats)[2]), nrow = 2, byrow = TRUE)
      matrix(c(longs[1, 1], lats[1, 1], longs[174, 174], lats[174, 174]), nrow = 2, byrow = TRUE)
      
      # Esto tiene que dar (0, 0). Es transformar el centro de la grilla de latlong a la pdef
      sp::coordinates(spTransform(
        SpatialPoints(
          coords = matrix(
            c(round((longs[trunc(ctl$pDef$iRef), trunc(ctl$pDef$jRef)] + longs[ceiling(ctl$pDef$iRef), ceiling(ctl$pDef$jRef)]) / 2, 3),
              round((lats[trunc(ctl$pDef$iRef), trunc(ctl$pDef$jRef)] + lats[ceiling(ctl$pDef$iRef), ceiling(ctl$pDef$jRef)]) / 2, 3)), nrow = 1, byrow = TRUE),
          proj4string = CRS(projargs = prjLongLatWRF, SRS_string = SRS_string)),
        grillaXY@proj4string))
      
      # Estos dos tienen que dar iguales, es el inverso, convertir las coordenadas 0,0 en pdef a latlong
      sp::coordinates(spTransform(
        SpatialPoints(
          coords = matrix(c(0, 0), nrow = 1, byrow = TRUE), 
          proj4string = grillaXY@proj4string),
        grillaXY@proj4string))
      c(ctl$pDef$lonRef, ctl$pDef$latRef)
      
      # Esto tiene que dar
      # (-870000, -870000)
      # (870000, 870000)
      ctl$pDef
      (ctl$pDef$dx * ctl$pDef$iSize) / 2
      (ctl$pDef$dy * ctl$pDef$jSize) / 2
      sp::coordinates(spTransform(
        SpatialPoints(coords = matrix(c(longs[1, 1], lats[1, 1], longs[174, 174], lats[174, 174]), nrow = 2, byrow = TRUE), 
                      proj4string = CRS(projargs = prjLongLatWRF, SRS_string = SRS_string)),
        grillaXY@proj4string))
      
      proj4string(grillaXY)
      prjLongLatWRF
      
      coords <- spTransform(
        SpatialPoints(
          coords = matrix(c(as.numeric(longs), as.numeric(lats)), ncol = 2), 
          proj4string = CRS(projargs = prjLongLatWRF, SRS_string = SRS_string)), 
        grillaXY@proj4string)
      matrix(c(as.numeric(longs), as.numeric(lats)), ncol = 2) - sp::coordinates(coords)
      lala <- SpatialGridDataFrame(
        grid = grillaXY, 
        data = data.frame(value=matrix(c(as.numeric(longs), as.numeric(lats)), ncol = 2) - sp::coordinates(coords)))
      mapearGrillaGGPlot(lala, continuo = T, dibujar = F)
      
      lala <- sp::coordinates(spTransform(
        SpatialPoints(
          coords = matrix(c(longs, lats), ncol = 2), 
          proj4string = CRS(projargs = prjLongLatWRF, SRS_string = SRS_string)),
        grillaXY@proj4string))
      
      diff(sort(unique(lala[,1])))
      diff(sort(unique(lala[,2])))
    }
    
    # Obtengo la grilla LatLong y sus centroides
    grillaLatLong <- getGrillaRectilineaLatLong(ctl)
    puntosGrillaLatLong <- SpatialPointsDataFrame(
      coords = sp::coordinates(grillaLatLong), 
      proj4string = grillaLatLong@proj4string,
      data = data.frame(value=rep(NA_real_, length(grillaLatLong))))
    
    # Proyecto los centroides de la grilla LatLong a la proyección de la grilla nativa
    puntosGrillaLatLong <- spTransform(puntosGrillaLatLong, grillaXY@proj4string)
    
    # Obtengo los valores interpolados en la grilla de LatLong
    vals <- interpp(x= puntosGrillaNativa, z='value', xo = puntosGrillaLatLong, linear = TRUE)
    
    grillaLatLongDF <- SpatialGridDataFrame(grid = grillaLatLong, data = data.frame(value=vals$z))
    length(grillaXY)
    length(grillaLatLong)
    
    escala <- crearEscalaTemperaturaWRF()
    rc <- crop(raster(grillaLatLongDF), extent(-60, -52, -36, -29))
    g <- as(rc, 'SpatialGridDataFrame')
    shpBase <- cargarSHP('C:/mch/ArchivosProcesosLocales/MapaUruguayVacio/uruguay_departamentos.shp')
    shpBase <- spTransform(shpBase, g@proj4string)
    mapearGrillaGGPlot(grilla = g, shpBase = shpBase, escala = escala, dibujar=F, 
                       titulo = paste0('GLL: ', proj4string(g)), 
                       subtitulo = paste0('GN: ', proj4string(grillaXY)))
    
    mapearGrillaGGPlot(grillaLatLongDF, shpBase = shpBase, escala = escala, dibujar=F)
    mapearGrillaGGPlot(grillaLatLongDF, shpBase = shpBase, continuo = T, isolineas = T, dibujar=F)
    
    shpBase <- spTransform(shpBase, res@proj4string)
    mapearGrillaGGPlot(res, shpBase = shpBase, escala = escala, dibujar=F)
    mapearGrillaGGPlot(res, escala = escala, dibujar=F)
  }
}

inicializarCTL <- function(ctl) {
  # GrADS views gridded data sets as multi-dimensional arrays varying in longitude, latitude, vertical level, 
  # variable, and time. 
  # In version 2.0, a fifth grid dimension was added. The fifth grid dimension is assumed to be used for ensembles, 
  # and so it is given the name E (or ens), but because it is generally implemented,  it may be used for other grid 
  # dimensions such as EOFs. The default size of the E dimension is 1 -- if no E dimension exists, it is not 
  # necessary to explicity declare it the descriptor file.
  
  # bytesPorHiperCuboT es la cantidad de bytes entre un paso de tiempo y otro. En la iteración se va incrementando
  # de a ctl$bytesPorGrillaLonLat[i] para construir ctl$bytesHastaVariableINivelZ
  ctl$bytesPorHiperCuboT <- ctl$tHeader
  # bytesPorGrillaLonLat es la cantidad de bytes en una grilla XY de un nivel fijo de elevación
  ctl$bytesPorGrillaLonLat <- integer(length=ctl$nVars)
  # bytesHastaVariableINivelZ[[i]][j] es la cantidad de bytes hasta el nivel de elevación j de la variable i desde
  # el comienzo de un hipercubo de tiempo
  ctl$bytesHastaVariableINivelZ <- list()
  i <- 1
  for (i in seq_along(ctl$vars)) {
    # Cada celda de la variable i ocupa bytesPorCeldaBytes
    if (ctl$vars[[i]]$structure == 40) { bytesPorCelda <- as.integer(ctl$vars[[i]]$arg1)
    } else { bytesPorCelda <- 4L }
    
    # Si hay pDef, cada grilla tiene ctl$pDef$iSize celdas en x por ctl$pDef$jSize en y, total = ctl$pDef$iSize * ctl$pDef$jSize
    # Si no hay pDef, tiene ctl$xdef$n en x y ctl$ydef$n en y, total = ctl$xdef$n * ctl$ydef$n
    if (!is.null(ctl$pDef)) { ctl$bytesPorGrillaLonLat[i] <- ctl$xyHeader + bytesPorCelda * ctl$pDef$iSize * ctl$pDef$jSize
    } else { ctl$bytesPorGrillaLonLat[i] <- ctl$xyHeader + bytesPorCelda * ctl$xdef$n * ctl$ydef$n }
    
    ctl$bytesHastaVariableINivelZ[[i]] <- integer(length=max(ctl$vars[[i]]$levs, 1))
    for (j in seq_along(ctl$bytesHastaVariableINivelZ[[i]])) {
      ctl$bytesHastaVariableINivelZ[[i]][j] <- ctl$bytesPorHiperCuboT
      ctl$bytesPorHiperCuboT <- ctl$bytesPorHiperCuboT + ctl$bytesPorGrillaLonLat[i]
    }
  }
  
  ctl$bytesPorHiperCuboT = ctl$bytesPorHiperCuboT + ctl$trailerBytes
  ctl$bytesHastaHiperCuboT <- ctl$fileHeader + (seq_along(ctl$tdef$vals) - 1) * ctl$bytesPorHiperCuboT
  
  return(ctl)
}

parseCTL <- function(ctlFile, convert360to180=FALSE, verbose=FALSE) {
  # verbose = TRUE
  # ctlFile <- 'D:/testsMCH/GrADS/wrfout_d02_2017-05-04_000000.ctl'
  ctl <- readChar(con = ctlFile, nchars = file.info(ctlFile)$size)
  
  atributosGrADS <- 'dset|chsub|dtype|index|stnmap|title|undef|unpack|fileheader|xyheader|xytrailer|theader|headerbytes|trailerbytes|xvar|yvar|zvar|stid|tvar|toffvar|cachesize|options|pdef|xdef|ydef|zdef|tdef|edef|vectorpairs|vars|endvars|\\@|\\*'
  expresionReg <- paste(paste0('(\n|^)', unlist(strsplit(atributosGrADS, split = '|', fixed = T))), collapse = '|')
  
  matches <- gregexpr(pattern = expresionReg, text = ctl, ignore.case = T)
  textoAtributos <- gsub(unlist(regmatches(x = ctl, m = matches)), pattern = '\n', replacement = '')
  largoMatches <- attr(x = matches[[1]], which = 'match.length')
  
  valoresAtributos <- character(length(textoAtributos))
  i <- 1
  for (i in seq.int(from = 1, to = length(matches[[1]]) - 1, by = 1)) {
    valoresAtributos[i] <- trimws(substr(ctl, start = matches[[1]][i] + largoMatches[i], stop = matches[[1]][i + 1] - 1)) }
  i <- length(matches[[1]])
  valoresAtributos[i] <- trimws(substr(ctl, start = matches[[1]][i] + largoMatches[i], stop = nchar(ctl)))
  
  rm(matches, largoMatches)
  
  # not all fields will be assigned to res, only those in the ctl file.
  # when using it you must check for for example res$dset != NULL
  res <- list(ctlFile=ctlFile)
  # these fields are an exception to the rule since they must always have a value to be able to read the bin
  res$fileHeader <- 0L
  res$tHeader <- 0L
  res$trailerBytes <- 0L
  res$xyHeader <- 0L
  
  # extensions
  res$convert360to180 <- convert360to180
  
  # partial support, more options need to be coded but this handles TRMM-3B42-RT and WRF ctl files 
  # most of the code comes from MCH's source code and is converted from Delphi to R 
  i <- 1
  # i <- 5
  # i <- i + 1
  # verbose<-T
  while (i <= length(textoAtributos)) {
    # lower case to check without case sensitivity
    param <- tolower(textoAtributos[i])
    if (verbose) print(paste(i, ': ', param, ' - ', valoresAtributos[i]))
    
    # ignore comments for now
    if (!startsWith(x = param, prefix = '*')) {
      if (param == 'dset') {
        res$dset <- valoresAtributos[i]
      } else if (param == 'title') {
        res$title <- valoresAtributos[i]
      } else if (param == 'undef') {
        res$undef <- as.numeric(valoresAtributos[i])
      } else if (param == 'fileheader') {
        res$fileHeader <- as.integer(valoresAtributos[i])
      } else if (param == 'xyheader') {
        res$xyHeader <- as.integer(valoresAtributos[i])
      } else if (param == 'headerbytes' || param == 'theader') {
        res$tHeader <- as.integer(valoresAtributos[i])
      } else if (param == 'trailerbytes') {
        res$trailerBytes <- as.integer(valoresAtributos[i])
      } else if (param == 'options') {
        res$options <- unlist(strsplit(valoresAtributos[i], split = '[[:space:]]+'))
      } else if (param == 'pdef') {
        res$pDef <- parsePDef(PDef = c('pdef', unlist(strsplit(x = valoresAtributos[i], split = '[[:space:]]+'))))
      } else if (param == 'xdef') {
        res$xdef <- parseDimDef(dimDef = c('xdef', unlist(strsplit(x = valoresAtributos[i], split = '[[:space:]]+'))))
        if (res$convert360to180) {
          if (res$xdef$type == 'linear') {
            res$xdef$from <- res$xdef$from - 180
            res$xdef$vals <- res$xdef$vals - 180
          } else {
            stop(paste0(
              'readGrADS.parseCTL: convert360to180 not implemented for xdef$type==', res$xdef$type))
          }
        }
      } else if (param == 'ydef') {
        res$ydef <- parseDimDef(dimDef = c('ydef', unlist(strsplit(x = valoresAtributos[i], split = '[[:space:]]+'))), allowGaussian = TRUE)
      } else if (param == 'zdef') {
        res$zdef <- parseDimDef(dimDef = c('zdef', unlist(strsplit(x = valoresAtributos[i], split = '[[:space:]]+'))))
      } else if (param == 'tdef') {
        res$tdef <- parseTDef(TDef = c('tdef', unlist(strsplit(x = valoresAtributos[i], split = '[[:space:]]+'))))
      } else if (param == 'vars') {
        varsDef <- unlist(strsplit(x = valoresAtributos[i], split = '\n+'))
        res$nVars <- as.integer(varsDef[1])
        varsDef <- varsDef[-1]
        res$vars <- vector(mode = "list", length = res$nVars)
        j <- 1
        for (j in seq_along(res$vars))
          res$vars[[j]] <- parseVarDef(varDef = unlist(strsplit(varsDef[[j]], split = '[[:space:]]+')))
        names(res$vars) <- lapply(res$vars, FUN = function(x) { return(x$varname) })
        # +1 for endvars
        i <- i + 1
      } else if (param == '*') {
        # no processing done for now
        res$comments[[length(res$comments)+1]] <- valoresAtributos[i]
      } else if (param == '@') {
        # no processing done for now
        res$attributes[[length(res$attributes)+1]] <- valoresAtributos[i]
      }else {
        stop(paste0('parseCTL: Error unknown GrADS Data Descriptor File Component in line ', i, ': "', ctl[[i]], '"'))
      }
    }
    i <- i + 1
  }
  
  res <- inicializarCTL(ctl = res)
  
  if ('big_endian' %in% res$options) { res$endian <- 'big' 
  } else if ('little_endian' %in% res$options) { res$endian <- 'little' 
  } else if ('byteswapped' %in% res$options) { res$endian <- 'big' 
  } else { res$endian <- 'little' }
  
  if (!is.null(res$pDef)) {
    res$numRows <- res$pDef$iSize
    res$numCols <- res$pDef$jSize
  } else {
    res$numRows <- res$xdef$n
    res$numCols <- res$ydef$n
  }
  res$numCells <- res$numRows * res$numCols
  
  return (res)
}

lambert_from_spherical <- function(xlon, ylat, a=6370, lat0, lon0, lat1, lat2, std_lon) {
  # Code partially based on the TERRAIN preprocessor for MM5 v2.0,
  # developed by Yong-Run Guo and Sue Chen, National Center for
  # Atmospheric Research, and Pennsylvania State University
  # 10/21/1993
  conv <- 57.29578
  
  # This part may be precomputed
  if (lat0 < 0) { sign = -1
  } else { sign = 1}
  
  pole <- 90
  if (abs(lat1) > 90) {
    lat1 <- 60
    lat2 <- 30
    lat1 <- sign * lat1
    lat2 <- sign * lat2
  }
  if (lat1 == lat2) { xn <- sin(lat2 / conv)
  } else {
    xn <- log10(cos(lat1 / conv)) - log10(cos(lat2 / conv))
    xn <- xn / (log10(tan((45 - sign * lat1 / 2) / conv)) - log10(tan((45 - sign * lat2 / 2) / conv)))
  }
  psi1 <- 90 - sign * lat1
  psi1 <- psi1 / conv
  if (lat0 < 0) {
    psi1 <- -psi1
    pole <- -pole
  }
  psi0 <- (pole - lat0) / conv
  xc <- 0
  yc <- -a / xn * sin(psi1) * (tan(psi0 / 2) / tan(psi1 / 2)) ^ xn
  
  1 / (tan(psi0 / 2) / tan(psi1 / 2) ^ abs(xn))
  
  # Actual computation for the specified location
  ylon <- xlon - std_lon
  if (ylon > 180) ylon <- ylon - 360;
  if (ylon < -180.) ylon <- ylon + 360;
  flp <- xn * ylon / conv
  psx <- (pole - ylat) / conv
  r <- -a / xn * sin(psi1) * (tan(psx / 2) / tan(psi1 / 2)) ^ xn
  if (lat0 < 0) {
    xloc <- r * sin(flp)
    yloc <- r * cos(flp)
  } else {
    xloc <- -r * sin(flp)
    yloc <- r * cos(flp)
  }
  
  xloc <- xloc - xc
  yloc <- yloc - yc
  ret <- c(xloc, yloc)
  return(ret)
}

