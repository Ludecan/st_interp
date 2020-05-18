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
if (iFrame >= 3) { script.dir.interpolarYMapearEx <- sys.frame(iFrame - 3)$ofile 
} else { script.dir.interpolarYMapearEx <- NULL }
while ((is.null(script.dir.interpolarYMapearEx) || is.na(regexpr('interpolarYMapearEx.r', script.dir.interpolarYMapearEx, fixed=T)[1])) && iFrame >= 0) {
  script.dir.interpolarYMapearEx <- sys.frame(iFrame)$ofile
  iFrame <- iFrame - 1
}
if (is.null(script.dir.interpolarYMapearEx)) { script.dir.interpolarYMapearEx <- ''
} else { script.dir.interpolarYMapearEx <- paste0(dirname(script.dir.interpolarYMapearEx), '/') }
source(paste0(script.dir.interpolarYMapearEx, 'interpolarEx.r'), encoding = 'WINDOWS-1252')
source(paste0(script.dir.interpolarYMapearEx, "mapearEx.r"), encoding = 'WINDOWS-1252')
source(paste0(script.dir.interpolarYMapearEx, 'funcionesAuxiliares.r'), encoding = 'WINDOWS-1252')
source(paste0(script.dir.interpolarYMapearEx, 'leerEscalas.r'), encoding = 'WINDOWS-1252')
source(paste0(script.dir.interpolarYMapearEx, '../sysutils/sysutils.r'), encoding = 'WINDOWS-1252')

createDefaultListaMapas <- function(
    paramsIyM, fechasObservaciones, nObservacionesTemporales=length(fechasObservaciones), 
    dibujarObservacionesEscalaFija=F, dibujarEscalaFija=T, 
    dibujarObservacionesEscalaAdaptada=F, dibujarEscalaAdaptada=F, 
    generarThumbnailFija=F, generarThumbnailAdaptada=F, 
    incluirIsolineaFija=F, incluirIsolineaAdaptada=F,
    dibujarPuntosObservacionesFija=F, dibujarPuntosObservacionesAdaptada=F,
    salvarGeoTiff=!is.null(paramsIyM$tlagsAR), salvarBin=F, titulo='', incluirSubtitulo=F, 
    recalcularSiYaExiste=T) {
  if (!is.null(fechasObservaciones)) {
    if (any(format(fechasObservaciones, format = "%H:%M") != "00:00")) { 
      formatoFechasArchivo <- "%Y_%m_%d_%H_%M"
      formatoFechasTitulo <- "%Y-%m-%d %H:%M"
    } else { 
      formatoFechasArchivo <- "%Y_%m_%d" 
      formatoFechasTitulo <- "%Y-%m-%d"
    }
    
    nombreArchivo <- pathParaGuardadoDeArchivos(filename=paste(paramsIyM$baseNomArchResultados, format(fechasObservaciones, format=formatoFechasArchivo), '.png', sep=''), pathProceso=paramsIyM$pathProceso)
  } else {
    nombreArchivo <- pathParaGuardadoDeArchivos(filename=paste(paramsIyM$baseNomArchResultados, 1:nObservacionesTemporales, '.png', sep=''), pathProceso=paramsIyM$pathProceso)
  }
    
  dibujarEscalaFija <- rep(dibujarEscalaFija, nObservacionesTemporales)
  dibujarEscalaAdaptada <- rep(dibujarEscalaAdaptada, nObservacionesTemporales)
  generarThumbnailFija <- rep(generarThumbnailFija, nObservacionesTemporales)
  generarThumbnailAdaptada <- rep(generarThumbnailAdaptada, nObservacionesTemporales)
  incluirIsolineasFija <- rep(incluirIsolineaFija, nObservacionesTemporales)
  incluirIsolineasAdaptada <- rep(incluirIsolineaAdaptada, nObservacionesTemporales)
  dibujarPuntosObservacionesFija <- rep(dibujarPuntosObservacionesFija, nObservacionesTemporales)
  dibujarPuntosObservacionesAdaptada <- rep(dibujarPuntosObservacionesAdaptada, nObservacionesTemporales)
  salvarGeoTiff <- rep(salvarGeoTiff, nObservacionesTemporales)
  salvarBin <- rep(salvarBin, nObservacionesTemporales)
  titulo <- paste(titulo, format(fechasObservaciones, format=formatoFechasTitulo))
  incluirSubtitulo <- rep(incluirSubtitulo, nObservacionesTemporales)
  recalcularSiYaExiste <- rep(recalcularSiYaExiste, nObservacionesTemporales)
  dibujarObservacionesEscalaFija <- rep(dibujarObservacionesEscalaFija, nObservacionesTemporales) 
  dibujarObservacionesEscalaAdaptada <- rep(dibujarObservacionesEscalaAdaptada, nObservacionesTemporales)
  
  return(data.frame(nombreArchivo, dibujarObservacionesEscalaFija, dibujarEscalaFija, dibujarObservacionesEscalaAdaptada,
                    dibujarEscalaAdaptada, generarThumbnailFija, generarThumbnailAdaptada, incluirIsolineasFija, 
                    incluirIsolineasAdaptada, dibujarPuntosObservacionesFija, dibujarPuntosObservacionesAdaptada,
                    salvarGeoTiff, salvarBin, titulo, incluirSubtitulo, recalcularSiYaExiste,
                    stringsAsFactors = FALSE))
}

mapearI <- function(ti, coordsObservaciones, fechasObservaciones, valoresObservaciones, pathsRegresores=NULL, valoresRegresoresSobreObservaciones=NULL, coordsAInterpolar, 
                    interpolacion, paramsIyM, shpMask, xyLims, listaMapas, espEscalaFija=NULL, espEscalaAdaptada=NULL,
                    puntosAResaltar=NULL) {
  if (listaMapas$incluirSubtitulo[ti]) {
    if (!is.null(pathsRegresores) && ncol(pathsRegresores) > 0) {
      subtitulo <- paste(paste(sum(!is.na(valoresObservaciones[ti,])), 'Observaciones'), paste(names(pathsRegresores), collapse = ' + '), sep=' + ')
    } else { subtitulo <- paste(sum(!is.na(valoresObservaciones[ti,])), 'Observaciones') }
  } else { subtitulo <- '' }
  
  coordsObservaciones$value <- as.numeric(valoresObservaciones[ti, ])
  
  nomArch <- listaMapas$nombreArchivo[ti]
  if (listaMapas$dibujarObservacionesEscalaFija[ti] || listaMapas$dibujarEscalaFija[ti] || listaMapas$generarThumbnailFija[ti]) {
    if (is.null(espEscalaFija)) {
      archEscala <- paste(paramsIyM$pathProceso, changeFileExt(basename(nomArch), '_escala.json'), sep='')
      if (file.exists(archEscala)) { espEscalaFija <- fromJSON(txt = archEscala)
      } else { espEscalaFija <- crearEspecificacionEscalaRelativaAlMinimoYMaximo(continuo = F, nDigitos = 1) }
    }
    escala <- darEscala(especificacion = espEscalaFija, valores = interpolacion$predictions@data[,1])
    
    if (listaMapas$dibujarObservacionesEscalaFija[ti]) {
      nombreObs <- appendToFileName(nomArch, '_obs')
      mapearPuntosGGPlot(puntos = coordsObservaciones, shpBase = shpMask$shp, escala = escala, nomArchResultados = nombreObs, 
                         xyLims = xyLims, dibujar=interactive(), titulo=listaMapas$titulo[ti], subtitulo=subtitulo)
    }
    
    if (listaMapas$dibujarEscalaFija[ti]) {
      mapearGrillaGGPlot(grilla=interpolacion$predictions, shpBase=shpMask$shp, escala=escala, nomArchResultados=nomArch,
                         xyLims=xyLims, dibujar=interactive(), titulo=listaMapas$titulo[ti], subtitulo=subtitulo,
                         isolineas = listaMapas$incluirIsolineasFija[ti], dibujarPuntosObservaciones = listaMapas$dibujarPuntosObservacionesFija[ti], 
                         coordsObservaciones = coordsObservaciones, puntosAResaltar = puntosAResaltar)
      
      #espEscalaFija <- crearEspecificacionEscalaEquiespaciadaParaLluvia(interpolacion$predictions@data[,1])
      #escala <- darEscala(especificacion = espEscalaFija, valores = interpolacion$predictions@data[,1])
      #mapearGrillaGGPlot(grilla=interpolacion$predictions, shpBase=shpMask$shp, escala=escala, nomArchResultados=nomArch,
      #                   xyLims=xyLims, dibujar=interactive(), titulo='Precipitación Acumulada 24 al 25 de mayo del 2017',
      #                   isolineas = F, dibujarPuntosObservaciones = T, 
      #                   coordsObservaciones=coordsObservaciones, puntosAResaltar = puntosAResaltar)
    }
    
    if (listaMapas$generarThumbnailFija[ti]) {
      nombreThumb <- appendToFileName(agregarCarpetaAlFinal(nomArch, carpeta='procesadas'), '_thumb')
      mapearGrillaGGPlot(grilla=interpolacion$predictions, shpBase=shpMask$shp, escala=escala, nomArchResultados=nombreThumb,
                         xyLims = xyLims, dibujarEscala=F, dibujarEjes=F)
    }
  }
  
  if (listaMapas$dibujarObservacionesEscalaAdaptada[ti] || listaMapas$dibujarEscalaAdaptada[ti] || listaMapas$generarThumbnailAdaptada[ti]) {
    if (is.null(espEscalaAdaptada)) {
      archEscala <- paste(paramsIyM$pathProceso, changeFileExt(basename(nomArch), '_escalaAdaptada.json'), sep='')
      if (file.exists(archEscala)) { espEscalaAdaptada <- fromJSON(txt = archEscala)
      } else { espEscalaAdaptada <- crearEspecificacionEscalaRelativaAlMinimoYMaximo(continuo = F, nDigitos = 1) }
    }
    escala <- darEscala(especificacion = espEscalaAdaptada, valores = interpolacion$predictions@data[,1])

    if (listaMapas$dibujarObservacionesEscalaAdaptada[ti]) {
      nombreObs <- appendToFileName(nomArch, '_adaptada_obs')
      mapearPuntosGGPlot(puntos = coordsObservaciones, shpBase = shpMask$shp, escala = escala, nomArchResultados = nombreObs, 
                         xyLims = xyLims, dibujar=interactive(), titulo=listaMapas$titulo[ti], subtitulo=subtitulo)
    }
    
    if (listaMapas$dibujarEscalaAdaptada[ti])
      mapearGrillaGGPlot(grilla=interpolacion$predictions, shpBase=shpMask$shp, escala=escala, nomArchResultados=appendToFileName(nomArch, '_adaptada'),
                         xyLims=xyLims, dibujar=interactive(), titulo=listaMapas$titulo[ti], subtitulo=subtitulo,
                         isolineas = listaMapas$incluirIsolineasAdaptada[ti], dibujarPuntosObservaciones = listaMapas$dibujarPuntosObservacionesAdaptada[ti], 
                         coordsObservaciones=coordsObservaciones, puntosAResaltar = puntosAResaltar)
    
    if (listaMapas$generarThumbnailAdaptada[ti]) {
      nombreThumb <- appendToFileName(agregarCarpetaAlFinal(nomArch, carpeta='procesadas'), '_thumbAdaptada')
      mapearGrillaGGPlot(grilla=interpolacion$predictions, shpBase=shpMask$shp, escala=escala, nomArchResultados=nombreThumb,
                         xyLims = xyLims, dibujarEscala=F, dibujarEjes=F)
    }
  }  
  
  # Comento código viejo. No se usaba para nada pero esto permitía hacer el mapa de los puntos de observaciones directamente o
  # el mapa de varianza de predicción. Si alguna vez interesa, habría que actualizar este código al nuevo formato 
  # de los parámetros de mapeo
  #if (length(especificacionesEscalasObservaciones) > 0) {
  #  escalas <- darEscalas(especificaciones = especificacionesEscalasObservaciones, valores = interpolacion$predictions@data[,1])
  #  if (paramsIyM$graficarObservaciones) {
  #    # iEscala <- 1
  #    if (paramsIyM$titulo != '') { titulo <- paste('Observaciones de', paramsIyM$titulo, fechasObservaciones[ti])
  #    } else { titulo <- ''}
  #    for (iEscala in seq_along(along.with = escalas)) {
  #      nomArch <- paste(paramsIyM$baseNomArchResultados, '_', format(fechasObservaciones[ti], format="%Y_%m_%d"), '_Obs.png', sep='')
  #      coordsObservaciones$value <- as.numeric(t(valoresObservaciones[ti,]))
  #      mapearPuntosGGPlot(puntos = coordsObservaciones, shpBase = shpMask$shp, nomArchResultados = nomArch, nDigitos = 1,
  #                         xyLims=xyLims, escala = escalas[[iEscala]], dibujar = interactive(), titulo=titulo, dibujarTexto=T, 
  #                         subtitulo=paste(sum(!is.na(valoresObservaciones[ti,])), 'Observaciones'))
  #    }
  #  }
  #  
  #  if (paramsIyM$graficarPrediccion) {
  #    # iEscala <- 1
  #    if (paramsIyM$titulo != '') { titulo <- paste(paramsIyM$titulo, fechasObservaciones[ti])
  #    } else { titulo <- ''}
  #    for (iEscala in seq_along(along.with = escalas)) {
  #      nomArch <- paste(paramsIyM$baseNomArchResultados, '_', format(fechasObservaciones[ti], format="%Y_%m_%d"), '.png', sep='')
  #      mapearGrillaGGPlot(grilla=interpolacion$predictions, shpBase=shpMask$shp, escala=escalas[[iEscala]], 
  #                         nomArchResultados=nomArch, xyLims=xyLims, dibujar=interactive(), titulo=titulo, 
  #                         subtitulo=subtitulo)
  #    }
  #  }
  #}

  #if (length(especificacionesEscalasVarianza) > 0) {
  #  if (paramsIyM$graficarVarianza && !is.null(interpolacion$campoVarianza)) {
  #    escalas <- darEscalas(especificaciones = especificacionesEscalasVarianza, interpolacion$predictions@data[,1])
  #    iEscala<-1
  #    if (paramsIyM$titulo != '') { titulo <- paste('Varianza de Interpolación de', paramsIyM$titulo, seriesfechasObservaciones[ti])
  #    } else { titulo <- ''}
  #    for (iEscala in seq_along(along.with = escalas)) {
  #      nomArch <- paste(paramsIyM$baseNomArchResultados, '_', format(fechasObservaciones[ti], format="%Y_%m_%d"), '_VarPredictor.png', sep='')
  #      mapearGrillaGGPlot(grilla = interpolacion$predictions, shpBase = shpMask$shp, escala = escalas[[iEscala]], nomArchResultados = nomArch, 
  #                         xyLims=xyLims, dibujar = interactive(), titulo=titulo, subtitulo=subtitulo)
  #    }
  #  }
  #}
}

interpolarYMapearI <- function(
    iTi, tsAInterpolar=1:nrow(valoresObservaciones), coordsObservaciones, fechasObservaciones, 
    valoresObservaciones, pathsRegresores=NULL, valoresRegresoresSobreObservaciones=NULL, 
    coordsAInterpolar, paramsIyM, shpMask, xyLims, listaMapas, returnInterpolacion=TRUE, 
    paramsParaRellenoRegresores=NULL, pathsRegresoresParaRellenoRegresores=NULL, espEscalaFija=NULL, 
    espEscalaAdaptada=NULL) {
  # tsAInterpolar=1:nrow(valoresObservaciones)
  # iTi <- 90
  # iTi <- which(as.character(fechasObservaciones[tsAInterpolar]) == '2018-10-24')
  ti <- tsAInterpolar[iTi]
  print(paste(ti, ': ', fechasObservaciones[ti], sep=''))

  nomArchGeoTiff <- changeFileExt(listaMapas$nombreArchivo[ti], '.tif')
  dir.create(dirname(listaMapas$nombreArchivo[ti]), showWarnings = F, recursive = T)
  
  if (!listaMapas$recalcularSiYaExiste[ti] && file.exists(nomArchGeoTiff) && file.info(nomArchGeoTiff)$size > 0) {
    existia <- evaluarConReintentos(interpolacion <- list(predictions = readGDAL(nomArchGeoTiff, silent = T)), segundosEntreIntentos = 1)
  } else { existia <- F }
  
  if (!existia) {
    # paramsIyM$modoDiagnostico <- T
    if (paramsIyM$modoDiagnostico) {
      if (is.POSIXct(fechasObservaciones)) { 
        paramsIyM$carpetaParaModoDiagnostico <- paste(dirname(listaMapas$nombreArchivo[ti]), '/', format(x = fechasObservaciones[ti], format = '%Y%m%d_%H%M'), '/', sep='')
        paramsIyM$strFecha <- format(fechasObservaciones[ti], '%Y-%m-%d %H:%M')
      } else { 
        paramsIyM$carpetaParaModoDiagnostico <- paste(dirname(listaMapas$nombreArchivo[ti]), '/', fechasObservaciones[ti], '/', sep='') 
        paramsIyM$strFecha <- fechasObservaciones[ti]
      }
      dir.create(paramsIyM$carpetaParaModoDiagnostico, recursive = T, showWarnings = F)
      paramsIyM$nombreModelo <- nombreModelo(params = paramsIyM, pathsRegresores = pathsRegresores)
      
      if (is.null(paramsIyM$especEscalaDiagnostico)) {
        paramsIyM$especEscalaDiagnostico <- crearEspecificacionEscalaRelativaAlMinimoYMaximo()
      }
    }
    # params <- paramsIyM
    # objParameters <- NULL
    interpolacion <- universalGridding(
      ti=ti, coordsObservaciones = coordsObservaciones, fechasObservaciones=fechasObservaciones, 
      valoresObservaciones = valoresObservaciones, pathsRegresores = pathsRegresores, 
      valoresRegresoresSobreObservaciones = valoresRegresoresSobreObservaciones, 
      coordsAInterpolar = coordsAInterpolar, params = paramsIyM, shpMask = shpMask, 
      paramsParaRellenoRegresores=paramsParaRellenoRegresores, 
      pathsRegresoresParaRellenoRegresores=pathsRegresoresParaRellenoRegresores)
    # mapearGrillaGGPlot(grilla=interpolacion$predictions, shpBase = shpMask$shp, continuo = T, dibujar=F)
    # mapearPuntosGGPlot(coordsObservaciones, shpBase = shpMask$shp, continuo = T, dibujar=F)
  }
  #setwd('C:/mch')
  #plotKML(interpolacion$predictions)
  
  nomArch <- listaMapas$nombreArchivo[ti]
  if (listaMapas$salvarBin[ti]) { 
    salvarInterpolacion(baseNomArchResultados=nomArch, interpolacion, formatoSalida='binary', 
                        salvarPrediccion=TRUE, salvarVarianza=FALSE) }
  
  if (!existia && listaMapas$salvarGeoTiff[ti]) { 
    writeGDAL(dataset = interpolacion$predictions, fname = nomArchGeoTiff, 
              options = c('COMPRESS=DEFLATE', 'PREDICTOR=2', 'ZLEVEL=9')) 
  }
  
  if (listaMapas$dibujarObservacionesEscalaFija[ti] || listaMapas$dibujarEscalaFija[ti] || listaMapas$generarThumbnailFija[ti] || 
      listaMapas$dibujarObservacionesEscalaAdaptada[ti] || listaMapas$dibujarEscalaAdaptada[ti] || listaMapas$generarThumbnailAdaptada[ti]) {
    # puntosAResaltar = paramsIyM$puntosAResaltar[[iTi]]
    mapearI(ti=ti, coordsObservaciones=coordsObservaciones, fechasObservaciones=fechasObservaciones,
            valoresObservaciones=valoresObservaciones, pathsRegresores=pathsRegresores, 
            valoresRegresoresSobreObservaciones=valoresRegresoresSobreObservaciones, 
            coordsAInterpolar=coordsAInterpolar, interpolacion=interpolacion, paramsIyM=paramsIyM, 
            shpMask=shpMask, xyLims=xyLims, listaMapas=listaMapas, espEscalaFija=espEscalaFija, 
            espEscalaAdaptada=espEscalaAdaptada, puntosAResaltar = paramsIyM$puntosAResaltar[[iTi]])
  }
  if (returnInterpolacion) { return(interpolacion)
  } else { return(NULL) }
}

interpolarYMapear <- function(coordsObservaciones, fechasObservaciones, valoresObservaciones, pathsRegresores=NULL, 
                              coordsAInterpolar, paramsIyM, shpMask, xyLims, listaMapas=NULL, returnInterpolacion=TRUE, 
                              paramsParaRellenoRegresores=NULL, pathsRegresoresParaRellenoRegresores=NULL, 
                              espEscalaFija=NULL, espEscalaAdaptada=NULL, tsAInterpolar=1:nrow(valoresObservaciones)) {
  # tsAInterpolar=1:nrow(valoresObservaciones)
  if (is.null(listaMapas)) listaMapas <- createDefaultListaMapas(paramsIyM = paramsIyM, fechasObservaciones = fechasObservaciones)
  
  # Preprocesamiento de los datos
  if (length(coordsObservaciones) <= 500 | !(paramsIyM$incorporarCoordenadas | paramsIyM$incorporarDistanciaAlAgua | paramsIyM$incorporarAltitud | paramsIyM$incorporarTiempo)) {
    # Elimino estaciones que no tengan ninguna observación para mejorar performance, sobre todo de la CV
    # Aplico esto solo si hay menos de 500 observaciones, sino asumo que es un relleno de regresores en donde no conviene
    # sacar los NA porque hay que volver a cachear los regresores estáticos
    iConDatos <- apply(valoresObservaciones, MARGIN = 2, FUN = function(x) { return( !all(is.na(x))) })
    coordsObservaciones <- coordsObservaciones[iConDatos,]
    valoresObservaciones <- valoresObservaciones[,iConDatos, drop=F]
    rm(iConDatos)
  }
  
  if (paramsIyM$difMaxFiltradoDeOutliersRLM > 0) {
    # params = paramsIyM
    # factorMADHaciaAbajo = paramsIyM$difMaxFiltradoDeOutliersRLM
    outliersRLM <- deteccionOutliersRLM(coordsObservaciones = coordsObservaciones, fechasObservaciones = fechasObservaciones, 
                                        valoresObservaciones = valoresObservaciones, params = paramsIyM, 
                                        pathsRegresores = pathsRegresores, listaMapas = listaMapas, 
                                        factorMADHaciaAbajo = paramsIyM$difMaxFiltradoDeOutliersRLM)
  } else { outliersRLM <- NULL }
  
  # Eliminación de outliers
  if (paramsIyM$difMaxFiltradoDeOutliersCV > 0) {
    # params = paramsIyM
    # maxOutlyingness = paramsIyM$difMaxFiltradoDeOutliersCV
    outliersCV <- deteccionOutliersUniversalGriddingCV(coordsObservaciones = coordsObservaciones, fechasObservaciones = fechasObservaciones, 
                                                       valoresObservaciones=valoresObservaciones, pathsRegresores = pathsRegresores, 
                                                       params=paramsIyM, maxOutlyingness=paramsIyM$difMaxFiltradoDeOutliersCV)
  } else { outliersCV <- NULL }
  
  if (!is.null(outliersRLM) | !is.null(outliersCV)) {
    outliers <- rbind(outliersRLM, outliersCV)
    outliers <- outliers[!duplicated(outliers$iOutlier), ]
    
    # Solo guardo los mapas de las fechas que tengan algún outlier interno al país
    if (!is.null(shpMask)) { iEstacionesInternas <- !is.na(over(x = geometry(coordsObservaciones), y = geometry(shpMask$shp)))
    } else { iEstacionesInternas <- rep(T, length(coordsObservaciones)) }
    
    outliersDeFechas <- list()
    fechasConOutliersInternos <- integer(0)
    i <- unique(outliers$iFecha)[1]
    for (i in unique(outliers$iFecha)) {
      iEstacionesOutlierDeFecha <- outliers$iEstacion[outliers$iFecha == i]
      if (any(iEstacionesInternas[iEstacionesOutlierDeFecha])) {
        fechasConOutliersInternos <- c(fechasConOutliersInternos, i)
        outliersDeFechas[[length(outliersDeFechas) + 1]] <- iEstacionesOutlierDeFecha
      }
    }
    
    if (length(fechasConOutliersInternos) > 0 & 
        any(listaMapas$dibujarEscalaFija[fechasConOutliersInternos] | listaMapas$dibujarEscalaAdaptada[fechasConOutliersInternos])) {
      listaMapasAux <- listaMapas
      listaMapasAux$nombreArchivo <- agregarCarpetaAlFinal(listaMapas$nombreArchivo, 'SinFiltrar')
      listaMapasAux$generarThumbnailFija <- F
      listaMapasAux$generarThumbnailAdaptada <- F
      listaMapasAux$salvarBin <- F
      listaMapasAux$salvarGeoTiff <- F
      
      paramsAux <- paramsIyM
      paramsAux$puntosAResaltar <- lapply(outliersDeFechas, FUN = function(x, coordsObservaciones) { return(coordsObservaciones[x,]) }, coordsObservaciones)
      paramsAux$difMaxFiltradoDeOutliersRLM <- 0
      paramsAux$difMaxFiltradoDeOutliersCV <- 0
      
      # paramsAux$nCoresAUsar <- 0
      # paramsIyM = paramsAux
      # listaMapas = listaMapasAux
      # tsAInterpolar = fechasConOutliersInternos
      # returnInterpolacion = F
      
      interpolarYMapear(coordsObservaciones = coordsObservaciones, fechasObservaciones = fechasObservaciones, 
                        valoresObservaciones = valoresObservaciones, pathsRegresores = pathsRegresores, 
                        coordsAInterpolar = coordsAInterpolar, paramsIyM = paramsAux, shpMask = shpMask, 
                        xyLims = xyLims, listaMapas = listaMapasAux,  returnInterpolacion = F, 
                        paramsParaRellenoRegresores = paramsParaRellenoRegresores, 
                        pathsRegresoresParaRellenoRegresores = pathsRegresoresParaRellenoRegresores, 
                        espEscalaFija = espEscalaFija, espEscalaAdaptada = espEscalaAdaptada, 
                        tsAInterpolar=fechasConOutliersInternos)
    }
    
    write.table(outliers, paste(paramsIyM$pathProceso, 'dfOutliers.txt', sep=''), sep = '\t', dec = '.', row.names = F)
    valoresObservaciones[outliers$iOutlier] <- NA
    rm(outliers)
  }
  
  # Reducción de series
  if (paramsIyM$radioReduccionSeriesKm > 0) {
    # radioReduccionSeriesKm = paramsIyM$radioReduccionSeriesKm
    # funcionReduccionSeries = paramsIyM$funcionReduccionSeries
    obsReducidas <- reducirSpatialPointsDataFrameYMatrizObservaciones(coordsObservaciones, valoresObservaciones, radioReduccionSeriesKm = paramsIyM$radioReduccionSeriesKm, 
                                                                      funcionReduccionSeries = paramsIyM$funcionReduccionSeries)
    coordsObservaciones <- obsReducidas$coordsObservaciones
    valoresObservaciones <- obsReducidas$valoresObservaciones
    rm(obsReducidas)
  }
  
  if (is.na(paramsIyM$cutoff)) paramsIyM$cutoff <- getDefaultSpatialCutoff(coordsObservaciones, params = paramsIyM)
  
  # Cargamos las observaciones de los regresores en las coordenadas de las estaciones
  if (!is.null(pathsRegresores)) {
    valoresRegresoresSobreObservaciones <- extraerValoresRegresoresSobreSP(
      objSP = coordsObservaciones, pathsRegresores = pathsRegresores, 
      nCoresAUsar = paramsIyM$nCoresAUsar)
  } else { valoresRegresoresSobreObservaciones <- NULL }
  
  if (paramsIyM$incorporarCoordenadas | paramsIyM$incorporarTiempo | paramsIyM$incorporarDistanciaAlAgua | paramsIyM$incorporarAltitud) {
    cachearRegresoresEstaticos(coordsObservaciones = geometry(coordsObservaciones), coordsAInterpolar = geometry(coordsAInterpolar), nCoresAUsar = paramsIyM$nCoresAUsar)
  }
  
  if (length(paramsIyM$tlagsAR) <= 0) {
    if (paramsIyM$nCoresAUsar <= 0) {
      nCoresDisponibles <- getAvailableCores(maxCoresPerGB = 1)
      if (length(tsAInterpolar) >= 1) { 
        nCoresAUsar <- min(nCoresDisponibles, length(tsAInterpolar))
        paramsIyM$nCoresAUsar <- ceiling(nCoresDisponibles / nCoresAUsar)
      } else { 
        nCoresAUsar <- 1
        paramsIyM$nCoresAUsar <- nCoresDisponibles
      }
    } else {
      nCoresAUsar <- paramsIyM$nCoresAUsar
      paramsIyM$nCoresAUsar <- 1
    }
    # paramsParaRellenoRegresores <- NULL
    # pathsRegresoresParaRellenoRegresores <- NULL
    if (!is.null(paramsParaRellenoRegresores)) paramsParaRellenoRegresores$nCoresAUsar <- paramsIyM$nCoresAUsar
    
    #returnInterpolacion <- T
    if (nCoresAUsar > 1) {
      cl <- makeCluster(getOption("cl.cores", nCoresAUsar))
      clusterExport(cl, varlist = c('script.dir.interpolarYMapearEx'))
      clusterEvalQ(cl, {
        set.seed(31)
        source(paste0(script.dir.interpolarYMapearEx, 'interpolarEx.r'), encoding = 'WINDOWS-1252')
        source(paste0(script.dir.interpolarYMapearEx, 'interpolarYMapearEx.r'), encoding = 'WINDOWS-1252')
      })
      if (exists(x = 'setMKLthreads')) { clusterEvalQ(cl = cl, expr = setMKLthreads(1)) }
      
      res <- parLapplyLB(
        cl=cl, X=seq(along.with = tsAInterpolar), fun=interpolarYMapearI, tsAInterpolar=tsAInterpolar,
        coordsObservaciones=coordsObservaciones, fechasObservaciones=fechasObservaciones, 
        valoresObservaciones=valoresObservaciones, pathsRegresores=pathsRegresores, 
        valoresRegresoresSobreObservaciones=valoresRegresoresSobreObservaciones, 
        coordsAInterpolar=coordsAInterpolar, paramsIyM=paramsIyM, shpMask=shpMask, xyLims=xyLims, 
        listaMapas=listaMapas, returnInterpolacion=returnInterpolacion, 
        paramsParaRellenoRegresores=paramsParaRellenoRegresores, 
        pathsRegresoresParaRellenoRegresores=pathsRegresoresParaRellenoRegresores, 
        espEscalaFija=espEscalaFija, espEscalaAdaptada=espEscalaAdaptada)
      stopCluster(cl)
    } else {
      set.seed(31)
      res <- lapply(
        X=seq(along.with = tsAInterpolar), FUN=interpolarYMapearI, tsAInterpolar=tsAInterpolar,
        coordsObservaciones=coordsObservaciones, fechasObservaciones=fechasObservaciones,
        valoresObservaciones=valoresObservaciones, pathsRegresores=pathsRegresores, 
        valoresRegresoresSobreObservaciones=valoresRegresoresSobreObservaciones, 
        coordsAInterpolar=coordsAInterpolar, paramsIyM=paramsIyM, shpMask=shpMask, xyLims=xyLims, 
        listaMapas=listaMapas, returnInterpolacion=returnInterpolacion, 
        paramsParaRellenoRegresores=paramsParaRellenoRegresores, 
        pathsRegresoresParaRellenoRegresores=pathsRegresoresParaRellenoRegresores, 
        espEscalaFija=espEscalaFija, espEscalaAdaptada=espEscalaAdaptada)
    }
  } else {
    nT <- max(tsAInterpolar)
    tsAInterpolar = 1:nT
    # Si hay componente autorregresivo no hay independecia temporal por lo que solo podemos aprovechar el paralelismo en forma espacial
    if (paramsIyM$nCoresAUsar <= 0) { paramsIyM$nCoresAUsar <- detectCores(T, T)
    } else { paramsIyM$nCoresAUsar <- 1 }
    
    tUltimoValorSinAR <- min(2 * paramsIyM$ventanaIgualacionDistribuciones + max(paramsIyM$tlagsAR) - 1, nT)
    nOrig <- ncol(pathsRegresores)
    
    vRegsSobreObservacionesTLags <- vector(mode = "list", length = length(paramsIyM$tlagsAR))
    for (i in 1:length(paramsIyM$tlagsAR))
      vRegsSobreObservacionesTLags[[i]] <- matrix(NA, nrow = nrow(valoresObservaciones), ncol=ncol(valoresObservaciones))
    names(vRegsSobreObservacionesTLags) <- paste('T_', paramsIyM$tlagsAR[i], sep='')

    if (returnInterpolacion) { res <- vector(mode = "list", length = length(valoresObservaciones)) }
    
    # Primeros valores sin AR para inicializar la ventana
    # ti <- 2
    for (iTi in 1:tUltimoValorSinAR) {
      interpolacionTi <- interpolarYMapearI(iTi = iTi, tsAInterpolar=tsAInterpolar, coordsObservaciones=coordsObservaciones, 
                                            fechasObservaciones=fechasObservaciones, valoresObservaciones=valoresObservaciones, 
                                            pathsRegresores=pathsRegresores, valoresRegresoresSobreObservaciones=valoresRegresoresSobreObservaciones, 
                                            coordsAInterpolar=coordsAInterpolar, paramsIyM=paramsIyM, shpMask=shpMask, 
                                            xyLims=xyLims, listaMapas=listaMapas, paramsParaRellenoRegresores=paramsParaRellenoRegresores, 
                                            pathsRegresoresParaRellenoRegresores=pathsRegresoresParaRellenoRegresores, 
                                            espEscalaFija=espEscalaFija, espEscalaAdaptada=espEscalaAdaptada)
      
      interpolacionTiSobreObservaciones <- as.numeric(t(over(x = coordsObservaciones, y = interpolacionTi$predictions)[,1]))
      for (i in 1:length(paramsIyM$tlagsAR)) {
        if (iTi + paramsIyM$tlagsAR[i] < nrow(valoresObservaciones)) {
          vRegsSobreObservacionesTLags[[i]][iTi + paramsIyM$tlagsAR[i], ] <- interpolacionTiSobreObservaciones
        }
      }
      
      if (returnInterpolacion) res[[iTi]] <- interpolacionTi
    }
    
    valoresRegresoresSobreObservaciones <- c(valoresRegresoresSobreObservaciones, vRegsSobreObservacionesTLags)
    rm(vRegsSobreObservacionesTLags)
    
    # Extiendo la matriz de paths de regresores para incluir los lags temporales 
    i <- 1
    for (i in 1:length(paramsIyM$tlagsAR)) {
      # Cargo los valores de los pasos de inicialización en el vector de paths de cada lag temporal
      pathsRegresores <- cbind(pathsRegresores, rep(NA_character_, nrow(pathsRegresores)))
      
      colnames(pathsRegresores)[nOrig + i] <- paste('T_', paramsIyM$tlagsAR[i], sep='')
      for (iTi in 1:nT) {
        if (iTi + paramsIyM$tlagsAR[i] <= nrow(pathsRegresores)) {
          pathsRegresores[iTi + paramsIyM$tlagsAR[i], nOrig + i] <- changeFileExt(listaMapas$nombreArchivo[iTi], '.tif')
        }
      }
    }

    # A partir de acá es con AR
    # iTi <- tUltimoValorSinAR + 1
    #iTi <- 632
    #for (iTi in (tUltimoValorSinAR+1):631) {
    for (iTi in (tUltimoValorSinAR+1):nT) {
      interpolacionTi <- interpolarYMapearI(iTi = iTi, tsAInterpolar = tsAInterpolar, coordsObservaciones=coordsObservaciones, fechasObservaciones=fechasObservaciones, valoresObservaciones=valoresObservaciones, 
                                            pathsRegresores=pathsRegresores, valoresRegresoresSobreObservaciones=valoresRegresoresSobreObservaciones, 
                                            coordsAInterpolar=coordsAInterpolar, paramsIyM=paramsIyM, shpMask=shpMask, xyLims=xyLims, listaMapas=listaMapas, 
                                            paramsParaRellenoRegresores=paramsParaRellenoRegresores, 
                                            pathsRegresoresParaRellenoRegresores=pathsRegresoresParaRellenoRegresores, 
                                            espEscalaFija=espEscalaFija, espEscalaAdaptada=espEscalaAdaptada)
      
      interpolacionTiSobreObservaciones <- as.numeric(t(over(x = coordsObservaciones, y = interpolacionTi$predictions)[1]))
      i<-1
      for (i in 1:length(paramsIyM$tlagsAR)) {
        if (iTi + paramsIyM$tlagsAR[i] <= nT) {
          valoresRegresoresSobreObservaciones[[nOrig + i]][iTi + paramsIyM$tlagsAR[i], ] <- interpolacionTiSobreObservaciones
        }
      }
      
      if (returnInterpolacion) res[[iTi]] <- interpolacionTi
    }
  }
  
  return(res)
}
