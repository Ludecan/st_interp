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
if (iFrame >= 3) { script.dir.parsearParamsInterpolar <- sys.frame(iFrame - 3)$ofile 
} else { script.dir.parsearParamsInterpolar <- NULL }
while ((is.null(script.dir.parsearParamsInterpolar) || is.na(regexpr('parsearParamsInterpolar.r', script.dir.parsearParamsInterpolar, fixed=T)[1])) && iFrame >= 0) {
  script.dir.parsearParamsInterpolar <- sys.frame(iFrame)$ofile
  iFrame <- iFrame - 1
}
if (is.null(script.dir.parsearParamsInterpolar)) { script.dir.parsearParamsInterpolar <- ''
} else { script.dir.parsearParamsInterpolar <- paste(dirname(script.dir.parsearParamsInterpolar), '/', sep='') }

source(paste(script.dir.parsearParamsInterpolar, '../parsearParams/parsearParamsUtils.r', sep=''))

# ----------------------------------------- Interpolar ----------------------------------------------------------
createParamsInterpolar <- function(pathEjecucion='./',
                                   pathProceso='',
                                   nomArchResultados, 
                                   coordsAInterpolarSonGrilla,
                                   proj4StringObservaciones,
                                   proj4StringAInterpolar,
                                   interpolationMethod='automap', # one of 'none', 'idw' (inverse distance weighting), 'automap' (kriging) or 'copula'
                                   metodoIgualacionDistribuciones='ninguna', # 'ninguna', 'regresionLineal', 'regresionLinealRobusta', 'regresionLinealConEliminacionDeOutliers', 'CDFMatching' o 'GLS'
                                   ventanaIgualacionDistribuciones=5,
                                   mLimitarValoresInterpolados='NoLimitar', #one of 'NoLimitar', 'LimitarMinimo', 'LimitarMaximo', 'LimitarMinimoyMaximo', 'UsarPromDesvEst' or 'UsarPromDesvEstYMinimoYMaximo'
                                   minimoLVI=NA, maximoLVI=NA, factorDesvEstLVI=3.5,
                                   block=NA,
                                   salvarPrediccion=FALSE, 
                                   salvarVarianza=FALSE,
                                   formatoSalida='binary', # one of 'binary' or 'netCDF'
                                   graficarObservaciones=TRUE,
                                   graficarPrediccion=TRUE,
                                   graficarVarianza=FALSE,
                                   graficarSD=TRUE,
                                   usarNugget=FALSE,
                                   pathSHPMapaBase='',
                                   resXImagenes=640,
                                   resYImagenes=480,
                                   dibujarEscala=TRUE,
                                   numDecimales=1,
                                   debugeando=FALSE,
                                   nmin=0,
                                   nmax=Inf,
                                   maxdist=Inf,
                                   inverseDistancePower=NA,
                                   nCoresAUsar=0,
                                   incorporarCoordenadas=FALSE, 
                                   formulaCoordenadas='x + y',
                                   incorporarTiempo=FALSE, 
                                   formulaTiempo='t',
                                   incorporarDistanciaAlAgua=FALSE, 
                                   formulaDistanciaAlAgua='I(dist^0.125)',
                                   incorporarAltitud=FALSE, 
                                   formulaAltitud='alt',
                                   descartarCoordenadasNoSignificativas=FALSE,
                                   umbralMascaraCeros=0,
                                   metodoRemocionDeSesgo='ninguno', # 'ninguno', 'IDW_ResiduosNegativos', 'IDW_ResiduosPositivos', 'IDW_ResiduosNegativosYPositivos'
                                   modelosVariograma=c('Exp', 'Sph', 'Pow', 'Cir', 'Pen'),
                                   cutoff=NA,
                                   tryFixNugget=FALSE, 
                                   nPuntosIniciales=2,
                                   usarFitVariogramGLS='auto',
                                   modelosVariogramaST=c('Separable', 'ProductSum', 'Metric', 'SimpleSumMetric', 'SumMetric'),
                                   fit.methodST=6,
                                   verbose=FALSE,
                                   radioReduccionSeriesKm=0,
                                   funcionReduccionSeries='mean',
                                   modoDiagnostico=FALSE,
                                   imitarSurfer=FALSE) {
  # funcion auxiliar para saber todos los parámetros que se tienen que pasar en params y valores
  # por defecto de algunos parámetros
  return(list(pathEjecucion=pathEjecucion,
              pathProceso=pathProceso,
              nomArchResultados=nomArchResultados, 
              coordsAInterpolarSonGrilla=coordsAInterpolarSonGrilla,
              proj4StringObservaciones=proj4StringObservaciones,
              proj4StringAInterpolar=proj4StringAInterpolar,
              interpolationMethod=interpolationMethod,
              metodoIgualacionDistribuciones=metodoIgualacionDistribuciones,
              ventanaIgualacionDistribuciones=ventanaIgualacionDistribuciones,
              mLimitarValoresInterpolados=mLimitarValoresInterpolados, 
              minimoLVI=minimoLVI, maximoLVI=maximoLVI, factorDesvEstLVI=factorDesvEstLVI,              
              block=block,
              salvarPrediccion=salvarPrediccion, 
              salvarVarianza=salvarVarianza, 
              formatoSalida=formatoSalida,
              graficarObservaciones=graficarObservaciones, 
              graficarPrediccion=graficarPrediccion, 
              graficarVarianza=graficarVarianza,
              graficarSD=graficarSD,
              usarNugget=usarNugget,
              pathSHPMapaBase=pathSHPMapaBase,
              resXImagenes=resXImagenes,
              resYImagenes=resYImagenes,
              dibujarEscala=dibujarEscala,
              numDecimales=numDecimales,
              debugeando=debugeando,
              nmin=nmin,
              nmax=nmax,
              maxdist=maxdist,
              inverseDistancePower=inverseDistancePower,
              nCoresAUsar=nCoresAUsar,
              incorporarCoordenadas=incorporarCoordenadas,
              formulaCoordenadas=formulaCoordenadas,
              incorporarTiempo=incorporarTiempo,
              incorporarDistanciaAlAgua=incorporarDistanciaAlAgua,
              formulaDistanciaAlAgua=formulaDistanciaAlAgua,
              incorporarAltitud=incorporarAltitud,
              formulaAltitud=formulaAltitud,
              descartarCoordenadasNoSignificativas=descartarCoordenadasNoSignificativas,
              umbralMascaraCeros=umbralMascaraCeros,
              metodoRemocionDeSesgo=metodoRemocionDeSesgo,
              modelosVariograma=modelosVariograma,
              cutoff=cutoff, tryFixNugget=tryFixNugget, 
              nPuntosIniciales=nPuntosIniciales, 
              usarFitVariogramGLS=usarFitVariogramGLS,
              modelosVariogramaST=modelosVariogramaST,
              fit.methodST=fit.methodST,
              verbose=verbose,
              radioReduccionSeriesKm=radioReduccionSeriesKm,
              funcionReduccionSeries=funcionReduccionSeries,
              modoDiagnostico=modoDiagnostico,
              imitarSurfer=imitarSurfer))
}

parsearParamsInterpolar <- function(params) {
  return(getParamValuesFromConstructorParams(params, funcCrearParams=createParamsInterpolar))
}

# --------------------------------- InterpolarFaltantesMultiObservaciones ---------------------------------------
createParamsInterpolarFaltantesMultiObservaciones <- function(pathEjecucion,
                                                              pathProceso='',
                                                              incluyeNombresObservaciones=FALSE,
                                                              nomArchResultados, 
                                                              proj4StringObservaciones,
                                                              proj4StringAInterpolar,
                                                              interpolationMethod = 'automap', # one of 'idw', 'automap' or 'copula'
                                                              metodoIgualacionDistribuciones='ninguna', # 'ninguna', 'regresionLineal', 'regresionLinealRobusta', 'regresionLinealConEliminacionDeOutliers' o 'CDFMatching'
                                                              ventanaIgualacionDistribuciones=5,
                                                              mLimitarValoresInterpolados='NoLimitar', #one of 'NoLimitar', 'LimitarMinimo', 'LimitarMaximo', 'LimitarMinimoyMaximo', 'UsarPromDesvEst' or 'UsarPromDesvEstYMinimoYMaximo'
                                                              minimoLVI=NA, maximoLVI=NA, factorDesvEstLVI=3.5,
                                                              block=NA,
                                                              formatoSalida='binary', # one of 'binary' or 'netCDF'
                                                              usarNugget=FALSE,
                                                              nmin=0,
                                                              nmax=Inf,
                                                              maxdist=Inf,
                                                              inverseDistancePower=NA,
                                                              incorporarCoordenadas=FALSE, 
                                                              formulaCoordenadas='x + y',
                                                              incorporarTiempo=FALSE, 
                                                              formulaTiempo='t',
                                                              incorporarDistanciaAlAgua=FALSE, 
                                                              formulaDistanciaAlAgua='I(dist^0.125)',
                                                              incorporarAltitud=FALSE, 
                                                              formulaAltitud='alt',
                                                              umbralMascaraCeros=0,
                                                              metodoRemocionDeSesgo='ninguno', # 'ninguno', 'IDW_ResiduosNegativos', 'IDW_ResiduosPositivos', 'IDW_ResiduosNegativosYPositivos'
                                                              modelosVariograma=c('Exp', 'Sph', 'Pow', 'Cir', 'Pen'),
                                                              cutoff=NA,
                                                              tryFixNugget=FALSE, 
                                                              nPuntosIniciales=2,
                                                              usarFitVariogramGLS=FALSE,
                                                              modelosVariogramaST=c('Separable', 'ProductSum', 'Metric', 'SimpleSumMetric', 'SumMetric'),
                                                              fit.methodST=6,
                                                              verbose=FALSE,
                                                              radioReduccionSeriesKm=0,
                                                              funcionReduccionSeries='mean',
                                                              difMaxFiltradoDeOutliersRLM=0,
                                                              difMaxFiltradoDeOutliersCV=0,
                                                              modoDiagnostico=FALSE,
                                                              imitarSurfer=FALSE) {
  return(list(pathEjecucion=pathEjecucion,
              pathProceso=pathProceso,
              incluyeNombresObservaciones=incluyeNombresObservaciones,
              nomArchResultados=nomArchResultados, 
              proj4StringObservaciones=proj4StringObservaciones,
              proj4StringAInterpolar=proj4StringAInterpolar,
              interpolationMethod=interpolationMethod,
              metodoIgualacionDistribuciones=metodoIgualacionDistribuciones,
              ventanaIgualacionDistribuciones=ventanaIgualacionDistribuciones,
              mLimitarValoresInterpolados=mLimitarValoresInterpolados,
              minimoLVI=minimoLVI, maximoLVI=maximoLVI, factorDesvEstLVI=factorDesvEstLVI,              
              block=block,
              formatoSalida=formatoSalida,
              usarNugget=usarNugget,
              nmin=nmin,
              nmax=nmax,
              maxdist=maxdist,
              inverseDistancePower=inverseDistancePower,
              incorporarCoordenadas=incorporarCoordenadas,
              formulaCoordenadas=formulaCoordenadas,
              incorporarTiempo=incorporarTiempo,
              formulaTiempo=formulaTiempo,
              incorporarDistanciaAlAgua=incorporarDistanciaAlAgua,
              formulaDistanciaAlAgua=formulaDistanciaAlAgua,
              incorporarAltitud=incorporarAltitud,
              formulaAltitud=formulaAltitud,
              umbralMascaraCeros=umbralMascaraCeros,
              metodoRemocionDeSesgo=metodoRemocionDeSesgo,
              modelosVariograma=modelosVariograma,
              cutoff=cutoff,
              tryFixNugget=tryFixNugget, 
              nPuntosIniciales=nPuntosIniciales,
              usarFitVariogramGLS=usarFitVariogramGLS,
              modelosVariogramaST=modelosVariogramaST,
              fit.methodST=fit.methodST,
              verbose=verbose,
              radioReduccionSeriesKm=radioReduccionSeriesKm,
              funcionReduccionSeries=funcionReduccionSeries,
              difMaxFiltradoDeOutliersRLM=difMaxFiltradoDeOutliersRLM,
              difMaxFiltradoDeOutliersCV=difMaxFiltradoDeOutliersCV,
              modoDiagnostico=modoDiagnostico,
              imitarSurfer=imitarSurfer))
}

parsearParamsInterpolarFaltantesMultiObservaciones <- function(params) {
  return(getParamValuesFromConstructorParams(params, funcCrearParams=createParamsInterpolarFaltantesMultiObservaciones))
}

# ------------------------------------------------- Mapear -----------------------------------------------------
createParamsMapear <- function(pathEjecucion,
                               pathProceso='',
                               incluyeNombresPuntos=FALSE,
                               coordsAMapearSonGrilla,
                               proj4String,
                               interpolationMethod = 'automap', # one of 'idw', 'automap' or 'copula'
                               pathSHPMapaBase='',
                               resXImagenes=640,
                               resYImagenes=480,
                               resXThumbnails=64,
                               resYThumbnails=64,
                               dibujarEscala=TRUE,
                               numDecimales=1) {
  # funcion auxiliar para saber todos los parámetros que se tienen que pasar en params y valores
  # por defecto de algunos parámetros
  return(list(pathEjecucion=pathEjecucion,
              pathProceso=pathProceso,
              incluyeNombresPuntos=incluyeNombresPuntos,
              coordsAMapearSonGrilla=coordsAMapearSonGrilla,
              proj4String=proj4String,
              pathSHPMapaBase=pathSHPMapaBase,
              resXImagenes=resXImagenes,
              resYImagenes=resYImagenes,
              resXThumbnails=resXThumbnails,
              resYThumbnails=resYThumbnails,
              dibujarEscala=dibujarEscala,
              numDecimales=numDecimales))
}

parsearParamsMapear <- function(params) {
  return(getParamValuesFromConstructorParams(params, funcCrearParams=createParamsMapear))
}