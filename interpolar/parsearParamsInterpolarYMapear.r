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
if (iFrame >= 3) { script.dir.parsearParamsInterpolarYMapear <- sys.frame(iFrame - 3)$ofile 
} else { script.dir.parsearParamsInterpolarYMapear <- NULL }
while ((is.null(script.dir.parsearParamsInterpolarYMapear) || is.na(regexpr('parsearParamsInterpolarYMapear.r', script.dir.parsearParamsInterpolarYMapear, fixed=T)[1])) && iFrame >= 0) {
  script.dir.parsearParamsInterpolarYMapear <- sys.frame(iFrame)$ofile
  iFrame <- iFrame - 1
}
if (is.null(script.dir.parsearParamsInterpolarYMapear)) { script.dir.parsearParamsInterpolarYMapear <- ''
} else { script.dir.parsearParamsInterpolarYMapear <- paste(dirname(script.dir.parsearParamsInterpolarYMapear), '/', sep='') }
source(paste(script.dir.parsearParamsInterpolarYMapear, '../parsearParams/parsearParamsUtils.r', sep=''))
# require('jsonlite')

createParamsInterpolarYMapear <- function(baseNomArchResultados='',
                                          pathEjecucion='', 
                                          pathProceso='',
                                          proj4StringObservaciones='+proj=longlat +datum=WGS84',
                                          proj4StringAInterpolar='+proj=utm +zone=21 +south +ellps=WGS84 +datum=WGS84 +towgs84=0,0,0 +units=m +no_defs',
                                          coordsAInterpolarSonGrilla=TRUE, 
                                          interpolationMethod='automap', # 'none'(sin interpolaci�n), 'idw' (inverse distance weighting), 'automap' (kriging), 'copula' o 'stUniversalKriging' (Kriging Espacio - Temporal)
                                          mLimitarValoresInterpolados='UsarPromDesvEst', #'NoLimitar', 'LimitarMinimo', 'LimitarMaximo', 'LimitarMinimoyMaximo', 'UsarPromDesvEst' o 'UsarPromDesvEstYMinimoYMaximo'
                                          minimoLVI=NA, maximoLVI=NA,
                                          factorDesvEstLVI=3.5,
                                          metodoIgualacionDistribuciones='regresionLineal', # 'ninguna', 'regresionLineal', 'regresionLinealRobusta', 'regresionLinealConEliminacionDeOutliers', 'CDFMatching' o 'GLS'
                                          formulaRegresion='',
                                          ventanaIgualacionDistribuciones=1,
                                          incorporarCoordenadas=FALSE, 
                                          formulaCoordenadas='x + y',
                                          incorporarTiempo=FALSE, 
                                          formulaTiempo='t',
                                          incorporarDistanciaAlAgua=FALSE, 
                                          formulaDistanciaAlAgua='I(dist^0.125)',
                                          incorporarAltitud=FALSE, 
                                          formulaAltitud='alt',
                                          descartarCoordenadasNoSignificativas=FALSE,
                                          rellenarRegresores=FALSE,
                                          invertirAjusteRegresores=FALSE,
                                          usarNugget=FALSE,
                                          block=NA,
                                          nmin=0,
                                          nmax=Inf,
                                          maxdist=Inf,
                                          inverseDistancePower=NA,
                                          umbralMascaraCeros=0,
                                          metodoRemocionDeSesgo='ninguno', # 'ninguno', 'IDW_ResiduosNegativos', 'IDW_ResiduosPositivos', 'IDW_ResiduosNegativosYPositivos'
                                          modelosVariograma=c('Exp', 'Sph', 'Pow', 'Cir', 'Pen'),
                                          cutoff=Inf,
                                          tlags=0:5, 
                                          nTsST=5,#max(tlags), 
                                          tlagsAR=NULL,
                                          tryFixNugget=FALSE,
                                          nPuntosIniciales=2,
                                          usarFitVariogramGLS='auto',
                                          modelosVariogramaST=c('Separable', 'ProductSum', 'Metric', 'SimpleSumMetric', 'SumMetric'),
                                          fit.methodST=6, 
                                          verbose=FALSE,
                                          pathSHPMapaBase='',
                                          nCoresAUsar=0,
                                          radioReduccionSeriesKm=0,
                                          funcionReduccionSeries='mean',
                                          difMaxFiltradoDeOutliersRLM=0,
                                          difMaxFiltradoDeOutliersCV=0,
                                          modoDiagnostico=FALSE,
                                          imitarSurfer=FALSE) {
  res <- list(pathEjecucion=pathEjecucion, pathProceso=pathProceso, baseNomArchResultados=baseNomArchResultados, proj4StringObservaciones=proj4StringObservaciones, 
              proj4StringAInterpolar=proj4StringAInterpolar, coordsAInterpolarSonGrilla=coordsAInterpolarSonGrilla, interpolationMethod=interpolationMethod, 
              mLimitarValoresInterpolados=mLimitarValoresInterpolados,minimoLVI=minimoLVI, maximoLVI=maximoLVI, factorDesvEstLVI=factorDesvEstLVI, 
              metodoIgualacionDistribuciones=metodoIgualacionDistribuciones, formulaRegresion=formulaRegresion, 
              ventanaIgualacionDistribuciones=ventanaIgualacionDistribuciones, incorporarCoordenadas=incorporarCoordenadas,
              formulaCoordenadas=formulaCoordenadas, incorporarTiempo=incorporarTiempo, formulaTiempo=formulaTiempo, 
              incorporarDistanciaAlAgua=incorporarDistanciaAlAgua, formulaDistanciaAlAgua=formulaDistanciaAlAgua, 
              incorporarAltitud=incorporarAltitud, formulaAltitud=formulaAltitud, descartarCoordenadasNoSignificativas=descartarCoordenadasNoSignificativas, 
              rellenarRegresores=rellenarRegresores, invertirAjusteRegresores=invertirAjusteRegresores, usarNugget=usarNugget, block=block, nmin=nmin, nmax=nmax, maxdist=maxdist,
              inverseDistancePower=inverseDistancePower,umbralMascaraCeros=umbralMascaraCeros, metodoRemocionDeSesgo=metodoRemocionDeSesgo, 
              modelosVariograma=modelosVariograma, cutoff=cutoff, tlags=tlags, nTsST=nTsST, tlagsAR=tlagsAR, tryFixNugget=tryFixNugget, 
              nPuntosIniciales=nPuntosIniciales, usarFitVariogramGLS=usarFitVariogramGLS, modelosVariogramaST= modelosVariogramaST, fit.methodST=fit.methodST, verbose=verbose,
              pathSHPMapaBase=pathSHPMapaBase, nCoresAUsar=nCoresAUsar,radioReduccionSeriesKm=radioReduccionSeriesKm,
              funcionReduccionSeries=funcionReduccionSeries, difMaxFiltradoDeOutliersRLM=difMaxFiltradoDeOutliersRLM,
              difMaxFiltradoDeOutliersCV=difMaxFiltradoDeOutliersCV, modoDiagnostico=modoDiagnostico, 
              imitarSurfer=imitarSurfer)
  if (res$usarNugget) { res$fixNugget <- NA
  } else { res$fixNugget <- 0 }
  return(res)
}

parsearParamsInterpolarYMapear <- function(params) {
  return(getParamValuesFromConstructorParams(params, funcCrearParams=createParamsInterpolarYMapear))
}

createParamsInterpolarYMapearPInterpolarFaltantesMultiObservaciones <- function(paramsMultiObservaciones) {
  return(createParamsInterpolarYMapear(
    baseNomArchResultados=paramsMultiObservaciones$nomArchResultados,
    pathEjecucion=paramsMultiObservaciones$pathEjecucion, 
    pathProceso=paramsMultiObservaciones$pathProceso,
    proj4StringObservaciones=paramsMultiObservaciones$proj4StringObservaciones,
    proj4StringAInterpolar=paramsMultiObservaciones$proj4StringAInterpolar,
    coordsAInterpolarSonGrilla=FALSE, 
    interpolationMethod=paramsMultiObservaciones$interpolationMethod,
    mLimitarValoresInterpolados=paramsMultiObservaciones$mLimitarValoresInterpolados,
    minimoLVI=paramsMultiObservaciones$minimoLVI, maximoLVI=paramsMultiObservaciones$maximoLVI,
    factorDesvEstLVI=paramsMultiObservaciones$factorDesvEstLVI,
    metodoIgualacionDistribuciones=paramsMultiObservaciones$metodoIgualacionDistribuciones,
    formulaRegresion='',
    ventanaIgualacionDistribuciones=paramsMultiObservaciones$ventanaIgualacionDistribuciones,
    incorporarCoordenadas=paramsMultiObservaciones$incorporarCoordenadas,
    formulaCoordenadas = paramsMultiObservaciones$formulaCoordenadas,
    incorporarTiempo = paramsMultiObservaciones$incorporarTiempo,
    formulaTiempo = paramsMultiObservaciones$formulaTiempo,
    incorporarDistanciaAlAgua=paramsMultiObservaciones$incorporarDistanciaAlAgua,
    formulaDistanciaAlAgua = paramsMultiObservaciones$formulaDistanciaAlAgua,
    incorporarAltitud=paramsMultiObservaciones$incorporarAltitud,
    formulaAltitud = paramsMultiObservaciones$formulaAltitud,
    rellenarRegresores=FALSE,
    invertirAjusteRegresores=FALSE,
    usarNugget=paramsMultiObservaciones$usarNugget,
    block=paramsMultiObservaciones$block,
    nmin=paramsMultiObservaciones$nmin,
    nmax=paramsMultiObservaciones$nmax,
    maxdist=paramsMultiObservaciones$maxdist,
    inverseDistancePower=paramsMultiObservaciones$inverseDistancePower,
    umbralMascaraCeros=paramsMultiObservaciones$umbralMascaraCeros,
    metodoRemocionDeSesgo = paramsMultiObservaciones$metodoRemocionDeSesgo,
    modelosVariograma=paramsMultiObservaciones$modelosVariograma,
    cutoff=paramsMultiObservaciones$cutoff,
    tlags=0:5, 
    nTsST=5,
    tlagsAR=NULL,
    tryFixNugget=paramsMultiObservaciones$tryFixNugget,
    nPuntosIniciales=paramsMultiObservaciones$nPuntosIniciales,
    usarFitVariogramGLS=paramsMultiObservaciones$usarFitVariogramGLS,
    modelosVariogramaST=paramsMultiObservaciones$modelosVariogramaST,
    fit.methodST=paramsMultiObservaciones$fit.methodST, 
    verbose=paramsMultiObservaciones$verbose,
    pathSHPMapaBase='',
    nCoresAUsar=0,
    radioReduccionSeriesKm=paramsMultiObservaciones$radioReduccionSeriesKm,
    funcionReduccionSeries=paramsMultiObservaciones$funcionReduccionSeries,
    difMaxFiltradoDeOutliersRLM=paramsMultiObservaciones$difMaxFiltradoDeOutliersRLM,
    difMaxFiltradoDeOutliersCV=paramsMultiObservaciones$difMaxFiltradoDeOutliersCV,
    imitarSurfer=paramsMultiObservaciones$imitarSurfer))
}
