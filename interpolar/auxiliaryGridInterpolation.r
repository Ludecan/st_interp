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

# All methods in this file (excepting getObservedValues) receive spatial objects as inputs, 
# so observations should be SpatialPointsDataFrames, grids should be SpatialPixelDataFrames

# ---------------------------------------------- Data selection -----------------------------------------------------
getTimeWindow <- function(timeStep, observations, nDaysBack=5, nDaysForth=nDaysBack) {
  # observations is a 2D matrix so that observations[i, j] is the observed value at time i and location j
  # timeStep is the time step
  # if I'm close to the end for every two days forward I can't get, take one day from the past
  if (timeStep + nDaysForth > nrow(observations)) {
    nDaysBack <- nDaysBack + round((timeStep + nDaysForth - nrow(observations)) / 2)
    tMax <- nrow(observations)
  } else {
    tMax <- timeStep + nDaysForth
  }
  
  if (timeStep - nDaysBack < 1) {
    tMax <- min(nrow(observations), tMax + round((1 - (timeStep - nDaysBack)) / 2))
    tMin <- 1
  } else {
    tMin <- timeStep - nDaysBack
  }

  return (list(timeStep=timeStep, tMin=tMin, tMax=tMax))
}

# ----------------------------------------- Bias removal methods -------------------------------------------------
simpleBiasAdjustment <- function(observationT, spObservations, gridOverObservationT, spSatelliteGrid, satelliteGridT, interpolationParams) {
  residuals <- observationT - gridOverObservationT
  spObservations$value <- residuals

  residualInterpolation <- interpolarEx(spObservations, spSatelliteGrid, interpolationParams, longitudesEnColumnas=F)
  residualInterpolation$predictions@data[,residualInterpolation$campoMedia] <- residualInterpolation$predictions@data[,residualInterpolation$campoMedia] + as.vector(satelliteGridT)

  #image.plot(satGridLons$vals, satGridLats$vals, matrix(residualInterpolation$predictions@data[,residualInterpolation$campoMedia], nrow=nrow(satelliteGridT), ncol=ncol(satelliteGridT)))
  return (matrix(residualInterpolation$predictions@data[,residualInterpolation$campoMedia], nrow=nrow(satelliteGridT), ncol=ncol(satelliteGridT)))
}

preCDFBiasAdjustment <- function(observationsTs, gridOverObservationsTs, satelliteGridT) {
  # NAs get removed
  #sobs <- sort(as.vector(observationsTs))
  #sestim <- sort(as.vector(gridOverObservationsTs))
  #ajuste <- satelliteGridT
  
  sestim <- ecdf(gridOverObservationsTs)
  probs <- sestim(satelliteGridT)
  ajuste <- quantile(x=observationsTs, probs=probs, na.rm=T)
  ajuste <- matrix(ajuste, nrow=nrow(satelliteGridT), ncol=ncol(satelliteGridT))
  
  #image.plot(satGridLons$vals, satGridLats$vals, satelliteGridT)
  #image.plot(satGridLons$vals, satGridLats$vals, ajuste)
  #image.plot(satGridLons$vals, satGridLats$vals, satelliteGridT-ajuste)
  
  return(ajuste)
  #for (i in 1:length(satelliteGridT)) {
  #  if (satelliteGridT[i] <= sestim[1]) {
  #    ajuste[i] <- sobs[1]
  #  } else if (satelliteGridT[i] >= sestim[length(sestim)]) {
  #    ajuste[i] <- sobs[length(sobs)]
  #  } else {
  #    j <- which.max(sestim < satelliteGridT[i])
  #    ajuste[i] <- sobs[j] + (sobs[j+1]-sobs[j]) / (sestim[j+1]-sestim[j]) * (satelliteGridT[i]-satelliteGridT[j]) # Ajuste lineal
  #  }
  #}
  #return(ajuste)
  
  #return (CDFmatching(observationsTs, gridOverObservationsTs, satelliteGridT))
}

# ----------------------------------------------- Regression -------------------------------------------------------

glmRegressGrid <- function(observationT, satelliteGridT, gridIndexes, family=gaussian(link=identity)) {
  if (family$family == 'Gamma' & min(observationT, na.rm=T) <= 0) {
    data <- data.frame(observations=as.vector(observationT)+1E-8, gridOverObservations=as.vector(satelliteGridT[gridIndexes]))
  } else { 
    data <- data.frame(observations=as.vector(observationT), gridOverObservations=as.vector(satelliteGridT[gridIndexes])) 
  }
  model <- glm(formula=observations ~ gridOverObservations + 1, data=data, family=family)
  newData <- data.frame(gridOverObservations=as.vector(satelliteGridT))
  
  #image.plot(satGridLons$vals, satGridLats$vals, matrix(predict(object=model, newdata=newData), nrow=nrow(satelliteGridT), ncol=ncol(satelliteGridT)))
  return(matrix(predict.glm(object=model, newdata=newData), nrow=nrow(satelliteGridT), ncol=ncol(satelliteGridT)))
}

# ----------------------------------------------- RnR Masks --------------------------------------------------------

getRnRMaskStationOnly <- function(observationT, spObservations, spSatelliteGrid, interpolationParams, excedanceProbability=0.95) {
  spObservations$value <- observationT
  binaryObservations <- crearObservacionesBinarias(spObservations)
  binaryInterpolation <- interpolarEx(binaryObservations, spSatelliteGrid, interpolationParams, longitudesEnColumnas=F)
  # lazy way. should look into it from observation, create series of binaryInterpolation at
  # locations where observaciones$value == 0 and get value with excedanceProbability excedance
  # probability
  threshold <- 0.3
  RnRmask <- getPredictionMatrix(binaryInterpolation) >= threshold
  mode(RnRmask) <- "integer"

  return (RnRmask)
}

getQuantileRnRMaskStationOnly <- function(observationT, spObservations, spSatelliteGrid, interpolationParams, excedanceProbability=0.95) {
  spObservations$value <- observationT
  binaryObservations <- crearObservacionesBinarias(spObservations)
  binaryInterpolation <- interpolarEx(binaryObservations, spSatelliteGrid, interpolationParams, longitudesEnColumnas=T)
  # lazy way. should look into it from observation, create series of binaryInterpolation at
  # locations where observaciones$value == 0 and get value with excedanceProbability excedance
  # probability
  threshold <- 0.5
  RnRmask <- getPredictionMatrix(binaryInterpolation) >= threshold
  mode(RnRmask) <- "integer"
  
  return (RnRmask)
}

getRnRMaskSatelliteOnly <- function(satelliteGridT, threshold=0.1) {
  RnRmask <- satelliteGridT >= threshold
  mode(RnRmask) <- "integer"
  
  return (RnRmask)
}

getRnRMaskStationAndSatellite <- function(observationT, spObservations, satelliteGridT, spSatelliteGrid, interpolationParams, excedanceProbability=0.95, satelliteThreshold=0.1, maxDistanceToStation=0.75) {
  stationRnRMask <- getRnRMaskStationOnly(observationT, spObservations, spSatelliteGrid, interpolationParams, excedanceProbability)
  satelliteRnRMask <- getRnRMaskSatelliteOnly(satelliteGridT, satelliteThreshold)
  
  distMatrix <- rdist(spSatelliteGrid@coords, spObservations@coords)
  distToNearestObs <- apply(distMatrix, 1, function(x) min(x))
  
  ix <- which(distToNearestObs > maxDistanceToStation)
  stationRnRMask[ix] <- satelliteRnRMask[ix]
  
  return (stationRnRMask)
}