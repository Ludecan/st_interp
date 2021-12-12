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

require("methods")

getTokens <- function(params, separadorTokens=';') {
  paramsSplit <- unlist(strsplit(params, separadorTokens))
  tokens <- array(data=NA, dim=c(length(paramsSplit), 2))
  for (i in seq_along(paramsSplit)) {
    posIgual <- regexpr('=', paramsSplit[i], fixed=T)[1]
    tokens[i, 1] <- substr(x=paramsSplit[i], start=1, stop=posIgual - 1)
    tokens[i, 2] <- substr(x=paramsSplit[i], start=posIgual+1, stop=nchar(paramsSplit[i]))
  }
  
  return(tokens)
}

getParamValue <- function(tokens, paramName, classOfValue=NA, obligatorio=FALSE, paramDefaultValue=NA) {
  i <- match(paramName, tokens[,1])
  if (!is.na(i)) {
    if (tokens[i, 2] != 'NA') {
      if (!is.na(classOfValue)) { return (as(tokens[i, 2], classOfValue))
      } else { return (type.convert(tokens[i, 2], as.is=T)) }
    } else { return(NA) }
  } else if (!obligatorio) { 
    if (!is.na(classOfValue)) { return (as(paramDefaultValue, classOfValue))
    } else { return(paramDefaultValue)
    }
  } else { stop(paste('parsearParamsUtils.getParamValue: No se encuentra el valor de ', paramName, sep='')) }
}

getParamValues <- function(
    params, separadorTokens=';', paramNames, classOfValues=rep(NA_character_, length(paramNames)), 
    obligatorios=rep(FALSE, length(paramNames)), paramDefaultValues=rep(NA, length(paramNames))) {
  tokens <- getTokens(params, separadorTokens=separadorTokens)
  res <- vector(mode = "list", length = length(paramNames))
  names(res) <- paramNames
  
  i <- 1
  i <- i - 1
  for (i in seq(along=paramNames)) {
    res[[paramNames[i]]] <- getParamValue(
        tokens=tokens, paramName=paramNames[i], classOfValue=classOfValues[i], 
        obligatorio=obligatorios[i], paramDefaultValue=paramDefaultValues[[i]])
    i <- i + 1
  }
  return(res)
}

getParamValuesFromConstructorParams <- function(params, separadorTokens=';', funcCrearParams) {
  parametrosFunc <- formals(funcCrearParams)
  paramNames <- names(parametrosFunc)
  classOfValues <- rep(NA_character_, length(parametrosFunc))
  obligatorios <- rep(FALSE, length(parametrosFunc))
  paramDefaultValues <- list()
  length(paramDefaultValues) <- length(parametrosFunc)

  i <- 2
  for (i in seq_along(parametrosFunc)) {
    if (class(parametrosFunc[[i]]) == 'NULL') {
      obligatorios[i] <- FALSE
      paramDefaultValues[[i]] <- NULL
      classOfValues[i] <- 'NULL'
    } else if (class(parametrosFunc[[i]]) == 'name') {
      if (parametrosFunc[[i]] == 'F') {
        obligatorios[i] <- FALSE
        paramDefaultValues[[i]] <- FALSE
        classOfValues[i] <- 'logical'
      } else if (parametrosFunc[[i]] == 'T') {
        obligatorios[i] <- FALSE
        paramDefaultValues[[i]] <- TRUE
        classOfValues[i] <- 'logical'
      } else {
        obligatorios[i] <- T
      }
    } else if (class(parametrosFunc[[i]]) == 'call') {
      obligatorios[i] <- FALSE
      paramDefaultValues[[i]] <- eval(parametrosFunc[[i]])
      classOfValues[i] <- class(paramDefaultValues[[i]])
    } else if (is.null(parametrosFunc[[i]])) {
      obligatorios[i] <- FALSE
      paramDefaultValues[[i]] <- '<NULL>'
    } else if (is.na(parametrosFunc[[i]])) {
      obligatorios[i] <- FALSE
      paramDefaultValues[[i]] <- NA
    } else {
      obligatorios[i] <- FALSE
      if (!obligatorios[i]) {
        paramDefaultValues[[i]] <- parametrosFunc[[i]]
        classOfValues[i] <- class(x=parametrosFunc[[i]])
      }
    }
  }
  # cbind(paramNames, obligatorios, paramDefaultValues, classOfValues)
  return (getParamValues(params, separadorTokens, paramNames, classOfValues, obligatorios, paramDefaultValues))
}
