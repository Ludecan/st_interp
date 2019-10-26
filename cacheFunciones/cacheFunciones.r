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

script.dir.cacheFunciones <- dirname((function() { attr(body(sys.function()), "srcfile") })()$filename)
source(paste(script.dir.cacheFunciones, '/../instalarPaquetes/instant_pkgs.r', sep=''))
instant_pkgs(c('digest'))

getPathCache <- function(objParametros, dirEjecucion='', prefijoNombreArchivoCache='') {
  return (paste(dirEjecucion, 'RCache/', prefijoNombreArchivoCache, digest(objParametros), '.rds', sep=''))
}

guardarCache <- function(pathCache, obj) {
  dir.create(dirname(pathCache), recursive=T, showWarnings=F)
  saveRDS(object = obj, file=pathCache)
}

cargarCache <- function(pathCache) {
  return(readRDS(file = pathCache))
}
