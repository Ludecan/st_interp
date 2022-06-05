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
if (iFrame >= 3) { script.dir.cacheFunciones <- sys.frame(iFrame - 3)$ofile 
} else { script.dir.cacheFunciones <- NULL }
while ((is.null(script.dir.cacheFunciones) || is.na(regexpr('mapearEx.r', script.dir.cacheFunciones, fixed=T)[1])) && iFrame >= 0) {
  script.dir.cacheFunciones <- sys.frame(iFrame)$ofile
  iFrame <- iFrame - 1
}
if (is.null(script.dir.cacheFunciones)) { script.dir.cacheFunciones <- ''
} else { script.dir.cacheFunciones <- paste0(dirname(script.dir.cacheFunciones), '/') }

source(paste0(script.dir.cacheFunciones, '../instalarPaquetes/instant_pkgs.r'))
instant_pkgs(pkgs=c('digest'))

getPathCache <- function(objParametros, dirEjecucion='', prefijoNombreArchivoCache='') {
  return (paste0(dirEjecucion, 'RCache_', .Platform$OS.type, '/', prefijoNombreArchivoCache, 
                 digest(objParametros), '.rds'))
}

guardarCache <- function(pathCache, obj) {
  dir.create(dirname(pathCache), recursive=T, showWarnings=F)
  saveRDS(object=obj, file=pathCache)
}

cargarCache <- function(pathCache) {
  return(readRDS(file = pathCache))
}
