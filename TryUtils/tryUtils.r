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

evaluarConReintentos <- function(expr, maxNIntentos=5, segundosEntreIntentos=5) {
  hecho <- FALSE
  nIntentos <- 0
  while (!hecho & nIntentos < maxNIntentos) {
    result = tryCatch({
      eval.parent(expr)
    }, error = function(e) {
      print(paste0("MY_ERROR:  ", e))
      result <- list(error=e)
      class(result) <- 'try-error'
      return(result)
    })

    nIntentos <- nIntentos + 1    
    if (class(result) == 'try-error') {
      if (nIntentos < maxNIntentos) { Sys.sleep(segundosEntreIntentos) }
    } else {
      hecho <- TRUE
    }
  }
  return(hecho)
}

tryExpr <- function(expr, silent=F) {
  result = tryCatch({
    eval.parent(expr)
  }, error = function(e) {
    if (!silent) print(e)
    class(e) <- c(class(e), 'try-error')
    return(e)
  })
  
  return (!'try-error' %in% class(result))
}