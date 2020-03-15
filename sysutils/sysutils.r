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

require('parallel')

getTotalRAM_GB <- function() {
  if (.Platform$OS.type == "windows") {
    memtot_gb <- memory.size(max=NA) / 1024
  } else {
    memtot_gb <- as.numeric(system("awk '/MemTot/ {print $2}' /proc/meminfo", intern=TRUE)) / 1024**2
  }
  return(memtot_gb)
}

getAvailableCores <- function(logicalCores=TRUE, maxCoresPerGB=Inf) {
  return(min(detectCores(T, logical=logicalCores), getTotalRAM_GB() * maxCoresPerGB))
}
