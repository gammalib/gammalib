/***************************************************************************
 *                          latlib - SWIG file                             *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010 by Jurgen Knodlseder                                *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 * ----------------------------------------------------------------------- *
 * Usage:                                                                  *
 * swig -c++ -python -Wall -includeall latlib.i                            *
 ***************************************************************************/
%module latlib

/* __ LAT ________________________________________________________________ */
#include "GLATObservation.i"
#include "GLATEventCube.i"
#include "GLATEventBin.i"
#include "GLATResponse.i"
#include "GLATInstDir.i"
#include "GLATRoi.i"
#include "GLATLtCube.i"



