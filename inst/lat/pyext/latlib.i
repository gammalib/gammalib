/***************************************************************************
 *                          latlib - SWIG file                             *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2011 by Jurgen Knodlseder                           *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *  This program is free software: you can redistribute it and/or modify   *
 *  it under the terms of the GNU General Public License as published by   *
 *  the Free Software Foundation, either version 3 of the License, or      *
 *  (at your option) any later version.                                    *
 *                                                                         *
 *  This program is distributed in the hope that it will be useful,        *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of         *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          *
 *  GNU General Public License for more details.                           *
 *                                                                         *
 *  You should have received a copy of the GNU General Public License      *
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.  *
 *                                                                         *
 * ----------------------------------------------------------------------- *
 * Usage:                                                                  *
 * swig -c++ -python -Wall latlib.i                                        *
 ***************************************************************************/
%module latlib

/* __ LAT ________________________________________________________________ */
%include "GLATObservation.i"
%include "GLATEventCube.i"
%include "GLATEventBin.i"
%include "GLATResponse.i"
%include "GLATInstDir.i"
%include "GLATRoi.i"
%include "GLATLtCube.i"
%include "GLATMeanPsf.i"



