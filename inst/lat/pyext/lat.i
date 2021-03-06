/***************************************************************************
 *                       lat.i - Fermi/LAT module                          *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2011-2017 by Juergen Knoedlseder                         *
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
 * swig -c++ -python -Wall lat.i                                           *
 ***************************************************************************/
/**
 * @file lat.i
 * @brief Fermi/LAT module 
 * @author Juergen Knoedlseder
 */
%module lat
%feature("autodoc", "1");

/* __ Headers needed for compilation _____________________________________ */
%{
#include <stddef.h>
#include "GException.hpp"
#include "GTools.hpp"
%}

/* __ Include standard typemaps for vectors and strings __________________ */
%include stl.i

/* __ Include GammaLib typemaps __________________________________________ */
%include typemap_GFilename.i
%include typemap_slices.i

/* __ Include interface classes __________________________________________ */
%import(module="gammalib.base") "GBase.i";

/* __ Make sure that exceptions are catched ______________________________ */
%import(module="gammalib.support") "GException.i";

/* __ Inform about base classes __________________________________________ */
%import(module="gammalib.obs") "GObservation.i";
%import(module="gammalib.obs") "GEvent.i";
%import(module="gammalib.obs") "GEventAtom.i";
%import(module="gammalib.obs") "GEventBin.i";
%import(module="gammalib.obs") "GEvents.i";
%import(module="gammalib.obs") "GEventList.i";
%import(module="gammalib.obs") "GEventCube.i";
%import(module="gammalib.obs") "GResponse.i";
%import(module="gammalib.obs") "GInstDir.i";
%import(module="gammalib.obs") "GRoi.i";

/* __ LAT ________________________________________________________________ */
%include "GLATAeff.i"
%include "GLATEdisp.i"
%include "GLATEventAtom.i"
%include "GLATEventBin.i"
%include "GLATEventCube.i"
%include "GLATEventList.i"
%include "GLATInstDir.i"
%include "GLATLtCube.i"
%include "GLATMeanPsf.i"
%include "GLATObservation.i"
%include "GLATPsf.i"
%include "GLATResponse.i"
%include "GLATRoi.i"



