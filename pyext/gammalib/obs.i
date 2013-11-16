/***************************************************************************
 *                       obs.i - Observation module                        *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2011-2013 by Juergen Knoedlseder                         *
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
 * swig -c++ -python -Wall obs.i                                           *
 ***************************************************************************/
/**
 * @file obs.i
 * @brief Observation module
 * @author Juergen Knoedlseder
 */
%module obs
%feature("autodoc", "1");

/* __ Headers needed for compilation _____________________________________ */
%{
#include <stddef.h>
#include "GException.hpp"
#include "GTools.hpp"
%}

/* __ Include standard typemaps for vectors and strings __________________ */
%include stl.i

/* __ Include interface classes __________________________________________ */
%import(module="gammalib.base") "GBase.i";
%import(module="gammalib.base") "GContainer.i";
%import(module="gammalib.base") "GRegistry.i";

/* __ Make sure that exceptions are catched ______________________________ */
%import(module="gammalib.support") "GException.i";

/* __ Observation handling _______________________________________________ */
%include "GEnergy.i"
%include "GEnergies.i"
%include "GTime.i"
%include "GTimes.i"
%include "GTimeReference.i"
%include "GEbounds.i"
%include "GGti.i"
%include "GCaldb.i"
%include "GObservations.i"
%include "GObservation.i"
%include "GObservationRegistry.i"
%include "GEvents.i"
%include "GEventList.i"
%include "GEventCube.i"
%include "GEvent.i"
%include "GEventAtom.i"
%include "GEventBin.i"
%include "GInstDir.i"
%include "GRoi.i"
%include "GPointing.i"
%include "GResponse.i"
%include "GPhoton.i"
%include "GPhotons.i"
%include "GSource.i"
