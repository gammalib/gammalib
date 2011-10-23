/***************************************************************************
 *            mwl  -  Multi-wavelength support python bindings             *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2011 by Jurgen Knodlseder                                *
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
 * swig -c++ -python -Wall mwl.i                                           *
 ***************************************************************************/
%module mwl
%feature("autodoc", "1");

/* __ Include standard typemaps for vectors and strings __________________ */
%include stl.i

/* __ Make sure that exceptions are catched ______________________________ */
%import(module="support") "GException.i";

/* __ Inform about base classes __________________________________________ */
%import(module="obs") "GObservation.i";
%import(module="obs") "GEvent.i";
%import(module="obs") "GEventBin.i";
%import(module="obs") "GEvents.i";
%import(module="obs") "GEventCube.i";

/* __ MWL ________________________________________________________________ */
%include "GMWLObservation.i"
%include "GMWLSpectrum.i"
%include "GMWLDatum.i"



