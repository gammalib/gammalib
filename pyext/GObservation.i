/***************************************************************************
 *    GObservation.i  -  Observation abstract base class SWIG interface    *
 * ----------------------------------------------------------------------- *
 *  copyright            : (C) 2008 by Jurgen Knodlseder                   *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GObservation.i
 * @brief GObservation class SWIG file.
 * @author J. Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GObservation.hpp"
%}


/***********************************************************************//**
 * @class GObservation
 *
 * @brief Abstract interface for the observation classes.
 ***************************************************************************/
class GObservation {
public:
    // Constructors and destructors
    GObservation();
    GObservation(const GObservation& obs);
    virtual ~GObservation();

    // Methods
};
