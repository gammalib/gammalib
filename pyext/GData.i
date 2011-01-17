/***************************************************************************
 *                   GData.i  -  Data class SWIG interface                 *
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
 * @file GData.i
 * @brief GData class SWIG file.
 * @author J. Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GData.hpp"
%}


/***********************************************************************//**
 * @class GData
 *
 * @brief Interface for the data class.
 ***************************************************************************/
class GData {
public:
    // Constructors and destructors
    GData();
    GData(const GData& data);
    virtual ~GData();

    // Methods
};
