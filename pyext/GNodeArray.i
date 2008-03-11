/***************************************************************************
 *          GNodeArray.i  -  Array of nodes class SWIG definition          *
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
 * @file GNodeArray.i
 * @brief GNodeArray class SWIG file.
 * @author J. Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GNodeArray.hpp"
%}


/***********************************************************************//**
 * @class GNodeArray
 *
 * @brief SWIG interface for the node array class.
 ***************************************************************************/
class GNodeArray {
public:
    // Constructors and destructors
    GNodeArray();
    GNodeArray(const GNodeArray& array);
    ~GNodeArray();

    // Methods
    void   nodes(const int& num, const double* array);
    void   nodes(const GVector& vector);
    void   set_value(const double& value);
    int    inx_left(void);
    int    inx_right(void);
    double wgt_left(void);
    double wgt_right(void);
};
