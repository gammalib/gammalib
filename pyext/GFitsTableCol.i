/***************************************************************************
 *    GFitsTableCol.i  - FITS table column abstract base class SWIG file   *
 * ----------------------------------------------------------------------- *
 *  copyright : (C) 2008-2010 by Jurgen Knodlseder                         *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GFitsTableCol.i
 * @brief GFitsTableCol class SWIG file.
 * @author J. Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GFitsTableCol.hpp"
%}


/***********************************************************************//**
 * @class GFitsTableCol
 *
 * @brief Abstract SWIG interface for FITS table column
 ***************************************************************************/
class GFitsTableCol {
public:
    // Constructors and destructors
    GFitsTableCol(void);
    GFitsTableCol(const std::string& name, const int& length,
                  const int& number,       const int& width);
    GFitsTableCol(const GFitsTableCol& column);
    virtual ~GFitsTableCol(void);

    // Virtual Methods
    virtual std::string string(const int& row, const int& inx = 0) = 0;
    virtual double      real(const int& row, const int& inx = 0) = 0;
    virtual int         integer(const int& row, const int& inx = 0) = 0;

    // Base class Methods
    void        name(const std::string& name);
    std::string name(void);
    int         colnum(void);
    int         type(void);
    int         repeat(void);
    int         width(void);
    int         number(void);
    int         length(void);
};
