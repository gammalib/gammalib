/***************************************************************************
 *        GFitsTableCol.hpp  - FITS table column abstract base class       *
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
 * @file GFitsTableCol.hpp
 * @brief GFitsTableCol class definition.
 * @author J. Knodlseder
 */

#ifndef GFITSTABLECOL_HPP
#define GFITSTABLECOL_HPP

/* __ Includes ___________________________________________________________ */
#include "GFitsCfitsio.hpp"

/* __ Namespaces _________________________________________________________ */


/***********************************************************************//**
 * @class GFitsTableCol
 *
 * @brief Abstract interface for FITS table column
 *
 * This class implements a FITS table column. Vector columns are supported.
 ***************************************************************************/
class GFitsTableCol {

    // Friend classes
    friend class GFitsTable;

    // I/O friends
    friend ostream& operator<< (ostream& os, const GFitsTableCol& column);

public:
    // Constructors and destructors
    GFitsTableCol();
    GFitsTableCol(const int& length, const int& size = 1);
    GFitsTableCol(const GFitsTableCol& column);
    virtual ~GFitsTableCol();

    // Operators
    GFitsTableCol& operator= (const GFitsTableCol& column);

    // Virtual Methods
    virtual void           save(void) = 0;
    virtual std::string    string(const int& row, const int& inx = 0) = 0;
    virtual double         real(const int& row, const int& inx = 0) = 0;
    virtual int            integer(const int& row, const int& inx = 0) = 0;
    virtual GFitsTableCol* clone(void) const = 0;

    // Base class Methods
    void        name(const std::string& name);
    std::string name(void);
    int         colnum(void);
    int         type(void);
    int         repeat(void);
    int         width(void);
    int         length(void);

protected:
    // Protected data area
    std::string m_name;        //!< Column name
    std::string m_format;      //!< Column format
    std::string m_unit;        //!< Column unit
    int         m_colnum;      //!< Column number in FITS file (starting from 1)
    int         m_type;        //!< Column type
    int         m_repeat;      //!< Repeat value of column
    int         m_width;       //!< Width of column
    int         m_length;      //!< Length of column
    __fitsfile  m_fitsfile;    //!< FITS file associated with column

private:
    // Private methods
    void init_members(void);
    void copy_members(const GFitsTableCol& column);
    void free_members(void);
    void connect(__fitsfile* fptr);
};

#endif /* GFITSTABLECOL_HPP */
