/***************************************************************************
 *        GFitsTableCol.hpp  - FITS table column abstract base class       *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2008 by Jurgen Knodlseder                                *
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

public:
    // Constructors and destructors
    GFitsTableCol(void);
    GFitsTableCol(const std::string& name, const int& length,
                  const int& number,       const int& width);
    GFitsTableCol(const GFitsTableCol& column);
    virtual ~GFitsTableCol(void);

    // Operators
    GFitsTableCol& operator= (const GFitsTableCol& column);

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

protected:
    // Protected data area
    std::string m_name;      //!< Column name
    std::string m_unit;      //!< Column unit
    int         m_colnum;    //!< @brief Column number (starting from 1).
                             //!< This parameter is used to signal if a
                             //!< table column corresponds to a FITS file
                             //!< column. If it is set to 0 there is no
                             //!< correspondance.
    int         m_type;      //!< Column type
    int         m_repeat;    //!< Repeat value of column
    int         m_width;     //!< Width of single column element
    int         m_number;    //!< @brief Number of elements in column.
                             //!< m_number = m_repeat / m_width
    int         m_length;    //!< Length of column
    int         m_size;      //!< Size of allocated data area (0 if not loaded)
    int         m_anynul;    //!< Number of NULLs encountered
    __fitsfile  m_fitsfile;  //!< FITS file associated with column

    // Protected virtual methods
    virtual GFitsTableCol* clone(void) const = 0;
    virtual std::string    ascii_format(void) const = 0;
    virtual std::string    binary_format(void) const = 0;
    virtual void           alloc_data(void) = 0;
    virtual void           init_data(void) = 0;
    virtual void*          ptr_data(void) = 0;
    virtual void*          ptr_nulval(void) = 0;

    // Protected methods
    virtual void save(void);
    virtual void fetch_data(void);
    virtual void load_column(void);
    virtual void save_column(void);
    virtual void dump_column(std::ostream& os, void* data) const;
    virtual int  offset(const int& row, const int& inx) const; 

private:
    // Private methods
    void init_members(void);
    void copy_members(const GFitsTableCol& column);
    void free_members(void);
    void connect(__fitsfile* fptr);
};

#endif /* GFITSTABLECOL_HPP */
