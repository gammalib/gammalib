/***************************************************************************
 *      GFitsTableUStCol.hpp  - FITS table unsigned short column class     *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010 by Jurgen Knodlseder                                *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GFitsTableUStCol.hpp
 * @brief GFitsTableUStCol class definition.
 * @author J. Knodlseder
 */

#ifndef GFITSTABLEUSTCOL_HPP
#define GFITSTABLEUSTCOL_HPP

/* __ Includes ___________________________________________________________ */
#include "GFitsTableCol.hpp"


/***********************************************************************//**
 * @class GFitsTableUStCol
 *
 * @brief Interface for FITS table unsigned long integer column
 *
 * This class implements a FITS table unsigned long integer column.
 ***************************************************************************/
class GFitsTableUStCol : public GFitsTableCol {

    // I/O friends
    friend std::ostream& operator<< (std::ostream& os, const GFitsTableUStCol& column);

public:
    // Constructors and destructors
    GFitsTableUStCol();
    GFitsTableUStCol(const std::string& name, const int& length,
                     const int& size = 1);
    GFitsTableUStCol(const GFitsTableUStCol& column);
    virtual ~GFitsTableUStCol();

    // Operators
    GFitsTableUStCol&     operator= (const GFitsTableUStCol& column);
    unsigned short&       operator() (const int& row, const int& inx = 0);
    const unsigned short& operator() (const int& row, const int& inx = 0) const;

    // Methods
    std::string     string(const int& row, const int& col = 0);
    double          real(const int& row, const int& col = 0);
    int             integer(const int& row, const int& col = 0);
    unsigned short* data(void);
    void            set_nullval(const unsigned short* value);

private:
    // Private methods
    void              init_members(void);
    void              copy_members(const GFitsTableUStCol& column);
    void              free_members(void);
    void              save(void);
    GFitsTableUStCol* clone(void) const;
    std::string       ascii_format(void) const;
    std::string       binary_format(void) const;
    void              alloc_data(void);
    void              init_data(void);
    void              fetch_data(void);
    void*             ptr_data(void) { return m_data; }
    void*             ptr_nulval(void) { return m_nulval; }

    // Private data area
    unsigned short* m_data;          //!< Data area
    unsigned short* m_nulval;        //!< NULL value
};

#endif /* GFITSTABLEUSTCOL_HPP */
