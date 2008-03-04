/***************************************************************************
 *           GFitsTableLngCol.hpp  - FITS table long column class          *
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
 * @file GFitsTableLngCol.hpp
 * @brief GFitsTableLngCol class definition.
 * @author J. Knodlseder
 */

#ifndef GFITSTABLELNGCOL_HPP
#define GFITSTABLELNGCOL_HPP

/* __ Includes ___________________________________________________________ */
#include "GFitsCfitsio.hpp"
#include "GFitsTableCol.hpp"

/* __ Namespaces _________________________________________________________ */


/* __ Structures _________________________________________________________ */


/***********************************************************************//**
 * @class GFitsTableLngCol
 *
 * @brief Interface for FITS table long integer column
 *
 * This class implements a FITS table long integer column.
 ***************************************************************************/
class GFitsTableLngCol : public GFitsTableCol {

public:
    // Constructors and destructors
    GFitsTableLngCol();
    GFitsTableLngCol(const GFitsTableLngCol& column);
    virtual ~GFitsTableLngCol();

    // Operators
    GFitsTableLngCol& operator= (const GFitsTableLngCol& column);
    long&             operator() (const int& row, const int& inx = 0);
    const long&       operator() (const int& row, const int& inx = 0) const;

    // Methods
    std::string       string(const int& row, const int& col = 0);
    double            real(const int& row, const int& col = 0);
    int               integer(const int& row, const int& col = 0);
    GFitsTableLngCol* clone(void) const;
    long*             data(void);
    void              set_nullval(const long* value);

private:
    // Private methods
    void init_members(void);
    void copy_members(const GFitsTableLngCol& column);
    void free_members(void);
    void fetch_data(void);

    // Private data area
    int   m_size;          //!< Size of data area
    long* m_data;          //!< Data area
    long* m_nulval;        //!< NULL value
    int   m_anynul;        //!< Number of NULLs encountered
};

#endif /* GFITSTABLELNGCOL_HPP */
