/***************************************************************************
 *           GFitsTableShtCol.hpp  - FITS table short column class         *
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
 * @file GFitsTableShtCol.hpp
 * @brief GFitsTableShtCol class definition.
 * @author J. Knodlseder
 */

#ifndef GFITSTABLESHTCOL_HPP
#define GFITSTABLESHTCOL_HPP

/* __ Includes ___________________________________________________________ */
#include "GFitsCfitsio.hpp"
#include "GFitsTableCol.hpp"

/* __ Namespaces _________________________________________________________ */


/* __ Structures _________________________________________________________ */


/***********************************************************************//**
 * @class GFitsTableShtCol
 *
 * @brief Interface for FITS table short integer column
 *
 * This class implements a FITS table short integer column.
 ***************************************************************************/
class GFitsTableShtCol : public GFitsTableCol {

    // I/O friends
    friend ostream& operator<< (ostream& os, const GFitsTableShtCol& column);

public:
    // Constructors and destructors
    GFitsTableShtCol();
    GFitsTableShtCol(const std::string& name, const int& length,
                     const int& size = 1);
    GFitsTableShtCol(const GFitsTableShtCol& column);
    virtual ~GFitsTableShtCol();

    // Operators
    GFitsTableShtCol& operator= (const GFitsTableShtCol& column);
    short&            operator() (const int& row, const int& inx = 0);
    const short&      operator() (const int& row, const int& inx = 0) const;

    // Methods
    void              save(void);
    std::string       string(const int& row, const int& col = 0);
    double            real(const int& row, const int& col = 0);
    int               integer(const int& row, const int& col = 0);
    GFitsTableShtCol* clone(void) const;
    short*            data(void);
    void              set_nullval(const short* value);

private:
    // Private methods
    void        init_members(void);
    void        copy_members(const GFitsTableShtCol& column);
    void        free_members(void);
    std::string ascii_format(void) const;
    std::string binary_format(void) const;
    void        alloc_data(void);
    void        init_data(void);
    void        fetch_data(void);
    void*       ptr_data(void) { return m_data; }
    void*       ptr_nulval(void) { return m_nulval; }

    // Private data area
    short* m_data;       //!< Data area
    short* m_nulval;     //!< NULL value
};

#endif /* GFITSTABLESHTCOL_HPP */
