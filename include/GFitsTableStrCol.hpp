/***************************************************************************
 *          GFitsTableStrCol.hpp  - FITS table string column class         *
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
 * @file GFitsTableStrCol.hpp
 * @brief GFitsTableStrCol class definition.
 * @author J. Knodlseder
 */

#ifndef GFITSTABLESTRCOL_HPP
#define GFITSTABLESTRCOL_HPP

/* __ Includes ___________________________________________________________ */
#include "GFitsCfitsio.hpp"
#include "GFitsTableCol.hpp"

/* __ Namespaces _________________________________________________________ */


/* __ Structures _________________________________________________________ */


/***********************************************************************//**
 * @class GFitsTableStrCol
 *
 * @brief Interface for FITS table string column
 *
 * This class implements a FITS table string column.
 ***************************************************************************/
class GFitsTableStrCol : public GFitsTableCol {

    // I/O friends
    friend ostream& operator<< (ostream& os, const GFitsTableStrCol& column);

public:
    // Constructors and destructors
    GFitsTableStrCol();
    GFitsTableStrCol(const std::string& name, const int& length,
                     const int& width, const int& size = 1);
    GFitsTableStrCol(const GFitsTableStrCol& column);
    ~GFitsTableStrCol();

    // Operators
    GFitsTableStrCol&  operator= (const GFitsTableStrCol& column);
    std::string&       operator() (const int& row, const int& inx = 0);
    const std::string& operator() (const int& row, const int& inx = 0) const;

    // Methods
    std::string  string(const int& row, const int& col = 0);
    double       real(const int& row, const int& col = 0);
    int          integer(const int& row, const int& col = 0);
    std::string* data(void);
    void         set_nullval(const std::string string);

private:
    // Private methods
    void              init_members(void);
    void              copy_members(const GFitsTableStrCol& column);
    void              free_members(void);
    void              save(void);
    GFitsTableStrCol* clone(void) const;
    std::string       ascii_format(void) const;
    std::string       binary_format(void) const;
    void              alloc_data(void);
    void              init_data(void);
    void              fetch_data(void);
    void*             ptr_data(void) { return m_buffer; }
    void*             ptr_nulval(void) { return m_nulval; }
    void              alloc_buffer(void);
    void              free_buffer(void);

    // Private data area
    std::string* m_data;    //!< Data area
    char**       m_buffer;  //!< Data area for CFITSIO transfer
    char*        m_nulval;  //!< NULL string
};

#endif /* GFITSTABLESTRCOL_HPP */
