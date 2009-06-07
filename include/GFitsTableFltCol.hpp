/***************************************************************************
 *           GFitsTableFltCol.hpp  - FITS table float column class         *
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
 * @file GFitsTableFltCol.hpp
 * @brief GFitsTableFltCol class definition.
 * @author J. Knodlseder
 */

#ifndef GFITSTABLEFLTCOL_HPP
#define GFITSTABLEFLTCOL_HPP

/* __ Includes ___________________________________________________________ */
#include "GFitsCfitsio.hpp"
#include "GFitsTableCol.hpp"


/***********************************************************************//**
 * @class GFitsTableFltCol
 *
 * @brief Interface for FITS table floating point column
 *
 * This class implements a FITS table floating point column.
 ***************************************************************************/
class GFitsTableFltCol : public GFitsTableCol {

    // I/O friends
    friend std::ostream& operator<< (std::ostream& os, const GFitsTableFltCol& column);

public:
    // Constructors and destructors
    GFitsTableFltCol();
    GFitsTableFltCol(const std::string& name, const int& length,
                     const int& size = 1);
    GFitsTableFltCol(const GFitsTableFltCol& column);
    virtual ~GFitsTableFltCol();

    // Operators
    GFitsTableFltCol& operator= (const GFitsTableFltCol& column);
    float&            operator() (const int& row, const int& inx = 0);
    const float&      operator() (const int& row, const int& inx = 0) const;

    // Methods
    std::string string(const int& row, const int& col = 0);
    double      real(const int& row, const int& col = 0);
    int         integer(const int& row, const int& col = 0);
    float*      data(void);
    void        set_nullval(const float* value);

private:
    // Private methods
    void              init_members(void);
    void              copy_members(const GFitsTableFltCol& column);
    void              free_members(void);
    void              save(void);
    GFitsTableFltCol* clone(void) const;
    std::string       ascii_format(void) const;
    std::string       binary_format(void) const;
    void              alloc_data(void);
    void              init_data(void);
    void              fetch_data(void);
    void*             ptr_data(void) { return m_data; }
    void*             ptr_nulval(void) { return m_nulval; }

    // Private data area
    float* m_data;       //!< Data area
    float* m_nulval;     //!< NULL value
};

#endif /* GFITSTABLEFLTCOL_HPP */
