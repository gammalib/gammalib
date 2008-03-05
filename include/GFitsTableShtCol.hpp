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

#ifndef GFITSTABLESHTCOL_HPP
#define GFITSTABLESHTCOL_HPP

/* __ Includes ___________________________________________________________ */
#include "GFitsCfitsio.hpp"
#include "GFitsTableCol.hpp"

/* __ Namespaces _________________________________________________________ */


/* __ Structures _________________________________________________________ */


/***************************************************************************
 *                    GFitsTableShtCol class definition                    *
 ***************************************************************************/
class GFitsTableShtCol : public GFitsTableCol {

public:
    // Constructors and destructors
    GFitsTableShtCol();
    GFitsTableShtCol(const GFitsTableShtCol& column);
    virtual ~GFitsTableShtCol();

    // Operators
    GFitsTableShtCol& operator= (const GFitsTableShtCol& column);

    // Methods
    void              save(void);
    std::string       string(const int& row, const int& col = 0);
    double            real(const int& row, const int& col = 0);
    int               integer(const int& row, const int& col = 0);
    GFitsTableShtCol* clone(void) const;
    double*           ptr_double(void);
    float*            ptr_float(void);
    short*            ptr_short(void);
    long*             ptr_long(void);
    int*              ptr_int(void);
    void              set_nullval(const short* value);

private:
    // Private methods
    void init_members(void);
    void copy_members(const GFitsTableShtCol& column);
    void free_members(void);
    void load(void);

    // Private data area
    int    m_size;          // Size of data area
    short* m_data;          // Data area
    short* m_nulval;        // NULL value
    int    m_anynul;        // Number of NULLs encountered
};


/***************************************************************************
 *                              Inline methods                             *
 ***************************************************************************/
inline short* GFitsTableShtCol::ptr_short(void) { return m_data; }
inline
GFitsTableShtCol* GFitsTableShtCol::clone(void) const
{
    return new GFitsTableShtCol(*this);
}

#endif /* GFITSTABLESHTCOL_HPP */
