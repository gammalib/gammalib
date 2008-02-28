/***************************************************************************
 *          GFitsTableDblCol.hpp  - FITS table double column class         *
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

#ifndef GFITSTABLEDBLCOL_HPP
#define GFITSTABLEDBLCOL_HPP

/* __ Includes ___________________________________________________________ */
#include "GFitsCfitsio.hpp"
#include "GFitsTableCol.hpp"

/* __ Namespaces _________________________________________________________ */


/* __ Structures _________________________________________________________ */


/***************************************************************************
 *                    GFitsTableDblCol class definition                    *
 ***************************************************************************/
class GFitsTableDblCol : public GFitsTableCol {

public:
    // Constructors and destructors
    GFitsTableDblCol();
    GFitsTableDblCol(const GFitsTableDblCol& column);
    virtual ~GFitsTableDblCol();

    // Operators
    GFitsTableDblCol& operator= (const GFitsTableDblCol& column);

    // Methods
    std::string       string(const int& row, const int& col = 0);
    double            real(const int& row, const int& col = 0);
    int               integer(const int& row, const int& col = 0);
    GFitsTableDblCol* clone(void) const;
    std::string*      ptr_string(void);
    double*           ptr_double(void);
    float*            ptr_float(void);
    short*            ptr_short(void);
    long*             ptr_long(void);
    int*              ptr_int(void);
    
private:
    // Private methods
    void init_members(void);
    void copy_members(const GFitsTableDblCol& column);
    void free_members(void);
    void load(void);

    // Private data area
    int     m_size;          // Size of data area
    double* m_data;          // Data area
    double  m_nulval;        // NULL value
};


/***************************************************************************
 *                              Inline methods                             *
 ***************************************************************************/
inline double* GFitsTableDblCol::ptr_double(void) { return m_data; }
inline 
GFitsTableDblCol* GFitsTableDblCol::clone(void) const 
{
    return new GFitsTableDblCol(*this);
}

#endif /* GFITSTABLEDBLCOL_HPP */
