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

#ifndef GFITSTABLESTRCOL_HPP
#define GFITSTABLESTRCOL_HPP

/* __ Includes ___________________________________________________________ */
#include "GFitsCfitsio.hpp"
#include "GFitsTableCol.hpp"

/* __ Namespaces _________________________________________________________ */


/* __ Structures _________________________________________________________ */


/***************************************************************************
 *                    GFitsTableStrCol class definition                    *
 ***************************************************************************/
class GFitsTableStrCol : public GFitsTableCol {

public:
    // Constructors and destructors
    GFitsTableStrCol();
    GFitsTableStrCol(const GFitsTableStrCol& column);
    ~GFitsTableStrCol();

    // Operators
    GFitsTableStrCol& operator= (const GFitsTableStrCol& column);

    // Methods
    std::string       string(const int& row, const int& col = 0);
    double            real(const int& row, const int& col = 0);
    int               integer(const int& row, const int& col = 0);
    GFitsTableStrCol* clone(void) const;
    float*            ptr_float(void);
    double*           ptr_double(void);
    short*            ptr_short(void);
    long*             ptr_long(void);
    int*              ptr_int(void);
    void              set_nullstr(const std::string string);

private:
    // Private methods
    void init_members(void);
    void copy_members(const GFitsTableStrCol& column);
    void free_members(void);
    void load(void);

    // Private data area
    int    m_size;       // Total number of strings
    int    m_num_subs;   // Number of substrings
    char** m_data;       // Data area
    char*  m_nulstr;     // NULL string
    int    m_anynul;     // Number of NULLs encountered
};


/***************************************************************************
 *                              Inline methods                             *
 ***************************************************************************/
inline
GFitsTableStrCol* GFitsTableStrCol::clone(void) const
{
    return new GFitsTableStrCol(*this);
}

#endif /* GFITSTABLESTRCOL_HPP */
