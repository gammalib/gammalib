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

/* __ Namespaces _________________________________________________________ */


/* __ Structures _________________________________________________________ */


/***********************************************************************//**
 * @class GFitsTableFltCol
 *
 * @brief Interface for FITS table floating point column
 *
 * This class implements a FITS table floating point column.
 ***************************************************************************/
class GFitsTableFltCol : public GFitsTableCol {

public:
    // Constructors and destructors
    GFitsTableFltCol();
    GFitsTableFltCol(const GFitsTableFltCol& column);
    virtual ~GFitsTableFltCol();

    // Operators
    GFitsTableFltCol& operator= (const GFitsTableFltCol& column);

    // Methods
    std::string       string(const int& row, const int& col = 0);
    double            real(const int& row, const int& col = 0);
    int               integer(const int& row, const int& col = 0);
    GFitsTableFltCol* clone(void) const;
    float*            data(void);
    void              set_nullval(const float* value);

private:
    // Private methods
    void init_members(void);
    void copy_members(const GFitsTableFltCol& column);
    void free_members(void);
    void load(void);

    // Private data area
    int    m_size;          //!< Size of data area
    float* m_data;          //!< Data area
    float* m_nulval;        //!< NULL value
    int    m_anynul;        //!< Number of NULLs encountered
};


/***************************************************************************
 *                              Inline methods                             *
 ***************************************************************************/

#endif /* GFITSTABLEFLTCOL_HPP */
