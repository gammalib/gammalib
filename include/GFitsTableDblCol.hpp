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
/**
 * @file GFitsTableDblCol.hpp
 * @brief GFitsTableDblCol class definition.
 * @author J. Knodlseder
 */

#ifndef GFITSTABLEDBLCOL_HPP
#define GFITSTABLEDBLCOL_HPP

/* __ Includes ___________________________________________________________ */
#include "GFitsCfitsio.hpp"
#include "GFitsTableCol.hpp"

/* __ Namespaces _________________________________________________________ */


/* __ Structures _________________________________________________________ */


/***********************************************************************//**
 * @class GFitsTableDblCol
 *
 * @brief Interface for FITS table double precision column
 *
 * This class implements a FITS table double precision column.
 ***************************************************************************/
class GFitsTableDblCol : public GFitsTableCol {

public:
    // Constructors and destructors
    GFitsTableDblCol();
    GFitsTableDblCol(const int& length, const int& size = 1);
    GFitsTableDblCol(const GFitsTableDblCol& column);
    virtual ~GFitsTableDblCol();

    // Operators
    GFitsTableDblCol& operator= (const GFitsTableDblCol& column);
    double&           operator() (const int& row, const int& inx = 0);
    const double&     operator() (const int& row, const int& inx = 0) const;

    // Methods
    std::string       string(const int& row, const int& inx = 0);
    double            real(const int& row, const int& inx = 0);
    int               integer(const int& row, const int& inx = 0);
    GFitsTableDblCol* clone(void) const;
    double*           data(void);
    void              set_nullval(const double* value);

private:
    // Private methods
    void init_members(void);
    void copy_members(const GFitsTableDblCol& column);
    void free_members(void);
    void fetch_data(void);

    // Private data area
    int     m_size;          //!< Size of data area
    double* m_data;          //!< Data area
    double* m_nulval;        //!< NULL value
    int     m_anynul;        //!< Number of NULLs encountered
};

#endif /* GFITSTABLEDBLCOL_HPP */
