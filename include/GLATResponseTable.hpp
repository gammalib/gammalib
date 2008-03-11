/***************************************************************************
 *         GLATResponseTable.hpp  -  GLAST LAT Response table class        *
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
 * @file GLATResponseTable.hpp
 * @brief GLATResponseTable class definition.
 * @author J. Knodlseder
 */

#ifndef GLATRESPONSETABLE_HPP
#define GLATRESPONSETABLE_HPP

/* __ Includes ___________________________________________________________ */
#include "GFitsHDU.hpp"

/* __ Namespaces _________________________________________________________ */


/***********************************************************************//**
 * @class GLATResponseTable
 *
 * @brief Interface for the GLAST LAT instrument response table classe.
 ***************************************************************************/
class GLATResponseTable {

public:
    // Constructors and destructors
    GLATResponseTable();
    GLATResponseTable(const GLATResponseTable& table);
    ~GLATResponseTable();

    // Operators
    GLATResponseTable& operator= (const GLATResponseTable & table);

    // Methods
    void   load(const GFitsHDU* hdu);
    void   save(GFitsHDU* hdu) const;
    int    index(const int& ie, const int& ic) const;
    double energy(const int& ie) const;
    int    num_energy(void) const { return m_energy_num; }
    int    num_ctheta(void) const { return m_ctheta_num; }

private:
    // Methods
    void init_members(void);
    void copy_members(const GLATResponseTable& table);
    void free_members(void);
    
    // Data
    int     m_energy_num;   //!< Number of energy bins in table
    int     m_ctheta_num;   //!< Number of cos theta bins in table
    double* m_energy_lo;    //!< Energy bins lower boundary (MeV)
    double* m_energy_hi;    //!< Energy bins upper boundary (MeV)
    double* m_ctheta_lo;    //!< cos(theta) bins lower boundary
    double* m_ctheta_hi;    //!< cos(theta) bins upper boundary
};

#endif /* GLATRESPONSETABLE_HPP */
