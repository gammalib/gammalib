/***************************************************************************
 *                GEbounds.hpp  -  Energy boundary class                   *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2009-2010 by Jurgen Knodlseder                           *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GEbounds.hpp
 * @brief Energy boundary class interface definition.
 * @author J. Knodlseder
 */

#ifndef GBOUNDS_HPP
#define GBOUNDS_HPP

/* __ Includes ___________________________________________________________ */
#include <iostream>
#include "GFits.hpp"
#include "GFitsBinTable.hpp"
#include "GEnergy.hpp"


/***********************************************************************//**
 * @class GEbounds
 *
 * @brief Interface for the GEbounds class.
 *
 * This class holds a list of energy intervals that are used for science
 * analysis.
 ***************************************************************************/
class GEbounds {

    // I/O friends
    friend std::ostream& operator<< (std::ostream& os, const GEbounds& ebds);

public:
    // Constructors and destructors
    GEbounds(void);
    GEbounds(const GEbounds& ebds);
    virtual ~GEbounds(void);

    // Operators
    GEbounds& operator= (const GEbounds& ebds);

    // Methods
    void    clear(void);
    void    append(const GEnergy& emin, const GEnergy& emax);
    void    insert(const GEnergy& emin, const GEnergy& emax);
    void    setlin(const GEnergy& emin, const GEnergy& emax, const int& num);
    void    setlog(const GEnergy& emin, const GEnergy& emax, const int& num);
	void    load(const std::string& filename,
                 const std::string& extname = "EBOUNDS");
	void    save(const std::string& filename, bool clobber,
                 const std::string& extname = "EBOUNDS");
    void    read(GFitsBinTable* hdu);
    void    write(GFits* file, const std::string& extname = "EBOUNDS");
    int     index(const GEnergy& eng) const;
    int     size(void) const { return m_num; }
    GEnergy emin(void) const { return m_emin; }
    GEnergy emax(void) const { return m_emax; }
    GEnergy emin(int inx) const;
    GEnergy emax(int inx) const;
    GEnergy emean(int inx) const;
    GEnergy elogmean(int inx) const;
    bool    isin(const GEnergy& eng) const;
    
protected:
    // Protected methods
    void      init_members(void);
    void      copy_members(const GEbounds& ebds);
    void      free_members(void);
    void      set_attributes(void);
    GEbounds* clone(void) const;
    void      insert_eng(int inx, const GEnergy& emin, const GEnergy& emax);
    void      merge_engs(void);

    // Protected data area
	int      m_num;         //!< Number of energy boundaries
    GEnergy  m_emin;        //!< Minimum energy covered
    GEnergy  m_emax;        //!< Maximum energy covered
    GEnergy* m_min;         //!< Energy bin minima
    GEnergy* m_max;         //!< Energy bin maxima
};

#endif /* GBOUNDS_HPP */
