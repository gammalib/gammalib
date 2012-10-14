/***************************************************************************
 *               GEbounds.hpp  -  Energy boundaries class                  *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2009-2012 by Juergen Knoedlseder                         *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *  This program is free software: you can redistribute it and/or modify   *
 *  it under the terms of the GNU General Public License as published by   *
 *  the Free Software Foundation, either version 3 of the License, or      *
 *  (at your option) any later version.                                    *
 *                                                                         *
 *  This program is distributed in the hope that it will be useful,        *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of         *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          *
 *  GNU General Public License for more details.                           *
 *                                                                         *
 *  You should have received a copy of the GNU General Public License      *
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.  *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GEbounds.hpp
 * @brief Energy boundaries class interface definition
 * @author Juergen Knoedlseder
 */

#ifndef GEBOUNDS_HPP
#define GEBOUNDS_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GBase.hpp"
#include "GFits.hpp"
#include "GFitsTable.hpp"
#include "GEnergy.hpp"


/***********************************************************************//**
 * @class GEbounds
 *
 * @brief Interface for the energy boundaries class
 *
 * This class holds a list of energy intervals.
 ***************************************************************************/
class GEbounds : public GBase {

public:
    // Constructors and destructors
    GEbounds(void);
    GEbounds(const GEbounds& ebds);
    explicit GEbounds(const std::string& filename,
                      const std::string& extname = "EBOUNDS");
    virtual ~GEbounds(void);

    // Operators
    GEbounds& operator= (const GEbounds& ebds);

    // Methods
    void        clear(void);
    GEbounds*   clone(void) const;
    void        append(const GEnergy& emin, const GEnergy& emax);
    void        insert(const GEnergy& emin, const GEnergy& emax);
    void        setlin(const GEnergy& emin, const GEnergy& emax, const int& num);
    void        setlog(const GEnergy& emin, const GEnergy& emax, const int& num);
    void        load(const std::string& filename,
                     const std::string& extname = "EBOUNDS");
    void        save(const std::string& filename, bool clobber,
                     const std::string& extname = "EBOUNDS") const;
    void        read(GFitsTable* hdu);
    void        write(GFits* file, const std::string& extname = "EBOUNDS") const;
    int         index(const GEnergy& eng) const;
    int         size(void) const { return m_num; }
    GEnergy     emin(void) const { return m_emin; }
    GEnergy     emax(void) const { return m_emax; }
    GEnergy     emin(int inx) const;
    GEnergy     emax(int inx) const;
    GEnergy     emean(int inx) const;
    GEnergy     elogmean(int inx) const;
    bool        isin(const GEnergy& eng) const;
    std::string print(void) const;


protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GEbounds& ebds);
    void free_members(void);
    void set_attributes(void);
    void insert_eng(int inx, const GEnergy& emin, const GEnergy& emax);
    void merge_engs(void);

    // Protected data area
    int      m_num;         //!< Number of energy boundaries
    GEnergy  m_emin;        //!< Minimum energy covered
    GEnergy  m_emax;        //!< Maximum energy covered
    GEnergy* m_min;         //!< Energy bin minima
    GEnergy* m_max;         //!< Energy bin maxima
};

#endif /* GEBOUNDS_HPP */
