/***************************************************************************
 *                GEbounds.hpp - Energy boundaries class                   *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2009-2020 by Juergen Knoedlseder                         *
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
#include "GContainer.hpp"
#include "GEnergy.hpp"

/* __ Forward declarations _______________________________________________ */
class GXmlElement;
class GFits;
class GFitsTable;
class GFilename;
class GEnergies;

/* __ Constants __________________________________________________________ */
namespace gammalib {
    const std::string extname_ebounds = "EBOUNDS";
}


/***********************************************************************//**
 * @class GEbounds
 *
 * @brief Energy boundaries container class
 *
 * This class holds a list of energy intervals that are defined by a minimum
 * and maximum energy. Energies are implement using the GEnergy class which
 * holds unit independent energy information.
 *
 * The class has no method for sorting of the energy boundaries; it is
 * expected that the energy boundaries are correctly set by the client.
 ***************************************************************************/
class GEbounds : public GContainer {

    // Operator friends
    friend bool operator==(const GEbounds& a, const GEbounds& b);
    friend bool operator!=(const GEbounds& a, const GEbounds& b);

public:
    // Constructors and destructors
    GEbounds(void);
    explicit GEbounds(const GFilename& filename);
    explicit GEbounds(const GXmlElement& xml);
    explicit GEbounds(const GEnergies& energies);
    GEbounds(const GEbounds& ebds);
    GEbounds(const GEnergy& emin, const GEnergy& emax);
    GEbounds(const int&         num,
             const GEnergy&     emin,
             const GEnergy&     emax,
             const std::string& method = "LOG",
             const double&      gamma  = 1.0);
    virtual ~GEbounds(void);

    // Operators
    GEbounds& operator=(const GEbounds& ebds);

    // Methods
    void           clear(void);
    GEbounds*      clone(void) const;
    std::string    classname(void) const;
    int            size(void) const;
    bool           is_empty(void) const;
    void           append(const GEnergy& emin, const GEnergy& emax);
    void           insert(const GEnergy& emin, const GEnergy& emax);
    void           merge(void);
    void           merge(const GEnergy& emin, const GEnergy& emax);
    void           remove(const int& index);
    void           reserve(const int& num);
    void           extend(const GEbounds& ebds);
    void           set(const GEnergies& energies);
    void           set(const int&         num,
                       const GEnergy&     emin,
                       const GEnergy&     emax,
                       const std::string& method = "LOG",
                       const double&      gamma  = 1.0);
    void           load(const GFilename& filename);
    void           save(const GFilename& filename,
                        const bool& clobber = false,
                        const std::string& unit = "keV") const;
    void           read(const GFitsTable& table);
    void           write(GFits& file,
                         const std::string& extname = gammalib::extname_ebounds,
                         const std::string& unit = "keV") const;
    void           read(const GXmlElement& xml);
    void           write(GXmlElement& xml) const;
    int            index(const GEnergy& eng) const;
    const GEnergy& emin(void) const;
    const GEnergy& emax(void) const;
    void           emin(const int& index, const GEnergy& energy);
    void           emax(const int& index, const GEnergy& energy);
    GEnergy        emin(const int& index) const;
    GEnergy        emax(const int& index) const;
    GEnergy        emean(const int& index) const;
    GEnergy        elogmean(const int& index) const;
    GEnergy        ewidth(const int& index) const;
    bool           contains(const GEnergy& eng) const;
    bool           contains(const GEnergy& emin, const GEnergy& emax) const;
    std::string    print(const GChatter& chatter = NORMAL) const;

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GEbounds& ebds);
    void free_members(void);
    void set_attributes(void);
    void insert_eng(const int& index, const GEnergy& emin, const GEnergy& emax);

    // Protected data area
    int      m_num;         //!< Number of energy boundaries
    GEnergy  m_emin;        //!< Minimum energy of all intervals
    GEnergy  m_emax;        //!< Maximum energy of all intervals
    GEnergy* m_min;         //!< Array of interval minimum energies
    GEnergy* m_max;         //!< Array of interval maximum energies
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GEbounds").
 ***************************************************************************/
inline
std::string GEbounds::classname(void) const
{
    return ("GEbounds");
}


/***********************************************************************//**
 * @brief Return number of energy boundaries
 *
 * @return Number of energy boundaries.
 ***************************************************************************/
inline
int GEbounds::size(void) const
{
    return m_num;
}


/***********************************************************************//**
 * @brief Signal if there are no energy boundaries
 *
 * @return True if there are no energy boundaries.
 ***************************************************************************/
inline
bool GEbounds::is_empty(void) const
{
    return (m_num == 0);
}


/***********************************************************************//**
 * @brief Return minimum energy of all intervals
 *
 * @return Minimum energy of all intervals.
 ***************************************************************************/
inline
const GEnergy& GEbounds::emin(void) const
{
    return m_emin;
}


/***********************************************************************//**
 * @brief Return maximum energy of all intervals
 *
 * @return Maximum energy of all intervals.
 ***************************************************************************/
inline
const GEnergy& GEbounds::emax(void) const
{
    return m_emax;
}


/***********************************************************************//**
 * @brief Energy boundaries inequality operator friend
 *
 * @param[in] a First energy boundaries.
 * @param[in] b Second energy boundaries.
 * @return True if both energy boundaries are different.
 ***************************************************************************/
inline
bool operator!=(const GEbounds& a, const GEbounds& b)
{
    return (!(a == b));
}

#endif /* GEBOUNDS_HPP */
