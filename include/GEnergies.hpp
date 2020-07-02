/***************************************************************************
 *                   GEnergies.hpp - Energy container class                *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2013-2020 by Juergen Knoedlseder                         *
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
 * @file GEnergies.hpp
 * @brief Energy container class definition
 * @author Juergen Knoedlseder
 */

#ifndef GENERGIES_HPP
#define GENERGIES_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include <vector>
#include "GContainer.hpp"
#include "GEnergy.hpp"
#include "GFits.hpp"
#include "GFitsTable.hpp"

/* __ Forward declarations _______________________________________________ */
class GFilename;
class GEbounds;

/* __ Constants __________________________________________________________ */
namespace gammalib {
    const std::string extname_energies = "ENERGIES";
}


/***********************************************************************//**
 * @class GEnergies
 *
 * @brief Energy container class.
 *
 * This class is a container for energies. Energies are implemented by the
 * GEnergy class which stores energies in a unit independent way.
 *
 * An energy container can be constructed from a FITS file or by defining
 * a number of energies within a given interval, spaced either linearly or
 * logarithmically.
 ***************************************************************************/
class GEnergies : public GContainer {

public:
    // Constructors and destructors
    GEnergies(void);
    explicit  GEnergies(const GFilename& filename);
    explicit  GEnergies(const GEbounds& ebounds);
    GEnergies(const GEnergies& energies);
    GEnergies(const int&         num,
              const GEnergy&     emin,
              const GEnergy&     emax,
              const std::string& method = "LOG",
              const double&      gamma  = 1.0);
    virtual ~GEnergies(void);
 
    // Operators
    GEnergies&     operator=(const GEnergies& energies);
    GEnergy&       operator[](const int& index);
    const GEnergy& operator[](const int& index) const;

    // Methods
    void           clear(void);
    GEnergies*     clone(void) const;
    std::string    classname(void) const;
    GEnergy&       at(const int& index);
    const GEnergy& at(const int& index) const;
    int            size(void) const;
    bool           is_empty(void) const;
    GEnergy&       append(const GEnergy& energy);
    GEnergy&       insert(const int& index, const GEnergy& energy);
    void           remove(const int& index);
    void           reserve(const int& num);
    void           extend(const GEnergies& energies);
    void           set(const GEbounds& ebounds);
    void           set(const int&         num,
                       const GEnergy&     emin,
                       const GEnergy&     emax,
                       const std::string& method = "LOG",
                       const double&      gamma  = 1.0);
    void           load(const GFilename& filename);
    void           save(const GFilename& filename,
                        const bool&      clobber = false) const;
    void           read(const GFitsTable& table);
    void           write(GFits& file,
                         const std::string& extname = gammalib::extname_energies) const;
    std::string    print(const GChatter& chatter = NORMAL) const;
  
protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GEnergies& energies);
    void free_members(void);
    void set_lin(const int& num, const GEnergy& emin, const GEnergy& emax);
    void set_log(const int& num, const GEnergy& emin, const GEnergy& emax);
    void set_pow(const int& num, const GEnergy& emin, const GEnergy& emax,
                 const double& gamma);

    // Protected data members
    std::vector<GEnergy> m_energies;  //!< List of energies
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GEnergies").
 ***************************************************************************/
inline
std::string GEnergies::classname(void) const
{
    return ("GEnergies");
}


/***********************************************************************//**
 * @brief Return reference to energy
 *
 * @param[in] index Energy index [0,...,size()-1].
 *
 * Returns a reference to the energy with the specified @p index.
 ***************************************************************************/
inline
GEnergy& GEnergies::operator[](const int& index)
{
    return (m_energies[index]);
}


/***********************************************************************//**
 * @brief Return reference to energy (const version)
 *
 * @param[in] index Energy index [0,...,size()-1].
 *
 * Returns a reference to the energy with the specified @p index.
 ***************************************************************************/
inline
const GEnergy& GEnergies::operator[](const int& index) const
{
    return (m_energies[index]);
}


/***********************************************************************//**
 * @brief Return number of energies in container
 *
 * @return Number of energies in container.
 *
 * Returns the number of energies in the energy container.
 ***************************************************************************/
inline
int GEnergies::size(void) const
{
    return (int)m_energies.size();
}


/***********************************************************************//**
 * @brief Signals if there are no energies in container
 *
 * @return True if container is empty, false otherwise.
 *
 * Signals if the energy container does not contain any energy.
 ***************************************************************************/
inline
bool GEnergies::is_empty(void) const
{
    return (m_energies.empty());
}


/***********************************************************************//**
 * @brief Reserves space for energies in container
 *
 * @param[in] num Number of energies.
 *
 * Reserves space for @p num energies in the container.
 ***************************************************************************/
inline
void GEnergies::reserve(const int& num)
{
    m_energies.reserve(num);
    return;
}

#endif /* GENERGIES_HPP */
