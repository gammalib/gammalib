/***************************************************************************
 *                   GPhotons.hpp - Photon container class                 *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2012 by Juergen Knoedlseder                              *
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
 * @file GPhotons.hpp
 * @brief Photon container class definition
 * @author Juergen Knoedlseder
 */

#ifndef GPHOTONS_HPP
#define GPHOTONS_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GBase.hpp"
#include "GPhoton.hpp"


/***********************************************************************//**
 * @class GPhotons
 *
 * @brief Container class for photons.
 ***************************************************************************/
class GPhotons : public GBase {

public:
    // Constructors and destructors
    GPhotons(void);
    GPhotons(const GPhotons& photons);
    virtual ~GPhotons(void);
 
    // Operators
    GPhotons&      operator=(const GPhotons& photons);
    GPhoton&       operator[](const int& index);
    const GPhoton& operator[](const int& index) const;

    // Methods
    void           clear(void);
    GPhotons*      clone(void) const;
    int            size(void) const { return m_photons.size(); }
    void           append(const GPhoton& photon);
    void           reserve(const int& number);
    std::string    print(void) const;
  
protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GPhotons& photons);
    void free_members(void);

    // Protected data members
    std::vector<GPhoton> m_photons;      //!< Photon list
};
#endif /* GPHOTONS_HPP */
