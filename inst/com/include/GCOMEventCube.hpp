/***************************************************************************
 *          GCOMEventCube.hpp - COMPTEL event bin container class          *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2012-2016 by Juergen Knoedlseder                         *
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
 * @file GCOMEventCube.hpp
 * @brief COMPTEL event bin container class interface definition
 * @author Juergen Knoedlseder
 */

#ifndef GCOMEVENTCUBE_HPP
#define GCOMEVENTCUBE_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include <vector>
#include "GEventCube.hpp"
#include "GCOMEventBin.hpp"
#include "GCOMInstDir.hpp"
#include "GSkyDir.hpp"
#include "GSkyMap.hpp"
#include "GEnergy.hpp"
#include "GTime.hpp"

/* __ Forward declarations _______________________________________________ */
class GFilename;
class GEbounds;
class GGti;


/***********************************************************************//**
 * @class GCOMEventCube
 *
 * @brief COMPTEL event bin container class
 *
 * This class is a container class for COMPTEL event bins.
 ***************************************************************************/
class GCOMEventCube : public GEventCube {

public:
    // Constructors and destructors
    GCOMEventCube(void);
    explicit GCOMEventCube(const GFilename& filename);
    explicit GCOMEventCube(const GSkyMap& map, const GEbounds& ebds,
                           const GGti& gti, const double& phimin,
                           const double& dphi);
    GCOMEventCube(const GCOMEventCube& cube);
    virtual ~GCOMEventCube(void);

    // Operators
    virtual GCOMEventCube&      operator=(const GCOMEventCube& cube);
    virtual GCOMEventBin*       operator[](const int& index);
    virtual const GCOMEventBin* operator[](const int& index) const;

    // Implemented pure virtual base class methods
    virtual void           clear(void);
    virtual GCOMEventCube* clone(void) const;
    virtual std::string    classname(void) const;
    virtual int            size(void) const;
    virtual int            dim(void) const;
    virtual int            naxis(const int& axis) const;
    virtual void           load(const GFilename& filename);
    virtual void           save(const GFilename& filename,
                                const bool&      clobber = false) const;
    virtual void           read(const GFits& file);
    virtual void           write(GFits& file) const;
    virtual int            number(void) const;
    virtual std::string    print(const GChatter& chatter = NORMAL) const;

    // Other methods
    void                   map(const GSkyMap& map, const double& phimin,
                               const double& dphi);
    const GSkyMap&         map(void) const;
    int                    nchi(void) const;
    int                    npsi(void) const;
    int                    nphi(void) const;
    int                    npix(void) const;

protected:
    // Protected methods
    void         init_members(void);
    void         copy_members(const GCOMEventCube& cube);
    void         free_members(void);
    void         set_scatter_directions(void);
    void         set_scatter_angles(const double& phimin, const double& dphi);
    virtual void set_energies(void);
    virtual void set_times(void);
    void         init_bin(void);
    void         set_bin(const int& index);

    // Protected members
    GCOMEventBin         m_bin;        //!< Actual event bin
    GCOMInstDir          m_dir;        //!< Actual event direction
    GSkyMap              m_map;        //!< Counts map stored as sky map
    GTime                m_time;       //!< Event cube mean time
    double               m_ontime;     //!< Event cube ontime (sec)
    GEnergy              m_energy;     //!< Event cube mean energy
    GEnergy              m_ewidth;     //!< Event cube energy bin width
    std::vector<GSkyDir> m_dirs;       //!< Array of event scatter directions
    std::vector<double>  m_solidangle; //!< Array of solid angles (sr)
    std::vector<double>  m_phi;        //!< Array of event scatter angles
    std::vector<double>  m_dphi;       //!< Array of event scatter angles widths
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GCOMEventCube").
 ***************************************************************************/
inline
std::string GCOMEventCube::classname(void) const
{
    return ("GCOMEventCube");
}


/***********************************************************************//**
 * @brief Return event cube sky map
 *
 * @return Sky map containing event cube.
 ***************************************************************************/
inline
const GSkyMap& GCOMEventCube::map(void) const
{
    return (m_map);
}


/***********************************************************************//**
 * @brief Return number of Chi bins
 *
 * @return Number of Chi bins.
 ***************************************************************************/
inline
int GCOMEventCube::nchi(void) const
{
    return (m_map.nx());
}


/***********************************************************************//**
 * @brief Return number of Psi bins
 *
 * @return Number of Psi bins.
 ***************************************************************************/
inline
int GCOMEventCube::npsi(void) const
{
    return (m_map.ny());
}


/***********************************************************************//**
 * @brief Return number of Phibar bins
 *
 * @return Number of Phibar bins.
 ***************************************************************************/
inline
int GCOMEventCube::nphi(void) const
{
    return (m_map.nmaps());
}


/***********************************************************************//**
 * @brief Return number of pixels
 *
 * @return Number of pixels.
 ***************************************************************************/
inline
int GCOMEventCube::npix(void) const
{
    return (m_map.npix());
}

#endif /* GCOMEVENTCUBE_HPP */
