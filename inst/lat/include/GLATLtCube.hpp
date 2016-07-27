/***************************************************************************
 *               GLATLtCube.hpp - Fermi/LAT livetime cube class            *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2016 by Juergen Knoedlseder                         *
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
 * @file GLATLtCube.hpp
 * @brief Fermi/LAT livetime cube class definition
 * @author Juergen Knoedlseder
 */

#ifndef GLATLTCUBE_HPP
#define GLATLTCUBE_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GBase.hpp"
#include "GGti.hpp"
#include "GLATLtCubeMap.hpp"

/* __ Forward declarations _______________________________________________ */
class GFilename;
class GSkyDir;
class GEnergy;
class GLATAeff;
class GLATPsf;


/***********************************************************************//**
 * @class GLATLtCube
 *
 * @brief Interface for the Fermi LAT livetime cube
 *
 * The livetime cube holds the livetime as function and zenith and azimuth
 * angles for a given observation. The azimuth dependence is optional. 
 ***************************************************************************/
class GLATLtCube : public GBase {

public:
    // Constructors and destructors
    GLATLtCube(void);
    explicit GLATLtCube(const GFilename& filename);
    GLATLtCube(const GLATLtCube& cube);
    virtual ~GLATLtCube(void);

    // Operators
    GLATLtCube& operator= (const GLATLtCube& cube);
    double      operator()(const GSkyDir& dir, const GEnergy& energy,
                           _ltcube_ctheta fct) const;
    double      operator()(const GSkyDir& dir, const GEnergy& energy,
                           _ltcube_ctheta_phi fct) const;
    double      operator()(const GSkyDir& dir, const GEnergy& energy,
                           const GLATAeff& aeff) const;
    double      operator()(const GSkyDir& dir, const GEnergy& energy,
                           const double& offset, const GLATPsf& psf,
                           const GLATAeff& aeff) const;

    // Methods
    void        clear(void);
    GLATLtCube* clone(void) const;
    std::string classname(void) const;
    void        load(const GFilename& filename);
    void        save(const GFilename& filename,
                     const bool&      clobber = false) const;
    std::string print(const GChatter& chatter = NORMAL) const;

private:
    // Methods
    void init_members(void);
    void copy_members(const GLATLtCube& cube);
    void free_members(void);
    
    // Protected members
    GLATLtCubeMap m_exposure;
    GLATLtCubeMap m_weighted_exposure;
    GGti          m_gti;
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GLATLtCube").
 ***************************************************************************/
inline
std::string GLATLtCube::classname(void) const
{
    return ("GLATLtCube");
}

#endif /* GLATLTCUBE_HPP */
