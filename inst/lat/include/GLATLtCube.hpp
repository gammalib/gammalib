/***************************************************************************
 *                GLATLtCube.hpp  -  Fermi/LAT livetime cube               *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2012 by Juergen Knoedlseder                         *
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
 * @author J. Knoedlseder
 */

#ifndef GLATLTCUBE_HPP
#define GLATLTCUBE_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include <iostream>
#include "GLog.hpp"
#include "GLATLtCubeMap.hpp"
#include "GLATAeff.hpp"
#include "GLATPsf.hpp"
#include "GGti.hpp"
#include "GSkyDir.hpp"
#include "GEnergy.hpp"


/***********************************************************************//**
 * @class GLATLtCube
 *
 * @brief Interface for the Fermi LAT livetime cube
 *
 * The livetime cube holds the livetime as function and zenith and azimuth
 * angles for a given observation. The azimuth dependence is optional. 
 ***************************************************************************/
class GLATLtCube {

    // I/O friends
    friend std::ostream& operator<< (std::ostream& os, const GLATLtCube& cube);
    friend GLog&         operator<< (GLog& log, const GLATLtCube& cube);

public:
    // Constructors and destructors
    GLATLtCube(void);
    GLATLtCube(const std::string& filename);
    GLATLtCube(const GLATLtCube& cube);
    virtual ~GLATLtCube(void);

    // Operators
    GLATLtCube& operator= (const GLATLtCube& cube);
    double      operator() (const GSkyDir& dir, const GEnergy& energy,
                            _ltcube_ctheta fct);
    double      operator() (const GSkyDir& dir, const GEnergy& energy,
                            _ltcube_ctheta_phi fct);
    double      operator() (const GSkyDir& dir, const GEnergy& energy,
                            const GLATAeff& aeff);
    double      operator() (const GSkyDir& dir, const GEnergy& energy,
                            const double& offset, const GLATPsf& psf,
                            const GLATAeff& aeff);

    // Methods
    void        clear(void);
    GLATLtCube* clone(void) const;
    void        load(const std::string& filename);
    void        save(const std::string& filename, bool clobber=false) const;
    std::string print(void) const;

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

#endif /* GLATLTCUBE_HPP */
