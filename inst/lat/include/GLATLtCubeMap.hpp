/***************************************************************************
 *          GLATLtCubeMap.hpp - Fermi LAT livetime cube map class          *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2013 by Juergen Knoedlseder                         *
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
 * @file GLATLtCubeMap.hpp
 * @brief Fermi LAT livetime cube map class definition
 * @author Juergen Knoedlseder
 */

#ifndef GLATLTCUBEMAP_HPP
#define GLATLTCUBEMAP_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GBase.hpp"
#include "GFitsTable.hpp"
#include "GSkymap.hpp"
#include "GSkyDir.hpp"
#include "GLATAeff.hpp"
#include "GLATPsf.hpp"

/* __ Typedefs ___________________________________________________________ */
typedef double (*_ltcube_ctheta)(const double& costheta);
typedef double (*_ltcube_ctheta_phi)(const double& costheta, const double& phi);


/***********************************************************************//**
 * @class GLATLtCubeMap
 *
 * @brief Fermi LAT livetime cube map class
 *
 * A livetime cube map holds a set of HEALPix skymaps that are a function
 * of the cosine of the zenith angle and (optionally) of the azimuth angle.
 ***************************************************************************/
class GLATLtCubeMap : public GBase {

public:
    // Constructors and destructors
    GLATLtCubeMap(void);
    GLATLtCubeMap(const GLATLtCubeMap& map);
    virtual ~GLATLtCubeMap(void);

    // Operators
    GLATLtCubeMap& operator=(const GLATLtCubeMap& cube);
    double         operator()(const GSkyDir& dir, _ltcube_ctheta fct) const;
    double         operator()(const GSkyDir& dir, _ltcube_ctheta_phi fct) const;
    double         operator()(const GSkyDir& dir, const GEnergy& energy,
                              const GLATAeff& aeff) const;
    double         operator()(const GSkyDir& dir, const GEnergy& energy,
                              const double& offset, const GLATPsf& psf,
                              const GLATAeff& aeff) const;

    // Methods
    void           clear(void);
    GLATLtCubeMap* clone(void) const;
    std::string    classname(void) const;
    void           read(const GFitsTable& table);
    void           write(GFits& file) const;
    const int&     ncostheta(void) const;
    const int&     nphi(void) const;
    bool           has_phi(void) const;
    double         costheta(const int& index) const;
    double         phi(const int& index) const;
    const double&  costhetamin(void) const;
    std::string    costhetabin(void) const;
    std::string    print(const GChatter& chatter = NORMAL) const;

private:
    // Methods
    void init_members(void);
    void copy_members(const GLATLtCubeMap& cube);
    void free_members(void);
    
    // Protected members
    GSkymap m_map;          //!< Lifetime cube map
    int     m_num_ctheta;   //!< Number of bins in cos theta
    int     m_num_phi;      //!< Number of bins in phi
    double  m_min_ctheta;   //!< Minimum cos theta value
    bool    m_sqrt_bin;     //!< Square root binning?
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GLATLtCubeMap").
 ***************************************************************************/
inline
std::string GLATLtCubeMap::classname(void) const
{
    return ("GLATLtCubeMap");
}


/***********************************************************************//**
 * @brief Return number of cosine theta bins
 *
 * @return Number of cosine theta bins.
 ***************************************************************************/
inline
const int& GLATLtCubeMap::ncostheta(void) const
{
    return m_num_ctheta;
}


/***********************************************************************//**
 * @brief Return number of phi bins
 *
 * @return Number of phi bins.
 ***************************************************************************/
inline
const int& GLATLtCubeMap::nphi(void) const
{
    return m_num_phi;
}


/***********************************************************************//**
 * @brief Signal if livetime cube map has phi dependence
 *
 * @return True if livetime cube map has phi dependence.
 ***************************************************************************/
inline
bool GLATLtCubeMap::has_phi(void) const
{
    return (m_num_phi != 0);
}


/***********************************************************************//**
 * @brief Return minimum of cosine theta
 *
 * @return Minimum of cosine theta.
 ***************************************************************************/
inline
const double& GLATLtCubeMap::costhetamin(void) const
{
    return m_min_ctheta;
}

#endif /* GLATLTCUBEMAP_HPP */
