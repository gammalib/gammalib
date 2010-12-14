/***************************************************************************
 *            GLATLtCubeMap.hpp  -  Fermi LAT lifetime cube map            *
 * ----------------------------------------------------------------------- *
 *  copyright : (C) 2010 by Jurgen Knodlseder                              *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GLATLtCubeMap.hpp
 * @brief GLATLtCubeMap class definition.
 * @author J. Knodlseder
 */

#ifndef GLATLTCUBEMAP_HPP
#define GLATLTCUBEMAP_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include <iostream>
#include "GLog.hpp"
#include "GFitsTable.hpp"
#include "GSkymap.hpp"
#include "GSkyDir.hpp"

/* __ Typedefs ___________________________________________________________ */
typedef double (*_ltcube_ctheta)(const double& costheta);
typedef double (*_ltcube_ctheta_phi)(const double& costheta, const double& phi);


/***********************************************************************//**
 * @class GLATLtCubeMap
 *
 * @brief Interface for the Fermi LAT lifetime cube map.
 *
 * A lifetime cube map holds a set of HEALPix skymaps that are a function
 * of the cosine of the zenith angle and (optionally) of the azimuth angle.
 ***************************************************************************/
class GLATLtCubeMap {

    // I/O friends
    friend std::ostream& operator<< (std::ostream& os, const GLATLtCubeMap& map);
    friend GLog&         operator<< (GLog& log, const GLATLtCubeMap& map);

public:
    // Constructors and destructors
    GLATLtCubeMap(void);
    GLATLtCubeMap(const GLATLtCubeMap& map);
    virtual ~GLATLtCubeMap(void);

    // Operators
    GLATLtCubeMap& operator= (const GLATLtCubeMap& cube);
    double         operator() (const GSkyDir& dir, _ltcube_ctheta fct);
    double         operator() (const GSkyDir& dir, _ltcube_ctheta_phi fct);

    // Methods
    void           clear(void);
    GLATLtCubeMap* clone(void) const;
    void           read(const GFitsTable* hdu);
    void           write(GFits* file) const;
    int            ncostheta(void) const { return m_num_ctheta; }
    int            nphi(void) const { return m_num_phi; }
    double         costheta(const int& index) const;
    double         phi(const int& index) const;
    double         costhetamin(void) const { return m_min_ctheta; }
    std::string    costhetabin(void) const;
    std::string    print(void) const;

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

#endif /* GLATLTCUBEMAP_HPP */
