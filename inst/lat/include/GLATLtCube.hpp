/***************************************************************************
 *                GLATLtCube.hpp  -  Fermi LAT lifetime cube               *
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
 * @file GLATLtCube.hpp
 * @brief GLATLtCube class definition.
 * @author J. Knodlseder
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
 * @brief Interface for the Fermi LAT lifetime cube.
 *
 * The livetime cube holds the lifetime as function and zenith and azimuth
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
                            const double& offset, const GLATPsf& psf);

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
    bool          m_livetime_correct;
};

#endif /* GLATLTCUBE_HPP */
