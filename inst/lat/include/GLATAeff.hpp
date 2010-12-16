/***************************************************************************
 *                 GLATAeff.hpp  -  Fermi LAT effective area               *
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
 * @file GLATAeff.hpp
 * @brief Fermi LAT effective area class definition.
 * @author J. Knodlseder
 */

#ifndef GLATAEFF_HPP
#define GLATAEFF_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
//#include <vector>
#include <iostream>
#include "GLog.hpp"
#include "GLATPointing.hpp"
#include "GLATResponseTable.hpp"
#include "GFits.hpp"
#include "GFitsTable.hpp"
#include "GSkyDir.hpp"
#include "GEnergy.hpp"
#include "GTime.hpp"


/***********************************************************************//**
 * @class GLATAeff
 *
 * @brief Interface for the Fermi LAT effective area.
 ***************************************************************************/
class GLATAeff {

    // I/O friends
    friend std::ostream& operator<< (std::ostream& os, const GLATAeff& aeff);
    friend GLog&         operator<< (GLog& log, const GLATAeff& aeff);

public:
    // Constructors and destructors
    GLATAeff(void);
    GLATAeff(const std::string filename);
    GLATAeff(const GLATAeff& aeff);
    virtual ~GLATAeff(void);

    // Operators
    GLATAeff& operator= (const GLATAeff& aeff);
    double    operator() (const double& logE, const double& ctheta);
    double    operator() (const GSkyDir& srcDir, const GEnergy& srcEng,
                          const GTime& srcTime, const GLATPointing& pnt);

    // Methods
    void         clear(void);
    GLATAeff*    clone(void) const;
    void         load(const std::string filename);
    void         read(const GFits* file);
    int          size(void) const { return nenergies()*ncostheta(); }
    int          nenergies(void) const { return m_aeff_bins.nenergies(); }
    int          ncostheta(void) const { return m_aeff_bins.ncostheta(); }
    double       costhetamin(void) const { return m_min_ctheta; }
    void         costhetamin(const double& ctheta);
    void         ltcube_energy(const GEnergy& energy);
    double       ltcube_ctheta(const double& costheta);
    double       ltcube_ctheta_phi(const double& costheta, const double& phi);
    bool         hasphi(void) const { return false; }
    std::string  print(void) const;

private:
    // Methods
    void init_members(void);
    void copy_members(const GLATAeff& aeff);
    void free_members(void);
    void read_aeff(const GFitsTable* hdu);
    
    // Protected members
    GLATResponseTable m_aeff_bins;      //!< Aeff energy and cos theta binning
    double*           m_aeff;           //!< Aeff array
    double            m_min_ctheta;     //!< Minimum valid cos(theta)
    double            m_ltcube_logE;    //!< log10 energy for ltcube methods
};

#endif /* GLATAEFF_HPP */
