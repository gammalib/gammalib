/***************************************************************************
 *             GLATPsf.hpp  -  Fermi LAT point spread function             *
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
 * @file GLATPsf.hpp
 * @brief Fermi LAT point spread function class definition.
 * @author J. Knodlseder
 */

#ifndef GLATPSF_HPP
#define GLATPSF_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
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
 * @class GLATPsf
 *
 * @brief Interface for the Fermi LAT point spread function.
 ***************************************************************************/
class GLATPsf {

    // I/O friends
    friend std::ostream& operator<< (std::ostream& os, const GLATPsf& psf);
    friend GLog&         operator<< (GLog& log, const GLATPsf& psf);

public:
    // Constructors and destructors
    GLATPsf(void);
    GLATPsf(const std::string filename);
    GLATPsf(const GLATPsf& psf);
    virtual ~GLATPsf(void);

    // Operators
    GLATPsf& operator= (const GLATPsf& psf);
    //double    operator() (const double& logE, const double& ctheta);
    //double    operator() (const GSkyDir& srcDir, const GEnergy& srcEng,
    //                      const GTime& srcTime, const GLATPointing& pnt);

    // Methods
    void         clear(void);
    GLATPsf*     clone(void) const;
    void         load(const std::string filename);
    void         save(const std::string filename, bool clobber = false);
    void         read(const GFits& file);
    void         write(GFits& file) const;
    int          size(void) const { return nenergies()*ncostheta(); }
    int          nenergies(void) const { return m_rpsf_bins.nenergies(); }
    int          ncostheta(void) const { return m_rpsf_bins.ncostheta(); }
    //double       costhetamin(void) const { return m_min_ctheta; }
    //void         costhetamin(const double& ctheta);
    //void         ltcube_energy(const GEnergy& energy);
    //double       ltcube_ctheta(const double& costheta);
    //double       ltcube_ctheta_phi(const double& costheta, const double& phi);
    bool         hasphi(void) const { return false; }
    std::string  print(void) const;

private:
    // Methods
    void init_members(void);
    void copy_members(const GLATPsf& psf);
    void free_members(void);
    void read_psf(const GFitsTable* hdu);
    void write_psf(GFits& file) const;
    
    // Protected members
    GLATResponseTable   m_rpsf_bins;    //!< PSF energy and cos theta binning
    std::vector<double> m_ncore;        //!< PSF ncore parameter
    std::vector<double> m_sigma;        //!< PSF sigma parameter
    std::vector<double> m_gcore;        //!< PSF gcore parameter
    std::vector<double> m_gtail;        //!< PSF gtail parameter
    std::vector<double> m_scale;        //!< PSF scaling parameters
    double              m_ltcube_logE;  //!< log10 energy for ltcube methods
};

#endif /* GLATPSF_HPP */
