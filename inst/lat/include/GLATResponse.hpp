/***************************************************************************
 *               GLATResponse.hpp  -  GLAST LAT Response class             *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2008-2010 by Jurgen Knodlseder                           *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GLATResponse.hpp
 * @brief GLATResponse class interface definition.
 * @author J. Knodlseder
 */

#ifndef GLATRESPONSE_HPP
#define GLATRESPONSE_HPP

/* __ Includes ___________________________________________________________ */
#include <vector>
#include "GLATAeff.hpp"
#include "GLATResponseTable.hpp"
#include "GLATEventAtom.hpp"
#include "GLATEventBin.hpp"
#include "GEvent.hpp"
#include "GResponse.hpp"
#include "GPointing.hpp"
#include "GInstDir.hpp"
#include "GRoi.hpp"
#include "GGti.hpp"
#include "GEbounds.hpp"
#include "GSkyDir.hpp"
#include "GEnergy.hpp"
#include "GTime.hpp"
#include "GNodeArray.hpp"
#include "GVector.hpp"
#include "GModel.hpp"
#include "GFits.hpp"
#include "GFitsTable.hpp"


/***********************************************************************//**
 * @class GLATResponse
 *
 * @brief Interface for the Fermi LAT instrument response function.
 ***************************************************************************/
class GLATResponse : public GResponse {

public:
    // Constructors and destructors
    GLATResponse(void);
    GLATResponse(const GLATResponse& rsp);
    virtual ~GLATResponse(void);

    // Operators
    GLATResponse& operator= (const GLATResponse & rsp);

    // Implement pure virtual base class methods
    void          clear(void);
    GLATResponse* clone(void) const;
    void          load(const std::string& rspname);
    bool          hasedisp(void) const { return false; }
    bool          hastdisp(void) const { return false; }
    std::string   print(void) const;
    
    // Implemented response computation methods
    double irf(const GInstDir& obsDir, const GEnergy& obsEng, const GTime& obsTime,
               const GSkyDir&  srcDir, const GEnergy& srcEng, const GTime& srcTime,
               const GPointing& pnt) const;
    double nirf(const GSkyDir&  srcDir, const GEnergy& srcEng, const GTime& srcTime,
                const GPointing& pnt, const GRoi& roi, const GEbounds& ebds,
                const GGti& gti) const;


    double diffrsp(const GEvent& event, const GModel& model,
                   const GEnergy& srcEng, const GTime& srcTime,
                   const GPointing& pnt) const;
    double live(const GSkyDir&  srcDir, const GEnergy& srcEng,
                const GTime& srcTime, const GPointing& pnt) const;
    //double aeff(const GSkyDir&  srcDir, const GEnergy& srcEng, const GTime& srcTime,
    //            const GPointing& pnt) const;
    double psf(const GInstDir& obsDir,
               const GSkyDir&  srcDir, const GEnergy& srcEng, const GTime& srcTime,
               const GPointing& pnt) const;
    double edisp(const GEnergy& obsEng,
                 const GSkyDir&  srcDir, const GEnergy& srcEng, const GTime& srcTime,
                 const GPointing& pnt) const;
    double tdisp(const GTime& obsTime,
                 const GSkyDir&  srcDir, const GEnergy& srcEng, const GTime& srcTime,
                 const GPointing& pnt) const;
    double npsf(const GSkyDir&  srcDir, const GEnergy& srcEng,
                const GTime& srcTime, const GPointing& pnt, const GRoi& roi) const;
    double nedisp(const GSkyDir&  srcDir, const GEnergy& srcEng,
                  const GTime& srcTime, const GPointing& pnt, const GEbounds& ebds) const;
    double ntdisp(const GSkyDir&  srcDir, const GEnergy& srcEng,
                  const GTime& srcTime, const GPointing& pnt, const GGti& gti) const;

    // Other Methods
    int       size(void) const { return m_aeff.size(); }
    GLATAeff* aeff(const int& index) const;
    void      save(const std::string& rspname) const;
    
    //double  aeff(const double& logE, const double& ctheta);
    //void    aeff_ctheta_min(const double& ctheta);
    //double  aeff_ctheta_min(void) const;
    double  psf(const double& delta, const double& logE, const double& ctheta);
    GVector psf(const GVector& delta, const double& logE, const double& ctheta);

private:
    // Private PSF methods
    void    psf_init(std::string section);
    void    psf_append(GFits& file) const;
    double  psf_scale(const double& energy) const;
    double  psf_scale_log(const double& logE) const;
    double  psf_base_value(const double& u, const double& gamma) const;
    GVector psf_base_vector(const GVector& u, const double& gamma) const;
    void    psf_init_members(void);
    void    psf_copy_members(const GLATResponse& rsp);
    void    psf_free_members(void);

    // Private energy dissipation methods
    void    edisp_init(std::string section);
    void    edisp_append(GFits& file) const;
    void    edisp_init_members(void);
    void    edisp_copy_members(const GLATResponse& rsp);
    void    edisp_free_members(void);

    // Other private methods
    void    init_members(void);
    void    copy_members(const GLATResponse& rsp);
    void    free_members(void);
    GVector get_fits_vector(const GFitsTable* hdu, const std::string& colname, int row = 0);
    double  diffrsp_atom(const GLATEventAtom& event, const GModel& model,
                         const GEnergy& srcEng, const GTime& srcTime,
                         const GPointing& pnt) const;
    double  diffrsp_bin(const GLATEventBin& event, const GModel& model,
                        const GEnergy& srcEng, const GTime& srcTime,
                        const GPointing& pnt) const;

    // Private PSF data
    GLATResponseTable m_psf_bins;        //!< PSF energy and cos theta binning
    GNodeArray        m_psf_angle;       //!< PSF vector binning
    double            m_psf_angle_max;   //!< PSF maximum angular distance covered
    double*           m_psf;             //!< PSF vector array
    double*           m_norm;            //!< PSF normalization array
    double*           m_sigma;           //!< PSF sigma array

    // Private Edisp data
    GLATResponseTable m_edisp_bins;      //!< Edisp energy and cos theta binning

    // Private members
    bool                   m_hasfront;    //!< Front IRF loaded
    bool                   m_hasback;     //!< Back IRG loaded
    std::vector<GLATAeff*> m_aeff;        //!< Effective areas
    std::vector<GLATAeff*> m_psf1;         //!< Point spread functions
    std::vector<GLATAeff*> m_edisp;       //!< Energy dispersions
};

#endif /* GLATRESPONSE_HPP */
