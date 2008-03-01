/***************************************************************************
 *               GLATResponse.hpp  -  GLAST LAT Response class             *
 * ----------------------------------------------------------------------- *
 *  copyright            : (C) 2008 by Jurgen Knodlseder                   *
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
 * @brief GLATResponse class definition.
 * @author J. Knodlseder
 */

#ifndef GLATRESPONSE_HPP
#define GLATRESPONSE_HPP

/* __ Includes ___________________________________________________________ */
#include "GVector.hpp"
#include "GSkyDir.hpp"
#include "GResponse.hpp"
#include "GFitsHDU.hpp"

/* __ Namespaces _________________________________________________________ */


/***************************************************************************
 *                        GLATResponse class definition                    *
 ***************************************************************************/
/**
 * @class GLATResponse
 *
 * @brief Interface for the GLAST LAT instrument response function classes.
 *
 * @author J. Knodlseder
 */
class GLATResponse : public GResponse {

public:
    // Constructors and destructors
    GLATResponse();
    GLATResponse(const GLATResponse& rsp);
    ~GLATResponse();

    // Operators
    GLATResponse& operator= (const GLATResponse & rsp);

    // Methods
    double irf(const GSkyDir& obsDir, const double& obsEng,
               const GSkyDir& srcDir, const double& srcEng,
               const GSkyDir& instPntDir, const double& instPosAng,
               const double& time);
    double psf(const GSkyDir& obsDir, const double& obsEng,
               const GSkyDir& srcDir, const double& srcEng,
               const GSkyDir& instPntDir, const double& instPosAng,
               const double& time);
    double aeff(const GSkyDir& obsDir, const double& obsEng,
                const GSkyDir& srcDir, const double& srcEng,
                const GSkyDir& instPntDir, const double& instPosAng,
                const double& time);
    double edisp(const GSkyDir& obsDir, const double& obsEng,
                 const GSkyDir& srcDir, const double& srcEng,
                 const GSkyDir& instPntDir, const double& instPosAng,
                 const double& time);

    void set_caldb(std::string caldb);
    void load(std::string rspname, std::string rsptype);
    void save(std::string rspname) const;
    GLATResponse* clone(void) const;

private:
    // Private Effective Area methods
    void aeff_init(void);

    // Private PSF methods
    void    psf_init(void);
    double  psf_get(const double& delta, const double& logE, const double& ctheta);
    double  psf_scale(const double& energy);
    double  psf_base_value(const double& u, const double& gamma);
    GVector psf_base_vector(const GVector& u, const double& gamma);
    
    // Private PSF data
    int     m_psf_angle_num;
    double  m_psf_angle_min;
    double  m_psf_angle_max;
    double  m_psf_angle_bin;
    int     m_psf_energy_num;
    int     m_psf_ctheta_num;
    double* m_psf;
    double* m_norm;
    double* m_sigma;

    // Provate energy dissipation methods
    void edisp_init(void);

    // Private methods
    void    init_members(void);
    void    copy_members(const GLATResponse& rsp);
    void    free_members(void);
    GVector get_fits_vector(const GFitsHDU* hdu, const std::string& colname, int row = 0);

    // Private data area
    std::string m_rsptype;
    int         m_section;   // 0=front, 1=back

};


/***************************************************************************
 *                              Inline methods                             *
 ***************************************************************************/
inline 
GLATResponse* GLATResponse::clone(void) const
{
    return new GLATResponse(*this);
}

#endif /* GLATRESPONSE_HPP */
