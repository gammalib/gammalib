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
    void load(std::string rspname);
    void save(std::string rspname) const;
    GLATResponse* clone(void) const;

private:
    // Private methods
    void    init_members(void);
    void    copy_members(const GLATResponse& rsp);
    void    free_members(void);
    void    init_aeff(void);
    void    init_psf(void);
    void    init_edisp(void);
    GVector get_fits_vector(const GFitsHDU* hdu, const std::string& colname, int row = 0);

    // Private data area

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
