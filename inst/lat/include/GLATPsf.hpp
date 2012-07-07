/***************************************************************************
 *             GLATPsf.hpp  -  Fermi LAT point spread function             *
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
 * @file GLATPsf.hpp
 * @brief Fermi LAT point spread function class definition.
 * @author J. Knoedlseder
 */

#ifndef GLATPSF_HPP
#define GLATPSF_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include <iostream>
#include "GLog.hpp"
#include "GLATInstDir.hpp"
#include "GLATPointing.hpp"
#include "GLATResponseTable.hpp"
#include "GFits.hpp"
#include "GFitsTable.hpp"
#include "GSkyDir.hpp"
#include "GEnergy.hpp"
#include "GTime.hpp"
#include "GIntegrand.hpp"


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
    double   operator() (const double& offset, const double& logE,
                         const double& ctheta);
    double   operator() (const GLATInstDir& obsDir, const GSkyDir&  srcDir,
                         const GEnergy& srcEng, const GTime& srcTime,
                         const GLATPointing& pnt);

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
    double       costhetamin(void) const { return m_min_ctheta; }
    void         costhetamin(const double& ctheta);
    bool         hasphi(void) const { return false; }
    bool         isfront(void) const { return m_front; }
    bool         isback(void) const { return !m_front; }
    int          version(void) const { return m_version; }
    std::string  print(void) const;

private:
    // Methods
    void   init_members(void);
    void   copy_members(const GLATPsf& psf);
    void   free_members(void);
    void   read_scale(const GFitsTable* hdu);
    void   write_scale(GFits& file) const;
    double scale_factor(const double& energy) const;

    // PSF version 1 methods and classes
    void          read_psf_v1(const GFitsTable* hdu);
    void          write_psf_v1(GFits& file) const;
    double        psf_v1(const double& offset, const double& logE,
                         const double& ctheta);
    static double base_fct_v1(const double& u, const double& gamma);
    static double base_int_v1(const double& u, const double& gamma);
    class base_integrand_v1 : public GIntegrand {
    public:
        base_integrand_v1(double ncore, double ntail, double sigma,
                          double gcore, double gtail) :
                          m_ncore(ncore), m_ntail(ntail), m_sigma(sigma),
                          m_gcore(gcore), m_gtail(gtail) { }
        double eval(double x) {
            double r = x / m_sigma;
            double u = 0.5 * r * r;
            double f = m_ncore * base_fct_v1(u, m_gcore) +
                       m_ntail * base_fct_v1(u, m_gtail);
            return (f*sin(x));
        }
    private:
        double m_ncore;
        double m_ntail;
        double m_sigma;
        double m_gcore;
        double m_gtail;
    };


    // PSF version 3 methods and classes
    void          read_psf_v3(const GFitsTable* hdu);
    void          write_psf_v3(GFits& file) const;
    
    // Protected members
    int                 m_version;      //!< PSF version (starting from 1)
    bool                m_front;        //!< PSF is for front section?
    GLATResponseTable   m_rpsf_bins;    //!< PSF energy and cos theta binning
    double              m_scale_par1;   //!< PSF scaling parameter 1
    double              m_scale_par2;   //!< PSF scaling parameter 2
    double              m_scale_index;  //!< PSF scaling index
    double              m_min_ctheta;   //!< Minimum valid cos(theta)

    // PSF version 1 protected members
    std::vector<double> m_ncore;        //!< PSF ncore parameter
    std::vector<double> m_sigma;        //!< PSF sigma parameter
    std::vector<double> m_gcore;        //!< PSF gcore parameter
    std::vector<double> m_gtail;        //!< PSF gtail parameter

    // Additional PSF version 3 protected members
    //std::vector<double> m_ncore;        //!< PSF ncore parameter
    std::vector<double> m_ntail;        //!< PSF ntail parameter
    std::vector<double> m_score;        //!< PSF score parameter
    std::vector<double> m_stail;        //!< PSF stail parameter
    //std::vector<double> m_gcore;        //!< PSF gcore parameter
    //std::vector<double> m_gtail;        //!< PSF gtail parameter
};

#endif /* GLATPSF_HPP */
