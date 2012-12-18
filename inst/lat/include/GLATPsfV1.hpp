/***************************************************************************
 *       GLATPsfV1.hpp  -  Fermi:LAT point spread function version 1       *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2012 by Juergen Knoedlseder                              *
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
 * @file GLATPsfV1.hpp
 * @brief Fermi:LAT point spread function version 1 class definition
 * @author Juergen Knoedlseder
 */

#ifndef GLATPSFV1_HPP
#define GLATPSFV1_HPP

/* __ Includes ___________________________________________________________ */
#include <vector>
#include <cmath>
#include "GLATPsfBase.hpp"
#include "GFits.hpp"
#include "GFitsTable.hpp"
#include "GFunction.hpp"


/***********************************************************************//**
 * @class GLATPsfV1
 *
 * @brief Interface for the Fermi LAT point spread function version 1
 *
 * This class has been inspired by code from the Fermi/LAT ScienceTools.
 * For comparison check the file irfs/latResponse/src/Psf.h
 ***************************************************************************/
class GLATPsfV1 : public GLATPsfBase {

public:
    // Constructors and destructors
    GLATPsfV1(void);
    GLATPsfV1(const GLATPsfV1& psf);
    virtual ~GLATPsfV1(void);

    // Operators
    GLATPsfV1& operator= (const GLATPsfV1& psf);

    // Methods
    void        clear(void);
    GLATPsfV1*  clone(void) const;
    void        read(const GFitsTable* hdu);
    void        write(GFits& file) const;
    double      psf(const double& offset, const double& logE,
                    const double& ctheta);
    int         version(void) const { return 1; }
    std::string print(void) const;

private:
    // Methods
    void          init_members(void);
    void          copy_members(const GLATPsfV1& psf);
    void          free_members(void);
    static double base_fct(const double& u, const double& gamma);
    static double base_int(const double& u, const double& gamma);

    // Integrand class. This class is used to perform the radial
    // integration of the PSF that assures the proper normalization
    // of the PSF at low energies.
    class base_integrand : public GFunction {
    public:
        base_integrand(double ncore, double ntail, double sigma,
                       double gcore, double gtail) :
                       m_ncore(ncore), m_ntail(ntail), m_sigma(sigma),
                       m_gcore(gcore), m_gtail(gtail) { }
        double eval(double x) {
            double r = x / m_sigma;
            double u = 0.5 * r * r;
            double f = m_ncore * base_fct(u, m_gcore) +
                       m_ntail * base_fct(u, m_gtail);
            return (f*std::sin(x));
        }
    private:
        double m_ncore;
        double m_ntail;
        double m_sigma;
        double m_gcore;
        double m_gtail;
    };

    // Protected members
    std::vector<double> m_ncore;        //!< PSF ncore parameter
    std::vector<double> m_sigma;        //!< PSF sigma parameter
    std::vector<double> m_gcore;        //!< PSF gcore parameter
    std::vector<double> m_gtail;        //!< PSF gtail parameter
};

#endif /* GLATPSFV1_HPP */
