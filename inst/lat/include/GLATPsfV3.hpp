/***************************************************************************
 *       GLATPsfV3.hpp  -  Fermi/LAT point spread function version 3       *
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
 * @file GLATPsfV3.hpp
 * @brief Fermi/LAT point spread function version 3 class definition
 * @author J. Knoedlseder
 */

#ifndef GLATPSFV3_HPP
#define GLATPSFV3_HPP

/* __ Includes ___________________________________________________________ */
#include <vector>
#include "GLATPsfBase.hpp"
#include "GFits.hpp"
#include "GFitsTable.hpp"
#include "GIntegrand.hpp"


/***********************************************************************//**
 * @class GLATPsfV3
 *
 * @brief Interface for the Fermi LAT point spread function version 3
 *
 * Version 3 of the Fermi/LAT PSF is the sum of two King model functions.
 * In contrast to version 2, this class interpolates the distributions
 * rather than the parameters.
 *
 * This class has been inspired by code from the Fermi/LAT ScienceTools.
 * For comparison check the file irfs/latResponse/src/Psf3.h
 ***************************************************************************/
class GLATPsfV3 : public GLATPsfBase {

public:
    // Constructors and destructors
    GLATPsfV3(void);
    GLATPsfV3(const GLATPsfV3& psf);
    virtual ~GLATPsfV3(void);

    // Operators
    GLATPsfV3& operator= (const GLATPsfV3& psf);

    // Methods
    void       clear(void);
    GLATPsfV3* clone(void) const;
    void       read(const GFitsTable* hdu);
    void       write(GFits& file) const;
    double     psf(const double& offset, const double& logE,
                   const double& ctheta);
    int        version(void) const { return 3; }

private:
    // Methods
    void          init_members(void);
    void          copy_members(const GLATPsfV3& psf);
    void          free_members(void);
    static double base_fct(const double& u, const double& gamma);
    static double base_int(const double& u, const double& gamma);
    double        eval_psf(const double& offset, const double& energy,
                           const int& index);

    // Protected members
    std::vector<double> m_ncore;        //!< PSF ncore parameter
    std::vector<double> m_ntail;        //!< PSF ntail parameter
    std::vector<double> m_score;        //!< PSF score parameter
    std::vector<double> m_stail;        //!< PSF stail parameter
    std::vector<double> m_gcore;        //!< PSF gcore parameter
    std::vector<double> m_gtail;        //!< PSF gtail parameter
};

#endif /* GLATPSFV3_HPP */
