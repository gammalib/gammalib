/***************************************************************************
 *                  GCTASupport.hpp - CTA support functions                *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2011-2018 by Juergen Knoedlseder                         *
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
 * @file GCTASupport.hpp
 * @brief Definition of support function used by CTA classes
 * @author Juergen Knoedlseder
 */

#ifndef GCTASUPPORT_HPP
#define GCTASUPPORT_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include <vector>

/* __ Namespaces _________________________________________________________ */

/* __ Constants __________________________________________________________ */

/* __ Forward declarations _______________________________________________ */
class GFits;
class GFitsHDU;
class GCTARoi;
class GEbounds;
class GPhases;
class GObservation;
class GResponse;
class GEvent;
class GCTAObservation;
class GCTAPointing;
class GCTAResponse;
class GCTAResponseIrf;
class GCTAResponseCube;
class GCTAAeff;
class GCTABackground;
class GCTAEventList;
class GCTAInstDir;

/* __ Prototypes _________________________________________________________ */
namespace gammalib {
    const GCTAObservation&  cta_obs(const std::string&  origin,
                                    const GObservation& obs);
    const GCTAPointing&     cta_pnt(const std::string&  origin,
                                    const GObservation& obs);
    const GCTAResponse*     cta_rsp(const std::string& origin,
                                    const GResponse&   rsp);
    const GCTAResponseIrf&  cta_rsp_irf(const std::string&  origin,
                                        const GObservation& obs);
    const GCTAResponseCube& cta_rsp_cube(const std::string&  origin,
                                         const GObservation& obs);
    const GCTAAeff&         cta_rsp_aeff(const std::string&  origin,
                                         const GObservation& obs);
    const GCTABackground&   cta_rsp_bkg(const std::string&  origin,
                                        const GObservation& obs);
    const GCTAEventList&    cta_event_list(const std::string&  origin,
                                           const GObservation& obs);
    const GCTAInstDir&      cta_dir(const std::string&  origin,
                                    const GEvent&       event);
    double                  cta_roi_arclength(const double& rad,
                                              const double& dist,
                                              const double& cosdist,
                                              const double& sindist,
                                              const double& roi,
                                              const double& cosroi);
    GCTARoi                 read_ds_roi(const GFitsHDU& hdu);
    GEbounds                read_ds_ebounds(const GFitsHDU& hdu);
    GPhases                 read_ds_phase(const GFitsHDU& hdu);
    std::string             read_ds_gti_extname(const GFitsHDU& hdu);
    std::string             gadf_hduclas4(const GFits&       fits,
                                          const std::string& hduclas4);
}

#endif /* GCTASUPPORT_HPP */
