/***************************************************************************
 *                  GCTASupport.cpp - CTA support functions                *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2011-2020 by Juergen Knoedlseder                         *
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
 * @file GCTASupport.cpp
 * @brief Implementation of support function used by CTA classes
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <typeinfo>
#include <cmath>
#include <iostream>
#include <string>
#include "GCTASupport.hpp"
#include "GTools.hpp"
#include "GMath.hpp"
#include "GFitsHDU.hpp"
#include "GEbounds.hpp"
#include "GPhases.hpp"
#include "GResponse.hpp"
#include "GModelData.hpp"
#include "GModelSpectral.hpp"
#include "GModelTemporal.hpp"
#include "GCTAObservation.hpp"
#include "GCTAPointing.hpp"
#include "GCTAResponse.hpp"
#include "GCTAResponseIrf.hpp"
#include "GCTAResponseCube.hpp"
#include "GCTAAeff.hpp"
#include "GCTABackground.hpp"
#include "GCTAEventList.hpp"
#include "GCTAEventCube.hpp"
#include "GCTARoi.hpp"
#include "GCTAInstDir.hpp"
#include "GCTAModelBackground.hpp"
#include "GCTAModelAeffBackground.hpp"
#include "GCTAModelIrfBackground.hpp"
#include "GCTAModelCubeBackground.hpp"
#include "GCTAModelRadialAcceptance.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_READ_DS_ROI                      "gammalib::read_ds_roi(GFitsHDU&)"
#define G_READ_DS_EBOUNDS              "gammalib::read_ds_ebounds(GFitsHDU&)"
#define G_READ_DS_PHASE                  "gammalib::read_ds_phase(GFitsHDU&)"
#define G_READ_DS_GTI              "gammalib::read_ds_gti_extname(GFitsHDU&)"
#define G_CTA_MODEL_NROI           "gammalib::cta_model_nroi(GObservation*, "\
                                             "GModelData*, GEnergy&, GTime&)"

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */


/***********************************************************************//**
 * @brief Extract ROI from data sub-space keywords
 *
 * @param[in] hdu FITS HDU
 *
 * @exception GException::invalid_value
 *            Invalid ROI data sub-space encountered
 *
 * Reads the ROI data sub-space keywords by searching for a DSTYPx keyword
 * named "POS(RA,DEC)". The data sub-space information is expected to be in
 * the format "CIRCLE(267.0208,-24.78,4.5)", where the 3 arguments are Right
 * Ascension, Declination and radius in units of degrees. No detailed syntax
 * checking is performed.
 *
 * If no ROI information has been found, an GCTARoi object with initial
 * values will be returned. 
 ***************************************************************************/
GCTARoi gammalib::read_ds_roi(const GFitsHDU& hdu)
{
    // Initialise ROI
    GCTARoi roi;

    // Get number of data sub-space keywords (default to 0 if keyword is
    // not found)
    int ndskeys = (hdu.has_card("NDSKEYS")) ? hdu.integer("NDSKEYS") : 0;

    // Loop over all data selection keys
    for (int i = 1; i <= ndskeys; ++i) {

        // Set data sub-space key strings
        std::string type_key  = "DSTYP"+gammalib::str(i);
        //std::string unit_key  = "DSUNI"+gammalib::str(i);
        std::string value_key = "DSVAL"+gammalib::str(i);

        // Continue only if type_key is found and if this key is POS(RA,DEC)
        if (hdu.has_card(type_key) && hdu.string(type_key) == "POS(RA,DEC)") {

            // ...
            //std::string unit              = gammalib::toupper(hdu.string(unit_key));
            std::string value             = hdu.string(value_key);
            std::string value_proc        = gammalib::strip_chars(value, "CIRCLE(");
            value_proc                    = gammalib::strip_chars(value_proc, ")");
            std::vector<std::string> args = gammalib::split(value_proc, ",");
            if (args.size() == 3) {

                // Extract sky direction and radius
                double ra  = gammalib::todouble(args[0]);
                double dec = gammalib::todouble(args[1]);
                double rad = gammalib::todouble(args[2]);
                GSkyDir     skydir;
                skydir.radec_deg(ra, dec);

                // If the header contains pointing information then set a
                // full instrument direction, including DETX and DETY
                // coordinates.
                if (hdu.has_card("RA_PNT") && hdu.has_card("DEC_PNT")) {
                    double ra_pnt  = hdu.real("RA_PNT");
                    double dec_pnt = hdu.real("DEC_PNT");
                    GSkyDir dir_pnt;
                    dir_pnt.radec_deg(ra_pnt, dec_pnt);
                    GCTAPointing pnt(dir_pnt);
                    GCTAInstDir  dir = pnt.instdir(skydir);
                    roi.centre(dir);
                }

                // ... otherwise just set the sky direction
                else {
                    GCTAInstDir dir(skydir);
                    roi.centre(dir);
                }

                // Set RoI radius
                roi.radius(rad);
            }
            else {
                std::string msg = "Invalid acceptance cone value \""+value+
                                  "\" encountered in data sub-space "
                                  "key \""+value_key+"\".";
                throw GException::invalid_value(G_READ_DS_ROI, msg);
            }

        } // endif: POS(RA,DEC) type found

    } // endfor: looped over data sub-space keys

    // Return roi
    return roi;
}


/***********************************************************************//**
 * @brief Read energy boundary data sub-space keywords
 *
 * @param[in] hdu FITS HDU
 *
 * @exception GException::invalid_value
 *            Invalid energy data sub-space encountered
 *
 * Reads the energy boundary data sub-space keywords by searching for a 
 * DSTYPx keyword named "ENERGY". The data sub-space information is expected
 * to be in the format "200:50000", where the 2 arguments are the minimum
 * and maximum energy. The energy unit is given by the keyword DSUNIx, which
 * supports keV, MeV, GeV and TeV (case independent). No detailed syntax
 * checking is performed.
 ***************************************************************************/
GEbounds gammalib::read_ds_ebounds(const GFitsHDU& hdu)
{
    // Initialise energy boundaries
    GEbounds ebounds;

    // Get number of data sub-space keywords (default to 0 if keyword is
    // not found)
    int ndskeys = (hdu.has_card("NDSKEYS")) ? hdu.integer("NDSKEYS") : 0;

    // Loop over all data sub-space keys
    for (int i = 1; i <= ndskeys; ++i) {

        // Set data sub-space key strings
        std::string type_key  = "DSTYP"+gammalib::str(i);
        std::string unit_key  = "DSUNI"+gammalib::str(i);
        std::string value_key = "DSVAL"+gammalib::str(i);

        // Continue only if type_key is found and if this key is ENERGY
        if (hdu.has_card(type_key) && hdu.string(type_key) == "ENERGY") {

            // Extract energy boundaries
            std::string value               = hdu.string(value_key);
            std::string unit                = hdu.string(unit_key);
            std::vector<std::string> values = gammalib::split(value, ":");
            if (values.size() == 2) {
                double  emin = gammalib::todouble(values[0]);
                double  emax = gammalib::todouble(values[1]);
                GEnergy e_min(emin, unit);
                GEnergy e_max(emax, unit);
                ebounds.append(e_min, e_max);
            }
            else {
                std::string msg = "Invalid energy value \""+value+
                                  "\" encountered in data selection key \""+
                                  value_key+"\"";
                throw GException::invalid_value(G_READ_DS_EBOUNDS, msg);
            }

        } // endif: ENERGY type_key found

    } // endfor: looped over data selection keys

    // Return
    return ebounds;
}


/***********************************************************************//**
 * @brief Read phase boundary data sub-space keywords
 *
 * @param[in] hdu FITS HDU
 * @return Phase intervals
 *
 * @exception GException::invalid_value
 *            Invalid phase data sub-space encountered
 *
 * Reads the phase boundary data sub-space keywords by searching for a 
 * DSTYPx keyword named "PHASE". The data sub-space information is expected
 * to be in the format "0.1:0.3,0.5:0.7", where each subset of numbers 
 * represents the minimum and maximum phase. No detailed syntax
 * checking is performed.
 ***************************************************************************/
GPhases gammalib::read_ds_phase(const GFitsHDU& hdu)
{
    // Initialise phase intervals
    GPhases phases;

    // Get number of data sub-space keywords (default to 0 if keyword is
    // not found)
    int ndskeys = (hdu.has_card("NDSKEYS")) ? hdu.integer("NDSKEYS") : 0;

    // Loop over all data sub-space keys
    for (int i = 1; i <= ndskeys; ++i) {

        // Set data sub-space key strings
        std::string type_key  = "DSTYP"+gammalib::str(i);
        //std::string unit_key  = "DSUNI"+gammalib::str(i);
        std::string value_key = "DSVAL"+gammalib::str(i);

        // Continue only if type_key is found and if this key is PHASE
        if (hdu.has_card(type_key) && hdu.string(type_key) == "PHASE") {

            // Extract phase boundaries
            std::string value                  = hdu.string(value_key);
            std::vector<std::string> intervals = gammalib::split(value, ",");
            for (int j = 0; j < intervals.size(); ++j) {
                std::vector<std::string> values = gammalib::split(intervals[j], ":");
                if (values.size() == 2) {
                    double phmin = gammalib::todouble(values[0]);
                    double phmax = gammalib::todouble(values[1]);
                    phases.append(phmin, phmax);
                }
                else {
                    std::string msg = "Invalid phase value \""+value+
                                      "\" encountered in data selection key \""+
                                      value_key+"\"";
                    throw GException::invalid_value(G_READ_DS_PHASE, msg);
                }

            } //endfor: looped over phase intervals

        } // endif: PHASE type_key found

    } // endfor: looped over data selection keys

    // Return
    return phases;
} 


/***********************************************************************//**
 * @brief Return Good Time Intervals extension name from data sub-space
 *        keywords
 *
 * @param[in] hdu FITS HDU
 * @return Good Time Interval extension name
 *
 * @exception GException::invalid_value
 *            Invalid Good Time Intervals data sub-space encountered
 *
 * Returns the name of the FITS extension that contains the Good Time
 * Intervals by screening the data sub-space keywords that are present in
 * the FITS header. The method searches for a DSTYPx keyword named "TIME"
 * and a corresponding DSVALx keyword named "TABLE", and the extension name
 * is extracted from the corresponding DSREFx keyword. Note that by
 * convention the extension name is preceeded by a colon, which is stripped
 * by this method.
 ***************************************************************************/
std::string gammalib::read_ds_gti_extname(const GFitsHDU& hdu)
{
    // Initialise extension name
    std::string extname;

    // Get number of data sub-space keys (default to 0 if "NDSKEYS" keyword
    // is not found)
    int ndskeys = (hdu.has_card("NDSKEYS")) ? hdu.integer("NDSKEYS") : 0;

    // Loop over all data sub-space keys
    for (int i = 1; i <= ndskeys; ++i) {

        // Set data sub-space key strings
        std::string type_key  = "DSTYP"+gammalib::str(i);
        //std::string unit_key  = "DSUNI"+gammalib::str(i);
        std::string value_key = "DSVAL"+gammalib::str(i);
        std::string ref_key   = "DSREF"+gammalib::str(i);

        // Skip if DSTYPi keyword does not exist, or if it is not "TIME"
        if (!hdu.has_card(type_key) || hdu.string(type_key) != "TIME") {
            continue;
        }

        // If DSVALi keyword does not exist or if it's value is not
        // "TABLE" then throw an exception
        if (!hdu.has_card(value_key)) {
            std::string msg = "Keyword \""+value_key+"\" missing for data "
                              "sub-space selection of type "+type_key+
                              "=\"TIME\". Please correct FITS header.";
            throw GException::invalid_value(G_READ_DS_GTI, msg);
        }
        if (hdu.string(value_key) != "TABLE") {
            std::string msg = "Cannot interpret keyword \""+value_key+"\" "
                              "value \""+hdu.string(value_key)+"\". Only "
                              "the value \"TABLE\" is supported.";
            throw GException::invalid_value(G_READ_DS_GTI, msg);
        }

        // If DSREFi keyword does not exist then throw an exception
        if (!hdu.has_card(ref_key)) {
            std::string msg = "Keyword \""+ref_key+"\" missing for data "
                              "sub-space selection of type "+type_key+
                              "=\"TIME\". Please correct FITS header.";
            throw GException::invalid_value(G_READ_DS_GTI, msg);
        }

        // Get extension name (strip any leading and trailing colons)
        extname = gammalib::strip_whitespace(gammalib::strip_chars(hdu.string(ref_key), ":"));

    } // endfor: looped over data sub-space keywords

    // Return
    return extname;
}


/***********************************************************************//**
 * @brief Return extension name for GADF response table of given HDU class 4
 *
 * @param[in] fits FITS file.
 * @param[in] hduclas4 HDU class 4.
 * @return Extension name.
 *
 * Returns the extension name for GADF response table of given HDU class 4.
 * If the response table is not found, an empty extension name is returned.
 ***************************************************************************/
std::string gammalib::gadf_hduclas4(const GFits&       fits,
                                    const std::string& hduclas4)
{
    // Initialise extension name
    std::string extname;

    // Loop over all FITS HDUs
    for (int i = 0; i < fits.size(); ++i) {

        // Get FITS HDU
        const GFitsHDU* hdu = fits[i];

        // Skip HDU if it has no "HDUCLASS", "HDUCLAS1" and "HDUCLAS4" keywords
        if (!(hdu->has_card("HDUCLASS") &&
              hdu->has_card("HDUCLAS1") &&
              hdu->has_card("HDUCLAS4"))) {
            continue;
        }

        // Skip HDU if "HDUCLASS" is not "GADF", "HDUCLAS1" is not response
        // and "HDUCLAS4" is not hduclas4
        if (!((gammalib::strip_whitespace(hdu->string("HDUCLASS")) == "GADF") &&
              (gammalib::strip_whitespace(hdu->string("HDUCLAS1")) == "RESPONSE") &&
              (gammalib::strip_whitespace(hdu->string("HDUCLAS4")) == hduclas4))) {
            continue;
        }

        // Extract extension name and break
        extname = hdu->extname();

    } // endfor: loop over headers

    // Return
    return extname;
}


/***********************************************************************//**
 * @brief Determine number of radial Romberg iterations
 *
 * @param[in] rho_max Maximum radial offset (radians).
 * @param[in] resolution Requested angular resolution (radians).
 * @param[in] iter_min Minimum number of iterations.
 * @param[in] iter_max Maximum number of iterations.
 * @return Number of radial Romberg iterations.
 *
 * Determines the number of radial Romberg iterations using the formula
 *
 * \f[
 *    iter = \log_2 \left( \frac{\rho_{\rm max}}{resolution} \right) + 1
 * \f]
 *
 * where
 * \f$\rho_{\rm max}\f$ is the maximum radial offset and
 * \f$resolution\f$ is the required angular resolution.
 *
 * The result will be constrained to the interval [@p iter_min,@p iter_max].
 ***************************************************************************/
int gammalib::iter_rho(const double& rho_max,
                       const double& resolution,
                       const int&    iter_min,
                       const int&    iter_max)
{
    // Initialise number of iterations
    int iter = iter_min;

    // Continue only if radial offset and resolution are positive
    if ((rho_max > 0.0) && (resolution > 0.0)) {
        double arg = rho_max / resolution;
        iter       = int(std::log(arg) * gammalib::inv_ln2) + 1;
        if (iter < iter_min) {
            iter = iter_min;
        }
        else if (iter > iter_max) {
            iter = iter_max;
        }
    }

    // Return number of iterations
    return iter;
}


/***********************************************************************//**
 * @brief Determine number of azimuthal Romberg iterations
 *
 * @param[in] rho Radial offset (radians).
 * @param[in] resolution Requested angular resolution (radians).
 * @param[in] iter_min Minimum number of iterations.
 * @param[in] iter_max Maximum number of iterations.
 * @return Number of azimuthal Romberg iterations.
 *
 * Determines the number of azimuthal Romberg iterations using the formula
 *
 * \f[
 *    iter = \log_2 \left( \frac{2\pi \rho}{resolution} \right) + 1
 * \f]
 *
 * where
 * \f$\rho\f$ is the radial offset and
 * \f$resolution\f$ is the required angular resolution.
 *
 * The result will be constrained to the interval [@p iter_min,@p iter_max].
 ***************************************************************************/
int gammalib::iter_phi(const double& rho,
                       const double& resolution,
                       const int&    iter_min,
                       const int&    iter_max)
{
    // Initialise number of iterations
    int iter = iter_min;

    // Continue only if radial offset and resolution are positive
    if ((rho > 0.0) && (resolution > 0.0)) {
        double arg = gammalib::twopi * rho / resolution;
        iter       = int(std::log(arg) * gammalib::inv_ln2) + 1;
        if (iter < iter_min) {
            iter = iter_min;
        }
        else if (iter > iter_max) {
            iter = iter_max;
        }
    }

    // Return number of iterations
    return iter;
}


/***********************************************************************//**
 * @brief Retrieve CTA observation from generic observation
 *
 * @param[in] origin Method asking for pointer retrieval.
 * @param[in] obs Generic observation.
 * @return Reference to CTA observation.
 *
 * @exception GException::invalid_argument
 *            Observation @p obs is not a CTA observations.
 *
 * Dynamically casts generic observation into a CTA observation. If the
 * generic observation is not a CTA observation, an exception is thrown.
 ***************************************************************************/
const GCTAObservation& gammalib::cta_obs(const std::string&  origin,
                                         const GObservation& obs)
{
    // Get pointer on CTA observation
    const GCTAObservation* cta = dynamic_cast<const GCTAObservation*>(&obs);

    // If pointer is not valid then throw an exception
    if (cta == NULL) {
        std::string cls = std::string(typeid(&obs).name());
        std::string msg = "Invalid observation type \""+cls+"\" provided with "
                          "name \""+obs.name()+"\" (ID="+obs.id()+"). "
                          "Please specify a \"GCTAObservation\" instance as "
                          "argument.";
        throw GException::invalid_argument(origin, msg);
    }

    // Return reference
    return (*cta);
}


/***********************************************************************//**
 * @brief Retrieve CTA pointing from generic observation
 *
 * @param[in] origin Method asking for pointer retrieval.
 * @param[in] obs Generic observation.
 * @return Reference to CTA pointing.
 *
 * Extract CTA pointing from a CTA observation.
 ***************************************************************************/
const GCTAPointing& gammalib::cta_pnt(const std::string&  origin,
                                      const GObservation& obs)
{
    // Retrieve CTA observation
    const GCTAObservation& cta = gammalib::cta_obs(origin, obs);

    // Return CTA pointing
    return (cta.pointing());
}


/***********************************************************************//**
 * @brief Retrieve CTA response from generic observation
 *
 * @param[in] origin Method asking for pointer retrieval.
 * @param[in] rsp Generic response.
 * @return Pointer to CTA response.
 *
 * @exception GException::invalid_argument
 *            Response @p rsp is not a CTA response.
 *
 * Extract CTA response from a CTA observation.
 ***************************************************************************/
const GCTAResponse* gammalib::cta_rsp(const std::string& origin,
                                      const GResponse&   rsp)
{
    // Retrieve CTA response
    const GCTAResponse* ptr = dynamic_cast<const GCTAResponse*>(&rsp);

    // If pointer is not valid then throw an exception
    if (ptr == NULL) {
        std::string cls = std::string(typeid(&rsp).name());
        std::string msg = "Invalid response type \""+cls+"\" provided. Please "
                          "specify a \"GCTAResponse\" as argument.";
        throw GException::invalid_argument(origin, msg);
    }

    // Return CTA response
    return (ptr);
}


/***********************************************************************//**
 * @brief Retrieve CTA IRF response from generic observation
 *
 * @param[in] origin Method asking for pointer retrieval.
 * @param[in] obs Generic observation.
 * @return Reference to CTA IRF response.
 *
 * @exception GException::invalid_argument
 *            Observation @p obs does not contain a CTA IRF response.
 *
 * Extract CTA IRF response from a CTA observation.
 ***************************************************************************/
const GCTAResponseIrf& gammalib::cta_rsp_irf(const std::string&  origin,
                                             const GObservation& obs)
{
    // Retrieve CTA observation
    const GCTAObservation& cta = gammalib::cta_obs(origin, obs);

    // Get pointer on CTA IRF response
    const GCTAResponseIrf* rsp = dynamic_cast<const GCTAResponseIrf*>(cta.response());

    // If pointer is not valid then throw an exception
    if (rsp == NULL) {
        std::string cls = std::string(typeid(cta.response()).name());
        std::string msg = "Invalid response type \""+cls+"\" provided in "
                          "CTA observation \""+obs.name()+"\" (ID="+obs.id()+"). "
                          "Please specify a CTA observation containing an IRF "
                          "response as argument.";
        throw GException::invalid_argument(origin, msg);
    }

    // Return CTA IRF response
    return (*rsp);
}


/***********************************************************************//**
 * @brief Retrieve CTA cube response from generic observation
 *
 * @param[in] origin Method asking for pointer retrieval.
 * @param[in] obs Generic observation.
 * @return Reference to CTA cube response.
 *
 * @exception GException::invalid_argument
 *            Observation @p obs does not contain a CTA cube response.
 *
 * Extract CTA cube response from a CTA observation.
 ***************************************************************************/
const GCTAResponseCube& gammalib::cta_rsp_cube(const std::string&  origin,
                                               const GObservation& obs)
{
    // Retrieve CTA observation
    const GCTAObservation& cta = gammalib::cta_obs(origin, obs);

    // Get pointer on CTA response cube
    const GCTAResponseCube* rsp = dynamic_cast<const GCTAResponseCube*>(cta.response());

    // If pointer is not valid then throw an exception
    if (rsp == NULL) {
        std::string cls = std::string(typeid(cta.response()).name());
        std::string msg = "Invalid response type \""+cls+"\" provided. Please "
                          "specify a CTA observation containing a response cube "
                          "as argument.";
        throw GException::invalid_argument(origin, msg);
    }

    // Return CTA response cube
    return (*rsp);
}


/***********************************************************************//**
 * @brief Retrieve CTA effective area response from generic observation
 *
 * @param[in] origin Method asking for pointer retrieval.
 * @param[in] obs Generic observation.
 * @return Reference to CTA effective area response.
 *
 * @exception GException::invalid_argument
 *            Observation @p obs does not contain a CTA effective area
 *            response.
 *
 * Extract CTA effective area response from a CTA observation.
 ***************************************************************************/
const GCTAAeff& gammalib::cta_rsp_aeff(const std::string&  origin,
                                       const GObservation& obs)
{
    // Retrieve CTA IRF response
    const GCTAResponseIrf& rsp = gammalib::cta_rsp_irf(origin, obs);

    // Get pointer on CTA effective area
    const GCTAAeff* aeff = rsp.aeff();
    if (aeff == NULL) {
        std::string msg = "Specified observation does not contain a valid "
                          "effective area response. Please specify a CTA "
                          "observation containing an effective area response.";
        throw GException::invalid_argument(origin, msg);
    }

    // Return CTA effective area response
    return (*aeff);
}


/***********************************************************************//**
 * @brief Retrieve CTA background response from generic observation
 *
 * @param[in] origin Method asking for pointer retrieval.
 * @param[in] obs Generic observation.
 * @return Reference to CTA background response.
 *
 * @exception GException::invalid_argument
 *            Observation @p obs does not contain a CTA background response.
 *
 * Extract CTA background response from a CTA observation.
 ***************************************************************************/
const GCTABackground& gammalib::cta_rsp_bkg(const std::string&  origin,
                                            const GObservation& obs)
{
    // Retrieve CTA IRF response
    const GCTAResponseIrf& rsp = gammalib::cta_rsp_irf(origin, obs);

    // Get pointer on CTA background response
    const GCTABackground* bkg = rsp.background();
    if (bkg == NULL) {
        std::string msg = "Specified observation does not contain a valid "
                          "background response. Please specify a CTA "
                          "observation containing a background response.";
        throw GException::invalid_argument(origin, msg);
    }

    // Return CTA background response
    return (*bkg);
}


/***********************************************************************//**
 * @brief Retrieve CTA event list from generic observation
 *
 * @param[in] origin Method asking for pointer retrieval.
 * @param[in] obs Generic observation.
 * @return Reference to CTA event list.
 *
 * @exception GException::invalid_argument
 *            Observation @p obs does not contain a CTA event list.
 *
 * Extract CTA event list from a CTA observation.
 ***************************************************************************/
const GCTAEventList& gammalib::cta_event_list(const std::string&  origin,
                                              const GObservation& obs)
{
    // Retrieve CTA observation
    const GCTAObservation& cta = gammalib::cta_obs(origin, obs);

    // Get pointer on CTA events list
    const GCTAEventList* events = dynamic_cast<const GCTAEventList*>(cta.events());

    // If pointer is not valid then throw an exception
    if (events == NULL) {
        std::string msg = "Specified observation does not contain a CTA event "
                          "list. Please specify a CTA observation containing "
                          "a CTA event list as argument.";
        throw GException::invalid_argument(origin, msg);
    }

    // Return CTA event list
    return (*events);
}


/***********************************************************************//**
 * @brief Retrieve CTA event cube from generic observation
 *
 * @param[in] origin Method asking for pointer retrieval.
 * @param[in] obs Generic observation.
 * @return Reference to CTA event cube.
 *
 * @exception GException::invalid_argument
 *            Observation @p obs does not contain a CTA event cube.
 *
 * Extract CTA event cube from a CTA observation.
 ***************************************************************************/
const GCTAEventCube& gammalib::cta_event_cube(const std::string&  origin,
                                              const GObservation& obs)
{
    // Retrieve CTA observation
    const GCTAObservation& cta = gammalib::cta_obs(origin, obs);

    // Get pointer on CTA event cube
    const GCTAEventCube* events = dynamic_cast<const GCTAEventCube*>(cta.events());

    // If pointer is not valid then throw an exception
    if (events == NULL) {
        std::string msg = "Specified observation does not contain a CTA event "
                          "cube. Please specify a CTA observation containing "
                          "a CTA event cube as argument.";
        throw GException::invalid_argument(origin, msg);
    }

    // Return CTA event cube
    return (*events);
}


/***********************************************************************//**
 * @brief Retrieve CTA instrument direction from generic event
 *
 * @param[in] origin Method asking for pointer retrieval.
 * @param[in] event Generic event.
 * @return Reference to CTA RoI.
 *
 * @exception GException::invalid_argument
 *            @p event does not contain a CTA instrument direction.
 *
 * Extract CTA Instrument Direction from an event.
 ***************************************************************************/
const GCTAInstDir& gammalib::cta_dir(const std::string& origin,
                                     const GEvent&      event)
{
    // Get pointer on CTA instrument direction
    const GCTAInstDir* dir = dynamic_cast<const GCTAInstDir*>(&(event.dir()));

    // If pointer is not valid then throw an exception
    if (dir == NULL) {
        std::string cls = std::string(typeid(&event).name());
        std::string msg = "Invalid event type \""+cls+"\" provided. Please "
                          "specify a \"GCTAEventAtom\" or \"GCTAEventBin\" "
                          "instance as argument.";
        throw GException::invalid_argument(origin, msg);
    }

    // Return reference
    return (*dir);
}


/***********************************************************************//**
 * @brief Retrieve spectral component from CTA background model
 *
 * @param[in] model Data model.
 * @return Pointer to spectral component.
 *
 * Retrieves the spectral component from a CTA background model. If the
 * background model has no spectral component, a NULL pointer will be
 * returned.
 ***************************************************************************/
const GModelSpectral* gammalib::cta_model_spectral(const GModelData& model)
{
    // Initialise pointer on spectral component
    GModelSpectral* spectral = NULL;

    // Cast pointer to possible CTA background models
    const GCTAModelBackground*       cta  = dynamic_cast<const GCTAModelBackground*>(&model);
    const GCTAModelAeffBackground*   aeff = dynamic_cast<const GCTAModelAeffBackground*>(&model);
    const GCTAModelIrfBackground*    irf  = dynamic_cast<const GCTAModelIrfBackground*>(&model);
    const GCTAModelCubeBackground*   cube = dynamic_cast<const GCTAModelCubeBackground*>(&model);
    const GCTAModelRadialAcceptance* rad  = dynamic_cast<const GCTAModelRadialAcceptance*>(&model);

    // Extract spectral component
    if (cta != NULL) {
        spectral = cta->spectral();
    }
    else if (aeff != NULL) {
        spectral = aeff->spectral();
    }
    else if (irf != NULL) {
        spectral = irf->spectral();
    }
    else if (cube != NULL) {
        spectral = cube->spectral();
    }
    else if (rad != NULL) {
        spectral = rad->spectral();
    }

    // Return pointer to spectral component
    return spectral;
}


/***********************************************************************//**
 * @brief Retrieve temporal component from CTA background model
 *
 * @param[in] model Data model.
 * @return Pointer to temporal component.
 *
 * Retrieves the temporal component from a CTA background model. If the
 * background model has no temporal component, a NULL pointer will be
 * returned.
 ***************************************************************************/
const GModelTemporal* gammalib::cta_model_temporal(const GModelData& model)
{
    // Initialise pointer on temporal component
    GModelTemporal* temporal = NULL;

    // Cast pointer to possible CTA background models
    const GCTAModelBackground*       cta  = dynamic_cast<const GCTAModelBackground*>(&model);
    const GCTAModelAeffBackground*   aeff = dynamic_cast<const GCTAModelAeffBackground*>(&model);
    const GCTAModelIrfBackground*    irf  = dynamic_cast<const GCTAModelIrfBackground*>(&model);
    const GCTAModelCubeBackground*   cube = dynamic_cast<const GCTAModelCubeBackground*>(&model);
    const GCTAModelRadialAcceptance* rad  = dynamic_cast<const GCTAModelRadialAcceptance*>(&model);

    // Extract temporal component
    if (cta != NULL) {
        temporal = cta->temporal();
    }
    else if (aeff != NULL) {
        temporal = aeff->temporal();
    }
    else if (irf != NULL) {
        temporal = irf->temporal();
    }
    else if (cube != NULL) {
        temporal = cube->temporal();
    }
    else if (rad != NULL) {
        temporal = rad->temporal();
    }

    // Return pointer to temporal component
    return temporal;
}
