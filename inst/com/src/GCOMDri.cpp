/***************************************************************************
 *                  GCOMDri.cpp - COMPTEL Data Space class                 *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2017 by Juergen Knoedlseder                              *
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
 * @file GCOMDri.cpp
 * @brief COMPTEL Data Space class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GWcs.hpp"
#include "GFits.hpp"
#include "GFitsImage.hpp"
#include "GCOMDri.hpp"
#include "GCOMSupport.hpp"

/* __ Method name definitions ____________________________________________ */

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */



/*==========================================================================
 =                                                                         =
 =                         Constructors/destructors                        =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Void constructor
 ***************************************************************************/
GCOMDri::GCOMDri(void)
{
    // Initialise class members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief File name constructor
 *
 * @param[in] filename COMPTEL Data Space FITS file name.
 ***************************************************************************/
GCOMDri::GCOMDri(const GFilename& filename)
{
    // Initialise class members
    init_members();

    // Load DRI FITS file
    load(filename);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] dri COMPTEL Data Space.
 ***************************************************************************/
GCOMDri::GCOMDri(const GCOMDri& dri)
{
    // Initialise class members
    init_members();

    // Copy members
    copy_members(dri);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GCOMDri::~GCOMDri(void)
{
    // Free members
    free_members();

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                               Operators                                 =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Assignment operator
 *
 * @param[in] dri COMPTEL Data Space.
 * @return COMPTEL Data Space.
 ***************************************************************************/
GCOMDri& GCOMDri::operator=(const GCOMDri& dri)
{
    // Execute only if object is not identical
    if (this != &dri) {

        // Free members
        free_members();

        // Initialise private members
        init_members();

        // Copy members
        copy_members(dri);

    } // endif: object was not identical

    // Return this object
    return *this;
}


/*==========================================================================
 =                                                                         =
 =                             Public methods                              =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Clear COMPTEL Data Space
 ***************************************************************************/
void GCOMDri::clear(void)
{
    // Free members
    free_members();

    // Initialise private members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone COMPTEL Data Space
 *
 * @return Pointer to deep copy of COMPTEL Data Space.
 ***************************************************************************/
GCOMDri* GCOMDri::clone(void) const
{
    return new GCOMDri(*this);
}


/***********************************************************************//**
 * @brief Load COMPTEL Data Space from DRI FITS file
 *
 * @param[in] filename DRI FITS file name.
 ***************************************************************************/
void GCOMDri::load(const GFilename& filename)
{
    // Open FITS file
    GFits fits(filename);

    // Get HDU (pointer is always valid)
    const GFitsImage& hdu = *fits.image(0);

    // Read DRI file
    read(hdu);

    // Close FITS file
    fits.close();
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Save COMPTEL Data Space into DRI FITS file
 *
 * @param[in] filename DRI FITS file name.
 ***************************************************************************/
void GCOMDri::save(const GFilename& filename, const bool& clobber) const
{
    // Create FITS file
    GFits fits;

    // Write data space into FITS file
    write(fits, filename.extname(gammalib::extname_dri));

    // Save FITS file
    fits.saveto(filename, clobber);

    // Close FITS file
    fits.close();
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Read COMPTEL Data Space from DRI FITS image
 *
 * @param[in] image DRI FITS image.
 ***************************************************************************/
void GCOMDri::read(const GFitsImage& image)
{
    // Clear
    clear();

    // Read sky map
    m_dri.read(image);

    // Correct WCS projection (HEASARC data format kluge)
    com_wcs_mer2car(m_dri);

    // Read attributes
    read_attributes(&image);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write COMPTEL Data Space into FITS image
 *
 * @param[in] fits FITS file.
 * @param[in] extname Extension name.
 ***************************************************************************/
void GCOMDri::write(GFits& fits, const std::string& extname) const
{
    // Write sky map into FITS file
    GFitsHDU *image = m_dri.write(fits, extname);

    // Write DRI attributes
    write_attributes(image);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Print COMPTEL Data Space
 *
 * @param[in] chatter Chattiness.
 * @return String containing COMPTEL Data Space information.
 *
 * @todo Implement method.
 ***************************************************************************/
std::string GCOMDri::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Compute scatter angle dimensions
        double      chimin = 0.0;
        double      chimax = 0.0;
        double      chibin = 0.0;
        double      psimin = 0.0;
        double      psimax = 0.0;
        double      psibin = 0.0;
        const GWcs* wcs    = dynamic_cast<const GWcs*>(m_dri.projection());
        if (wcs != NULL) {
            chibin = wcs->cdelt(0);
            chimin = wcs->crval(0) - (wcs->crpix(0)-0.5) * chibin;
            chimax = chimin + m_dri.nx() * chibin;
            psibin = wcs->cdelt(1);
            psimin = wcs->crval(1) - (wcs->crpix(1)-0.5) * psibin;
            psimax = psimin + m_dri.ny() * psibin;
        }

        // Compute Phibar maximum
        double phimax = m_phimin + m_dri.nmaps() * m_phibin;

        // Append header
        result.append("=== GCOMDri ===");

        // Append Phibar information
        result.append("\n"+gammalib::parformat("Chi range"));
        result.append(gammalib::str(chimin)+" - ");
        result.append(gammalib::str(chimax)+" deg");
        result.append("\n"+gammalib::parformat("Chi bin size"));
        result.append(gammalib::str(chibin)+" deg");
        result.append("\n"+gammalib::parformat("Psi range"));
        result.append(gammalib::str(psimin)+" - ");
        result.append(gammalib::str(psimax)+" deg");
        result.append("\n"+gammalib::parformat("Psi bin size"));
        result.append(gammalib::str(psibin)+" deg");
        result.append("\n"+gammalib::parformat("Phibar range"));
        result.append(gammalib::str(m_phimin)+" - ");
        result.append(gammalib::str(phimax)+" deg");
        result.append("\n"+gammalib::parformat("Phibar bin size"));
        result.append(gammalib::str(m_phibin)+" deg");

        // Append energy boundaries
        result.append("\n"+m_ebounds.print(gammalib::reduce(chatter)));

        // Append GTI
        result.append("\n"+m_gti.print(gammalib::reduce(chatter)));

        // EXPLICIT: Append sky map
        if (chatter >= EXPLICIT) {
            result.append("\n"+m_dri.print(gammalib::reduce(chatter)));
        }

    } // endif: chatter was not silent

    // Return result
    return result;
}


/*==========================================================================
 =                                                                         =
 =                             Private methods                             =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void GCOMDri::init_members(void)
{
    // Initialise members
    m_dri.clear();
    m_ebounds.clear();
    m_gti.clear();
    m_phimin = 0.0;
    m_phibin = 0.0;
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] dri COMPTEL Data Space.
 ***************************************************************************/
void GCOMDri::copy_members(const GCOMDri& dri)
{
    // Copy members
    m_dri     = dri.m_dri;
    m_ebounds = dri.m_ebounds;
    m_gti     = dri.m_gti;
    m_phimin  = dri.m_phimin;
    m_phibin  = dri.m_phibin;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GCOMDri::free_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Read DRI attributes from FITS HDU
 *
 * @param[in] hdu FITS HDU pointer.
 ***************************************************************************/
void GCOMDri::read_attributes(const GFitsHDU* hdu)
{
    // Get phibar attributes
    m_phibin = hdu->real("CDELT3");
    m_phimin = hdu->real("CRVAL3") - (hdu->real("CRPIX3")-0.5) * m_phibin;
    
    // Get time attributes
    GTime tstart = com_time(hdu->integer("VISDAY"), hdu->integer("VISTIM"));
    GTime tstop  = com_time(hdu->integer("VIEDAY"), hdu->integer("VIETIM"));

    // Set Good Time Intervals
    m_gti = GGti(tstart, tstop);

    // Optionally read energy attributes
    if (hdu->has_card("E_MIN") && hdu->has_card("E_MAX")) {

        // Get energy attributes
        GEnergy emin(hdu->real("E_MIN"), "MeV");
        GEnergy emax(hdu->real("E_MAX"), "MeV");

        // Set energy boundaries
        m_ebounds = GEbounds(emin, emax);

    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write DRI attributes into FITS HDU
 *
 * @param[in] hdu FITS HDU pointer.
 ***************************************************************************/
void GCOMDri::write_attributes(GFitsHDU* hdu) const
{
    // Set Phibar keywords
    double crval3 = m_phimin + 0.5 * m_phibin;

    // Write Phibar keywords
    hdu->card("CTYPE3", "Phibar", "Compton scatter angle");
    hdu->card("CRPIX3", 1.0, "Pixel coordinate of reference point (starting from 1)");
    hdu->card("CRVAL3", crval3, "[deg] Coordinate value at reference point");
    hdu->card("CDELT3", m_phibin, "[deg] Coordinate increment at reference point");

    // Write OGIP time keywords
    m_gti.reference().write(*hdu);
    hdu->card("TSTART", m_gti.tstart().secs(), "[s] Start time");
    hdu->card("TSTOP",  m_gti.tstop().secs(), "[s] Stop time");
    hdu->card("DATE-OBS", m_gti.tstart().utc(), "Start of observation in UTC");
    hdu->card("DATE-END", m_gti.tstop().utc(), "Stop of observation in UTC");

    // Set time keywords
    int visday = com_tjd(m_gti.tstart());
    int vistim = com_tics(m_gti.tstart());
    int vieday = com_tjd(m_gti.tstop());
    int vietim = com_tics(m_gti.tstop());

    // Write COMPTEL time keywords
    hdu->card("VISDAY", visday, "[TJD] Data validity interval start day");
    hdu->card("VISTIM", vistim, "[tics] Data validity interval start time");
    hdu->card("VIEDAY", vieday, "[TJD] Data validity interval end day");
    hdu->card("VIETIM", vietim, "[tics] Data validity interval start time");

    // If there are energy boundaries then write them
    if (m_ebounds.size() > 0) {
        hdu->card("E_MIN", m_ebounds.emin().MeV(), "[MeV] Lower bound of energy range");
        hdu->card("E_MAX", m_ebounds.emax().MeV(), "[MeV] Upper bound of energy range");
    }

    // Return
    return;
}
