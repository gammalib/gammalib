/***************************************************************************
 *               GCOMEventList.cpp - COMPTEL event list class              *
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
 * @file GCOMEventList.hpp
 * @brief COMPTEL event list class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <typeinfo>
#include "GMath.hpp"
#include "GFits.hpp"
#include "GException.hpp"
#include "GCOMSupport.hpp"
#include "GCOMEventList.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_OPERATOR                          "GCOMEventList::operator[](int&)"
#define G_ROI                                     "GCOMEventList::roi(GRoi&)"

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
 *
 * Creates an empty COMPTEL event list.
 ***************************************************************************/
GCOMEventList::GCOMEventList(void) : GEventList()
{
    // Initialise class members for clean destruction
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief File name constructor
 *
 * @param[in] filename COMPTEL event list filename.
 *
 * Construct COMPTEL event list object by loading the events from a
 * FITS file.
 ***************************************************************************/
GCOMEventList::GCOMEventList(const GFilename& filename) : GEventList()
{
    // Initialise members
    init_members();

    // Load event list
    load(filename);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] list COMPTEL event list.
 ***************************************************************************/
GCOMEventList::GCOMEventList(const GCOMEventList& list) : GEventList(list)
{
    // Initialise members
    init_members();

    // Copy members
    copy_members(list);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GCOMEventList::~GCOMEventList(void)
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
 * @param[in] list COMPTEL event list.
 * @return COMPTEL event list.
 ***************************************************************************/
GCOMEventList& GCOMEventList::operator=(const GCOMEventList& list)
{
    // Execute only if object is not identical
    if (this != &list) {

        // Copy base class members
        this->GEventList::operator=(list);

        // Free members
        free_members();

        // Initialise members
        init_members();

        // Copy members
        copy_members(list);

    } // endif: object was not identical

    // Return this object
    return *this;
}


/***********************************************************************//**
 * @brief COMPTEL event atom access operator
 *
 * @param[in] index Event index [0,...,size()-1].
 * @return Pointer to COMPTEL event atom.
 *
 * @exception GException::out_of_range
 *            Event index outside valid range.
 *
 * Returns pointer to a COMPTEL event atom.
 ***************************************************************************/
GCOMEventAtom* GCOMEventList::operator[](const int& index)
{
    // Optionally check if the index is valid
    #if defined(G_RANGE_CHECK)
    if (index < 0 || index >= size()) {
        throw GException::out_of_range(G_OPERATOR, "Event index", index, size());
    }
    #endif

    // Return pointer
    return (&(m_events[index]));
}


/***********************************************************************//**
 * @brief COMPTEL event atom access operator
 *
 * @param[in] index Event index [0,...,size()-1].
 * @return Pointer to COMPTEL event atom.
 *
 * @exception GException::out_of_range
 *            Event index outside valid range.
 *
 * Returns pointer to a COMPTEL event atom.
 ***************************************************************************/
const GCOMEventAtom* GCOMEventList::operator[](const int& index) const
{
    // Optionally check if the index is valid
    #if defined(G_RANGE_CHECK)
    if (index < 0 || index >= size()) {
        throw GException::out_of_range(G_OPERATOR, "Event index", index, size());
    }
    #endif

    // Return pointer
    return (&(m_events[index]));
}


/*==========================================================================
 =                                                                         =
 =                             Public methods                              =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Clear COMPTEL event list
 *
 * Clears COMPTEL event list by resetting all class members to an
 * initial state. Any information that was present before will be lost.
 ***************************************************************************/
void GCOMEventList::clear(void)
{
    // Free class members (base and derived classes, derived class first)
    free_members();
    this->GEventList::free_members();
    this->GEvents::free_members();

    // Initialise members
    this->GEvents::init_members();
    this->GEventList::init_members();
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone event list
 *
 * @return Pointer to deep copy of COMPTEL event list.
 ***************************************************************************/
GCOMEventList* GCOMEventList::clone(void) const
{
    return new GCOMEventList(*this);
}


/***********************************************************************//**
 * @brief Load COMPTEL events from FITS file
 *
 * @param[in] filename COMPTEL event list FITS file name.
 *
 * Loads COMPTEL events from a FITS file into the event list.
 ***************************************************************************/
void GCOMEventList::load(const GFilename& filename)
{
    // Open FITS file
    GFits fits(filename);

    // Read event list from FITS file
    read(fits);

    // Close FITS file
    fits.close();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Save COMPTEL events
 *
 * @param[in] filename COMPTEL event list FITS file name.
 * @param[in] clobber Overwrite existing FITS file?
 *
 * Save COMPTEL events from a FITS file into the event list.
 ***************************************************************************/
void GCOMEventList::save(const GFilename& filename,
                         const bool&      clobber) const
{
    // Open FITS file
    GFits fits;

    // Write events
    write(fits);

    // Save events
    fits.saveto(filename, clobber);

    // Close FITS file
    fits.close();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Read COMPTEL events from FITS file.
 *
 * @param[in] file FITS file.
 *
 * Read the COMPTEL event list from a FITS file object.
 *
 * @todo Implement method.
 ***************************************************************************/
void GCOMEventList::read(const GFits& file)
{
    // Clear object
    clear();

    // Get HDU (pointer is always valid)
    const GFitsTable& hdu = *file.table(1);

    // Read event data
    read_events(hdu);

    // Read start and stop time
    GTime start(com_time(hdu.integer("VISDAY"), hdu.integer("VISTIM")));
    GTime stop(com_time(hdu.integer("VIEDAY"), hdu.integer("VIETIM")));
    
    // Append start and stop time as single time interval to GTI
    m_gti.append(start, stop);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write COMPTEL event list into FITS file.
 *
 * @param[in] file FITS file.
 *
 * Write the LCOMPTELAT event list into FITS file.
 *
 * @todo Implement method.
 ***************************************************************************/
void GCOMEventList::write(GFits& file) const
{
    // TODO: You need to implement an interface to the FITS file that so
    // that you can save your events.

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set region of interest
 *
 * @param[in] roi Region of interest.
 *
 * @exception GException::invalid_argument
 *            Specified RoI is not a COMPTEL RoI.
 *
 * Sets the region of interest for the observation.
 ***************************************************************************/
void GCOMEventList::roi(const GRoi& roi)
{
    // Get pointer on COMPTEL region of interest
    const GCOMRoi* comroi = dynamic_cast<const GCOMRoi*>(&roi);
    if (comroi == NULL) {
        std::string cls = std::string(typeid(&roi).name());
        std::string msg = "Region of interest of type \""+cls+"\" is "
                          "not a COMPTEL RoI. Please specify a "
                          "COMPTEL RoI as argument.";
        throw GException::invalid_argument(G_ROI, msg);
    }

    // Set RoI
    m_roi = *comroi;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Append event to event list
 *
 * @param[in] event Event.
 *
 * Appends an event to the end of the event list.
 ***************************************************************************/
void GCOMEventList::append(const GCOMEventAtom& event)
{
    // Append event
    m_events.push_back(event);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Remove events from event list
 *
 * @param[in] index Index from which on events should be removed.
 * @param[in] number Number of event to remove (default: 1).
 *
 * Removes events from the event list. This method does nothing if @p index
 * points beyond the event list. The method does also gently accept
 * @p number arguments where @p index + @p number reach beyond the event
 * list. In that case, all events from event @p index on will be removed.
 ***************************************************************************/
void GCOMEventList::remove(const int& index, const int& number)
{
    // Continue only if index is valid
    if (index < size()) {

        // Determine number of elements to remove
        int n_remove = (index + number > size()) ? size() - index : number;

        // Remove events
        m_events.erase(m_events.begin() + index,
                       m_events.begin() + index + n_remove);

    } // endif: index was valid

    // Return
    return;
}


/***********************************************************************//**
 * @brief Print COMPTEL event list information
 *
 * @param[in] chatter Chattiness.
 * @return String containing COMPTEL event list information.
 ***************************************************************************/
std::string GCOMEventList::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GCOMEventList ===");

        // Append information
        result.append("\n"+gammalib::parformat("Number of events") +
                      gammalib::str(number()));

        // Append time information
        if (gti().size() > 0) {
            result.append("\n"+gammalib::parformat("MJD interval"));
            result.append(gammalib::str(tstart().mjd())+" - ");
            result.append(gammalib::str(tstop().mjd())+" days");
            result.append("\n"+gammalib::parformat("UTC interval"));
            result.append(tstart().utc()+" - ");
            result.append(tstop().utc());
            result.append("\n"+gammalib::parformat("TJD interval"));
            result.append(gammalib::str(com_tjd(tstart()))+":");
            result.append(gammalib::str(com_tics(tstart()))+" - ");
            result.append(gammalib::str(com_tjd(tstop()))+":");
            result.append(gammalib::str(com_tics(tstop())));
        }
        else {
            result.append("\n"+gammalib::parformat("MJD interval"));
            result.append("not defined");
        }

        // Append energy information
        result.append("\n"+gammalib::parformat("Energy interval"));
        if (ebounds().size() > 0) {
            result.append(gammalib::str(ebounds().emin().MeV()));
            result.append(" - ");
            result.append(gammalib::str(ebounds().emax().MeV())+" MeV");
        }
        else {
            result.append("not defined");
        }

    } // endif: chatter was not silent

    // Return result
    return result;
}


/*==========================================================================
 =                                                                         =
 =                            Private methods                              =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void GCOMEventList::init_members(void)
{
    // Initialise members
    m_roi.clear();
    m_events.clear();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] list COMPTEL event list.
 ***************************************************************************/
void GCOMEventList::copy_members(const GCOMEventList& list)
{
    // Copy members
    m_roi    = list.m_roi;
    m_events = list.m_events;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GCOMEventList::free_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Read COMPTEL events from FITS table
 *
 * @param[in] table FITS table.
 *
 * Reads COMPTEL events from a FITS table. The method expects the data in the
 * format that is provided by the HEASARC archive.
 *
 * The method also sets the energy boundaries and Region of Interest of the
 * data set as this information is not provided in the FITS table header.
 * For that purpose it scans all data and determines the range that is
 * covered by the data.
 *
 * The method assumes that no events exist so far.
 ***************************************************************************/
void GCOMEventList::read_events(const GFitsTable& table)
{
    // Extract number of events in FITS file
    int num = table.nrows();

    // If there are events then load them
    if (num > 0) {

        // Use pointing direction as roi centre
        GSkyDir centre;
        centre.lb_deg(table.real("GLON_SCZ"), table.real("GLAT_SCZ"));

        // Initialise boundaries
        double radius_max = 0.0;
        double phibar_min = 0.0;
        double phibar_max = 0.0;

        // Reserve data
        m_events.reserve(num);

        // Get column pointers
        const GFitsTableCol* ptr_tjd        = table["TJD"];           // days
        const GFitsTableCol* ptr_tics       = table["TICS"];          // ticks
        const GFitsTableCol* ptr_glon_scat  = table["GLON_SCAT"];     // rad
        const GFitsTableCol* ptr_glat_scat  = table["GLAT_SCAT"];     // rad
        const GFitsTableCol* ptr_phi_scat   = table["AZIMUTH_SCAT"];  // rad
        const GFitsTableCol* ptr_theta_scat = table["ZENITH_SCAT"];   // rad
        const GFitsTableCol* ptr_phibar     = table["PHIBAR"];        // rad
        const GFitsTableCol* ptr_eha        = table["EARTH_HORIZON"]; // rad
        const GFitsTableCol* ptr_e1         = table["E_D1"];          // keV
        const GFitsTableCol* ptr_e2         = table["E_D2"];          // keV
        const GFitsTableCol* ptr_psd        = table["PSD"];           // channel
        const GFitsTableCol* ptr_tof        = table["TOF"];           // channel
        const GFitsTableCol* ptr_modcom     = table["MODCOM"];        // id
        const GFitsTableCol* ptr_reflag     = table["RC_REFLAG"];     //
        const GFitsTableCol* ptr_veto       = table["RC_VETO"];       //

        // Disable scaling of TOF and PSD values so that the original
        // channel values are recovered
        //ptr_psd->scale(1.0, 0.0);
        //ptr_tof->scale(1.0, 0.0);

        // Initialise boundaries
        GEnergy emin;
        GEnergy emax;

        // Copy data from columns into GCOMEventAtom objects
        for (int i = 0; i < num; ++i) {

            // Allocate event
            GCOMEventAtom event;

            // Set instrument direction. Note that contrary to what is
            // specified in the FITS file, angles are provided in radians.
            // Also, GLON and GLAT are inverted.
            GCOMInstDir inst_dir;
            GSkyDir     sky_dir;
            sky_dir.lb(ptr_glat_scat->real(i), ptr_glon_scat->real(i));
            inst_dir.dir(sky_dir);
            inst_dir.phibar(ptr_phibar->real(i) * gammalib::rad2deg);

            // Set total energy
            GEnergy etot;
            etot.keV(ptr_e1->real(i) + ptr_e2->real(i));

            // Set phibar value (deg)
            double phibar = ptr_phibar->real(i) * gammalib::rad2deg;

            // Set event information
            event.time(ptr_tjd->integer(i), ptr_tics->integer(i));
            event.energy(etot);
            event.dir(inst_dir);
            event.phi(ptr_phi_scat->real(i) * gammalib::rad2deg);     // rad -> deg
            event.theta(ptr_theta_scat->real(i) * gammalib::rad2deg); // rad -> deg
            event.phibar(phibar);                                     // rad -> deg
            event.eha(ptr_eha->real(i) * gammalib::rad2deg);          // rad -> deg
            event.e1(ptr_e1->real(i) * 1.0e-3);                       // keV -> MeV
            event.e2(ptr_e2->real(i) * 1.0e-3);                       // keV -> MeV
            event.psd(ptr_psd->integer(i));
            event.tof(ptr_tof->integer(i));
            event.modcom(ptr_modcom->integer(i));
            event.reflag(ptr_reflag->integer(i));
            event.veto(ptr_veto->integer(i));

            // Append event
            m_events.push_back(event);

            // Compute distance from pointing direction
            double radius = centre.dist_deg(sky_dir);

            // Update boundaries
            if (i == 0) {
                emin       = etot;
                emax       = etot;
                phibar_min = phibar;
                phibar_max = phibar;
                radius_max = radius;
            }
            else {
                if (etot < emin) {
                    emin = etot;
                }
                if (etot > emax) {
                    emax = etot;
                }
                if (phibar < phibar_min) {
                    phibar_min = phibar;
                }
                if (phibar > phibar_max) {
                    phibar_max = phibar;
                }
                if (radius > radius_max) {
                    radius_max = radius;
                }
            }

        } // endfor: looped over all events

        // Set region of interest
        double phibar = 0.5 * (phibar_min + phibar_max);
        m_roi         = GCOMRoi(GCOMInstDir(centre, phibar), radius_max,
                                phibar_min, phibar_max);

        // Set energy boundaries
        m_ebounds = GEbounds(emin, emax);

    } // endif: there were events

    // Return
    return;
}
