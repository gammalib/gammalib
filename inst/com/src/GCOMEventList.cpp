/***************************************************************************
 *               GCOMEventList.cpp - COMPTEL event list class              *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2017-2021 by Juergen Knoedlseder                         *
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
#define G_DEBUG_TOFCOR                            //!< Debug TOF correction


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
 ***************************************************************************/
void GCOMEventList::read(const GFits& file)
{
    // Clear object
    clear();

    // Get HDU (pointer is always valid)
    const GFitsTable& hdu = *file.table(1);

    // Read event data
    read_events(hdu);

    // If there are events then get start and stop time from the first
    // and last event and append them as single time interval to the GTI
    if (size() > 0) {
        GTime start = m_events[0].time();
        GTime stop  = m_events[size()-1].time();
        m_gti.append(start, stop);
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write COMPTEL event list into FITS file.
 *
 * @param[in] file FITS file.
 *
 * Write the COMPTEL event list into FITS file.
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
            result.append("\n"+gammalib::parformat("TJD interval"));
            result.append(gammalib::str(gammalib::com_tjd(tstart()))+":");
            result.append(gammalib::str(gammalib::com_tics(tstart()))+" - ");
            result.append(gammalib::str(gammalib::com_tjd(tstop()))+":");
            result.append(gammalib::str(gammalib::com_tics(tstop())));
            result.append("\n"+gammalib::parformat("MJD interval"));
            result.append(gammalib::str(tstart().mjd())+" - ");
            result.append(gammalib::str(tstop().mjd())+" days");
            result.append("\n"+gammalib::parformat("UTC interval"));
            result.append(tstart().utc()+" - ");
            result.append(tstop().utc());
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

        // Get data set representation version
        int verno = table.real("DSD_REP");
        #if defined(G_DEBUG_TOFCOR)
        std::cout << "GCOMEventList::read_events: DSD_REP = ";
        std::cout << verno << std::endl;
        #endif

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

        // Disable scaling of TOF and PSD values so that the original channel
        // values are recovered. The values are divided by 128.0 below. This
        // emulates the behaviour of the evpdal13.pevpsr.f function.
        ptr_psd->scale(1.0, 0.0);
        ptr_tof->scale(1.0, 0.0);

        // Initialise boundaries
        GEnergy emin;
        GEnergy emax;

        // Debug TOFCOR: initialise flag for first notice
        #if defined(G_DEBUG_TOFCOR)
        bool noticed = false;
        #endif

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

            // Set PSD and TOF values
            double psd = double(ptr_psd->integer(i))/128.0;
            double tof = double(ptr_tof->integer(i))/128.0;

            // Correct TOF if data representation flag is less than 3 and
            // rejection flag is >=4 (see dalaaa19.pevpsr.f)
            if (verno < 3) {
                if (ptr_reflag->integer(i) >= 4) {
                    #if defined(G_DEBUG_TOFCOR)
                    double tof_raw = tof;
                    #endif
                    tof = tofcor(ptr_e1->real(i), ptr_e2->real(i), tof);
                    #if defined(G_DEBUG_TOFCOR)
                    if (!noticed) {
                        std::cout << "GCOMEventList::read_events: ";
                        std::cout << "apply TOF correction (";
                        std::cout << tof_raw;
                        std::cout << " => ";
                        std::cout << tof;
                        std::cout << ")" << std::endl;
                        noticed = true;
                    }
                    #endif
                }
            }

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
            event.psd(psd);
            event.tof(tof);
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


/***********************************************************************//**
 * @brief Compute TOF correction
 *
 * @param[in] d1e D1 energy deposit (keV).
 * @param[in] d2e D2 energy deposit (keV).
 * @param[in] tof TOF value from EVP dataset (version < 3).
 *
 * Computes the TOF correction for EVP data sets with representation version
 * < 3.
 ***************************************************************************/
double GCOMEventList::tofcor(const double& d1e, const double& d2e, double tof) const
{
    // Begin and end of D1E-range for which fit is valid
    const double d1sta =    50.0;
    const double d1end = 20000.0;

    // Begin and end of D2E-range for which fit is valid
    const double d2sta =   600.0;
    const double d2end = 30000.0;

    // Degree of polynomial fit to ToF(D1E) for which fit is valid
    const int ndegd1 = 6;

    // D1-energy 'breaking' the D1E-fit range in 2 intervals
    const double d1ebrk = 2250.0;

    // Degree of polynomial fit to ToF(D2E) for which fit is valid
    const int ndegd2 = 4;

    // D2-energy 'breaking' the D2E-fit range in 3 intervals
    const double d2ebrk1 = 1400.0;
    const double d2ebrk2 = 5500.0;

    // Reference ToF: the ToF corrections are calculated with respect to this TOF
    const double tofref = 118.30;

    // ToF the corrected peak will be centered at
    const double tofcen = 120.0;

    // Polynomial coefficients of ToF(D1E) for interval 1: D1E < D1EBRK
    const double d11cof[] = {111.74858,
                              28.280247,
                             -45.024305,
                              35.183210,
                             -14.639463,
                               3.1342536,
                              -0.2711735};

    // Polynomial coefficients of ToF(D1E) for interval 2: D1E > D1EBRK
    const double d12cof[] = {116.25374,
                               0.500407092,
                               0.38182720,
                              -0.080145513,
                               0.0065569790,
                              -0.00024650067,
                               3.5077240e-6};

    // Polynomial coefficients of ToF(D2E) for interval 1: D2E < D2EBRK1
    const double d21cof[] = { 181.77024,
                             -252.41070,
                              371.09898,
                             -232.83985,
                               52.918785};

    // Polynomial coefficients of ToF(D2E) for interval 2: D2EBRK1 < D2E < D2EBRK2
    const double d22cof[] = { 120.91608,
                               -0.15048490,
                               -0.45526025,
                                0.11710009,
                               -0.0082172427};

    // Polynomial coefficients of ToF(D2E) for interval 3: D2EBRK2 < D2E
    const double d23cof[] = { 119.24278,
                               -0.43134699,
                                0.060183080,
                               -0.0026847790,
                                3.7720986e-5};

    // Make sure that D1 and D2 energies are within fit range
    double d1ener = d1e;
    double d2ener = d2e;
    if (d1ener < d1sta) {
        d1ener = d1sta;
    }
    if (d1ener > d1end) {
        d1ener = d1end;
    }
    if (d2ener < d2sta) {
        d2ener = d2sta;
    }
    if (d2ener > d2end) {
        d2ener = d2end;
    }

    // Compute D1 correction
    double d1ecor = -tofref;
    double d1emev = d1ener * 1.0e-3;
    double d1edum = 1.0;
    if (d1ener <= d1ebrk) {
        for (int i = 0; i <= ndegd1; ++i) {
            d1ecor += d11cof[i] * d1edum;
            d1edum *= d1emev;
        }
    }
    else {
        for (int i = 0; i <= ndegd1; ++i) {
            d1ecor += d12cof[i] * d1edum;
            d1edum *= d1emev;
        }
    }

    // Compute D2 correction
    double d2ecor = -tofref;
    double d2emev = d2ener * 1.0e-3;
    double d2edum = 1.0;
    if (d2ener <= d2ebrk1) {
        for (int i = 0; i <= ndegd2; ++i) {
            d2ecor += d21cof[i] * d2edum;
            d2edum *= d2emev;
        }
    }
    else if (d2ener > d2ebrk2) {
        for (int i = 0; i <= ndegd2; ++i) {
            d2ecor += d23cof[i] * d2edum;
            d2edum *= d2emev;
        }
    }
    else {
        for (int i = 0; i <= ndegd2; ++i) {
            d2ecor += d22cof[i] * d2edum;
            d2edum *= d2emev;
        }
    }

    // Compute total correction
    tof += tofcen - (tofref + d1ecor + d2ecor);

    // Return corrected TOF
    return tof;
}
