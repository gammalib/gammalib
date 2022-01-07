/***************************************************************************
 *             GCOMSelection.cpp - COMPTEL selection set class             *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2017-2022 by Juergen Knoedlseder                         *
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
 * @file GCOMSelection.cpp
 * @brief COMPTEL selection set class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GFitsHDU.hpp"
#include "GCOMTools.hpp"
#include "GCOMSupport.hpp"
#include "GCOMStatus.hpp"
#include "GCOMSelection.hpp"
#include "GCOMEventAtom.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_FPMTFLAG                            "GCOMSelection::fpmtflag(int&)"
#define G_USE_D1_GET                            "GCOMSelection::use_d1(int&)"
#define G_USE_D1_SET                     "GCOMSelection::use_d1(int&, bool&)"
#define G_USE_D2_GET                            "GCOMSelection::use_d2(int&)"
#define G_USE_D2_SET                     "GCOMSelection::use_d2(int&, bool&)"

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
GCOMSelection::GCOMSelection(void)
{
    // Initialise class members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] select COMPTEL selection set.
 ***************************************************************************/
GCOMSelection::GCOMSelection(const GCOMSelection& select)
{
    // Initialise class members
    init_members();

    // Copy members
    copy_members(select);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GCOMSelection::~GCOMSelection(void)
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
 * @param[in] select COMPTEL selection set.
 * @return COMPTEL selection set.
 ***************************************************************************/
GCOMSelection& GCOMSelection::operator=(const GCOMSelection& select)
{
    // Execute only if object is not identical
    if (this != &select) {

        // Free members
        free_members();

        // Initialise private members
        init_members();

        // Copy members
        copy_members(select);

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
 * @brief Clear COMPTEL selection set
 ***************************************************************************/
void GCOMSelection::clear(void)
{
    // Free members
    free_members();

    // Initialise private members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone COMPTEL selection set
 *
 * @return Pointer to deep copy of COMPTEL selection set.
 ***************************************************************************/
GCOMSelection* GCOMSelection::clone(void) const
{
    return new GCOMSelection(*this);
}


/***********************************************************************//**
 * @brief Initialise selection statistics
 ***************************************************************************/
void GCOMSelection::init_statistics(void) const
{
    // Initialise statistics members
    m_num_events_checked  = 0;
    m_num_events_used     = 0;
    m_num_events_rejected = 0;
    m_num_e1_min          = 0;
    m_num_e1_max          = 0;
    m_num_e2_min          = 0;
    m_num_e2_max          = 0;
    m_num_tof_min         = 0;
    m_num_tof_max         = 0;
    m_num_psd_min         = 0;
    m_num_psd_max         = 0;
    m_num_reflag_min      = 0;
    m_num_reflag_max      = 0;
    m_num_vetoflag_min    = 0;
    m_num_vetoflag_max    = 0;
    m_num_no_scatter      = 0;
    m_num_invalid_modcom  = 0;
    m_num_d1module_off    = 0;
    m_num_d2module_off    = 0;
    m_num_fpmt            = 0;
    for (int i = 0; i < 7; ++i) {
        m_num_d1[i] = 0;
    }
    for (int i = 0; i < 14; ++i) {
        m_num_d2[i] = 0;
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Check if event should be used
 *
 * @param[in] event Event.
 * @return True if event should be used, false otherwise.
 *
 * Checks if an event should be used for DRE binning.
 ***************************************************************************/
bool GCOMSelection::use_event(const GCOMEventAtom& event) const
{
    // Initialise D1 & D2 module status
    static const GCOMStatus status;

    // Initialise usage flag
    bool use = true;

    // Increment number of checked events
    m_num_events_checked++;

    // Check for bad minitelescopes
    if (event.modcom() < 1 || event.modcom() > 98) {
        m_num_invalid_modcom++;
        use = false;
    }

    // Check whether the event has a scatter angle determined. This is
    // signaled by a scatter angle -10e20 in radians.
    else if (event.theta() < -1.0e3) {
        m_num_no_scatter++;
        use = false;
    }

    // Apply event selection
    else if (event.e1() < m_e1_min) {
        m_num_e1_min++;
        use = false;
    }
    else if (event.e1() > m_e1_max) {
        m_num_e1_max++;
        use = false;
    }
    else if (event.e2() < m_e2_min) {
        m_num_e2_min++;
        use = false;
    }
    else if (event.e2() > m_e2_max) {
        m_num_e2_max++;
        use = false;
    }
    else if (event.tof() < m_tof_min) {
        m_num_tof_min++;
        use = false;
    }
    else if (event.tof() > m_tof_max) {
        m_num_tof_max++;
        use = false;
    }
    else if (event.psd() < m_psd_min) {
        m_num_psd_min++;
        use = false;
    }
    else if (event.psd() > m_psd_max) {
        m_num_psd_max++;
        use = false;
    }
    else if (event.reflag() < m_reflag_min) {
        m_num_reflag_min++;
        use = false;
    }
    else if (event.reflag() > m_reflag_max) {
        m_num_reflag_max++;
        use = false;
    }
    else if (event.veto() < m_vetoflag_min) {
        m_num_vetoflag_min++;
        use = false;
    }
    else if (event.veto() > m_vetoflag_max) {
        m_num_vetoflag_max++;
        use = false;
    }

    // If event should be used then handle now the module status and
    // decide whether it should be rejected
    if (use) {

        // Extract module IDs from MODCOM
        int id2 = (event.modcom()-1)/7 + 1;      // [1-14]
        int id1 =  event.modcom() - (id2-1) * 7; // [1-7]

        // Get event TJD
        int tjd = gammalib::com_tjd(event.time());

        // Set D1 and D2 status
        int d1status = status.d1status(tjd, id1);
        int d2status = status.d2status(tjd, id2);

        // Exclude event if the coresponding modules should not be used
        if (!m_use_d1[id1-1]) {
            m_num_d1module_off++;
            use = false;
        }
        else if (!m_use_d2[id2-1]) {
            m_num_d2module_off++;
            use = false;
        }

        // Exclude event if the corresponding modules are signalled
        // as off
        else if (d1status != 1) {
            m_num_d1module_off++;
            use = false;
        }
        else if (d2status < 1) {
            m_num_d2module_off++;
            use = false;
        }

        // Handle D2 modules with failed PMTs based on value of
        // failure flag
        else if (d2status > 1) {

            // If failure flag = 0 then reject event and signal that the
            // module was off
            if (m_fpmtflag == 0) {
                m_num_d2module_off++;
                use = false;
            }

            // ... otherwise if failure flag = 2 and the module has a
            // valid exclusion circle then reject event if it is
            // within the exclusion circle; if the module has no
            // exclusion region then reject the entire module and
            // signal that it was off.
            else if (m_fpmtflag == 2) {
                double r = gammalib::com_exd2r(id2);
                if (r >= 0.1) {
                    double dx  = gammalib::com_exd2x(id2) - 0.1 * event.x_d2();
                    double dy  = gammalib::com_exd2y(id2) - 0.1 * event.y_d2();
                    double dsq = dx*dx + dy*dy;
                    double rsq = r*r;
                    if (dsq <= rsq) {
                        m_num_fpmt++;
                        use = false;
                    }
                }
                else {
                    m_num_d2module_off++;
                    use = false;
                }
            }

        } // endif: handled D2 modules with failed PMTs

    } // endif: handled module status

    // Update acceptance and rejection statistics
    if (use) {

        // Extract module IDs from MODCOM
        int id2 = (event.modcom()-1)/7 + 1;      // [1-14]
        int id1 =  event.modcom() - (id2-1) * 7; // [1-7]

        // Update statistics
        m_num_events_used++;
        m_num_d1[id1-1]++;
        m_num_d2[id2-1]++;

    }
    else {
        m_num_events_rejected++;
    }

    // Return usage flag
    return use;
}


/***********************************************************************//**
 * @brief Set failed PMT flag for D2 modules
 *
 * @param[in] fpmtflag Failed PMT flag for D2 modules.
 *
 * Set the failed PMT flag for D2 modules. The following values can be
 * set:
 * - 0: excluded D2 modules with failed PMTs
 * - 1: include D2 modules with failed PMTs
 * - 2: include D2 modules with failed PMTs by excluding zone around failed
 *      PMTs
 ***************************************************************************/
void GCOMSelection::fpmtflag(const int& fpmtflag)
{
    // Check validify of values
    if (fpmtflag < 0 || fpmtflag > 2) {
        std::string msg = "Invalid value "+gammalib::str(fpmtflag)+
                          " specified for failed PMT flag for D2 modules. "
                          "Please specify 0, 1 or 2.";
        throw GException::invalid_argument(G_FPMTFLAG, msg);
    }

    // Set flag
    m_fpmtflag = fpmtflag;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Return D1 module usage flag
 *
 * @param[in] id1 D1 module identifier [0,...,6].
 * @return True if D1 module @p id1 should be used.
 ***************************************************************************/
const bool& GCOMSelection::use_d1(const int& id1) const
{
    // Check module validity range
    if (id1 < 0 or id1 > 6) {
        std::string msg = "Invalid D1 module identifier "+gammalib::str(id1)+
                          " specified. The D1 module identifier needs to be "
                          "comprised within 0 and 6.";
        throw GException::invalid_argument(G_USE_D1_GET, msg);
    }

    // Return usage flag
    return (m_use_d1[id1]);
}


/***********************************************************************//**
 * @brief Set D1 module usage flag
 *
 * @param[in] id1 D1 module identifier [0,...,6].
 * @param[in] use D1 module usage.
 ***************************************************************************/
void GCOMSelection::use_d1(const int& id1, const bool& use)
{
    // Check module validity range
    if (id1 < 0 or id1 > 6) {
        std::string msg = "Invalid D1 module identifier "+gammalib::str(id1)+
                          " specified. The D1 module identifier needs to be "
                          "comprised within 0 and 6.";
        throw GException::invalid_argument(G_USE_D1_SET, msg);
    }

    // Set usage flag
    m_use_d1[id1] = use;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Return D2 module usage flag
 *
 * @param[in] id2 D2 module identifier [0,...,13].
 * @return True if D2 module @p id2 should be used.
 ***************************************************************************/
const bool& GCOMSelection::use_d2(const int& id2) const
{
    // Check module validity range
    if (id2 < 0 or id2 > 13) {
        std::string msg = "Invalid D2 module identifier "+gammalib::str(id2)+
                          " specified. The D2 module identifier needs to be "
                          "comprised within 0 and 13.";
        throw GException::invalid_argument(G_USE_D2_GET, msg);
    }

    // Return usage flag
    return (m_use_d2[id2]);
}


/***********************************************************************//**
 * @brief Set D2 module usage flag
 *
 * @param[in] id2 D2 module identifier [0,...,13].
 * @param[in] use D2 module usage.
 ***************************************************************************/
void GCOMSelection::use_d2(const int& id2, const bool& use)
{
    // Check module validity range
    if (id2 < 0 or id2 > 13) {
        std::string msg = "Invalid D2 module identifier "+gammalib::str(id2)+
                          " specified. The D2 module identifier needs to be "
                          "comprised within 0 and 13.";
        throw GException::invalid_argument(G_USE_D2_SET, msg);
    }

    // Set usage flag
    m_use_d2[id2] = use;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Read selection set from FITS HDU keywords
 *
 * @param[in] hdu FITS HDU.
 ***************************************************************************/
void GCOMSelection::read(const GFitsHDU& hdu)
{
    // Read D1 energy selection
    if (hdu.has_card("D1EMIN") && hdu.has_card("D1EMAX")) {
        m_e1_min = hdu.real("D1EMIN");
        m_e1_max = hdu.real("D1EMAX");
    }

    // Read D2 energy selection
    if (hdu.has_card("D2EMIN") && hdu.has_card("D2EMAX")) {
        m_e2_min = hdu.real("D2EMIN");
        m_e2_max = hdu.real("D2EMAX");
    }

    // Read ToF selection
    if (hdu.has_card("TOFMIN") && hdu.has_card("TOFMAX")) {
        m_tof_min = hdu.real("TOFMIN");
        m_tof_max = hdu.real("TOFMAX");
    }

    // Read PSD selection
    if (hdu.has_card("PSDMIN") && hdu.has_card("PSDMAX")) {
        m_psd_min = hdu.real("PSDMIN");
        m_psd_max = hdu.real("PSDMAX");
    }

    // Read rejection flag selection
    if (hdu.has_card("RFLMIN") && hdu.has_card("RFLMAX")) {
        m_reflag_min = hdu.integer("RFLMIN");
        m_reflag_max = hdu.integer("RFLMAX");
    }

    // Read veto flag selection
    if (hdu.has_card("VFLMIN") && hdu.has_card("VFLMAX")) {
        m_vetoflag_min = hdu.integer("VFLMIN");
        m_vetoflag_max = hdu.integer("VFLMAX");
    }

    // Read D2 PMT failure handling flag
    if (hdu.has_card("D2FPMT")) {
        m_fpmtflag = hdu.integer("D2FPMT");
    }

    // Read D1 module usage flag
    if (hdu.has_card("D1USE")) {
        std::string use = gammalib::strip_whitespace(hdu.string("D1USE"));
        if (use.length() == 7) {
            for (int i = 0; i < 7; ++i) {
                m_use_d1[i] = (use[i] == '1');
            }
        }
    }

    // Read D2 module usage flag
    if (hdu.has_card("D2USE")) {
        std::string use = gammalib::strip_whitespace(hdu.string("D2USE"));
        if (use.length() == 14) {
            for (int i = 0; i < 14; ++i) {
                m_use_d2[i] = (use[i] == '1');
            }
        }
    }

    // Read selection statistics
    if (hdu.has_card("NEVCHK")) {
        m_num_events_checked = hdu.integer("NEVCHK");
    }
    if (hdu.has_card("NEVUSE")) {
        m_num_events_used = hdu.integer("NEVUSE");
    }
    if (hdu.has_card("NEVREJ")) {
        m_num_events_rejected = hdu.integer("NEVREJ");
    }
    if (hdu.has_card("ND1ELO") && hdu.has_card("ND1EHI")) {
        m_num_e1_min = hdu.integer("ND1ELO");
        m_num_e1_max = hdu.integer("ND1EHI");
    }
    if (hdu.has_card("ND2ELO") && hdu.has_card("ND2EHI")) {
        m_num_e2_min = hdu.integer("ND2ELO");
        m_num_e2_max = hdu.integer("ND2EHI");
    }
    if (hdu.has_card("NTOFLO") && hdu.has_card("NTOFHI")) {
        m_num_tof_min = hdu.integer("NTOFLO");
        m_num_tof_max = hdu.integer("NTOFHI");
    }
    if (hdu.has_card("NPSDLO") && hdu.has_card("NPSDHI")) {
        m_num_psd_min = hdu.integer("NPSDLO");
        m_num_psd_max = hdu.integer("NPSDHI");
    }
    if (hdu.has_card("NRFLLO") && hdu.has_card("NRFLHI")) {
        m_num_reflag_min = hdu.integer("NRFLLO");
        m_num_reflag_max = hdu.integer("NRFLHI");
    }
    if (hdu.has_card("NVFLLO") && hdu.has_card("NVFLHI")) {
        m_num_vetoflag_min = hdu.integer("NVFLLO");
        m_num_vetoflag_max = hdu.integer("NVFLHI");
    }
    if (hdu.has_card("NNOSCT")) {
        m_num_no_scatter = hdu.integer("NNOSCT");
    }
    for (int i = 0; i < 7; ++i) {
        std::string key = "ND1-"+gammalib::str(i+1,"%2.2d");
        if (hdu.has_card(key)) {
             m_num_d1[i] = hdu.integer(key);
        }
    }
    for (int i = 0; i < 14; ++i) {
        std::string key = "ND2-"+gammalib::str(i+1,"%2.2d");
        if (hdu.has_card(key)) {
             m_num_d2[i] = hdu.integer(key);
        }
    }
    if (hdu.has_card("ND1OFF")) {
        m_num_d1module_off = hdu.integer("ND1OFF");
    }
    if (hdu.has_card("ND2OFF")) {
        m_num_d2module_off = hdu.integer("ND2OFF");
    }
    if (hdu.has_card("ND2FPM")) {
        m_num_fpmt = hdu.integer("ND2FPM");
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write selection set keywords into FITS HDU
 *
 * @param[in] hdu FITS HDU.
 ***************************************************************************/
void GCOMSelection::write(GFitsHDU& hdu) const
{
    // Write D1 energy selection
    hdu.card("D1EMIN",  m_e1_min, "[MeV] Minimum D1 energy deposit");
    hdu.card("D1EMAX",  m_e1_max, "[MeV] Maximum D1 energy deposit");

    // Write D2 energy selection
    hdu.card("D2EMIN",  m_e2_min, "[MeV] Minimum D2 energy deposit");
    hdu.card("D2EMAX",  m_e2_max, "[MeV] Maximum D2 energy deposit");

    // Write ToF selection
    hdu.card("TOFMIN",  m_tof_min, "ToF selection minimum");
    hdu.card("TOFMAX",  m_tof_max, "ToF selection maximum");

    // Write PSD selection
    hdu.card("PSDMIN",  m_psd_min, "PSD selection minimum");
    hdu.card("PSDMAX",  m_psd_max, "PSD selection maximum");

    // Write rejection flag selection
    hdu.card("RFLMIN",  m_reflag_min, "Rejection flag minimum");
    hdu.card("RFLMAX",  m_reflag_max, "Rejection flag maximum");

    // Write veto flag selection
    hdu.card("VFLMIN",  m_vetoflag_min, "Veto flag minimum");
    hdu.card("VFLMAX",  m_vetoflag_max, "Veto flag maximum");

    // Write D2 PMT failure handling flag
    hdu.card("D2FPMT",  m_fpmtflag, "D2 PMT failure handling");

    // Write D1 module usage
    std::string d1use = "0000000";
    for (int i = 0; i < 7; ++i) {
        if (m_use_d1[i]) {
            d1use[i] = '1';
        }
    }
    hdu.card("D1USE", d1use, "D1 module usage");

    // Write D2 module usage
    std::string d2use = "00000000000000";
    for (int i = 0; i < 14; ++i) {
        if (m_use_d2[i]) {
            d2use[i] = '1';
        }
    }
    hdu.card("D2USE", d2use, "D2 module usage");

    // Write selection statistics
    hdu.card("NEVCHK",  m_num_events_checked,  "Number of checked events");
    hdu.card("NEVUSE",  m_num_events_used,     "Number of used events");
    hdu.card("NEVREJ",  m_num_events_rejected, "Number of rejected events");
    hdu.card("ND1ELO",  m_num_e1_min,          "Number of events < D1EMIN");
    hdu.card("ND1EHI",  m_num_e1_max,          "Number of events > D1EMAX");
    hdu.card("ND2ELO",  m_num_e2_min,          "Number of events < D2EMIN");
    hdu.card("ND2EHI",  m_num_e2_max,          "Number of events > D2EMAX");
    hdu.card("NTOFLO",  m_num_tof_min,         "Number of events < TOFMIN");
    hdu.card("NTOFHI",  m_num_tof_max,         "Number of events > TOFMIN");
    hdu.card("NPSDLO",  m_num_psd_min,         "Number of events < PSDMIN");
    hdu.card("NPSDHI",  m_num_psd_max,         "Number of events > PSDMAX");
    hdu.card("NRFLLO",  m_num_reflag_min,      "Number of events < RFLMIN");
    hdu.card("NRFLHI",  m_num_reflag_max,      "Number of events > RFLMAX");
    hdu.card("NVFLLO",  m_num_vetoflag_min,    "Number of events < VFLMIN");
    hdu.card("NVFLHI",  m_num_vetoflag_max,    "Number of events > VFLMAX");
    hdu.card("NNOSCT",  m_num_no_scatter,      "Number of events w/o scatter angle");
    hdu.card("NINVMT",  m_num_invalid_modcom,  "Number of events with invalid minitelescope");
    for (int i = 0; i < 7; ++i) {
        std::string key = "ND1-"+gammalib::str(i+1,"%2.2d");
        std::string com = "Number of events in D1 module "+gammalib::str(i+1);
        hdu.card(key, m_num_d1[i], com);
    }
    for (int i = 0; i < 14; ++i) {
        std::string key = "ND2-"+gammalib::str(i+1,"%2.2d");
        std::string com = "Number of events in D2 module "+gammalib::str(i+1);
        hdu.card(key, m_num_d2[i], com);
    }
    hdu.card("ND1OFF",  m_num_d1module_off,    "Number of events with D1 module off");
    hdu.card("ND2OFF",  m_num_d2module_off,    "Number of events with D2 module off");
    hdu.card("ND2FPM",  m_num_fpmt,            "Number of events excluded due to failed PMTs");

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set orbital period
 *
 * @param[in] period Orbital period (days).
 * @param[in] time Time of phase zero.
 *
 * Set the orbital phase for orbital phase selection.
 ***************************************************************************/
void GCOMSelection::orbital_period(const double& period, const GTime& time)
{
    // Clear temporal phase curve
    m_orbital_phase_curve.clear();

    // Get frequency of orbit in units of Hz
    double f0 = 1.0 / (period * gammalib::sec_in_day);

    // Set phase curve elements
    m_orbital_phase_curve.mjd(time);   // Reference time
    m_orbital_phase_curve.phase(0.0);  // Phase at reference time
    m_orbital_phase_curve.f0(f0);      // Frequency at reference time
    m_orbital_phase_curve.f1(0.0);     // No frequency derivative
    m_orbital_phase_curve.f2(0.0);     // No 2nd frequency derivative

    // Return
    return;
}


/***********************************************************************//**
 * @brief Print COMPTEL selection set
 *
 * @param[in] chatter Chattiness.
 * @return String containing COMPTEL selection set information.
 ***************************************************************************/
std::string GCOMSelection::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GCOMSelection ===");

        // Append number of checked events
        result.append("\n"+gammalib::parformat("Checked events"));
        result.append(gammalib::str(m_num_events_checked));
        result.append("\n"+gammalib::parformat("Accepted events"));
        result.append(gammalib::str(m_num_events_used));
        result.append("\n"+gammalib::parformat("Rejected events"));
        result.append(gammalib::str(m_num_events_rejected));

        // Append E1 selection statistics
        result.append("\n"+gammalib::parformat("E1 selection"));
        result.append(gammalib::str(m_num_e1_min)+" < [");
        result.append(gammalib::str(m_e1_min)+" - ");
        result.append(gammalib::str(m_e1_max)+" MeV] < ");
        result.append(gammalib::str(m_num_e1_max));

        // Append E2 selection statistics
        result.append("\n"+gammalib::parformat("E2 selection"));
        result.append(gammalib::str(m_num_e2_min)+" < [");
        result.append(gammalib::str(m_e2_min)+" - ");
        result.append(gammalib::str(m_e2_max)+" MeV] < ");
        result.append(gammalib::str(m_num_e2_max));

        // Append TOF selection statistics
        result.append("\n"+gammalib::parformat("TOF selection"));
        result.append(gammalib::str(m_num_tof_min)+" < [");
        result.append(gammalib::str(m_tof_min)+" - ");
        result.append(gammalib::str(m_tof_max)+"] < ");
        result.append(gammalib::str(m_num_tof_max));

        // Append PSD selection statistics
        result.append("\n"+gammalib::parformat("PSD selection"));
        result.append(gammalib::str(m_num_psd_min)+" < [");
        result.append(gammalib::str(m_psd_min)+" - ");
        result.append(gammalib::str(m_psd_max)+"] < ");
        result.append(gammalib::str(m_num_psd_max));

        // Append rejection flag selection statistics
        result.append("\n"+gammalib::parformat("Rejection flag selection"));
        result.append(gammalib::str(m_num_reflag_min)+" < [");
        result.append(gammalib::str(m_reflag_min)+" - ");
        result.append(gammalib::str(m_reflag_max)+"] < ");
        result.append(gammalib::str(m_num_reflag_max));

        // Append veto flag selection statistics
        result.append("\n"+gammalib::parformat("Veto flag selection"));
        result.append(gammalib::str(m_num_vetoflag_min)+" < [");
        result.append(gammalib::str(m_vetoflag_min)+" - ");
        result.append(gammalib::str(m_vetoflag_max)+"] < ");
        result.append(gammalib::str(m_num_vetoflag_max));

        // Append D1 module usage
        for (int i = 0; i < 7; ++i) {
            result.append("\n"+gammalib::parformat("D1 module "+gammalib::str(i+1)));
            result.append(gammalib::str(m_num_d1[i]));
            if (m_use_d1[i]) {
                result.append(" [use]");
            }
            else {
                result.append(" [don't use]");
            }
        }

        // Append D2 module usage
        for (int i = 0; i < 14; ++i) {
            result.append("\n"+gammalib::parformat("D2 module "+gammalib::str(i+1)));
            result.append(gammalib::str(m_num_d2[i]));
            if (m_use_d2[i]) {
                result.append(" [use]");
            }
            else {
                result.append(" [don't use]");
            }
        }

        // Append D1 module off
        result.append("\n"+gammalib::parformat("D1 modules off"));
        result.append(gammalib::str(m_num_d1module_off));

        // Append D2 module off
        result.append("\n"+gammalib::parformat("D2 modules off"));
        result.append(gammalib::str(m_num_d2module_off));

        // Append D2 PMT failure handling statistics
        result.append("\n"+gammalib::parformat("D2 PMT failure"));
        result.append(gammalib::str(m_num_fpmt)+" [");
        result.append(gammalib::str(m_fpmtflag)+"]");

        // Append other statistics
        result.append("\n"+gammalib::parformat("No scatter angle"));
        result.append(gammalib::str(m_num_no_scatter));
        result.append("\n"+gammalib::parformat("Invalid minitelescope"));
        result.append(gammalib::str(m_num_invalid_modcom));

        // Append orbital phase selection if it exists
        if (!m_orbital_phases.is_empty()) {
            result.append("\n"+m_orbital_phases.print());
            result.append("\n"+m_orbital_phase_curve.print());
        }

        // Append pulsar phase selection if it exists
        if (!m_pulsar_phases.is_empty()) {
            result.append("\n"+m_pulsar_phases.print());
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
void GCOMSelection::init_members(void)
{
    // Initialise members
    m_e1_min       = 0.070;  //!< Minimum D1 energy deposit (MeV)
    m_e1_max       =  20.0;  //!< Maximum D1 energy deposit (MeV)
    m_e2_min       = 0.650;  //!< Minimum D2 energy deposit (MeV)
    m_e2_max       =  30.0;  //!< Maximum D2 energy deposit (MeV)
    m_tof_min      =   115;  //!< Minimum TOF window
    m_tof_max      =   130;  //!< Maximum TOF window
    m_psd_min      =     0;  //!< Minimum PSD window
    m_psd_max      =   110;  //!< Maximum PSD window
    m_reflag_min   =     1;  //!< Minimum rejection flag
    m_reflag_max   =  1000;  //!< Maximum rejection flag
    m_vetoflag_min =     0;  //!< Minimum veto flag
    m_vetoflag_max =     0;  //!< Maximum veto flag
    m_fpmtflag     =     0;  //!< D2 PMT failure flag
    m_orbital_phases.clear();
    m_pulsar_phases.clear();
    m_orbital_phase_curve.clear();
    for (int i = 0; i < 7; ++i) {
        m_use_d1[i] = true;
    }
    for (int i = 0; i < 14; ++i) {
        m_use_d2[i] = true;
    }

    // Initialise statistics
    init_statistics();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] select COMPTEL selection set.
 ***************************************************************************/
void GCOMSelection::copy_members(const GCOMSelection& select)
{
    // Copy members
    m_e1_min              = select.m_e1_min;
    m_e1_max              = select.m_e1_max;
    m_e2_min              = select.m_e2_min;
    m_e2_max              = select.m_e2_max;
    m_tof_min             = select.m_tof_min;
    m_tof_max             = select.m_tof_max;
    m_psd_min             = select.m_psd_min;
    m_psd_max             = select.m_psd_max;
    m_reflag_min          = select.m_reflag_min;
    m_reflag_max          = select.m_reflag_max;
    m_vetoflag_min        = select.m_vetoflag_min;
    m_vetoflag_max        = select.m_vetoflag_max;
    m_fpmtflag            = select.m_fpmtflag;
    m_orbital_phases      = select.m_orbital_phases;
    m_pulsar_phases       = select.m_pulsar_phases;
    m_orbital_phase_curve = select.m_orbital_phase_curve;
    for (int i = 0; i < 7; ++i) {
        m_use_d1[i] = select.m_use_d1[i];
    }
    for (int i = 0; i < 14; ++i) {
        m_use_d2[i] = select.m_use_d2[i];
    }

    // Copy statistics
    m_num_events_checked  = select.m_num_events_checked;
    m_num_events_used     = select.m_num_events_used;
    m_num_events_rejected = select.m_num_events_rejected;
    m_num_e1_min          = select.m_num_e1_min;
    m_num_e1_max          = select.m_num_e1_max;
    m_num_e2_min          = select.m_num_e2_min;
    m_num_e2_max          = select.m_num_e2_max;
    m_num_tof_min         = select.m_num_tof_min;
    m_num_tof_max         = select.m_num_tof_max;
    m_num_psd_min         = select.m_num_psd_min;
    m_num_psd_max         = select.m_num_psd_max;
    m_num_reflag_min      = select.m_num_reflag_min;
    m_num_reflag_max      = select.m_num_reflag_max;
    m_num_vetoflag_min    = select.m_num_vetoflag_min;
    m_num_vetoflag_max    = select.m_num_vetoflag_max;
    m_num_no_scatter      = select.m_num_no_scatter;
    m_num_invalid_modcom  = select.m_num_invalid_modcom;
    m_num_d1module_off    = select.m_num_d1module_off;
    m_num_d2module_off    = select.m_num_d2module_off;
    m_num_fpmt            = select.m_num_fpmt;

    // Copy module usage
    for (int i = 0; i < 7; ++i) {
        m_num_d1[i] = select.m_num_d1[i];
    }
    for (int i = 0; i < 14; ++i) {
        m_num_d2[i] = select.m_num_d2[i];
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GCOMSelection::free_members(void)
{
    // Return
    return;
}
