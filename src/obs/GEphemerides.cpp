/***************************************************************************
 *                   GEphemerides.cpp - Ephemerides class                  *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2022 by Juergen Knoedlseder                              *
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
 * @file GEphemerides.cpp
 * @brief Ephemerides class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GEphemerides.hpp"
#include "GFits.hpp"
#include "GFitsTable.hpp"
#include "GFitsTableCol.hpp"
#include "GSkyDir.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_EPHEMERIS    "GEphemerides::ephemeris(GTime&, GVector*, GVector*, "\
                                                         "GVector*, double*)"
#define G_FETCH_DATA                             "GEphemerides::fetch_data()"

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
GEphemerides::GEphemerides(void)
{
    // Initialise class members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] ephemerides Ephemerides.
 ***************************************************************************/
GEphemerides::GEphemerides(const GEphemerides& ephemerides)
{
    // Initialise class members
    init_members();

    // Copy members
    copy_members(ephemerides);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GEphemerides::~GEphemerides(void)
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
 * @param[in] ephemerides Ephemerides.
 * @return Ephemerides.
 ***************************************************************************/
GEphemerides& GEphemerides::operator=(const GEphemerides& ephemerides)
{
    // Execute only if object is not identical
    if (this != &ephemerides) {

        // Free members
        free_members();

        // Initialise private members
        init_members();

        // Copy members
        copy_members(ephemerides);

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
 * @brief Clear Ephemerides
 ***************************************************************************/
void GEphemerides::clear(void)
{
    // Free members
    free_members();

    // Initialise private members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone Ephemerides
 *
 * @return Pointer to deep copy of Ephemerides.
 ***************************************************************************/
GEphemerides* GEphemerides::clone(void) const
{
    return new GEphemerides(*this);
}


/***********************************************************************//**
 * @brief Load Ephemerides
 *
 * @param[in] filename Ephemerides file name.
 *
 * Load ephemerides from FITS file.
 ***************************************************************************/
void GEphemerides::load(const GFilename& filename)
{
    // Clear members
    clear();

    // Open JPL ephemerides FITS file
    GFits fits(filename);

    // Get ephemerides table
    const GFitsTable* ephem = fits.table("EPHEM");

    // Extract number of events in FITS file
    int num = ephem->nrows();

    // Continue if there are ephemerides
    if (num > 0) {

        // Reserve data
        m_times.reserve(num);
        m_earth.reserve(num);
        m_earth_dt.reserve(num);
        m_earth_d2t.reserve(num);
        m_earth_d3t.reserve(num);
        m_sun.reserve(num);
        m_tdb2tt.reserve(num);

        // Get validity interval
        double tstart = ephem->real("TSTART");
        m_tstart.jd(tstart);
        m_tstop.jd(ephem->real("TSTOP"));

        // Get column pointers
        const GFitsTableCol* ptr_earth  = (*ephem)["EARTH"];    // light seconds
        const GFitsTableCol* ptr_sun    = (*ephem)["SUN"];      // light seconds
        const GFitsTableCol* ptr_tbd2tt = (*ephem)["TIMEDIFF"]; // sec

        // Extract data
        for (int i = 0; i < num; ++i) {

            // Set time
            GTime time;
            time.jd(tstart + double(i));

            // Set vectors
            GVector earth(ptr_earth->real(i,0),
                          ptr_earth->real(i,1),
                          ptr_earth->real(i,2));
            GVector earth_dt(ptr_earth->real(i,3),
                             ptr_earth->real(i,4),
                             ptr_earth->real(i,5));
            GVector earth_d2t(ptr_earth->real(i,6),
                              ptr_earth->real(i,7),
                              ptr_earth->real(i,8));
            GVector earth_d3t(ptr_earth->real(i,9),
                              ptr_earth->real(i,10),
                              ptr_earth->real(i,11));
            GVector sun(ptr_sun->real(i,0),
                        ptr_sun->real(i,1),
                        ptr_sun->real(i,2));

            // Set tbd2tt
            double tbd2tt = ptr_tbd2tt->real(i);

            // Push back data
            m_times.push_back(time);
            m_earth.push_back(earth);
            m_earth_dt.push_back(earth_dt);
            m_earth_d2t.push_back(earth_d2t);
            m_earth_d3t.push_back(earth_d3t);
            m_sun.push_back(sun);
            m_tdb2tt.push_back(tbd2tt);

        } // endfor: looped over rows

        // Store attributes
        m_name     = ephem->string("NAME");
        m_filename = filename;

    } // endif: there were ephemerides to load

    // Return
    return;
}


/***********************************************************************//**
 * @brief Get ephemeris vector and TBD->TT value for a given time
 *
 * @param[in] time Time.
 * @param[out] rce Pointer to vector from SSBC to Earth (light-s).
 * @param[out] rcs Pointer to vector from SSBC to Sun (light-s).
 * @param[out] vce Pointer to first time derivative of vector from SSBC to
 *                 Earth (light-s/s).
 * @param[out] etut Time difference TDB-TT (s)
 *
 * Get ephemeris vector and TBD-TT value for a given time. Information is
 * only returned for pointers that are not NULL.
 *
 * The code was inspired from the ftools routine xtereadeph.f and the COMPASS
 * routine bvceph.f.
 ***************************************************************************/
void GEphemerides::ephemeris(const GTime& time,
                             GVector*     rce,
                             GVector*     rcs,
                             GVector*     vce,
                             double*      etut) const
{
    // Fetch ephemerides
    const_cast<GEphemerides*>(this)->fetch_data();

    // Compute jd_utc as int
    double jd     = time.jd("UTC");
    int    jd_utc = int(jd);

    // Compute dt in TT as day fraction, between -0.5 and 0.5. While the
    // time system in the xtereadeph.f routine is not specified, the
    // bvceph.f routine explicitly adds the UTC->TT conversion coefficient
    // to the delta
    //double dt     = jd - double(jd_utc);
    double dt = time.jd("TT") - double(int(time.jd("TT")));
    if (dt > 0.5) {
        jd_utc += 1;
        dt     -= 1.0;
    }

    // Get Julian Day validity interval as int
    int jd0 = int(m_tstart.jd());
    int jd1 = int(m_tstop.jd());

    // If date is outside range of ephemeris then throw an exception
    if (jd_utc < jd0) {
        std::string msg = "Time JD "+gammalib::str(jd_utc)+ " is before "
                          "validity start JD "+gammalib::str(jd0)+" of "
                          "ephemerides. Please specify a time later than "
                          "the validity start.";
        throw GException::invalid_argument(G_EPHEMERIS, msg);
    }
    else if (jd_utc > jd1) {
        std::string msg = "Time JD "+gammalib::str(jd_utc)+ " is after "
                          "validity end JD "+gammalib::str(jd1)+" of "
                          "ephemerides. Please specify a time earlier than "
                          "the validity end.";
        throw GException::invalid_argument(G_EPHEMERIS, msg);
    }

    // Compute ephemeris index
    int index = jd_utc - jd0;

    // Make rough computation of time derivative of TDB-TDT
    double tbd2tt_dot = (index+1 < size()) ? m_tdb2tt[index+1] - m_tdb2tt[index]
                                           : m_tdb2tt[index]   - m_tdb2tt[index-1];

    // Get Earth position and velocity vector at the desired time using
    // Taylor expansion
    const double c1 = 1.0/2.0;
    const double c2 = 1.0/6.0;
    if (rce != NULL) {
        *rce = m_earth[index] + (m_earth_dt[index] + (m_earth_d2t[index] *
               c1 + m_earth_d3t[index] * dt * c2) * dt) * dt;
    }
    if (vce != NULL) {
        *vce = (m_earth_dt[index] + (m_earth_d2t[index] + m_earth_d3t[index] *
               dt * c1)) * dt * gammalib::sec2day;
    }

    // Get Sun vector
    if (rcs != NULL) {
        *rcs = m_sun[index];
    }

    // Get TBD-TT (seconds)
    if (etut != NULL) {
        *etut = m_tdb2tt[index] + dt * tbd2tt_dot;
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Get time difference between geocentric and SSB (seconds)
 *
 * @param[in] srcdir Source direction.
 * @param[in] time Geocentric time.
 * @return Time difference in seconds.
 ***************************************************************************/
double GEphemerides::geo2ssb(const GSkyDir& srcdir, const GTime& time) const
{
    // Set constants
    const double t_sun = 4.92549089483e-6;                  // s
    const double r_sun = 2.315;                             // light s
    const double inv_c = 1000.0 / gammalib::speed_of_light; // s/km

    // Get ephemerides
    GVector rce(3);
    GVector rcs(3);
    GVector vce(3);
    double  etut;
    ephemeris(time, &rce, &rcs, &vce, &etut);

    // Get sky direction as celestial vector
    GVector dir = srcdir.celvector();

    // Calculate light-travel time to barycenter in Euclidean space
    double light_travel_time = dir * rce;

    // Now calculate the time delay due to the gravitational field of the
    // Sun (I.I. Shapiro, Phys. Rev. Lett. 13, 789 (1964))
    GVector rsa       = rce - rcs;
    double  sundis    = norm(rsa);
    double  sunsiz    = r_sun/sundis;
    double  cos_theta = (dir * rsa) / sundis;

    // Special handling if sky direction is inside the Sun
    /*
    if (cth + 1.0 < 0.5 * sunsiz * sunsiz) {

    }
    */

    // Compute Shapiro delay
    double shapiro_delay = -2.0 * t_sun * std::log(1.0 + cos_theta);

    // Compute time difference between geocentric and Solar System Barycentric
    // frame
    double geo2ssb = light_travel_time - shapiro_delay + etut;

    // Return barycentric correction (seconds)
    return geo2ssb;
}


/***********************************************************************//**
 * @brief Get time difference between observatory and SSB (seconds)
 *
 * @param[in] srcdir Source direction.
 * @param[in] time Geocentric time.
 * @param[in] obs Observatory vector (km).
 * @return Time difference in seconds.
 *
 * Computes the time difference between the location of an observatory,
 * specified by an observatory position vector @p obs in units of km, and the
 * Solar System Barycentre.
 ***************************************************************************/
double GEphemerides::geo2ssb(const GSkyDir& srcdir,
                             const GTime&   time,
                             const GVector& obs) const
{
    // Set constants
    const double t_sun = 4.92549089483e-6;                  // s
    const double r_sun = 2.315;                             // light s
    const double inv_c = 1000.0 / gammalib::speed_of_light; // s/km

    // Get ephemerides
    GVector rce(3);
    GVector rcs(3);
    GVector vce(3);
    double  etut;
    ephemeris(time, &rce, &rcs, &vce, &etut);

    // Compute vector from Solar System Barycentre to observatory
    GVector rca = rce + obs * inv_c;

    // Get sky direction as celestial vector
    GVector dir = srcdir.celvector();

    // Calculate light-travel time to barycenter in Euclidean space
    double light_travel_time = dir * rca;

    // Now calculate the time delay due to the gravitational field of the
    // Sun (I.I. Shapiro, Phys. Rev. Lett. 13, 789 (1964))
    GVector rsa       = rca - rcs;
    double  sundis    = norm(rsa);
    double  sunsiz    = r_sun       / sundis;
    double  cos_theta = (dir * rsa) / sundis;

    // Special handling if sky direction is inside the Sun
    /*
    if (cth + 1.0 < 0.5 * sunsiz * sunsiz) {

    }
    */

    // Compute Shapiro delay
    double shapiro_delay = -2.0 * t_sun * std::log(1.0 + cos_theta);

    // Compute time difference between geocentric and Solar System Barycentric
    // frame
    double geo2ssb = light_travel_time - shapiro_delay + etut;

    // Return barycentric correction (seconds)
    return geo2ssb;
}


/***********************************************************************//**
 * @brief Print Ephemerides
 *
 * @param[in] chatter Chattiness.
 * @return String containing Ephemerides information.
 ***************************************************************************/
std::string GEphemerides::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GEphemerides ===");

        // Append information
        result.append("\n"+gammalib::parformat("Ephemerides name"));
        result.append(m_name);
        result.append("\n"+gammalib::parformat("File name"));
        result.append(m_filename.url());
        if (!is_empty()) {
            result.append("\n"+gammalib::parformat("Validity MJD range"));
            result.append(gammalib::str(m_tstart.mjd()));
            result.append(" - ");
            result.append(gammalib::str(m_tstop.mjd()));
            result.append("\n"+gammalib::parformat("Validity UTC range"));
            result.append(m_tstart.utc());
            result.append(" - ");
            result.append(m_tstop.utc());
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
void GEphemerides::init_members(void)
{
    // Initialise members
    m_name.clear();
    m_filename.clear();
    m_tstart.clear();
    m_tstop.clear();
    m_times.clear();
    m_earth.clear();
    m_earth_dt.clear();
    m_earth_d2t.clear();
    m_earth_d3t.clear();
    m_sun.clear();
    m_tdb2tt.clear();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] ephemerides Ephemerides.
 ***************************************************************************/
void GEphemerides::copy_members(const GEphemerides& ephemerides)
{
    // Copy members
    m_name      = ephemerides.m_name;
    m_filename  = ephemerides.m_filename;
    m_tstart    = ephemerides.m_tstart;
    m_tstop     = ephemerides.m_tstop;
    m_times     = ephemerides.m_times;
    m_earth     = ephemerides.m_earth;
    m_earth_dt  = ephemerides.m_earth_dt;
    m_earth_d2t = ephemerides.m_earth_d2t;
    m_earth_d3t = ephemerides.m_earth_d3t;
    m_sun       = ephemerides.m_sun;
    m_tdb2tt    = ephemerides.m_tdb2tt;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GEphemerides::free_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Fetch ephemerides data
 *
 * @exception GException::file_error
 *            Ephemerides data file not found.
 *
 * This method loads the DE200 ephemerides data from the file provided in
 * the reference data repository if no ephemerides are loaded.
 ***************************************************************************/
void GEphemerides::fetch_data(void)
{
    // If no ephemerides data are loaded then load the JPL DE200 ephemerides
    if (is_empty()) {

        // Initialise filename
        GFilename filename("$GAMMALIB/share/refdata/ephem_jpl_de200.fits");

        // If file exists then load ephemerides
        if (filename.is_fits()) {
            load(filename);
        }

        // ... otherwise the file may not yet be installed, possibly because
        // a unit test is performed. Therefore try to get file in the source
        // directory.
        else {
        
            // Initialise filename based on source directory
            GFilename src_filename("$TEST_SRCDIR/refdata/ephem_jpl_de200.fits");

            // If file exists then load ephemerides
            if (src_filename.is_fits()) {
                load(src_filename);
            }

            // ... otherwise throw an exception
            else {
                std::string msg = "Could not find ephemerides file in \""+
                                  filename.url()+"\" or in \""+
                                  src_filename.url()+"\". Unable to return "
                                  "information based on ephemerides.";
                throw GException::file_error(G_FETCH_DATA, msg);
            }

        } // endelse: check for file in source directory

    } // endif: no data were loaded

    // Return
    return;
}
