/***************************************************************************
 *              GPulsarEphemeris.cpp - Pulsar ephemeris class              *
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
 * @file GPulsarEphemeris.cpp
 * @brief Pulsar ephemeris class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GPulsarEphemeris.hpp"
#include "GTools.hpp"

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
GPulsarEphemeris::GPulsarEphemeris(void)
{
    // Initialise class members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] ephemeris Pulsar ephemeris.
 ***************************************************************************/
GPulsarEphemeris::GPulsarEphemeris(const GPulsarEphemeris& ephemeris)
{
    // Initialise class members
    init_members();

    // Copy members
    copy_members(ephemeris);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GPulsarEphemeris::~GPulsarEphemeris(void)
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
 * @param[in] ephemeris Pulsar ephemeris.
 * @return Pulsar ephemeris.
 ***************************************************************************/
GPulsarEphemeris& GPulsarEphemeris::operator=(const GPulsarEphemeris& ephemeris)
{
    // Execute only if object is not identical
    if (this != &ephemeris) {

        // Free members
        free_members();

        // Initialise private members
        init_members();

        // Copy members
        copy_members(ephemeris);

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
 * @brief Clear Pulsar ephemeris
 ***************************************************************************/
void GPulsarEphemeris::clear(void)
{
    // Free members
    free_members();

    // Initialise private members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone Pulsar ephemeris
 *
 * @return Pointer to deep copy of Pulsar ephemeris.
 ***************************************************************************/
GPulsarEphemeris* GPulsarEphemeris::clone(void) const
{
    return new GPulsarEphemeris(*this);
}


/***********************************************************************//**
 * @brief Returns pulsar phase
 *
 * @param[in] time Time.
 * @param[in] timesys Time system for evaluation of time.
 * @return Pulsar phase.
 *
 * Computes
 *
 * \f[
 *    \Phi(t) = \Phi_0 + f(t-t_0) + \frac{1}{2}\dot{f} (t-t_0)^2 +
 *                                  \frac{1}{6}\ddot{f} (t-t_0)^3
 * \f]
 *
 * where
 * \f$t\f$ is the specified @p time in seconds,
 * \f$t_0\f$ is a reference time in seconds,
 * \f$\Phi_0\f$ is the phase at the reference time,
 * \f$f\f$ is the variation frequency at the reference time,
 * \f$\dot{f}\f$ is the first derivative of the variation frequency at the
 * reference time, and
 * \f$\dot{f}\f$ is the second derivative of the variation frequency at the
 * reference time.
 *
 * The phase \f$\Phi(t)\f$ is in the interval [0,1].
 ***************************************************************************/
double GPulsarEphemeris::phase(const GTime&       time,
                               const std::string& timesys) const
{
    // Set constants
    const double c1 = 0.5;
    const double c2 = 1.0 / 6.0;

    // Compute time since reference time in seconds
    double dt = time.secs(timesys) - m_t0.secs(m_timesys);

    // Computes phase
    double phase = m_phase + dt * (m_f0 + dt * (c1 * m_f1 + c2 * m_f2 * dt));

    // Put phase into interval [0,1]
    phase -= floor(phase);

    // Return phase
    return phase;
}


/***********************************************************************//**
 * @brief Print Pulsar ephemeris
 *
 * @param[in] chatter Chattiness.
 * @return String containing Pulsar ephemeris information.
 ***************************************************************************/
std::string GPulsarEphemeris::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GPulsarEphemeris ===");

        // Append pulsar name and direction
        result.append("\n"+gammalib::parformat("Pulsar name"));
        result.append(m_name);
        result.append("\n"+gammalib::parformat("Pulsar Right Ascension"));
        result.append(gammalib::str(m_dir.ra_deg()));
        result.append(" deg");
        result.append("\n"+gammalib::parformat("Pulsar Declination"));
        result.append(gammalib::str(m_dir.dec_deg()));
        result.append(" deg");

        // Get phase information
        double f0 = this->f0();

        // Append phase information
        result.append("\n"+gammalib::parformat("Time system of ephemeris"));
        result.append(m_timesys);
        result.append("\n"+gammalib::parformat("Reference epoch"));
        result.append("MJD ");
        result.append(gammalib::str(t0().mjd(m_timesys)));
        result.append("\n"+gammalib::parformat("Phase"));
        result.append(gammalib::str(phase()));
        result.append("\n"+gammalib::parformat("Frequency"));
        result.append(gammalib::str(f0));
        result.append(" Hz");
        result.append("\n"+gammalib::parformat("Frequency derivative"));
        result.append(gammalib::str(f1()));
        result.append(" Hz/s");
        result.append("\n"+gammalib::parformat("2nd frequency derivative"));
        result.append(gammalib::str(f2()));
        result.append(" Hz/s^2");
        result.append("\n"+gammalib::parformat("Period"));
        if (f0 != 0.0) {
            double p0 = 1.0 / f0;
            if (p0 < 1.0) {
                result.append(gammalib::str(p0*1000.0));
                result.append(" ms");
            }
            else {
                result.append(gammalib::str(p0));
                result.append(" s");
            }
        }
        else {
            result.append("infinity");
        }

        // Append validity information
        result.append("\n"+gammalib::parformat("Validity MJD range"));
        result.append(gammalib::str(tstart().mjd(m_timesys)));
        result.append(" - ");
        result.append(gammalib::str(tstop().mjd(m_timesys)));
        result.append("\n"+gammalib::parformat("Validity UTC range"));
        result.append(tstart().utc());
        result.append(" - ");
        result.append(tstop().utc());

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
void GPulsarEphemeris::init_members(void)
{
    // Initialise members
    m_name.clear();
    m_tstart.clear();
    m_tstop.clear();
    m_dir.clear();
    m_t0.clear();
    m_timesys = "TT";
    m_phase   = 0.0;
    m_f0      = 0.0;
    m_f1      = 0.0;
    m_f2      = 0.0;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] ephemeris Pulsar ephemeris.
 ***************************************************************************/
void GPulsarEphemeris::copy_members(const GPulsarEphemeris& ephemeris)
{
    // Copy members
    m_name    = ephemeris.m_name;
    m_tstart  = ephemeris.m_tstart;
    m_tstop   = ephemeris.m_tstop;
    m_dir     = ephemeris.m_dir;
    m_timesys = ephemeris.m_timesys;
    m_t0      = ephemeris.m_t0;
    m_phase   = ephemeris.m_phase;
    m_f0      = ephemeris.m_f0;
    m_f1      = ephemeris.m_f1;
    m_f2      = ephemeris.m_f2;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GPulsarEphemeris::free_members(void)
{
    // Return
    return;
}
