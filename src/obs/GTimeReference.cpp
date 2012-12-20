/***************************************************************************
 *                 GTimeReference.cpp - Time reference class               *
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
 * @file GTimeReference.cpp
 * @brief Time reference class interface implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GTimeReference.hpp"
#include "GTools.hpp"
#include "GException.hpp"

/* __ Constants __________________________________________________________ */

/* __ Method name definitions ____________________________________________ */
#define G_READ                              "GTimeReference::read(GFitsHDU*)"
#define G_SET     "GTimeReference::set(double&, std::string&, std::string&, "\
                                                              "std::string&)"

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
GTimeReference::GTimeReference(void)
{
    // Initialise private members
    init_members();
  
    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] ref Time reference.
 ***************************************************************************/
GTimeReference::GTimeReference(const GTimeReference& ref)
{ 
    // Initialise private members
    init_members();

    // Copy members
    copy_members(ref);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Time reference constructor
 *
 * @param[in] mjdref Reference MJD (days).
 * @param[in] timeunit Time unit ("sec(s)", "day(s)").
 * @param[in] timesys Time system (ignored so far).
 * @param[in] timeref Time reference (ignored so far).
 *
 * Sets the time reference from a MJD reference day, a time unit, a time
 * system and a time reference.
 ***************************************************************************/
GTimeReference::GTimeReference(const double&      mjdref,
                               const std::string& timeunit,
                               const std::string& timesys,
                               const std::string& timeref)
{
    // Initialise private members
    init_members();

    // Set time reference
    set(mjdref, timeunit, timesys, timeref);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Time reference constructor
 *
 * @param[in] mjdrefi Integer part of reference MJD (days).
 * @param[in] mjdreff Fractional part of reference MJD (days).
 * @param[in] timeunit Time unit (sec, days).
 * @param[in] timesys Time system (TT).
 * @param[in] timeref Local time reference.
 *
 * Sets the time reference from a MJD reference day (specified by an integer
 * and a fractional part), a time unit, a time system and a time reference.
 ***************************************************************************/
GTimeReference::GTimeReference(const int&         mjdrefi,
                               const double&      mjdreff,
                               const std::string& timeunit,
                               const std::string& timesys,
                               const std::string& timeref)
{
    // Initialise private members
    init_members();

    // Set time reference
    set(mjdrefi, mjdreff, timeunit, timesys, timeref);

    // Return
    return;
}


/***********************************************************************//**
 * @brief FITS header constructor
 *
 * @param[in] hdu FITS extension.
 *
 * Constructs time reference from the information found in the FITS header.
 * See GTimeReference::read for more information on the expected format.
 ***************************************************************************/
GTimeReference::GTimeReference(const GFitsHDU* hdu)
{
    // Initialise private members
    init_members();

    // Read reference from FITS header
    read(hdu);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GTimeReference::~GTimeReference(void)
{
    // Free members
    free_members();

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                                Operators                                =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Assignment operator
 *
 * @param[in] ref Time reference.
 ***************************************************************************/
GTimeReference& GTimeReference::operator= (const GTimeReference& ref)
{ 
    // Execute only if object is not identical
    if (this != &ref) {

        // Free members
        free_members();

        // Initialise private members
        init_members();

        // Copy members
        copy_members(ref);

    } // endif: object was not identical
  
    // Return
    return *this;
}


/*==========================================================================
 =                                                                         =
 =                             Public methods                              =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Clear time reference
 ***************************************************************************/
void GTimeReference::clear(void)
{
    // Free members
    free_members();

    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone object
 *
 * @return Pointer to deep copy of time reference.
 ***************************************************************************/
GTimeReference* GTimeReference::clone(void) const
{
    // Clone this image
    return new GTimeReference(*this);
}


/***********************************************************************//**
 * @brief Read time reference from FITS header
 *
 * @param[in] hdu FITS extension.
 *
 * GException::no_valid_time_ref
 *             No valid reference MJD found in header.
 *
 * Reads the time reference information from a FITS header. The method
 * requires either the keyword "MJDREF" or the pair of keywords "MJDREFI"
 * and "MJDREFF" to be set. The following keywords are optional (the
 * assumed default values in absent of the keywords is given in parentheses):
 *
 *     TIMEUNIT ("s")
 *     TIMESYS  ("TT")
 *     TIMEREF  ("LOCAL")
 *
 * Nothing is done if the HDU pointer is NULL.
 ***************************************************************************/
void GTimeReference::read(const GFitsHDU* hdu)
{
    // Continue only if HDU is valid
    if (hdu != NULL) {

        // Get reference MJD
        double mjdref  = (hdu->hascard("MJDREF"))  ? hdu->real("MJDREF") : 0.0;
        int    mjdrefi = (hdu->hascard("MJDREFI")) ? hdu->integer("MJDREFI") : 0;
        double mjdreff = (hdu->hascard("MJDREFF")) ? hdu->real("MJDREFF") : 0.0;

        // Get remaining keywords. To accept a large variety of FITS headers,
        // all keywords are optionally.
        std::string timeunit = (hdu->hascard("TIMEUNIT")) ? hdu->string("TIMEUNIT") : "s";
        std::string timesys  = (hdu->hascard("TIMESYS"))  ? hdu->string("TIMESYS")  : "TT";
        std::string timeref  = (hdu->hascard("TIMEREF"))  ? hdu->string("TIMEREF")  : "LOCAL";

        // Set time reference
        if (hdu->hascard("MJDREF")) {
            set(mjdref, timeunit, timesys, timeref);
        }
        else if (hdu->hascard("MJDREFI") && hdu->hascard("MJDREFF")) {
            set(mjdrefi, mjdreff, timeunit, timesys, timeref);
        }
        else {
            throw GException::no_valid_time_ref(G_READ,
                  "Require either keyword \"MJDREF\" or keyword pair"
                  " \"MJDREFI\" and \"MJDREFF\".");
        }

    } // endif: HDU was valid

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write time reference into FITS header
 *
 * @param[in] hdu FITS extension.
 *
 * Writes or updates the time reference information in a FITS header.
 * Depending of whether the keyword "MJDREF" or the pair of keywords "MJDREFI"
 * and "MJDREFF" exist already in the header, the method either writes the
 * reference MJD as floating point value, or split into an integer and a
 * fractional part. If nothing has been written yet, splitting into an
 * integer and fractional part will be used as this preserves the highest
 * possible accuracy.
 *
 * The following additional keywords are written:
 *     TIMEUNIT
 *     TIMESYS
 *     TIMEREF
 *
 * Nothing is done if the HDU pointer is NULL.
 ***************************************************************************/
void GTimeReference::write(GFitsHDU* hdu) const
{
    // Continue only if HDU is valid
    if (hdu != NULL) {

        // Case A: use floating point reference MJD
        if (hdu->hascard("MJDREF")) {
            hdu->card("MJDREF",   mjdref(),   "[days] Time reference MJD");
            hdu->card("TIMEUNIT", timeunit(), "Time unit");
            hdu->card("TIMESYS",  timesys(),  "Time system");
            hdu->card("TIMEREF",  timeref(),  "Time reference");
        }

        // Case B: use fractional reference MJD
        else {
            hdu->card("MJDREFI",  mjdrefi(),  "[days] Integer part of time reference MJD");
            hdu->card("MJDREFF",  mjdreff(),  "[days] Fractional part of time reference MJD");
            hdu->card("TIMEUNIT", timeunit(), "Time unit");
            hdu->card("TIMESYS",  timesys(),  "Time system");
            hdu->card("TIMEREF",  timeref(),  "Time reference");
        }


    } // endif: HDU was valid

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set time reference
 *
 * @param[in] mjdref Reference MJD (days).
 * @param[in] timeunit Time unit ("s", "d", "sec(s)", "day(s)").
 * @param[in] timesys Time system.
 * @param[in] timeref Time reference.
 *
 * @exception GException::time_invalid_unit
 *            Invalid time unit specified.
 *
 * Sets the time reference from a MJD reference day, a time unit, a time
 * system and a time reference.
 *
 * @todo Implement checking of "timesys" and "timeref" parameters.
 ***************************************************************************/
void GTimeReference::set(const double&      mjdref,
                         const std::string& timeunit,
                         const std::string& timesys,
                         const std::string& timeref)
{
    // Check timeunit string
    std::string ltimeunit = tolower(timeunit);
    if (ltimeunit == "d" || ltimeunit == "day" || ltimeunit == "days") {
        m_unit_sec = false;
    }
    else if (ltimeunit == "s" || ltimeunit == "sec" || ltimeunit == "secs") {
        m_unit_sec = true;
    }
    else {
        throw GException::time_invalid_unit(G_SET, timeunit,
              "Valid timeunit values are: \"d\", \"day\", \"days\","
              " \"s\", \"sec\" or \"secs\"");
    }

    // Set members
    m_mjdref   = mjdref;
    m_timeunit = ltimeunit;
    m_timesys  = toupper(timesys);
    m_timeref  = toupper(timeref);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set time reference
 *
 * @param[in] mjdrefi Integer part of reference MJD (days).
 * @param[in] mjdreff Fractional part of reference MJD (days).
 * @param[in] timeunit Time unit (sec, days).
 * @param[in] timesys Time system (TT).
 * @param[in] timeref Local time reference.
 *
 * Sets the time reference from a MJD reference day (specified by an integer
 * and a fractional part), a time unit, a time system and a time reference.
 ***************************************************************************/
void GTimeReference::set(const int&         mjdrefi,
                         const double&      mjdreff,
                         const std::string& timeunit,
                         const std::string& timesys,
                         const std::string& timeref)
{
    // Compute reference MJD
    double mjdref = double(mjdrefi) + mjdreff;
    
    // Set time
    set(mjdref, timeunit, timesys, timeref);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Return MJD reference (units: days)
 *
 * @return Modified Julian reference day (days).
 *
 * Returns the Modified Julian reference day.
 ***************************************************************************/
const double& GTimeReference::mjdref(void) const
{
    // Return MDJ reference
    return m_mjdref;
}


/***********************************************************************//**
 * @brief Returns integer part of MJD reference (units: days)
 *
 * @return Integer part of Modified Julian reference day (days).
 *
 * Returns the integer part of the Modified Julian reference day.
 ***************************************************************************/
int GTimeReference::mjdrefi(void) const
{
    // Return integer part of MJD reference
    return int(m_mjdref);
}


/***********************************************************************//**
 * @brief Returns fractional part of MJD reference (units: days)
 *
 * @return Fractional part of Modified Julian reference day (days).
 *
 * Returns the fractional part of the Modified Julian reference day.
 ***************************************************************************/
double GTimeReference::mjdreff(void) const
{
    // Return fractional part of MJD reference
    return (m_mjdref-double(mjdrefi()));
}


/***********************************************************************//**
 * @brief Return time unit
 *
 * @return Time unit.
 *
 * Returns the reference time unit.
 ***************************************************************************/
const std::string& GTimeReference::timeunit(void) const
{
    // Return time unit
    return m_timeunit;
}


/***********************************************************************//**
 * @brief Return time system
 *
 * @return Time system.
 *
 * Returns the reference time system.
 ***************************************************************************/
const std::string& GTimeReference::timesys(void) const
{
    // Return time system
    return m_timesys;
}


/***********************************************************************//**
 * @brief Return time reference
 *
 * @return Time reference.
 *
 * Returns the reference time reference.
 ***************************************************************************/
const std::string& GTimeReference::timeref(void) const
{
    // Return time reference
    return m_timeref;
}


/***********************************************************************//**
 * @brief Return the time unit in seconds
 *
 * @return Time unit in seconds.
 *
 * Returns 1 if the time using is in seconds and 86400 if the time unit is
 * in days.
 ***************************************************************************/
double GTimeReference::unitseconds(void) const
{
    // Set time unit in seconds
    double unit = (m_unit_sec) ? 1.0 : sec_in_day;

    // Return time unit
    return unit;
}


/***********************************************************************//**
 * @brief Print time reference
 *
 * @return String containing the time reference.
 ***************************************************************************/
std::string GTimeReference::print(void) const
{
    // Initialise result string
    std::string result;

    // Append header
    result.append("=== GTimeReference ===\n");

    // Append information
    result.append(parformat("MJD reference time")+str(mjdref())+"\n");
    result.append(parformat("Time unit")+timeunit()+"\n");
    result.append(parformat("Time system")+timesys()+"\n");
    result.append(parformat("Time reference")+timeref()+"\n");

    // Return
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
void GTimeReference::init_members(void)
{
    // Initialise members
    m_mjdref   = 0.0;
    m_timeunit = "secs";
    m_timesys  = "TT";
    m_timeref  = "local";
    m_unit_sec = true;
  
    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] ref Time reference.
 ***************************************************************************/
void GTimeReference::copy_members(const GTimeReference& ref)
{
    // Copy members
    m_mjdref   = ref.m_mjdref;
    m_timeunit = ref.m_timeunit;
    m_timesys  = ref.m_timesys;
    m_timeref  = ref.m_timeref;
    m_unit_sec = ref.m_unit_sec;
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GTimeReference::free_members(void)
{
    // Return
    return;
}
