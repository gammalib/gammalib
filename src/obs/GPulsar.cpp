/***************************************************************************
 *                        GPulsar.cpp - Pulsar class                       *
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
 * @file GPulsar.cpp
 * @brief Pulsar class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GPulsar.hpp"
#include "GFits.hpp"
#include "GGti.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_AT                            "GPulsarEphemeris& GPulsar::at(int&)"
#define G_EPHEMERIS                              "GPulsar::ephemeris(GTime&)"
#define G_LOAD                      "GPulsar::load(GFilename&, std::string&)"
#define G_LOAD_FITS            "GPulsar::load_fits(GFilename&, std::string&)"
#define G_LOAD_INTEGRAL   "GPulsar::load_integral(GFitsTable*, std::string&)"
#define G_LOAD_FERMI         "GPulsar::load_fermi(GFitsTable*, std::string&)"
#define G_LOAD_PSRTIME      "GPulsar::load_psrtime(GFilename&, std::string&)"
#define G_LOAD_PARFILE                    "GPulsar::load_parfile(GFilename&)"

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
GPulsar::GPulsar(void)
{
    // Initialise class members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Filename constructor
 *
 * @param[in] filename File name of pulsar ephemerides
 * @param[in] name Pulsar name
 *
 * Constructs a pulsar from an ephemerides file. In case that several pulsars
 * are present in the specified ephemerides file the name of the pulsar that
 * should be constructed is to be specified.
 ***************************************************************************/
GPulsar::GPulsar(const GFilename& filename, const std::string& name)
{
    // Initialise class members
    init_members();

    // Load pulsar ephemerides
    load(filename, name);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] pulsar Pulsar.
 ***************************************************************************/
GPulsar::GPulsar(const GPulsar& pulsar)
{
    // Initialise class members
    init_members();

    // Copy members
    copy_members(pulsar);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GPulsar::~GPulsar(void)
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
 * @param[in] pulsar Pulsar.
 * @return Pulsar.
 ***************************************************************************/
GPulsar& GPulsar::operator=(const GPulsar& pulsar)
{
    // Execute only if object is not identical
    if (this != &pulsar) {

        // Free members
        free_members();

        // Initialise private members
        init_members();

        // Copy members
        copy_members(pulsar);

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
 * @brief Clear Pulsar
 ***************************************************************************/
void GPulsar::clear(void)
{
    // Free members
    free_members();

    // Initialise private members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone Pulsar
 *
 * @return Pointer to deep copy of Pulsar.
 ***************************************************************************/
GPulsar* GPulsar::clone(void) const
{
    return new GPulsar(*this);
}

/***********************************************************************//**
 * @brief Return reference to ephemeris
 *
 * @param[in] index Ephemeris index [0,...,size()-1].
 *
 * @exception GException::out_of_range
 *            Ephemeris index is out of range.
 *
 * Returns a reference to the ephemeris with the specified @p index.
 ***************************************************************************/
GPulsarEphemeris& GPulsar::at(const int& index)
{
    // Compile option: raise an exception if index is out of range
    #if defined(G_RANGE_CHECK)
    if (index < 0 || index >= size()) {
        throw GException::out_of_range(G_AT, "Ephemeris index", index, size());
    }
    #endif

    // Return reference
    return (m_ephemerides[index]);
}


/***********************************************************************//**
 * @brief Return reference to ephemeris (const version)
 *
 * @param[in] index Ephemeris index [0,...,size()-1].
 *
 * @exception GException::out_of_range
 *            Ephemeris index is out of range.
 *
 * Returns a const reference to the ephemeris with the specified @p index.
 ***************************************************************************/
const GPulsarEphemeris& GPulsar::at(const int& index) const
{
    // Compile option: raise an exception if index is out of range
    #if defined(G_RANGE_CHECK)
    if (index < 0 || index >= size()) {
        throw GException::out_of_range(G_AT, "Ephemeris index", index, size());
    }
    #endif

    // Return reference
    return (m_ephemerides[index]);
}


/***********************************************************************//**
 * @brief Return pulsar ephemeris
 *
 * param[in] time Time.
 * @return Pulsar ephemeris.
 *
 * @exception GException::invalid_argument
 *            No valid ephemeris found for specified @p time.
 *
 * Returns the pulsar ephemeris as function of time.
 ***************************************************************************/
const GPulsarEphemeris& GPulsar::ephemeris(const GTime& time) const
{
    // Find appropriate ephemerides
    int ifound = -1;
    for (int i = 0; i < size(); ++i) {
        if (time >= m_ephemerides[i].tstart() &&
            time <= m_ephemerides[i].tstop()) {
            ifound = i;
            break;
        }
    }

    // Throw an expection if no valid ephemerides were found
    if (ifound == -1) {
        std::string msg = "No valid ephemeris found for MJD "+
                          gammalib::str(time.mjd())+". Please specify "
                          "ephemerides that comprise the time.";
        throw GException::invalid_argument(G_EPHEMERIS, msg);
    }

    // Return ephemerides
    return (m_ephemerides[ifound]);
}


/***********************************************************************//**
 * @brief Return validity intervals of pulsar ephemerides
 *
 * @return Validity intervals of pulsar ephemerides.
 *
 * Returns the validity intervals of pulsar ephemerides as Good Time
 * Intervals.
 ***************************************************************************/
GGti GPulsar::validity(void) const
{
    // Initialise Good Time Intervals
    GGti gti;

    // Loop over ephemerides
    for (int i = 0; i < size(); ++i) {
        gti.append(m_ephemerides[i].tstart(), m_ephemerides[i].tstop());
    }

    // Return Good Time Intervals
    return gti;
}


/***********************************************************************//**
 * @brief Load Pulsar from ephemerides file
 *
 * @param[in] filename File name of pulsar ephemerides
 * @param[in] name Pulsar name
 *
 * @exception GException::file_error
 *            Could not open specified file as an ASCII file
 *
 * Load a pulsar from an ephemerides file. In case that several pulsars are
 * present in the specified ephemerides file the name of the pulsar that
 * should be constructed is to be specified.
 ***************************************************************************/
void GPulsar::load(const GFilename& filename, const std::string& name)
{
    // Clear pulsar
    clear();

    // If filename is a FITS file then load ephemerides from a FITS file
    if (filename.is_fits()) {
        load_fits(filename, name);
    }

    // ... otherwise we expect that file is an ASCII file and we dispatch
    // according to ASCII file content
    else {

        // Allocate line buffer
        const int n = 1000;
        char  line[n];

        // Open file as ASCII file
        FILE* fptr = std::fopen(filename.url().c_str(), "r");
        if (fptr == NULL) {
            std::string msg = "Pulsar ephemerides file \""+filename.url()+
                              "\" not found or readable. Please specify a "
                              "valid and readable ephemerides file.";
            throw GException::file_error(G_LOAD, msg);
        }

        // Read first line in buffer
        if (std::fgets(line, n, fptr) != NULL) {

            // Close file as we don't need it anymore
            std::fclose(fptr);

            // Split line in elements
            std::vector<std::string> elements = gammalib::split(line, " ");

            // If there are many elements in line we deal with a file in
            // psrtime format
            if (elements.size() > 10) {
                load_psrtime(filename, name);
            }

            // ... otherwise we deal with a tempo2 par file
            else {
                load_parfile(filename);
            }

        } // endif: there was a line in the file

        // ... otherwise the file was empty and we simply close it
        else {
            std::fclose(fptr);
        }
    
    } // endelse: handled ASCII files

    // Return
    return;
}


/***********************************************************************//**
 * @brief Print Pulsar
 *
 * @param[in] chatter Chattiness.
 * @return String containing Pulsar information.
 ***************************************************************************/
std::string GPulsar::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Get ephemerides validity interval
        GGti gti = validity();

        // Append header
        result.append("=== GPulsar ===");

        // Append information
        result.append("\n"+gammalib::parformat("Pulsar name"));
        result.append(m_name);
        result.append("\n"+gammalib::parformat("Number of ephemerides"));
        result.append(gammalib::str(m_ephemerides.size()));
        if (!is_empty()) {
            result.append("\n"+gammalib::parformat("Validity MJD range"));
            result.append(gammalib::str(gti.tstart().mjd()));
            result.append(" - ");
            result.append(gammalib::str(gti.tstop().mjd()));
            result.append("\n"+gammalib::parformat("Validity UTC range"));
            result.append(gti.tstart().utc());
            result.append(" - ");
            result.append(gti.tstop().utc());
        }

        // EXPLICIT: Append ephemerides
        if (chatter >= EXPLICIT) {
            for (int i = 0; i < m_ephemerides.size(); ++i) {
                result.append("\n"+m_ephemerides[i].print());
            }
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
void GPulsar::init_members(void)
{
    // Initialise members
    m_name.clear();
    m_ephemerides.clear();
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] pulsar Pulsar.
 ***************************************************************************/
void GPulsar::copy_members(const GPulsar& pulsar)
{
    // Copy members
    m_name        = pulsar.m_name;
    m_ephemerides = pulsar.m_ephemerides;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GPulsar::free_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Load Pulsar from ephemerides FITS file
 *
 * @param[in] filename Name of pulsar ephemerides FITS file
 * @param[in] name Pulsar name
 *
 * @exception GException::file_error
 *            FITS file format not recognised
 *
 * Load a pulsar from an ephemerides FITS file. In case that several pulsars
 * are present in the specified ephemerides file the name of the pulsar that
 * should be constructed is to be specified.
 ***************************************************************************/
void GPulsar::load_fits(const GFilename& filename, const std::string& name)
{
    // Open FITS file
    GFits fits(filename);

    // If FITS file is an INTEGRAL file then load ephemerides from INTEGRAL
    // FITS table
    if (fits.contains("GNRL-EPHE-CAT")) {
        const GFitsTable* table = fits.table("GNRL-EPHE-CAT");
        load_integral(table, name);
    }

    // ... otherwise if FITS file is a Fermi D4 file then load ephemerides
    // from Fermi D4 table
    else if (fits.contains("SPIN_PARAMETERS")) {
        const GFitsTable* table = fits.table("SPIN_PARAMETERS");
        load_fermi(table, name);
    }

    // ... otherwise signal that the FITS file was not recognised
    else {
        std::string msg = "Could not recognised format of FITS file \""+
                          filename.url()+"\". Please specify either an "
                          "INTEGRAL or a D4 Fermi file.";
        throw GException::file_error(G_LOAD_FITS, msg);
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Load Pulsar from INTEGRAL ephemerides FITS table
 *
 * @param[in] table INTEGRAL ephemerides FITS table
 * @param[in] name Pulsar name
 *
 * @exception GException::invalid_argument
 *            Specified pulsar name not found in FITS table
 *
 * Load a pulsar from INTEGRAL ephemerides FITS table. In the INTEGRAL format
 * only a single pulsar will be present in the FITS file, specified by
 * the SOURCEID keyword.
 ***************************************************************************/
void GPulsar::load_integral(const GFitsTable* table, const std::string& name)
{
    // Get pulsar attributes from FITS table
    m_name     = table->string("SOURCEID");
    double ra  = table->real("RA_OBJ");
    double dec = table->real("DEC_OBJ");

    // Check name if provided
    if (!name.empty()) {
        if (name != m_name) {
            std::string msg = "Pulsar \""+m_name+"\" found in ephemerides "
                              "file yet pulsar \""+name+"\" was requested. "
                              "Please check the requested pulsar name.";
            throw GException::invalid_argument(G_LOAD_INTEGRAL, msg);
        }
    }

    // Set sky direction
    GSkyDir dir;
    dir.radec_deg(ra, dec);

    // Get pointers to table columns
    const GFitsTableCol* mjdstart = (*table)["MJDSTART"];
    const GFitsTableCol* mjdstop  = (*table)["MJDSTOP"];
    const GFitsTableCol* tref0    = (*table)["TREF0"];
    const GFitsTableCol* t2peak   = (*table)["T2PEAK"];
    const GFitsTableCol* f0       = (*table)["F0"];
    const GFitsTableCol* f1       = (*table)["FDOT0"];
    const GFitsTableCol* f2       = (*table)["FDOTDOT0"];

    // Extract ephemerides
    for (int i = 0; i < mjdstart->nrows(); ++i) {

        // Set validity interval
        GTime tstart;
        GTime tstop;
        tstart.mjd(mjdstart->real(i));
        tstop.mjd(mjdstop->real(i));

        // Set t0
        GTime t0;
        t0.mjd(tref0->real(i));
        t0 += t2peak->real(i);

        // Allocate ephemeris
        GPulsarEphemeris ephemeris;

        // Set attributes
        ephemeris.name(m_name);
        ephemeris.dir(dir);
        ephemeris.tstart(tstart);
        ephemeris.tstop(tstop);
        ephemeris.t0(t0);
        ephemeris.f0(f0->real(i));
        ephemeris.f1(f1->real(i));
        ephemeris.f2(f2->real(i));

        // Push back ephemeris
        m_ephemerides.push_back(ephemeris);

    } // endfor: looped over ephemerides

    // Return
    return;
}


/***********************************************************************//**
 * @brief Load Pulsar from Fermi ephemerides FITS table
 *
 * @param[in] table Fermi ephemerides FITS table
 * @param[in] name Pulsar name
 *
 * @exception GException::invalid_argument
 *            Specified pulsar name not found in FITS table
 *
 * Load a pulsar from Fermi ephemerides FITS table.
 ***************************************************************************/
void GPulsar::load_fermi(const GFitsTable* table, const std::string& name)
{
    // Get pointers to table columns
    const GFitsTableCol* psrname      = (*table)["PSRNAME"];
    const GFitsTableCol* ra           = (*table)["RA"];
    const GFitsTableCol* dec          = (*table)["DEC"];
    const GFitsTableCol* valid_since  = (*table)["VALID_SINCE"];
    const GFitsTableCol* valid_until  = (*table)["VALID_UNTIL"];
    const GFitsTableCol* toabary_int  = (*table)["TOABARY_INT"];
    const GFitsTableCol* toabary_frac = (*table)["TOABARY_FRAC"];
    const GFitsTableCol* f0           = (*table)["F0"];
    const GFitsTableCol* f1           = (*table)["F1"];
    const GFitsTableCol* f2           = (*table)["F2"];

    // Extract ephemerides
    for (int i = 0; i < psrname->nrows(); ++i) {

        // Skip if pulsar name does not correspond to specified name
        if (psrname->string(i) != name) {
            continue;
        }

        // Save the pulsar name
        m_name = psrname->string(i);

        // Set sky direction
        GSkyDir dir;
        dir.radec_deg(ra->real(i), dec->real(i));

        // Set validity interval
        GTime tstart;
        GTime tstop;
        tstart.mjd(valid_since->real(i));
        tstop.mjd(valid_until->real(i));

        // Set t0
        GTime t0;
        t0.mjd(toabary_int->integer(i)+toabary_frac->real(i));

        // Allocate ephemeris
        GPulsarEphemeris ephemeris;

        // Set attributes
        ephemeris.name(m_name);
        ephemeris.dir(dir);
        ephemeris.tstart(tstart);
        ephemeris.tstop(tstop);
        ephemeris.t0(t0);
        ephemeris.f0(f0->real(i));
        ephemeris.f1(f1->real(i));
        ephemeris.f2(f2->real(i));

        // Push back ephemeris
        m_ephemerides.push_back(ephemeris);

    } // endfor: looped over ephemerides

    // Signal if no pulsar was found
    if (m_ephemerides.empty()) {
        if (name.empty()) {
            std::string msg = "No pulsar name was specified. Please specify "
                              "the name of the pulsar for which ephemerides "
                              "should be extracted from the file.";
            throw GException::invalid_argument(G_LOAD_FERMI, msg);
        }
        else {
            std::string msg = "No pulsar \""+name+"\" found in ephemerides "
                              "file. Please check the specified file or the "
                              "requested pulsar name.";
            throw GException::invalid_argument(G_LOAD_FERMI, msg);
        }
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Load Pulsar from ephemerides psrtime file
 *
 * @param[in] filename Name of pulsar ephemerides psrtime file
 * @param[in] name Pulsar name
 *
 * @exception GException::file_error
 *            Could not open ephemerides ASCII file
 * @exception GException::invalid_argument
 *            No pulsar name specified for a file that contains multiple
 *            pulsars
 *
 * Load a pulsar from an ephemerides psrtime file. In case that several
 * pulsars are present in the specified ephemerides file the name of the
 * pulsar that should be constructed is to be specified.
 ***************************************************************************/
void GPulsar::load_psrtime(const GFilename& filename, const std::string& name)
{
    // Allocate line buffer
    const int n = 1000;
    char  line[n];

    // Open file
    FILE* fptr = std::fopen(filename.url().c_str(), "r");
    if (fptr == NULL) {
        std::string msg = "Pulsar ephemerides file \""+filename.url()+
                          "\" not found or readable. Please specify a "
                          "valid and readable ephemerides file.";
        throw GException::file_error(G_LOAD_PSRTIME, msg);
    }

    // Initialise used pulsar name
    std::string psrname;

    // Read lines
    while (std::fgets(line, n, fptr) != NULL) {

        // Split line in elements
        std::vector<std::string> elements = gammalib::split(line, " ");

        // Get pulsar name
        std::string psrb = elements[0];

        // Skip if pulsar name does not correspond to specified name
        if (!name.empty()) {
            std::string psrb_alt = "PSR B" + psrb;
            if ((psrb != name) && (psrb_alt != name)) {
                continue;
            }
        }

        // ... otherwise if no pulsar name was specified then check
        // that the file only contains a single pulsar name
        else {
            if (!psrname.empty() && (psrb != psrname)) {
                std::string msg = "Pulsar ephemerides file \""+filename.url()+
                                  "\" contains several pulsar but not pulsar "
                                  "name was specified. Please specify for which "
                                  "pulsar you want to load ephemerides.";
                throw GException::invalid_argument(G_LOAD_PSRTIME, msg);
            }
        }
        
        // Save the pulsar name
        psrname = psrb;
        m_name  = "PSR B" + psrb;

        // Convert elements in attributes
        double ra_h     = gammalib::todouble(elements[1]);
        double ra_m     = gammalib::todouble(elements[2]);
        double ra_s     = gammalib::todouble(elements[3]);
        double dec_d    = gammalib::todouble(elements[4]);
        double dec_m    = gammalib::todouble(elements[5]);
        double dec_s    = gammalib::todouble(elements[6]);
        double mjdstart = gammalib::todouble(elements[7]);
        double mjdstop  = gammalib::todouble(elements[8]);
        double t0geo    = gammalib::todouble(elements[9]);
        double f0       = gammalib::todouble(gammalib::replace_segment(elements[10],"D","E"));
        double f1       = gammalib::todouble(gammalib::replace_segment(elements[11],"D","E"));
        double f2       = gammalib::todouble(gammalib::replace_segment(elements[12],"D","E"));

        // Compute Right Ascension and Declination
        double ra  = (ra_h + ra_m/60.0 + ra_s/3600.0) * 15.0;
        double dec = (dec_d < 0.0) ? -(-dec_d + dec_m/60.0 + dec_s/3600.0)
                                   : (dec_d + dec_m/60.0 + dec_s/3600.0);
       
        // Set sky direction
        GSkyDir dir;
        dir.radec_deg(ra, dec);

        // Set validity interval
        GTime tstart;
        GTime tstop;
        tstart.mjd(mjdstart);
        tstop.mjd(mjdstop);

        // Set t0
        GTime t0;
        t0.mjd(t0geo); // TODO: conversion to Bary ???

        // Allocate ephemeris
        GPulsarEphemeris ephemeris;

        // Set attributes
        ephemeris.name(m_name);
        ephemeris.dir(dir);
        ephemeris.tstart(tstart);
        ephemeris.tstop(tstop);
        ephemeris.t0(t0);
        ephemeris.f0(f0);
        ephemeris.f1(f1);
        ephemeris.f2(f2);

        // Push back ephemeris
        m_ephemerides.push_back(ephemeris);

    } // endwhile: looped over lines

    // Close file
    std::fclose(fptr);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Load Pulsar from ephemeris par file
 *
 * @param[in] filename Name of pulsar ephemeris par file
 * @param[in] name Pulsar name
 *
 * @exception GException::file_error
 *            Could not open ephemeris ASCII file
 * @exception GException::invalid_argument
 *            No valid Right Ascension and Declination found
 *
 * Load a pulsar from an ephemeris par file.
 ***************************************************************************/
void GPulsar::load_parfile(const GFilename& filename)
{
    // Allocate line buffer
    const int n = 1000;
    char  line[n];

    // Open file
    FILE* fptr = std::fopen(filename.url().c_str(), "r");
    if (fptr == NULL) {
        std::string msg = "Pulsar ephemeris file \""+filename.url()+
                          "\" not found or readable. Please specify a "
                          "valid and readable ephemeris file.";
        throw GException::file_error(G_LOAD_PARFILE, msg);
    }

    // Allocate ephemeris
    GPulsarEphemeris ephemeris;
    double           ra  = 9999.0;
    double           dec = 9999.0;

    // Read lines
    while (std::fgets(line, n, fptr) != NULL) {

        // Split line in elements
        std::vector<std::string> elements = gammalib::split(line, " ");

        // Extract keywords
        if (elements[0] == "PSRJ") {
            m_name = "PSR " + gammalib::rstrip_chars(elements[1], "\n");
            ephemeris.name(m_name);
        }
        else if (elements[0] == "RAJ") {
            std::vector<std::string> ras = gammalib::split(elements[1], ":");
            double ra_h                  = gammalib::todouble(ras[0]);
            double ra_m                  = gammalib::todouble(ras[1]);
            double ra_s                  = gammalib::todouble(ras[2]);
            ra = (ra_h + ra_m/60.0 + ra_s/3600.0) * 15.0;
        }
        else if (elements[0] == "DECJ") {
            std::vector<std::string> decs = gammalib::split(elements[1], ":");
            double dec_d                  = gammalib::todouble(decs[0]);
            double dec_m                  = gammalib::todouble(decs[1]);
            double dec_s                  = gammalib::todouble(decs[2]);
            dec = (dec_d < 0.0) ? -(-dec_d + dec_m/60.0 + dec_s/3600.0)
                                : (dec_d + dec_m/60.0 + dec_s/3600.0);
        }
        else if (elements[0] == "START") {
            GTime tstart;
            tstart.mjd(gammalib::todouble(elements[1]));
            ephemeris.tstart(tstart);
        }
        else if (elements[0] == "FINISH") {
            GTime tstop;
            tstop.mjd(gammalib::todouble(elements[1]));
            ephemeris.tstop(tstop);
        }
        else if (elements[0] == "TZRMJD") {  // TODO: Are you sure???
            GTime t0;
            t0.mjd(gammalib::todouble(elements[1]));
            ephemeris.t0(t0);
        }
        else if (elements[0] == "F0") {
            ephemeris.f0(gammalib::todouble(elements[1]));
        }
        else if (elements[0] == "F1") {
            ephemeris.f1(gammalib::todouble(elements[1]));
        }
        else if (elements[0] == "F2") {
            ephemeris.f2(gammalib::todouble(elements[1]));
        }

    } // endwhile: looped over lines

    // If Right Ascension and Declination is valid then set sky direction
    // and push back ephemeris
    if (ra >= 0.0 && ra <= 360.0 && dec >= -90.0 && dec <= 90.0) {
        GSkyDir dir;
        dir.radec_deg(ra, dec);
        ephemeris.dir(dir);
        m_ephemerides.push_back(ephemeris);
    }

    // ... otherwise throw an exception
    else {
        std::string msg = "No valid Right Ascension or Declination found in "
                          "pulsar ephemeris file \""+filename.url()+"\". "
                          "Please verify the ephemeris file.";
        throw GException::invalid_argument(G_LOAD_PARFILE, msg);
    }

    // Close file
    std::fclose(fptr);

    // Return
    return;
}
