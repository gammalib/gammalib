/***************************************************************************
 *                 GCaldb.cpp - Calibration database class                 *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2011-2014 by Juergen Knoedlseder                         *
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
 * @file GCaldb.cpp
 * @brief Calibration database class implementation.
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <unistd.h>        // access() function
#include <cstdlib>         // std::getenv() function
#include "GException.hpp"
#include "GTools.hpp"
#include "GCaldb.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_CALDB1                                                   "GCaldb()"
#define G_CALDB2                                        "GCaldb(std::string)"
#define G_ROOTDIR                                         "GCaldb::rootdir()"
#define G_SET_DATABASE                    "GCaldb::set_database(std::string)"
#define G_PATH                       "GCaldb::path(std::string, std::string)"
#define G_CIF                         "GCaldb::cif(std::string, std::string)"

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
 * The void constructor determines the calibration database root directory
 * from the GAMMALIB_CALDB or CALDB environment variables.
 ***************************************************************************/
GCaldb::GCaldb(void)
{
    // Initialise class members
    init_members();
    
    // Set calibration database root directory
    set_database(rootdir());

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] caldb Calibration database.
 ***************************************************************************/
GCaldb::GCaldb(const GCaldb& caldb)
{
    // Initialise class members for clean destruction
    init_members();

    // Copy members
    copy_members(caldb);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Initalisation constructor
 *
 * @param[in] pathname Calibration database root directory.
 *
 * This constructor sets the calibration database using the specified root
 * directory. Unless the specified pathname is not empty, any existing CALDB
 * or GAMMALIB_CALDB environment variables will be ignored.
 *
 * If the pathname is empty, however, the constructor attempts to detemine
 * the CALDB root directory from the GAMMALIB_CALDB environment variable. If
 * GAMMALIB_CALDB is not found, the method will check for the CALDB
 * environment variable. If this fails, the constructor will throw an
 * exception.
 *
 * Any environment variable contained in the pathname will be expanded.
 ***************************************************************************/
GCaldb::GCaldb(const std::string& pathname)
{
    // Initialise members
    init_members();

    // If pathname is not empty then use it as root pathname
    if (pathname.length() > 0) {
        set_database(gammalib::expand_env(pathname));
    }

    // ... otherwise try to determine root pathname from GAMMALIB_CALDB
    // or CALDB environment variables.
    else {
        set_database(rootdir());
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Calibration database constructor
 *
 * @param[in] mission Mission name (case insensitive).
 * @param[in] instrument Instrument name (case insensitive).
 *
 * Constructs a calibration database instance by opening the calibration
 * database for the specified @p mission and @p instrument. Opening consists
 * of loading the Calibration Index File (CIF) in memory. Once opened,
 * calibration information can be accessed.
 *
 * For more information about opening the database, refer to GCaldb::open.
 ***************************************************************************/
GCaldb::GCaldb(const std::string& mission, const std::string& instrument)
{
    // Initialise members
    init_members();

    // Set calibration database root directory
    set_database(rootdir());

    // Open database
    open(mission, instrument);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GCaldb::~GCaldb(void)
{
    // Free members
    free_members();

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                              Operators                                  =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Assignment operator
 *
 * @param[in] caldb Calibration database.
 ***************************************************************************/
GCaldb& GCaldb::operator= (const GCaldb& caldb)
{
    // Execute only if object is not identical
    if (this != &caldb) {

        // Free members
        free_members();

        // Initialise private members
        init_members();

        // Copy members
        copy_members(caldb);

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
 * @brief Clear calibration database
 ***************************************************************************/
void GCaldb::clear(void)
{
    // Free members
    free_members();

    // Initialise private members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone calibration database
 *
 * Cloning provides a copy of the calibration database.
 ***************************************************************************/
GCaldb* GCaldb::clone(void) const
{
    // Clone this object
    return (new GCaldb(*this));
}


/***********************************************************************//**
 * @brief Returns number of entries in calibration database
 *
 * This method returns the number of entries that were found in an opened
 * calibration database. An entry corresponds to one line in the Calibration
 * Index File (CIF). If the index file is empty, or if no calibration
 * database has been opened, the method returns 0.
 ***************************************************************************/
int GCaldb::size(void) const
{
    // Initialise size
    int size = 0;

    // If CIF is open then determine number of entries
    if (m_cif != NULL) {
        size = m_cif->nrows();
    }

    // Return size
    return size;
}


/***********************************************************************//**
 * @brief Set calibration database root directory
 *
 * @param[in] pathname Calibration database root directory.
 *
 * Set a new root directory for the calibration database. Any previously
 * opened database will be closed.
 * Any environment variable present in the pathname will be expanded.
 ***************************************************************************/
void GCaldb::dir(const std::string& pathname)
{
    // Clear any existing database
    clear();

    // Set calibration database
    set_database(gammalib::expand_env(pathname));

    // Return
    return;
}


/***********************************************************************//**
 * @brief Open calibration database
 *
 * @param[in] mission Mission name (case insensitive).
 * @param[in] instrument Instrument name (case insensitive; optional).
 *
 * Opens the calibration database for a given mission and instrument. Opening
 * consists of loading the Calibration Index File (CIF) in memory. Once
 * opened, calibration information can be accessed.
 ***************************************************************************/
void GCaldb::open(const std::string& mission, const std::string& instrument)
{
    // Close any open database
    close();

    // Set full CIF filename (this also validates the mission and
    // instrument names)
    m_cifname = cifname(mission, instrument);

    // Store mission and instrument
    m_mission    = mission;
    m_instrument = instrument;

    // Open CIF FITS file
    m_fits.open(m_cifname);

    // Store pointer to first extension which holds the CIF table
    m_cif = m_fits.table(1);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Close calibration database
 *
 * Close any previously opened calibration database.
 ***************************************************************************/
void GCaldb::close(void)
{
    // Reset database parameters
    m_mission.clear();
    m_instrument.clear();
    m_cifname.clear();
    
    // Close CIF FITS file
    m_fits.close();

    // Reset CIF table pointer
    m_cif = NULL;
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Return calibration file name based on selection parameters
 *
 * @param[in] detector Detector name (not used if empty).
 * @param[in] filter Filter name (not used if empty).
 * @param[in] codename Code name (not used if empty).
 * @param[in] date Date in yyyy-mm-dd format (not used if empty).
 * @param[in] time Time in hh:mm:ss format (not used if empty).
 * @param[in] expr Boolean selection expression (not used if empty).
 * @return Calibration filename (empty if not found).
 *
 * Returns a calibration file name based on selection parameters. If more
 * files satisfy the specified selection parameters, the first file will
 * be returned.
 *
 * @todo data should support "now" and probably be implemented using < condition.
 * @todo time should support "now" and probably be implemented using < condition.
 * @todo expr should support arbitrary Boolean expressions.
 ***************************************************************************/
std::string GCaldb::filename(const std::string& detector,
                             const std::string& filter,
                             const std::string& codename,
                             const std::string& date,
                             const std::string& time,
                             const std::string& expr)
{
    // Initialise empty filename
    std::string filename;

    // Continue only if CIF is opened
    if (m_cif != NULL) {

        // Loop over all CIF entries
        for (int i = 0; i < size(); ++i) {

            // Initialise selection flags
            bool match_detector = false;
            bool match_filter   = false;
            bool match_codename = false;
            bool match_date     = false;
            bool match_time     = false;
            bool match_expr     = false;

            // Check for detector
            if (detector.length() > 0) {
                const GFitsTableCol* column = (*m_cif)["DETNAM"];
                if (detector == column->string(i)) {
                    match_detector = true;
                }
            }
            else {
                match_detector = true;
            }

            // Check for filter
            if (filter.length() > 0) {
                const GFitsTableCol* column = (*m_cif)["FILTER"];
                if (filter == column->string(i)) {
                    match_filter = true;
                }
            }
            else {
                match_filter = true;
            }

            // Check for code name
            if (codename.length() > 0) {
                const GFitsTableCol* column = (*m_cif)["CAL_CNAM"];
                if (codename == column->string(i)) {
                    match_codename = true;
                }
            }
            else {
                match_codename = true;
            }

            // Check for date
            if (date.length() > 0) {
                const GFitsTableCol* column = (*m_cif)["CAL_VSD"];
                if (date == column->string(i)) {
                    match_date = true;
                }
            }
            else {
                match_date = true;
            }

            // Check for time
            if (time.length() > 0) {
                const GFitsTableCol* column = (*m_cif)["CAL_VST"];
                if (time == column->string(i)) {
                    match_time = true;
                }
            }
            else {
                match_time = true;
            }

            // Check for expression
            if (expr.length() > 0) {
                const GFitsTableCol* column = (*m_cif)["CAL_CBD"];
                int num_columns = column->elements(i);
                for (int k = 0; k < num_columns; ++k) {
                    if (expr == column->string(i,k)) {
                        match_expr = true;
                        break;
                    }
                }
            }
            else {
                match_expr = true;
            }

            // Select file if all selection match
            if (match_detector && match_filter && match_codename &&
                match_date && match_time && match_expr) {
                const GFitsTableCol* cal_dir  = (*m_cif)["CAL_DIR"];
                const GFitsTableCol* cal_file = (*m_cif)["CAL_FILE"];
                filename = m_caldb + "/" +
                           cal_dir->string(i) + "/" +
                           cal_file->string(i);
                break;
            }

        } // endfor: looped over all CIF entries
        
    } // endif: CIF has been opened

    // Return filename
    return filename;
}


/***********************************************************************//**
 * @brief Print calibration database information
 *
 * @param[in] chatter Chattiness (defaults to NORMAL).
 * @return String containing calibration database information.
 ***************************************************************************/
std::string GCaldb::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append Header
        result.append("=== GCaldb ===");

        // Append information
        result.append("\n"+gammalib::parformat("Database root")+m_caldb);

        // Append information about opened database
        if (m_cif != NULL) {
            result.append("\n"+gammalib::parformat("Selected Mission"));
            result.append(gammalib::toupper(m_mission));
            result.append("\n"+gammalib::parformat("Selected Instrument"));
            result.append(gammalib::toupper(m_instrument));
            result.append("\n"+gammalib::parformat("Calibration Index File"));
            result.append(m_cifname);
            result.append("\n"+gammalib::parformat("Number of entries"));
            result.append(gammalib::str(size()));
        }

    } // endif: chatter was not silent
    
    // Return
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
void GCaldb::init_members(void)
{
    // Initialise members
    m_caldb.clear();
    m_mission.clear();
    m_instrument.clear();
    m_cifname.clear();
    m_fits.clear();
    m_cif = NULL;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] caldb Calibration database.
 ***************************************************************************/
void GCaldb::copy_members(const GCaldb& caldb)
{
    // Copy attributes
    m_caldb      = caldb.m_caldb;
    m_mission    = caldb.m_mission;
    m_instrument = caldb.m_instrument;
    m_cifname    = caldb.m_cifname;
    m_fits       = caldb.m_fits;

    // Set CIF table pointer. We do not copy over the pointer as the FITS
    // file has been copied here.
    if (caldb.m_cif != NULL) {
        m_cif = m_fits.table(1);
    }
    else {
        m_cif = NULL;
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GCaldb::free_members(void)
{
    // Close any open database
    close();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Return path to CALDB root directory
 *
 * @exception GException::directory_not_found
 *            Calibration directory not found.
 * @exception GException::directory_not_accessible
 *            No read permission granted to calibration directory.
 *
 * The calibration directory path is given by one of the following
 *
 *     $CALDB/data/<mission>
 *     $CALDB/data/<mission>/<instrument>
 *     $GAMMALIB_CALDB/data/<mission>
 *     $GAMMALIB_CALDB/data/<mission>/<instrument>
 *
 * where \<mission\> is the name of the mission and \<instrument\> is the
 * optional instrument name (all lower case). The arguments provided to the
 * method are transformed to lower case.
 ***************************************************************************/
std::string GCaldb::rootdir(void) const
{
    // Get root directory from GAMMALIB_CALDB and CALDB environment variables
    char* ptr1 = std::getenv("GAMMALIB_CALDB");
    char* ptr2 = std::getenv("CALDB");

    // Throw an exception if non of the environment variables was set
    if (ptr1 == NULL && ptr2 == NULL) {
        throw GException::env_not_found(G_ROOTDIR,
              "GAMMALIB_CALDB or CALDB",
              "Please set the GAMMALIB_CALDB or the CALDB environment"
              " variable to a valid calibration database root directory.");
    }

    // Initialise root directory string
    std::string rootdir;

    // If GAMMALIB_CALDB was set then check whether the directory exists
    if (ptr1 != NULL) {
        if (access(ptr1, F_OK) == 0) {
            rootdir = std::string(ptr1);
        }
    }

    // ... otherwise, if CALDB was set then check whether the directory
    // exists
    else if (ptr2 != NULL) {
        if (access(ptr2, F_OK) == 0) {
            rootdir = std::string(ptr2);
        }
    }
    
    // Return root directory
    return rootdir;
}


/***********************************************************************//**
 * @brief Set calibration database
 *
 * @param[in] pathname Calibration database root directory.
 *
 * @exception GException::directory_not_found
 *            Calibration database root directory not found.
 * @exception GException::directory_not_accessible
 *            No read permission granted for calibration database root
 *            directory.
 *
 * Set calibration database to a specific root directory. The method verifies
 * the existence of the root directory prior to setting.
 ***************************************************************************/
void GCaldb::set_database(const std::string& pathname)
{
    // Verify that specified root directory exists
    if (access(pathname.c_str(), F_OK) != 0) {
        throw GException::directory_not_found(G_SET_DATABASE, pathname);
    }

    // Verify that specified root directory allows read access
    if (access(pathname.c_str(), R_OK) != 0) {
        throw GException::directory_not_accessible(G_SET_DATABASE, pathname,
              "Requested read permission not granted.");
    }

    // Set calibration database root directory
    m_caldb = pathname;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Return path to calibration directory
 *
 * @param[in] mission Mission name (case insensitive).
 * @param[in] instrument Instrument name (case insensitive; optional).
 *
 * @exception GException::directory_not_found
 *            Calibration directory not found.
 * @exception GException::directory_not_accessible
 *            No read permission granted to calibration directory.
 *
 * The calibration directory path is given by one of the following
 *
 *     $CALDB/data/<mission>
 *     $CALDB/data/<mission>/<instrument>
 *     $GAMMALIB_CALDB/data/<mission>
 *     $GAMMALIB_CALDB/data/<mission>/<instrument>
 *
 * where \<mission\> is the name of the mission and \<instrument\> is the
 * optional instrument name (all lower case). The arguments provided to the
 * method are transformed to lower case.
 ***************************************************************************/
std::string GCaldb::path(const std::string& mission, const std::string& instrument)
{
    // Verify that mission name is valid and directory is accessible
    std::string path = m_caldb + "/data/" + gammalib::tolower(mission);
    if (access(path.c_str(), F_OK) != 0) {
        throw GException::directory_not_found(G_PATH, path,
              "Requested mission \""+gammalib::toupper(mission)+"\" not found in"
              " calibration database.");
    }
    if (access(path.c_str(), R_OK) != 0) {
        throw GException::directory_not_accessible(G_PATH, path,
              "Requested read permission not granted for mission \""+
              gammalib::toupper(mission)+"\".");
    }

    // If an instrument has been specified, verify that instrument name is
    // valid and directory is accessible
    if (instrument.length() > 0) {

        // Add instrument to path
        path += "/" + gammalib::tolower(instrument);

        // Verify path
        if (access(path.c_str(), F_OK) != 0) {
            throw GException::directory_not_found(G_PATH, path,
                  "Requested instrument \""+gammalib::toupper(instrument)+"\" on"
                  " mission \""+gammalib::toupper(mission)+"\" not found in"
                  " calibration database.");
        }
        if (access(path.c_str(), R_OK) != 0) {
            throw GException::directory_not_accessible(G_PATH, path,
                "Requested read permission not granted for instrument \""+
                gammalib::toupper(instrument)+"\" on mission \""+gammalib::toupper(mission)+
                "\".");
        }
        
    } // endif: instrument has been specified
    
    // Return path
    return path;
}


/***********************************************************************//**
 * @brief Return absolute CIF filename
 *
 * @param[in] mission Mission name (case insensitive).
 * @param[in] instrument Instrument name (case insensitive; optional).
 *
 * @exception GException::file_not_found
 *            CIF not found.
 *
 * The calibration directory path is given by one of the following
 *
 *     $CALDB/data/<mission>/caldb.indx
 *     $CALDB/data/<mission>/<instrument>/caldb.indx
 *     $GAMMALIB_CALDB/data/<mission>/caldb.indx
 *     $GAMMALIB_CALDB/data/<mission>/<instrument>/caldb.indx
 *
 * where \<mission\> is the name of the mission and \<instrument\> is the
 * optional instrument name (all lower case). The arguments provided to the
 * method are transformed to lower case.
 ***************************************************************************/
std::string GCaldb::cifname(const std::string& mission, const std::string& instrument)
{
    // Get path to calibration directory
    std::string cif = path(mission, instrument);

    // Add calibration index filename
    cif += "/caldb.indx";

    // Verify that CIF exists
    if (access(cif.c_str(), F_OK) != 0) {
        throw GException::file_not_found(G_CIF, cif,
              "Calibration Index File (CIF) not found.");
    }

    // Return cif
    return cif;
}
