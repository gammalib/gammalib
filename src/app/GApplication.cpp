/***************************************************************************
 *             GApplication.cpp - GammaLib application base class          *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2013 by Juergen Knoedlseder                         *
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
 * @file GApplication.cpp
 * @brief GammaLib application base class
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GApplication.hpp"
#include "GTools.hpp"

/* __ Method name definitions ____________________________________________ */

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */
const int header_width = 80;                                //!< Header width

/* __ Debug definitions __________________________________________________ */


/*==========================================================================
 =                                                                         =
 =                         Constructors/destructors                        =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Void constructor
 ***************************************************************************/
GApplication::GApplication(void)
{
    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Application constructor
 *
 * @param[in] name Application name.
 * @param[in] version Application version.
 *
 * Constructs an application from an application @p name and @p version. The
 * application parameters will be loaded from the parameter file. No log
 * file will be opened.
 *
 * This constructor should be used for Python scripts.
 ***************************************************************************/
GApplication::GApplication(const std::string& name, const std::string& version)
{
    // Initialise members
    init_members();

    // Set application name and version
    m_name    = name;
    m_version = version;

    // Initialise application parameters
    m_pars.load(par_filename());

    // Return
    return;
}
    

/***********************************************************************//**
 * @brief Application constructor
 *
 * @param[in] name Application name.
 * @param[in] version Application version.
 * @param[in] argc Number of command line arguments
 * @param[in] argv Command line arguments
 *
 * Constructs an application from an application @p name, @p version and a
 * number @p argc of command line arguments @p argv. The application
 * parameters will be loaded from the parameter file and the log file will
 * be opened.
 *
 * This constructor should be used for C++ applications.
 ***************************************************************************/
GApplication::GApplication(const std::string& name, const std::string& version,
                           int argc, char *argv[])
{
    // Initialise members
    init_members();

    // Set application name and version
    m_name    = name;
    m_version = version;

    // Save arguments as vector of strings
    for (int i = 0; i < argc; ++i) {
        m_args.push_back(gammalib::strip_whitespace(argv[i]));
    }

    // Initialise application parameters
    m_pars.load(par_filename(), m_args);

    // Initialise the application logger
    logFileOpen();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] app Application.
 ***************************************************************************/
GApplication::GApplication(const GApplication& app)
{
    // Initialise members
    init_members();

    // Copy members
    copy_members(app);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GApplication::~GApplication(void)
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
 * @param[in] app Application.
 * @return Application.
 ***************************************************************************/
GApplication& GApplication::operator=(const GApplication& app)
{
    // Execute only if object is not identical
    if (this != &app) {

        // Free members
        free_members();

        // Initialise private members for clean destruction
        init_members();

        // Copy members
        copy_members(app);

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
 * @brief Clear application
 ***************************************************************************/
void GApplication::clear(void)
{
    // Free members
    free_members();

    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone application
 *
 * @return Pointer to deep copy of application.
 ***************************************************************************/
GApplication* GApplication::clone(void) const
{
    // Clone application
    return new GApplication(*this);
}


/***********************************************************************//**
 * @brief Return application elapsed time in calendar seconds
 *
 * @return Application elapsed time in calendar seconds.
 ***************************************************************************/
double GApplication::telapse(void) const
{
    // Get actual time
    std::time_t acttime;
    std::time(&acttime);

    // Compute elapsed time
    double telapse = difftime(acttime, m_tstart);

    // Return elapsed time
    return telapse;
}


/***********************************************************************//**
 * @brief Return application elapsed time in CPU seconds
 *
 * @return Application elapsed time in CPU seconds.
 ***************************************************************************/
double GApplication::celapse(void) const
{
    // Compute elapsed CPU clock time
    double celapse = ((double) (clock() - m_cstart)) / CLOCKS_PER_SEC;

    // Return elapsed time
    return celapse;
}


/***********************************************************************//**
 * @brief Open log file
 *
 * @param[in] clobber (defaults to true)
 *
 * Opens application log file and sets the logger chattiness dependent on
 * the chatter parameter of the application.
 ***************************************************************************/
void GApplication::logFileOpen(const bool& clobber)
{
    // Initialise the application logger
    log.open(log_filename(), clobber);

    // Set logger chattiness
    set_log_chatter();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Signal terse logging
 *
 * @return True if at least terse logging is requested.
 * 
 * The terse level is used for the most crucial application information that
 * should be logged in all cases (chatter >= 1).
 ***************************************************************************/
bool GApplication::logTerse(void) const
{
    // Get chatter level (circumvent const correctness)
    int chatter = ((GPar*)&m_pars["chatter"])->integer();

    // Return terse logging condition
    return (chatter > 0);
}


/***********************************************************************//**
 * @brief Signal normal logging
 *
 * @return True if at least normal logging is requested.
 * 
 * The normal level is used for standard application information that should
 * be logged in normal operations (chatter >= 2).
 ***************************************************************************/
bool GApplication::logNormal(void) const
{
    // Get chatter level (circumvent const correctness)
    int chatter = ((GPar*)&m_pars["chatter"])->integer();

    // Return normal logging condition
    return (chatter > 1);
}


/***********************************************************************//**
 * @brief Signal explicit logging
 *
 * @return True if at least explicit logging is requested.
 * 
 * The explicit level is used for detailed application information that
 * should be logged in detailed studies (chatter >= 3).
 ***************************************************************************/
bool GApplication::logExplicit(void) const
{
    // Get chatter level (circumvent const correctness)
    int chatter = ((GPar*)&m_pars["chatter"])->integer();

    // Return explicit logging condition
    return (chatter > 2);
}


/***********************************************************************//**
 * @brief Signal verbose logging
 *
 * @return True if verbose logging is requested.
 * 
 * The verbose level is used for full application information (chatter >= 4).
 ***************************************************************************/
bool GApplication::logVerbose(void) const
{
    // Get chatter level (circumvent const correctness)
    int chatter = ((GPar*)&m_pars["chatter"])->integer();

    // Return verbose logging condition
    return (chatter > 3);
}


/***********************************************************************//**
 * @brief Signal debug logging
 *
 * @return True if the debugging mode has been selected.
 * 
 * The debug level is used for application debugging.
 ***************************************************************************/
bool GApplication::logDebug(void) const
{
    // Get debug condition (circumvent const correctness)
    bool debug = ((GPar*)&m_pars["debug"])->boolean();

    // Return debug condition
    return (debug);
}


/***********************************************************************//**
 * @brief Return clobber
 *
 * @return True if the clobber parameter is true.
 * 
 * The clobber indicates if existing files should be overwritten.
 ***************************************************************************/
bool GApplication::clobber(void) const
{
    // Get clobber condition (circumvent const correctness)
    bool clobber = ((GPar*)&m_pars["clobber"])->boolean();

    // Return clobber condition
    return (clobber);
}


/***********************************************************************//**
 * @brief Signal if specified parameter exists
 *
 * @param[in] name Parameter name.
 * @return True if an application parameter with the specified name exists.
 ***************************************************************************/
bool GApplication::haspar(const std::string& name) const
{
    // Return test result
    return (m_pars.haspar(name));
}


/***********************************************************************//**
 * @brief Returns parameter filename
 *
 * @return Parameter filename.
 *
 * The parameter filename is given by the task name to which the suffix
 * '.par' is added.
 ***************************************************************************/
std::string GApplication::par_filename(void) const
{
    // Return
    return (m_name+".par");
}


/***********************************************************************//**
 * @brief Returns log filename
 *
 * @return Log filename.
 *
 * The log filename is given by the task name to which the suffix ".log" is
 * added.
 ***************************************************************************/
std::string GApplication::log_filename(void) const
{
    // Return
    return (m_name+".log");
}


/***********************************************************************//**
 * @brief Write application header in log file
 *
 * Dump the application header into the log file. The header is composed of
 * a fixed width block delimined by "*" characters that contains information
 * about the application name and version.
 ***************************************************************************/
void GApplication::log_header(void)
{
    // Reset any indentation
    log.indent(0);

    // Dump header
    log << gammalib::fill("*", header_width) << std::endl;
    log << "*" << gammalib::center(m_name, header_width-2) << "*" << std::endl;
    log << "* " << gammalib::fill("-", header_width-4) << " *" << std::endl;
    log << "* Version: " << gammalib::left(m_version, header_width-12) << "*" << std::endl;
    log << gammalib::fill("*", header_width) << std::endl;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write application trailer in log file
 *
 * The application trailer gives the total number of elapsed calendar and
 * CPU seconds.
 ***************************************************************************/
void GApplication::log_trailer(void)
{
    // Reset any indentation
    log.indent(0);

    // Dump trailer
    log << "Application \"" << m_name << "\" terminated after ";
    log << telapse() << " wall clock seconds, consuming ";
    log << celapse() << " seconds of CPU time.";
    if (telapse() < 0.1) {
         log << " This was rather quick!";
    }
    if (telapse() > 86400.0) {
         log << " Hope it was worth waiting ...";
    }
    log << std::endl;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write application parameters in log file
 *
 * Writes all application parameters in the log file.
 ***************************************************************************/
void GApplication::log_parameters(void)
{
    // Write header
    log.header1("Parameters");

    // Write parameters in logger
    for (int i = 0; i < m_pars.size(); ++i) {

        // Add line feed
        if (i > 0) {
            log << std::endl;
        }

        // Set parameter name
        std::string name = " " + m_pars.m_pars[i].name() + " ";
        name = name + gammalib::fill(".", 28-name.length()) + ": ";

        // Write parameter
        log << name << m_pars.m_pars[i].m_value;

    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Print application
 *
 * @param[in] chatter Chattiness (defaults to NORMAL).
 * @return String containing application information.
 ***************************************************************************/
std::string GApplication::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GApplication ===");

        // Append application name, version and arguments
        result.append("\n"+gammalib::parformat("Name")+name());
        result.append("\n"+gammalib::parformat("Version")+version());
        for (int i = 0; i < m_args.size(); ++i) {
            if (i == 0) {
                result.append("\n"+gammalib::parformat("Command")+m_args[i]);
            }
            else if (i == 1) {
                result.append("\n"+gammalib::parformat("Arguments")+m_args[i]);
            }
            else {
                result.append("\n"+gammalib::parformat("                           ")+m_args[i]);
            }
        }

        // Append parameters
        for (int i = 0; i < m_pars.size(); ++i) {
            result.append("\n"+gammalib::parformat(m_pars.m_pars[i].name()) +
                          m_pars.m_pars[i].m_value);
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
void GApplication::init_members(void)
{
    // Initialise public members
    log.clear();

    // Initialise protected members
    m_name.clear();
    m_version.clear();
    m_args.clear();
    m_pars.clear();

    // Save the execution calendar start time
    std::time(&m_tstart);

    // Save the execution start clock
    m_cstart = std::clock();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] app Application.
 ***************************************************************************/
void GApplication::copy_members(const GApplication& app)
{
    // Copy public attributes
    log = app.log;

    // Copy protected attributes
    m_name    = app.m_name;
    m_version = app.m_version;
    m_args    = app.m_args;
    m_tstart  = app.m_tstart;
    m_cstart  = app.m_cstart;
    m_pars    = app.m_pars;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GApplication::free_members(void)
{
    // Save application parameters
    m_pars.save(par_filename());

    // Put trailer in log file
    log_trailer();

    // Close log file
    log.close();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set chattiness of logger
 ***************************************************************************/
void GApplication::set_log_chatter(void)
{
    // Get chattiness from application parameter
    int chatter = (const_cast<GPar*>(&m_pars["chatter"]))->integer();

    // Make sure that chatter is within valid range
    if (chatter < 0) {
        chatter = 0;
    }
    else if (chatter > 4) {
        chatter = 4;
    }

    // Set logger chatter
    log.chatter((GChatter)chatter);

    // Return
    return;
}
