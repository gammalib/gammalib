/***************************************************************************
 *             GApplication.cpp - GammaLib application base class          *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2018 by Juergen Knoedlseder                         *
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
 *
 * Constructs an empty application.
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
 * constructor will set the parameter filename to "<name>.par" and the log
 * filename to "<name>".log. The parameters will be loaded from the parameter
 * file.
 *
 * No log file will be opened. To open the log file an explicit call to the
 * logFileOpen() method is required.
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

    // Set default parfile and logfile name
    m_parfile = name+".par";
    m_logfile = name+".log";

    // Initialise application parameters
    m_pars.load(par_filename());

    // Signal that application parameters have been loaded
    m_pars_loaded = true;

    // Set log filename and chattiness
    set_log_filename();
    set_log_chatter();

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
 * number @p argc of command line arguments @p argv. The constructor will
 * set the parameter filename to "<name>.par" and the log filename to
 * "<name>".log. The parameters will be loaded from the parameter file.
 *
 * No log file will be opened. To open the log file an explicit call to the
 * logFileOpen() method is required.
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

    // Set default parfile and logfile name
    m_parfile = name+".par";
    m_logfile = name+".log";

    // Save arguments as vector of strings
    for (int i = 0; i < argc; ++i) {
        m_args.push_back(gammalib::strip_whitespace(argv[i]));
    }

    // Catch --help option before doing anything else. The actual action
    // needs to be done by the client, but we skip the loading of the
    // parameter file and the opening of the logger if the help option
    // was specified.
    m_need_help = false;
    if (m_args.size() >= 2) {
        if (m_args[1] == "--help") {
            m_need_help = true;
        }
    }

    // If no --help option has been specified then proceed with loading
    // the parameter file and opening the logger
    if (!m_need_help) {

        // Initialise application parameters
        m_pars.load(par_filename(), m_args);

        // Signal that application parameters have been loaded
        m_pars_loaded = true;

        // Set log filename and chattiness
        set_log_filename();
        set_log_chatter();

    } // endif: no --help option specified

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
    // Save application name and version as we need them to reconstruct
    // the application after freeing and initialising its members
    std::string name    = m_name;
    std::string version = m_version;

    // Free members
    free_members();

    // Initialise members
    init_members();

    // Recover saved application name and version
    m_name    = name;
    m_version = version;

    // Set default parfile and logfile name
    m_parfile = name+".par";
    m_logfile = name+".log";

    // Initialise application parameters
    m_pars.load(par_filename());

    // Signal that application parameters have been loaded
    m_pars_loaded = true;

    // Set log filename and chattiness
    set_log_filename();
    set_log_chatter();

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
 * @param[in] clobber Overwrite (true) or append (false) to existing file
 *                    (defaults to true)
 *
 * Opens application log file and sets the logger chattiness dependent on
 * the chatter parameter of the application.
 ***************************************************************************/
void GApplication::logFileOpen(const bool& clobber)
{
    // Open logger only if filename is non-empty
    if (!log_filename().empty()) {

        // Initialise the application logger
        log.open(log_filename(), clobber);

        // Write header into log file
        log_header();

    } // endif: file name was not empty

    // Set logger chattiness
    set_log_chatter();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Close log file
 *
 * Close the log file. In case that some characters have been written
 * through the logger a trailer will be appended to the logger before
 * closing the log file. The trailer informs about the computation time
 * used by the application.
 ***************************************************************************/
void GApplication::logFileClose(void)
{
    // Continue only if log file is open and if something has been written
    // through this logger. This avoid writing trailers for logger than
    // have not been used.
    if (log.is_open() && (log.written_size() > 0)) {

        // Write line feed before trailer into logger
        log << std::endl;
            
        // Write trailer into logger
        log_trailer();

        // Close log file
        log.close();

    } // endif: log file was open

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
    int chatter = const_cast<GApplicationPar*>(&m_pars["chatter"])->integer();

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
    int chatter = const_cast<GApplicationPar*>(&m_pars["chatter"])->integer();

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
    int chatter = const_cast<GApplicationPar*>(&m_pars["chatter"])->integer();

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
    int chatter = const_cast<GApplicationPar*>(&m_pars["chatter"])->integer();

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
    bool debug = const_cast<GApplicationPar*>(&m_pars["debug"])->boolean();

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
    bool clobber = const_cast<GApplicationPar*>(&m_pars["clobber"])->boolean();

    // Return clobber condition
    return (clobber);
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
    log << "*" << gammalib::centre(m_name, header_width-2) << "*" << std::endl;
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
 * @brief Write string in log file
 *
 * @param[in] chatter Minimum required chattiness
 * @param[in] string String
 * @param[in] linefeed Terminate string with linefeed?
 *
 * Writes a string into the log file if chattiness is at least @p chatter.
 * If @p linefeed is true the string is terminated with a linefeed.
 ***************************************************************************/
void GApplication::log_string(const GChatter&    chatter,
                              const std::string& string,
                              const bool&        linefeed)
{
    // Get chattiness of application
    GChatter chattiness = static_cast<GChatter>((&m_pars["chatter"])->integer());

    // Only write message if chattiness is at least equal to the minimum
    // required chattiness
    if (chattiness >= chatter) {
        log << string;
        if (linefeed) {
            log << std::endl;
        }
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write parameter value in log file
 *
 * @param[in] chatter Minimum required chattiness
 * @param[in] name Parameter name string
 * @param[in] value Value string
 *
 * Writes a parameter value into the log file if chattiness is at least
 * @p chatter.
 ***************************************************************************/
void GApplication::log_value(const GChatter&    chatter,
                             const std::string& name,
                             const std::string& value)
{
    // Get chattiness of application
    GChatter chattiness = static_cast<GChatter>((&m_pars["chatter"])->integer());

    // Only write message if chattiness is at least equal to the minimum
    // required chattiness
    if (chattiness >= chatter) {
        log << gammalib::parformat(name);
        log << value << std::endl;
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write parameter value in log file
 *
 * @param[in] chatter Minimum required chattiness
 * @param[in] name Parameter name string
 * @param[in] value Integer value
 *
 * Writes a parameter value into the log file if chattiness is at least
 * @p chatter.
 ***************************************************************************/
void GApplication::log_value(const GChatter&    chatter,
                             const std::string& name,
                             const int&         value)
{
    // Get chattiness of application
    GChatter chattiness = static_cast<GChatter>((&m_pars["chatter"])->integer());

    // Only write message if chattiness is at least equal to the minimum
    // required chattiness
    if (chattiness >= chatter) {
        log << gammalib::parformat(name);
        log << value << std::endl;
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write parameter value in log file
 *
 * @param[in] chatter Minimum required chattiness
 * @param[in] name Parameter name string
 * @param[in] value Floating point value
 *
 * Writes a parameter value into the log file if chattiness is at least
 * @p chatter.
 ***************************************************************************/
void GApplication::log_value(const GChatter&    chatter,
                             const std::string& name,
                             const double&      value)
{
    // Get chattiness of application
    GChatter chattiness = static_cast<GChatter>((&m_pars["chatter"])->integer());

    // Only write message if chattiness is at least equal to the minimum
    // required chattiness
    if (chattiness >= chatter) {
        log << gammalib::parformat(name);
        log << value << std::endl;
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write header 1 in log file
 *
 * @param[in] chatter Minimum required chattiness
 * @param[in] header Header string
 *
 * Writes a header of level 1 into the log file if chattiness is at least
 * @p chatter.
 ***************************************************************************/
void GApplication::log_header1(const GChatter&    chatter,
                               const std::string& header)
{
    // Get chattiness of application
    GChatter chattiness = static_cast<GChatter>((&m_pars["chatter"])->integer());

    // Only write message if chattiness is at least equal to the minimum
    // required chattiness
    if (chattiness >= chatter) {
        log << std::endl;
        log.header1(header);
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write header 2 in log file
 *
 * @param[in] chatter Minimum required chattiness
 * @param[in] header Header string
 *
 * Writes a header of level 2 into the log file if chattiness is at least
 * @p chatter.
 ***************************************************************************/
void GApplication::log_header2(const GChatter&    chatter,
                               const std::string& header)
{
    // Get chattiness of application
    GChatter chattiness = static_cast<GChatter>((&m_pars["chatter"])->integer());

    // Only write message if chattiness is at least equal to the minimum
    // required chattiness
    if (chattiness >= chatter) {
        log.header2(header);
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write header 3 in log file
 *
 * @param[in] chatter Minimum required chattiness
 * @param[in] header Header string
 *
 * Writes a header of level 3 into the log file if chattiness is at least
 * @p chatter.
 ***************************************************************************/
void GApplication::log_header3(const GChatter&    chatter,
                               const std::string& header)
{
    // Get chattiness of application
    GChatter chattiness = static_cast<GChatter>((&m_pars["chatter"])->integer());

    // Only write message if chattiness is at least equal to the minimum
    // required chattiness
    if (chattiness >= chatter) {
        log.header3(header);
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write application parameters in log file
 *
 * @param[in] chatter Minimum required chattiness
 *
 * Writes all application parameters in the log file. For parameters that
 * have not yet been queried the method does not write the current value
 * but signals [not queried].
 ***************************************************************************/
void GApplication::log_parameters(const GChatter& chatter)
{
    // Get chattiness of application
    GChatter chattiness = static_cast<GChatter>((&m_pars["chatter"])->integer());

    // Only write parameters if the chattiness is at least equal to the
    // minimum required chattiness
    if (chattiness >= chatter) {

        // Write header
        log.header1("Parameters");

        // Write parameters in logger
        for (int i = 0; i < m_pars.size(); ++i) {

            // Set parameter name
            std::string name = " " + m_pars.m_pars[i].name() + " ";
            name = name + gammalib::fill(".", 28-name.length()) + ": ";

            // Write parameter
            if (m_pars.m_pars[i].is_query() &&
                !m_pars.m_pars[i].was_queried()) {
                log << name << "[not queried]" << std::endl;
            }
            else {
                log << name << m_pars.m_pars[i].m_value << std::endl;
            }

        } // endfor: looped over all parameters

    } // endif: chattiness satisfied minimum required level

    // Return
    return;
}


/***********************************************************************//**
 * @brief Print application
 *
 * @param[in] chatter Chattiness.
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
    m_parfile.clear();
    m_logfile.clear();
    m_args.clear();
    m_pars.clear();
    m_pars_loaded = false;
    m_need_help   = false;

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
 *
 * Copies all class members.
 ***************************************************************************/
void GApplication::copy_members(const GApplication& app)
{
    // Copy public attributes
    log = app.log;

    // Copy protected attributes
    m_name        = app.m_name;
    m_version     = app.m_version;
    m_parfile     = app.m_parfile;
    m_logfile     = app.m_logfile;
    m_args        = app.m_args;
    m_tstart      = app.m_tstart;
    m_cstart      = app.m_cstart;
    m_pars        = app.m_pars;
    m_pars_loaded = app.m_pars_loaded;
    m_need_help   = app.m_need_help;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GApplication::free_members(void)
{
    // Save application parameters if they have been loaded
    if (m_pars_loaded) {
        m_pars.save(par_filename());
    }

    // Close log file
    logFileClose();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set chattiness of logger
 ***************************************************************************/
void GApplication::set_log_chatter(void)
{
    // Get chattiness from application parameter
    int chatter = m_pars["chatter"].integer();

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


/***********************************************************************//**
 * @brief Set log filename from "logfile" parameter (if it exists)
 ***************************************************************************/
void GApplication::set_log_filename(void)
{
    // Continue only if logfile parameter exists
    if (m_pars.contains("logfile")) {

        // Get log filename from application parameter
        m_logfile = m_pars["logfile"].filename();

    }

    // Return
    return;
}
