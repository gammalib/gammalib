/***************************************************************************
 *             GApplication.cpp - GammaLib application base class          *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2011 by Jurgen Knodlseder                           *
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
 * @author Jurgen Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GApplication.hpp"
#include "GException.hpp"
#include "GTools.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_OPERATOR                   "GApplication::operator[](std::string&)"

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
 ***************************************************************************/
GApplication::GApplication(const std::string& name, const std::string& version)
{
    // Initialise members
    init_members();

    // Set application name and version
    m_name    = name;
    m_version = version;

    // Initialise the application logger
    //logFileOpen();

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
 ***************************************************************************/
GApplication::GApplication(const std::string& name, const std::string& version,
                           int argc, char *argv[])
{
    // Initialise members
    init_members();

    // Set application name and version
    m_name    = name;
    m_version = version;

    // Initialise the application logger
    logFileOpen();

    // Save arguments as vector of strings
    for (int i = 0; i < argc; ++i)
        m_args.push_back(strip_whitespace(argv[i]));

    // Initialise application parameters
    m_pars.load(par_filename(), m_args);

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
 ***************************************************************************/
GApplication& GApplication::operator= (const GApplication& app)
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


/***********************************************************************//**
 * @brief Parameter access operator
 *
 * @param[in] name Parameter name.
 ***************************************************************************/
GPar& GApplication::operator[](const std::string& name)
{
    // Return parameter
    return (m_pars[name]);
}


/***********************************************************************//**
 * @brief Parameter access operator (const version)
 *
 * @param[in] name Parameter name.
 ***************************************************************************/
const GPar& GApplication::operator[](const std::string& name) const
{
    // Return parameter
    return (m_pars[name]);
}


/*==========================================================================
 =                                                                         =
 =                             Public methods                              =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Return application name
 ***************************************************************************/
std::string GApplication::name(void) const
{
    // Return name
    return m_name;
}


/***********************************************************************//**
 * @brief Return application version
 ***************************************************************************/
std::string GApplication::version(void) const
{
    // Return version
    return m_version;
}


/***********************************************************************//**
 * @brief Return application elapsed time in calendar seconds
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
 * @param[in] clobber (default: true)
 *
 * Opens application log file.
 ***************************************************************************/
void GApplication::logFileOpen(bool clobber)
{
    // Initialise the application logger
    log.open(log_filename(), clobber);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Signal terse logging
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
 ***************************************************************************/
bool GApplication::haspar(const std::string& name) const
{
    // Return test result
    return (m_pars.haspar(name));
}


/***********************************************************************//**
 * @brief Print application
 ***************************************************************************/
std::string GApplication::print(void) const
{
    // Initialise result string
    std::string result;

    // Append header
    result.append("=== GApplication ===");

    // Append application name, version and arguments
    result.append("\n"+parformat("Name")+name());
    result.append("\n"+parformat("Version")+version());
    for (int i = 0; i < m_args.size(); ++i) {
        if (i == 0)
            result.append("\n"+parformat("Command")+m_args[i]);
        else if (i == 1)
            result.append("\n"+parformat("Arguments")+m_args[i]);
        else
            result.append("\n"+parformat("                           ")+m_args[i]);
    }

    // Append parameters
    for (int i = 0; i < m_pars.size(); ++i) {
        result.append("\n"+parformat(m_pars.m_pars[i].name()) +
                           m_pars.m_pars[i].m_value);
    }

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
 * @brief Returns parameter filename
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
    log << fill("*", header_width) << std::endl;
    log << "*" << center(m_name, header_width-2) << "*" << std::endl;
    log << "* " << fill("-", header_width-4) << " *" << std::endl;
    log << "* Version: " << left(m_version, header_width-12) << "*" << std::endl;
    log << fill("*", header_width) << std::endl;

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
    if (telapse() < 0.1)
         log << " This was rather quick!";
    if (telapse() > 86400.0)
         log << " We hope it was worth waiting ...";
    log << std::endl;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write application parameters in log file
 *
 * The application trailer gives the total number of elapsed CPU time.
 ***************************************************************************/
void GApplication::log_parameters(void)
{
    // Write header
    log.header1("Parameters");

    // Write parameters in logger
    for (int i = 0; i < m_pars.size(); ++i) {

        // Add line feed
        if (i > 0)
            log << std::endl;

        // Set parameter name
        std::string name = " " + m_pars.m_pars[i].name() + " ";
        name = name + fill(".", 28-name.length()) + ": ";

        // Write parameter
        log << name << m_pars.m_pars[i].m_value;

    }

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                                 Friends                                 =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Output operator
 *
 * @param[in] os Output stream.
 * @param[in] app Application.
 ***************************************************************************/
std::ostream& operator<< (std::ostream& os, const GApplication& app)
{
     // Write application in output stream
    os << app.print();

    // Return output stream
    return os;
}


/***********************************************************************//**
 * @brief Log operator
 *
 * @param[in] log Logger.
 * @param[in] app Application.
 ***************************************************************************/
GLog& operator<< (GLog& log, const GApplication& app)
{
    // Write application into logger
    log << app.print();

    // Return logger
    return log;
}
