/***************************************************************************
 *             GApplication.cpp - GammaLib application base class          *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010 by Jurgen Knodlseder                                *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
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
#define G_PAR                                "GApplication::par(std::string)"

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
    // Initialise private members for clean destruction
    init_members();

    // Save the execution start time
    time(&m_tstart);

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
    // Initialise private members for clean destruction
    init_members();

    // Set application name and version
    m_name    = name;
    m_version = version;

    // Save the execution start time
    time(&m_tstart);

    // Initialise the application logger
    log.open(log_filename(), true);

    // Initialise application parameters
    m_pars.load(par_filename());

    // Get standard parameters from parameter file
    get_par_standard();

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
    // Initialise private members for clean destruction
    init_members();

    // Set application name and version
    m_name    = name;
    m_version = version;

    // Save the execution start time
    time(&m_tstart);

    // Initialise the application logger
    log.open(log_filename(), true);

    // Save arguments as vector of strings
    for (int i = 0; i < argc; ++i)
        m_args.push_back(strip_whitespace(argv[i]));

    // Initialise application parameters
    m_pars.load(par_filename(), m_args);

    // Get standard parameters from parameter file
    get_par_standard();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] app Object from which the instance should be built.
 ***************************************************************************/
GApplication::GApplication(const GApplication& app)
{
    // Initialise private members for clean destruction
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
 * @param[in] app Object which should be assigned.
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
 * @brief Return application elapsed time in CPU seconds
 ***************************************************************************/
double GApplication::telapse(void) const
{
    // Get actual time
    time_t acttime;
    time(&acttime);

    // Compute elapsed time
    double telapse = difftime(acttime, m_tstart);

    // Return elapsed time
    return telapse;
}


/***********************************************************************//**
 * @brief Return pointer to application parameter
 *
 * @param[in] name Parameter to be returned.
 *
 * @exception GException::par_error
 *            Parameter not found..
 ***************************************************************************/
GPar* GApplication::par(const std::string& name)
{
    // Get pointer to application parameter
    GPar* ptr = m_pars.par(name);

    // Throw exception if parameter does not exist
    if (ptr == NULL) {
        throw GException::par_error(G_PAR, "Parameter '"+name+"' not found "
                                    "in parameter file.");
    }

    // Return pointer
    return ptr;
}


/***********************************************************************//**
 * @brief Signal terse logging
 *
 * The terse level is used for the most crucial application information that
 * should be logged in all cases (chatter >= 1).
 ***************************************************************************/
bool GApplication::logTerse(void) const
{
    // Return terse logging condition
    return (m_chatter > 0);
}


/***********************************************************************//**
 * @brief Signal normal logging
 *
 * The normal level is used for standard application information that should
 * be logged in normal operations (chatter >= 2).
 ***************************************************************************/
bool GApplication::logNormal(void) const
{
    // Return normal logging condition
    return (m_chatter > 1);
}


/***********************************************************************//**
 * @brief Signal explicit logging
 *
 * The explicit level is used for detailed application information that
 * should be logged in detailed studies (chatter >= 3).
 ***************************************************************************/
bool GApplication::logExplicit(void) const
{
    // Return explicit logging condition
    return (m_chatter > 2);
}


/***********************************************************************//**
 * @brief Signal verbose logging
 *
 * The verbose level is used for full application information (chatter >= 4).
 ***************************************************************************/
bool GApplication::logVerbose(void) const
{
    // Return verbose logging condition
    return (m_chatter > 3);
}


/***********************************************************************//**
 * @brief Signal debug logging
 *
 * The debug level is used for application debugging.
 ***************************************************************************/
bool GApplication::logDebug(void) const
{
    // Return debug condition
    return (m_debug);
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
    m_tstart = 0;
    m_name.clear();
    m_version.clear();
    m_args.clear();
    m_pars.clear();
    m_chatter = 1;
    m_clobber = true;
    m_debug   = false;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] app Object from which members which should be copied.
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
    m_pars    = app.m_pars;
    m_chatter = app.m_chatter;
    m_clobber = app.m_clobber;
    m_debug   = app.m_debug;

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
 * The log filename is given by the task name to which the suffix '.log' is
 * added.
 ***************************************************************************/
std::string GApplication::log_filename(void) const
{
    // Return
    return (m_name+".log");
}


/***********************************************************************//**
 * @brief Get standard parameters from parameter file
 ***************************************************************************/
void GApplication::get_par_standard(void)
{
    // Get parameter pointers (NULL if not found)
    GPar* par_chatter = par("chatter");
    GPar* par_clobber = par("clobber");
    GPar* par_debug   = par("debug");

    // Extract parameters if they exist
    if (par_chatter != NULL)
        m_chatter = par_chatter->integer();
    if (par_clobber != NULL)
        m_clobber = par_clobber->boolean();
    if (par_debug != NULL)
        m_debug = par_debug->boolean();
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Write application header in log file
 *
 * Dump the application header into the log file. The header is composed of
 * a fixed width block delimined by '*' characters that contains information
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
 * The application trailer gives the total number of elapsed CPU time.
 ***************************************************************************/
void GApplication::log_trailer(void)
{
    // Get elapsed time
    double t = telapse();

    // Reset any indentation
    log.indent(0);

    // Dump trailer
    log << "Application \"" << m_name << "\" terminated after";
    log << " consuming " << telapse() << " seconds of CPU time.";
    if (t < 0.1)
         log << " This was rather quick!";
    if (t > 86400.0)
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
 * @brief Put object in output stream
 *
 * @param[in] os Output stream into which the model will be dumped
 * @param[in] app Object to be dumped
 ***************************************************************************/
std::ostream& operator<< (std::ostream& os, const GApplication& app)
{
    // Put object in stream
    os << "=== GApplication ===" << std::endl;
    os << " Name ......................: " << app.name() << std::endl;
    os << " Version ...................: " << app.version();
    for (int i=0; i < app.m_args.size(); ++i) {
        if (i == 0)
            os << std::endl << " Command ...................: " << app.m_args[i];
        else if (i == 1)
            os << std::endl << " Arguments .................: " << app.m_args[i];
        else
            os << std::endl << "                              " << app.m_args[i];
    }

    // Return output stream
    return os;
}


