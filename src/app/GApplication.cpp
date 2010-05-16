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

    // Initialise application parameters
    m_pars.load(parfilename());
    
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
    
    // Save arguments as vector of strings
    for (int i=0; i < argc; ++i)
        m_args.push_back(strip_whitespace(argv[i]));

    // Initialise application parameters
    m_pars.load(parfilename(), m_args);
    
    // DUMMY FOR TESTING
    std::cout << m_pars << std::endl;
    std::string value = m_pars.par("chatter")->value();
    std::cout << value << std::endl;
    value = m_pars.par("clobber")->value();
    std::cout << value << std::endl;

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
    // Save application parameters
    m_pars.save(parfilename());

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
    // Initialise members
    m_name.clear();
    m_version.clear();
    m_args.clear();
    m_tstart = 0;
  
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
    // Copy attributes
    m_name    = app.m_name;
    m_version = app.m_version;
    m_args    = app.m_args;
    m_tstart  = app.m_tstart;
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GApplication::free_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Returns parameter filename
 ***************************************************************************/
std::string GApplication::parfilename(void) const
{
    // Return
    return (m_name+".par");
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


