/***************************************************************************
 *             GApplication.hpp - GammaLib application base class          *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2011 by Jurgen Knodlseder                           *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GApplication.hpp
 * @brief GammaLib application base class
 * @author Jurgen Knodlseder
 */

#ifndef GAPPLICATION_HPP
#define GAPPLICATION_HPP

/* __ Includes ___________________________________________________________ */
#include <ctime>
#include <vector>
#include <string>
#include <iostream>
#include "GLog.hpp"
#include "GPars.hpp"


/***********************************************************************//**
 * @class GApplication
 *
 * @brief GammaLib application interface defintion.
 ***************************************************************************/
class GApplication {

    // I/O friends
    friend std::ostream& operator<< (std::ostream& os, const GApplication& app);

public:
    // Constructors and destructors
    GApplication(void);
    GApplication(const std::string& name, const std::string& version);
    GApplication(const std::string& name, const std::string& version,
                 int argc, char* argv[]);
    GApplication(const GApplication& app);
    ~GApplication(void);

    // Operators
    GApplication& operator= (const GApplication& app);

    // Methods
    std::string name(void) const;
    std::string version(void) const;
    double      telapse(void) const;
    GPar*       par(const std::string& name);
    bool        logTerse(void) const;
    bool        logNormal(void) const;
    bool        logExplicit(void) const;
    bool        logVerbose(void) const;
    bool        logDebug(void) const;
    bool        clobber(void) const { return m_clobber; }

    // Public members
    GLog        log;                       //!< Application logger
    
protected:
    // Protected methods
    void        init_members(void);
    void        copy_members(const GApplication& app);
    void        free_members(void);
    std::string par_filename(void) const;
    std::string log_filename(void) const;
    void        get_par_standard(void);
    void        log_header(void);
    void        log_trailer(void);
    void        log_parameters(void);

    // Protected data members
    std::string              m_name;       //!< Application name
    std::string              m_version;    //!< Application version
    std::vector<std::string> m_args;       //!< Command line arguments
    std::time_t              m_tstart;     //!< Start time of execution
    GPars                    m_pars;       //!< Application parameters
    int                      m_chatter;    //!< Chatter level
    bool                     m_clobber;    //!< Clobber flag
    bool                     m_debug;      //!< Debugging mode activated
};

#endif /* GAPPLICATION_HPP */
