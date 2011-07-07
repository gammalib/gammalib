/***************************************************************************
 *             GApplication.hpp - GammaLib application base class          *
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
 * @brief GammaLib application interface defintion
 *
 * This class provides the base class for ftools-like executables based on
 * GammaLib. The ftools-like executables will be implemented as derived
 * classes, and automatically have access to task parameters and the
 * application logger.
 ***************************************************************************/
class GApplication {

    // I/O friends
    friend std::ostream& operator<< (std::ostream& os, const GApplication& app);
    friend GLog&         operator<< (GLog&        log, const GApplication& app);

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
    GPar&         operator[](const std::string& name);
    const GPar&   operator[](const std::string& name) const;

    // Methods
    std::string name(void) const;
    std::string version(void) const;
    double      telapse(void) const;
    double      celapse(void) const;
    void        logFileOpen(bool clobber = true);
    bool        logTerse(void) const;
    bool        logNormal(void) const;
    bool        logExplicit(void) const;
    bool        logVerbose(void) const;
    bool        logDebug(void) const;
    bool        clobber(void) const;
    bool        haspar(const std::string& name) const;
    std::string print(void) const;

    // Public members
    GLog        log;                       //!< Application logger
    
protected:
    // Protected methods
    void        init_members(void);
    void        copy_members(const GApplication& app);
    void        free_members(void);
    std::string par_filename(void) const;
    std::string log_filename(void) const;
    void        log_header(void);
    void        log_trailer(void);
    void        log_parameters(void);

    // Protected data members
    std::string              m_name;       //!< Application name
    std::string              m_version;    //!< Application version
    std::vector<std::string> m_args;       //!< Command line arguments
    std::time_t              m_tstart;     //!< Calendar start time of execution
    std::clock_t             m_cstart;     //!< Clock start time of execution
    GPars                    m_pars;       //!< Application parameters
};

#endif /* GAPPLICATION_HPP */
