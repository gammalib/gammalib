/***************************************************************************
 *             GApplication.hpp - GammaLib application base class          *
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
 * @file GApplication.hpp
 * @brief GammaLib application base class
 * @author Juergen Knoedlseder
 */

#ifndef GAPPLICATION_HPP
#define GAPPLICATION_HPP

/* __ Includes ___________________________________________________________ */
#include <ctime>
#include <vector>
#include <string>
#include "GBase.hpp"
#include "GLog.hpp"
#include "GPars.hpp"


/***********************************************************************//**
 * @class GApplication
 *
 * @brief GammaLib application interface definition
 *
 * This class provides the base class for ftools-like executables based on
 * GammaLib. The ftools-like executables will be implemented as derived
 * classes, and automatically have access to task parameters and the
 * application logger.
 ***************************************************************************/
class GApplication : public GBase {

public:
    // Constructors and destructors
    GApplication(void);
    GApplication(const std::string& name, const std::string& version);
    GApplication(const std::string& name, const std::string& version,
                 int argc, char* argv[]);
    GApplication(const GApplication& app);
    ~GApplication(void);

    // Operators
    GApplication& operator=(const GApplication& app);
    GPar&         operator[](const std::string& name);
    const GPar&   operator[](const std::string& name) const;

    // Methods
    void               clear(void);
    GApplication*      clone(void) const;
    const std::string& name(void) const;
    const std::string& version(void) const;
    double             telapse(void) const;
    double             celapse(void) const;
    void               logFileOpen(const bool& clobber = true);
    bool               logTerse(void) const;
    bool               logNormal(void) const;
    bool               logExplicit(void) const;
    bool               logVerbose(void) const;
    bool               logDebug(void) const;
    bool               clobber(void) const;
    bool               has_par(const std::string& name) const;
    const std::string& par_filename(void) const;
    const std::string& log_filename(void) const;
    void               log_header(void);
    void               log_trailer(void);
    void               log_parameters(void);
    std::string        print(const GChatter& chatter = NORMAL) const;

    // Public members
    GLog log;   //!< Application logger

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GApplication& app);
    void free_members(void);
    void set_log_chatter(void);
    void set_log_filename(void);

    // Protected data members
    std::string              m_name;       //!< Application name
    std::string              m_version;    //!< Application version
    std::string              m_parfile;    //!< Parameter filename
    std::string              m_logfile;    //!< Log filename
    std::vector<std::string> m_args;       //!< Command line arguments
    std::time_t              m_tstart;     //!< Calendar start time of execution
    std::clock_t             m_cstart;     //!< Clock start time of execution
    GPars                    m_pars;       //!< Application parameters
};


/***********************************************************************//**
 * @brief Parameter access operator
 *
 * @param[in] name Parameter name.
 ***************************************************************************/
inline
GPar& GApplication::operator[](const std::string& name)
{
    return (m_pars[name]);
}


/***********************************************************************//**
 * @brief Parameter access operator (const version)
 *
 * @param[in] name Parameter name.
 ***************************************************************************/
inline
const GPar& GApplication::operator[](const std::string& name) const
{
    return (m_pars[name]);
}


/***********************************************************************//**
 * @brief Return application name
 *
 * @return Application name.
 ***************************************************************************/
inline
const std::string& GApplication::name(void) const
{
    return m_name;
}


/***********************************************************************//**
 * @brief Return application version
 *
 * @return Application version.
 ***************************************************************************/
inline
const std::string& GApplication::version(void) const
{
    return m_version;
}


/***********************************************************************//**
 * @brief Signal if specified parameter exists
 *
 * @param[in] name Parameter name.
 * @return True if an application parameter with the specified name exists.
 ***************************************************************************/
inline
bool GApplication::has_par(const std::string& name) const
{
    return (m_pars.contains(name));
}


/***********************************************************************//**
 * @brief Returns parameter filename
 *
 * @return Parameter filename.
 ***************************************************************************/
inline
const std::string& GApplication::par_filename(void) const
{
    // Return
    return (m_parfile);
}


/***********************************************************************//**
 * @brief Returns log filename
 *
 * @return Log filename.
 ***************************************************************************/
inline
const std::string& GApplication::log_filename(void) const
{
    // Return
    return (m_logfile);
}

#endif /* GAPPLICATION_HPP */
