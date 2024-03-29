/***************************************************************************
 *             GApplication.hpp - GammaLib application base class          *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2023 by Juergen Knoedlseder                         *
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
#include "GApplicationPars.hpp"

/* __ Forward declarations _______________________________________________ */
class GFitsHDU;
class GFits;
class GFilename;


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
    GApplication(const std::string& name,
                 const std::string& version);
    GApplication(const std::string&      name,
                 const std::string&      version,
                 const GApplicationPars& pars);
    GApplication(const std::string& name,
                 const std::string& version,
                 int                argc,
                 char*              argv[]);
    GApplication(const GApplication& app);
    ~GApplication(void);

    // Operators
    GApplication&          operator=(const GApplication& app);
    GApplicationPar&       operator[](const int& index);
    const GApplicationPar& operator[](const int& index) const;
    GApplicationPar&       operator[](const std::string& name);
    const GApplicationPar& operator[](const std::string& name) const;

    // Public methods
    void                            clear(void);
    GApplication*                   clone(void) const;
    std::string                     classname(void) const;
    const std::string&              name(void) const;
    const std::string&              version(void) const;
    double                          telapse(void) const;
    double                          celapse(void) const;
    double                          gCO2e(const std::string& country) const;
    void                            add_celapse(const double& celapse);
    void                            logFileOpen(const bool& clobber = true);
    void                            logFileClose(void);
    bool                            logTerse(void) const;
    bool                            logNormal(void) const;
    bool                            logExplicit(void) const;
    bool                            logVerbose(void) const;
    bool                            logDebug(void) const;
    bool                            clobber(void) const;
    bool                            has_par(const std::string& name) const;
    const std::string&              par_filename(void) const;
    const std::string&              log_filename(void) const;
    void                            log_header(void);
    void                            log_trailer(void);
    void                            log_string(const GChatter& chatter,
                                               const std::string& string,
                                               const bool&        linefeed = true);
    void                            log_value(const GChatter&    chatter,
                                              const std::string& name,
                                              const std::string& value,
                                              const std::string& unit = "");
    void                            log_value(const GChatter&    chatter,
                                              const std::string& name,
                                              const int&         value,
                                              const std::string& unit = "");
    void                            log_value(const GChatter&    chatter,
                                              const std::string& name,
                                              const double&      value,
                                              const std::string& unit = "");
    void                            log_header1(const GChatter&    chatter,
                                                const std::string& header);
    void                            log_header2(const GChatter&    chatter,
                                                const std::string& header);
    void                            log_header3(const GChatter&    chatter,
                                                const std::string& header);
    void                            log_parameters(const GChatter& chatter);
    const bool&                     need_help(void) const;
    void                            statistics(const bool& statistics);
    const bool&                     statistics(void) const;
    const GApplicationPars&         pars(void) const;
    void                            pars(const GApplicationPars& pars);
    const std::vector<std::string>& args(void) const;
    void                            stamp(GFitsHDU& hdu) const;
    void                            stamp(GFits& fits) const;
    void                            stamp(const GFilename& filename) const;
    std::string                     print(const GChatter& chatter = NORMAL) const;

    // Public members
    GLog log;   //!< Application logger

#ifndef SWIG
protected:
#endif
    // Returns the number of running instance of this method. The number
    // of running instances has been implement as static method to avoid
    // the static initialization order fiasco of static members; using
    // static methods we follow the "construct on first use idiom")
    static int& running() {
        static int m_running = 0;
        return m_running;
    }

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GApplication& app);
    void free_members(void);
    void set_statistics(void);
    void set_log_chatter(void);
    void set_log_filename(void);
    void write_statistics(void);

    // Protected data members
    std::string              m_name;        //!< Application name
    std::string              m_version;     //!< Application version
    std::string              m_parfile;     //!< Parameter filename
    std::string              m_logfile;     //!< Log filename
    std::vector<std::string> m_args;        //!< Command line arguments
    std::time_t              m_tstart;      //!< Calendar start time of execution
    double                   m_cstart;      //!< Clock start time of execution
    double                   m_celapse;     //!< Internal CPU seconds counter
    GApplicationPars         m_pars;        //!< Application parameters
    bool                     m_pars_loaded; //!< Application parameters loaded
    bool                     m_need_help;   //!< --help specified
    bool                     m_statistics;  //!< Enable writing of statistics
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GApplication").
 ***************************************************************************/
inline
std::string GApplication::classname(void) const
{
    return ("GApplication");
}


/***********************************************************************//**
 * @brief Parameter access operator
 *
 * @param[in] name Parameter name.
 * @return Reference to application parameter.
 *
 * Returns a reference to the application parameter with the given @p name.
 ***************************************************************************/
inline
GApplicationPar& GApplication::operator[](const std::string& name)
{
    return (m_pars[name]);
}


/***********************************************************************//**
 * @brief Parameter access operator (const version)
 *
 * @param[in] name Parameter name.
 * @return Constant reference to application parameter
 *
 * Returns a const reference to the application parameter with the given
 * @p name.
 ***************************************************************************/
inline
const GApplicationPar& GApplication::operator[](const std::string& name) const
{
    return (m_pars[name]);
}


/***********************************************************************//**
 * @brief Parameter access operator
 *
 * @param[in] index Parameter index [0,...,pars().size()-1].
 * @return Reference to application parameter
 *
 * Returns a reference to the application parameter with the given @p index.
 * No range checking is performed for the index.
 ***************************************************************************/
inline
GApplicationPar& GApplication::operator[](const int& index)
{
    return (m_pars[index]);
}


/***********************************************************************//**
 * @brief Parameter access operator (const version)
 *
 * @param[in] index Parameter index [0,...,pars().size()-1].
 * @return Constant reference to application parameter
 *
 * Returns a const reference to the application parameter with the given
 * @p index. No range checking is performed for the index.
 ***************************************************************************/
inline
const GApplicationPar& GApplication::operator[](const int& index) const
{
    return (m_pars[index]);
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
 * @brief Add CPU seconds to internal counter
 *
 * @param[in] celapse CPU seconds.
 *
 * Adds the @p celapse CPU seconds to the internal CPU seconds counter.
 ***************************************************************************/
inline
void GApplication::add_celapse(const double& celapse)
{
    m_celapse += celapse;
    return;
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
    // Set logfile name from parameters (in case it has changed)
    const_cast<GApplication*>((this))->set_log_filename();

    // Return
    return (m_logfile);
}


/***********************************************************************//**
 * @brief Signals if --help option has been specified
 *
 * @return True if --help option has been specified.
 ***************************************************************************/
inline
const bool& GApplication::need_help(void) const
{
    // Return
    return (m_need_help);
}


/***********************************************************************//**
 * @brief Set if statistics will be written
 *
 * @param[in] statistics Write statistics?
 ***************************************************************************/
inline
void GApplication::statistics(const bool& statistics)
{
    m_statistics = statistics;
    return;
}


/***********************************************************************//**
 * @brief Signals if statistics will be written
 *
 * @return True if statistics will be written.
 ***************************************************************************/
inline
const bool& GApplication::statistics(void) const
{
    // Return
    return (m_statistics);
}


/***********************************************************************//**
 * @brief Return application parameters
 *
 * @return Application parameters.
 ***************************************************************************/
inline
const GApplicationPars& GApplication::pars(void) const
{
    return m_pars;
}


/***********************************************************************//**
 * @brief Set application parameters
 *
 * @param[in] pars Application parameters.
 ***************************************************************************/
inline
void GApplication::pars(const GApplicationPars& pars)
{
    m_pars = pars;
    return;
}


/***********************************************************************//**
 * @brief Return command line arguments
 *
 * @return Command line arguments.
 ***************************************************************************/
inline
const std::vector<std::string>& GApplication::args(void) const
{
    return m_args;
}

#endif /* GAPPLICATION_HPP */
