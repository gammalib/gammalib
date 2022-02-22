/***************************************************************************
 *                        GDaemon.hpp - Daemon class                       *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2022 by Juergen Knoedlseder                              *
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
 * @file GDaemon.hpp
 * @brief Daemon class definition
 * @author Juergen Knoedlseder
 */

#ifndef GDAEMON_HPP
#define GDAEMON_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GBase.hpp"
#include "GTools.hpp"
#include "GFilename.hpp"
#include "GLog.hpp"

/* __ Forward declarations _______________________________________________ */
class GCsv;
class GXml;

/* __ Constants __________________________________________________________ */


/***********************************************************************//**
 * @class GDaemon
 *
 * @brief Daemon class
 *
 * @todo Add class description.
 ***************************************************************************/
class GDaemon : public GBase {

public:
    // Constructors and destructors
    GDaemon(void);
    GDaemon(const GDaemon& daemon);
    virtual ~GDaemon(void);

    // Operators
    GDaemon& operator=(const GDaemon& daemon);

    // Implemented pure virtual base class methods
    virtual void        clear(void);
    virtual GDaemon*    clone(void) const;
    virtual std::string classname(void) const;
    virtual std::string print(const GChatter& chatter = NORMAL) const;

    // Other methods
    void start(void);
    bool alive(void) const;

protected:
    // Protected methods
    void      init_members(void);
    void      copy_members(const GDaemon& daemon);
    void      free_members(void);
    void      create_lock_file(void);
    void      delete_lock_file(void);
    void      write_heartbeat(void);
    GFilename lock_filename(void) const;
    GFilename heartbeat_filename(void) const;
    GFilename statistics_filename(void) const;
    pid_t     lock_pid(void) const;
    void      update_statistics(void);
    void      update_xml(const GCsv& statistics);
    void      create_xml(const GFilename& filename);
    void      update_dates(GXml& xml, const GCsv& statistics);
    void      update_countries_header(GXml& xml, const GCsv& statistics);
    void      update_countries_data(GXml& xml, const GCsv& statistics);
    void      update_versions_data(GXml& xml, const GCsv& statistics);
    void      update_daily(GXml& xml, const GCsv& statistics);

    // Protected members
    pid_t    m_pid;        //!< Process ID
    int      m_period;     //!< Wake-up period in seconds
    int      m_heartbeat;  //!< Heartbeat period in seconds
    GLog     m_log;        //!< Logger
    GChatter m_chatter;    //!< Chattiness of logger
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GDaemon").
 ***************************************************************************/
inline
std::string GDaemon::classname(void) const
{
    return ("GDaemon");
}


/***********************************************************************//**
 * @brief Returns name of daemon lock file
 *
 * @return Daemon lock filename.
 ***************************************************************************/
inline
GFilename GDaemon::lock_filename(void) const
{
    return (gammalib::gamma_filename("daemon.lock"));
}


/***********************************************************************//**
 * @brief Returns name of daemon heartbeat file
 *
 * @return Daemon heartbeat filename.
 ***************************************************************************/
inline
GFilename GDaemon::heartbeat_filename(void) const
{
    return (gammalib::gamma_filename("daemon.heartbeat"));
}

#endif /* GDAEMON_HPP */
