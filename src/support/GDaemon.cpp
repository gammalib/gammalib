/***************************************************************************
 *                        GDaemon.cpp - Daemon class                       *
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
 * @file GDaemon.cpp
 * @brief Daemon class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <cstdlib>         // std::getenv()
#include <cstdio>          // std::fopen()
#include <unistd.h>        // sleep()
#include <fcntl.h>         // for file locking
#include <signal.h>        // kill()
#if defined(__linux__) || defined(__linux) || defined(linux) || defined(__gnu_linux__)
#include <sys/prctl.h>     // prctl(), PR_SET_NAME
#endif
#include "GDaemon.hpp"
#include "GCsv.hpp"
#include "GTime.hpp"
#include "GXml.hpp"
#include "GXmlNode.hpp"

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
GDaemon::GDaemon(void)
{
    // Initialise class members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] daemon Daemon.
 ***************************************************************************/
GDaemon::GDaemon(const GDaemon& daemon)
{
    // Initialise class members
    init_members();

    // Copy members
    copy_members(daemon);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GDaemon::~GDaemon(void)
{
    // Free members
    free_members();

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                               Operators                                 =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Assignment operator
 *
 * @param[in] daemon Daemon.
 * @return Daemon.
 ***************************************************************************/
GDaemon& GDaemon::operator=(const GDaemon& daemon)
{
    // Execute only if object is not identical
    if (this != &daemon) {

        // Free members
        free_members();

        // Initialise private members
        init_members();

        // Copy members
        copy_members(daemon);

    } // endif: object was not identical

    // Return this object
    return *this;
}


/*==========================================================================
 =                                                                         =
 =                             Public methods                              =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Clear Daemon
 ***************************************************************************/
void GDaemon::clear(void)
{
    // Free members
    free_members();

    // Initialise private members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone Daemon
 *
 * @return Pointer to deep copy of Daemon.
 ***************************************************************************/
GDaemon* GDaemon::clone(void) const
{
    return new GDaemon(*this);
}


/***********************************************************************//**
 * @brief Starts the daemon
 *
 * Starts up the daemon and entered the never-ending daemon event loop.
 ***************************************************************************/
void GDaemon::start(void)
{
    // Open logger
    m_log.open(gammalib::gamma_filename("daemon.log"), true);
    m_log.chatter((GChatter)m_chatter);
    m_log.date(true);

    // Log start up of daemon
    m_log << "GammaLib daemon start up" << std::endl;

    // Set daemon name
    #if defined(__linux__) || defined(__linux) || defined(linux) || defined(__gnu_linux__)
    prctl(PR_SET_NAME, "gammalibd");
    #endif

    // Create lock file
    create_lock_file();

    // Main event handling loop
    while (1) {

        // Update application statistics
        update_statistics();

        // Force logger flushing
        m_log.flush(true);

        // Wait for some period
        sleep(m_period);

    } // endwhile: main event loop

    // Delete lock file
    delete_lock_file();

    // Close log file
    m_log.close();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Check if daemon is alive
 *
 * @return True if daemon is alive
 *
 * Checks if the daemon is alive. The daemon is alive if a lock file exists
 * and if the process with the PID that is found in the lock file is alive.
 *
 * This method can be used by a client to check whether a daemon is alive.
 ***************************************************************************/
bool GDaemon::alive(void) const
{
    // Initialise flag
    bool alive = false;

    // Get process ID in lock file
    pid_t pid = lock_pid();

    // If process ID is positive then check if process is still alive
    if (pid > 0) {

        // Check if process is alive
        if (0 == kill(pid, 0)) {
            alive = true;
        }
    
    } // endif: process ID was positive

    // Return flag
    return alive;
}


/***********************************************************************//**
 * @brief Print Daemon
 *
 * @param[in] chatter Chattiness.
 * @return String containing Daemon information.
 ***************************************************************************/
std::string GDaemon::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GDaemon ===");

        // Append information
        result.append("\n"+gammalib::parformat("Process ID"));
        result.append(gammalib::str(m_pid));
        result.append("\n"+gammalib::parformat("Wake up period (s)"));
        result.append(gammalib::str(m_period));
        result.append("\n"+gammalib::parformat("Logger chatter level"));
        result.append(gammalib::str(m_chatter));

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
void GDaemon::init_members(void)
{
    // Initialise members
    m_pid     = 0;       //!< No process ID
    m_period  = 3600.0;  //!< Wake-up daemon every hour
    m_log.clear();
    m_chatter = NORMAL;  //!< NORMAL chatter level
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] daemon Daemon.
 ***************************************************************************/
void GDaemon::copy_members(const GDaemon& daemon)
{
    // Copy members
    m_pid     = daemon.m_pid;
    m_period  = daemon.m_period;
    m_log     = daemon.m_log;
    m_chatter = daemon.m_chatter;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GDaemon::free_members(void)
{
    // Delete lock file
    delete_lock_file();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Create the lock file
 *
 * Creates the daemon lock file and write process ID into lock file.
 ***************************************************************************/
void GDaemon::create_lock_file(void)
{
    // Get lock filename
    GFilename lockfile = lock_filename();

    // Open lock file. Continue only if opening was successful.
    FILE* fptr = fopen(lockfile.url().c_str(), "w");

    // Successful?
    if (fptr != NULL) {

        // Get process ID
        m_pid = getpid();

        // Write process ID into lockfile
        fprintf(fptr, "%d\n", m_pid);

        // Close lockfile
        fclose(fptr);

        // Log creation of lock file
        m_log << "Created lock file \"" << lockfile.url();
        m_log << "\" for PID " << m_pid << std::endl;

    } // endif: Lock file opened successfully

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete daemon lock file
 *
 * Deletes the daemon lock file on disk if the process ID in the lock file
 * corresponds to the process ID of the instance.
 ***************************************************************************/
void GDaemon::delete_lock_file(void) const
{
    // If process ID in lock file is the ID of the current process then
    // delete the lock file
    if (lock_pid() == getpid()) {
    
        // Get lock filename
        GFilename lockfile = lock_filename();

        // Delete lock file
        std::remove(lockfile.url().c_str());

    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Returns process ID in lock file
 *
 * @return Process ID in lock file.
 *
 * Returns the process ID that is found in the lock file. If no lock file
 * exists the method returns 0.
 ***************************************************************************/
pid_t GDaemon::lock_pid(void) const
{
    // Initialis process ID
    pid_t pid = 0;

    // Get lock filename
    GFilename lockfile = lock_filename();

    // Open lock file. Continue only if opening was successful.
    FILE* fptr = fopen(lockfile.url().c_str(), "r");
    if (fptr != NULL) {

        // Allocate line buffer
        const int n = 1000;
        char      line[n];

        // Read line
        fgets(line, n-1, fptr);

        // Close file
        fclose(fptr);

        // Extract process ID
        pid = gammalib::toint(std::string(line));

    } // endif: lock file opened successfully

    // Return process ID
    return pid;
}


/***********************************************************************//**
 * @brief Update application statistics
 *
 * Update application statistics in the @p statistics.xml file by scanning
 * the low-level file @p statistics.csv.
 *
 * @todo Delete low-level file if scan was successful.
 ***************************************************************************/
void GDaemon::update_statistics(void)
{
    // Set low-level statistics filename
    GFilename filename = gammalib::gamma_filename("statistics.csv");

    // Continue only if file exists
    if (filename.exists()) {

        // OpenMP critical zone to write statitics
        #pragma omp critical(GDaemon_update_statistics)
        {

            // Get file lock. Continue only in case of success.
            struct flock lock;
            lock.l_type   = F_RDLCK;  // Want a read lock
            lock.l_whence = SEEK_SET; // Want beginning of file
            lock.l_start  = 0;        // No offset, lock entire file ...
            lock.l_len    = 0;        // ... to the end
            lock.l_pid    = getpid(); // Current process ID
            int fd        = open(filename.url().c_str(), O_RDONLY);
            if (fd != -1) {

                // Lock file
                fcntl(fd, F_SETLKW, &lock);

                // Load CSV file
                GCsv statistics(filename, ",");

                // Update high-level statistics
                update_high_level_statistics(statistics);

                // Remove low-level statistics file
                std::remove(filename.url().c_str());

                // Unlock file
                lock.l_type = F_UNLCK;
                fcntl(fd, F_SETLK, &lock);

                // Close file
                close(fd);

            } // endif: file locking successful
        
        } // end of OMP critial zone

    } // endif: there was a statistics ASCII file

    // Return
    return;
}


/***********************************************************************//**
 * @brief Update high-level statistics
 *
 * Update high-level statistics by scanning the low-level statistics file.
 * The high-level statistics is written out as XML file. If no XML file
 * exists a new one is created, otherwise an existing file is updated.
 ***************************************************************************/
void GDaemon::update_high_level_statistics(const GCsv& statistics)
{
    // Log updating of high-level statistics
    m_log << "Update high-level statistics" << std::endl;

    // Set current time
    GTime now;
    now.now();

    // Set high-level statistics filename
    GFilename filename = gammalib::gamma_filename("statistics.xml");

    // Initialise high-level statistics
    GXml xml;

    // If a high-level statistics file exists already then load information
    // from the file
    if (filename.exists()) {
        xml.load(filename);
    }

    // ... otherwise create empty high-level statistics file
    else {

        // Append base
        GXmlNode* base = xml.append(GXmlElement("statistics title=\"High-level statistic\""));

        // Append header
        GXmlNode* header = base->append("header");
        GXmlNode* dates  = header->append("dates");
        GXmlNode* d1     = dates->append("creation");
        d1->append(GXmlText(now.utc()));
        GXmlNode* d2     = dates->append("modified");
        d2->append(GXmlText(now.utc()));
        GXmlNode* d3     = dates->append("start");
        d3->append(GXmlText(""));
        GXmlNode* d4     = dates->append("stop");
        d4->append(GXmlText(""));

        // Append data
        GXmlNode* data  = base->append("data");
        GXmlNode* daily = data->append("daily");

        // Log creation of high-level statistics file
        m_log << "Created high-level statistics XML file \"";
        m_log << filename.url() << "\"" << std::endl;

    } // endelse: created empty high-level statistics file

    // Get some useful pointers
    GXmlNode* base   = xml.element("statistics",0);
    GXmlNode* header = base->element("header",0);
    GXmlNode* data   = base->element("data",0);
    GXmlNode* daily  = data->element("daily",0);

    // Loop over statistics
    for (int i = 0; i < statistics.nrows(); ++i) {

        // Skip header row
        if (i == 0) {
            continue;
        }

        // Extract relevant attributes
        std::string date    = statistics.string(i,0).substr(0,10);
        std::string tool    = statistics.string(i,3);
        std::string version = statistics.string(i,4);
        double      wall    = statistics.real(i,5);
        double      cpu     = statistics.real(i,6);
        double      gCO2e   = statistics.real(i,7);

        // Skip tools without name
        if (tool.empty()) {
            continue;
        }

        // Get pointer to relevant date element. If no date element exists
        // that corresponds to the date of the statistics then append an
        // element
        GXmlNode* n_date = NULL;
        int ndates = daily->elements("date");
        if (ndates == 0) {
            n_date = daily->append("date");
            n_date->append(GXmlElement("value", date));
        }
        else {
            for (int k = 0; k < ndates; ++k) {
                GXmlElement* n = daily->element("date", k);
                if (daily->element("date", k)->element("value",0)->value() == date) {
                    n_date = n;
                    break;
                }
            }
            if (n_date == NULL) {
                n_date = daily->append("date");
                n_date->append(GXmlElement("value", date));
            }
        }

        // Get pointer to "tools" node. If no such node exists then append
        // one
        GXmlNode* n_tools = NULL;
        if (n_date->elements("tools") == 0) {
            n_tools = n_date->append("tools");
        }
        else {
            n_tools = n_date->element("tools",0);
        }

        // Get pointer to relevant tool element. If no tool element exists
        // that corresponds to the relevant tool then append an element
        GXmlElement* n_tool = NULL;
        int ntools = n_tools->elements();
        if (ntools == 0) {
            n_tool = n_tools->append(tool);
        }
        else {
            for (int k = 0; k < ntools; ++k) {
                GXmlElement* n = n_tools->element(k);
                if (n->name() == tool) {
                    n_tool = n;
                    break;
                }
            }
            if (n_tool == NULL) {
                n_tool = n_tools->append(tool);
            }
        }

        // Set or update attributes
        int calls = 1;
        if (!n_tool->has_attribute("version")) {
            n_tool->attribute("version", version);
        }
        if (n_tool->has_attribute("calls")) {
            calls += gammalib::toint(n_tool->attribute("calls"));
        }
        if (n_tool->has_attribute("wall")) {
            wall += gammalib::todouble(n_tool->attribute("wall"));
        }
        if (n_tool->has_attribute("cpu")) {
            cpu += gammalib::todouble(n_tool->attribute("cpu"));
        }
        if (n_tool->has_attribute("gCO2e")) {
            gCO2e += gammalib::todouble(n_tool->attribute("gCO2e"));
        }
        n_tool->attribute("calls",  gammalib::str(calls));
        n_tool->attribute("wall",   gammalib::str(wall));
        n_tool->attribute("cpu",    gammalib::str(cpu));
        n_tool->attribute("gCO2e",  gammalib::str(gCO2e));

    } // endfor: looped over statistics

    // Update modified time in header
    GXmlNode* modified = header->element("dates > modified");
    static_cast<GXmlText*>((*modified)[0])->text(now.utc());

    // Update first and last time in header
    std::string start;
    std::string stop;
    if (statistics.nrows() > 1) {
        start                     = statistics.string(1,0);
        stop                      = statistics.string(statistics.nrows()-1,0);
        GXmlNode*   n_start       = header->element("dates > start");
        GXmlNode*   n_stop        = header->element("dates > stop");
        std::string current_start = static_cast<GXmlText*>((*n_start)[0])->text();
        std::string current_stop  = static_cast<GXmlText*>((*n_stop)[0])->text();
        if (current_start.empty() || (start < current_start)) {
            static_cast<GXmlText*>((*n_start)[0])->text(start);
        }
        if (current_stop.empty() || (stop > current_stop)) {
            static_cast<GXmlText*>((*n_stop)[0])->text(stop);
        }
    }

    // Save high-level statistics file
    xml.save(filename);

    // Return
    return;
}
