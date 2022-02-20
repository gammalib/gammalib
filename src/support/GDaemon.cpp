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
        result.append(gammalib::str((int)(m_pid)));
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
    m_period  = 3600;    //!< Wake-up daemon every hour
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
        fprintf(fptr, "%d\n", (int)(m_pid));

        // Close lockfile
        fclose(fptr);

        // Log creation of lock file
        m_log << "Created lock file \"" << lockfile.url();
        m_log << "\" for PID " << (int)(m_pid) << std::endl;

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
                update_xml(statistics);

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
void GDaemon::update_xml(const GCsv& statistics)
{
    // Log updating of high-level statistics
    m_log << "Update high-level statistics" << std::endl;

    // Set high-level statistics filename
    GFilename filename = gammalib::gamma_filename("statistics.xml");

    // If no high-level statistics file exists then create one
    if (!filename.exists()) {
        create_xml(filename);
    }

    // Load high-level statistics file. If an exception occurs during this
    // process the input file is probably corrupt, and a failure message
    // is written into the log file. Processing only continues if no
    // exception occured.
    GXml xml;
    try {
        xml.load(filename);
    }
    catch (...) {
        xml.clear();
        m_log << "*** Failure occured in loading high-level statistics ";
        m_log << "XML file \"" << filename.url() << "\"" << std::endl;
    }

    // Update file only if XML document is not empty
    if (!xml.is_empty()) {

        // Update dates in header
        update_dates(xml, statistics);
        
        // Update countries in header and data
        update_countries_header(xml, statistics);
        update_countries_data(xml, statistics);

        // Update versions in data
        update_versions_data(xml, statistics);

        // Update daily information
        update_daily(xml, statistics);

        // Save updated high-level statistics XML file
        xml.save(filename);

    } // endif: XML document was not empty

    // Return
    return;
}


/***********************************************************************//**
 * @brief Create high-level statistics XML file
 *
 * @param[in] filename File name.
 *
 * Creates an empty high-level statistics XML file with the specified
 * @p filename. This method overwrites any existing XML file.
 ***************************************************************************/
void GDaemon::create_xml(const GFilename& filename)
{
    // Initialise high-level statistics XML object
    GXml xml;

    // Append base node
    GXmlNode* base = xml.append(GXmlElement("statistics title=\"High-level statistics\""));

    // Append header
    GXmlNode* header    = base->append("header");
    GXmlNode* dates     = header->append("dates");
    GXmlNode* countries = header->append("countries");
    GXmlNode* d1        = dates->append("creation");
    d1->append(GXmlText(""));
    GXmlNode* d2 = dates->append("modified");
    d2->append(GXmlText(""));
    GXmlNode* d3 = dates->append("start");
    d3->append(GXmlText(""));
    GXmlNode* d4 = dates->append("stop");
    d4->append(GXmlText(""));

    // Append data
    GXmlNode* data  = base->append("data");
    data->append("countries");
    data->append("versions");
    data->append("daily");

    // Log creation of high-level statistics file
    m_log << "Created high-level statistics XML file \"";
    m_log << filename.url() << "\"" << std::endl;

    // Save high-level statistics XML file
    xml.save(filename);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Update dates in high-level statistics header.
 *
 * @param[in,out] xml High-level statistics XML object.
 * @param[in] statistics Low-level statistics CSV object.
 *
 * Updates dates in the high-level statistics header.
 ***************************************************************************/
void GDaemon::update_dates(GXml& xml, const GCsv& statistics)
{
    // Set current time
    GTime now;
    now.now();

    // Get pointer to "header" node
    GXmlNode* base   = xml.element("statistics",0);
    GXmlNode* header = base->element("header",0);

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

    // Return
    return;
}


/***********************************************************************//**
 * @brief Update countries in high-level statistics header.
 *
 * @param[in,out] xml High-level statistics XML object.
 * @param[in] statistics Low-level statistics CSV object.
 *
 * Update country statistics in header block.
 ***************************************************************************/
void GDaemon::update_countries_header(GXml& xml, const GCsv& statistics)
{
    // Get pointer to "countries" node in "header" node. In case that no
    // "countries" node exists in "header" node then append one.
    GXmlNode* base   = xml.element("statistics",0);
    GXmlNode* header = base->element("header",0);
    if (header->elements("countries") == 0) {
        header->append("countries");
    }
    GXmlNode* countries = header->element("countries",0);

    // Loop over statistics
    for (int i = 1; i < statistics.nrows(); ++i) {

        // Extract relevant attributes
        std::string country = statistics.string(i,2);

        // Update list of countries in header
        int ncountries = countries->elements("country");
        if (ncountries == 0) {
            GXmlNode* n = countries->append("country");
            n->append(GXmlText(country));
        }
        else {
            bool found = false;
            for (int k = 0; k < ncountries; ++k) {
                if (countries->element("country", k)->value() == country) {
                    found = true;
                    break;
                }
            }
            if (!found) {
                GXmlNode* n = countries->append("country");
                n->append(GXmlText(country));
            }
        }

    } // endfor: looped over statistics

    // Return
    return;
}


/***********************************************************************//**
 * @brief Update countries in high-level statistics data.
 *
 * @param[in,out] xml High-level statistics XML object.
 * @param[in] statistics Low-level statistics CSV object.
 *
 * Update country statistics in data block.
 ***************************************************************************/
void GDaemon::update_countries_data(GXml& xml, const GCsv& statistics)
{
    // Get pointer to "countries" node. In case that no "countries" node
    // exists in "data" node then append one.
    GXmlNode* base = xml.element("statistics",0);
    GXmlNode* data = base->element("data",0);
    if (data->elements("countries") == 0) {
        data->append("countries");
    }
    GXmlNode* countries = data->element("countries",0);

    // Loop over statistics
    for (int i = 1; i < statistics.nrows(); ++i) {

        // Extract relevant attributes
        std::string country = statistics.string(i,2);
        double      wall    = statistics.real(i,5);
        double      cpu     = statistics.real(i,6);
        double      gCO2e   = statistics.real(i,7);

        // Get pointer to relevant "country" element. If no "country"
        // element exists that corresponds to the country of the tool then
        // append a new element
        GXmlElement* country_element = NULL;
        int ncountries = countries->elements();
        if (ncountries == 0) {
            country_element = countries->append(country);
        }
        else {
            for (int k = 0; k < ncountries; ++k) {
                GXmlElement* element = countries->element(k);
                if (element->name() == country) {
                    country_element = element;
                    break;
                }
            }
            if (country_element == NULL) {
                country_element = countries->append(country);
            }
        }

        // Now we have the relevant "country" element and we can update the
        // statistics
        int calls = 1;
        if (country_element->has_attribute("calls")) {
            calls += gammalib::toint(country_element->attribute("calls"));
        }
        if (country_element->has_attribute("wall")) {
            wall += gammalib::todouble(country_element->attribute("wall"));
        }
        if (country_element->has_attribute("cpu")) {
            cpu += gammalib::todouble(country_element->attribute("cpu"));
        }
        if (country_element->has_attribute("gCO2e")) {
            gCO2e += gammalib::todouble(country_element->attribute("gCO2e"));
        }
        country_element->attribute("calls",  gammalib::str(calls));
        country_element->attribute("wall",   gammalib::str(wall));
        country_element->attribute("cpu",    gammalib::str(cpu));
        country_element->attribute("gCO2e",  gammalib::str(gCO2e));

    } // endfor: looped over statistics

    // Return
    return;
}


/***********************************************************************//**
 * @brief Update versions in high-level statistics data.
 *
 * @param[in,out] xml High-level statistics XML object.
 * @param[in] statistics Low-level statistics CSV object.
 *
 * Update version statistics in data block.
 ***************************************************************************/
void GDaemon::update_versions_data(GXml& xml, const GCsv& statistics)
{
    // Get pointer to "versions" node. In case that no "versions" node exists
    // in "data" node then append one.
    GXmlNode* base = xml.element("statistics",0);
    GXmlNode* data = base->element("data",0);
    if (data->elements("versions") == 0) {
        data->append("versions");
    }
    GXmlNode* versions = data->element("versions",0);

    // Loop over statistics
    for (int i = 1; i < statistics.nrows(); ++i) {

        // Extract relevant attributes
        std::string version = statistics.string(i,4);
        double      wall    = statistics.real(i,5);
        double      cpu     = statistics.real(i,6);
        double      gCO2e   = statistics.real(i,7);

        // Get pointer to relevant "version" element. If no "version"
        // element exists that corresponds to the version of the tool then
        // append a new element
        GXmlElement* version_element = NULL;
        int nversion = versions->elements();
        if (nversion == 0) {
            version_element = versions->append(version);
        }
        else {
            for (int k = 0; k < nversion; ++k) {
                GXmlElement* element = versions->element(k);
                if (element->name() == version) {
                    version_element = element;
                    break;
                }
            }
            if (version_element == NULL) {
                version_element = versions->append(version);
            }
        }

        // Now we have the relevant "version" element and we can update the
        // statistics
        int calls = 1;
        if (version_element->has_attribute("calls")) {
            calls += gammalib::toint(version_element->attribute("calls"));
        }
        if (version_element->has_attribute("wall")) {
            wall += gammalib::todouble(version_element->attribute("wall"));
        }
        if (version_element->has_attribute("cpu")) {
            cpu += gammalib::todouble(version_element->attribute("cpu"));
        }
        if (version_element->has_attribute("gCO2e")) {
            gCO2e += gammalib::todouble(version_element->attribute("gCO2e"));
        }
        version_element->attribute("calls",  gammalib::str(calls));
        version_element->attribute("wall",   gammalib::str(wall));
        version_element->attribute("cpu",    gammalib::str(cpu));
        version_element->attribute("gCO2e",  gammalib::str(gCO2e));

    } // endfor: looped over statistics

    // Return
    return;
}


/***********************************************************************//**
 * @brief Update daily statistics
 *
 * @param[in,out] xml High-level statistics XML object.
 * @param[in] statistics Low-level statistics CSV object.
 *
 * Update daily statistics in data block.
 ***************************************************************************/
void GDaemon::update_daily(GXml& xml, const GCsv& statistics)
{
    // Get pointer to "daily" node. In case that no "daily" node exists
    // in "data" node then append one.
    GXmlNode* base = xml.element("statistics",0);
    GXmlNode* data = base->element("data",0);
    if (data->elements("daily") == 0) {
        data->append("daily");
    }
    GXmlNode* daily = data->element("daily",0);

    // Loop over statistics
    for (int i = 1; i < statistics.nrows(); ++i) {

        // Extract relevant attributes
        std::string date    = statistics.string(i,0).substr(0,10);
        std::string country = statistics.string(i,2);
        std::string tool    = statistics.string(i,3);
        std::string version = statistics.string(i,4);
        double      wall    = statistics.real(i,5);
        double      cpu     = statistics.real(i,6);
        double      gCO2e   = statistics.real(i,7);

        // If tool name is empty then set tool name to "unknown"
        if (tool.empty()) {
            tool = "unknown";
        }

        // Get pointer to relevant "date" node. If no "date" node exists
        // that corresponds to the date of the statistics then append an
        // "date" node for the required date.
        GXmlNode* date_node = NULL;
        int ndates = daily->elements("date");
        if (ndates == 0) {
            date_node = daily->append("date");
            date_node->append(GXmlElement("value", date));
        }
        else {
            for (int k = 0; k < ndates; ++k) {
                GXmlElement* element = daily->element("date", k);
                if (element->element("value",0)->value() == date) {
                    date_node = element;
                    break;
                }
            }
            if (date_node == NULL) {
                date_node = daily->append("date");
                date_node->append(GXmlElement("value", date));
            }
        }

        // Get pointer to "tools" node. If no such node exists then append
        // one
        GXmlNode* tools = NULL;
        if (date_node->elements("tools") == 0) {
            tools = date_node->append("tools");
        }
        else {
            tools = date_node->element("tools",0);
        }

        // Get pointer to relevant "tool" element. If no "tool" element
        // exists that corresponds to the name of the tool then append
        // a new element
        GXmlElement* tool_element = NULL;
        int ntools = tools->elements();
        if (ntools == 0) {
            tool_element = tools->append(tool);
        }
        else {
            for (int k = 0; k < ntools; ++k) {
                GXmlElement* element = tools->element(k);
                if (element->name() == tool) {
                    tool_element = element;
                    break;
                }
            }
            if (tool_element == NULL) {
                tool_element = tools->append(tool);
            }
        }

        // Now we have the relevant "tool" element and we can update the
        // statistics
        int calls = 1;
        if (!tool_element->has_attribute("version")) {
            tool_element->attribute("version", version);
        }
        if (!tool_element->has_attribute("country")) {
            tool_element->attribute("country", country);
        }
        if (tool_element->has_attribute("calls")) {
            calls += gammalib::toint(tool_element->attribute("calls"));
        }
        if (tool_element->has_attribute("wall")) {
            wall += gammalib::todouble(tool_element->attribute("wall"));
        }
        if (tool_element->has_attribute("cpu")) {
            cpu += gammalib::todouble(tool_element->attribute("cpu"));
        }
        if (tool_element->has_attribute("gCO2e")) {
            gCO2e += gammalib::todouble(tool_element->attribute("gCO2e"));
        }
        tool_element->attribute("calls",  gammalib::str(calls));
        tool_element->attribute("wall",   gammalib::str(wall));
        tool_element->attribute("cpu",    gammalib::str(cpu));
        tool_element->attribute("gCO2e",  gammalib::str(gCO2e));

    } // endfor: looped over statistics

    // Return
    return;
}
