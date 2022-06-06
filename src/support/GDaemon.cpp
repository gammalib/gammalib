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
#include <signal.h>        // sigaction
#include <string.h>        // memset
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
#define G_APPEND_LOG //!< Append output to any existing log file

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
    // Get process ID
    m_pid = getpid();

    // Open logger
    #if defined(G_APPEND_LOG)
    m_log.open(gammalib::gamma_filename("daemon.log"), false);
    #else
    m_log.open(gammalib::gamma_filename("daemon.log"), true);
    #endif
    m_log.chatter((GChatter)m_chatter);
    m_log.date(true);

    // Log start up of daemon
    m_log << "[" << (int)(m_pid) << "] ";
    m_log << "GammaLib daemon start up" << std::endl;

    // Set daemon name
    #if defined(__linux__) || defined(__linux) || defined(linux) || defined(__gnu_linux__)
    prctl(PR_SET_NAME, "gammalibd");
    #endif

    // Create lock file
    create_lock_file();

    // Create heartbeat file
    write_heartbeat();

    // Main event handling loop
    while (true) {

        // If another daemon is running then exit this daemon
        pid_t pid = lock_pid();
        if ((pid != 0) && (pid != m_pid)) {
            m_log << "[" << (int)(m_pid) << "] ";
            m_log << "Another daemon with PID " << (int)(pid) << " is running, ";
            m_log << "exit this daemon" << std::endl;
            break;
        }

        // Put all activities into try-catch block
        try {

            // Update host-country
            update_host_country();

            // Update application statistics
            update_statistics();

        }
        catch (const std::exception &except) {
            m_log << "[" << (int)(m_pid) << "] ";
            m_log << "Exception catched by daemon" << std::endl;
            m_log << except.what() << std::endl;
        }

        // Force logger flushing
        m_log.flush(true);

        // Determine number of heartbeat cycles before wakeup
        int ncycle = int(double(m_period) / double(m_heartbeat) + 0.5);
        if (ncycle < 1) {
            ncycle = 1;
        }

        // Loop over heartbeat cycles and wait for heartbeats
        for (int i = 0; i < ncycle; ++i) {
            write_heartbeat();
            sleep(m_heartbeat);
        }

    } // endwhile: main event loop

    // Log termination of event loop
    m_log << "[" << (int)(m_pid) << "] ";
    m_log << "Terminated event loop" << std::endl;

    // Delete lock file
    delete_lock_file();

    // Force logger flushing
    m_log.flush(true);

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
 * Checks if the daemon is alive by comparing the time found in the heartbeat
 * file to the current time. If the heartbeat file is not older than twice
 * the heartbeat period the daemon is considered alive.
 *
 * This method can be used by a client to check whether a daemon is alive.
 * The method does not rely on any process ID and works even over multiple
 * machines that share the same disk space.
 ***************************************************************************/
bool GDaemon::alive(void) const
{
    // Initialise flag
    bool alive = false;

    // Set current time
    GTime now;
    now.now();

    // Get heartbeat filename
    GFilename filename = heartbeat_filename();

    // OpenMP critical zone to listen for heartbeats
    #pragma omp critical(GDaemon_alive)
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

            // Open heartbeat file. Continue only if opening was successful.
            FILE* fptr = fopen(filename.url().c_str(), "r");
            if (fptr != NULL) {

                // Allocate line buffer
                const int n = 1000;
                char      line[n];

                // Read line
                fgets(line, n-1, fptr);

                // Close file
                fclose(fptr);

                // Extract time. Put this in a try-catch block to catch
                // any ill-conditioned time string
                try {
                    GTime heartbeat;
                    heartbeat.utc(std::string(line));
                    if ((now - heartbeat) < 2.0 * double(m_heartbeat)) {
                        alive = true;
                    }
                }
                catch (const std::exception &except) {
                    ;
                }

            } // endif: lock file opened successfully

            // Unlock file
            lock.l_type = F_UNLCK;
            fcntl(fd, F_SETLK, &lock);

            // Close file
            close(fd);

        } // endif: file locking successful

    } // end of OMP critial zone

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
        result.append("\n"+gammalib::parformat("Heartbeat period (s)"));
        result.append(gammalib::str(m_heartbeat));
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
    m_pid       = 0;     //!< No process ID
    m_period    = 3600;  //!< Wake-up daemon every hour
    m_heartbeat = 60;    //!< One heartbeat per minute
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
    m_pid       = daemon.m_pid;
    m_period    = daemon.m_period;
    m_heartbeat = daemon.m_heartbeat;
    m_log       = daemon.m_log;
    m_chatter   = daemon.m_chatter;

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
    if (fptr != NULL) {

        // Write process ID into lockfile
        fprintf(fptr, "%d\n", (int)(m_pid));

        // Close lockfile
        fclose(fptr);

        // Log creation of lock file
        m_log << "[" << (int)(m_pid) << "] ";
        m_log << "Created lock file \"" << lockfile.url();
        m_log << "\"" << std::endl;

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
void GDaemon::delete_lock_file(void)
{
    // If process ID in lock file is the ID of the current process then
    // delete the lock file
    if (lock_pid() == getpid()) {

        // Get lock filename
        GFilename lockfile = lock_filename();

        // Delete lock file
        std::remove(lockfile.url().c_str());

        // Log removal of lock file
        m_log << "[" << (int)(m_pid) << "] ";
        m_log << "Removed lock file \"" << lockfile.url();
        m_log << "\"" << std::endl;

    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write heartbeat file
 *
 * Creates and updates the heartbeat file that contains the current UTC time.
 ***************************************************************************/
void GDaemon::write_heartbeat(void)
{
    // Set current time
    GTime now;
    now.now();

    // Get heartbeat filename
    GFilename filename = heartbeat_filename();

    // Open heartbeat file. Continue only if opening was successful.
    FILE* fptr = fopen(filename.url().c_str(), "w");
    if (fptr != NULL) {

        // Write current time into heartbeat file
        fprintf(fptr, "%s\n", now.utc().c_str());

        // Close lockfile
        fclose(fptr);

    } // endif: Heartbeat file opened successfully

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

        // OpenMP critical zone to read statitics
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

                // Log updating of high-level statistics
                m_log << "[" << (int)(m_pid) << "] ";
                m_log << "Update high-level statistics" << std::endl;

                // Lock file
                fcntl(fd, F_SETLKW, &lock);

                // Put all activities into try-catch block
                try {

                    // Load CSV file
                    GCsv statistics(filename, ",");

                    // Recover valid "statistics.xml" file
                    recover_valid_xml();

                    // Set high-level statistics filename
                    GFilename filename_work = gammalib::gamma_filename("statistics.xml");
                    GFilename filename_copy = gammalib::gamma_filename("statistics.xml~");

                    // Load high-level statistics file
                    GXml xml;
                    xml.load(filename_work);

                    // Update file only if XML document is not empty
                    if (!xml.is_empty()) {

                        // Save XML file into a copy
                        xml.save(filename_copy);

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
                        xml.save(filename_work);

                        // Now get rid of the copy
                        std::remove(filename_copy.url().c_str());

                        // And finally get rid of the low-level statistics file
                        std::remove(filename.url().c_str());

                    } // endif: XML document was not empty

                }
                catch (const std::exception &except) {
                    m_log << "[" << (int)(m_pid) << "] ";
                    m_log << "Exception catched by daemon" << std::endl;
                    m_log << except.what() << std::endl;
                }

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
 * @brief Update host country
 *
 * Update host country code in $HOME/.gamma/host-country file.
 ***************************************************************************/
void GDaemon::update_host_country(void)
{
    // Get host country code by forcing a query
    std::string country = gammalib::host_country(true);

    // Continue only if host country is a two-digit code
    if (country.length() == 2) {

        // Set host country filename
        GFilename filename = gammalib::gamma_filename("host-country");

        // OpenMP critical zone to write host country code
        #pragma omp critical(GDaemon_update_host_country)
        {

            // Open host-country file, and in case of success, write
            // country
            FILE* fptr = fopen(filename.url().c_str(), "w");
            if (fptr != NULL) {
                fprintf(fptr, "%s\n", country.c_str());
                fclose(fptr);
            }

        } // end of OMP critial zone

    } // endif: we had a two-digit country code

    // Return
    return;
}


/***********************************************************************//**
 * @brief Recovers a valid XML file
 *
 * This method recovers a valid XML file in @p statistics.xml. Several cases
 * are covered.
 *
 * If none of the files @p statistics.xml and @p statistics.xml~ exists
 * the method will create a new XML file using the create_xml() method.
 *
 * If only the copy @p statistics.xml~ exists there was a problem during
 * writing the working file, hence the @p statistics.xml~ file is copied
 * into @p statistics.xml. In case that @p statistics.xml~ is corrupt a
 * new XML file is created.
 *
 * If both files exist, the integrity of both files is checked. If only
 * @p statistics.xml is corrupted, @p statistics.xml~ will be copied
 * into @p statistics.xml. If only @p statistics.xml~ is corrupted it
 * is ignored. If both files are corrupted they are secured and a new XML
 * file will be created using the create_xml() method. If both files are okay,
 * the file @p statistics.xml is secured and @p statistics.xml~ is copied
 * into @p statistics.xml.
 *
 * Before existing, any @p statistics.xml.copy file is removed.
 ***************************************************************************/
void GDaemon::recover_valid_xml(void)
{
    // Set filenames
    GFilename filename_work = gammalib::gamma_filename("statistics.xml");
    GFilename filename_copy = gammalib::gamma_filename("statistics.xml~");

    // If none of the files exist then create a new working file
    if (!filename_work.exists() && !filename_copy.exists()) {
        create_xml(filename_work);
    }

    // ... else if only the copy exists and if the copy is okay, then use it.
    // Otherwise create a new XML file
    else if (!filename_work.exists() && filename_copy.exists()) {
        GXml xml_copy;
        try {
            xml_copy.load(filename_copy);
            xml_copy.save(filename_work);
            m_log << "[" << (int)(m_pid) << "] ";
            m_log << "No \"statistics.xml\" file found, use copy ";
            m_log << "\"statistics.xml~\"" << std::endl;
        }
        catch (const std::exception &except) {
            m_log << "[" << (int)(m_pid) << "] ";
            m_log << "No \"statistics.xml\" file found and corrupt ";
            m_log << "\"statistics.xml~\" file encountered, secure ";
            m_log << "it and create new XML file" << std::endl;
            GTime now;
            now.now();
            std::string newcopy = filename_work.url() + "." + now.utc() + "~";
            std::rename(filename_copy.url().c_str(), newcopy.c_str());
            create_xml(filename_work);
        }
    }

    // ... else if both of the files exist then an issue occured in writing
    // the working file (otherwise the copy would have been removed)
    else if (filename_work.exists() && filename_copy.exists()) {

        // Check integrity of working file and copy
        bool integrity_work = true;
        bool integrity_copy = true;
        GXml xml_work;
        GXml xml_copy;
        try {
            xml_work.load(filename_work);
            integrity_work = true;
        }
        catch (const std::exception &except) {
            integrity_work = false;
        }
        try {
            xml_copy.load(filename_copy);
            integrity_copy = true;
        }
        catch (const std::exception &except) {
            integrity_copy = false;
        }

        // If only the copy is okay but the working file is corrupt the
        // situation is clear: the working file got corrupted. In that
        // case save copy as working file.
        if (integrity_copy && !integrity_work) {
            m_log << "[" << (int)(m_pid) << "] ";
            m_log << "Corrupt \"statistics.xml\" file encountered, ";
            m_log << "use copy \"statistics.xml~\"" << std::endl;
            xml_copy.save(filename_work);
        }

        // ... otherwise if both files are corrupted then rename the
        // working file and create a new fresh high-level statistics.
        // This case should actually never happen!
        else if (!integrity_copy && !integrity_work) {
            m_log << "[" << (int)(m_pid) << "] ";
            m_log << "Corrupt \"statistics.xml\" and \"statistics.xml~\" ";
            m_log << "files encountered, secure them and create new ";
            m_log << "file" << std::endl;
            GTime now;
            now.now();
            std::string newwork = filename_work.url() + "." + now.utc();
            std::string newcopy = filename_work.url() + "." + now.utc() + "~";
            std::rename(filename_work.url().c_str(), newwork.c_str());
            std::rename(filename_copy.url().c_str(), newcopy.c_str());
            create_xml(filename_work);
        }

        // ... otherwise if both files are okay, then the working file is
        // secured and the copy is used, as in normal workings it is not
        // expected that a copy exists, so some issue happend during the
        // last update. This case should actually never happen!
        else if (integrity_copy && integrity_work) {
            m_log << "[" << (int)(m_pid) << "] ";
            m_log << "Unexpected \"statistics.xml\" file encountered, ";
            m_log << "secure file and use copy \"statistics.xml~\"";
            m_log << std::endl;
            GTime now;
            now.now();
            std::string newwork = filename_work.url() + "." + now.utc();
            std::rename(filename_work.url().c_str(), newwork.c_str());
            xml_copy.save(filename_work);
        }

    } // endelse: both files existed

    // Make sure we get rid of any copy
    std::remove(filename_copy.url().c_str());

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
    // Set current time
    GTime now;
    now.now();

    // Initialise high-level statistics XML object
    GXml xml;

    // Append base node
    GXmlNode* base = xml.append(GXmlElement("statistics title=\"High-level statistics\""));

    // Append header
    GXmlNode* header    = base->append("header");
    GXmlNode* dates     = header->append("dates");
    GXmlNode* countries = header->append("countries");
    GXmlNode* d1        = dates->append("creation");
    d1->append(GXmlText(GXmlText(now.utc())));
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
    m_log << "[" << (int)(m_pid) << "] ";
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
    if (modified->is_empty()) {
        modified->append(GXmlText(now.utc()));
    }
    else {
        static_cast<GXmlText*>((*modified)[0])->text(now.utc());
    }

    // Update first and last time in header
    std::string start;
    std::string stop;
    if (statistics.nrows() > 1) {
        start                = statistics.string(1,0);
        stop                 = statistics.string(statistics.nrows()-1,0);
        GXmlNode* start_node = header->element("dates > start");
        GXmlNode* stop_node  = header->element("dates > stop");
        if (start_node->is_empty()) {
            start_node->append(GXmlText(start));
        }
        else {
            std::string current_start = static_cast<GXmlText*>((*start_node)[0])->text();
            if (current_start.empty() || (start < current_start)) {
                static_cast<GXmlText*>((*start_node)[0])->text(start);
            }
        }
        if (stop_node->is_empty()) {
            stop_node->append(GXmlText(stop));
        }
        else {
            std::string current_stop  = static_cast<GXmlText*>((*stop_node)[0])->text();
            if (current_stop.empty() || (stop > current_stop)) {
                static_cast<GXmlText*>((*stop_node)[0])->text(stop);
            }
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

        // If country is empty then set country to "unknown"
        if (country.empty()) {
            country = "unknown";
        }

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

        // If country is empty then set country to "unknown"
        if (country.empty()) {
            country = "unknown";
        }

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

        // If version is empty then set version to "unknown"
        if (version.empty()) {
            version = "unknown";
        }

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
