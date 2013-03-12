/***************************************************************************
 *                     GVOClient.hpp - VO client class                     *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2013 by Juergen Knoedlseder                              *
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
 * @file GVOClient.cpp
 * @brief VO client class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <cstdlib>         // std::getenv() function
#include "GVOClient.hpp"
#include "GException.hpp"
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
GVOClient::GVOClient(void)
{
    // Initialise members
    init_members();

    // Find Hub
    find_hub();
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] client VO client.
 ***************************************************************************/
GVOClient::GVOClient(const GVOClient& client)
{
    // Initialise members
    init_members();

    // Copy members
    copy_members(client);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GVOClient::~GVOClient(void)
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
 * @param[in] client VO client.
 * @return VO client.
 ***************************************************************************/
GVOClient& GVOClient::operator= (const GVOClient& client)
{
    // Execute only if object is not identical
    if (this != &client) {

        // Free members
        free_members();

        // Initialise members
        init_members();

        // Copy members
        copy_members(client);

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
 * @brief Clear object.
 *
 * Reset object to a clean initial state.
 ***************************************************************************/
void GVOClient::clear(void)
{
    // Free memory and initialise members
    free_members();
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone object
 ***************************************************************************/
GVOClient* GVOClient::clone(void) const
{
    // Clone client
    return new GVOClient(*this);
}


/***********************************************************************//**
 * @brief Register client to SAMP Hub
 *
 * Connects the VO client to the Hub. A socket to the Hub is opened and the
 * client is registered at the Hub.
 *
 * @todo Implement method
 ***************************************************************************/
void GVOClient::connect(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Unregister client from SAMP Hub
 *
 * Disconnects the VO client from the Hub. The client is unregistered and the
 * socket to the Hub is closed.
 *
 * @todo Implement method
 ***************************************************************************/
void GVOClient::disconnect(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Signals if client has Hub information
 *
 * @return True if client has Hub information, false otherwise.
 *
 * Checks if all mandatory Hub tokens are non-empty.
 ***************************************************************************/
bool GVOClient::hashub(void) const
{
    // Return Hub information status
    return (!m_secret.empty() && !m_hub_url.empty() && !m_version.empty());
}


/***********************************************************************//**
 * @brief Print VO client information
 *
 * @return String containing VO client information
 ***************************************************************************/
std::string GVOClient::print(void) const
{
    // Initialise result string
    std::string result;

    // Append header
    result.append("=== GVOClient ===");

    // Append client information
    result.append("\n"+parformat("Hub key")+m_secret);
    result.append("\n"+parformat("Hub URL")+m_hub_url);
    result.append("\n"+parformat("SAMP protocal version")+m_version);
    result.append("\n"+parformat("Client key")+m_client_key);

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
void GVOClient::init_members(void)
{
    // Initialise members
    m_secret.clear();
    m_hub_url.clear();
    m_version.clear();
    m_client_key.clear();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] client VO client.
 ***************************************************************************/
void GVOClient::copy_members(const GVOClient& client)
{
    // Copy members
    m_secret     = client.m_secret;
    m_hub_url    = client.m_hub_url;
    m_version    = client.m_version;
    m_client_key = client.m_client_key;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GVOClient::free_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Find SAMP Hub
 *
 * @return True of SAMP Hub has been found, false otherwise.
 *
 * Search a valid SAMP Hub and retrieve all mandatory token for this Hub.
 * The manadtory tokens are
 *
 *     samp.secret           Opaque text string required for Hub registration
 *     samp.hub.xmlrpc.url   XML-RPC endpoint for communication
 *     samp.profile.version  Version of SAMP profile
 *
 * Implements IVOA standard REC-SAMP-1.3-20120411.
 *
 * @todo If no valid Hub was found the method should attempt to start its
 * own Hub.
 ***************************************************************************/
bool GVOClient::find_hub(void)
{
    // Allocate line buffer
    const int n = 1000; 
    char      line[n];

    // Initialise find flag to false
    bool found = false;

    // Get lockfile URL
    std::string lockurl = hub_lockfile();

    // Continue only if a URL has been found
    if (!lockurl.empty()) {

        // If we have a file:// prefix then strip it now. This is a kluge
        // and should be remplaced by a method that allows opening any kind
        // of URL
        if (lockurl.compare(0, 7, "file://") == 0) {
            lockurl = lockurl.substr(7, std::string::npos);
        }
    
        // Open SAMP lockfile. Continue only if opening was successful
        FILE* fptr = fopen(lockurl.c_str(), "r");
        if (fptr != NULL) {

            // Parse lockfile and search for mandatory tokens
            while (fgets(line, n, fptr) != NULL) {

                // Convert line to C++ string
                std::string cline = std::string(line);

                // Check for secret key
                if (cline.compare(0, 12, "samp.secret=") == 0) {
                    m_secret = strip_chars(cline.substr(12, std::string::npos), "\r\n");
                }

                // Check for Hub URL
                else if (cline.compare(0, 20, "samp.hub.xmlrpc.url=") == 0) {
                    m_hub_url = strip_chars(cline.substr(20, std::string::npos), "\r\n");
                }

                // Check for profile version
                else if (cline.compare(0, 21, "samp.profile.version=") == 0) {
                    m_version = strip_chars(cline.substr(21, std::string::npos), "\r\n");
                }

            }

            // Close SAMP lockfile
            fclose(fptr);

            // Check for existence of mandatory tokens
            found = hashub();
        }

    } // endif: URL has been found

    //TODO: Create GammaLib own Hub in case that no Hub has been found
    //if (!found) {
    //}

    // If no Hub has been found, clear all Hub related members
    if (!found) {
        m_secret.clear();
        m_hub_url.clear();
        m_version.clear();
    }

    // Return find flag
    return found;
}


/***********************************************************************//**
 * @brief Returns SAMP Hub lockfile URL
 *
 * @return SAMP Hub lockfile URL (empty if no lockfile was found).
 *
 * Implements IVOA standard REC-SAMP-1.3-20120411.
 ***************************************************************************/
std::string GVOClient::hub_lockfile(void) const
{
    // Initialise result
    std::string url = "";

    // Check for existence of the SAMP_HUB environment variable first
    char* hub_ptr = std::getenv("SAMP_HUB");
    if (hub_ptr != NULL) {

        // Check for mandatory std-lockurl: prefix (no other prefixe is
        // supported so far)
        std::string lockurl = std::string(hub_ptr);
        if (lockurl.compare(0, 12, "std-lockurl:") == 0) {

            // Extract URL
            url = lockurl.substr(12, std::string::npos);

        } // endif: std-lockurl: prefix found
  
    }

    // ... otherwise the lockfile should be $HOME/.samp
    else {

        // Get user's HOME directory path as the prefix of the full
        // path. If the HOME environment variable is not set we
        // expect that .samp is in the local directory. This is non
        // standard, but prevents for creating an exception here.
        std::string prefix = "";
        char* home_ptr = std::getenv("HOME");
        if (home_ptr != NULL) {
            prefix = std::string(home_ptr) + "/";
        }

        // Set filename
        url = prefix + ".samp";

    } // endelse: no SAMP_HUB environment variable found

    // Return URL
    return url;
}
