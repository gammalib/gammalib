/***************************************************************************
 *                     GVOHub.hpp - VO SAMP Hub class                      *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2014 by Thierry Louge                                    *
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
 * @file GVOHub.cpp
 * @brief VO SAMP Hub class implementation
 * @author Thierry Louge
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <cstdlib>         // std::getenv() function
#include <cstdio>          // std::fopen(), etc. functions
#include <cstring>         // std::memset() function
#include <unistd.h>        // close() function
#include <netdb.h>         // getaddrinfo() function
#include <sys/socket.h>    // socket(), connect() functions
#include "GVOHub.hpp"
#include "GVOClient.hpp"
#include "GException.hpp"
#include "GTools.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_HANDLE_REQUEST                  "GVOHub::handle_request(socklen_t)"
#define G_START_HUB                                     "GVOHub::start_hub()"

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
GVOHub::GVOHub(void)
{
    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] hub VO Hub.
 ***************************************************************************/
GVOHub::GVOHub(const GVOHub& hub)
{
    // Initialise members
    init_members();

    // Copy members
    copy_members(hub);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GVOHub::~GVOHub(void)
{
    
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
 * @param[in] hub VO hub.
 * @return VO hub.
 ***************************************************************************/
GVOHub& GVOHub::operator=(const GVOHub& hub)
{
    // Execute only if object is not identical
    if (this != &hub) {

        // Free members
        free_members();

        // Initialise members
        init_members();

        // Copy members
        copy_members(hub);

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
void GVOHub::clear(void)
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
GVOHub* GVOHub::clone(void) const
{
    // Clone client
    return new GVOHub(*this);
}


/***********************************************************************//**
 * @brief Start Hub
 ***************************************************************************/
void GVOHub::start(void)
{
    // Creates Hub
    printf("Creating the hub\n");
    init_hub();
    
    // Starts the hub
    printf("Starting the hub\n");
    start_hub();

    // Return
    return;
}


/***********************************************************************//**
 * @brief 
 *
 *
 ***************************************************************************/
void GVOHub::register_service(void)
{

    // Return
    return;
}


/***********************************************************************//**
 * @brief 
 *
 *
 ***************************************************************************/
void GVOHub::unregister(void)
{

    // Return
    return;
}


/***********************************************************************//**
 * @brief Print VO hub information
 *
 * @param[in] chatter Chattiness (defaults to NORMAL).
 * @return String containing VO hub information
 ***************************************************************************/
std::string GVOHub::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GVOHub ===");

        // Append client information
        result.append("\n"+gammalib::parformat("Name")+m_name);

        // Append Hub information
        result.append("\n"+gammalib::parformat("Hub key")+m_secret);
        result.append("\n"+gammalib::parformat("Hub URL")+m_hub_url);
        result.append("\n"+gammalib::parformat("Hub host (port)"));
        result.append(m_hub_host+" ("+m_hub_port+")");
        result.append("\n"+gammalib::parformat("SAMP protocol version"));
        result.append(m_version);
        result.append("\n"+gammalib::parformat("Hub connection"));
        if (m_socket == -1) {
            result.append("no");
        }
        else {
            if (!m_client_key.empty()) {
                result.append("registered as \""+m_client_id);
                result.append("\" on Hub \""+m_hub_id);
                result.append("\" via socket "+gammalib::str(m_socket));
            }
            else {
                result.append("established on socket ");
                result.append(gammalib::str(m_socket));
            }
        }

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
void GVOHub::init_members(void)
{
    // Initialise members
    m_name     = "GammaLib";
    m_secret   = "mysupersecret#0032557sentence";
    m_hub_url  = "http://localhost:8001";
    m_hub_host = "127.0.0.1";
    m_hub_port = "8001";
    m_version  = "1.3";
    m_client_key.clear();
    m_hub_id   = "b79884e0";
    m_client_id.clear();
    m_socket   = -1;         // Signals no socket

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] hub VO hub.
 ***************************************************************************/
void GVOHub::copy_members(const GVOHub& hub)
{
    // Copy members
    m_name       = hub.m_name;
    m_secret     = hub.m_secret;
    m_hub_url    = hub.m_hub_url;
    m_hub_host   = hub.m_hub_host;
    m_hub_port   = hub.m_hub_port;
    m_version    = hub.m_version;
    m_client_key = hub.m_client_key;
    m_hub_id     = hub.m_hub_id;
    m_client_id  = hub.m_client_id;
    m_socket     = hub.m_socket;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GVOHub::free_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Creates the lockfile, fill it
 *
 *
 * Implements IVOA standard REC-SAMP-1.3-20120411.
 ***************************************************************************/
void GVOHub::init_hub(void)
{
    // Get lockfile URL
    std::string lockurl = get_hub_lockfile();

    // Open SAMP lockfile. Continue only if opening was successful
    FILE* fptr = fopen(lockurl.c_str(), "w");
    if (fptr != NULL) {

        // If successfull writes basic hub configuration to the lockfile
        fputs ("# SAMP lockfile \n# Required keys:\n",fptr);
        char sampsecret[256] = "samp.secret= ";
        strcat(sampsecret,m_secret.c_str());
        strcat(sampsecret,"\n");
        fputs (sampsecret,fptr);
        char sampurl[256] = "samp.hub.xmlrpc.url=";
        strcat(sampurl,m_hub_url.c_str());
        strcat(sampurl,"\n");
        fputs (sampurl,fptr);
        char samprofile[256] = "samp.profile.version=";
        strcat(samprofile,m_version.c_str());
        strcat(samprofile,"\n");
        fputs (samprofile,fptr);
        fputs ("# Info stored by hub for some private reason:\n",fptr);
        char sampid[256] = "gammalib.hubid=";
        strcat(sampid,m_hub_id.c_str());
        strcat(sampid,"\n");
        fputs (sampid,fptr);
    }
    
    // Close SAMP lockfile
    fclose(fptr);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Returns SAMP Hub lockfile URL
 *
 * @return SAMP Hub lockfile URL (empty if no lockfile was found).
 *
 * Implements IVOA standard REC-SAMP-1.3-20120411.
 ***************************************************************************/
std::string GVOHub::get_hub_lockfile(void) const
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


/***********************************************************************//**
 * @brief Reads the client message and runs appropriate function
 *
 * @exception GExpection::invalid_value
 *            Buffer could not be read or written
 ***************************************************************************/
void GVOHub::handle_request(const socklen_t& sock)
{
    // Initialize buffer
    char buffer[256];
    bzero(buffer,256);

    // Read buffer
    int n = read(sock,buffer,255);
    if (n < 0) {
        std::string msg = "Could not read socket buffer.";
        throw GException::invalid_value(G_HANDLE_REQUEST, msg);
    }
    printf("Here is the message: %s\n",buffer);

    // Here, obviously not simple string matching but xml parsing to be done with 
    // help of GXml object. This is only for architectural purpose.
    if (strcmp(buffer,"register") == 0) {
        register_service();
    }
    
    // Same for writing, we'll use GXml
    n = write(sock,"I got your message",18);
    if (n < 0) {
        std::string msg = "Could not write socket buffer.";
        throw GException::invalid_value(G_HANDLE_REQUEST, msg);
    }
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Starts the SAMP hub socket and listens on it
 *
 * @return 1 if error comes out, nothing if everything is right (infinite listening loop)
 *
 ***************************************************************************/
void GVOHub::start_hub(void)
{
    struct sockaddr_in serv_addr,cli_addr;
    int portno,pid;
    socklen_t hub_socket,newsocket,sockfd,clilen;

    // Prepare TCP/IP structure
    bzero((char *) &serv_addr, sizeof(serv_addr));
    portno                    = 8001;
    serv_addr.sin_family      = AF_INET;
    serv_addr.sin_addr.s_addr = INADDR_ANY;
    serv_addr.sin_port        = htons(portno);
    hub_socket                = socket(AF_INET, SOCK_STREAM, 0);
    if (hub_socket < 0) {
        std::string msg = "Socket could not be created.";
        throw GException::invalid_value(G_START_HUB, msg);
    }
    
    // Server socket is opened. Now, bind it to the port, with family etc.
    if (bind(hub_socket, (struct sockaddr *) &serv_addr,sizeof(serv_addr)) < 0) {
        std::string msg = "Could not bind to server socket.";
        throw GException::invalid_value(G_START_HUB, msg);
    }

    // Now start listening for the clients: 5 requests simultaneously maximum
    listen(hub_socket,5);
    clilen = sizeof(cli_addr);
    while (1) {

    	// Accept actual connection from the client 
//std::cout << "Now accept connection the socket" << std::endl;
    	newsocket = accept(hub_socket, (struct sockaddr *)&cli_addr, &clilen);
//std::cout << "Connection accepted" << std::endl;
    	if (newsocket < 0) {
            std::string msg = "Client connection to socket not accepted.";
            throw GException::invalid_value(G_START_HUB, msg);
    	}

    	// Create child process to handle the request
        pid = fork();
//std::cout << "Called fork " << pid << std::endl;
        if (pid < 0) {
            std::string msg = "Not child process created.";
            throw GException::invalid_value(G_START_HUB, msg);
        }
        if (pid == 0) {
            // Child process: client process
            handle_request(newsocket);
            exit(0);
        }
        else {
//std::cout << "Now close the socket" << std::endl;
            close(newsocket);
//std::cout << "Socket closed" << std::endl;
        }

    }

    // Return
    return;
}
