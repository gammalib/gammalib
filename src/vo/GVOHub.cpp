/***************************************************************************
 *                     GVOHub.cpp - VO SAMP Hub class                      *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2014-2020 by Thierry Louge                               *
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
#include <csignal>         // signal() function
#include <cerrno>          // errno
#include <unistd.h>        // close() function
#include <netdb.h>         // getaddrinfo() function
#include <netinet/in.h>    // sockaddr_in, INADDR_ANY, htons
#include <fstream>
#include <sys/shm.h>
#include <sys/socket.h>    // socket(), connect() functions
#include <sys/wait.h>      // waitpid() function
#include <arpa/inet.h>     // inet_addr() function
#include "GTools.hpp"
#include "GException.hpp"
#include "GXml.hpp"
#include "GXmlNode.hpp"
#include "GVOHub.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_START_HUB                                     "GVOHub::start_hub()"
#define G_GET_SOCKET                                   "GVOHub::get_socket()"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */
//#define G_CONSOLE_DUMP               //!< Show method headers
//#define G_CONSOLE_ERRORS             //!< Show error messages
//#define G_SHOW_MESSAGE               //!< Show posted and received messages


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
    // Free members
    free_members();

    // Initialise members
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
    // Start Hub
    start_hub();

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

        // Append Hub information
        result.append("\n"+gammalib::parformat("Hub identifier")+m_hub_id);
        result.append("\n"+gammalib::parformat("Hub key")+m_secret);
        result.append("\n"+gammalib::parformat("Hub URL")+hub_url());
        result.append("\n"+gammalib::parformat("Hub host")+m_hub_host);
        result.append("\n"+gammalib::parformat("Hub port")+m_hub_port);
        result.append("\n"+gammalib::parformat("Hub path")+m_hub_path);
        result.append("\n"+gammalib::parformat("SAMP protocol version"));
        result.append(m_version);

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
    m_secret        = random_string(15);
    m_hub_host      = "127.0.0.1";
    m_hub_port      = "2526";
    m_hub_path      = "xmlrpc";
    m_version       = "1.3";
    m_hub_id        = "gammalib_hub";
    m_socket        = -1;        // Signals no socket
    m_shutdown      = false;
    m_clients.clear();

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
    m_secret   = hub.m_secret;
    m_hub_host = hub.m_hub_host;
    m_hub_port = hub.m_hub_port;
    m_hub_path = hub.m_hub_path;
    m_version  = hub.m_version;
    m_hub_id   = hub.m_hub_id;
    m_socket   = hub.m_socket;
    m_shutdown = hub.m_shutdown;
    m_clients  = hub.m_clients;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GVOHub::free_members(void)
{ 
    // Close sockets
    if (m_socket != -1) {
        close(m_socket);
        m_socket = -1;
    }

    // Remove lockfile
    delete_samp_file();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Starts the SAMP hub socket and listens on it
 *
 * @exception GException::runtime_error
 *            Problem with creating of or listening on socket
 *
 * This is the main Hub event loop.
 ***************************************************************************/
void GVOHub::start_hub(void)
{
    // Get socket
    m_socket = get_socket();

    // Create SAMP file
    create_samp_file();

    // Start listening for the clients: 5 requests simultaneously pending
    // at maximum
    if (listen(m_socket, 5) < 0) {
        std::string msg = "Unable to start listening on Hub socket. Errno="+
                          gammalib::str(errno);
        throw GException::runtime_error(G_START_HUB, msg);
    }

    // ...
    struct sockaddr_in cli_addr;
    socklen_t clilen = sizeof(cli_addr);

    // Main event handling loop
    while (true) {

        // Accept connection from the client 
    	int socket = accept(m_socket, (struct sockaddr *)&cli_addr, &clilen);
    	if (socket < 0) {
            std::string msg = "Client connection to socket not accepted.";
            throw GException::runtime_error(G_START_HUB, msg);
    	}

        // Handle request
        handle_request(socket);

        // Close socket
        close(socket);

        // If we have received a shutdown request then exit the event
        // loop
        if (m_shutdown) {
            break;
        }

    } // endwhile: main event loop

    // Close socket
    close(m_socket);
    m_socket = -1;

    // Delete SAMP file
    delete_samp_file();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Reads the client message and runs appropriate function
 *
 * @param[in] sock Socket to client.
 ***************************************************************************/
void GVOHub::handle_request(const socklen_t& sock)
{
    // Header
    #if defined(G_CONSOLE_DUMP)
    std::cout << "GVOHub::handle_request" << std::endl;
    #endif

    // Initialize buffer
    char buffer[1001];

    // Initialise message
    std::string message = "";

    // Read from socket until nothing is received anymore.
    int timeout = 2000; // Initial timeout is 2 sec
    int n       = 0;
    do {
        n = gammalib::recv(sock, buffer, 1000, 0, timeout);
        if (n > 0) {
            buffer[n] = '\0';
            message.append(std::string(buffer));
        }
        timeout = 10; // The timeout now is 0.01 sec 
    } while (n > 0);

    // Dump the buffer
    #if defined(G_SHOW_MESSAGE)
    std::cout << std::endl;
    std::cout << "GVOHub has received the following message:" << std::endl;
    std::cout << "==========================================" << std::endl;
    std::cout << message << std::endl;
    #endif

    // Extract response into an XML object
    GXml xml;
    size_t start = message.find("<?xml");
    if (start != std::string::npos) {
        xml = GXml(message.substr(start, std::string::npos));
    }

    // Get methodName value
    std::string method_called;
    const GXmlNode* node = xml.element("methodCall > methodName");
    if (node != NULL) {
        const GXmlText* text = static_cast<const GXmlText*>((*node)[0]);
        if (text != NULL) {
            method_called = text->text();
        }
    }

    // Dump the method that is called
    #if defined(G_CONSOLE_DUMP)
    std::cout << "GVOHub has received the message type \"";
    std::cout << method_called << "\"" << std::endl;
    #endif

    // Dispatch according to method
    if (method_called.compare("samp.hub.ping") == 0) {
        request_ping(sock);
    }
    else if (method_called.compare("samp.hub.register") == 0) {
        request_register(xml, sock);
    }
    else if (method_called.compare("samp.hub.unregister") == 0) {
        request_unregister(xml, sock);
    }
    else if (method_called.compare("samp.hub.declareMetadata") == 0) {
        request_declare_metadata(xml, sock);
    }
    else if (method_called.compare("samp.hub.declareSubscriptions") == 0) {
        request_declare_subscriptions(xml, sock);
    }
    else if (method_called.compare("samp.hub.setXmlrpcCallback") == 0) {
        request_set_xml_rpc_callback(xml, sock);
    }
    else if (method_called.compare("samp.hub.getSubscriptions") == 0) {
        request_get_subscriptions(xml, sock);
    }
    else if (method_called.compare("samp.hub.getRegisteredClients") == 0) {
        request_get_registered_clients(xml, sock);
    }
    else if (method_called.compare("samp.hub.getSubscribedClients") == 0) {
        request_get_subscribed_clients(xml, sock);
    }
    else if (method_called.compare("samp.hub.getMetadata") == 0) {
        request_get_metadata(xml, sock);
    }
    else if (method_called.compare("samp.hub.notify") == 0) {
        request_ping(sock);
    }
    else if (method_called.compare("samp.hub.notifyAll") == 0) {
        request_notify_all(xml, sock);
    }
    else if (method_called.compare("samp.hub.call") == 0) {
        request_ping(sock);
    }
    else if (method_called.compare("samp.hub.callAll") == 0) {
        request_ping(sock);
    }
    else if (method_called.compare("samp.hub.callAndWait") == 0) {
        request_ping(sock);
    }
    else if (method_called.compare("samp.hub.reply") == 0) {
        request_ping(sock);
    }
    else if (method_called.compare("samp.hub.shutdown") == 0) {
        request_shutdown(sock);
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Handles ping requests
 *
 * @param[in] sock Socket.
 *
 * Handles all incoming ping requests.
 ***************************************************************************/
void GVOHub::request_ping(const socklen_t& sock)
{
    // Header
    #if defined(G_CONSOLE_DUMP)
    std::cout << "GVOHub::request_ping" << std::endl;
    #endif

    // Post void message
    post_samp_void(sock);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Handles registration requests
 *
 * @param[in] xml XML message sent by client.
 * @param[in] sock Socket.
 *
 * Handles all incoming registration requests.
 ***************************************************************************/
void GVOHub::request_register(const GXml& xml, const socklen_t& sock)
{
    // Header
    #if defined(G_CONSOLE_DUMP)
    std::cout << "GVOHub::request_register" << std::endl;
    #endif

    // Set client reference
    int         counter = 0;
    std::string reference;
    do {
        reference = "c"+gammalib::str(counter);
        int i = 0;
        for (; i < m_clients.size(); ++i) {
            if (m_clients[i].reference == reference) {
                counter++;
                break;
            }
        }
        if (i >= m_clients.size()) {
            break;
        }
    } while (true);

    // Create a new client
    struct client voclient;
    voclient.reference   = reference;
    voclient.private_key = random_string(15);

    // Attached client
    m_clients.push_back(voclient);

    // Declare response
    std::string response = "";

    // Set response
    response.append("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
    response.append("<methodResponse>\n");
    response.append("<params>\n");
    response.append("  <param><value><struct>\n");
    response.append("    <member>\n");
    response.append("      <name>samp.private-key</name>\n");
    response.append("      <value>"+voclient.private_key+"</value>\n");
    response.append("    </member>\n");
    response.append("    <member>\n");
    response.append("      <name>samp.hub-id</name>\n");
    response.append("      <value>"+m_hub_id+"</value>\n");
    response.append("    </member>\n");
    response.append("    <member>\n");
    response.append("      <name>samp.self-id</name>\n");
    response.append("      <value>"+voclient.reference+"</value>\n");
    response.append("    </member>\n");
    response.append("    <member>\n");
    response.append("      <name>samp.status</name>\n");
    response.append("      <value>samp.ok</value>\n");
    response.append("    </member>\n");
    response.append("    <member>\n");
    response.append("      <name>samp.url-translator</name>\n");
    response.append("      <value>"+hub_url()+"</value>\n");
    response.append("    </member>\n");
    response.append("  </struct></value></param>\n");
    response.append("</params>\n");
    response.append("</methodResponse>\n");

    // Post response
    post_string(response, sock);

    // Notify all clients about the registering
    for (int i = 0; i < m_clients.size(); ++i) {
        for (int j = 0; j < m_clients[i].subscriptions.size(); ++j) {
            if ("samp.hub.event.register" == m_clients[i].subscriptions[j]) {
                notify_register(m_clients[i], voclient.reference);
                break;
            }
        }
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Handles unregistration requests
 *
 * @param[in] xml XML message sent by client.
 * @param[in] sock Socket.
 *
 * Handles all incoming unregistration requests.
 ***************************************************************************/
void GVOHub::request_unregister(const GXml& xml, const socklen_t& sock)
{
    // Optionally dump header
    #if defined(G_CONSOLE_DUMP)
    std::cout << "GVOHub::request_unregister" << std::endl;
    #endif 

    // Get client index
    int i = get_client_index(xml);

    // Continue only if index is valid
    if (i != -1) {

        // Post void message
        post_samp_void(sock);

        // Notify all clients about the registering
        for (int k = 0; k < m_clients.size(); ++k) {
            if (i != k) {
                for (int j = 0; j < m_clients[k].subscriptions.size(); ++j) {
                    if ("samp.hub.event.unregister" == m_clients[k].subscriptions[j]) {
                        notify_unregister(m_clients[k], m_clients[i].reference);
                    }
                }
            }
        }

        // Remove client
        m_clients.erase(m_clients.begin()+i);

    } // endif: valid client found

    // Signal if the client was unknown
    #if defined(G_CONSOLE_ERRORS)
    else {
        std::cout << "*** ERROR: GVOHub::request_unregister: ";
        std::cout << "Unable to find client \"" << get_client_key(xml);
        std::cout << "\"" << std::endl;
    }
    #endif

    // Return
    return;
}


/***********************************************************************//**
 * @brief Handles metadata declaration requests
 *
 * @param[in] xml XML message sent by client.
 * @param[in] sock Socket.
 *
 * Handles all incoming metadata declaration requests.
 ***************************************************************************/
void GVOHub::request_declare_metadata(const GXml& xml, const socklen_t& sock)
{
    // Header
    #if defined(G_CONSOLE_DUMP)
    std::cout << "GVOHub::request_declare_metadata" << std::endl;
    #endif

    // Get client index
    int i = get_client_index(xml);

    // Continue only if index is valid
    if (i != -1) {

        // Get node
        const GXmlNode* node = xml.element("methodCall > params > param[1] > value > struct");

        // Store metadata
        m_clients[i].name          = get_response_value(node, "samp.name");;
        m_clients[i].description   = get_response_value(node, "samp.description.text");
        m_clients[i].icon          = get_response_value(node, "samp.icon.url");
        m_clients[i].documentation = get_response_value(node, "samp.documentation.url");
        m_clients[i].affiliation   = get_response_value(node, "author.affiliation");
        m_clients[i].email         = get_response_value(node, "author.email");
        m_clients[i].author_name   = get_response_value(node, "author.name");
        m_clients[i].homepage      = get_response_value(node, "home.page");

        // Post void message
        post_samp_void(sock);

        // Notify all clients about the registering
        for (int k = 0; k < m_clients.size(); ++k) {
            if (i != k) {
                for (int j = 0; j < m_clients[k].subscriptions.size(); ++j) {
                    if ("samp.hub.event.metadata" == m_clients[k].subscriptions[j]) {
                        notify_metadata(m_clients[k], m_clients[i].reference);
                    }
                }
            }
        }

    } // endif: client found

    // Signal if the client was unknown
    #if defined(G_CONSOLE_ERRORS)
    else {
        std::cout << "*** ERROR: GVOHub::request_declare_metadata: ";
        std::cout << "Unable to find client \"" << get_client_key(xml);
        std::cout << "\"" << std::endl;
    }
    #endif

    // Return
    return;
}


/***********************************************************************//**
 * @brief Handles subscriptions declaration requests
 *
 * @param[in] xml XML message sent by client.
 * @param[in] sock Socket.
 *
 * Handles all incoming subscription declaration requests.
 ***************************************************************************/
void GVOHub::request_declare_subscriptions(const GXml& xml,const socklen_t& sock)
{
    // Header
    #if defined(G_CONSOLE_DUMP)
    std::cout << "GVOHub::request_declare_subscriptions" << std::endl;
    #endif

    // Get client index
    int i = get_client_index(xml);

    // Continue only if index is valid
    if (i != -1) {

        // Get subscriptions
        std::vector<std::string> subscriptions = get_subscriptions(xml);

        // Store subscriptions
        for (int k = 0; k < subscriptions.size(); ++k) {
            m_clients[i].subscriptions.push_back(subscriptions[k]);
        }

        // Post SAMP ok
        post_samp_ok(sock);

    } // endif: valid client found

    // Signal if the client was unknown
    #if defined(G_CONSOLE_ERRORS)
    else {
        std::cout << "*** ERROR: GVOHub::request_declare_subscriptions: ";
        std::cout << "Unable to find client \"" << get_client_key(xml);
        std::cout << "\"" << std::endl;
    }
    #endif

    // Return
    return;
}


/***********************************************************************//**
 * @brief Handles XML-RPC callback setting requests
 *
 * @param[in] xml XML message sent by client.
 * @param[in] sock Socket.
 *
 * Handles all incoming XML-RPC callback setting requests.
 ***************************************************************************/
void GVOHub::request_set_xml_rpc_callback(const GXml&      xml,
                                          const socklen_t& sock)
{
    // Header
    #if defined(G_CONSOLE_DUMP)
    std::cout << "GVOHub::request_set_xml_rpc_callback" << std::endl;
    #endif

    // Get client index
    int i = get_client_index(xml);

    // Continue only if index is valid
    if (i != -1) {

        // Set callback URL
        m_clients[i].url = get_callback_url(xml);

        // Post SAMP ok
        post_samp_ok(sock);

    } // endif: client found

    // Signal if the client was unknown
    #if defined(G_CONSOLE_ERRORS)
    else {
        std::cout << "*** ERROR: GVOHub::request_set_xml_rpc_callback: ";
        std::cout << "Unable to find client \"" << get_client_key(xml);
        std::cout << "\"" << std::endl;
    }
    #endif

    // Return
    return;
}


/***********************************************************************//**
 * @brief Handles subscriptions getting requests
 *
 * @param[in] xml XML message sent by client.
 * @param[in] sock Socket.
 *
 * Handles all incoming subscriptions getting requests.
 ***************************************************************************/
void GVOHub::request_get_subscriptions(const GXml& xml, const socklen_t& sock)
{
    // Header
    #if defined(G_CONSOLE_DUMP)
    std::cout << "GVOHub::request_get_subscriptions" << std::endl;
    #endif

    // Get the client name
    std::string client_name = get_client_name(xml);

    // Declare response
    std::string response = "";

    // Set response
    response.append("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
    response.append("<methodResponse>\n");
    response.append("<params>\n");
    response.append("  <param><value><struct>\n");

    // If key is "gammalib_hub" then return some special information
    if (client_name == "gammalib_hub") {
        // Do nothing
    }

    // ... handle all other clients
    else {

        // Get client index
        int i = get_client_index(client_name);

        // Continue only if index is valid
        if (i != -1) {

            // Append all subscriptions
            for (int k = 0; k < m_clients[i].subscriptions.size(); ++k) {
                response.append("    <member>\n");
                response.append("      <name>");
                response.append(m_clients[i].subscriptions[k]);
                response.append("</name>\n");
                response.append("      <value></value>\n");
                response.append("    </member>\n");
            }

        }

        // Signal if the client was unknown
        #if defined(G_CONSOLE_ERRORS)
        else {
            std::cout << "*** ERROR: GVOHub::request_get_subscriptions: ";
            std::cout << "Unable to find client \"" << get_client_key(xml);
            std::cout << "\"" << std::endl;
        }
        #endif

    }

    // Finish response
    response.append("  </struct></value></param>\n");
    response.append("</params>\n");
    response.append("</methodResponse>\n");

    // Post response
    post_string(response, sock);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Handles registered client information requests
 *
 * @param[in] xml XML message sent by client.
 * @param[in] sock Socket.
 *
 * Handles all incoming registered client information requests.
 ***************************************************************************/
void GVOHub::request_get_registered_clients(const GXml& xml, const socklen_t& sock)
{
    // Header
    #if defined(G_CONSOLE_DUMP)
    std::cout << "GVOHub::request_get_registered_clients" << std::endl;
    #endif

    // Get client key
    std::string key = get_client_key(xml);

    // Declare message
    std::string msg = "";

    // Set response
    msg.append("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
    msg.append("<methodResponse>\n");
    msg.append("<params>\n");
    msg.append("  <param><value><array><data>\n");
    msg.append("    <value>gammalib_hub</value>\n");

    // Loop over all clients. Do not send back current client's registration
    for (int i = 0; i < m_clients.size(); ++i) {
        if (m_clients[i].private_key != key) {
            msg.append("    <value>");
            msg.append(m_clients[i].reference);
            msg.append("</value>\n");
        }
    }

    // Finish response
    msg.append("  </data></array></value></param>\n");
    msg.append("</params>\n");
    msg.append("</methodResponse>\n");

    // Post response
    post_string(msg, sock);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Handles subscribed client information requests
 *
 * @param[in] xml XML message sent by client.
 * @param[in] sock Socket.
 *
 * Handles all incoming subscribed client information requests.
 ***************************************************************************/
void GVOHub::request_get_subscribed_clients(const GXml& xml, const socklen_t& sock)
{
    // Header
    #if defined(G_CONSOLE_DUMP)
    std::cout << "GVOHub::request_get_subscribed_clients" << std::endl;
    #endif

    // Declare message
    std::string msg = "";    

    msg.append("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
    msg.append("<methodResponse>\n");
    msg.append("<params>\n");
    msg.append("  <param><value><struct>\n");
    msg.append("    <member>\n");
    msg.append("      <name>c0</name>\n");
    msg.append("      <value><struct></struct></value>\n");
    msg.append("    </member>\n");
    msg.append("    <member>\n");
    msg.append("      <name>gammalib_hub</name>\n");
    msg.append("      <value><struct></struct></value>\n");
    msg.append("    </member>\n");
    msg.append("  </struct></value></param>\n");
    msg.append("</params>\n");
    msg.append("</methodResponse>");

    // Post response
    post_string(msg, sock);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Handles a metadata requests
 *
 * @param[in] xml XML message sent by client.
 * @param[in] sock Socket.
 *
 * Handles all metadata requests.
 ***************************************************************************/
void GVOHub::request_get_metadata(const GXml& xml, const socklen_t& sock)
{
    // Header
    #if defined(G_CONSOLE_DUMP)
    std::cout << "GVOHub::request_get_metadata" << std::endl;
    #endif

    // Get the client name
    std::string client_name = get_client_name(xml);

    // Declare response
    std::string response = "";

    // Set response
    response.append("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
    response.append("<methodResponse>\n");
    response.append("<params>\n");
    response.append("  <param><value><struct>\n");
    response.append("    <member>\n");
    response.append("      <name>samp.status</name>\n");
    response.append("      <value>samp.ok</value>\n");
    response.append("    </member>\n");

    // If key is "gammalib_hub" then return some special information
    if (client_name == "gammalib_hub") {
		response.append("    <member>\n");
        response.append("      <name>samp.name</name>\n");
        response.append("      <value>gammalib_hub</value>\n");
        response.append("    </member>\n");
   		response.append("    <member>\n");
        response.append("      <name>samp.description.text</name>\n");
        response.append("      <value>GammaLib VO Hub</value>\n");
        response.append("    </member>\n");
		response.append("    <member>\n");
        response.append("      <name>samp.icon.url</name>\n");
        response.append("      <value>http://a.fsdn.com/allura/p/gammalib/icon</value>\n");
        response.append("    </member>\n");
		response.append("    <member>\n");
        response.append("      <name>samp.documentation.url</name>\n");
        response.append("      <value>http://cta.irap.omp.eu/gammalib/user_manual/modules/vo.html</value>\n");
        response.append("    </member>\n");
		response.append("    <member>\n");
        response.append("      <name>author.affiliation</name>\n");
        response.append("      <value>IRAP, Toulouse, France</value>\n");
        response.append("    </member>\n");
		response.append("    <member>\n");
        response.append("      <name>author.email</name>\n");
        response.append("      <value>jurgen.knodlseder@irap.omp.eu</value>\n");
        response.append("    </member>\n");
		response.append("    <member>\n");
        response.append("      <name>author.name</name>\n");
        response.append("      <value>J. Knoedlseder, T. Louge</value>\n");
        response.append("    </member>\n");
		response.append("    <member>\n");
        response.append("      <name>home.page</name>\n");
        response.append("      <value>http://cta.irap.omp.eu/gammalib/</value>\n");
        response.append("    </member>\n");
    }

    // ... otherwise get the index
    else {

        // Get client index
        int i = get_client_index(client_name);

        // Continue only if index is valid
        if (i != -1) {

            // Append response
            response.append("    <member>\n");
            response.append("      <name>samp.name</name>\n");
            response.append("      <value>"+m_clients[i].name+"</value>\n");
            response.append("    </member>\n");
            response.append("    <member>\n");
            response.append("      <name>samp.description.text</name>\n");
            response.append("      <value>"+m_clients[i].description+"</value>\n");
            response.append("    </member>\n");
            response.append("    <member>\n");
            response.append("      <name>samp.icon.url</name>\n");
            response.append("      <value>"+m_clients[i].icon+"</value>\n");
            response.append("    </member>\n");
            response.append("    <member>\n");
            response.append("      <name>samp.documentation.url</name>\n");
            response.append("      <value>"+m_clients[i].documentation+"</value>\n");
            response.append("    </member>\n");
            response.append("    <member>\n");
            response.append("      <name>author.affiliation</name>\n");
            response.append("      <value>"+m_clients[i].affiliation+"</value>\n");
            response.append("    </member>\n");
            response.append("    <member>\n");
            response.append("      <name>author.email</name>\n");
            response.append("      <value>"+m_clients[i].email+"</value>\n");
            response.append("    </member>\n");
            response.append("    <member>\n");
            response.append("      <name>author.name</name>\n");
            response.append("      <value>"+m_clients[i].author_name+"</value>\n");
            response.append("    </member>\n");
            response.append("    <member>\n");
            response.append("      <name>home.page</name>\n");
            response.append("      <value>"+m_clients[i].homepage+"</value>\n");
            response.append("    </member>\n");

        } // endif: index was valid

        // Signal if the client was unknown
        #if defined(G_CONSOLE_ERRORS)
        else {
            std::cout << "*** ERROR: GVOHub::request_get_metadata: ";
            std::cout << "Unable to find client \"" << get_client_key(xml);
            std::cout << "\"" << std::endl;
        }
        #endif

    } // endelse: client was not the hub

    // Finish response
    response.append("  </struct></value></param>\n");
    response.append("</params>\n");
    response.append("</methodResponse>\n");

    // Post response
    post_string(response, sock);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Handle request to notify all clients
 *
 * @param[in] xml XML message sent by client.
 * @param[in] sock Socket.
 *
 * Handles requests to notify all clients. The method will loop over all
 * clients to check the subscriptions of each client, and send a notification
 * to all clients that are subscribed to the requested message type. 
 ***************************************************************************/
void GVOHub::request_notify_all(const GXml& xml, const socklen_t& sock)
{
    // Header
    #if defined(G_CONSOLE_DUMP)
    std::cout << "GVOHub::request_notify_all entrance" << std::endl;
    #endif

    // Get message type
    std::string mtype = get_mtype(xml);

    #if defined(G_CONSOLE_DUMP)
    std::cout << "mtype:"+mtype << std::endl;
    #endif

    // Post void message
    post_samp_void(sock);

    // Loop over all clients
    for (int i = 0; i < m_clients.size(); ++i) {

        // Loop over all subscriptions
        for (int j = 0; j < m_clients[i].subscriptions.size(); ++j) {

            // Continue only if message type matches client's subscription
            if (mtype == m_clients[i].subscriptions[j]) {

                // Image loading?
                if (mtype == "image.load.fits") {
                    notify_image_load(m_clients[i], xml);
		        }

                // ... or VO table loading?
                else if (mtype == "table.load.votable") {
                    notify_table_load(m_clients[i], xml);
		        }

            } // endif: message type matched client subscription

        } // endfor: looped over all message types

    } // endfor: looped over all clients

    // Return
    return;

}


/***********************************************************************//**
 * @brief Handles Hub shutdown requests
 *
 * @param[in] sock Socket.
 *
 * Handles all incoming Hub shutdown requests.
 ***************************************************************************/
void GVOHub::request_shutdown(const socklen_t& sock)
{
    // Header
    #if defined(G_CONSOLE_DUMP)
    std::cout << "GVOHub::request_shutdown" << std::endl;
    #endif

    // Set shutdown flag
    m_shutdown = true;

    // Post SAMP ok
    post_samp_ok(sock);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Extract client key from XML request
 *
 * @param[in] xml XML message sent by client.
 *
 * Extracts the clients key from the XML request.
 ***************************************************************************/
std::string GVOHub::get_client_key(const GXml& xml) const
{
    // Header
    #if defined(G_CONSOLE_DUMP)
    std::cout << "GVOHub::get_client_key" << std::endl;
    #endif

    // Initialise response
    std::string client_key = "";

    // Get the client's private key
    const GXmlNode* node = xml.element("methodCall > params > param > value");
    if (node != NULL) {
        const GXmlText* text = static_cast<const GXmlText*>((*node)[0]);
        if (text != NULL) {
            client_key = text->text();
        }
    }

    // Return key
    return client_key;
}


/***********************************************************************//**
 * @brief Extract client index in shared memory from XML request
 *
 * @param[in] xml XML message sent by client.
 * @return Client index (-1 if not found).
 *
 * Extracts the client index from the XML request. The method returns -1 if
 * no client was found.
 ***************************************************************************/
int GVOHub::get_client_index(const GXml& xml) const
{
    // Header
    #if defined(G_CONSOLE_DUMP)
    std::cout << "GVOHub::get_client_index" << std::endl;
    #endif

    // Initialise index
    int index = -1;

    // Get the client's private key
    std::string key = get_client_key(xml);

    // Continue only if string is not empty
    if (!key.empty()) {

        // Searches for key in clients
        for (int i = 0; i < m_clients.size(); ++i) {
            if (m_clients[i].private_key == key) {
                index = i;
                break;
            }
        }

    } // endif: key was not empty

    // Return index
    return index;
}


/***********************************************************************//**
 * @brief Extract client index in shared memory from XML request
 *
 * @param[in] reference Client reference.
 * @return Client index (-1 if not found).
 *
 * Extracts the client index from the client reference. The method returns -1
 * if no client was found.
 ***************************************************************************/
int GVOHub::get_client_index(const std::string& reference) const
{
    // Header
    #if defined(G_CONSOLE_DUMP)
    std::cout << "GVOHub::get_client_index" << std::endl;
    #endif

    // Initialise index
    int index = -1;

    // Search for index
    for (int i = 0; i < m_clients.size(); ++i) {
        if (m_clients[i].reference == reference) {
            index = i;
            break;
        }
    }

    // Return index
    return index;
}


/***********************************************************************//**
 * @brief Returns value for a named parameter
 *
 * @param[in] node XML node containing the response values.
 * @param[in] name Parameter name.
 * @return Parameter value.
 *
 * Returns the @p value for the parameter @p name in a list of @p member
 * XML nodes with the structure
 *
 *     <member>
 *       <name>name</name>
 *       <value>value</value>
 *     </member>
 *
 * If the specified parameter @p name was not found, or the the XML @p node
 * is NULL, an empty string is returned.
 ***************************************************************************/
std::string GVOHub::get_response_value(const GXmlNode*    node,
                                       const std::string& name) const
{
    // Declare empty value
    std::string value = "";

    // Search for parameter name, and if found, return its value
    if (node != NULL) {
        #if defined(G_CONSOLE_DUMP)
        std::cout << "GVOHub::get_response_value parsing non null node" << std::endl;
        #endif
        int num = node->elements("member");            
        for (int i = 0; i < num; ++i) {
            const GXmlNode* member = node->element("member", i);
            std::string one_name;
            std::string one_value;
            gammalib::xml_get_name_value_pair(member, one_name, one_value);
            #if defined(G_CONSOLE_DUMP)
            std::cout << "GVOHub::get_response_value parsing "+one_name+ " " + one_value << std::endl;
            #endif
            if (one_name == name) {
                value = one_value;
                break;
            }
        }

    } else {
        #if defined(G_CONSOLE_DUMP)
        std::cout << "GVOHub::get_response_value parsing NULL node" << std::endl;
        #endif
    }
    // Return value
    return value;
}


/***********************************************************************//**
 * @brief Returns subscriptions from XML document
 *
 * @param[in] xml client query XML document.
 * @return List of subscriptions.
 ***************************************************************************/
std::vector<std::string> GVOHub::get_subscriptions(const GXml& xml) const                                         
{
    // Declare subscriptions
    std::vector<std::string> subscriptions;

    // Get the client's private key
    const GXmlNode* node = xml.element("methodCall > params > param[1] > "
                                       "value > struct");
    if (node != NULL) {
        int num = node->elements("member");
        for (int i = 0; i < num; ++i) {
            const GXmlNode* member = node->element("member", i);
            if (member != NULL) {
                const GXmlNode* name = member->element("name", 0);
                if (name != NULL) {
                    const GXmlText* text = static_cast<const GXmlText*>((*name)[0]);
                    if (text != NULL) {
                        subscriptions.push_back(text->text());
                    }
                }
            }
        }
    }

    // Return subscriptions
    return subscriptions;
}


/***********************************************************************//**
 * @brief Returns callback URL of client
 *
 * @param[in] xml client query XML document.
 * @return Callback URL of client.
 ***************************************************************************/
std::string GVOHub::get_callback_url(const GXml& xml) const
{
    // Declare URL
    std::string url;

    // Get the client's private key
    const GXmlNode* node = xml.element("methodCall > params > param[1] > "
                                       "value");
    if (node != NULL) {
		const GXmlText* text = static_cast<const GXmlText*>((*node)[0]);
        if (text != NULL) {
		    url.append(text->text());
        }
    }

    // Return URL
    return url;
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
 * @brief Create the lockfile
 *
 * Creates the SAMP lockfile and fill it with information. Implements IVOA
 * standard REC-SAMP-1.3-20120411.
 ***************************************************************************/
void GVOHub::create_samp_file(void) const
{
    // Get lockfile URL
    std::string lockurl = get_hub_lockfile();

    // Open SAMP lockfile. Continue only if opening was successful
    FILE* fptr = fopen(lockurl.c_str(), "w");

    // Successful?
    if (fptr != NULL) {

        // Write lockfile
        fprintf(fptr, "# SAMP lockfile\n");
        fprintf(fptr, "# Required keys:\n");
        fprintf(fptr, "samp.secret=%s\n", m_secret.c_str());
        fprintf(fptr, "samp.hub.xmlrpc.url=%s\n", hub_url().c_str());
        fprintf(fptr, "samp.profile.version=%s\n", m_version.c_str());
        fprintf(fptr, "# Info stored by hub for some private reason:\n");
        fprintf(fptr, "gammalib.hubid=%s\n", m_hub_id.c_str());

        // Close lockfile
        fclose(fptr);

    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete the lockfile
 *
 * Deletes the SAMP lockfile on disk.
 ***************************************************************************/
void GVOHub::delete_samp_file(void) const
{
    // Get lockfile URL
    std::string lockurl = get_hub_lockfile();

    // Delete lockfile
    std::remove(lockurl.c_str());

    // Return
    return;
}


/***********************************************************************//**
 * @brief Get Hub socket
 *
 * @return Hub socket (-1 if no socket was found).
 *
 * Returns Hub socket. Starting from an initial port of 2526 the method
 * searches for the next free port and binds it to a socket.
 ***************************************************************************/
int GVOHub::get_socket(void)
{
    // Initialise socket, port and socket address structure
    int    sock = -1;
    int    port = 2526;
    struct sockaddr_in serv_addr;

    // Loop over 100 ports at most
    for (int i = 0; i < 100; ++i) {

        // Clean TCP/IP structure
        std::memset(&serv_addr, 0, sizeof(serv_addr));

        // Set TCP/IP structure
        serv_addr.sin_family      = AF_INET;
        serv_addr.sin_addr.s_addr = inet_addr(m_hub_host.c_str());
        serv_addr.sin_port        = htons(port);

        // Create Hub socket
        sock = socket(AF_INET, SOCK_STREAM, 0);

        // Creation of hub main socket
        if (sock < 0) {
            std::string msg = "Unable to create Hub socket. Errno="+
                              gammalib::str(errno);
            throw GException::runtime_error(G_GET_SOCKET, msg);
        }

        // Now bind socket to the address and port
        if (bind(sock, (struct sockaddr*)&serv_addr, sizeof(serv_addr)) < 0) {
            if (errno == EADDRINUSE) {
                port++;
                close(sock);
            }
            else {
                std::string msg = "Unable to bind Hub socket to address "+
                                  m_hub_host+":"+gammalib::str(port)+
                                  ". Errno="+gammalib::str(errno);
                throw GException::runtime_error(G_GET_SOCKET, msg);
            }
        }
        else {
            break;
        }

    } // endfor: looped over ports

    // Store port as string
    m_hub_port = gammalib::str(port);

    // Return socket
    return sock;
}


/***********************************************************************//**
 * @brief Post string content to client
 *
 * @param[in] content String content to post.
 * @param[in] sock Socket.
 *
 * Posts the content of a string to a client.
 ***************************************************************************/
void GVOHub::post_string(const std::string& content, const socklen_t& sock) const
{
    // Continue only if socket is valid
    if (sock != -1) {

        // Determine content length
        int length = content.length();

        // Set prefix (note that the response has no POST!)
        std::string prefix = "HTTP/1.1 200 OK\n"
                             "Connection: close\n"
                             "Content-Type: text/xml\n"
                             "Content-Length: "+gammalib::str(length)+"\n\n";


        // Build post string
        std::string post = prefix + content;

        // Dump message to post
        #if defined(G_SHOW_MESSAGE)
        std::cout << std::endl;
        std::cout << "GVOHub response:" << std::endl;
        std::cout << "================" << std::endl;
        std::cout << post << std::endl;
        #endif

        // Send content to socket
        bool done = false;
        do {
            int length      = post.length();
            int sent_length = send(sock, post.c_str(), length, 0);
            if (sent_length < length) {
                post = post.substr(sent_length, std::string::npos);
            }
            else {
                done = true;
            }
        } while (!done);

    } // endif: socket was valid

    // Return
    return;
}


/***********************************************************************//**
 * @brief Post SAMP ok massage to client
 *
 * @param[in] sock Socket.
 *
 * Posts a SAMP OK message to a client.
 ***************************************************************************/
void GVOHub::post_samp_ok(const socklen_t& sock) const
{
    // Declare response
    std::string response = "";

    // Compose response
    response.append("<?xml version='1.0' encoding='UTF-8'?>\n");
    response.append("<methodResponse>\n");
    response.append("<params>\n");
    response.append("  <param><value><struct>\n");
    response.append("    <member>\n");
    response.append("      <name>samp.status</name>\n");
    response.append("      <value>samp.ok</value>\n");
    response.append("    </member>\n");
    response.append("  </struct></value></param>\n");
    response.append("</params>\n");
    response.append("</methodResponse>\n");

    // Post string
    post_string(response, sock);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Post SAMP void massage to client
 *
 * @param[in] sock Socket.
 *
 * Posts a void SAMP message to a client. Void messages can typically be
 * discared by the client.
 ***************************************************************************/
void GVOHub::post_samp_void(const socklen_t& sock) const
{
    // Declare response
    std::string response = "";

    // Compose response
    response.append("<?xml version='1.0' encoding='UTF-8'?>\n");
    response.append("<methodResponse>\n");
    response.append("  <params>\n");
    response.append("    <param>\n");
    response.append("      <value/>\n");
    response.append("    </param>\n");
    response.append("  </params>\n");
    response.append("</methodResponse>\n");

    // Post string
    post_string(response, sock);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Send notification to client
 *
 * @param[in] url URL.
 * @param[in] notification Notification.
 *
 * Send a @p notification to a @p url.
 ***************************************************************************/
void GVOHub::notify(const std::string& url,
                    const std::string& notification) const
{
    // Initalise socket
    int socket = -1;

    // Extract host, port and path from URL
    std::string host = "127.0.0.1";
    std::string port = "2525";
    std::string path = "xmlrpc";
    if (url.compare(0, 7, "http://") == 0) {
        size_t length;
        size_t start = 7;
        size_t stop  = url.find(":", start);
        size_t end   = std::string::npos;
        if (stop != std::string::npos) {
            length = stop - start;
        }
        else {
            length = std::string::npos;
        }
        host = url.substr(start, length);
        if (stop != std::string::npos) {
            stop++;
            end  = url.find("/", stop);
            if (end != std::string::npos) {
                length = end - stop;
            }
            else {
                length = std::string::npos;
            }
            port = url.substr(stop, length);
        }
        if (end != std::string::npos) {
            end++;
            length = url.length() - end;
            path   = url.substr(end, length);
        }
    }

    // Set hints
    struct addrinfo hints;
    std::memset(&hints, 0, sizeof(hints));
    hints.ai_family   = AF_INET;
    hints.ai_socktype = SOCK_STREAM;

    // Get server information
    struct addrinfo* servinfo;
    if (getaddrinfo(host.c_str(), port.c_str(), &hints, &servinfo) == 0) {

        // Loop through all the results and connect to the first we can
        for (struct addrinfo* ptr = servinfo; ptr != NULL; ptr = ptr->ai_next) {

            // Create socket
            socket = ::socket(ptr->ai_family,
                              ptr->ai_socktype,
                              ptr->ai_protocol);

            // Connect to socket if socket is valid
            if (socket != -1) {
                if (::connect(socket,
                              ptr->ai_addr,
                              ptr->ai_addrlen) == -1) {
                    close(socket);
                    socket = -1;
                }
                else {
                    // Socket connection successful, break now
                    break;
                }
            }

        } // endfor: looped through all results

    } // endif: server information was valid

    // Continue only if socket is valid
    if (socket != -1) {

        // Determine notification length
        int length = notification.length();

        // Set prefix
        std::string prefix = "POST /"+path+" HTTP/1.0\n"
                             "Connection: close\n"
                             "User-Agent: GammaLib\n"
                             "Content-Type: text/xml\n"
                             "Content-Length: "+gammalib::str(length)+"\n\n";

        // Build post string
        std::string post = prefix + notification;

        // Dump message to post
        #if defined(G_SHOW_MESSAGE)
        std::cout << std::endl;
        std::cout << "GVOHub sends following notification:" << std::endl;
        std::cout << "====================================" << std::endl;
        std::cout << post << std::endl;
        #endif

        // Send content to socket
        bool done = false;
        do {
            int length      = post.length() + 1; // +1 for terminating 0
            int sent_length = send(socket, post.c_str(), length, 0);
            if (sent_length < length) {
                post = post.substr(sent_length, std::string::npos);
            }
            else {
                done = true;
            }
        } while (!done);

        // Close socket
        close(socket);

    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Notify client about another client registering at Hub
 *
 * @param[in] client Client to notify.
 * @param[in] reference Reference to client which registered.
 *
 * Notify a client about another client that has registered at the Hub. The
 * method sends a "samp.hub.event.register" message to the @p client.
 ***************************************************************************/
void GVOHub::notify_register(const GVOHub::client& client,
                             const std::string&    reference)
{
    // Declare notification
    std::string msg = "";

    // Set response for samp.hub.event.register
    msg.append("<?xml version=\"1.0\"?>\n");
    msg.append("<methodCall>\n");
    msg.append("<methodName>samp.client.receiveNotification</methodName>\n");
    msg.append("<params>\n");
    msg.append("  <param><value>"+client.private_key+"</value></param>\n");
    msg.append("  <param><value>"+m_hub_id+"</value></param>\n");
    msg.append("  <param><value><struct>\n");
    msg.append("    <member>\n");
    msg.append("      <name>samp.mtype</name>\n");
    msg.append("      <value>samp.hub.event.register</value>\n");
    msg.append("    </member>\n");
    msg.append("    <member>\n");
    msg.append("      <name>samp.params</name>\n");
    msg.append("      <value><struct>\n");
    msg.append("        <member>\n");
    msg.append("          <name>id</name>\n");
    msg.append("          <value>"+reference+"</value>\n");
    msg.append("        </member>\n");
    msg.append("      </struct></value>\n");
    msg.append("    </member>\n");
    msg.append("  </struct></value></param>\n");
    msg.append("</params>\n");
    msg.append("</methodCall>\n");

    // Post notification
    notify(client.url, msg);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Notify client about another client unregistering at Hub
 *
 * @param[in] client Client to notify.
 * @param[in] reference Reference to client who unregistered.
 *
 * Notify a client about another client that has unregistered at the Hub. The
 * method sends a "samp.hub.event.unregister" message to the @p client.
 ***************************************************************************/
void GVOHub::notify_unregister(const GVOHub::client& client,
                               const std::string&    reference)
{
    // Declare notification
    std::string msg = "";

    // Set response for samp.hub.event.unregister
    msg.append("<?xml version=\"1.0\"?>\n");
    msg.append("<methodCall>\n");
    msg.append("<methodName>samp.client.receiveNotification</methodName>\n");
    msg.append("<params>\n");
    msg.append("  <param><value>"+client.private_key+"</value></param>\n");
    msg.append("  <param><value>"+m_hub_id+"</value></param>\n");
    msg.append("  <param><value><struct>\n");
    msg.append("    <member>\n");
    msg.append("      <name>samp.mtype</name>\n");
    msg.append("      <value>samp.hub.event.unregister</value>\n");
    msg.append("    </member>\n");
    msg.append("    <member>\n");
    msg.append("      <name>samp.params</name>\n");
    msg.append("      <value><struct>\n");
    msg.append("        <member>\n");
    msg.append("          <name>id</name>\n");
    msg.append("          <value>"+reference+"</value>\n");
    msg.append("        </member>\n");
    msg.append("      </struct></value>\n");
    msg.append("    </member>\n");
    msg.append("  </struct></value></param>\n");
    msg.append("</params>\n");
    msg.append("</methodCall>\n");

    // Post notification
    notify(client.url, msg);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Notify client about the metadata of another client
 *
 * @param[in] client Client to notify.
 * @param[in] reference Reference of client for which metadata should be sent.
 *
 * Notify a client about another the metadata of another client. The method
 * sends a "samp.hub.event.metadata" message to the @p client.
 ***************************************************************************/
void GVOHub::notify_metadata(const GVOHub::client& client,
                             const std::string&    reference)
{
    // Get index of reference client
    int k;
    for (k = 0; k < m_clients.size(); ++k) {
        if (m_clients[k].reference == reference) {
            break;
        }
    }

    // Continue only if client index is valid
    if (k < m_clients.size()) {

        // Declare notification
        std::string msg = "";

        // Set response for samp.hub.event.metadata
        msg.append("<?xml version=\"1.0\"?>\n");
        msg.append("<methodCall><methodName>samp.client.receiveNotification</methodName>\n");
        msg.append("<params>\n");
        msg.append("  <param><value>"+client.private_key+"</value></param>\n");
        msg.append("  <param><value>"+m_hub_id+"</value></param>\n");
        msg.append("  <param><value><struct>\n");
        msg.append("    <member>\n");
        msg.append("      <name>samp.mtype</name>\n");
        msg.append("      <value>samp.hub.event.metadata</value>\n");
        msg.append("    </member>\n");
        msg.append("    <member>\n");
        msg.append("      <name>samp.params</name>\n");
        msg.append("      <value><struct>\n");
        msg.append("        <member>\n");
        msg.append("          <name>id</name>\n");
        msg.append("          <value>"+m_clients[k].reference+"</value>\n");
        msg.append("        </member>\n");
        msg.append("        <member>\n");
        msg.append("          <name>metadata</name>\n");
        msg.append("          <value><struct>\n");
        msg.append("            <member>\n");
        msg.append("              <name>samp.name</name>\n");
        msg.append("              <value>"+m_clients[k].name+"</value>\n");
        msg.append("            </member>\n");
        msg.append("            <member>\n");
        msg.append("              <name>samp.description.text</name>\n");
        msg.append("              <value>"+m_clients[k].description+"</value>\n");
        msg.append("            </member>\n");
        msg.append("            <member>\n");
        msg.append("              <name>samp.icon.url</name>\n");
        msg.append("              <value>"+m_clients[k].icon+"</value>\n");
        msg.append("            </member>\n");
        msg.append("            <member>\n");
        msg.append("              <name>samp.documentation.url</name>\n");
        msg.append("              <value>"+m_clients[k].documentation+"</value>\n");
        msg.append("            </member>\n");
        msg.append("            <member>\n");
        msg.append("              <name>author.affiliation</name>\n");
        msg.append("              <value>"+m_clients[k].affiliation+"</value>\n");
        msg.append("            </member>\n");
        msg.append("            <member>\n");
        msg.append("              <name>author.name</name>\n");
        msg.append("              <value>"+m_clients[k].author_name+"</value>\n");
        msg.append("            </member>\n");
        msg.append("            <member>\n");
        msg.append("              <name>author.email</name>\n");
        msg.append("              <value>"+m_clients[k].email+"</value>\n");
        msg.append("            </member>\n");
        msg.append("            <member>\n");
        msg.append("              <name>home.page</name>\n");
        msg.append("              <value>"+m_clients[k].homepage+"</value>\n");
        msg.append("            </member>\n");
        msg.append("          </struct></value>\n");
        msg.append("        </member>\n");				
        msg.append("      </struct></value>\n");
        msg.append("    </member>\n");
        msg.append("  </struct></value></param>\n");
        msg.append("</params>\n");
        msg.append("</methodCall>\n");

        // Post notification
        notify(client.url, msg);

    } // endif: reference client was valid

    // Return
    return;
}


/***********************************************************************//**
 * @brief Notify client about image loading
 *
 * @param[in] client Client.
 * @param[in] xml Message from which image information is extracted.
 *
 * Notify a client that a FITS image is avialable for loading. The method
 * sends a "image.load.fits" message to the @p client. Information about the
 * image is extracted from the @p xml document.
 ***************************************************************************/
void GVOHub::notify_image_load(const GVOHub::client& client,
                               const GXml&           xml)
{
	// Get client response
    const GXmlNode* node = xml.element("methodCall > params > param[1] > "
                                       "value > struct > member > value > "
                                       "struct");

    // Extract image information from client response
	std::string name = get_response_value(node, "name");
	std::string url  = get_response_value(node, "url");
	std::string id   = get_response_value(node, "image-id");

    // Declare notification
    std::string msg = "";

	// Set notification for image.load.fits
    msg.append("<?xml version=\"1.0\"?>\n");
    msg.append("<methodCall>\n");
    msg.append("<methodName>samp.client.receiveNotification</methodName>\n");
    msg.append("<params>\n");
    msg.append("  <param><value>"+client.private_key+"</value></param>\n");
    msg.append("  <param><value>"+m_hub_id+"</value></param>\n");
    msg.append("  <param><value><struct>\n");
    msg.append("    <member>\n");
    msg.append("      <name>samp.mtype</name>\n");
    msg.append("      <value>image.load.fits</value>\n");
    msg.append("    </member>\n");
    msg.append("    <member>\n");
    msg.append("      <name>samp.params</name>\n");
    msg.append("      <value><struct>\n");
    msg.append("        <member>\n");
    msg.append("          <name>name</name>\n");
    msg.append("          <value>"+name+"</value>\n");
    msg.append("        </member>\n");
    msg.append("        <member>\n");
    msg.append("          <name>image-id</name>\n");
    msg.append("          <value>"+id+"</value>\n");
    msg.append("        </member>\n");
    msg.append("        <member>\n");
    msg.append("          <name>url</name>\n");
    msg.append("          <value>"+url+"</value>\n");
    msg.append("        </member>\n");
    msg.append("      </struct></value>\n");
    msg.append("    </member>\n");
    msg.append("  </struct></value></param>\n");
    msg.append("</params>\n");
    msg.append("</methodCall>\n");

    // Post notification
    notify(client.url, msg);

    // Return
    return;
}

/***********************************************************************//**
 * @brief Notify client about VO table loading
 *
 * @param[in] client Client.
 * @param[in] xml Message from which table information is extracted.
 *
 * Notify a client that a VO table is avialable for loading. The method sends
 * sends a "table.load.votable" message to the @p client. Information about
 * the VO table is extracted from the @p xml document.
 ***************************************************************************/
void GVOHub::notify_table_load(const GVOHub::client& client,
                               const GXml&           xml)
{
	// Get client response
    const GXmlNode* node = xml.element("methodCall > params > param[1] > "
                                       "value > struct > member > value > "
                                       "struct");

    // Extract image information from client response
	std::string name = get_response_value(node, "name");
	std::string url  = get_response_value(node, "url");
	std::string id   = get_response_value(node, "table-id");

    // Declare notification
    std::string msg = "";

	// Set notification for table.load.votable
    msg.append("<?xml version=\"1.0\"?>\n");
    msg.append("<methodCall>\n");
    msg.append("<methodName>samp.client.receiveNotification</methodName>\n");
    msg.append("<params>\n");
    msg.append("  <param><value>"+client.private_key+"</value></param>\n");
    msg.append("  <param><value>"+m_hub_id+"</value></param>\n");
    msg.append("  <param><value><struct>\n");
    msg.append("    <member>\n");
    msg.append("      <name>samp.mtype</name>\n");
    msg.append("      <value>table.load.votable</value>\n");
    msg.append("    </member>\n");
    msg.append("    <member>\n");
    msg.append("      <name>samp.params</name>\n");
    msg.append("      <value><struct>\n");
    msg.append("        <member>\n");
    msg.append("          <name>name</name>\n");
    msg.append("          <value>"+name+"</value>\n");
    msg.append("        </member>\n");
    msg.append("        <member>\n");
    msg.append("          <name>table-id</name>\n");
    msg.append("          <value>"+id+"</value>\n");
    msg.append("        </member>\n");
    msg.append("        <member>\n");
    msg.append("          <name>url</name>\n");
    msg.append("          <value>"+url+"</value>\n");
    msg.append("        </member>\n");
    msg.append("      </struct></value>\n");
    msg.append("    </member>\n");
    msg.append("  </struct></value></param>\n");
    msg.append("</params>\n");
    msg.append("</methodCall>\n");

    // Post notification
    notify(client.url, msg);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Generates random string of characters
 *
 * @param[in] length Length of random string
 ***************************************************************************/
std::string GVOHub::random_string(const size_t& length) const
{
    srand(time(0));
    std::string str = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmno"
                      "pqrstuvwxyz";
    while(str.size() != length) {
        int pos = ((rand() % (str.size() - 1)));
        str.erase(pos, 1);
    }

    // Return string
    return str;
}


/***********************************************************************//**
 * @brief Extract mtype from XML request
 *
 * @param[in] xml XML message sent by client.
 *
 * Extracts the mtype identifying calling message from the XML request.
 ***************************************************************************/
std::string GVOHub::get_mtype(const GXml& xml) const
{
    // Header
    #if defined(G_CONSOLE_DUMP)
    std::cout << "In GVOHub::get_mtype" << std::endl;
    #endif

    // Get the XML node containing the client's private key
    const GXmlNode* node = xml.element("methodCall > params > param[1] > "
                                       "value > struct ");

    // Get the client's private key
    std::string client_key = get_response_value(node, "samp.mtype");

    // Return key
    return client_key;
}


/***********************************************************************//**
 * @brief Return Hub URL
 *
 * @return Hub URL.
 *
 * Returns the XML-RPC endpoint for communication with the hub.
 ***************************************************************************/
std::string GVOHub::hub_url(void) const
{
    // Set Hub URL
    std::string hub_url = "http://"+m_hub_host+":"+m_hub_port+"/"+m_hub_path;

    // Return Hub URL
    return (hub_url);
}


/***********************************************************************//**
 * @brief Return client name from XML message
 *
 * @param[in] xml XML message sent by client.
 *
 * Returns client name from XML message.
 ***************************************************************************/
std::string GVOHub::get_client_name(const GXml& xml) const
{
    // Initialise client name
    std::string client_name = "";

    // Get the client name
    const GXmlNode* node = xml.element("methodCall > params > param[1] > value");
    if (node != NULL) {
        const GXmlText* text = static_cast<const GXmlText*>((*node)[0]);
        if (text != NULL) {
            client_name = text->text();
        }
    }

    // Return client name
    return client_name;
}
