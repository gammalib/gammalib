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
#include <cstring>         // std::memset() function
#include <fcntl.h>         // close() function
#include <netdb.h>         // getaddrinfo() function
#include <sys/socket.h>    // socket(), connect() functions
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
 * client is registered at the Hub. Metadata are sent.
 ***************************************************************************/
void GVOClient::connect(void)
{
    // Register to Hub
    register_to_hub();

    // Send meta data
    send_metadata();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Unregister client from SAMP Hub
 *
 * Disconnects the VO client from the Hub. The client is unregistered and the
 * socket to the Hub is closed.
 ***************************************************************************/
void GVOClient::disconnect(void)
{
    // Unregister from Hub
    unregister_from_hub();

    // Close socket
    if (m_socket != -1) {
        close(m_socket);
        m_socket = -1;
    }

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
 * @brief Signals if client is connected to Hub
 *
 * @return True if client is connected to Hub, false otherwise.
 *
 * Checks if a client is connected to Hub.
 ***************************************************************************/
bool GVOClient::isconnected(void) const
{
    // Return
    return false;
}


/***********************************************************************//**
 * @brief Returns Hub response a XML object
 *
 * @return Hub response.
 ***************************************************************************/
GXml GVOClient::response(void) const
{
    // Declare empty XML document
    GXml xml;

    // Continue only if connection has been established
    if (m_socket != -1) {

        // Receive a string
        std::string response = receive_string();

        // Find start of XML text
        size_t start = response.find("<?xml");

        // If found then convert text into XML document
        if (start != std::string::npos) {
            xml = GXml(response.substr(start, std::string::npos));
        }

    } // endif: connection has been established

    // Return XML document
    return xml;
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
    result.append("\n"+parformat("Name")+m_name);

    // Append Hub information
    result.append("\n"+parformat("Hub key")+m_secret);
    result.append("\n"+parformat("Hub URL")+m_hub_url);
    result.append("\n"+parformat("Hub host (port)")+m_hub_host+" ("+m_hub_port+")");
    result.append("\n"+parformat("SAMP protocol version")+m_version);

    // Append connection information
    result.append("\n"+parformat("Hub connection"));
    if (m_socket == -1) {
        result.append("no");
    }
    else {
        if (!m_client_key.empty()) {
            result.append("registered as \""+m_client_id);
            result.append("\" on Hub \""+m_hub_id);
            result.append("\" via socket "+str(m_socket));
        }
        else {
            result.append("established on socket "+str(m_socket));
        }
    }

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
    m_name = "GammaLib";
    m_secret.clear();
    m_hub_url.clear();
    m_hub_host.clear();
    m_hub_port.clear();
    m_version.clear();
    m_client_key.clear();
    m_hub_id.clear();
    m_client_id.clear();
    m_socket = -1;         // Signals no socket

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
    m_name       = client.m_name;
    m_secret     = client.m_secret;
    m_hub_url    = client.m_hub_url;
    m_hub_host   = client.m_hub_host;
    m_hub_port   = client.m_hub_port;
    m_version    = client.m_version;
    m_client_key = client.m_client_key;
    m_hub_id     = client.m_hub_id;
    m_client_id  = client.m_client_id;
    m_socket     = client.m_socket;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GVOClient::free_members(void)
{
    // Unregister from Hub
    unregister_from_hub();

    // Close socket
    if (m_socket != -1) {
        close(m_socket);
        m_socket = -1;
    }

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
    std::string lockurl = get_hub_lockfile();

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

            // Extract host and port from Hub URL
            if (m_hub_url.compare(0, 7, "http://") == 0) {
                size_t length;
                size_t start = 7;
                size_t stop  = m_hub_url.find(":", start);
                if (stop != std::string::npos) {
                    length = stop - start;
                }
                else {
                    length = std::string::npos;
                }
                m_hub_host = m_hub_url.substr(start, length);
                if (stop != std::string::npos) {
                    stop = stop + 1;
                    size_t end = m_hub_url.find("/", stop);
                    if (end != std::string::npos) {
                        length = end - stop;
                    }
                    else {
                        length = std::string::npos;
                    }
                    m_hub_port = m_hub_url.substr(stop, length);
                }
            }

            // Check for existence of mandatory tokens
            found = hashub();

        } // endif: SAMP lockfile opened

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
 * @brief Connect to SAMP Hub
 *
 * Connects to Hub by creating a socket and connecting to this socket. The
 * method expects that a Hub has already been found. If no Hub has been found
 * (i.e. hashub() returns false), the method does nothing.
 *
 * The method extracts the hostname and the port from the Hub endpoint URL.
 * If Hub connection fails, any created socket is closed.
 *
 * If the method is successful, m_socket will contain on exit a non-negative
 * number. If any failure occurs, m_socket will be set to -1.
 ***************************************************************************/
void GVOClient::connect_to_hub(void)
{
    // Close any existing socket
    if (m_socket != -1) {
        close(m_socket);
        m_socket = -1;
    }

    // Continue only if host and port information was found
    if (!m_hub_host.empty() && !m_hub_port.empty()) {

        // Set hints
        struct addrinfo hints;
        std::memset(&hints, 0, sizeof(hints));
        hints.ai_family   = AF_INET;
        hints.ai_socktype = SOCK_STREAM;

        // Get server information
        struct addrinfo* servinfo;
        if (getaddrinfo(m_hub_host.c_str(), m_hub_port.c_str(),
                        &hints, &servinfo) == 0) {

            // Loop through all the results and connect to the first
            // we can
            for (struct addrinfo* ptr = servinfo; ptr != NULL; ptr = ptr->ai_next) {

                // Create socket
                m_socket = socket(ptr->ai_family,
                                  ptr->ai_socktype,
                                  ptr->ai_protocol);

                // Connect to socket if socket is valid
                if (m_socket != -1) {
                    if (::connect(m_socket,
                                  ptr->ai_addr,
                                  ptr->ai_addrlen) == -1) {
                        close(m_socket);
                        m_socket = -1;
                    }
                }

            } // endfor: looped through all results

        } // endif: server information was valid

    } // endif: host and port information found

    // Return
    return;
}


/***********************************************************************//**
 * @brief Register client at SAMP Hub
 ***************************************************************************/
void GVOClient::register_to_hub(void)
{
    // Connect to Hub
    connect_to_hub();

    // Continue only if Hub connection has been established
    if (m_socket != -1) {

        // Declare message
        std::string msg = "";

        // Set metadata header
        msg.append("<?xml version=\"1.0\"?>\n");
        msg.append("<methodCall>\n");
        msg.append("<methodName>samp.hub.register</methodName>\n");
        msg.append("<params>\n");
        msg.append("<param><value><string>"+m_secret+"</string></value></param>\n");
        msg.append("</params>\n");
        msg.append("</methodCall>\n");

        // Post message
        post_string(msg);

        // Get Hub response
        GXml xml = response();

        // Extract Hub and client identifiers
        m_client_key = get_response_value(xml, "samp.private-key");
        m_hub_id     = get_response_value(xml, "samp.hub-id");
        m_client_id  = get_response_value(xml, "samp.self-id");

    } // endif: Hub connection has been established

    // Return
    return;
}


/***********************************************************************//**
 * @brief Unregister client from SAMP Hub
 ***************************************************************************/
void GVOClient::unregister_from_hub(void)
{
    // Connect to Hub
    connect_to_hub();

    // Continue only if Hub connection has been established
    if (m_socket != -1) {

        // Declare message
        std::string msg = "";

        // Set metadata header
        msg.append("<?xml version=\"1.0\"?>\n");
        msg.append("<methodCall>\n");
        msg.append("<methodName>samp.hub.unregister</methodName>\n");
        msg.append("<params>\n");
        msg.append("<param><value><string>"+m_client_key+"</string></value></param>\n");
        msg.append("</params>\n");
        msg.append("</methodCall>\n");

        // Post message
        post_string(msg);

        // Get Hub response
        GXml xml = response();

        // Reset Hub and client identifiers
        m_client_key.clear();
        m_hub_id.clear();
        m_client_id.clear();

    } // endif: Hub connection has been established

    // Return
    return;
}


/***********************************************************************//**
 * @brief Send client metadata to Hub
 ***************************************************************************/
void GVOClient::send_metadata(void)
{
    // Connect to Hub. This was needed to post the metadata to the Hub.
    // Apparently, the connection is lost after each write/read cycle.
    // Maybe this is socket standard?
    connect_to_hub();

    // Continue only if Hub connection has been established
    if (m_socket != -1) {

        // Declare message
        std::string msg = "";

        // Set metadata header
        msg.append("<?xml version=\"1.0\"?>\n");
        msg.append("<methodCall>\n");
        msg.append("<methodName>samp.hub.declareMetadata</methodName>\n");
        msg.append("<params>\n");
        msg.append("<param><value><string>"+m_client_key+"</string></value></param>\n");
        msg.append("<param><value><struct>\n");

        // Set SAMP name
        msg.append("<member>\n");
        msg.append("<name>samp.name</name>\n");
        msg.append("<value><string>"+m_name+"</string></value>\n");
        msg.append("</member>\n");

        // Set SAMP description text
        msg.append("<member>\n");
        msg.append("<name>samp.description.text</name>\n");
        msg.append("<value><string>GammaLib client</string></value>\n");
        msg.append("</member>\n");

        // Set SAMP icon URL
        msg.append("<member>\n");
        msg.append("<name>samp.icon.url</name>\n");
        msg.append("<value><string>http://a.fsdn.com/allura/p/gammalib/icon</string></value>\n");
        msg.append("</member>\n");

        // Set author affiliation
        msg.append("<member>\n");
        msg.append("<name>author.affiliation</name>\n");
        msg.append("<value><string>IRAP, Toulouse, France</string></value>\n");
        msg.append("</member>\n");

        // Set author e-mail
        msg.append("<member>\n");
        msg.append("<name>author.email</name>\n");
        msg.append("<value><string>jurgen.knodlseder@irap.omp.eu</string></value>\n");
        msg.append("</member>\n");

        // Set author name
        msg.append("<member>\n");
        msg.append("<name>author.name</name>\n");
        msg.append("<value><string>Juergen Knoedlseder</string></value>\n");
        msg.append("</member>\n");

        // Set metadata trailer
        msg.append("</struct></value></param>\n");
        msg.append("</params>\n");
        msg.append("</methodCall>\n");

        // Post message
        post_string(msg);

        // Get Hub response
        GXml xml = response();

    } // endif: Hub connection has been established

    // Return
    return;
}


/***********************************************************************//**
 * @brief Returns value for a SAMP Hub response parameter
 *
 * @param[in] xml Hub response XML document.
 * @param[in] name Parameter name.
 * @return Parameter value.
 *
 * Returns value for a SAMP Hub response parameter. If the specified
 * parameter was not found or if the response structure is not compliant,
 * an empty string is returned.
 ***************************************************************************/
std::string GVOClient::get_response_value(const GXml&        xml,
                                          const std::string& name) const
{
    // Declare value
    std::string value = "";

    // Search for value of specified member
    GXmlNode* node = xml.element("methodResponse", 0);
    if (node != NULL) {
        node = node->element("params", 0);
        if (node != NULL) {
            node = node->element("param", 0);
            if (node != NULL) {
                node = node->element("value", 0);
                if (node != NULL) {
                    node = node->element("struct", 0);
                    if (node != NULL) {
                        int num = node->elements("member");
                        for (int i = 0; i < num; ++i) {
                            GXmlNode* member = node->element("member", i);
                            std::string one_name;
                            std::string one_value;
                            get_name_value_pair(member, one_name, one_value);
                            if (one_name == name) {
                                value = one_value;
                                break;
                            }
                        }
                    }
                }
            }
        }
    }

    // Return value
    return value;
}


/***********************************************************************//**
 * @brief Extract name / value pair from XML node
 *
 * @param[in] node Pointer to XML node.
 * @param[out] name Name string.
 * @param[out] value Value string.
 *
 * Extracts a name / value pair from a XML node. If the XML node pointer is
 * NULL, the name and value strings will be empty.
 ***************************************************************************/
void GVOClient::get_name_value_pair(const GXmlNode* node,
                                    std::string&    name,
                                    std::string&    value) const
{
    // Clear name and value strings
    name.clear();
    value.clear();

    // Continue only if node is valid
    if (node != NULL) {

        // Get name node and extract text content
        GXmlNode* ptr = node->element("name", 0);
        if (ptr != NULL) {
            GXmlText* text = static_cast<GXmlText*>(ptr->child(0));
            if (text != NULL) {
                name = text->text();
            }
        }

        // Get value node and extract text content
        ptr = node->element("value", 0);
        if (ptr != NULL) {
            GXmlText* text = static_cast<GXmlText*>(ptr->child(0));
            if (text != NULL) {
                value = text->text();
            }
        }
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Post string content to Hub
 *
 * @param[in] content String content to post
 *
 * Posts the content of a string to the Hub.
 *
 * The method does nothing if no Hub connection has been established.
 ***************************************************************************/
void GVOClient::post_string(const std::string& content) const
{
    // Continue only if Hub connection has been established
    if (m_socket != -1) {

        // Determine content length
        int length = content.length();

        // Set prefix
        std::string prefix = "POST /xmlrpc HTTP/1.0\n"
                             "User-Agent: GammaLib\n"
                             "Content-Type: text/xml\n"
                             "Content-Length: "+str(length)+"\n\n";

        // Build post string
        std::string post = prefix + content;

        // Send content to socket
        bool done = false;
        do {
            int length      = post.length();
            int sent_length = send(m_socket, post.c_str(), length, 0);
            if (sent_length < length) {
                post = post.substr(sent_length, std::string::npos);
            }
            else {
                done = true;
            }
        } while (!done);

    } // endif: Hub connection had been established
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Receive string content from Hub
 *
 * @return String received from Hub.
 *
 * Reads information sent by Hub into a string.
 * 
 * The method does nothing if no Hub connection has been established.
 ***************************************************************************/
std::string GVOClient::receive_string(void) const
{
    // Initialise empty string
    std::string result = "";

    // Continue only if Hub connection has been established
    if (m_socket != -1) {

        // Define buffer
        char buffer[1001];

        // Read buffer until it is empty
        int n = 0;
        do {
            n = recv(m_socket, buffer, 1000, 0);
            if (n > 0) {
                buffer[n+1] = '\0';
                result.append(std::string(buffer));
            }
        } while (n > 0);

    } // endif: Hub connection had been established

    // Return result
    return result;
}


/***********************************************************************//**
 * @brief Returns SAMP Hub lockfile URL
 *
 * @return SAMP Hub lockfile URL (empty if no lockfile was found).
 *
 * Implements IVOA standard REC-SAMP-1.3-20120411.
 ***************************************************************************/
std::string GVOClient::get_hub_lockfile(void) const
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
