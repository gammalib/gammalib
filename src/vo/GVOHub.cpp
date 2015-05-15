/***************************************************************************
 *                     GVOHub.cpp - VO SAMP Hub class                      *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2014-2015 by Thierry Louge                               *
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
#include "GVOHub.hpp"
#include "GVOClient.hpp"
#include "GException.hpp"
#include "GTools.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_INIT_MEMBERS			                     "GVOHub::init_members()"
#define G_START_HUB                                     "GVOHub::start_hub()"
#define G_HANDLE_REQUEST                 "GVOHub::handle_request(socklen_t&)"
#define G_REGISTER_SERVICE      "GVOHub::register_service(GXml&, socklen_t&)"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */
#define G_KEY 12345

/* __ Debug definitions __________________________________________________ */
#define G_CONSOLE_DUMP


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

    // Start Hub
    start_hub();

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

    // Start Hub
    start_hub();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GVOHub::~GVOHub(void)
{
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

        // Start Hub
        start_hub();

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

    // Start Hub
    start_hub();

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
    m_name          = "GammaLib";
    m_secret        = random_string(15);
    m_hub_url       = "http://localhost:8001/xmlrpc";
    m_hub_host      = "127.0.0.1";
    m_hub_port      = "8001";
    m_version       = "1.3";
    m_hub_id        = "gammalib_hub";
    m_socket        = -1;        // Signals no socket
    m_cback_socket  = -1;
    m_nb_clients    = 0;     //  No registered clients at initialization
    m_clients       = NULL;
    m_shmem         = NULL;
    m_shmem_handler = 0;

    // Get shared memory for handling of clients
    if ((m_shmem_handler = shmget(G_KEY, sizeof(clients), 0666 | IPC_CREAT)) < 0) {
        std::string msg = "Could not open shared memory.";
        throw GException::invalid_value(G_INIT_MEMBERS, msg);
	}

    // Get pointer on shared memory
    if ((m_shmem = shmat(m_shmem_handler, NULL, 0)) == (void*) -1) {
        std::string msg = "Could not attach shared memory.";
        throw GException::invalid_value(G_INIT_MEMBERS, msg);
	}

    // Initialise shared memory
    m_clients             = static_cast<clients*>(m_shmem);
    m_clients->registered = 0;
    for (int i = 0; i < GVOHUB_NB_CLIENTS; ++i) {
        strcpy(m_clients->metadata[i].reference,    "NoRef");
        strcpy(m_clients->metadata[i].private_key,  "Unknown");
        strcpy(m_clients->metadata[i].name,         "Unknown");
        strcpy(m_clients->metadata[i].icon,         "Unknown");
        strcpy(m_clients->metadata[i].documentation,"Unknown");
        strcpy(m_clients->metadata[i].affiliation,  "Unknown");
        strcpy(m_clients->metadata[i].author_name,  "Unknown");
        strcpy(m_clients->metadata[i].email,        "Unknown");
        strcpy(m_clients->metadata[i].homepage,     "Unknown");
        strcpy(m_clients->metadata[i].port,         "Unknown");
        strcpy(m_clients->metadata[i].description,  "Unknown");
    }
    
    // Create SAMP file
    create_samp_file();

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
    m_name         = hub.m_name;
    m_secret       = hub.m_secret;
    m_hub_url      = hub.m_hub_url;
    m_hub_host     = hub.m_hub_host;
    m_hub_port     = hub.m_hub_port;
    m_version      = hub.m_version;
    m_hub_id       = hub.m_hub_id;
    m_socket       = hub.m_socket;
    m_cback_socket = hub.m_cback_socket;
    m_nb_clients   = hub.m_nb_clients;

    //TODO: What to do with shared memory? Copy over all data?

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GVOHub::free_members(void)
{ 

    // Remove lockfile
    std::string lockurl = get_hub_lockfile();
    std::remove(lockurl.c_str());

    // Release shared memory
    struct shmid_ds *buf;
    shmdt(m_shmem);
    shmctl(m_shmem_handler, IPC_RMID, buf);

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
    // Set maximum number of clients
    int max_clients = 10;

    // Prepare TCP/IP structure
    struct sockaddr_in serv_addr;
    std::memset(&serv_addr, 0, sizeof(serv_addr));
    serv_addr.sin_family      = AF_INET;
    serv_addr.sin_addr.s_addr = htonl(INADDR_ANY);
    serv_addr.sin_port        = htons(8001);
    
    // Initialise all client_socket[] to 0 so not checked
    socklen_t client_socket[max_clients];
    for (int i = 0; i < max_clients; ++i) {
        client_socket[i] = 0;
    }

    // Create Hub socket
    socklen_t hub_socket = socket(AF_INET, SOCK_STREAM, 0);
std::cout << "Created Hub socket " << hub_socket << std::endl;
    
    // Creation of hub main socket
    if (hub_socket < 0) {
        std::string msg = "Unable to create Hub socket. Errno="+
                          gammalib::str(errno);
        throw GException::runtime_error(G_START_HUB, msg);
    }
    
    // Set hub main socket to allow multiple connections
    int opt = 1;
    if (setsockopt(hub_socket, SOL_SOCKET, SO_REUSEADDR, (char*)&opt, sizeof(opt)) < 0) {
        std::string msg = "Unable to set Hub socket to multiple connections."
                          " Errno="+gammalib::str(errno);
        throw GException::runtime_error(G_START_HUB, msg);
    }
std::cout << "Allow for multiple connections on " << hub_socket << std::endl;
    
    // Server socket is opened. Now, bind it to the port, with family etc.
    if (bind(hub_socket, (struct sockaddr*)&serv_addr, sizeof(serv_addr)) < 0) {
        std::string msg = "Unable to bind Hub socket to server socket. Errno="+
                          gammalib::str(errno);
        throw GException::runtime_error(G_START_HUB, msg);
    }
std::cout << "Bind to socket " << hub_socket << std::endl;

    // Now start listening for the clients: 5 requests simultaneously pending
    // maximum
    if (listen(hub_socket, 5) < 0) {
        std::string msg = "Unable to start listening on Hub socket. Errno="+
                          gammalib::str(errno);
        throw GException::runtime_error(G_START_HUB, msg);
    }

    // ...
    struct sockaddr_in cli_addr;
    socklen_t clilen = sizeof(cli_addr);

    // Main event handling loop
    while (1) {
        
    	// Accept connection from the client 
    	socklen_t newsocket = accept(hub_socket, (struct sockaddr *)&cli_addr, &clilen);
    	if (newsocket < 0) {
            std::string msg = "Client connection to socket not accepted.";
            throw GException::runtime_error(G_START_HUB, msg);
    	}

    	// Create child process to handle the request
        int pid = fork();
        if (pid < 0) {
            std::string msg = "No child process created.";
            throw GException::runtime_error(G_START_HUB, msg);
        }

        // If we have a PID of 0 we are in the child process. In this case
        // we have to handle the incoming requests ...
        if (pid == 0) {
            handle_request(newsocket);
            close(newsocket);
            exit(0);
        }

        // ... otherwise we are in the parent process. Set SIG_IGN to child
	    // so that the kernel can close it without sending it to zombie state
        // <defunct>
        else {
            int status = 0;
            waitpid(pid, &status, 0);
            std::signal(pid, SIG_IGN);	
        }

    } // endwhile: main event loop
    
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

    // Initialise response
    std::string response = "";

    // Read from socket until nothing is received anymore.
    int timeout = 5000; // Initial timeout is 5 sec
    int n       = 0;
    do {
        n = gammalib::recv(sock, buffer, 1000, 0, timeout);
        if (n > 0) {
            buffer[n] = '\0';
            response.append(std::string(buffer));
        }
        timeout = 10; // The timeout now is 0.01 sec 
    } while (n > 0);
    
    // Dump the buffer
    #if defined(G_CONSOLE_DUMP)
    std::cout << "Hub has received the following message:" << std::endl;
    std::cout << response << std::endl;
    #endif
    
    // Extract response into an XML object
    GXml   xml;
    size_t start = response.find("<?xml");
    if (start != std::string::npos) {
        xml = GXml(response.substr(start, std::string::npos));
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

    // Dump the buffer
    #if defined(G_CONSOLE_DUMP)
    std::cout << "Method called: " << method_called << std::endl;
    #endif

    // Dispatch according to method
    if (method_called.compare("samp.hub.ping") == 0) {
        request_ping(sock);
    }
    else if (method_called.compare("samp.hub.register") == 0) {
        request_register(xml, sock);
    }
    else if (method_called.compare("samp.hub.unregister") == 0) {
        request_unregister(xml,sock);
    }
    else if (method_called.compare("samp.hub.declareMetadata") == 0) {
        request_declare_metadata(xml,sock);
    }
    else if (method_called.compare("samp.hub.declareSubscriptions") == 0) {
        request_declare_subscriptions(xml,sock);
    }
    else if (method_called.compare("samp.hub.setXmlrpcCallback") == 0) {
        request_set_xml_rpc_callback(xml,sock);
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
        request_ping(sock);
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
    
    // Declare message
    std::string msg = "";

    // Set response
    msg.append("<?xml version=\"1.0\"?>\n");
    msg.append("<methodResponse>\n");
    msg.append("<params>\n");
    msg.append("  <param><value/></param>\n");
    msg.append("</params>\n");
    msg.append("</methodResponse>\n");
    
    // Post response
    post_string(msg,sock);

    //Return
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

    // Initialize values
    char        cl_id[31];
    char        selfid[124];
    char        privatekey[124];
    std::string secret     = random_string(15);
    std::string translator = "http://localhost:8001/xmlrpc";
    
    // Search for a new free slot
    bool found = false;
    for (int i = 0; i < GVOHUB_NB_CLIENTS; ++i) {
        if (strcmp(m_clients->metadata[i].reference, "NoRef") == 0) {
            sprintf(cl_id, "c%d",i);
            strcpy(m_clients->metadata[i].reference, cl_id);
            sprintf(m_clients->metadata[i].private_key, "%s", secret.c_str());
            found = true;
            break;
        } 
    }

    // Throw an exception if no free slot is available
    if (!found) {
        std::string msg = "No free slot found for client.";
        throw GException::invalid_value(G_REGISTER_SERVICE, msg);
    }

    // Set self ID and private key
    sprintf(selfid,     "    <member><name>samp.self-id</name><value>%s</value></member>\n",cl_id);
    sprintf(privatekey, "    <member><name>samp.private-key</name><value>%s</value></member>\n",secret.c_str());

    // Declare response
    std::string msg = "";
    msg.append("<?xml version=\"1.0\" ?>\n");
    msg.append("<methodResponse>\n");
    msg.append("<params>\n");
    msg.append("  <param><value><struct>\n");
    msg.append("    <member><name>samp.hub-id</name><value>" + m_hub_id + "</value></member>\n");
    msg.append(selfid);
    msg.append(privatekey);
    msg.append("    <member><name>samp.status</name><value>samp.ok</value></member>\n");
    msg.append("    <member><name>samp.url-translator</name><value>" + translator + "</value></member>\n");
    msg.append("  </struct></value></param>\n");
    msg.append("</params>\n");
    msg.append("</methodResponse>\n");

    // Post response
    post_string(msg,sock);

    // Activate callbacks
    activate_callbacks("samp.hub.event.register", cl_id);	

    // Increment number of registered clients
    m_clients->registered ++;

    //Return
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

        // Declare message
        std::string msg = "";

        // Set response
        msg.append("<?xml version=\"1.0\"?>\n");
        msg.append("<methodResponse>\n");
        msg.append("<params>\n");
        msg.append("  <param><value/></param>\n");
        msg.append("</params>\n");
        msg.append("</methodResponse>\n");
            
        // Post response
        post_string(msg, sock);

        // Activate callbacks
        activate_callbacks("samp.hub.event.unregister",
                           m_clients->metadata[i].reference);

        // Reset client entry
        strcpy(m_clients->metadata[i].reference,      "NoRef");
        sprintf(m_clients->metadata[i].name,          "Unknown");
        sprintf(m_clients->metadata[i].description,   "Unknown");
        sprintf(m_clients->metadata[i].icon,          "Unknown");
        sprintf(m_clients->metadata[i].documentation, "Unknown");
        sprintf(m_clients->metadata[i].affiliation,   "Unknown");
        sprintf(m_clients->metadata[i].email,         "Unknown");
        sprintf(m_clients->metadata[i].author_name,   "Unknown");
        sprintf(m_clients->metadata[i].homepage,      "Unknown");
        
    } // endif: valid client found
    
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

        // Get metadata
        std::string client_name        = get_response_value(xml, "samp.name");
        std::string client_description = get_response_value(xml, "samp.description.text");
        std::string client_icon        = get_response_value(xml, "samp.icon.url");
        std::string client_doc         = get_response_value(xml, "samp.documentation.url");
        std::string client_affi        = get_response_value(xml, "author.affiliation");
        std::string client_mail        = get_response_value(xml, "author.email");
        std::string client_authname    = get_response_value(xml, "author.name");
        std::string client_page        = get_response_value(xml, "home.page");

        // Store metadata
        sprintf(m_clients->metadata[i].name,          "%s", client_name.c_str());
        sprintf(m_clients->metadata[i].description,   "%s", client_description.c_str());
        sprintf(m_clients->metadata[i].icon,          "%s", client_icon.c_str());
        sprintf(m_clients->metadata[i].documentation, "%s", client_doc.c_str());
        sprintf(m_clients->metadata[i].affiliation,   "%s", client_affi.c_str());
        sprintf(m_clients->metadata[i].email,         "%s", client_mail.c_str());
        sprintf(m_clients->metadata[i].author_name,   "%s", client_authname.c_str());
        sprintf(m_clients->metadata[i].homepage,      "%s", client_page.c_str());

        // Declare message
        std::string msg = "";

        // The hub response to Metadata declaration has no useful
        // content and can be discarded by the client.
        msg.append("<?xml version=\"1.0\"?>\n");
        msg.append("<methodResponse>\n");
        msg.append("<params><param><value/></param></params>\n");
        msg.append("</methodResponse>\n");

        // Post response
        post_string(msg,sock);
            
        // Activate callbacks
        activate_callbacks("samp.hub.event.metadata", m_clients->metadata[i].reference);

    } // endif: client found
        
    // Signal if the client was unknown
    #if defined(G_CONSOLE_DUMP)
    else {
        std::cout << " *** ERROR: Unable to find client " << get_client_key(xml) << std::endl;
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

        // Do stuff
        int j;
        std::list<std::string>::iterator p;
        std::ofstream myfile;
        std::string tempo;
        std::list<std::string> client_subscriptions = get_registrations(xml);
        for (j = 0, p=client_subscriptions.begin(); j<client_subscriptions.size()-1; j++,p++ ) {
            tempo = *p;
            sprintf(m_clients->metadata[i].registered_methods[j], "%s", tempo.c_str());
        }

        // Declare response
        std::string msg = "";
        msg.append("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
        msg.append("<methodResponse>\n");
        msg.append("<params>\n");
        msg.append("  <param><value><struct>\n");
        msg.append("    <member><name>samp.status</name><value>samp.ok</value></member>\n");
        msg.append("  </struct></value></param>\n");
        msg.append("</params>\n");
        msg.append("</methodResponse>\n");

        // Post string
        post_string(msg,sock);
        
    } // endif: valid client found
    
    //Return
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
void GVOHub::request_set_xml_rpc_callback(const GXml& xml,const socklen_t& sock)
{
    // Header
    #if defined(G_CONSOLE_DUMP)
    std::cout << "GVOHub::request_declare_subscriptions" << std::endl;
    #endif

    // Get client index
    int i = get_client_index(xml);

    // Continue only if index is valid
    if (i != -1) {

        // Get callback port
        std::string client_port = get_callback_port(xml);

        // Set callback port value
		sprintf(m_clients->metadata[i].port, "%s", client_port.c_str());

        // Declare message
        std::string msg = "";

        // The hub response to rpc callback declaration has no useful
        // content and can be discarded by the client.
        msg.append("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
        msg.append("<methodResponse>\n");
        msg.append("<params>\n");
        msg.append("  <param><value><struct>\n");
        msg.append("    <member><name>samp.status</name><value>samp.ok</value></member>\n");
        msg.append("  </struct></value></param>\n");
        msg.append("</params>\n");
        msg.append("</methodResponse>\n");

        // Post message
        post_string(msg, sock);

    } // endif: client found

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
void GVOHub::request_get_subscriptions(const GXml& xml,const socklen_t& sock)
{
    // Header
    #if defined(G_CONSOLE_DUMP)
    std::cout << "GVOHub::request_get_subscriptions" << std::endl;
    #endif

    // Get client index
    int i = get_client_index(xml);

    // Continue only if index is valid
    if (i != -1) {

        // Declare message
        std::string msg = "";

        // Set response
        msg.append("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
        msg.append("<methodResponse>\n");
        msg.append("<params>\n");
        msg.append("  <param><value><struct>\n");

        // Append all subscriptions
        char tempo[128];
		for (int j = 0; j < GVOHUB_NB_METHODS; ++j) {
			strcpy(tempo, m_clients->metadata[i].registered_methods[j]);
			if (strcmp("", tempo) != 0) {
				msg.append("    <member>\n");
                msg.append("      <name>");
				msg.append(tempo);
				msg.append("</name>\n");
                msg.append("      <value></value>\n)");
                msg.append("    </member>\n");
			}
		}

        // Finish response
        msg.append("  </struct></value></param>\n");
        msg.append("</params>\n");
        msg.append("</methodResponse>\n");
        
        // Post response
        post_string(msg,sock);

    } // endif: client index was valid
   
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
    msg.append("<?xml version=\'1.0\' encoding=\"UTF-8\"?>\n");
    msg.append("<methodResponse>\n");
    msg.append("<params>\n");
    msg.append("  <param><value><array><data>\n");
    msg.append("    <value>gammalib_hub</value>\n");

    // Loop over all clients. Do not send back current client's registration
    for (int i = 0; i < GVOHUB_NB_CLIENTS; ++i) {
        if (strcmp(m_clients->metadata[i].reference,"NoRef") != 0) {
            if (strcmp(m_clients->metadata[i].private_key, key.c_str()) != 0) {
                msg.append("    <value>");
                msg.append(m_clients->metadata[i].reference);
                msg.append("</value>\n");
            }
        }
    }
    
    // Finish response
    msg.append("  </data></array></value></param>\n");
    msg.append("</params>\n");
    msg.append("</methodResponse>\n");

    // Post response
    post_string(msg, sock);
    
    //Return
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

/*
    std::string buffer;
    const int n = 1000; 
    char      line[n];
    FILE* fptr = fopen("SAMP-registered.txt", "r");
    if (fptr != NULL) {
        while (fgets(line, n, fptr) != NULL) {

                // Convert line to C++ string
                std::string cline = std::string(line);

                // Check for secret key
                if (cline.compare(0, 10, "Registered") == 0) {
                    buffer = gammalib::strip_chars(cline.substr(10, std::string::npos), "\r\n");
                }
    }
    
    // Closes file
    fclose(fptr);
*/

    // Declare message
    std::string msg = "";

    // Set response
    msg.append("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
    msg.append("<methodResponse>\n");
    msg.append("<params>\n");
    msg.append("  <param><value/></param>\n");
    msg.append("</params>\n");
    msg.append("</methodResponse>\n");

    // Post message
    post_string(msg,sock);
    
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
void GVOHub::request_get_metadata(const GXml& xml,const socklen_t& sock)
{
    // Header
    #if defined(G_CONSOLE_DUMP)
    std::cout << "GVOHub::request_get_metadata" << std::endl;
    #endif

    // Get client key
    std::string key = get_client_key(xml);

    // Declare message
    std::string msg = "";

    // Set response
    msg.append("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
    msg.append("<methodResponse>\n");
    msg.append("<params>\n");
    msg.append("  <param><value><struct>\n");
    msg.append("    <member>\n");
    msg.append("      <name>samp.status</name>\n");
    msg.append("      <value>samp.ok</value>\n");
    msg.append("    </member>\n");

    // If key is "gammalib_hub" then return some special information
    if (key == "gammalib_hub") {
		msg.append("    <member><name>samp.name</name><value>gammalib_hub</value></member>\n");
   		msg.append("    <member><name>samp.description.text</name><value>GammaLib VO Hub</value></member>\n");
		msg.append("    <member><name>samp.icon.url</name><value></value></member>\n");
		msg.append("    <member><name>samp.documentation.url</name><value>gammalib github documentation</value></member>\n");
		msg.append("    <member><name>author.affiliation</name><value>IRAP, CNRS Toulouse</value></member>\n");
		msg.append("    <member><name>author.email</name><value>dontcontactme@please</value></member>\n");
		msg.append("    <member><name>author.name</name><value>J. Knoedlseder, T. Louge</value></member>\n");
		msg.append("    <member><name>home.page</name><value>no dedicated page</value></member>\n");
    }

    // ... otherwise get the index
    else {
 
        // Get client index
        int i = get_client_index(xml);

        // Continue only if index is valid
        if (i != -1) {

            // Append response
            msg.append("    <member><name>samp.name</name><value>");
            msg.append(m_clients->metadata[i].name);
            msg.append("</value></member>\n");
            msg.append("    <member><name>samp.description.text</name><value>");
            msg.append(m_clients->metadata[i].description);
            msg.append("</value></member>\n");
            msg.append("    <member><name>samp.icon.url</name><value>");
            msg.append(m_clients->metadata[i].icon);
            msg.append("</value></member>\n");
            msg.append("    <member><name>samp.documentation.url</name><value>");
            msg.append(m_clients->metadata[i].documentation);
            msg.append("</value></member>\n");
            msg.append("    <member><name>author.affiliation</name><value>");
            msg.append(m_clients->metadata[i].affiliation);
            msg.append("</value></member>\n");
            msg.append("    <member><name>author.email</name><value>");
            msg.append(m_clients->metadata[i].email);
            msg.append("</value></member>\n");
            msg.append("    <member><name>author.name</name><value>");
            msg.append(m_clients->metadata[i].author_name);
            msg.append("</value></member>\n");
            msg.append("    <member><name>home.page</name><value>");
            msg.append(m_clients->metadata[i].homepage);
            msg.append("</value></member>\n");

        } // endif: index was valid
        
    } // endelse: client was not the hub

    // Finish response
    msg.append("  </struct></value></param>\n");
    msg.append("</params>\n");
    msg.append("</methodResponse>\n");
    
    // Post response
    post_string(msg, sock);
   
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
 *
 * Extracts the client index in shared memory from the XML request. The
 * method returns -1 if no client was found.
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

        // Searches for key in metadata
        for (int i = 0; i < GVOHUB_NB_CLIENTS; ++i) {
            if (strcmp(m_clients->metadata[i].private_key, key.c_str()) == 0) {
                index = i;
                break;
            }
        }

    } // endif: key was not empty

    // Return index
    return index;
}


/***********************************************************************//**
 * @brief Returns value for a SAMP client query parameter
 *
 * @param[in] xml client query XML document.
 * @param[in] name Parameter name.
 * @return Parameter value.
 *
 * Returns value for a SAMP client query parameter. If the specified
 * parameter was not found or if the response structure is not compliant,
 * an empty string is returned.
 ***************************************************************************/
std::string GVOHub::get_response_value(const GXml&        xml,
                                       const std::string& name) const
{
    // Declare value
    std::string value = "";

    // Get the client's private key
    const GXmlNode* node = xml.element("methodCall > params > param > value > struct");
    if (node != NULL) {
        int num = node->elements("member");            
        for (int i = 0; i < num; ++i) {
            const GXmlNode* member = node->element("member", i);
            std::string one_name;
            std::string one_value;
            get_name_value_pair(member, one_name, one_value);
            if (one_name == name) {
                value = one_value;
                break;
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
void GVOHub::get_name_value_pair(const GXmlNode* node,
                                 std::string&    name,
                                 std::string&    value) const
{
    // Clear name and value strings
    name.clear();
    value.clear();

    // Continue only if node is valid
    if (node != NULL) {

        // Get name node and extract text content
        const GXmlNode* ptr = node->element("name", 0);
        if (ptr != NULL) {
            const GXmlText* text = static_cast<const GXmlText*>((*ptr)[0]);
            if (text != NULL) {
                name = text->text();
            }
        }

        // Get value node and extract text content
        ptr = node->element("value", 0);
        if (ptr != NULL) {
            const GXmlText* text = static_cast<const GXmlText*>((*ptr)[0]);
            if (text != NULL) {
                value = text->text();
            }
        }
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Returns registrations from XML document
 *
 * @param[in] xml client query XML document.
 * @return List of registrations.
 ***************************************************************************/
std::list<std::string> GVOHub::get_registrations(const GXml& xml) const                                         
{
    // Declare value
    std::list<std::string> value;

    // Get the client's private key
    const GXmlNode* node = xml.element("methodCall > params > param > value > struct");
    if (node != NULL) {
        int num = node->elements("member");
        for (int i = 0; i < num; ++i) {
            const GXmlNode* member = node->element("member", i);
            if (member != NULL) {
                const GXmlNode* name = member->element("name", 0);
                if (name != NULL) {
                    const GXmlText* text = static_cast<const GXmlText*>((*name)[0]);
                    if (text != NULL) {
                        value.push_front(text->text());
                    }
                }
            }
        }
    }

    // Return value
    return value;
}


/***********************************************************************//**
 * @brief Returns callback port
 *
 * @param[in] xml client query XML document.
 * @return Callback port value.
 ***************************************************************************/
std::string GVOHub::get_callback_port(const GXml& xml) const                                         
{
    // Declare value
    std::string value;

    // Get the client's private key
    const GXmlNode* node = xml.element("methodCall > params > param[1] > value");
    if (node != NULL) {
		const GXmlText* text = static_cast<const GXmlText*>((*node)[0]);
        if (text != NULL) {
		    value.append(text->text());
        }
    }

    // Return value
    return value;
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
 * @brief Test if any registered client needs to know that the method 
 *		  has been called
 *
 * @param[in] method: name of the method
 * @param[in] cl_id id of calling client
 ***************************************************************************/
void GVOHub::activate_callbacks(std::string method,char cl_id[31])
{
    int k;
    char port[16];
    char pkey[16];
    std::string spkey;

    // Loop over all clients
    for (int i = 0; i < GVOHUB_NB_CLIENTS; ++i) {
    
        // Loop over all methods
        for (int j = 0; j < GVOHUB_NB_METHODS; ++j) {
    
            // Continue only if method is valid
            if (strcmp("",m_clients->metadata[i].registered_methods[j]) != 0) {

                //
                if (strcmp(method.c_str(), m_clients->metadata[i].registered_methods[j]) == 0) {
                    strcpy(port,m_clients->metadata[i].port);
                    std::string port_string = std::string(port);
                    //Eliminate http://127.0.0.1/ from port string -----> improve this code is necessary
                    port_string.replace(0,17,"");
                    std::size_t pos = port_string.find("/");
                    std::string port_only = port_string.substr(0,pos);

                    // Set hints
                    struct addrinfo hints;
                    std::memset(&hints, 0, sizeof(hints));
                    hints.ai_family   = AF_INET;
                    hints.ai_socktype = SOCK_STREAM;

                    // Get server information
                    struct addrinfo* servinfo;
                    if (getaddrinfo("127.0.0.1", port_only.c_str(),
						&hints, &servinfo) == 0) {

                        // Loop through all the results and connect to the first
                        // we can
                        for (struct addrinfo* ptr = servinfo; ptr != NULL; ptr = ptr->ai_next) {

                            // Create socket
                            m_cback_socket = socket(ptr->ai_family,
                                                    ptr->ai_socktype,
                                                    ptr->ai_protocol);

                            // Connect to socket if socket is valid
                            if (m_cback_socket != -1) {
                                if (::connect(m_cback_socket,
                                              ptr->ai_addr,
                                              ptr->ai_addrlen) == -1) {
                                    close(m_cback_socket);
                                    m_cback_socket = -1;
                                }
                            }

                        } // endfor: looped through all results

                    } // endif: server information was valid

                    //
                    std::string msg = "";

                    //
                    if (strcmp(method.c_str(),"samp.hub.event.register") == 0) {

                        // Set methodResponse elts for samp.hub.event.register
                        msg.append("<?xml version=\"1.0\"?>\n");
                        msg.append("<methodCall><methodName>samp.client.receiveNotification</methodName>\n");
                        msg.append("\t<params>\n");
                        msg.append("\t\t<param><value>");
                        msg.append(m_clients->metadata[i].private_key);
                        msg.append("</value></param>\n");
                        msg.append("\t\t<param><value>"+m_hub_id+"</value></param>\n");
                        msg.append("\t\t<param>\n\t\t\t<value>\n\t\t\t\t<struct>\n");
                        msg.append("\t\t\t\t\t<member><name>samp.mtype</name><value>samp.hub.event.register</value></member>\n");
                        msg.append("\t\t\t\t\t<member><name>samp.params</name><value>\n\t\t\t\t\t\t<struct>\n\t\t\t\t\t\t\t<member>");
                        msg.append("<name>id</name><value>");
                        msg.append(cl_id);
                        msg.append("</value></member>\n\t\t\t\t\t\t</struct>\n\t\t\t\t\t</value></member>\n");
                        msg.append("\t\t\t\t</struct>\n\t\t\t</value>\n\t\t</param>\n");
                        msg.append("\t</params>\n");
                        msg.append("</methodCall>\n");
                    
                        // Post response
                        post_string_callback(msg);
                    
                    }
                    
                    //
                    else if (strcmp(method.c_str(),"samp.hub.event.unregister") == 0) {

                        // Set methodResponse elts for samp.hub.event.unregister
                        msg.append("<?xml version=\"1.0\"?>\n");
                        msg.append("<methodCall><methodName>samp.client.receiveNotification</methodName>\n");
                        msg.append("\t<params>\n");
                        msg.append("\t\t<param><value>");
                        msg.append(m_clients->metadata[i].private_key);
                        msg.append("</value></param>\n");
                        msg.append("\t\t<param><value>"+m_hub_id+"</value></param>\n");
                        msg.append("\t\t<param>\n\t\t\t<value>\n\t\t\t\t<struct>\n");
	    				msg.append("\t\t\t\t\t<member><name>samp.mtype</name><value>samp.hub.event.unregister</value></member>\n");
                        msg.append("\t\t\t\t\t<member><name>samp.params</name><value>\n\t\t\t\t\t\t<struct>\n\t\t\t\t\t\t\t<member>");
                        msg.append("<name>id</name><value>");
                        msg.append(cl_id);
                        msg.append("</value></member>\n\t\t\t\t\t\t</struct>\n\t\t\t\t\t</value></member>\n");
	    				msg.append("\t\t\t\t</struct>\n\t\t\t</value>\n\t\t</param>\n");
                        msg.append("\t</params>\n");
                        msg.append("</methodCall>\n");
                        
                        // Post response
                        post_string_callback(msg);
					
                    }
                    
                    //
			        else if (strcmp(method.c_str(),"samp.hub.event.metadata") == 0) {

                        // Set methodResponse elts for samp.hub.event.unregister
                        msg.append("<?xml version=\"1.0\"?>\n");
                        msg.append("<methodCall><methodName>samp.client.receiveNotification</methodName>\n");
                        msg.append("\t<params>\n");
                        msg.append("\t\t<param><value>");
                        msg.append(m_clients->metadata[i].private_key);
                        msg.append("</value></param>\n");
                        msg.append("\t\t<param><value>"+m_hub_id+"</value></param>\n");
                        msg.append("\t\t<param>\n\t\t\t<value>\n\t\t\t\t<struct>\n");
	    				msg.append("\t\t\t\t\t<member><name>samp.mtype</name><value>samp.hub.event.metadata</value></member>\n");
                        msg.append("\t\t\t\t\t<member><name>samp.params</name><value>\n\t\t\t\t\t\t<struct>\n\t\t\t\t\t\t\t");
                        msg.append("\t\t\t\t\t<member><name>id</name><value>");
                        for (k = 0; k < GVOHUB_NB_CLIENTS; k++) {
                            if (strcmp(cl_id,m_clients->metadata[k].reference) == 0) {
                                break;					
                            }
                        }
                        msg.append(m_clients->metadata[k].reference);
                        msg.append("</value>\n\t\t\t\t\t</member>\n");
                        msg.append("\t\t\t\t\t<member><name>metadata</name><value>\n\t\t\t\t\t\t<struct>\n\t\t\t\t\t\t\t");
                        msg.append("\t\t\t\t\t<member>\n\t\t\t\t\t\t<name>samp.name</name>\n\t\t\t\t\t\t<value>");
                        msg.append(m_clients->metadata[k].name);
                        msg.append("</value>\n\t\t\t\t\t</member>\n");
                        msg.append("\t\t\t\t\t<member>\n\t\t\t\t\t\t<name>samp.description.text</name>\n\t\t\t\t\t\t<value>");
                        msg.append(m_clients->metadata[k].description);
                        msg.append("</value>\n\t\t\t\t\t</member>\n");
                        msg.append("\t\t\t\t\t<member>\n\t\t\t\t\t\t<name>samp.icon.url</name>\n\t\t\t\t\t\t<value>");
                        msg.append(m_clients->metadata[k].icon);
                        msg.append("</value>\n\t\t\t\t\t</member>\n");
                        msg.append("\t\t\t\t\t<member>\n\t\t\t\t\t\t<name>samp.documentation.url</name>\n\t\t\t\t\t\t<value>");
                        msg.append(m_clients->metadata[k].documentation);
                        msg.append("</value>\n\t\t\t\t\t</member>\n");
                        msg.append("\t\t\t\t\t<member>\n\t\t\t\t\t\t<name>author.affiliation</name>\n\t\t\t\t\t\t<value>");
                        msg.append(m_clients->metadata[k].affiliation);
                        msg.append("</value>\n\t\t\t\t\t</member>\n");
                        msg.append("\t\t\t\t\t<member>\n\t\t\t\t\t\t<name>author.email</name>\n\t\t\t\t\t\t<value>");
                        msg.append(m_clients->metadata[k].email);
                        msg.append("</value>\n\t\t\t\t\t</member>\n");
                        msg.append("\t\t\t\t\t<member>\n\t\t\t\t\t\t<name>author.name</name>\n\t\t\t\t\t\t<value>");
                        msg.append(m_clients->metadata[k].author_name);
                        msg.append("</value>\n\t\t\t\t\t</member>\n");
                        msg.append("\t\t\t\t\t<member>\n\t\t\t\t\t\t<name>home.page</name>\n\t\t\t\t\t\t<value>");
                        msg.append(m_clients->metadata[k].homepage);
                        msg.append("</value>\n\t\t\t\t\t</member>\n");	
                        msg.append("</struct></value></member>\n");				
                        msg.append("</struct>\n\t\t\t\t\t</value></member>\n");
	    				msg.append("\t\t\t\t</struct>\n\t\t\t</value>\n\t\t</param>\n");
                        msg.append("\t</params>\n");
                        msg.append("</methodCall>\n");
                        
                        // Post response
                        post_string_callback(msg);
					
                    }
                }
            } // endif: method was valid
        } // endfor: looped over all methods
    } // endfor: looped over all clients
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Creates the lockfile, fill it
 *
 * Implements IVOA standard REC-SAMP-1.3-20120411.
 ***************************************************************************/
void GVOHub::create_samp_file(void)
{
    // Get lockfile URL
    std::string lockurl = get_hub_lockfile();

    // Open SAMP lockfile. Continue only if opening was successful
    FILE* fptr = fopen(lockurl.c_str(), "w");
    if (fptr != NULL) {

        // If successful writes basic hub configuration to the lockfile
        fputs("# SAMP lockfile \n# Required keys:\n",fptr);
        char sampsecret[256] = "samp.secret=";
        strcat(sampsecret,m_secret.c_str());
        strcat(sampsecret,"\n");
        fputs (sampsecret,fptr);
        char sampurl[256] = "samp.hub.xmlrpc.url=";
        strcat(sampurl,m_hub_url.c_str());
        strcat(sampurl,"\n");
        fputs(sampurl,fptr);
        char samprofile[256] = "samp.profile.version=";
        strcat(samprofile,m_version.c_str());
        strcat(samprofile,"\n");
        fputs(samprofile,fptr);
        fputs("# Info stored by hub for some private reason:\n",fptr);
        char sampid[256] = "gammalib.hubid=";
        strcat(sampid,m_hub_id.c_str());
        strcat(sampid,"\n");
        fputs(sampid,fptr);
    }
    
    // Close SAMP lockfile
    fclose(fptr);

    // Return
    return;
}



/***********************************************************************//**
 * @brief Post string content to client
 *
 * @param[in] content String content to post.
 * @param[in] sock Socket.
 *
 * Posts the content of a string to a client.
 ***************************************************************************/
void GVOHub::post_string(const std::string& content,const socklen_t& sock) const
{
    // Continue only if socket is valid
    if (sock != -1) {

        // Determine content length
        int length = content.length();

        // Set prefix
        std::string prefix = "HTTP/1.1 200 OK\n"
                             "Connection: close\n"
                             "Content-Type: text/xml\n"
                             "Content-Length: "+gammalib::str(length)+"\n\n";

        // Build post string
        std::string post = prefix + content;

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
 * @brief Post string content to client for callback
 *
 * @param[in] content String content to post
 *
 * Posts the content of a string to the client for a callback.
 *
 * The method does nothing if no callback socket is open.
 ***************************************************************************/
void GVOHub::post_string_callback(const std::string& content) const
{
    // Continue only if Hub connection has been established
    if (m_cback_socket != -1) {

        // Determine content length
        int length = content.length();

        // Set prefix
        std::string prefix = "POST /xmlrpc HTTP/1.1\n"
                             "Connection: keep-alive\n"
                             "User-Agent: GammaLib\n"
                             "Host: localhost:"+gammalib::str(m_cback_socket)+"\n"
                             "Content-Type: text/xml\n"
                             "Content-Length: "+gammalib::str(length)+"\n\n";

        // Build post string
        std::string post = prefix + content;

        // Send content to socket
        bool done = false;
        do {
            int length      = post.length();
            int sent_length = send(m_cback_socket, post.c_str(), length, 0);
            if (sent_length < length) {
                post = post.substr(sent_length, std::string::npos);
            }
            else {
                done = true;
            }
        } while (!done);

        // Read buffer until it is empty
        char buffer[1001];
        int n = 0;
        do {
            n = gammalib::recv(m_cback_socket, buffer, 1000, 0, 100);
            if (n > 0) {
                buffer[n+1] = '\0';
            }
        } while (n > 0);

    } // endif: Hub connection had been established
	
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
    std::string str = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz";
    int pos;
    while(str.size() != length) {
        pos = ((rand() % (str.size() - 1)));
        str.erase (pos, 1);
    }
   
    // Return string
    return str;
}

