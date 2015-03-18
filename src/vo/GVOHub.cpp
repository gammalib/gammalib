/***************************************************************************
 *                     GVOHub.hpp - VO SAMP Hub class                      *
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
#include <errno.h>
#endif
#include <cstdlib>         // std::getenv() function
#include <cstdio>          // std::fopen(), etc. functions
#include <cstring>         // std::memset() function
#include <unistd.h>        // close() function
#include <netdb.h>         // getaddrinfo() function
#include <sys/socket.h>    // socket(), connect() functions
#include <fstream>
#include "GVOHub.hpp"
#include "GVOClient.hpp"
#include "GException.hpp"
#include "GTools.hpp"
#include "GVOApp.hpp"

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
void GVOHub::ping_service(const socklen_t& sock)
{
    printf("ping call received\n");
    // Declare message
    std::string msg = "";

    // Set methodResponse elts
    msg.append("<?xml version=\"1.0\"?>\n");
    msg.append("<methodResponse>\n");
    msg.append("</methodResponse>\n");
    int length = msg.length();
    
    int n = write(sock,msg.c_str(),length);
    if (n < 0) {
	printf("Could not write socket buffer.\n");
        std::string msg = "Could not write socket buffer.";
        throw GException::invalid_value(G_HANDLE_REQUEST, msg);
    }

    printf("Response message sent\n");
    //Return
    return;
}

/***********************************************************************//**
 * @brief 
 *
 *
 ***************************************************************************/
void GVOHub::register_service(const GXml& xml, const socklen_t& sock)
{
    
    printf("Registration request called\n");
    
    // Declare message
    std::string msg = "";
    char buffer[124];
    char privatekey[124];
    char cl_id[256];
    int i;
    std::string content;
    std::string translator = "http://localhost:8001/";
    printf("Message composition\n");
    
    FILE* fptr = fopen("SAMP-registered.txt", "w");
    if (fptr != NULL) {

        // If successfull writes application definition
        fputs ("# SAMP registered applications and Methods \n",fptr);
	sprintf(cl_id,"Registered client-key:%d\n", m_nb_clients);
        fputs (cl_id,fptr);
    }
    
    // Closes file
    fclose(fptr);

    printf("Message composition 2 \n");
    sprintf(buffer,"<member><name>samp.self-id</name><value><string>client-id:%d</string></value></member>\n",m_nb_clients);
    sprintf(privatekey,"<member><name>samp.private-key</name><value><string>client-key:%d</string></value></member>\n",m_nb_clients);
    printf(buffer);
    m_nb_clients ++;
    
    // Set methodResponse elts
    msg.append("<?xml version=\"1.0\"?>\n");
    msg.append("<methodResponse>\n");
    msg.append("\t<params>\n");
    msg.append("\t\t<param>\n\t\t\t<value>\n\t\t\t\t<struct>\n");
    //msg.append("<member><name>samp.private-key</name><value><string>client-key:" + m_hub_id + "</string></value></member>\n");
    msg.append(privatekey);
    msg.append("<member><name>samp.hub-id</name><value><string>client-id:" + m_hub_id + "</string></value></member>\n");
    msg.append(buffer);
    msg.append("<member><name>samp.url-translator</name><value><string>" + translator + "</string></value></member>\n");
    msg.append("</struct></value></param>\n");
    msg.append("</params>\n");
    msg.append("</methodResponse>\n");

    int length = msg.length();
    printf("Sending message %s\n",msg.c_str());
    fflush(stdout);
    int n = write(sock,msg.c_str(),length);
    if (n < 0) {
	printf("Could not write socket buffer.\n");
        std::string msg = "Could not write socket buffer.";
        throw GException::invalid_value(G_HANDLE_REQUEST, msg);
    }

    printf("Response message sent\n");
    //Return
    return;
}


/***********************************************************************//**
 * @brief 
 *
 *
 ***************************************************************************/
void GVOHub::unregister(const socklen_t& sock)
{
    /* Code to create new file with replacement to erase registered app. File output should be renamed to the original file name
    string strReplace = "HELLO";
    string strNew = "GOODBYE";
    ifstream filein("filein.txt"); //File to read from
    ofstream fileout("fileout.txt"); //Temporary file
    if(!filein || !fileout)
    {
        cout << "Error opening files!" << endl;
        return 1;
    }

    string strTemp;
    //bool found = false;
    while(filein >> strTemp)
    {
        if(strTemp == strReplace){
            strTemp = strNew;
            //found = true;
        }
        strTemp += "\n";
        fileout << strTemp;
        //if(found) break;
    }*/
    close(sock);
    return;
}

/***********************************************************************//**
 * @brief 
 *
 *
 ***************************************************************************/
void GVOHub::register_metadata(const GXml& xml,const socklen_t& sock)
{
    char buffer[124];
    printf("Client Metadata received\n");
    // Search for metadata values to store
    std::string client_calling = "";
    // Search for value of methodName
    const GXmlNode* node = xml.element("methodCall", 0);
    node = node->element("params", 0);
    node = node->element("param", 0);
    node = node->element("value", 0);
    if (node != NULL) {
        const GXmlText* text = static_cast<const GXmlText*>((*node)[0]);
        if (text != NULL) {
            client_calling = text->text();
        }
    }
    std::cout << client_calling << std::endl;
    std::string client_name = get_response_value(xml, "samp.name");
    sprintf(buffer,"\n%s metadata : Name %s\n",client_calling.c_str(),client_name.c_str());
    
    printf(buffer);
    fflush(stdout);
    FILE* fptr = fopen("SAMP-registered.txt", "a");
    if (fptr != NULL) {
        fputs (buffer,fptr);
    }
    
    // Closes file
    fclose(fptr);
    // Declare message
    std::string msg = "";

    //The hub response to Metadata declaration has no usefull content and can be discarded by the client.
    // Set methodResponse elts
    msg.append("<?xml version=\"1.0\"?>\n");
    msg.append("<methodResponse>\n");
    msg.append("</methodResponse>\n");

    int length = msg.length();
    
    int n = write(sock,msg.c_str(),length);
    if (n < 0) {
	printf("Could not write socket buffer.\n");
        std::string msg = "Could not write socket buffer.";
        throw GException::invalid_value(G_HANDLE_REQUEST, msg);
    }

    printf("Response message sent\n");
    //Return
    return;
}

/***********************************************************************//**
 * @brief 
 *
 *
 ***************************************************************************/
void GVOHub::set_xml_rpc_callback(const GXml& xml,const socklen_t& sock)
{
    char buffer[124];
    printf("Client callback registration\n");
    // Search for metadata values to store
    std::string client_subscriptions = get_callback_port(xml, "samp.hub.setXmlrpcCallback");
    std::cout << client_subscriptions << std::endl;
    
    FILE* fptr = fopen("SAMP-registered.txt", "a");
    if (fptr != NULL) {
        fputs (client_subscriptions.c_str(),fptr);
    }
    
    // Closes file
    fclose(fptr);
    // Declare message
    std::string msg = "";

    //The hub response to rpc callback declaration has no usefull content and can be discarded by the client.
    // Set methodResponse elts
    msg.append("<?xml version=\"1.0\"?>\n");
    msg.append("<methodResponse>\n");
    msg.append("</methodResponse>\n");

    int length = msg.length();
    
    int n = write(sock,msg.c_str(),length);
    if (n < 0) {
	printf("Could not write socket buffer.\n");
        std::string msg = "Could not write socket buffer.";
        throw GException::invalid_value(G_HANDLE_REQUEST, msg);
    }

    printf("Response message sent\n");
    //Return
    return;
}
/***********************************************************************//**
 * @brief 
 *
 *
 ***************************************************************************/
void GVOHub::get_registered_clients(const GXml& xml,const socklen_t& sock)
{
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
    // Declare message
    std::string msg = "";

    //The hub response to rpc callback declaration has no usefull content and can be discarded by the client.
    // Set methodResponse elts
    msg.append("<?xml version=\"1.0\"?>\n");
    msg.append("<methodResponse>\n");
    msg.append("</methodResponse>\n");

    int length = msg.length();
    
    int n = write(sock,msg.c_str(),length);
    if (n < 0) {
	printf("Could not write socket buffer.\n");
        std::string msg = "Could not write socket buffer.";
        throw GException::invalid_value(G_HANDLE_REQUEST, msg);
    }

    printf("Response message sent\n");
    //Return
    return;
  }
}
/***********************************************************************//**
 * @brief 
 *
 *
 ***************************************************************************/
void GVOHub::declare_subscriptions(const GXml& xml,const socklen_t& sock)
{
    std::list<std::string>::iterator p;
    int i;
    std::ofstream myfile;
    std::string content;
    printf("Declare Subscriptions\n");
    
    std::string client_calling = "";
    // Search for value of methodName
    const GXmlNode* node = xml.element("methodCall", 0);
    node = node->element("params", 0);
    node = node->element("param", 0);
    node = node->element("value", 0);
    if (node != NULL) {
        const GXmlText* text = static_cast<const GXmlText*>((*node)[0]);
        if (text != NULL) {
            client_calling = text->text();
        }
    }
    std::cout << client_calling << std::endl;
    
    std::list<std::string> client_subscriptions = get_registrations(xml, "samp.hub.declareSubscriptions");
    std::string temp;
    temp.append("\n");
    FILE* fptr = fopen("SAMP-registered.txt", "a");
    if (fptr != NULL) {
      for( i=0,p=client_subscriptions.begin(); i<client_subscriptions.size()-1; i++,p++ ) {
	 std::cout << "Client registered to following methods calls:" << std::endl;
         std::cout << *p << std::endl;
	 temp.append(*p);
	 temp.append(" registered for ");
	 temp.append(client_calling);
	 temp.append("\n");
      }
      fputs (temp.c_str(),fptr);
    }
    
    // Declare message
    std::string msg = "";
    // Set methodResponse elts
    msg.append("<?xml version=\"1.0\"?>\n");
    msg.append("<methodResponse>\n");
    msg.append("</methodResponse>\n");

    int length = msg.length();
    
    int n = write(sock,msg.c_str(),length);
    if (n < 0) {
	printf("Could not write socket buffer.\n");
        std::string msg = "Could not write socket buffer.";
        throw GException::invalid_value(G_HANDLE_REQUEST, msg);
    }

    printf("Response message sent\n");
    //Return
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
                //result.append("registered as \""+m_client_id);
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
    m_name       = "GammaLib";
    m_secret     = "mysupersecret#0032557sentence";
    m_hub_url    = "http://localhost:8001";
    m_hub_host   = "127.0.0.1";
    m_hub_port   = "8001";
    m_version    = "1.3";
    m_client_key.clear();
    m_hub_id     = "b79884e0";
    //m_client_id.clear();
    m_socket     = -1;        // Signals no socket
    m_nb_clients = 0;	      //  No registered clients at initialization
    connected().clear();
    
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
    m_name           = hub.m_name;
    m_secret         = hub.m_secret;
    m_hub_url        = hub.m_hub_url;
    m_hub_host       = hub.m_hub_host;
    m_hub_port       = hub.m_hub_port;
    m_version        = hub.m_version;
    m_client_key     = hub.m_client_key;
    m_hub_id         = hub.m_hub_id;
    //m_client_id      = hub.m_client_id;
    m_socket         = hub.m_socket;
    connected() = hub.connected();

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
    char buffer[4096];
    bzero(buffer,4096);
   
    // Declare empty XML document
    GXml xml;
    
    // Read buffer
    int n = read(sock,buffer,4096);
    if (n < 0) {
        std::string msg = "Could not read socket buffer.";
        throw GException::invalid_value(G_HANDLE_REQUEST, msg);
    }
    printf("Here is the message: %s\n",buffer);
    std::string response = buffer;
	
    // Find start of XML text
    size_t start = response.find("<?xml");
    // If found then convert text into XML document
    if (start != std::string::npos) {
        xml = GXml(response.substr(start, std::string::npos));
    } else {
      //No XML tag found.
      return;
      
    }
    
    std::string method_called = "";
    // Search for value of methodName
    const GXmlNode* node = xml.element("methodCall", 0);
    const GXmlNode* ptr = node->element("methodName", 0);
    if (ptr != NULL) {
        const GXmlText* text = static_cast<const GXmlText*>((*ptr)[0]);
        if (text != NULL) {
            method_called = text->text();
        }
    }
    // endif: connection has been established
    //std::string method_called = get_response_value(xml, "methodName");
    if (method_called.compare("samp.hub.register") == 0) {
        register_service(xml, sock);
    }
    if (method_called.compare("samp.hub.unregister") == 0) {
        unregister(sock);
    }
    if (method_called.compare("samp.hub.declareMetadata") == 0) {
        register_metadata(xml,sock);
    }
    if (method_called.compare("samp.hub.ping") == 0) {
        ping_service(sock);
    }
    if (method_called.compare("samp.hub.declareSubscriptions") == 0) {
        declare_subscriptions(xml,sock);
    }
    if (method_called.compare("samp.hub.getSubscriptions") == 0) {
        ping_service(sock);
    }
    if (method_called.compare("samp.hub.setXmlrpcCallback") == 0) {
        set_xml_rpc_callback(xml,sock);
    }
    if (method_called.compare("samp.hub.getRegisteredClients") == 0) {
        get_registered_clients(xml, sock);
    }
    if (method_called.compare("samp.hub.getSubscribedClients") == 0) {
        ping_service(sock);
    }
    if (method_called.compare("samp.hub.notify") == 0) {
        ping_service(sock);
    }
    if (method_called.compare("samp.hub.notifyAll") == 0) {
        ping_service(sock);
    }
    if (method_called.compare("samp.hub.call") == 0) {
        ping_service(sock);
    }
    if (method_called.compare("samp.hub.callAll") == 0) {
        ping_service(sock);
    }
    if (method_called.compare("samp.hub.callAndWait") == 0) {
        ping_service(sock);
    }
    if (method_called.compare("samp.hub.reply") == 0) {
        ping_service(sock);
    }
    if (method_called.compare("samp.hub.getMetadata") == 0) {
        ping_service(sock);
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
    int max_sd, portno,pid,max_clients = 10, i,activity,addrlen;
    socklen_t sockfd, client_socket[10];
    int opt = 1;
    // Prepare TCP/IP structure
    bzero((char *) &serv_addr, sizeof(serv_addr));
    portno                    = 8001;
    serv_addr.sin_family      = AF_INET;
    serv_addr.sin_addr.s_addr = INADDR_ANY;
    serv_addr.sin_port        = htons(portno);
    
    //set of socket descriptors
    fd_set readfds;
    
    // Initialise all client_socket[] to 0 so not checked
    for (i = 0; i < max_clients; ++i) {
        client_socket[i] = 0;
    }

    // ...
    socklen_t hub_socket = socket(AF_INET, SOCK_STREAM, 0);
    
    // Creation of hub main socket
    if (hub_socket < 0) {
        std::string msg = "Hub socket could not be created.";
        throw GException::invalid_value(G_START_HUB, msg);
    }
    
    // Set hub main socket to allow multiple connections
    if (setsockopt(hub_socket, SOL_SOCKET, SO_REUSEADDR, (char *)&opt, sizeof(opt)) < 0) {
        std::string msg = "Hub socket could not be set to multiple connections.";
        throw GException::invalid_value(G_START_HUB, msg);
    }
    
    // Server socket is opened. Now, bind it to the port, with family etc.
    if (bind(hub_socket, (struct sockaddr *) &serv_addr,sizeof(serv_addr)) < 0) {
        std::string msg = "Hub socket could not bind to server socket.";
        throw GException::invalid_value(G_START_HUB, msg);
    }

    // Now start listening for the clients: 5 requests simultaneously pending maximum
    if (listen(hub_socket, 5) < 0) {
        std::string msg = "Hub socket could not start listening.";
        throw GException::invalid_value(G_START_HUB, msg);
    }

    // ...
    socklen_t clilen = sizeof(cli_addr);

    // Main event handling loop
    while (1) {
        
    	// Accept connection from the client 
    	socklen_t newsocket = accept(hub_socket, (struct sockaddr *)&cli_addr, &clilen);
    	if (newsocket < 0) {
            std::string msg = "Client connection to socket not accepted.";
            throw GException::invalid_value(G_START_HUB, msg);
    	}
	//handle_request(newsocket);
	
    	// Create child process to handle the request
        pid = fork();
        if (pid < 0) {
            std::string msg = "No child process created.";
            throw GException::invalid_value(G_START_HUB, msg);
        }
        if (pid == 0) {
            // Child process: client process
            handle_request(newsocket);
            exit(0);
        }
        else {
            close(newsocket);

        }
        

    } // endwhile

    // Return
    return;
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

    // Search for value of specified member
    const GXmlNode* node = xml.element("methodCall", 0);
    if (node != NULL) {
        node = node->element("params", 0);
        if (node != NULL) {
            node = node->element("param", 1);
            if (node != NULL) {
	      node = node->element("value", 0);
	      if (node != NULL) {
		node = node->element("struct", 0);
		if (node != NULL) {
		  int num = node->elements("member");            
                    for (int i = 0; i < num; ++i) {
			   /*std::cout << i << std::endl;
			   std::cout << *(node->element("member", i)) << std::cout;*/
                            const GXmlNode* member = node->element("member", i);
                            //const GXmlNode* member_name = member->element("name", 0);
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
std::list<std::string> GVOHub::get_registrations(const GXml&        xml,
                                          const std::string& name) const                                         
{
    // Declare value
    std::list<std::string> value;
    
    // Search for value of specified member
    const GXmlNode* node = xml.element("methodCall", 0);
    const GXmlNode* subnode;
    const GXmlText* subtext;
    if (node != NULL) {
        node = node->element("params", 0);
        if (node != NULL) {
            node = node->element("param", 1);
            if (node != NULL) {
	      node = node->element("value", 0);
	      if (node != NULL) {
		node = node->element("struct", 0);
		if (node != NULL) {
		  int num = node->elements("member");
                    for (int i = 0; i < num; ++i) {
			   /*std::cout << i << std::endl;
			   std::cout << "Affichage du membre\n";
			   std::cout << *(node->element("member", i)) << std::cout;*/
			   subnode = node->element("member", i);
			   subnode = subnode->element("name", 0);
			   std::cout << node->element("member",i)->element("name", 0)->print(NORMAL,0) << std::cout;
			   //value.push_front(node->element("member",i)->element("name", 0)->print(NORMAL,0));
			   const GXmlText* text = static_cast<const GXmlText*>((*subnode)[0]);
			   if (text != NULL) {
			      value.push_front(text->text());
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
std::string GVOHub::get_callback_port(const GXml&        xml,
                                          const std::string& name) const                                         
{
    // Declare value
    std::string value;
    
    // Search for value of specified member
    const GXmlNode* node = xml.element("methodCall", 0);
    const GXmlNode* subnode;
    const GXmlText* subtext;
    if (node != NULL) {
        node = node->element("params", 0);
        if (node != NULL) {
            subnode = node->element("param", 0);
            if (subnode != NULL) {
	      subnode = subnode->element("value", 0);
	      if (subnode != NULL) {
		const GXmlText* text = static_cast<const GXmlText*>((*subnode)[0]);
		  if (text != NULL) {
		    value = text->text();
		    value.append(" port : ");
		  }
	      }
	    }
	    subnode = node->element("param", 1);
            if (subnode != NULL) {
	      subnode = subnode->element("value", 0);
	      if (subnode != NULL) {
		const GXmlText* text = static_cast<const GXmlText*>((*subnode)[0]);
		  if (text != NULL) {
		    value.append(text->text());
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

