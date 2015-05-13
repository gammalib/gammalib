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
#include <signal.h>
#include <sys/types.h> 
#include <sys/wait.h>
#include <sys/ipc.h>
#include <sys/shm.h>


/* __ Method name definitions ____________________________________________ */
#define G_HANDLE_REQUEST                  "GVOHub::handle_request(socklen_t)"
#define G_START_HUB                                     "GVOHub::start_hub()"
#define G_INIT_MEMBERS			  "GVOHub::init_members()"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */
#define CLEF 12345

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
    msg.append("\t<params>\n");
    msg.append("\t\t<param>\n\t\t\t<value/>\n\t\t</param>\n");
    msg.append("\t</params>\n");
    msg.append("</methodResponse>\n");
    
    
    
    post_string(msg,sock);
    //printf("Response message sent\n");
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
    //sem_wait(&reg_clients->lock);
    // Declare message
    std::string msg = "";
    char selfid[124];
    char privatekey[124];
    char cl_id[31];
    std::string secret = random_string(15);
    std::string translator = "http://localhost:8001/xmlrpc";

    for (int i = 0; i < GVOHUB_NB_CLIENTS; i++) {
	printf("Reference of client %d on registration call: %s\n",i,reg_clients->metadatas[i].reference);
	if (strcmp(reg_clients->metadatas[i].reference,"NoRef") == 0) {
		printf("Client %d ready to be filled with new application\n",i);
		sprintf(cl_id,"c%d",i);
		strcpy(reg_clients->metadatas[i].reference,cl_id);
		sprintf(reg_clients->metadatas[i].private_key,secret.c_str());
		break;
	} else {
		printf("Reference of client %d : %s\n",i,reg_clients->metadatas[i].reference);	
	}
    }

    //printf("Message composition 2 \n");
    sprintf(selfid,"\t\t\t\t\t<member><name>samp.self-id</name><value>%s</value></member>\n",cl_id);
    sprintf(privatekey,"\t\t\t\t\t<member><name>samp.private-key</name><value>%s</value></member>\n",secret.c_str());
    reg_clients->registered ++;
    printf(" clients registered %d \n", reg_clients->registered);

    // Set methodResponse elts
    msg.append("<?xml version=\"1.0\" ?>\n");
    msg.append("<methodResponse>\n");
    msg.append("\t<params>\n");
    msg.append("\t\t<param>\n\t\t\t<value>\n\t\t\t\t<struct>\n");
    msg.append("\t\t\t\t\t<member><name>samp.hub-id</name><value>" + m_hub_id + "</value></member>\n");
    msg.append(selfid);
    msg.append(privatekey);
    msg.append("\t\t\t\t\t<member><name>samp.status</name><value>samp.ok</value></member>\n");
    msg.append("\t\t\t\t\t<member><name>samp.url-translator</name><value>" + translator + "</value></member>\n");
    msg.append("\t\t\t\t</struct>\n\t\t\t</value>\n\t\t</param>\n");
    msg.append("\t</params>\n");
    msg.append("</methodResponse>\n");
    post_string(msg,sock);
    
    activate_callbacks("samp.hub.event.register",cl_id);	
    //sem_post(&reg_clients->lock);
    //Return
    return;
}


/***********************************************************************//**
 * @brief 
 *
 *
 ***************************************************************************/
void GVOHub::unregister(const GXml& xml,const socklen_t& sock)
{
    //sem_wait(&reg_clients->lock);
    std::string client_calling = "";
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
    for (int i = 0; i < GVOHUB_NB_CLIENTS; i++) {
	printf("Client %s matching candidate:\n",reg_clients->metadatas[i].private_key);
	if (strcmp(reg_clients->metadatas[i].private_key,client_calling.c_str()) == 0) {
		printf("Client %s unregistration:\n",reg_clients->metadatas[i].reference);
		
		
	
	    // Declare message
	    std::string msg = "";

	    // Set methodResponse elts
	    msg.append("<?xml version=\"1.0\"?>\n");
	    msg.append("<methodResponse>\n");
	    msg.append("\t<params>\n");
	    msg.append("\t\t<param>\n\t\t\t<value/>\n\t\t</param>\n");
	    msg.append("\t</params>\n");
	    msg.append("</methodResponse>\n");
	    post_string(msg,sock);
	    activate_callbacks("samp.hub.event.unregister",reg_clients->metadatas[i].reference);
	    strcpy(reg_clients->metadatas[i].reference,"NoRef");
	    sprintf(reg_clients->metadatas[i].name,"NoName");
	    sprintf(reg_clients->metadatas[i].description,"NoDescription");
	    sprintf(reg_clients->metadatas[i].icon,"NoIcon");
	    sprintf(reg_clients->metadatas[i].documentation,"NoDoc");
	    sprintf(reg_clients->metadatas[i].affiliation,"NoAffiliation");
	    sprintf(reg_clients->metadatas[i].email,"NoMail");
	    sprintf(reg_clients->metadatas[i].author_name,"NoAuthor");
	    sprintf(reg_clients->metadatas[i].homepage,"NoPage");	
	}
    }
    //sem_post(&reg_clients->lock);	
    //Return
    return;
}

/***********************************************************************//**
 * @brief 
 *
 *
 ***************************************************************************/
void GVOHub::register_metadata(const GXml& xml,const socklen_t& sock)
{
    char buffer[4096];
    #if GVO_HUB_testing == 1
    printf("Client Metadata received\n");
    #endif
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

    for (int i = 0; i < GVOHUB_NB_CLIENTS; i++) {
	if (strcmp(reg_clients->metadatas[i].private_key,client_calling.c_str()) == 0) {
		printf("Client %s ready to be filled with metadata\n",reg_clients->metadatas[i].reference);
		std::string client_name = get_response_value(xml, "samp.name");
		sprintf(reg_clients->metadatas[i].name,client_name.c_str());
   		std::string client_description = get_response_value(xml, "samp.description.text");
		sprintf(reg_clients->metadatas[i].description,client_description.c_str());
    		std::string client_icon = get_response_value(xml, "samp.icon.url");
		sprintf(reg_clients->metadatas[i].icon,client_icon.c_str());
    		std::string client_doc = get_response_value(xml, "samp.documentation.url");
		sprintf(reg_clients->metadatas[i].documentation,client_doc.c_str());
    		std::string client_affi = get_response_value(xml, "author.affiliation");
		sprintf(reg_clients->metadatas[i].affiliation,client_affi.c_str());
    		std::string client_mail = get_response_value(xml, "author.email");
		sprintf(reg_clients->metadatas[i].email,client_mail.c_str());
    		std::string client_authname = get_response_value(xml, "author.name");
		sprintf(reg_clients->metadatas[i].author_name,client_authname.c_str());
   		std::string client_page = get_response_value(xml, "home.page");
		sprintf(reg_clients->metadatas[i].homepage,client_page.c_str());
		// Declare message
		std::string msg = "";

		//The hub response to Metadata declaration has no usefull content and can be discarded by the client.
		// Set methodResponse elts
		msg.append("<?xml version=\"1.0\"?>\n");
		msg.append("<methodResponse><params><param><value/></param></params>\n");
		msg.append("</methodResponse>\n");
		post_string(msg,sock);
		activate_callbacks("samp.hub.event.metadata",reg_clients->metadatas[i].reference);
		break;
	} else {
		if (i == GVOHUB_NB_CLIENTS - 1) {
			printf("OOPS... unable to find client: %s\n",client_calling.c_str());
		}
	}
    }
    

    
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
    int i;
    printf("Client callback registration\n");

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
    // Search for metadata values to store
    std::string client_port = get_callback_port(xml, "samp.hub.setXmlrpcCallback");
    std::cout << client_port << std::endl;
    
    for (i = 0; i < GVOHUB_NB_CLIENTS; i++) {
	//printf("comparing Reference of client %s with seeked %s\n",reg_clients->metadatas[i].private_key,client_calling.c_str());
	if (strcmp(reg_clients->metadatas[i].private_key,client_calling.c_str()) == 0) {
		printf("Client %s ready to be filled with port information\n",reg_clients->metadatas[i].reference);
		sprintf(reg_clients->metadatas[i].port,client_port.c_str());
		break;
	}
    }
    printf("%s\n",reg_clients->metadatas[i].port);
    // Declare message
    std::string msg = "";

    //The hub response to rpc callback declaration has no usefull content and can be discarded by the client.
    // Set methodResponse elts
    msg.append("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
    msg.append("<methodResponse>\n");
    msg.append("\t<params>\n");
    msg.append("\t\t<param>\n\t\t\t<value>\n\t\t\t\t<struct>\n");
    msg.append("\t\t\t\t\t<member>\n\t\t\t\t\t\t<name>samp.status</name>\n\t\t\t\t\t\t<value>samp.ok</value>\n\t\t\t\t\t</member>\n");
    msg.append("\t\t\t\t</struct>\n\t\t\t</value>\n\t\t</param>\n");
    msg.append("\t</params>\n");
    msg.append("</methodResponse>\n");
    post_string(msg,sock);
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
    char reference[64];
    std::string tofind;
    std::size_t found;
    std::string str2;
    // Declare message
    std::string msg = "";

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

    // Set methodResponse elts
    msg.append("<?xml version=\'1.0\' encoding=\"UTF-8\"?>\n");
    msg.append("<methodResponse>\n");
    msg.append("\t<params>\n");
    msg.append("\t\t<param>\n\t\t\t<value>\n\t\t\t\t<array>\n\t\t\t\t\t<data>\n");
    msg.append("\t\t\t\t\t\t<value>gammalib_hub</value>\n");


    for (int i = 0; i < GVOHUB_NB_CLIENTS; i++) {
	if (strcmp(reg_clients->metadatas[i].reference,"NoRef") == 0) {
		printf("Client %d not registered\n",i);
	} else {
		//Do not send back current client's registration
		if (strcmp(reg_clients->metadatas[i].private_key,client_calling.c_str()) != 0) {
			msg.append("\t\t\t\t\t\t<value>");
			msg.append(reg_clients->metadatas[i].reference);
			msg.append("</value>\n");
			printf("Reference of client %d : %s\n",i,reg_clients->metadatas[i].reference);
		}
	}
    }
    
    msg.append("\t\t\t\t\t</data>\n\t\t\t\t</array>\n\t\t\t</value>\n");
    msg.append("\t\t</param>\n");
    msg.append("\t</params>\n");
    msg.append("</methodResponse>\n");
    printf("Message ready to be sent \n");
    post_string(msg,sock);
    
    //Return
    return;
  
}
/***********************************************************************//**
 * @brief 
 *
 *
 ***************************************************************************/
void GVOHub::get_subscribed_clients(const GXml& xml,const socklen_t& sock)
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

    msg.append("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
    msg.append("<methodResponse>\n");
    msg.append("\t<params>\n");
    msg.append("\t\t<param>\n\t\t\t<value>\n\t\t\t\t<struct>\n");
    msg.append("\t\t\t\t\t<member>\n\t\t\t\t\t\t<name>c0</name>\n\t\t\t\t\t\t<value><struct></struct></value>\n\t\t\t\t\t</member>\n");
    msg.append("\t\t\t\t\t<member>\n\t\t\t\t\t\t<name>gammalib_hub</name>\n\t\t\t\t\t\t<value><struct></struct></value>\n\t\t\t\t\t</member>\n");
    msg.append("\t\t\t\t</struct>\n\t\t\t</value>\n\t\t</param>\n");
    msg.append("\t</params>\n");
    msg.append("</methodResponse>\n");


    //int length = msg.length();
    post_string(msg,sock);
    
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
    int i,j;
    std::ofstream myfile;
    //char content[GVOHUB_NB_METHODS];
    printf("Declare Subscriptions\n");
    
    std::string client_calling = "";
    std::string tempo;
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
    
    std::list<std::string> client_subscriptions = get_registrations(xml, "samp.hub.declareSubscriptions");
    std::cout << "Client registered to following methods calls:" << std::endl;
     
    for (int i = 0; i < GVOHUB_NB_CLIENTS; i++) {
	//printf("comparing Reference of client %s with seeked %s\n",reg_clients->metadatas[i].private_key,client_calling.c_str());
	if (strcmp(reg_clients->metadatas[i].private_key,client_calling.c_str()) == 0) {
      		for( j=0,p=client_subscriptions.begin(); j<client_subscriptions.size()-1; j++,p++ ) {
        		std::cout << *p << std::endl;
			tempo = *p;
			//content = *p;
			sprintf(reg_clients->metadatas[i].registered_methods[j],tempo.c_str());
			/*printf("Registered: %s\n",reg_clients->metadatas[i].registered_methods[j]);
			printf("Registration done\n");*/
		}
	}
    }
	
      
    
    // Declare message
    std::string msg = "";
    // Set methodResponse elts
    msg.append("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
    msg.append("<methodResponse>\n");
    msg.append("\t<params>\n");
    msg.append("\t\t<param>\n\t\t\t<value>\n\t\t\t\t<struct>\n");
    msg.append("\t\t\t\t\t<member>\n\t\t\t\t\t\t<name>samp.status</name>\n\t\t\t\t\t\t<value>samp.ok</value>\n\t\t\t\t\t</member>\n");
    msg.append("\t\t\t\t</struct>\n\t\t\t</value>\n\t\t</param>\n");
    msg.append("\t</params>\n");
    msg.append("</methodResponse>\n");

    post_string(msg,sock);
    //Return
    return;
}

/***********************************************************************//**
 * @brief 
 *
 *
 ***************************************************************************/
void GVOHub::get_metadata(const GXml& xml,const socklen_t& sock)
{
    std::list<std::string>::iterator p;
    int i;
    std::ofstream myfile;
    std::string content;
    printf("Retrieving metadata\n");
    
    std::string client_calling = "";
    // Search for value of methodName
    const GXmlNode* node = xml.element("methodCall", 0);
    node = node->element("params", 0);
    node = node->element("param", 1);
    node = node->element("value", 0);
    if (node != NULL) {
        const GXmlText* text = static_cast<const GXmlText*>((*node)[0]);
        if (text != NULL) {
            client_calling = text->text();
        }
    }
    // Declare message
    std::string msg = "";
    // Set methodResponse elts
    msg.append("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
    msg.append("<methodResponse>\n");
    msg.append("\t<params>\n");
    msg.append("\t\t<param>\n\t\t\t<value>\n\t\t\t\t<struct>\n");
    msg.append("\t\t\t\t\t<member>\n\t\t\t\t\t\t<name>samp.status</name>\n\t\t\t\t\t\t<value>samp.ok</value>\n\t\t\t\t\t</member>\n");
	
    for (int i = 0; i < GVOHUB_NB_CLIENTS; i++) {
	if (strcmp("gammalib_hub",client_calling.c_str()) == 0) {
		msg.append("\t\t\t\t\t<member>\n\t\t\t\t\t\t<name>samp.name</name>\n\t\t\t\t\t\t<value>gammalib_hub</value>\n\t\t\t\t\t</member>\n");
   		msg.append("\t\t\t\t\t<member>\n\t\t\t\t\t\t<name>samp.description.text</name>\n\t\t\t\t\t\t<value>Gamma lib hub.</value>\n\t\t\t\t\t</member>\n");
		msg.append("\t\t\t\t\t<member>\n\t\t\t\t\t\t<name>samp.icon.url</name>\n\t\t\t\t\t\t<value></value>\n\t\t\t\t\t</member>\n");
		msg.append("\t\t\t\t\t<member>\n\t\t\t\t\t\t<name>samp.documentation.url</name>\n\t\t\t\t\t\t<value>gammalib github documentation</value>\n\t\t\t\t\t</member>\n");
		msg.append("\t\t\t\t\t<member>\n\t\t\t\t\t\t<name>author.affiliation</name>\n\t\t\t\t\t\t<value>IRAP, CNRS Toulouse</value>\n\t\t\t\t\t</member>\n");
		msg.append("\t\t\t\t\t<member>\n\t\t\t\t\t\t<name>author.email</name>\n\t\t\t\t\t\t<value>dontcontactme@please</value>\n\t\t\t\t\t</member>\n");
		msg.append("\t\t\t\t\t<member>\n\t\t\t\t\t\t<name>author.name</name>\n\t\t\t\t\t\t<value>J. Knödlseder, T. Louge</value>\n\t\t\t\t\t</member>\n");
		msg.append("\t\t\t\t\t<member>\n\t\t\t\t\t\t<name>home.page</name>\n\t\t\t\t\t\t<value>no dedicated page</value>\n\t\t\t\t\t</member>\n");
		break;	
	}
	if (strcmp(reg_clients->metadatas[i].reference,client_calling.c_str()) == 0) {
		printf("Client %s metadata retrieval:\n",reg_clients->metadatas[i].reference);
		msg.append("\t\t\t\t\t<member>\n\t\t\t\t\t\t<name>samp.name</name>\n\t\t\t\t\t\t<value>");
		msg.append(reg_clients->metadatas[i].name);
		msg.append("</value>\n\t\t\t\t\t</member>\n");
   		msg.append("\t\t\t\t\t<member>\n\t\t\t\t\t\t<name>samp.description.text</name>\n\t\t\t\t\t\t<value>");
		msg.append(reg_clients->metadatas[i].description);
		msg.append("</value>\n\t\t\t\t\t</member>\n");
		msg.append("\t\t\t\t\t<member>\n\t\t\t\t\t\t<name>samp.icon.url</name>\n\t\t\t\t\t\t<value>");
		msg.append(reg_clients->metadatas[i].icon);
		msg.append("</value>\n\t\t\t\t\t</member>\n");
		msg.append("\t\t\t\t\t<member>\n\t\t\t\t\t\t<name>samp.documentation.url</name>\n\t\t\t\t\t\t<value>");
		msg.append(reg_clients->metadatas[i].documentation);
		msg.append("</value>\n\t\t\t\t\t</member>\n");
		msg.append("\t\t\t\t\t<member>\n\t\t\t\t\t\t<name>author.affiliation</name>\n\t\t\t\t\t\t<value>");
		msg.append(reg_clients->metadatas[i].affiliation);
		msg.append("</value>\n\t\t\t\t\t</member>\n");
		msg.append("\t\t\t\t\t<member>\n\t\t\t\t\t\t<name>author.email</name>\n\t\t\t\t\t\t<value>");
		msg.append(reg_clients->metadatas[i].email);
		msg.append("</value>\n\t\t\t\t\t</member>\n");
		msg.append("\t\t\t\t\t<member>\n\t\t\t\t\t\t<name>author.name</name>\n\t\t\t\t\t\t<value>");
		msg.append(reg_clients->metadatas[i].author_name);
		msg.append("</value>\n\t\t\t\t\t</member>\n");
		msg.append("\t\t\t\t\t<member>\n\t\t\t\t\t\t<name>home.page</name>\n\t\t\t\t\t\t<value>");
		msg.append(reg_clients->metadatas[i].homepage);
		msg.append("</value>\n\t\t\t\t\t</member>\n");
		break;	
	} else {
		if (i == GVOHUB_NB_CLIENTS - 1) {
			printf("OOPS... unable to find metadata for client: %s\n",client_calling.c_str());
		}
	}
    }
    
    msg.append("\t\t\t\t</struct>\n\t\t\t</value>\n\t\t</param>\n");
    msg.append("\t</params>\n");
    msg.append("</methodResponse>\n");
    printf("Metadata for %s : \n %s",client_calling.c_str(),msg.c_str());
    post_string(msg,sock);
   
    //Return
    return;
}
/***********************************************************************//**
 * @brief 
 *
 *
 ***************************************************************************/
void GVOHub::get_subscriptions(const GXml& xml,const socklen_t& sock)
{
    std::list<std::string>::iterator p;
    int i,j;
    std::ofstream myfile;
    std::string content;
    printf("Retrieving subscriptions\n");
    std::string client_calling = "";
    // Search for value of methodName
    const GXmlNode* node = xml.element("methodCall", 0);
    node = node->element("params", 0);
    node = node->element("param", 1);
    node = node->element("value", 0);
    if (node != NULL) {
        const GXmlText* text = static_cast<const GXmlText*>((*node)[0]);
        if (text != NULL) {
            client_calling = text->text();
        }
    }
    // Declare message
    std::string msg = "";
    // Set methodResponse elts
    msg.append("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
    msg.append("<methodResponse>\n");
    msg.append("\t<params>\n");
    msg.append("\t\t<param>\n\t\t\t<value>\n\t\t\t\t<struct>\n");
    //msg.append("\t\t\t\t\t<member>\n\t\t\t\t\t\t<name>samp.status</name>\n\t\t\t\t\t\t<value>samp.ok</value>\n\t\t\t\t\t</member>\n");
    char tempo[128];
    for (i = 0; i < GVOHUB_NB_CLIENTS; i++) {
	if (strcmp(reg_clients->metadatas[i].reference,client_calling.c_str()) == 0) {
		//printf("retrieving client %s\n",reg_clients->metadatas[i].reference);
		for (j = 0; j < GVOHUB_NB_METHODS; j++) {
			strcpy(tempo,reg_clients->metadatas[i].registered_methods[j]);
			if (strcmp("",tempo) == 0) {
				continue;
			} else {
				msg.append("\t\t\t\t\t<member>\n\t\t\t\t\t\t<name>");
				msg.append(tempo);
				msg.append("</name>\n\t\t\t\t\t\t<value></value>\n\t\t\t\t\t</member>\n");
			}
		}
	}
     }
    
    msg.append("\t\t\t\t</struct>\n\t\t\t</value>\n\t\t</param>\n");
    msg.append("\t</params>\n");
    msg.append("</methodResponse>\n");
    post_string(msg,sock);
   
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
    //m_secret is a private string to identify clients. Every client has his own m_secret
    m_secret     = "";
    m_hub_url    = "http://localhost:8001/xmlrpc";
    m_hub_host   = "127.0.0.1";
    m_hub_port   = "8001";
    m_version    = "1.3";
    m_client_key.clear();
    m_hub_id     = "gammalib_hub";
    //m_client_id.clear();
    m_socket     = -1;        // Signals no socket
    m_nb_clients = 0;	      //  No registered clients at initialization

    //puts the structure array of clients in shared memory
    if ((shm_handler = shmget(CLEF, sizeof(clients_descriptor), 0666 | IPC_CREAT)) < 0)	
	{
		std::string msg = "Could not open shared memory.";
        	throw GException::invalid_value(G_INIT_MEMBERS, msg);
	}
    if ((ptr_mem_partagee = shmat(shm_handler, NULL, 0)) == (void*) -1)
	{
		std::string msg = "Could not attach shared memory.";
        	throw GException::invalid_value(G_INIT_MEMBERS, msg);
	}
    reg_clients = (clients_descriptor *)ptr_mem_partagee;
    reg_clients->registered = 0;
    for (int i = 0; i < GVOHUB_NB_CLIENTS; i++) {
	strcpy(reg_clients->metadatas[i].reference,"NoRef");
	strcpy(reg_clients->metadatas[i].private_key,"Unknown");
	strcpy(reg_clients->metadatas[i].name,"Unknown");
	strcpy(reg_clients->metadatas[i].icon,"Unknown");
	strcpy(reg_clients->metadatas[i].documentation,"Unknown");
	strcpy(reg_clients->metadatas[i].affiliation,"Unknown");
	strcpy(reg_clients->metadatas[i].author_name,"Unknown");
	strcpy(reg_clients->metadatas[i].email,"Unknown");
	strcpy(reg_clients->metadatas[i].homepage,"Unknown");
	strcpy(reg_clients->metadatas[i].port,"Unknown");
	strcpy(reg_clients->metadatas[i].description,"Unknown");
	
    }
    for (int j = 0; j < GVOHUB_NB_CLIENTS; j++) {
	printf("%s\n",reg_clients->metadatas[j].reference);
    }
    //sem_init(&reg_clients->lock, 0, 1);
    // Creates Hub
    printf("Creating the hub\n");
    create_samp_file();
    
    // Starts the hub
    printf("Starting the hub\n");
    start_hub();
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
    m_hub_url        = hub.m_hub_url;
    m_hub_host       = hub.m_hub_host;
    m_hub_port       = hub.m_hub_port;
    m_version        = hub.m_version;
    m_client_key     = hub.m_client_key;
    m_hub_id         = hub.m_hub_id;
    m_socket         = hub.m_socket;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GVOHub::free_members(void)
{ 
    struct shmid_ds *buf;

    std::cout << "I am here" << std::endl;
    std::string lockurl = get_hub_lockfile();
    printf("Deleting file: %s\n",lockurl.c_str());
    fflush(stdout);
    std::remove(lockurl.c_str());
    shmdt(ptr_mem_partagee);
    shmctl(shm_handler, IPC_RMID, buf);
    // Return
    return;
}


/***********************************************************************//**
 * @brief Creates the lockfile, fill it
 *
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

        // If successfull writes basic hub configuration to the lockfile
        fputs ("# SAMP lockfile \n# Required keys:\n",fptr);
        char sampsecret[256] = "samp.secret=";
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
    char buffer[8192];
    bzero(buffer,8192);
    int i,j;
    // Declare empty XML document
    GXml xml;
    
    // Read buffer
    int n = read(sock,buffer,8192);
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
    if (method_called.compare("samp.hub.ping") == 0) {
        ping_service(sock);
    }
    if (method_called.compare("samp.hub.register") == 0) {
        register_service(xml, sock);
    }
    if (method_called.compare("samp.hub.unregister") == 0) {
        unregister(xml,sock);
    }
    if (method_called.compare("samp.hub.declareMetadata") == 0) {
        register_metadata(xml,sock);
    }
    if (method_called.compare("samp.hub.declareSubscriptions") == 0) {
        declare_subscriptions(xml,sock);
    }
    if (method_called.compare("samp.hub.setXmlrpcCallback") == 0) {
        set_xml_rpc_callback(xml,sock);
    }
    if (method_called.compare("samp.hub.getSubscriptions") == 0) {
        get_subscriptions(xml, sock);
    }
    if (method_called.compare("samp.hub.getRegisteredClients") == 0) {
        get_registered_clients(xml, sock);
    }
    if (method_called.compare("samp.hub.getSubscribedClients") == 0) {
        get_subscribed_clients(xml, sock);
    }
    if (method_called.compare("samp.hub.getMetadata") == 0) {
        get_metadata(xml, sock);
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

    //Catching ctrl+C to close SHM
    //signal(SIGINT, intHandler);

    // Main event handling loop
    
    while (1) {
        
    	// Accept connection from the client 
    	socklen_t newsocket = accept(hub_socket, (struct sockaddr *)&cli_addr, &clilen);
    	if (newsocket < 0) {
            std::string msg = "Client connection to socket not accepted.";
            throw GException::invalid_value(G_START_HUB, msg);
    	}
	/*
	handle_request(newsocket);
	close(newsocket);
	*/
    	// Create child process to handle the request
        
	pid = fork();
        if (pid < 0) {
            std::string msg = "No child process created.";
            throw GException::invalid_value(G_START_HUB, msg);
        }
        if (pid == 0) {
            // Child process: client process
            handle_request(newsocket);
	    close(newsocket);
            exit(0);
        }
        else {
	    //Father process code zone: set SIG_IGN to child
	    // so that kernel can close it without senndit it to zombie state <defunct>
	    int status = 0;
	    waitpid( pid, &status, 0 );
	    signal(pid, SIG_IGN);	
	}

    } 
    
    
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
			   //std::cout << node->element("member",i)->element("name", 0)->print(NORMAL,0) << std::cout;
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
		  /*if (text != NULL) {
		    value = text->text();
		    value.append(" port : ");
		  }*/
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

/***********************************************************************//**
 * @brief Post string content to Hub
 *
 * @param[in] content String content to post
 *
 * Posts the content of a string to the Hub.
 *
 * The method does nothing if no Hub connection has been established.
 ***************************************************************************/
void GVOHub::post_string(const std::string& content,const socklen_t& sock) const
{
    // Continue only if Hub connection has been established
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
	#if GVO_HUB_testing == 1
	printf(post.c_str());
	#endif
        // Send content to socket
        bool done = false;
        do {
            int length      = post.length();
            int sent_length = send(sock, post.c_str(), length, 0);
	    //printf("Sent through socket: \n %s",post.c_str());
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
std::string GVOHub::receive_string(const socklen_t& sock) const
{
    // Initialise empty string
    std::string result = "";

    // Continue only if Hub connection has been established
    printf("sock = %d\n",sock); 
    if (sock != -1) {
	
        // Define buffer
        char buffer[4097];

        // Read buffer until it is empty
        int n = 0;
	int m = 0;
        //do {
	   
	    printf("reading socket\n"); 
            n = recv(sock, buffer, 4096,0);
	    printf("%d octets received: message %s\n",n,buffer);
	    fflush(stdout);
            if (n > 0) {
                //buffer[n+1] = '\0';
                result.append(std::string(buffer));
            }
            
            printf("reading socket\n"); 
            m = recv(sock, buffer, 4096,0);
	    printf("%d octets received: message %s\n",n,buffer);
	    fflush(stdout);
            if (n > 0) {
                buffer[n+m+1] = '\0';
                result.append(std::string(buffer));
            }
        //} while (n > 0);

    } // endif: Hub connection had been established
    printf("receive_string now returning\n"); 
    fflush(stdout);
    // Return result
    return result;
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
std::list<std::string> GVOHub::get_clientid(const GXml&        xml,
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
	      //std::cout << node->element("member",i)->element("name", 0)->print(NORMAL,0) << std::cout;
	      const GXmlText* text = static_cast<const GXmlText*>((*node)[0]);
	      if (text != NULL) {
		value.push_front(text->text());
		#if GVO_HUB_testing == 1
		std::cout << "Client id for which information is seeked:" + text->text() + "\n";
		#endif
	      }
            }
	}
    }
		

    // Return value
    return value;
}
/***********************************************************************//**
 * @brief 
 *
 * @param[in] 
 *
 * 
 ***************************************************************************/
std::string GVOHub::random_string( size_t length )
{
   srand(time(0));
   std::string str = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz";
   int pos;
   while(str.size() != length) {
    pos = ((rand() % (str.size() - 1)));
    str.erase (pos, 1);
   }
   return str;
}

/***********************************************************************//**
 * @brief Test if any registered client needs to know that the method 
		has been called
 *
 * @param[in] method: name of the method, cl_id id of calling client
 *
 * 
 *
 * 
 ***************************************************************************/
void GVOHub::activate_callbacks(std::string method,char cl_id[31])
{
    int i,j,k;
    char port[16];
    char pkey[16];
    std::string spkey;
    
    for (i = 0; i < GVOHUB_NB_CLIENTS; i++) {
	for (j = 0; j < GVOHUB_NB_METHODS; j++) {
		if (strcmp("",reg_clients->metadatas[i].registered_methods[j]) == 0) {
			continue;
		} else {
			//printf("Comparing %s and %s \n",method.c_str(),reg_clients->metadatas[i].registered_methods[j]);
			if (strcmp(method.c_str(),reg_clients->metadatas[i].registered_methods[j]) == 0) {
				strcpy(port,reg_clients->metadatas[i].port);
				printf("Callback needed on port %s\n",port);
				std::string port_string = std::string(port);
				//Eliminate http://127.0.0.1/ from port string -----> improve this code is necessary
				port_string.replace(0,17,"");
				printf("Callback needed on port %s\n",port_string.c_str());
				std::size_t pos = port_string.find("/");
				std::string port_only = port_string.substr(0,pos);
				printf("%s\n",port_only.c_str());
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
					cback_socket = socket(ptr->ai_family,
						          ptr->ai_socktype,
						          ptr->ai_protocol);

					// Connect to socket if socket is valid
					if (cback_socket != -1) {
					    if (::connect(cback_socket,
						          ptr->ai_addr,
						          ptr->ai_addrlen) == -1) {
						close(cback_socket);
						cback_socket = -1;
					    }
					}

				    } // endfor: looped through all results

				} // endif: server information was valid
				std::string msg_cback = "";

				if (strcmp(method.c_str(),"samp.hub.event.register") == 0) {
					// Set methodResponse elts for samp.hub.event.register
					msg_cback.append("<?xml version=\"1.0\"?>\n");
				
					msg_cback.append("<methodCall><methodName>samp.client.receiveNotification</methodName>\n");
					msg_cback.append("\t<params>\n");
					msg_cback.append("\t\t<param><value>");
					msg_cback.append(reg_clients->metadatas[i].private_key);
					msg_cback.append("</value></param>\n");
					msg_cback.append("\t\t<param><value>"+m_hub_id+"</value></param>\n");
					msg_cback.append("\t\t<param>\n\t\t\t<value>\n\t\t\t\t<struct>\n");
	    				msg_cback.append("\t\t\t\t\t<member><name>samp.mtype</name><value>samp.hub.event.register</value></member>\n");
					msg_cback.append("\t\t\t\t\t<member><name>samp.params</name><value>\n\t\t\t\t\t\t<struct>\n\t\t\t\t\t\t\t<member>");
					msg_cback.append("<name>id</name><value>");
					msg_cback.append(cl_id);
					msg_cback.append("</value></member>\n\t\t\t\t\t\t</struct>\n\t\t\t\t\t</value></member>\n");
	    				msg_cback.append("\t\t\t\t</struct>\n\t\t\t</value>\n\t\t</param>\n");
					msg_cback.append("\t</params>\n");
					 
					msg_cback.append("</methodCall>\n");
					post_string_toclient(msg_cback);
				}
				if (strcmp(method.c_str(),"samp.hub.event.unregister") == 0) {
					// Set methodResponse elts for samp.hub.event.unregister
					msg_cback.append("<?xml version=\"1.0\"?>\n");
				
					msg_cback.append("<methodCall><methodName>samp.client.receiveNotification</methodName>\n");
					msg_cback.append("\t<params>\n");
					msg_cback.append("\t\t<param><value>");
					msg_cback.append(reg_clients->metadatas[i].private_key);
					msg_cback.append("</value></param>\n");
					msg_cback.append("\t\t<param><value>"+m_hub_id+"</value></param>\n");
					msg_cback.append("\t\t<param>\n\t\t\t<value>\n\t\t\t\t<struct>\n");
	    				msg_cback.append("\t\t\t\t\t<member><name>samp.mtype</name><value>samp.hub.event.unregister</value></member>\n");
					msg_cback.append("\t\t\t\t\t<member><name>samp.params</name><value>\n\t\t\t\t\t\t<struct>\n\t\t\t\t\t\t\t<member>");
					msg_cback.append("<name>id</name><value>");
					msg_cback.append(cl_id);
					msg_cback.append("</value></member>\n\t\t\t\t\t\t</struct>\n\t\t\t\t\t</value></member>\n");
	    				msg_cback.append("\t\t\t\t</struct>\n\t\t\t</value>\n\t\t</param>\n");
					msg_cback.append("\t</params>\n");
					 
					msg_cback.append("</methodCall>\n");
					post_string_toclient(msg_cback);
					
				}
			        if (strcmp(method.c_str(),"samp.hub.event.metadata") == 0) {
					// Set methodResponse elts for samp.hub.event.unregister
					msg_cback.append("<?xml version=\"1.0\"?>\n");
				
					msg_cback.append("<methodCall><methodName>samp.client.receiveNotification</methodName>\n");
					msg_cback.append("\t<params>\n");
					msg_cback.append("\t\t<param><value>");
					msg_cback.append(reg_clients->metadatas[i].private_key);
					msg_cback.append("</value></param>\n");
					msg_cback.append("\t\t<param><value>"+m_hub_id+"</value></param>\n");
					msg_cback.append("\t\t<param>\n\t\t\t<value>\n\t\t\t\t<struct>\n");
	    				msg_cback.append("\t\t\t\t\t<member><name>samp.mtype</name><value>samp.hub.event.metadata</value></member>\n");
					msg_cback.append("\t\t\t\t\t<member><name>samp.params</name><value>\n\t\t\t\t\t\t<struct>\n\t\t\t\t\t\t\t");
					msg_cback.append("\t\t\t\t\t<member><name>id</name><value>");
					for (k = 0; k < GVOHUB_NB_CLIENTS; k++) {
						if (strcmp(cl_id,reg_clients->metadatas[k].reference) == 0) {
							break;					
						}
					}
					msg_cback.append(reg_clients->metadatas[k].reference);
					msg_cback.append("</value>\n\t\t\t\t\t</member>\n");
					msg_cback.append("\t\t\t\t\t<member><name>metadata</name><value>\n\t\t\t\t\t\t<struct>\n\t\t\t\t\t\t\t");
					msg_cback.append("\t\t\t\t\t<member>\n\t\t\t\t\t\t<name>samp.name</name>\n\t\t\t\t\t\t<value>");
					msg_cback.append(reg_clients->metadatas[k].name);
					msg_cback.append("</value>\n\t\t\t\t\t</member>\n");
   					msg_cback.append("\t\t\t\t\t<member>\n\t\t\t\t\t\t<name>samp.description.text</name>\n\t\t\t\t\t\t<value>");
					msg_cback.append(reg_clients->metadatas[k].description);
					msg_cback.append("</value>\n\t\t\t\t\t</member>\n");
					msg_cback.append("\t\t\t\t\t<member>\n\t\t\t\t\t\t<name>samp.icon.url</name>\n\t\t\t\t\t\t<value>");
					msg_cback.append(reg_clients->metadatas[k].icon);
					msg_cback.append("</value>\n\t\t\t\t\t</member>\n");
					msg_cback.append("\t\t\t\t\t<member>\n\t\t\t\t\t\t<name>samp.documentation.url</name>\n\t\t\t\t\t\t<value>");
					msg_cback.append(reg_clients->metadatas[k].documentation);
					msg_cback.append("</value>\n\t\t\t\t\t</member>\n");
					msg_cback.append("\t\t\t\t\t<member>\n\t\t\t\t\t\t<name>author.affiliation</name>\n\t\t\t\t\t\t<value>");
					msg_cback.append(reg_clients->metadatas[k].affiliation);
					msg_cback.append("</value>\n\t\t\t\t\t</member>\n");
					msg_cback.append("\t\t\t\t\t<member>\n\t\t\t\t\t\t<name>author.email</name>\n\t\t\t\t\t\t<value>");
					msg_cback.append(reg_clients->metadatas[k].email);
					msg_cback.append("</value>\n\t\t\t\t\t</member>\n");
					msg_cback.append("\t\t\t\t\t<member>\n\t\t\t\t\t\t<name>author.name</name>\n\t\t\t\t\t\t<value>");
					msg_cback.append(reg_clients->metadatas[k].author_name);
					msg_cback.append("</value>\n\t\t\t\t\t</member>\n");
					msg_cback.append("\t\t\t\t\t<member>\n\t\t\t\t\t\t<name>home.page</name>\n\t\t\t\t\t\t<value>");
					msg_cback.append(reg_clients->metadatas[k].homepage);
					msg_cback.append("</value>\n\t\t\t\t\t</member>\n");	
					msg_cback.append("</struct></value></member>\n");				
					msg_cback.append("</struct>\n\t\t\t\t\t</value></member>\n");
	    				msg_cback.append("\t\t\t\t</struct>\n\t\t\t</value>\n\t\t</param>\n");
					msg_cback.append("\t</params>\n");
					 
					msg_cback.append("</methodCall>\n");
					post_string_toclient(msg_cback);
					
				}
				

				
			}
		}
	}
	
    }
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
void GVOHub::post_string_toclient(const std::string& content) const
{
    // Continue only if Hub connection has been established
    if (cback_socket != -1) {

        // Determine content length
        int length = content.length();

        // Set prefix
        std::string prefix = "POST /xmlrpc HTTP/1.1\n"
			     "Connection: keep-alive\n"
                             "User-Agent: GammaLib\n"
			     "Host: localhost:"+gammalib::str(cback_socket)+"\n"
                             "Content-Type: text/xml\n"
                             "Content-Length: "+gammalib::str(length)+"\n\n";

        // Build post string
        std::string post = prefix + content;

        // Send content to socket
        bool done = false;
        do {
            int length      = post.length();
            int sent_length = send(cback_socket, post.c_str(), length, 0);
	    #if GVO_HUB_testing == 1
	    printf("Send to client: %s\n",post.c_str());
	    #endif
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
	  printf("receive socket instruction from callback response\n");
	  fflush(stdout);
	  n = recv(cback_socket, buffer, 1000, 0);
	  printf("data received: %s\n",buffer);
	  fflush(stdout);
	  if (n > 0) {
	      buffer[n+1] = '\0';
	      
	  }
	} while (n > 0);
	/*printf("%s",buffer);
	// Declare message
	    std::string msg = "";

	    // Set methodResponse elts
	    msg.append("<?xml version=\"1.0\"?>\n");
	    msg.append("<methodResponse>\n");
	    msg.append("\t<params>\n");
	    msg.append("\t\t<param>\n\t\t\t<value/>\n\t\t</param>\n");
	    msg.append("\t</params>\n");
	    msg.append("</methodResponse>\n");
	    
	    
	    
	    post_string(msg,cback_socket);*/
    } // endif: Hub connection had been established
	
    // Return
    return;
}

