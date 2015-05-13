/***************************************************************************
 *                      GVOHub.hpp - VO SAMP Hub class                     *
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
 * @file GVOHub.hpp
 * @brief SAMP hub class interface definition
 * @author Thierry Louge
 */

#ifndef GVOHUB_HPP
#define GVOHUB_HPP
#define GVOHUB_NB_CLIENTS 5
#define GVOHUB_NB_METHODS 128
#define GVO_HUB_testing 1

/* __ Includes ___________________________________________________________ */
#include <string>
#include <list>
#include <sys/socket.h>
#include <semaphore.h>
#include "GBase.hpp"
#include "GXml.hpp"
#include "GXmlNode.hpp"
#include "GVOApp.hpp"


/***********************************************************************//**
 * @class GVOHub
 *
 * @brief VO SAMP Hub class
 *
 * This class implements a SAMP hub for exchanges 
 * through VO-compatible applications.
 ***************************************************************************/
class GVOHub : public GBase {

public:
    // Constructors and destructors
    GVOHub(void);
    GVOHub(const GVOHub& hub);
    virtual ~GVOHub(void);

    // Operators
    GVOHub& operator=(const GVOHub& hub);

    // Methods
    void        clear(void);
    GVOHub*     clone(void) const;
    std::string print(const GChatter& chatter = NORMAL) const;
    

protected:
    // Protected methods
    void        init_members(void);
    void        copy_members(const GVOHub& client);
    void        free_members(void);
    void        create_samp_file(void);
    void        start_hub(void);
    void        register_service(const GXml& xml,const socklen_t& sock);
    void        ping_service(const socklen_t& sock);
    void        unregister(const GXml& xml,const socklen_t& sock);
    void        handle_request(const socklen_t& sock);
    void        register_metadata(const GXml& xml,const socklen_t& sock);
    void        set_xml_rpc_callback(const GXml& xml,const socklen_t& sock);
    void        get_registered_clients(const GXml& xml,const socklen_t& sock);
    void        get_subscribed_clients(const GXml& xml,const socklen_t& sock);
    void        get_metadata(const GXml& xml,const socklen_t& sock);
    void        get_subscriptions(const GXml& xml,const socklen_t& sock);
    void        declare_subscriptions(const GXml& xml,const socklen_t& sock);
    void	post_string(const std::string& content,const socklen_t& sock) const;
    std::string	receive_string(const socklen_t& sock) const;

    // Low-level methods
    void                   post_string(const std::string& string) const;
    std::string            receive_string(void) const;
    std::string            get_response_value(const GXml& xml,
                                              const std::string& name) const;
    std::list<std::string> get_registrations(const GXml& xml,
                                             const std::string& name) const;
    std::list<std::string> get_clientid(const GXml& xml,
                                             const std::string& name) const;
    std::string get_callback_port(const GXml& xml,
                                             const std::string& name) const;
    void                   get_name_value_pair(const GXmlNode* node,
                                               std::string& name,
                                               std::string& value) const;
    std::string            get_hub_lockfile(void) const;
    std::string		   random_string( size_t length );
    void 		   activate_callbacks(std::string method,char cl_id[31]);
    void 		   post_string_toclient(const std::string& content) const;

    // Protected members
    std::string m_name;        //!< Client name
    std::string m_secret;      //!< Secret Hub key
    std::string m_hub_url;     //!< The XML-RPC endpoint for communication with the hub
    std::string m_hub_host;    //!< Hub host (extracted from XML-RPC endpoint)
    std::string m_hub_port;    //!< Hub port (extracted from XML-RPC endpoint)
    std::string m_version;     //!< The version of the SAMP Standard Profile implemented by the hub
    std::string m_client_key;  //!< Private client key
    std::string m_hub_id;      //!< Hub identifier used by the hub when it sends message itself rather than forwarding from others
    int         m_socket;      //!< Hub socket
    int		cback_socket;  //!< Hub socket to callback clients
    int         m_nb_clients;  //!< Number of already registered clients
    typedef struct {
	char private_key[16];
	char name[32];
	char reference[32];
	char description[256];
	char icon[32];
	char documentation[32];
	char affiliation[32];
	char author_name[32];
	char email[32];
	char homepage[32];
	char port[32];
	char registered_methods[GVOHUB_NB_METHODS][GVOHUB_NB_METHODS]; //128 methods of 128 cars max.
    } connected_shm;
    typedef struct {
    	int registered;
	sem_t lock;
        connected_shm metadatas[GVOHUB_NB_CLIENTS];
    } clients_descriptor;
    clients_descriptor *reg_clients;
    connected_shm *clients;
    void* ptr_mem_partagee; //pointer to shm address
    int shm_handler; //handler to shm 
    
};

#endif /* GVOHUB_HPP */
