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

/* __ Definitions ________________________________________________________ */
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


/***********************************************************************//**
 * @class GVOHub
 *
 * @brief VO SAMP Hub class
 *
 * This class implements a SAMP hub for exchanges through VO-compatible
 * applications.
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
    std::string classname(void) const;
    void        start(void);
    std::string print(const GChatter& chatter = NORMAL) const;
    

protected:
    // Protected methods
    void                   init_members(void);
    void                   copy_members(const GVOHub& client);
    void                   free_members(void);
    void                   start_hub(void);
    void                   handle_request(const socklen_t& sock);
    void                   request_ping(const socklen_t& sock);
    void                   request_register(const GXml& xml,
                                            const socklen_t& sock);
    void                   request_unregister(const GXml& xml,
                                              const socklen_t& sock);
    void                   request_declare_metadata(const GXml& xml,
                                                    const socklen_t& sock);
    void                   request_declare_subscriptions(const GXml& xml,
                                                         const socklen_t& sock);
    void                   request_set_xml_rpc_callback(const GXml& xml,
                                                        const socklen_t& sock);
    void                   request_get_subscriptions(const GXml& xml,
                                                     const socklen_t& sock);
    void                   request_get_registered_clients(const GXml& xml,
                                                          const socklen_t& sock);
    void                   request_get_subscribed_clients(const GXml& xml,
                                                          const socklen_t& sock);
    void                   request_get_metadata(const GXml& xml,
                                                const socklen_t& sock);
    void                   request_shutdown(const socklen_t& sock);
    std::string            get_client_key(const GXml& xml) const;
    int                    get_client_index(const GXml& xml) const;
    std::string            get_response_value(const GXml& xml,
                                              const std::string& name) const;
    void                   get_name_value_pair(const GXmlNode* node,
                                               std::string&    name,
                                               std::string&    value) const;
    std::list<std::string> get_registrations(const GXml& xml) const;
    std::string            get_callback_port(const GXml& xml) const;
    std::string            get_hub_lockfile(void) const;
    void 		           activate_callbacks(std::string method, char cl_id[31]);

    // Low-level methods
    void        create_samp_file(void) const;
    void	    post_string(const std::string& content, const socklen_t& sock) const;
    void 		post_string_callback(const std::string& content) const;
    std::string	random_string(const size_t& length) const;

    // Protected structures
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
    	int           registered;
        sem_t         lock;
        connected_shm metadata[GVOHUB_NB_CLIENTS];
    } clients;

    // Protected members
    std::string m_secret;        //!< Secret Hub key
    std::string m_hub_url;       //!< The XML-RPC endpoint for communication with the hub
    std::string m_hub_host;      //!< Hub host (extracted from XML-RPC endpoint)
    std::string m_hub_port;      //!< Hub port (extracted from XML-RPC endpoint)
    std::string m_version;       //!< The version of the SAMP Standard Profile implemented by the hub
    std::string m_hub_id;        //!< Hub identifier used by the hub when it sends message itself rather than forwarding from others
    socklen_t   m_socket;        //!< Hub socket
    int		    m_cback_socket;  //!< Hub socket to callback clients
    clients*    m_clients;       //!< Structure of registered clients
    void*       m_shmem;         //!< Pointer to shared memory address
    int         m_shmem_handler; //!< Handler to shared memory 
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GVOHub").
 ***************************************************************************/
inline
std::string GVOHub::classname(void) const
{
    return ("GVOHub");
}

#endif /* GVOHUB_HPP */
