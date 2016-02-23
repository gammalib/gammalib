/***************************************************************************
 *                      GVOHub.hpp - VO SAMP Hub class                     *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2014-2016 by Thierry Louge                               *
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

/* __ Includes ___________________________________________________________ */
#include <string>
#include <vector>
#include <sys/socket.h>
#include "GBase.hpp"

/* __ Forward declarations _______________________________________________ */
class GXml;
class GXmlNode;


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
    void                     init_members(void);
    void                     copy_members(const GVOHub& client);
    void                     free_members(void);
    void                     start_hub(void);
    void                     handle_request(const socklen_t& sock);
    void                     request_ping(const socklen_t& sock);
    void                     request_register(const GXml&      xml,
                                              const socklen_t& sock);
    void                     request_unregister(const GXml&      xml,
                                                const socklen_t& sock);
    void                     request_declare_metadata(const GXml&      xml,
                                                      const socklen_t& sock);
    void                     request_declare_subscriptions(const GXml&      xml,
                                                           const socklen_t& sock);
    void                     request_set_xml_rpc_callback(const GXml&      xml,
                                                          const socklen_t& sock);
    void                     request_get_subscriptions(const GXml&      xml,
                                                       const socklen_t& sock);
    void                     request_get_registered_clients(const GXml&      xml,
                                                            const socklen_t& sock);
    void                     request_get_subscribed_clients(const GXml&      xml,
                                                            const socklen_t& sock);
    void                     request_get_metadata(const GXml&      xml,
                                                  const socklen_t& sock);
    void                     request_notify_all(const GXml&      xml,
                                                const socklen_t& sock);
    void                     request_shutdown(const socklen_t& sock);
    std::string              get_client_key(const GXml& xml) const;
    int                      get_client_index(const GXml& xml) const;
    int                      get_client_index(const std::string& reference) const;
    std::string              get_response_value(const GXmlNode*    node,
                                                const std::string& name) const;
    void                     get_name_value_pair(const GXmlNode* node,
                                                 std::string&    name,
                                                 std::string&    value) const;
    std::vector<std::string> get_subscriptions(const GXml& xml) const;
    std::string              get_callback_url(const GXml& xml) const;
    std::string              get_hub_lockfile(void) const;
    std::string              get_mtype(const GXml& xml) const;
    std::string              get_destination(const GXml& xml) const;

    // Client structure
    struct client {
        std::string              name;
        std::string              private_key;
        std::string              reference;
        std::string              description;
        std::string              icon;
        std::string              documentation;
        std::string              affiliation;
        std::string              author_name;
        std::string              email;
        std::string              homepage;
        std::string              url;
        std::vector<std::string> subscriptions;
    };

    // Low-level methods
    void        create_samp_file(void) const;
    void        post_samp_ok(const socklen_t& sock) const;
    void        post_samp_void(const socklen_t& sock) const;
    void	    post_string(const std::string& content,
                            const socklen_t&   sock) const;
    void        notify(const std::string& url,
                       const std::string& notification) const;
    void        notify_register(const client&      client,
                                const std::string& reference);
    void        notify_unregister(const client&      client,
                                  const std::string& reference);
    void        notify_metadata(const client&      client,
                                const std::string& reference);
    void        notify_image_load(const client& client,
                                  const GXml&   xml);
    std::string	random_string(const size_t& length) const;

    // Protected members
    std::string         m_secret;    //!< Secret Hub key
    std::string         m_hub_url;   //!< The XML-RPC endpoint for communication with the hub
    std::string         m_hub_host;  //!< Hub host
    std::string         m_hub_port;  //!< Hub port
    std::string         m_hub_path;  //!< Hub path
    std::string         m_version;   //!< The version of the SAMP Standard Profile implemented by the hub
    std::string         m_hub_id;    //!< Hub identifier used by the hub when it sends message itself rather than forwarding from others
    socklen_t           m_socket;    //!< Hub socket
    bool                m_shutdown;  //!< Shutdown request
    std::vector<client> m_clients;   //!< Clients
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
