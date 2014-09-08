/***************************************************************************
 *                      GVOHub.hpp - VO SAMP Hub class                     *
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
 * @file GVOHub.hpp
 * @brief SAMP hub class interface definition
 * @author Thierry Louge
 */

#ifndef GVOHUB_HPP
#define GVOHUB_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include <sys/socket.h>
#include "GBase.hpp"
#include "GXml.hpp"
#include "GXmlNode.hpp"


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
    void        start(void);
    std::string print(const GChatter& chatter = NORMAL) const;
    //void        connect(void);
    //void        disconnect(void);
    //bool        has_hub(void) const;
    //bool        is_connected(void) const;
    //GXml        response(void) const;

protected:
    // Protected methods
    void        init_members(void);
    void        copy_members(const GVOHub& client);
    void        free_members(void);
    void	    init_hub(void);
    void        start_hub(void);
    //bool        find_hub(void);
    //void        connect_to_hub(void);
    void        register_service(void);
    void        unregister(void);
    void 	    handle_request(const socklen_t& sock);
    //void        send_metadata(void);
    
    // Low-level methods
    void        post_string(const std::string& string) const;
    std::string receive_string(void) const;
    std::string get_response_value(const GXml& xml, const std::string& name) const;
    void        get_name_value_pair(const GXmlNode* node, std::string& name, std::string& value) const;
    std::string get_hub_lockfile(void) const;

    // Protected data area
    std::string m_name;        //!< Client name
    std::string m_secret;      //!< Secret Hub key
    std::string m_hub_url;     //!< The XML-RPC endpoint for communication with the hub
    std::string m_hub_host;    //!< Hub host (extracted from XML-RPC endpoint)
    std::string m_hub_port;    //!< Hub port (extracted from XML-RPC endpoint)
    std::string m_version;     //!< The version of the SAMP Standard Profile implemented by the hub
    std::string m_client_key;  //!< Private client key
    std::string m_hub_id;      //!< Hub identifier
    std::string m_client_id;   //!< Client identifier
    int         m_socket;      //!< Hub socket
};

#endif /* GVOHUB_HPP */
