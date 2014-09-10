/***************************************************************************
 *                     GVOClient.hpp - VO client class                     *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2013-2014 by Juergen Knoedlseder                         *
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
 * @file GVOClient.hpp
 * @brief VO client class interface definition
 * @author Juergen Knoedlseder
 */

#ifndef GVOCLIENT_HPP
#define GVOCLIENT_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GBase.hpp"
#include "GXml.hpp"
#include "GXmlNode.hpp"


/***********************************************************************//**
 * @class GVOClient
 *
 * @brief VO client class
 *
 * This class implements a client for the Virtual Observatory. Upon
 * construction of an instance of the class, the client will search for a
 * SAMP Hub. If no such Hub exists, a GammaLib owned Hub will be started
 * using the GVOHub class. Once a Hub is running, the client will retrieve
 * its information and store it.
 *
 * The connect() and disconnect() methods exist to connect or disconnect
 * from the Hub. Connecting means opening a socket to the Hub and registering
 * the client to the Hub.
 ***************************************************************************/
class GVOClient : public GBase {

public:
    // Constructors and destructors
    GVOClient(void);
    GVOClient(const GVOClient& client);
    virtual ~GVOClient(void);

    // Operators
    GVOClient& operator=(const GVOClient& client);

    // Methods
    void        clear(void);
    GVOClient*  clone(void) const;
    std::string classname(void) const;
    void        connect(void);
    void        disconnect(void);
    bool        has_hub(void) const;
    bool        is_connected(void) const;
    GXml        response(void) const;
    std::string print(const GChatter& chatter = NORMAL) const;

protected:
    // Protected methods
    void        init_members(void);
    void        copy_members(const GVOClient& client);
    void        free_members(void);
    bool        find_hub(void);
    void        connect_to_hub(void);
    void        register_to_hub(void);
    void        unregister_from_hub(void);
    void        send_metadata(void);
    
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


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GVOClient").
 ***************************************************************************/
inline
std::string GVOClient::classname(void) const
{
    return ("GVOClient");
}

#endif /* GVOCLIENT_HPP */
