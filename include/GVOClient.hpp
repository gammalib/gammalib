/***************************************************************************
 *                     GVOClient.hpp - VO client class                     *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2013-2016 by Juergen Knoedlseder                         *
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

/* __ Forward declarations _______________________________________________ */
class GFitsHDU;
class GXml;
class GXmlNode;
class GVOTable;


/***********************************************************************//**
 * @class GVOClient
 *
 * @brief VO client class
 *
 * This class implements a client for the Virtual Observatory. Upon
 * construction of an instance of the class, the client will search for a
 * VO Hub. The has_hub() method signals whether VO Hub information was
 * found.
 *
 * The connect() method will connect the client to the VO Hub. If no VO
 * Hub has been found or if the VO Hub that was found is not alive the
 * client will start an own VO Hub using the GVOHub class. The
 * is_connected() method signals whether the client is connected to a
 * Hub.
 * 
 * The disconnect() method disconnects a client from the VO Hub.
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
    bool        ping_hub(void) const;
    void        shutdown_hub(void) const;
    GXml        execute(const std::string& request) const;
    void        publish(const GFitsHDU& hdu);
    void        publish(const GVOTable& votable);
    std::string print(const GChatter& chatter = NORMAL) const;
    
protected:
    // Protected methods
    void        init_members(void);
    void        copy_members(const GVOClient& client);
    void        free_members(void);
    bool        find_hub(void);
    bool        require_hub(void);
    void        register_to_hub(void);
    void        unregister_from_hub(void);
    void        send_metadata(void) const;
    
    // Low-level methods
    void        connect_to_hub(void) const;
    void        post_string(const std::string& string) const;
    std::string receive_string(void) const;
    std::string get_response_value(const GXml& xml, const std::string& name) const;
    void        get_name_value_pair(const GXmlNode* node, std::string& name, std::string& value) const;
    std::string get_hub_lockfile(void) const;
    bool        response_is_valid(const GXml& xml) const;
    int         response_error_code(const GXml& xml) const;
    std::string response_error_message(const GXml& xml) const;
    
    // Protected data area
    std::string m_name;        //!< Client name
    std::string m_secret;      //!< Secret Hub key
    std::string m_hub_url;     //!< The XML-RPC endpoint for communication with the hub
    std::string m_hub_host;    //!< Hub host (extracted from XML-RPC endpoint)
    std::string m_hub_port;    //!< Hub port (extracted from XML-RPC endpoint)
    std::string m_hub_path;    //!< Hub path (extracted from XML-RPC endpoint)
    std::string m_version;     //!< The version of the SAMP Standard Profile implemented by the hub
    std::string m_client_key;  //!< Private client key
    std::string m_hub_id;      //!< Hub identifier
    std::string m_client_id;   //!< Client identifier
    mutable int m_socket;      //!< Hub socket
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
