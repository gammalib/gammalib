/***************************************************************************
 *                      GVOApp.hpp - VO SAMP Hub class                     *
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
 * @file GVOApp.hpp
 * @brief SAMP hub class interface definition
 * @author Thierry Louge
 */

#ifndef GVOAPP_HPP
#define GVOAPP_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include <list>
#include "GBase.hpp"
#include "GXml.hpp"
#include "GXmlNode.hpp"


/***********************************************************************//**
 * @class GVOApp
 *
 * @brief VO SAMP Hub class
 *
 * This class implements a SAMP hub for exchanges 
 * through VO-compatible applications.
 ***************************************************************************/
class GVOApp : public GBase {

public:
    // Constructors and destructors
    GVOApp(void);
    GVOApp(const GVOApp& hub);
    virtual ~GVOApp(void);

    // Operators
    GVOApp& operator=(const GVOApp& hub);

    // Methods
    void        clear(void);
    GVOApp*     clone(void) const;
    void        start(void);
    std::string print(const GChatter& chatter = NORMAL) const;

protected:
    // Protected methods
    void        init_members(void);
    void        copy_members(const GVOApp& client);
    void        free_members(void);
    void        init_app(void);
    void        start_app(void);
    void        post_string(const std::string& string) const;
    std::string receive_string(void) const;
    std::string get_response_value(const GXml& xml, const std::string& name) const;
    void        get_name_value_pair(const GXmlNode* node, std::string& name, std::string& value) const;

    // Protected data area
    std::string             m_name;           //!< Client name
    std::string             m_id;             //!< Id assigned to client by the hub
    std::list <std::string> registered_calls; //!< Calls for which is registered the application.
    
    
};

#endif /* GVOAPP_HPP */
