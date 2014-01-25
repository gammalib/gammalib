/***************************************************************************
 *                      GVOClient.i - VO client class                      *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2014 by Juergen Knoedlseder                              *
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
 * @file GVOClient.i
 * @brief VO client class interface definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GVOClient.hpp"
%}


/***********************************************************************//**
 * @class GVOClient
 *
 * @brief VO client class
 ***************************************************************************/
class GVOClient : public GBase {
public:
    // Constructors and destructors
    GVOClient(void);
    GVOClient(const GVOClient& client);
    virtual ~GVOClient(void);

    // Methods
    void       clear(void);
    GVOClient* clone(void) const;
    void       connect(void);
    void       disconnect(void);
    bool       has_hub(void) const;
    bool       is_connected(void) const;
    GXml       response(void) const;
};


/***********************************************************************//**
 * @brief GVOClient class extension
 ***************************************************************************/
%extend GVOClient {
    GVOClient copy() {
        return (*self);
    }
};
