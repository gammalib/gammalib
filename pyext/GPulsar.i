/***************************************************************************
 *                         GPulsar.i - Pulsar class                        *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2022 by Juergen Knoedlseder                              *
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
 * @file GPulsar.i
 * @brief Pulsar class definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GPulsar.hpp"
%}


/***********************************************************************//**
 * @class GPulsar
 *
 * @brief Pulsar class
 ***************************************************************************/
class GPulsar : public GBase {

public:
    // Constructors and destructors
    GPulsar(void);
    GPulsar(const GFilename& filename, const std::string& name = "");
    GPulsar(const GPulsar& pulsar);
    virtual ~GPulsar(void);

    // Implemented pure virtual base class methods
    virtual void        clear(void);
    virtual GPulsar*    clone(void) const;
    virtual std::string classname(void) const;

    // Other methods
    int                     size(void) const;
    bool                    is_empty(void) const;
    const std::string&      name(void) const;
    void                    name(const std::string& name);
    const GPulsarEphemeris& ephemeris(const GTime& time) const;
    GGti                    validity(void) const;
    void                    load(const GFilename& filename,
                                 const std::string& name = "");
};


/***********************************************************************//**
 * @brief GPulsar class extension
 ***************************************************************************/
%extend GPulsar {
    GPulsar copy() {
        return (*self);
    }
};
