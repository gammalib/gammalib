/***************************************************************************
 *      GModelData.i  -  Abstract virtual data model class python I/F      *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2011 by Jurgen Knodlseder                                *
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
 * @file GModelData.i
 * @brief Abstract data model base class python interface
 * @author J. Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GModelData.hpp"
#include "GTools.hpp"
%}


/***********************************************************************//**
 * @class GModelData
 *
 * @brief Abstract data model base class python interface
 ***************************************************************************/
class GModelData : public GModel {
public:
    // Constructors and destructors
    GModelData(void);
    explicit GModelData(const GXmlElement& xml);
    GModelData(const GModelData& model);
    virtual ~GModelData(void);

    // Pure virtual methods
    virtual void        clear(void) = 0;
    virtual GModelData* clone(void) const = 0;
    virtual std::string type(void) const = 0;
    virtual double      eval(const GEvent& event,
                             const GObservation& obs) const = 0;
    virtual double      eval_gradients(const GEvent& event,
                                       const GObservation& obs) const = 0;
    virtual double      npred(const GEnergy& obsEng, const GTime& obsTime,
                              const GObservation& obs) const = 0;
    virtual void        read(const GXmlElement& xml) = 0;
    virtual void        write(GXmlElement& xml) const = 0;
};


/***********************************************************************//**
 * @brief GModelData class extension
 ***************************************************************************/
%extend GModelData {
};


/***********************************************************************//**
 * @brief GModelData type casts
 ***************************************************************************/
%inline %{
    GModelData* cast_GModelData(GModel* model) {
        return dynamic_cast<GModelData*>(model);
    }
%};

