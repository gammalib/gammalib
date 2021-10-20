/***************************************************************************
 *            GModelData.hpp - Abstract virtual data model class           *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2011-2021 by Juergen Knoedlseder                         *
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
 * @file GModelData.hpp
 * @brief Abstract data model base class interface definition
 * @author Juergen Knoedlseder
 */

#ifndef GMODELDATA_HPP
#define GMODELDATA_HPP

/* __ Includes ___________________________________________________________ */
#include <vector>
#include <string>
#include "GModel.hpp"

/* __ Forward declarations _______________________________________________ */
class GEnergy;
class GTime;
class GEvent;
class GEvents;
class GObservation;
class GRan;
class GXmlElement;
class GVector;
class GMatrixSparse;


/***********************************************************************//**
 * @class GModelData
 *
 * @brief Abstract data model class
 *
 * This abstract base class implements a model of the event distribution in
 * the instrument's data space.
 ***************************************************************************/
class GModelData : public GModel {

public:
    // Constructors and destructors
    GModelData(void);
    explicit GModelData(const GXmlElement& xml);
    GModelData(const GModelData& model);
    virtual ~GModelData(void);

    // Operators
    virtual GModelData& operator=(const GModelData& model);

    // Pure virtual methods
    virtual void        clear(void) = 0;
    virtual GModelData* clone(void) const = 0;
    virtual std::string classname(void) const = 0;
    virtual std::string type(void) const = 0;
    virtual bool        is_constant(void) const = 0;
    virtual double      eval(const GEvent&       event,
                             const GObservation& obs,
                             const bool&         gradients = false) const = 0;
    virtual double      npred(const GEnergy&       obsEng,
                              const GTime&         obsTime,
                              const GPolarization& obsPol,
                              const GObservation&  obs) const = 0;
    virtual GEvents*    mc(const GObservation& obs, GRan& ran) const = 0;
    virtual void        read(const GXmlElement& xml) = 0;
    virtual void        write(GXmlElement& xml) const = 0;
    virtual std::string print(const GChatter& chatter = NORMAL) const = 0;

    // Implemented pure virtual base class methods
    virtual GVector     eval(const GObservation& obs,
                             GMatrixSparse* gradients = NULL) const;

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GModelData& model);
    void free_members(void);
};

#endif /* GMODELDATA_HPP */
