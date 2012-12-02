/***************************************************************************
 *        GCOMModelDRBFitting.i  -  COMPTEL DRB model fitting class        *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2012 by Juergen Knoedlseder                              *
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
 * @file GCOMModelDRBFitting.hpp
 * @brief COMPTEL DRB model fitting class interface definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GCOMModelDRBFitting.hpp"
%}


/***********************************************************************//**
 * @class GCOMModelDRBFitting
 *
 * @brief COMPTEL DRB model fitting class
 ***************************************************************************/
class GCOMModelDRBFitting : public GModelData {
public:
    // Constructors and destructors
    GCOMModelDRBFitting(void);
    explicit GCOMModelDRBFitting(const GXmlElement& xml);
    GCOMModelDRBFitting(const GCOMModelDRBFitting& model);
    virtual ~GCOMModelDRBFitting(void);

    // Implemented pure virtual methods
    virtual void                 clear(void);
    virtual GCOMModelDRBFitting* clone(void) const;
    virtual std::string          type(void) const;
    virtual double               eval(const GEvent& event,
                                      const GObservation& obs) const;
    virtual double               eval_gradients(const GEvent& event,
                                                const GObservation& obs) const;
    virtual double               npred(const GEnergy& obsEng, const GTime& obsTime,
                                       const GObservation& obs) const;
    virtual GCOMEventCube*       mc(const GObservation& obs, GRan& ran) const;
    virtual void                 read(const GXmlElement& xml);
    virtual void                 write(GXmlElement& xml) const;
};


/***********************************************************************//**
 * @brief GCOMModelDRBFitting class extension
 ***************************************************************************/
%extend GCOMModelDRBFitting {
    GCOMModelDRBFitting copy() {
        return (*self);
    }
};
