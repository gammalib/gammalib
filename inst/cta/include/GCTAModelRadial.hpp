/***************************************************************************
 *          GCTAModelRadial.hpp - Radial model abstract base class         *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2011-2018 by Juergen Knoedlseder                         *
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
 * @file GCTAModelRadial.hpp
 * @brief Abstract radial acceptance model class interface definition
 * @author Juergen Knoedlseder
 */

#ifndef GCTAMODELRADIAL_HPP
#define GCTAMODELRADIAL_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GBase.hpp"
#include "GModelPar.hpp"
#include "GXmlElement.hpp"
#include "GRan.hpp"
#include "GCTAInstDir.hpp"


/***********************************************************************//**
 * @class GCTAModelRadial
 *
 * @brief Abstract radial acceptance model class
 *
 * This class implements the radial component of the CTA radial acceptance
 * model.
 ***************************************************************************/
class GCTAModelRadial : public GBase {

public:
    // Constructors and destructors
    GCTAModelRadial(void);
    GCTAModelRadial(const GCTAModelRadial& model);
    virtual ~GCTAModelRadial(void);

    // Operators
    virtual GCTAModelRadial& operator=(const GCTAModelRadial& model);
    virtual GModelPar&       operator[](const int& index);
    virtual const GModelPar& operator[](const int& index) const;
    virtual GModelPar&       operator[](const std::string& name);
    virtual const GModelPar& operator[](const std::string& name) const;

    // Pure virtual methods
    virtual void             clear(void) = 0;
    virtual GCTAModelRadial* clone(void) const = 0;
    virtual std::string      classname(void) const = 0;
    virtual std::string      type(void) const = 0;
    virtual double           eval(const double& offset,
                                  const bool&   gradients = false) const = 0;
    virtual GCTAInstDir      mc(const GCTAInstDir& dir, GRan& ran) const = 0;
    virtual double           omega(void) const = 0;
    virtual void             read(const GXmlElement& xml) = 0;
    virtual void             write(GXmlElement& xml) const = 0;
    virtual std::string      print(const GChatter& chatter = NORMAL) const = 0;

    // Methods
    int size(void) const;

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GCTAModelRadial& model);
    void free_members(void);

    // Proteced members
    std::vector<GModelPar*> m_pars;  //!< Parameter pointers
};


/***********************************************************************//**
 * @brief Return number of model parameters
 *
 * @return Number of model parameters.
 ***************************************************************************/
inline
int GCTAModelRadial::size(void) const
{
    return ((int)m_pars.size());
}

#endif /* GCTAMODELRADIAL_HPP */
