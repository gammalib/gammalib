/***************************************************************************
 *          GCTAModelRadialProfile.i - Radial Profile model class          *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2011-2013 by Juergen Knoedlseder                         *
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
 * @file GCTAModelRadialProfile.i
 * @brief Radial Profile model class Python interface definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GCTAModelRadialProfile.hpp"
#include "GTools.hpp"
%}


/***********************************************************************//**
 * @class GCTAModelRadialProfile
 *
 * @brief Radial Profile CTA model class
 ***************************************************************************/
class GCTAModelRadialProfile : public GCTAModelRadial {
public:
    // Constructors and destructors
    GCTAModelRadialProfile(void);
    explicit GCTAModelRadialProfile(const double& width, const double& core,
                                    const double& tail);
    explicit GCTAModelRadialProfile(const GXmlElement& xml);
    GCTAModelRadialProfile(const GCTAModelRadialProfile& model);
    virtual ~GCTAModelRadialProfile(void);

    // Implemented pure virtual methods
    virtual void                    clear(void);
    virtual GCTAModelRadialProfile* clone(void) const;
    virtual std::string             type(void) const;
    virtual double                  eval(const double& offset) const;
    virtual double                  eval_gradients(const double& offset) const;
    virtual GCTAInstDir             mc(const GCTAInstDir& dir, GRan& ran) const;
    virtual double                  omega(void) const;
    virtual void                    read(const GXmlElement& xml);
    virtual void                    write(GXmlElement& xml) const;

    // Other methods
    double width(void) const;
    double core(void) const;
    double tail(void) const;
    void   width(const double& width);
    void   core(const double& core);
    void   tail(const double& tail);
};


/***********************************************************************//**
 * @brief GCTAModelRadialProfile class extension
 ***************************************************************************/
%extend GCTAModelRadialProfile {
};
