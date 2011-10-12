/***************************************************************************
 *       GCTAModelRadialPolynom.i  -  Radial Polynom CTA model class       *
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
 * @file GCTAModelRadialPolynom.i
 * @brief Radial Polynom model class Python interface definition
 * @author J. Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GCTAModelRadialPolynom.hpp"
#include "GTools.hpp"
%}


/***********************************************************************//**
 * @class GCTAModelRadialPolynom
 *
 * @brief Radial Polynom CTA model class
 ***************************************************************************/
class GCTAModelRadialPolynom : public GCTAModelRadial {
public:
    // Constructors and destructors
    GCTAModelRadialPolynom(void);
    explicit GCTAModelRadialPolynom(const std::vector<double>& coeffs);
    explicit GCTAModelRadialPolynom(const GXmlElement& xml);
    GCTAModelRadialPolynom(const GCTAModelRadialPolynom& model);
    virtual ~GCTAModelRadialPolynom(void);

    // Implemented pure virtual methods
    virtual void                    clear(void);
    virtual GCTAModelRadialPolynom* clone(void) const;
    virtual std::string             type(void) const;
    virtual double                  eval(const double& offset) const;
    virtual double                  eval_gradients(const double& offset) const;
    virtual GCTAInstDir             mc(const GCTAInstDir& dir, GRan& ran) const;
    virtual double                  omega(void) const;
    virtual void                    read(const GXmlElement& xml);
    virtual void                    write(GXmlElement& xml) const;

    // Other methods
    int    size(void) const { return m_coeffs.size(); }
    //double coeff(void) const;
    //void   coeff(const double& value);
};


/***********************************************************************//**
 * @brief GCTAModelRadialPolynom class extension
 ***************************************************************************/
%extend GCTAModelRadialPolynom {
    char *__str__() {
        return tochar(self->print());
    }
};


/***********************************************************************//**
 * @brief GCTAModelRadialPolynom type casts
 ***************************************************************************/
%inline %{
    GCTAModelRadialPolynom* cast_GCTAModelRadialPolynom(GCTAModelRadial* model) {
        GCTAModelRadialPolynom* cast = dynamic_cast<GCTAModelRadialPolynom*>(model);
        if (cast == NULL) {
            throw GException::bad_type("cast_GCTAModelRadialPolynom(GCTAModelRadial* model)",
                                       "GCTAModelRadial not of type GCTAModelRadialPolynom");
        }
        return cast;
    }
%}
