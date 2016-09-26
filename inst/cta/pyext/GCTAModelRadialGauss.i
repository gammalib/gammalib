/***************************************************************************
 *            GCTAModelRadialGauss.i - Radial Gaussian model class         *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2011-2016 by Juergen Knoedlseder                         *
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
 * @file GCTAModelRadialGauss.i
 * @brief Radial Gaussian model class interface definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GCTAModelRadialGauss.hpp"
#include "GTools.hpp"
%}


/***********************************************************************//**
 * @class GCTAModelRadialGauss
 *
 * @brief Radial Gaussian CTA model class
 ***************************************************************************/
class GCTAModelRadialGauss : public GCTAModelRadial {
public:
    // Constructors and destructors
    GCTAModelRadialGauss(void);
    explicit GCTAModelRadialGauss(const double& sigma);
    explicit GCTAModelRadialGauss(const GXmlElement& xml);
    GCTAModelRadialGauss(const GCTAModelRadialGauss& model);
    virtual ~GCTAModelRadialGauss(void);

    // Implemented pure virtual methods
    virtual void                  clear(void);
    virtual GCTAModelRadialGauss* clone(void) const;
    virtual std::string           classname(void) const;
    virtual std::string           type(void) const;
    virtual double                eval(const double& offset,
                                       const bool& gradients = false) const;
    virtual GCTAInstDir           mc(const GCTAInstDir& dir, GRan& ran) const;
    virtual double                omega(void) const;
    virtual void                  read(const GXmlElement& xml);
    virtual void                  write(GXmlElement& xml) const;

    // Other methods
    double sigma(void) const;
    void   sigma(const double& sigma);
};


/***********************************************************************//**
 * @brief GCTAModelRadialGauss class extension
 ***************************************************************************/
%extend GCTAModelRadialGauss {
};
