/***************************************************************************
 *    GCTABackgroundPerfTable.i - CTA performance table background class   *
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
 * @file GCTABackgroundPerfTable.i
 * @brief CTA performance table background class definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GCTABackgroundPerfTable.hpp"
%}


/***********************************************************************//**
 * @class GCTABackgroundPerfTable
 *
 * @brief CTA performance table background class
 ***************************************************************************/
class GCTABackgroundPerfTable : public GCTABackground {
public:
    // Constructors and destructors
    GCTABackgroundPerfTable(void);
    explicit GCTABackgroundPerfTable(const std::string& filename);
    GCTABackgroundPerfTable(const GCTABackgroundPerfTable& bgd);
    virtual ~GCTABackgroundPerfTable(void);

    // Implemented pure virtual operators
    virtual double operator()(const double& logE, 
                              const double& detx, 
                              const double& dety,
                              const bool&   etrue = false) const;

    // Implemented pure virtual methods
    void                       clear(void);
    GCTABackgroundPerfTable*   clone(void) const;
    void                       load(const std::string& filename);
    std::string                filename(void) const;
    GCTAInstDir                mc(const GEnergy& energy,
                                  const GTime& time,
                                  GRan& ran) const;
    const GModelSpectralNodes& spectrum(void) const;
    std::string                print(const GChatter& chatter = NORMAL) const;

    // Methods
    int           size(void) const;
    void          sigma(const double& sigma);
    const double& sigma(void) const;
};


/***********************************************************************//**
 * @brief GCTABackgroundPerfTable class extension
 ***************************************************************************/
%extend GCTABackgroundPerfTable {
    GCTABackgroundPerfTable copy() {
        return (*self);
    }
};
