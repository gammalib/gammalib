/***************************************************************************
 *        GCOMInstChars.i - COMPTEL Instrument Characteristics class       *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2017-2021 by Juergen Knoedlseder                         *
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
 * @file GCOMInstChars.i
 * @brief COMPTEL Instrument Characteristics interface definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GCOMInstChars.hpp"
%}


/***********************************************************************//**
 * @class GCOMInstChars
 *
 * @brief Interface for the COMPTEL Instrument Characteristics class
 ***************************************************************************/
class GCOMInstChars : public GBase {

public:
    // Constructors and destructors
    GCOMInstChars(void);
    GCOMInstChars(const GCOMInstChars& ict);
    GCOMInstChars(const GCaldb& caldb, const std::string& ictname);
    ~GCOMInstChars(void);

    // Methods
    void           clear(void);
    GCOMInstChars* clone(void) const;
    std::string    classname(void) const;
    void           caldb(const GCaldb& caldb);
    const GCaldb&  caldb(void) const;
    void           load(const std::string& ictname);
    double         trans_D1(const double& energy) const;
    double         trans_V1(const double& energy) const;
    double         prob_D1inter(const double& energy) const;
    double         prob_no_multihit(const double& energy) const;
    double         prob_no_selfveto(const double& energy, const double& zenith) const;
    double         trans_D2(const double& energy, const double& phigeo) const;
    double         trans_V23(const double& energy, const double& phigeo) const;
    double         prob_D2inter(const double& energy, const double& phigeo) const;
    double         multi_scatter(const double& energy, const double& phigeo) const;
    double         psd_correction(const double& energy, const double& phigeo) const;
};


/***********************************************************************//**
 * @brief GCOMInstChars class extension
 ***************************************************************************/
%extend GCOMInstChars {
    GCOMInstChars copy() {
        return (*self);
    }
};
