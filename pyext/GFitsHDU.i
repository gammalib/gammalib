/***************************************************************************
 *                   GFitsHDU.i - FITS HDU handling class                  *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2008-2013 by Juergen Knoedlseder                         *
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
 * @file GFitsHDU.i
 * @brief FITS HDU class Python interface definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GFitsHDU.hpp"
#include "GTools.hpp"
%}


/***********************************************************************//**
 * @class GFitsHDU
 *
 * @brief FITS Header Data Unit (HDU) class
 *
 * The HDU is the basic unit of a FITS file. Each HDU consists of a header
 * and a data area. The header is composed of cards and is implemented by
 * the GFitsHeader class. The data are is either an image or a table and
 * is implemented by the abstract GFitsData base class.
 ***************************************************************************/
class GFitsHDU : public GBase {

public:
    // Constructors and destructors
    GFitsHDU(void);
    GFitsHDU(const GFitsHDU& hdu);
    virtual ~GFitsHDU(void);

    // Public enumerators
    enum HDUType {
        HT_IMAGE = 0,
        HT_ASCII_TABLE = 1,
        HT_BIN_TABLE = 2
    };

    // Pure virtual methods
    virtual void      clear(void) = 0;
    virtual GFitsHDU* clone(void) const = 0;
    virtual HDUType   exttype(void) const = 0;

    // Implemented methods
    int                size(void) const;
    const std::string& extname(void) const;
    void               extname(const std::string& extname);
    const int&         extno(void) const;
    void               extno(const int& extno);
    GFitsHeader&       header(void);
    bool               has_card(const int& cardno) const;
    bool               has_card(const std::string& keyname) const;
    GFitsHeaderCard&   card(const int& cardno);
    GFitsHeaderCard&   card(const std::string& keyname);
    void               card(const GFitsHeaderCard& card);
    void               card(const std::string& keyname,
                            const std::string& value,
                            const std::string& comment);
    void               card(const std::string& keyname,
                            const double& value,
                            const std::string& comment);
    void               card(const std::string& keyname,
                            const int& value,
                            const std::string& comment);
    std::string        string(const std::string& keyname) const;
    double             real(const std::string& keyname) const;
    int                integer(const std::string& keyname) const;
};


/***********************************************************************//**
 * @brief GFitsHDU class SWIG extension
 ***************************************************************************/
%extend GFitsHDU {
};
