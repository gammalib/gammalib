/***************************************************************************
 *              GTestCase.i - Test case class Python interface             *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2012 Jean-Baptiste Cayrou                                *
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
 * @file GTestCase.i
 * @brief Test case class Python interface defintion
 * @author Jean-Baptiste Cayrou
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GTestCase.hpp"
#include "GTools.hpp"
%}


/***********************************************************************//**
 * @class GTestCase
 *
 * @brief Test case Python interface defintion
 ***************************************************************************/
class GTestCase {

public:
    // public enumerators
    enum ErrorKind {
        FAIL_TEST,
        ERROR_TEST
    };

    // Constructors and destructors
    GTestCase(void);
    GTestCase(const GTestCase& test);
    GTestCase(ErrorKind kind, const std::string& name = "");
    virtual ~GTestCase(void);

    // Methods
    std::string name(void) const;
    void        name(const std::string& name);
    std::string message(void) const;
    void        message(const std::string& message);
    std::string type(void) const;
    void        type(const std::string& type);
    ErrorKind   kind(void) const;
    void        kind(ErrorKind kind);
    bool        passed(void) const;
    void        passed(const bool& passed);
    double      duration(void) const;
    void        duration(const double& duration);
};


/***********************************************************************//**
 * @brief GTestCase class extension
 ***************************************************************************/
%extend GTestCase {
    char *__str__() {
        return tochar(self->print());
    }
    GTestCase copy() {
        return (*self);
    }
};
