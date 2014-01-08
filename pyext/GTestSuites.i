/***************************************************************************
 *              GTestSuites.i - Test suite container class                 *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2012-2013 Jean-Baptiste Cayrou                           *
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
 * @file GTestSuites.i
 * @brief Test suite container class interface defintion
 * @author Jean-Baptiste Cayrou
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GTestSuites.hpp"
#include "GTools.hpp"
%}


/***********************************************************************//**
 * @class GTestSuites
 *
 * @brief Test suites container class
 ***************************************************************************/
class GTestSuites : public GContainer {
    public:

    // Constructors and destructors
    GTestSuites(void);
    GTestSuites(const GTestSuites& suites);
    GTestSuites(const std::string& name);
    virtual ~GTestSuites(void);

    // Methods
    void               clear(void);
    GTestSuites*       clone(void) const;
    GTestSuite*        at(const int& index);
    int                size(void) const;
    bool               is_empty(void) const;
    GTestSuite*        set(const int& index, const GTestSuite& suite);
    GTestSuite*        append(const GTestSuite& suite);
    GTestSuite*        insert(const int& index, const GTestSuite& suite);
    void               remove(const int& index);
    void               reserve(const int& num);
    void               extend(const GTestSuites& suites);
    const std::string& name(void) const;
    void               name(const std::string& name);
    void               cout(const bool& flag);
    int                errors(void) const;
    int                failures(void) const;
    int                tests(void) const;
    const time_t&      timestamp(void) const;
    bool               run(void);
    void               save(const std::string& filename) const;
};


/***********************************************************************//**
 * @brief GTestSuites class extension
 ***************************************************************************/
%extend GTestSuites {
    GTestSuite* __getitem__(const int& index) {
        if (index >= 0 && index < self->size()) {
            return (*self)[index];
        }
        else {
            throw GException::out_of_range("__getitem__(int)", index, self->size());
        }
    }
    void __setitem__(const int& index, const GTestSuite& suite) {
        if (index >= 0 && index < self->size()) {
            self->set(index, suite);
            return;
        }
        else {
            throw GException::out_of_range("__setitem__(int)", index, self->size());
        }
    }
    GTestSuites copy() {
        return (*self);
    }
};
