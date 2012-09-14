/***************************************************************************
 *                        GLog.i - Information logger                      *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2011 by Jurgen Knodlseder                           *
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
 * @file GLog.i
 * @brief Information logger python bindings
 * @author Jurgen Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GLog.hpp"
#include "GTools.hpp"
%}


/***********************************************************************//**
 * @class GLog
 *
 * @brief Information logger interface defintion.
 ***************************************************************************/
class GLog {
public:
    // Constructors and destructors
    GLog(void);
    GLog(const std::string& filename, bool clobber = false);
    GLog(const GLog& log);
    ~GLog(void);

    // Methods
    void        clear(void);
    int         size(void) const;
    void        open(const std::string& filename, bool clobber = false);
    void        close(void);
    void        date(bool flag);
    void        cout(bool flag);
    void        cerr(bool flag);
    void        name(const std::string& name);
    void        max_size(int size);
    void        indent(int indent);
    void        header0(const std::string& arg);
    void        header1(const std::string& arg);
    void        header2(const std::string& arg);
    void        header3(const std::string& arg);
    bool        date(void) const;
    bool        cout(void) const;
    bool        cerr(void) const;
    std::string name(void) const;
    int         max_size(void) const;
    int         indent(void) const;
    std::string filename(void) const;
};


/***********************************************************************//**
 * @brief GLog class SWIG extension
 ***************************************************************************/
%extend GLog {
    GLog copy() {
        return (*self);
    }
    void __call__(GLog& log) {
        (*self) << log;
    }
    void __call__(const std::string& arg) {
        (*self) << arg;
    }
    void __call__(const bool& value) {
        (*self) << value;
    }
    void __call__(const int& value) {
        (*self) << value;
    }
    void __call__(const double& value) {
        (*self) << value;
    }
    void parformat(const std::string& arg) {
        (*self) << parformat(arg);
    }
    void toupper(const std::string& arg) {
        (*self) << toupper(arg);
    }
    void tolower(const std::string& arg) {
        (*self) << tolower(arg);
    }
    void fill(const std::string& arg, int n) {
        (*self) << fill(arg, n);
    }
    void left(const std::string& arg, int n) {
        (*self) << left(arg, n);
    }
    void right(const std::string& arg, int n) {
        (*self) << right(arg, n);
    }
    void center(const std::string& arg, int n) {
        (*self) << center(arg, n);
    }
}
