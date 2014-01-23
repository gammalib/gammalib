/***************************************************************************
 *                        GLog.i - Information logger                      *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2013 by Juergen Knoedlseder                         *
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
 * @brief Information logger class definition
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
 * @brief Information logger interface definition.
 ***************************************************************************/
class GLog {
public:
    // Constructors and destructors
    GLog(void);
    GLog(const std::string& filename, const bool& clobber = false);
    GLog(const GLog& log);
    virtual ~GLog(void);

    // Methods
    void               clear(void);
    int                size(void) const;
    void               open(const std::string& filename, const bool& clobber = false);
    void               close(void);
    void               flush(const bool& force = false);
    void               date(const bool& flag);
    void               cout(const bool& flag);
    void               cerr(const bool& flag);
    void               name(const std::string& name);
    void               max_size(const int& size);
    void               indent(const int& indent);
    void               chatter(const GChatter& chatter);
    void               header0(const std::string& arg);
    void               header1(const std::string& arg);
    void               header2(const std::string& arg);
    void               header3(const std::string& arg);
    const bool&        date(void) const;
    const bool&        cout(void) const;
    const bool&        cerr(void) const;
    const std::string& name(void) const;
    const int&         max_size(void) const;
    const int&         indent(void) const;
    const GChatter&    chatter(void) const;
    const std::string& filename(void) const;
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
        (*self) << gammalib::parformat(arg);
    }
    void toupper(const std::string& arg) {
        (*self) << gammalib::toupper(arg);
    }
    void tolower(const std::string& arg) {
        (*self) << gammalib::tolower(arg);
    }
    void fill(const std::string& arg, const int& n) {
        (*self) << gammalib::fill(arg, n);
    }
    void left(const std::string& arg, const int& n) {
        (*self) << gammalib::left(arg, n);
    }
    void right(const std::string& arg, const int& n) {
        (*self) << gammalib::right(arg, n);
    }
    void centre(const std::string& arg, const int& n) {
        (*self) << gammalib::centre(arg, n);
    }
    void chatter(const int& chatter) {
        self->chatter(GChatter(chatter));
    }
}
