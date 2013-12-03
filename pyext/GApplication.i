/***************************************************************************
 *              GApplication.i - GammaLib application base class           *
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
 * @file GApplication.i
 * @brief GammaLib application base class Python interface defintion
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GApplication.hpp"
#include "GTools.hpp"
%}

// Include (int ARGC, char **ARGV) typemap to allow passing command line
// arguments to GApplication constructor
%include "argcargv.i"


/***********************************************************************//**
 * @class GApplication
 *
 * @brief GammaLib application Python interface defintion.
 ***************************************************************************/
class GApplication : public GBase {
public:
    // Constructors and destructors
    GApplication(void);
    GApplication(const std::string& name, const std::string& version);
    GApplication(const std::string& name, const std::string& version,
                 int ARGC, char **ARGV);
    GApplication(const GApplication& app);
    ~GApplication(void);

    // Methods
    void               clear(void);
    GApplication*      clone(void) const;
    const std::string& name(void) const;
    const std::string& version(void) const;
    double             telapse(void) const;
    double             celapse(void) const;
    void               logFileOpen(const bool& clobber = true);
    bool               logTerse(void) const;
    bool               logNormal(void) const;
    bool               logExplicit(void) const;
    bool               logVerbose(void) const;
    bool               logDebug(void) const;
    bool               clobber(void) const;
    bool               haspar(const std::string& name) const;
    const std::string& par_filename(void) const;
    const std::string& log_filename(void) const;
    void               log_header(void);
    void               log_trailer(void);
    void               log_parameters(void);

    // Public members
    GLog log;   //!< Application logger
};


/***********************************************************************//**
 * @brief GApplication class extension
 ***************************************************************************/
%extend GApplication {
    GPar& __getitem__(const std::string& name) {
        return (*self)[name];
    }
    void __setitem__(const std::string& name, const GPar& val) {
        (*self)[name] = val;
        return;
    }
    void __setitem__(const std::string& name, const int& val) {
        (*self)[name].integer(val);
        return;
    }
    void __setitem__(const std::string& name, const double& val) {
        (*self)[name].real(val);
        return;
    }
    void __setitem__(const std::string& name, const std::string& val) {
        (*self)[name].value(val);
        return;
    }
    GApplication copy() {
        return (*self);
    }
};
