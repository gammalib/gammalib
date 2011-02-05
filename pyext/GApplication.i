/***************************************************************************
 *   GApplication.i - GammaLib application base class Python interface     *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2011 by Jurgen Knodlseder                           *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GApplication.i
 * @brief GammaLib application base class Python interface defintion
 * @author Jurgen Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GApplication.hpp"
%}


/***********************************************************************//**
 * @class GApplication
 *
 * @brief GammaLib application Python interface defintion.
 ***************************************************************************/
class GApplication {
public:
    // Constructors and destructors
    GApplication(void);
    GApplication(const std::string& name, const std::string& version);
    GApplication(const std::string& name, const std::string& version,
                 int argc, char* argv[]);
    GApplication(const GApplication& app);
    ~GApplication(void);

    // Methods
    std::string name(void) const;
    std::string version(void) const;
    double      telapse(void) const;
    double      celapse(void) const;
    GPar*       par(const std::string& name);
    const GPar* par(const std::string& name) const;
    bool        logTerse(void) const;
    bool        logNormal(void) const;
    bool        logExplicit(void) const;
    bool        logVerbose(void) const;
    bool        logDebug(void) const;
    bool        clobber(void) const;
    bool        haspar(const std::string& name) const;
};


/***********************************************************************//**
 * @brief GApplication class extension
 ***************************************************************************/
%extend GApplication {
    char *__str__() {
        return tochar(self->print());
    }
    GPar& __getitem__(const std::string& name) {
        return (*self)(name);
    }
    void __setitem__(const std::string& name, const GPar& val) {
        (*self)(name) = val;
    }
    GApplication copy() {
        return (*self);
    }
}
