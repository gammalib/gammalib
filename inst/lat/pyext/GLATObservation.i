/***************************************************************************
 *        GLATObservation.i  -  LAT Observation class SWIG interface       *
 * ----------------------------------------------------------------------- *
 *  copyright : (C) 2008-2010 by Jurgen Knodlseder                         *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GLATObservation.i
 * @brief GLATObservation class SWIG file.
 * @author J. Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GLATObservation.hpp"
%}
%include stl.i

/***********************************************************************//**
 * @class GLATObservation
 *
 * @brief Interface for the LAT observation classes.
 ***************************************************************************/
class GLATObservation : public GObservation {
public:
    // Constructors and destructors
    GLATObservation();
    GLATObservation(const GLATObservation& obs);
    virtual ~GLATObservation();

    // Implement pure virtual methods
    void             clear(void);
    GLATObservation* clone(void) const;
    GLATResponse*    response(void) const;
    GLATPointing*    pointing(const GTime& time) const;
    std::string      instrument(void) const;

    // Other methods
    void        load_unbinned(const std::string& ft1name, const std::string& ft2name,
                              const std::string& ltcube_name);
    void        load_binned(const std::string& cntmap_name, const std::string& expmap_name,
                            const std::string& ltcube_name);
    void        response(const std::string& irfname, std::string caldb = "");
    GLATLtCube* ltcube(void) const;
};


/***********************************************************************//**
 * @brief GLATObservation class extension
 ***************************************************************************/
%extend GLATObservation {
    GLATObservation(const GObservation& obs) {
        GLATObservation* lat = (GLATObservation*)&obs;
        return lat;
    }
    char *__str__() {
        static std::string result = self->print();
        return ((char*)result.c_str());
    }
    GLATObservation copy() {
        return (*self);
    }
};
