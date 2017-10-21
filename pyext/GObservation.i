/***************************************************************************
 *             GObservation.i - Abstract observation base class            *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2008-2017 by Juergen Knoedlseder                         *
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
 * @file GObservation.i
 * @brief Abstract observation base class interface definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GObservation.hpp"
#include "GEvents.hpp"
#include "GEventList.hpp"
#include "GEventCube.hpp"
#include "GTools.hpp"
%}

/* __ Typemaps ___________________________________________________________ */
/*
 * Convert double output argument from and into a Python float object. This
 * typemap is needed for the likelihood() method.
 */
%typemap(in) double *npred (double temp) {
    temp = PyFloat_AsDouble($input);
    $1   = &temp;
}
%typemap(argout) double *npred {
    $result = PyFloat_FromDouble(*$1);
}


/***********************************************************************//**
 * @class GObservation
 *
 * @brief Abstract observation base class
 ***************************************************************************/
class GObservation : public GBase {
public:
    // Constructors and destructors
    GObservation(void);
    GObservation(const GObservation& obs);
    virtual ~GObservation(void);

    // Pure virtual methods
    virtual void             clear(void) = 0;
    virtual GObservation*    clone(void) const = 0;
    virtual std::string      classname(void) const = 0;
    virtual void             response(const GResponse& rsp) = 0;
    virtual const GResponse* response(void) const = 0;
    virtual std::string      instrument(void) const = 0;
    virtual double           ontime(void) const = 0;
    virtual double           livetime(void) const = 0;
    virtual double           deadc(const GTime& time) const = 0;
    virtual void             read(const GXmlElement& xml) = 0;
    virtual void             write(GXmlElement& xml) const = 0;

    // Virtual methods
    virtual GEvents*         events(void);
    virtual void             events(const GEvents& events);
    virtual double           likelihood(const GModels& models,
                                        GVector*       gradient,
                                        GMatrixSparse* curvature,
                                        double*        npred) const;
    virtual double           model(const GModels& models,
                                   const GEvent&  event,
                                   GVector*       gradient = NULL) const;
    virtual int              nobserved(void) const;
    virtual double           npred(const GModels& models,
                                   GVector*       gradient = NULL) const;
    virtual double           npred(const GModel& model) const;
    virtual double           model_grad(const GModel&    model,
                                        const GModelPar& par,
                                        const GEvent&    event) const;
    virtual double           npred_grad(const GModel&    model,
                                        const GModelPar& par) const;

    // Implemented methods
    void               name(const std::string& name);
    void               id(const std::string& id);
    void               statistic(const std::string& statistic);
    const std::string& name(void) const;
    const std::string& id(void) const;
    const std::string& statistic(void) const;
};


/***********************************************************************//**
 * @brief GObservation class extension
 ***************************************************************************/
%extend GObservation {
};
