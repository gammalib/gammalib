/***************************************************************************
 *           GCTAOnOffObservation.i - CTA on-off observation class         *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2013-2017 by Michael Mayer                               *
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
 * @file GCTAOnOffObservation.i
 * @brief CTA on-off observation class definition
 * @author Michael Mayer
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GCTAObservation.hpp"
#include "GCTAOnOffObservation.hpp"
%}


/* __ Typemaps ___________________________________________________________ */
/*
 * The following typemap allows to digest GObservation references for
 * arguments that require GCTAObservation references. This typemap is needed
 * to implement automatic type casting.
 */
%typemap(in) GCTAObservation& (GObservation *argp = NULL, int res = 0) {
  res = SWIG_ConvertPtr($input, (void**)&argp, $descriptor(GCTAObservation *), 0 | 0);
  if (!SWIG_IsOK(res)) {
    res = SWIG_ConvertPtr($input, (void**)&argp, $descriptor(GObservation *), 0 | 0);
    if (!SWIG_IsOK(res)) {
      SWIG_exception_fail(SWIG_ArgError(res), "in method '" "$symname" "', argument "
                         "$argnum"" of type '" "$type""'");
    }
  }
  $1 = dynamic_cast<$ltype>(argp);
  if ($1 == NULL) {
    SWIG_exception_fail(SWIG_ArgError(res), "in method '" "$symname" "', argument "
                        "$argnum"" not of type '" "$type""'");
  }
}


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
 * @class GCTAOnOffObservation
 *
 * @brief CTA On/Off observation class
 ***************************************************************************/
class GCTAOnOffObservation : public GObservation {

public:
    // Constructors and destructors
    GCTAOnOffObservation(void);
    GCTAOnOffObservation(const GObservations& obs);
    GCTAOnOffObservation(const GPha& pha_on,
                         const GPha& pha_off,
                         const GArf& arf,
                         const GRmf& rmf);
    GCTAOnOffObservation(const GCTAObservation& obs,
                         const GSkyDir&         srcdir,
                         const GEbounds&        etrue,
                         const GEbounds&        ereco,
                         const GSkyRegions&     on,
                         const GSkyRegions&     off);
    GCTAOnOffObservation(const GCTAOnOffObservation& obs);
    virtual ~GCTAOnOffObservation(void);
 
    // Implemented pure virtual methods
    virtual void                  clear(void);
    virtual GCTAOnOffObservation* clone(void) const;
	virtual std::string           classname(void) const;
    virtual void                  response(const GResponse& rsp);
    virtual const GResponse*      response(void) const;
    virtual std::string           instrument(void) const;
	virtual double                ontime(void) const;
	virtual double                livetime(void) const;
	virtual double                deadc(const GTime& time = GTime()) const;
    virtual void                  read(const GXmlElement& xml);
    virtual void                  write(GXmlElement& xml) const;

    // Overloaded virtual methods
	virtual double likelihood(const GModels& models,
                              GVector*       gradient,
                              GMatrixSparse* curvature,
                              double*        npred) const;
    virtual int    nobserved(void) const;

    // Other methods
    void        instrument(const std::string& instrument);
    const GPha& on_spec(void) const;
    const GPha& off_spec(void) const;
    const GArf& arf(void) const;
    const GRmf& rmf(void) const;
    GPha model_gamma(const GModels& models) const;
};


/***********************************************************************//**
 * @brief GCTOnOffAObservation class extension
 ***************************************************************************/
%extend GCTAOnOffObservation {
    GCTAOnOffObservation copy() {
        return (*self);
    }
}
