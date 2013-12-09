/***************************************************************************
 *           GCTAOnOffObservation.i - CTA on-off observation class         *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2013 by Michael Mayer                                    *
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


/***********************************************************************//**
 * @class GCTAOnOffObservation
 *
 * @brief CTA on-off observation class
 ***************************************************************************/
class GCTAOnOffObservation : public GBase {
public:
    // Constructors and destructors
    GCTAOnOffObservation(void);
    GCTAOnOffObservation(const GEbounds& ereco, const GSkyRegions& on,
                         const GSkyRegions& off);
    GCTAOnOffObservation(const GCTAOnOffObservation& obs);
    virtual ~GCTAOnOffObservation(void);
 
    // Methods
    void                  clear(void);
    GCTAOnOffObservation* clone(void) const;
    void                  name(const std::string& name);
    void                  instrument(const std::string& instrument);
    void                  id(const std::string& id);
    void                  on_regions(const GSkyRegions& regions);
    void                  off_regions(const GSkyRegions& regions);
    const std::string&    name(void) const;
    const std::string&    instrument(void) const;
    const std::string&    id(void) const;
    const GPha&           on_spec(void) const;
    const GPha&           off_spec(void) const;
    const GArf&           arf(void) const;
    const GRmf&           rmf(void) const;
    void                  fill(const GCTAObservation& obs);
    void                  compute_response(const GCTAObservation& obs,
                                           const GEbounds& etrue);
    void                  read(const GXmlElement& xml);
    void                  write(GXmlElement& xml) const;
};


/***********************************************************************//**
 * @brief GCTOnOffAObservation class extension
 ***************************************************************************/
%extend GCTAOnOffObservation {
    GCTAOnOffObservation copy() {
        return (*self);
    }
}
