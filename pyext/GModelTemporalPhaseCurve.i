/***************************************************************************
 *      GModelTemporalPhaseCurve.i - Temporal phase curve model class      *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2017 by Juergen Knoedlseder                              *
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
 * @file GModelTemporalPhaseCurve.i
 * @brief Temporal phase curve model class interface definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GModelTemporalPhaseCurve.hpp"
%}


/***********************************************************************//**
 * @class GModelTemporalPhaseCurve
 *
 * @brief Temporal phase curve model class
 ***************************************************************************/
class GModelTemporalPhaseCurve : public GModelTemporal {

public:
    // Constructors and destructors
    GModelTemporalPhaseCurve(void);
    explicit GModelTemporalPhaseCurve(const GXmlElement& xml);
    GModelTemporalPhaseCurve(const GFilename& filename,
                             const GTime&     mjd,
                             const double&    phase,
                             const double&    f0,
                             const double&    f1,
                             const double&    f2,
                             const double&    norm = 1.0,
                             const bool&      normalize = true);
    GModelTemporalPhaseCurve(const GModelTemporalPhaseCurve& model);
    virtual ~GModelTemporalPhaseCurve(void);

    // Implemented virtual base class methods
    virtual void                      clear(void);
    virtual GModelTemporalPhaseCurve* clone(void) const;
    virtual std::string               classname(void) const;
    virtual std::string               type(void) const;
    virtual double                    eval(const GTime& srcTime,
                                           const bool& gradients = false) const;
    virtual GTimes                    mc(const double& rate, const GTime& tmin,
                                         const GTime& tmax, GRan& ran) const;
    virtual void                      read(const GXmlElement& xml);
    virtual void                      write(GXmlElement& xml) const;

    // Other methods
    const GFilename& filename(void) const;
    void             filename(const GFilename& filename);
    double           norm(void) const;
    void             norm(const double& norm);
    GTime            mjd(void) const;
    void             mjd(const GTime& time);
    double           phase(void) const;
    double           phase(const GTime& time) const;
    void             phase(const double& phase);
    double           f0(void) const;
    void             f0(const double& f0);
    double           f1(void) const;
    void             f1(const double& f1);
    double           f2(void) const;
    void             f2(const double& f2);
    bool             normalize(void) const;
};


/***********************************************************************//**
 * @brief GModelTemporalPhaseCurve class extension
 ***************************************************************************/
%extend GModelTemporalPhaseCurve {
    GModelTemporalPhaseCurve copy() {
        return (*self);
    }
};
