/***************************************************************************
 *        GCTAResponseIrf.i - CTA instrument response function class       *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2018 by Juergen Knoedlseder                         *
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
 * @file GCTAResponseIrf.i
 * @brief CTA instrument response function class definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GCTAResponseIrf.hpp"
%}


/***********************************************************************//**
 * @class GCTAResponseIrf
 *
 * @brief CTA instrument response function class
 ***************************************************************************/
class GCTAResponseIrf : public GCTAResponse {
public:
    // Constructors and destructors
    GCTAResponseIrf(void);
    GCTAResponseIrf(const GCTAResponseIrf& rsp);
    explicit GCTAResponseIrf(const GXmlElement& xml);
    GCTAResponseIrf(const std::string& rspname, const GCaldb& caldb);
    virtual ~GCTAResponseIrf(void);

    // Implement pure virtual base class methods
    virtual void             clear(void);
    virtual GCTAResponseIrf* clone(void) const;
    virtual std::string      classname(void) const;
    virtual bool             use_edisp(void) const;
    virtual bool             use_tdisp(void) const;
    virtual double           irf(const GEvent&       event,
                                 const GPhoton&      photon,
                                 const GObservation& obs) const;
    virtual double           nroi(const GModelSky&    model,
                                  const GEnergy&      obsEng,
                                  const GTime&        obsTime,
                                  const GObservation& obs) const;
    virtual GEbounds         ebounds(const GEnergy& obsEnergy) const;
    virtual void             read(const GXmlElement& xml);
    virtual void             write(GXmlElement& xml) const;

    // Overload virtual base class methods

    // Other Methods
    bool                  apply_edisp(void) const;
    void                  apply_edisp(const bool& apply_edisp) const;
    GCTAEventAtom*        mc(const double& area, const GPhoton& photon,
                             const GObservation& obs, GRan& ran) const;
    void                  caldb(const GCaldb& caldb);
    const GCaldb&         caldb(void) const;
    const std::string&    rspname(void) const;
    void                  load(const std::string& rspname);
    void                  load_aeff(const GFilename& filename);
    void                  load_psf(const GFilename& filename);
    void                  load_edisp(const GFilename& filename);
    void                  load_background(const GFilename& filename);
    void                  offset_sigma(const double& sigma);
    double                offset_sigma(void) const;
    const GCTAAeff*       aeff(void) const;
    void                  aeff(GCTAAeff* aeff);
    const GCTAPsf*        psf(void) const;
    void                  psf(GCTAPsf* psf);
    const GCTAEdisp*      edisp(void) const;
    void                  edisp(GCTAEdisp* edisp);
    const GCTABackground* background(void) const;
    void                  background(GCTABackground* background);
    const double&         lo_save_thres(void) const;
    const double&         hi_save_thres(void) const;

    // Low-level response methods
    double aeff(const double& theta,
                const double& phi,
                const double& zenith,
                const double& azimuth,
                const double& srcLogEng) const;
    double psf(const double& delta,
               const double& theta,
               const double& phi,
               const double& zenith,
               const double& azimuth,
               const double& srcLogEng) const;
    double psf_delta_max(const double& theta,
                         const double& phi,
                         const double& zenith,
                         const double& azimuth,
                         const double& srcLogEng) const;
    double edisp(const GEnergy& obsEng,
                 const double&  theta,
                 const double&  phi,
                 const double&  zenith,
                 const double&  azimuth,
                 const double&  srcLogEng) const;
    double npsf(const GSkyDir&      srcDir,
                const double&       srcLogEng,
                const GTime&        srcTime,
                const GCTAPointing& pnt,
                const GCTARoi&      roi) const;
};


/***********************************************************************//**
 * @brief GCTAResponseIrf class extension
 ***************************************************************************/
%extend GCTAResponseIrf {
    GCTAResponseIrf copy() {
        return (*self);
    }
%pythoncode {
    def __getstate__(self):
        state = (self.rspname(), self.caldb())
        return state
    def __setstate__(self, state):
        if state[0]:
            self.__init__(state[0], state[1])
        else:
            self.__init__()
}
};
