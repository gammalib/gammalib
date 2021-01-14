/***************************************************************************
 *          GCTAObservation.i - CTA Observation class interface            *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2021 by Juergen Knoedlseder                         *
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
 * @file GCTAObservation.i
 * @brief CTA observation interface definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GCTAObservation.hpp"
%}


/***********************************************************************//**
 * @class GCTAObservation
 *
 * @brief CTA observation class
 ***************************************************************************/
class GCTAObservation : public GObservation {
public:
    // Constructors and destructors
    GCTAObservation(void);
    explicit GCTAObservation(const GFilename& filename);
    GCTAObservation(const GFilename& cntcube,
                    const GFilename& expcube,
                    const GFilename& psfcube,
                    const GFilename& bkgcube);
    GCTAObservation(const GFilename& cntcube,
                    const GFilename& expcube,
                    const GFilename& psfcube,
                    const GFilename& edispcube,
                    const GFilename& bkgcube);
    GCTAObservation(const GCTAObservation& obs);
    virtual ~GCTAObservation(void);

    // Implemented pure virtual base class methods
    virtual void                clear(void);
    virtual GCTAObservation*    clone(void) const;
    virtual std::string         classname(void) const;
    virtual void                response(const GResponse& rsp);
    virtual const GCTAResponse* response(void) const;
    virtual std::string         instrument(void) const;
    virtual double              ontime(void) const;
    virtual double              livetime(void) const;
    virtual double              deadc(const GTime& time = GTime()) const;
    virtual void                read(const GXmlElement& xml);
    virtual void                write(GXmlElement& xml) const;

    // Other methods
    void                instrument(const std::string& instrument);
    bool                has_response(void) const;
    bool                has_events(void) const;
    void                read(const GFits& fits);
    void                write(GFits& fits,
                              const std::string& evtname = "EVENTS",
                              const std::string& gtiname = "GTI") const;
    void                load(const GFilename& filename);
    void                load(const GFilename& cntcube,
                             const GFilename& expcube,
                             const GFilename& psfcube,
                             const GFilename& bkgcube);
    void                load(const GFilename& cntcube,
                             const GFilename& expcube,
                             const GFilename& psfcube,
                             const GFilename& edispcube,
                             const GFilename& bkgcube);
    void                save(const GFilename& filename,
                             const bool&      clobber = false) const;
    void                response(const std::string& rspname,
                                 const GCaldb&      caldb);
    void                response(const GCTACubeExposure& expcube,
                                 const GCTACubePsf&      psfcube,
                                 const GCTACubeBackground& bkgcube);
    void                response(const GCTACubeExposure&   expcube,
                                 const GCTACubePsf&        psfcube,
                                 const GCTACubeEdisp&      edispcube,
                                 const GCTACubeBackground& bkgcube);
    void                pointing(const GCTAPointing& pointing);
    const GCTAPointing& pointing(void) const;
    void                off_regions(const GSkyRegions& off_regions);
    const GSkyRegions&  off_regions(void) const;
    GCTARoi             roi(void) const;
    GGti                gti(void) const;
    GEbounds            ebounds(void) const;
    void                object(const std::string& object);
    const std::string&  object(void) const;
    void                ra_obj(const double& ra);
    const double&       ra_obj(void) const;
    void                dec_obj(const double& dec);
    const double&       dec_obj(void) const;
    void                ontime(const double& ontime);
    void                livetime(const double& livetime);
    void                deadc(const double& deadc);
    void                eventfile(const GFilename& filename);
    const GFilename&    eventfile(void) const;
    std::string         eventtype(void) const;
    void                dispose_events(void);
    void                lo_user_thres(const double& lo_user_thres);
    const double&       lo_user_thres(void) const;
    void                hi_user_thres(const double& hi_user_thres);
    const double&       hi_user_thres(void) const;
    void                n_tels(const int& tels);
    const int&          n_tels(void) const;
};


/***********************************************************************//**
 * @brief GCTAObservation class extension
 ***************************************************************************/
%extend GCTAObservation {
    GCTAObservation copy() {
        return (*self);
    }
%pythoncode {
    def __getstate__(self):
        if self.has_response():
            response = self.response()
        else:
            response = None
        state = (gammalib.GObservation.__getstate__(self),
                 self.instrument(), self.eventfile(), response,
                 self.pointing(), self.off_regions(), self.ontime(),
                 self.livetime(), self.deadc(), self.ra_obj(), self.dec_obj(),
                 self.lo_user_thres(), self.hi_user_thres(), self.n_tels())
        return state
    def __setstate__(self, state):
        self.__init__()
        gammalib.GObservation.__setstate__(self, state[0])
        self.instrument(state[1])
        self.eventfile(state[2])
        if state[3] != None:
            self.response(state[3])
        self.pointing(state[4])
        self.off_regions(state[5])
        self.ontime(state[6])
        self.livetime(state[7])
        self.deadc(state[8])
        self.ra_obj(state[9])
        self.dec_obj(state[10])
        self.lo_user_thres(state[11])
        self.hi_user_thres(state[12])
        self.n_tels(state[13])
}
};
