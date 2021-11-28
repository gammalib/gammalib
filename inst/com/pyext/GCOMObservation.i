/***************************************************************************
 *             GCOMObservation.i - COMPTEL observation class               *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2012-2021 by Juergen Knoedlseder                         *
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
 * @file GCOMObservation.i
 * @brief COMPTEL observation class interface definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GCOMObservation.hpp"
%}


/***********************************************************************//**
 * @class GCOMObservation
 *
 * @brief Interface class for COMPTEL observations
 ***************************************************************************/
class GCOMObservation : public GObservation {
public:
    // Constructors and destructors
    GCOMObservation(void);
    explicit GCOMObservation(const GXmlElement& xml);
    GCOMObservation(const GCOMDri& dre,
                    const GCOMDri& drb,
                    const GCOMDri& drg,
                    const GCOMDri& drx);
    GCOMObservation(const GFilename& drename,
                    const GFilename& drbname,
                    const GFilename& drgname,
                    const GFilename& drxname);
    GCOMObservation(const GFilename&              evpname,
                    const GFilename&              timname,
                    const std::vector<GFilename>& oadnames);
    GCOMObservation(const GCOMObservation& obs);
    virtual ~GCOMObservation(void);

    // Implement pure virtual methods
    virtual void                clear(void);
    virtual GCOMObservation*    clone(void) const;
    virtual std::string         classname(void) const;
    virtual void                response(const GResponse& rsp);
    virtual const GCOMResponse* response(void) const;
    virtual std::string         instrument(void) const;
    virtual double              ontime(void) const;
    virtual double              livetime(void) const;
    virtual double              deadc(const GTime& time = GTime()) const;
    virtual void                read(const GXmlElement& xml);
    virtual void                write(GXmlElement& xml) const;

    // Other methods
    bool             is_unbinned(void) const;
    bool             is_binned(void) const;
    void             load(const GFilename& drename,
                          const GFilename& drbname,
                          const GFilename& drgname,
                          const GFilename& drxname);
    void             load(const GFilename&              evpname,
                          const GFilename&              timname,
                          const std::vector<GFilename>& oadnames);
    void             response(const GCaldb& caldb, const std::string& rspname);
    void             response(const GCOMResponse& response);
    void             obs_id(const double& id);
    void             ontime(const double& ontime);
    void             livetime(const double& livetime);
    void             deadc(const double& deadc);
    void             ewidth(const double& ewidth);
    const double&    obs_id(void) const;
    const double&    ewidth(void) const;
    const GCOMDri&   drb(void) const;
    const GCOMDri&   drg(void) const;
    const GCOMDri&   drx(void) const;
    GCOMDri          drm(const GModels& models) const;
    const GCOMTim&   tim(void) const;
    const GCOMOads&  oads(void) const;
    const GFilename& drename(void) const;
    const GFilename& drbname(void) const;
    const GFilename& drgname(void) const;
    const GFilename& drxname(void) const;
    const int&       phi_first(void) const;
    const int&       phi_last(void) const;
    void             drename(const GFilename& drename);
    void             drbname(const GFilename& drbname);
    void             drgname(const GFilename& drgname);
    void             drxname(const GFilename& drxname);
    void             phi_first(const int& phi_first);
    void             phi_last(const int& phi_last);
    void             compute_drb(const std::string& method,
                                 const GCOMDri&     drm,
                                 const int&         nrunav = 3,
                                 const int&         navgr  = 3,
                                 const int&         nincl  = 13,
                                 const int&         nexcl  = 0);
};


/***********************************************************************//**
 * @brief GCOMObservation class extension
 ***************************************************************************/
%extend GCOMObservation {
    GCOMObservation copy() {
        return (*self);
    }
%pythoncode {
    def __getstate__(self):
        xml = gammalib.GXmlElement()
        self.write(xml)
        state = (xml,)
        return state
    def __setstate__(self, state):
        if state[0].elements() > 0:
            self.__init__(state[0])
        else:
            self.__init__()
}
};
