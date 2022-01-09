/***************************************************************************
 *                  GCOMDri.hpp - COMPTEL Data Space class                 *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2017-2022 by Juergen Knoedlseder                         *
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
 * @file GCOMDri.hpp
 * @brief COMPTEL Data Space class definition
 * @author Juergen Knoedlseder
 */

#ifndef GCOMDRI_HPP
#define GCOMDRI_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GBase.hpp"
#include "GSkyMap.hpp"
#include "GEbounds.hpp"
#include "GGti.hpp"
#include "GCOMSelection.hpp"

/* __ Forward declarations _______________________________________________ */
class GFilename;
class GFits;
class GFitsHDU;
class GFitsImage;
class GModel;
class GModelSky;
class GCOMOad;
class GCOMTim;
class GCOMStatus;
class GCOMObservation;

/* __ Constants __________________________________________________________ */
namespace gammalib {
    const std::string extname_dri = "DRI";
}


/***********************************************************************//**
 * @class GCOMDri
 *
 * @brief COMPTEL Data Space class
 ***************************************************************************/
class GCOMDri : public GBase {

public:
    // Constructors and destructors
    GCOMDri(void);
    explicit GCOMDri(const GFilename& filename);
    GCOMDri(const GSkyMap& map, const double& phimin = 0.0,
                                const double& phibin = 0.0,
                                const int&    nphibin = 0);
    GCOMDri(const GCOMDri& dri);
    virtual ~GCOMDri(void);

    // Operators
    GCOMDri&      operator=(const GCOMDri& dri);
    double&       operator[](const int& index);
    const double& operator[](const int& index) const;

    // Implemented pure virtual base class methods
    virtual void        clear(void);
    virtual GCOMDri*    clone(void) const;
    virtual std::string classname(void) const;
    virtual std::string print(const GChatter& chatter = NORMAL) const;

    // Other methods
    int                size(void) const;
    int                nchi(void) const;
    int                npsi(void) const;
    int                nphibar(void) const;
    const GSkyMap&     map(void) const;
    const std::string& name(void) const;
    void               name(const std::string& name);
    const GEbounds&    ebounds(void) const;
    void               ebounds(const GEbounds& ebounds);
    const GGti&        gti(void) const;
    void               gti(const GGti& gti);
    const double&      phimin(void) const;
    const double&      phibin(void) const;
    const double&      tof_correction(void) const;
    void               tof_correction(const double& tofcor);
    const double&      phase_correction(void) const;
    void               phase_correction(const double& phasecor);
    const int&         num_superpackets(void) const;
    void               num_superpackets(const int& number);
    const int&         num_used_superpackets(void) const;
    void               num_used_superpackets(const int& number);
    const int&         num_skipped_superpackets(void) const;
    void               num_skipped_superpackets(const int& number);
    void               compute_dre(const GCOMObservation& obs,
                                   const GCOMSelection&   select = GCOMSelection(),
                                   const double&          zetamin = 5.0);
    void               compute_drg(const GCOMObservation& obs,
                                   const GCOMSelection&   select = GCOMSelection(),
                                   const double&          zetamin = 5.0);
    void               compute_drx(const GCOMObservation& obs,
                                   const GCOMSelection&   select = GCOMSelection());
    void               compute_drm(const GCOMObservation& obs,
                                   const GModel&          model);
    double             cone_content(const GSkyDir& dir,
                                    const double&  armmin,
                                    const double&  armmax) const;
    void               load(const GFilename& filename);
    void               save(const GFilename& filename,
                            const bool&      clobber = false) const;
    void               read(const GFitsImage& image);
    void               write(GFits&             fits,
                             const std::string& extname = "") const;

protected:
    // Protected methods
    void   init_members(void);
    void   copy_members(const GCOMDri& dri);
    void   free_members(void);
    void   init_cube(void);
    void   init_statistics(void);
    bool   use_superpacket(const GCOMOad&       oad,
                           const GCOMTim&       tim,
                           const GCOMSelection& select);
    void   read_attributes(const GFitsHDU* hdu);
    void   write_attributes(GFitsHDU* hdu) const;
    double compute_geometry(const int& tjd, const double&        theta,
                                            const double&        phi,
                                            const GCOMSelection& select,
                                            const GCOMStatus&    status) const;
    double compute_surface(const double& x1, const double& y1, const double& r1,
                           const double& x2, const double& y2, const double& r2) const;
    double compute_overlap(const double& x1, const double& y1, const double& r1,
                           const double& x2, const double& y2, const double& r2,
                           const double& x3, const double& y3, const double& r3) const;
    void   compute_tof_correction(void);

    // Protected members
    std::string m_name;     //!< Data cube name
    GSkyMap     m_dri;      //!< Data cube
    GEbounds    m_ebounds;  //!< Energy boundaries of data cube
    GGti        m_gti;      //!< Good Time Intervals of data cube
    double      m_phimin;   //!< Phibar minimum (deg)
    double      m_phibin;   //!< Phibar binsize (deg)
    double      m_tofcor;   //!< ToF correction
    double      m_phasecor; //!< Pulsar phase correction

    // Computation statistics
    GTime m_tstart;                   //!< Selection start time
    GTime m_tstop;                    //!< Selection stop time
    int   m_num_superpackets;         //!< Number of superpackets
    int   m_num_used_superpackets;    //!< Number of used superpackets
    int   m_num_skipped_superpackets; //!< Number of skipped superpackets

    // Selection parameters
    bool          m_has_selection;    //!< Signal that selection was applied
    GCOMSelection m_selection;        //!< Selection parameters
    double        m_zetamin;          //!< Minimum zeta angle
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GCOMDri").
 ***************************************************************************/
inline
std::string GCOMDri::classname(void) const
{
    return ("GCOMDri");
}


/***********************************************************************//**
 * @brief DRI bin access operators
 *
 * @param[in] index DRI bin index [0,...,size()-1].
 * @return Reference to DRI bin.
 ***************************************************************************/
inline
double& GCOMDri::operator[](const int& index)
{
    return (const_cast<double&>((m_dri.pixels()[index])));
}


/***********************************************************************//**
 * @brief DRI bin access operators (const version)
 *
 * @param[in] index DRI bin index [0,...,size()-1].
 * @return Reference to DRI bin.
 ***************************************************************************/
inline
const double& GCOMDri::operator[](const int& index) const
{
    return (m_dri.pixels()[index]);
}


/***********************************************************************//**
 * @brief Return number of bins
 *
 * @return Number of bins.
 ***************************************************************************/
inline
int GCOMDri::size(void) const
{
    return (m_dri.npix()*m_dri.nmaps());
}


/***********************************************************************//**
 * @brief Return number of Chi bins
 *
 * @return Number of Chi bins.
 ***************************************************************************/
inline
int GCOMDri::nchi(void) const
{
    return (m_dri.nx());
}


/***********************************************************************//**
 * @brief Return number of Psi bins
 *
 * @return Number of Psi bins.
 ***************************************************************************/
inline
int GCOMDri::npsi(void) const
{
    return (m_dri.ny());
}


/***********************************************************************//**
 * @brief Return number of Phibar bins
 *
 * @return Number of Phibar bins.
 ***************************************************************************/
inline
int GCOMDri::nphibar(void) const
{
    return (m_dri.nmaps());
}


/***********************************************************************//**
 * @brief Return DRI sky map
 *
 * @return Sky map containing DRI data.
 ***************************************************************************/
inline
const GSkyMap& GCOMDri::map(void) const
{
    return (m_dri);
}


/***********************************************************************//**
 * @brief Return DRI cube name
 *
 * @return DRI cube name.
 ***************************************************************************/
inline
const std::string& GCOMDri::name(void) const
{
    return (m_name);
}


/***********************************************************************//**
 * @brief Set DRI cube name
 *
 * @param[in] name DRI cube name.
 *
 * Sets the name of the DRI cube.
 ***************************************************************************/
inline
void GCOMDri::name(const std::string& name)
{
    m_name = name;
    return;
}


/***********************************************************************//**
 * @brief Return energy boundaries of DRI cube
 *
 * @return Energy boundaries of DRI cube.
 ***************************************************************************/
inline
const GEbounds& GCOMDri::ebounds(void) const
{
    return (m_ebounds);
}


/***********************************************************************//**
 * @brief Set energy boundaries of DRI cube
 *
 * @param[in] ebounds Energy boundaries of DRI cube.
 *
 * Sets energy boundaries of DRI cube
 ***************************************************************************/
inline
void GCOMDri::ebounds(const GEbounds& ebounds)
{
    m_ebounds = ebounds;
    return;
}


/***********************************************************************//**
 * @brief Return Good Time Intervals of DRI cube
 *
 * @return Good Time Intervals of DRI cube.
 ***************************************************************************/
inline
const GGti& GCOMDri::gti(void) const
{
    return (m_gti);
}


/***********************************************************************//**
 * @brief Set Good Time Intervals of DRI cube
 *
 * @param[in] gti Good Time Intervals of DRI data.
 *
 * Sets the Good Time Intervals of DRI cube.
 ***************************************************************************/
inline
void GCOMDri::gti(const GGti& gti)
{
    m_gti = gti;
    return;
}


/***********************************************************************//**
 * @brief Return minimum Compton scatter angle of DRI cube
 *
 * @return Minimum Compton scatter angle of DRI cube (deg).
 ***************************************************************************/
inline
const double& GCOMDri::phimin(void) const
{
    return (m_phimin);
}


/***********************************************************************//**
 * @brief Return Compton scatter angle bin of DRI cube
 *
 * @return Compton scatter angle bin of DRI cube (deg).
 ***************************************************************************/
inline
const double& GCOMDri::phibin(void) const
{
    return (m_phibin);
}


/***********************************************************************//**
 * @brief Return ToF correction factor
 *
 * @return ToF correction factor.
 *
 * Returns the ToF correction factor that corrects for the event selection
 * in a ToF window.
 ***************************************************************************/
inline
const double& GCOMDri::tof_correction(void) const
{
    return (m_tofcor);
}


/***********************************************************************//**
 * @brief Set ToF correction factor
 *
 * @param[in] tofcor ToF correction factor.
 *
 * Set the ToF correction factor that corrects for the event selection
 * in a ToF window.
 ***************************************************************************/
inline
void GCOMDri::tof_correction(const double& tofcor)
{
    m_tofcor = tofcor;
    return;
}


/***********************************************************************//**
 * @brief Return pulsar phase correction factor
 *
 * @return Pulsar phase correction factor.
 *
 * Returns the pulsar phase correction factor that corrects for the phase
 * selection for pulsar analysis.
 ***************************************************************************/
inline
const double& GCOMDri::phase_correction(void) const
{
    return (m_phasecor);
}


/***********************************************************************//**
 * @brief Set pulsar phase correction factor
 *
 * @param[in] phasecor Pulsar phase correction factor.
 *
 * Set the pulsar phase correction factor that corrects for the phase
 * selection for pulsar analysis.
 ***************************************************************************/
inline
void GCOMDri::phase_correction(const double& phasecor)
{
    m_phasecor = phasecor;
    return;
}


/***********************************************************************//**
 * @brief Return number of superpackets read for DRI
 *
 * @return Number of superpackets read for DRI.
 *
 * Returns the number of superpackets read for DRI.
 ***************************************************************************/
inline
const int& GCOMDri::num_superpackets(void) const
{
    return (m_num_superpackets);
}


/***********************************************************************//**
 * @brief Set number of superpackets read for DRI
 *
 * @param[in] number Number of superpackets read for DRI.
 *
 * Set the number of superpackets read for DRI.
 ***************************************************************************/
inline
void GCOMDri::num_superpackets(const int& number)
{
    m_num_superpackets = number;
    return;
}


/***********************************************************************//**
 * @brief Return number of superpackets used for DRI
 *
 * @return Number of superpackets used for DRI.
 *
 * Returns the number of superpackets used for DRI.
 ***************************************************************************/
inline
const int& GCOMDri::num_used_superpackets(void) const
{
    return (m_num_used_superpackets);
}


/***********************************************************************//**
 * @brief Set number of superpackets used for DRI
 *
 * @param[in] number Number of superpackets used for DRI.
 *
 * Set the number of superpackets used for DRI.
 ***************************************************************************/
inline
void GCOMDri::num_used_superpackets(const int& number)
{
    m_num_used_superpackets = number;
    return;
}


/***********************************************************************//**
 * @brief Return number of superpackets skipped for DRI
 *
 * @return Number of superpackets skipped for DRI.
 *
 * Returns the number of superpackets skipped for DRI.
 ***************************************************************************/
inline
const int& GCOMDri::num_skipped_superpackets(void) const
{
    return (m_num_skipped_superpackets);
}


/***********************************************************************//**
 * @brief Set number of superpackets skipped for DRI
 *
 * @param[in] number Number of superpackets skipped for DRI.
 *
 * Set the number of superpackets skipped for DRI.
 ***************************************************************************/
inline
void GCOMDri::num_skipped_superpackets(const int& number)
{
    m_num_skipped_superpackets = number;
    return;
}

#endif /* GCOMDRI_HPP */
