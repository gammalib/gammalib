/***************************************************************************
 *         GSPIModelDataSpace.hpp - INTEGRAL/SPI data space model          *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2020 by Juergen Knoedlseder                              *
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
 * @file GSPIModelDataSpace.hpp
 * @brief INTEGRAL/SPI data space model interface definition
 * @author Juergen Knoedlseder
 */

#ifndef GSPIMODELDATASPACE_HPP
#define GSPIMODELDATASPACE_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include <vector>
#include "GModelData.hpp"
#include "GModelPar.hpp"
#include "GSPIEventCube.hpp"

/* __ Forward declarations _______________________________________________ */
class GEnergy;
class GTime;
class GSPIInstDir;
class GEvent;
class GObservation;
class GXmlElement;
class GSPIObservation;

/* __ Constants __________________________________________________________ */


/***********************************************************************//**
 * @class GSPIModelDataSpace
 *
 * @brief INTEGRAL/SPI data space model
 *
 * This class implements an INTEGRAL/SPI data space model that is based on
 * model values that are found in the INTEGRAL/SPI event cube.
 ***************************************************************************/
class GSPIModelDataSpace : public GModelData {

public:
    // Constructors and destructors
    GSPIModelDataSpace(void);
    GSPIModelDataSpace(const GSPIObservation& obs,
                       const std::string&     name,
                       const std::string&     method,
                       const int&             index);
    explicit GSPIModelDataSpace(const GXmlElement& xml);
    GSPIModelDataSpace(const GSPIModelDataSpace& model);
    virtual ~GSPIModelDataSpace(void);

    // Operators
    virtual GSPIModelDataSpace& operator=(const GSPIModelDataSpace& model);

    // Implemented pure virtual methods
    virtual void                clear(void);
    virtual GSPIModelDataSpace* clone(void) const;
    virtual std::string         classname(void) const;
    virtual std::string         type(void) const;
    virtual bool                is_constant(void) const;
    virtual double              eval(const GEvent&       event,
                                     const GObservation& obs,
                                     const bool&         gradients = false) const;
    virtual double              npred(const GEnergy&      obsEng,
                                      const GTime&        obsTime,
                                      const GObservation& obs) const;
    virtual GSPIEventCube*      mc(const GObservation& obs, GRan& ran) const;
    virtual void                read(const GXmlElement& xml);
    virtual void                write(GXmlElement& xml) const;
    virtual std::string         print(const GChatter& chatter = NORMAL) const;

protected:
    // Protected methods
    void   init_members(void);
    void   copy_members(const GSPIModelDataSpace& model);
    void   free_members(void);
    void   set_pointers(void);
    void   setup_model(const GObservation& obs) const;
    void   setup_pars(GSPIEventCube* cube);
    void   setup_pointing_indices(GSPIEventCube*            cube,
                                  std::vector<int>*         indices,
                                  std::vector<std::string>* names);
    void   setup_detector_indices(GSPIEventCube* cube,
                                  std::vector<int>*         indices,
                                  std::vector<std::string>* names);
    void   setup_energy_indices(GSPIEventCube* cube,
                                std::vector<int>*         indices,
                                std::vector<std::string>* names);
    void   setup_point(GSPIEventCube*            cube,
                       std::vector<int>*         indices,
                       std::vector<std::string>* names);
    void   setup_orbit(GSPIEventCube*            cube,
                       std::vector<int>*         indices,
                       std::vector<std::string>* names);
    void   setup_date(GSPIEventCube*            cube,
                      std::vector<int>*         indices,
                      std::vector<std::string>* names,
                      const double&             time);
    void   add_gedfail(GSPIEventCube*            cube,
                       std::vector<int>*         indices,
                       std::vector<std::string>* names);
    void   add_gedanneal(GSPIEventCube*            cube,
                         std::vector<int>*         indices,
                         std::vector<std::string>* names);
    void   setup_dete(GSPIEventCube*            cube,
                      std::vector<int>*         indices,
                      std::vector<std::string>* names);
    void   setup_evtclass(GSPIEventCube*            cube,
                          std::vector<int>*         indices,
                          std::vector<std::string>* names);
    void   setup_ebin(GSPIEventCube*            cube,
                      std::vector<int>*         indices,
                      std::vector<std::string>* names);
    double get_date_time(const std::string& method) const;
    void   split_pointing_indices(GSPIEventCube*            cube,
                                  std::vector<int>*         indices,
                                  std::vector<std::string>* names,
                                  const GTime&              time,
                                  const std::string&        reason) const;

    // Protected data members
    mutable GSPIObservation* m_obs;         //!< SPI observation
    std::string              m_method;      //!< Fitting method
    int                      m_index;       //!< Index of model in event bins
    int                      m_map_size;    //!< Size of parameter map
    int*                     m_map;         //!< Parameter map
    std::vector<GModelPar>   m_parameters;  //!< Model parameters
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GSPIModelDataSpace").
 ***************************************************************************/
inline
std::string GSPIModelDataSpace::classname(void) const
{
    return ("GSPIModelDataSpace");
}


/***********************************************************************//**
 * @brief Return model type
 *
 * @return Model type.
 *
 * Returns "DataSpace" as model type.
 ***************************************************************************/
inline
std::string GSPIModelDataSpace::type(void) const
{
    return ("DataSpace");
}


/***********************************************************************//**
 * @brief Signal if model is temporally constant
 *
 * @return True.
 *
 * Signals that the data space model is temporally constant.
 ***************************************************************************/
inline
bool GSPIModelDataSpace::is_constant(void) const
{
    return (true);
}

#endif /* GSPIMODELDATASPACE_HPP */
