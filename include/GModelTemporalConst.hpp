/***************************************************************************
 *        GModelTemporalConst.hpp  -  Temporal constant model class        *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2009-2010 by Jurgen Knodlseder                           *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GModelTemporalConst.hpp
 * @brief GModelTemporalConst class interface definition.
 * @author J. Knodlseder
 */

#ifndef GMODELTEMPORALCONST_HPP
#define GMODELTEMPORALCONST_HPP

/* __ Includes ___________________________________________________________ */
#include <iostream>
#include "GLog.hpp"
#include "GModelPar.hpp"
#include "GModelTemporal.hpp"
#include "GTime.hpp"


/***********************************************************************//**
 * @class GModelTemporalConst
 *
 * @brief Interface definition for the constant model class.
 *
 * This class implements the temporal component of the factorised model for
 * a model that is constant in time.
 ***************************************************************************/
class GModelTemporalConst  : public GModelTemporal {

    // I/O friends
    friend std::ostream& operator<< (std::ostream& os, const GModelTemporalConst& model);
    friend GLog&         operator<< (GLog& log, const GModelTemporalConst& model);

public:
    // Constructors and destructors
    GModelTemporalConst(void);
    GModelTemporalConst(const GModelTemporalConst& model);
    virtual ~GModelTemporalConst(void);

    // Operators
    GModelTemporalConst& operator= (const GModelTemporalConst& model);

    // Implemented virtual methods
    void                 clear(void);
    GModelTemporalConst* clone(void) const;
    int                  size(void) const { return m_npars; }
    std::string          type(void) const { return "Constant"; }
    double               eval(const GTime& srcTime);
    double               eval_gradients(const GTime& srcTime);
    void                 read(const GXmlElement& xml);
    void                 write(GXmlElement& xml) const;
    std::string          print(void) const;

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GModelTemporalConst& model);
    void free_members(void);

    // Implemented pure virtual methods
    GModelPar** par(void) { return m_par; }

    // Protected members
    int        m_npars;           //!< Number of parameters
    GModelPar* m_par[1];          //!< Pointers to parameters
    GModelPar  m_norm;            //!< Constant
};

#endif /* GMODELTEMPORALCONST_HPP */
