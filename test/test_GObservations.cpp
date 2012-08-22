/***************************************************************************
 *             test_GObservations.cpp  -  Test GObersations class          *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2012 by Jean-Baptiste Cayrou                             *
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
 * @file test_GObservations.cpp
 * @brief Test GObservations class with a test intrument ("testint/")
 * @author J.-B. Cayrou
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <iostream>
#include "GammaLib.hpp"
#include "testinst/GTestLib.hpp"

/* __ OpenMP section _____________________________________________________ */
#ifdef _OPENMP
#include <omp.h>
#ifdef __APPLE__
#ifdef __MACH__
#include <pthread.h>
pthread_attr_t gomp_thread_attr;
#endif
#endif
#endif

/* __ Coding definitions _________________________________________________ */
#define RATE      13.0        //!< Events per seconde. For events generation.
#define UN_BINNED 0
#define BINNED    1


#ifdef _OPENMP
/***********************************************************************//**
 * @brief Test observations optimizer.
 *
 * @param[in] mode Testing mode.
 * 
 * This method supports two testing modes: 0 = unbinned and 1 = binned.
 ***************************************************************************/
GModelPar& test_observations_optimizer(int mode=0)
{
    // Create Test Model
    GTestModelData model;
    
    // Create Models conteners
    GModels models;
    models.append(model);
    
    // Time iterval
    GTime tmin(0,0,   "sec");
    GTime tmax(1800,0,"sec");

    // Rate : events/sec
    double rate = RATE;
    
    // Create observations
    GObservations obs;
    
    // Add some observation
    for (int i=0; i<6; ++i) {

        // Random Generator
        GRan ran;
        ran.seed(i);
        
        GEvents *events;
        
        if (mode == UN_BINNED) {
            // Create a list of events
            events = model.generateList(rate,tmin,tmax,ran);
        }
        else {
            // Create a cube of events
            events = model.generateCube(rate,tmin,tmax,ran);
        }

        // Create an observation
        GTestObservation ob;
        
        // Add events to the observation
        ob.events(events);
        ob.ontime(tmax.met()-tmin.met());
        obs.append(ob);
    }

    // Add the model to the observation
    obs.models(models);
    
    // Create a GLog for show the interations of optimizer.
    GLog log;
  
    // Show in shell
    //log.cout(true);
    
    // Create an optimizer.
    GOptimizerLM opt(log);
    
    opt.max_stalls(50);
    
    // Optimize
    obs.optimize(opt);
    
    //std::cout << obs << std::endl;
    std::cout << opt << std::endl;
   
    std::cout << (obs.models())[0] << std::endl;
    
    if (opt.status() != 0) { //check if converged
        std::cout<<"ERROR : optimizer did not converge."<<std::endl;
        throw;
    }
    
    GModelPar& result = ((obs.models())[0])[0];
    
    if (fabs(result.value()-RATE) > result.error()*3) {
        std::cout<<"ERROR : Value is not precise enough."<<std::endl;
        throw;
    }
    
    return (((obs.models())[0])[0]);
}


/***********************************************************************//**
 * @brief Test optimizer with unbinned events.
 ***************************************************************************/
void test_observations_optimizer_unbinned()
{
    std::cout << "**** Unbinned Test ****" << std::endl <<std::endl;
    
    // Test with one thread
    double t_start = omp_get_wtime();
    
    std::cout << "* Unbinned : Test with 1 thread" << std::endl;
    omp_set_num_threads(1);
    GModelPar& result1 = test_observations_optimizer(UN_BINNED);
    
    double t_elapsed1 = omp_get_wtime()-t_start;
    
    // Test with 10 threads
    
    t_start = omp_get_wtime();
    
    std::cout << "* Unbinned : Test with 10 threads" << std::endl;
    omp_set_num_threads(10);
    GModelPar& result2 = test_observations_optimizer(UN_BINNED);
    
    double t_elapsed2 = omp_get_wtime()-t_start;

    //Compare times.
    std::cout << "Time with 1 thread : " << t_elapsed1 << " s" << std::endl;
    std::cout << "Time with 10 thread : " << t_elapsed2 << " s" << std::endl;
    
    std::cout<<std::endl;
}


/***********************************************************************//**
 * @brief Test optimizer with binned events.
 ***************************************************************************/
void test_observations_optimizer_binned()
{    
    std::cout << "**** Binned Test ****" << std::endl<<std::endl;
    
    // Test with one thread
    double t_start = omp_get_wtime();
    
    std::cout<<"* Binned : Test with 1 thread"<<std::endl;
    omp_set_num_threads(1);
    GModelPar& result1 = test_observations_optimizer(BINNED);
    
    double t_elapsed1 = omp_get_wtime()-t_start;
    
    // Test with 10 threads
    
    t_start = omp_get_wtime();
    std::cout<<"* Binned : Test with 10 threads"<<std::endl;
    omp_set_num_threads(10);
    GModelPar& result2 = test_observations_optimizer(BINNED);
    
    double t_elapsed2 = omp_get_wtime()-t_start;

    //Compare times.
    std::cout<<"Time with 1 thread : "<<t_elapsed1<<" s"<<std::endl;
    std::cout<<"Time with 10 threads : "<<t_elapsed2<<" s"<<std::endl;
}
#endif


/***********************************************************************//**
 * @brief Main test code.
 ***************************************************************************/
int main(void)
{
    // Dump header
    std::cout << std::endl;
    std::cout << "*******************************" << std::endl;
    std::cout << "* GObservations class testing *" << std::endl; 
    std::cout << "*******************************" << std::endl;
    
    std::cout<<"Test openMP results:"<<std::endl;

    //TODO: test gaussian.
    #ifdef _OPENMP
    test_observations_optimizer_unbinned();
    test_observations_optimizer_binned();
    #else
    std::cout << "GammaLib is not compiled with openmp option." << std::endl;
    #endif
    
    // Return
    return 0;
}
