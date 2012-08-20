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
#ifdef __APPLE__ & __MACH__
#include <pthread.h>
pthread_attr_t gomp_thread_attr;
#endif
#endif

/* __ Coding definitions _________________________________________________ */
#define RATE      13.0        //!< Events per seconde. For events generation.
#define UN_BINNED 0
#define BINNED    1


#ifdef _OPENMP

class TestGObservation : public GTestSuite
{
    public:
        TestGObservation(const std::string& name) : GTestSuite(name){ return; }

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

            // Create an optimizer.
            GOptimizerLM opt(log);

            opt.max_stalls(50);

            // Optimize
            obs.optimize(opt);

            //std::cout << obs << std::endl;
            //std::cout << opt << std::endl;

            // Get the result
            GModelPar& result = ((obs.models())[0])[0];

            //check if converged
            test_assert(opt.status()==0,"Check if converged","Optimizer did not convered"); 

            //check if value is correct
            test_assert(fabs(result.value()-RATE) < result.error()*3,"check if value is correct","Value is not precise enough."); 

            return (((obs.models())[0])[0]);
        }

        /***********************************************************************//**
        * @brief Test optimizer with unbinned events and 1 thread
        ***************************************************************************/
        void test_observations_optimizer_unbinned_1()
        {
             // Test with 1 thread
            omp_set_num_threads(1);
            test_observations_optimizer(UN_BINNED);

            return;
        }

        /***********************************************************************//**
         * @brief Test optimizer with unbinned events and 10 thread
         ***************************************************************************/
        void test_observations_optimizer_unbinned_10()
        {
             // Test with 10 thread
            omp_set_num_threads(10);
            test_observations_optimizer(UN_BINNED);

            return;
        }

        /***********************************************************************//**
         * @brief Test optimizer with binned events and 1 thread
         ***************************************************************************/
        void test_observations_optimizer_binned_1()
        {
             // Test with 1 thread
            omp_set_num_threads(1);
            test_observations_optimizer(BINNED);

            return;
        }

        /***********************************************************************//**
         * @brief Test optimizer with binned events and 10 thread
         ***************************************************************************/
        void test_observations_optimizer_binned_10()
        {
             // Test with 10 thread
            omp_set_num_threads(10);
            test_observations_optimizer(BINNED);

            return;
        }

        /***********************************************************************//**
         * @brief Set tests
        ***************************************************************************/
        void set(void)
        {
            //Unbinned
            add_test(static_cast<pfunction>(&TestGObservation::test_observations_optimizer_unbinned_1),"Unbinned 1 thread");
            add_test(static_cast<pfunction>(&TestGObservation::test_observations_optimizer_unbinned_10),"Unbinned 10 threads");

            //Binned
            add_test(static_cast<pfunction>(&TestGObservation::test_observations_optimizer_binned_1),"Binned 1 thread");
            add_test(static_cast<pfunction>(&TestGObservation::test_observations_optimizer_binned_10),"Binned 10 thread");

            return;
        }
};

#endif


/***********************************************************************//**
 * @brief Main test code.
 ***************************************************************************/
int main(void)
{
    GTestSuites testsuites("GObservations class testing");

    //TODO: test gaussian.
    #ifdef _OPENMP
    //Create a test suite
    TestGObservation test("Test OpenMP");
    //Append to the container
    testsuites.append(test);

    //Run
    bool was_successful=testsuites.run();

    #else
            std::cout<<"GammaLib is not compiled with openmp option."<<std::endl;
    #endif

    //save xml report
    testsuites.save("reports/GObservations.xml");

    // Return
    return was_successful;
}
