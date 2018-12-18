// INCLUDES
// Cosmosis includes
#include <cosmosis/datablock/datablock.hh>
#include <cosmosis/datablock/section_names.h>
// CppTransport includes
#include "gelaton_mpi.h"
#include <boost/math/interpolators/barycentric_rational.hpp>
#include <boost/range/adaptors.hpp>
#include <boost/math/tools/roots.hpp>
#include "transport-runtime/utilities/spline1d.h"
#include "transport-runtime/tasks/integration_detail/abstract.h"
#include "transport-runtime/models/advisory_classes.h"
#include <memory>
#include "transport-runtime/tasks/integration_detail/twopf_task.h"
#include "transport-runtime/enumerations.h"
#include "boost/filesystem.hpp"
// Include batcher file for integration
#include "sampling_integration_batcher.h"

// Conveniences for CppTransport
using DataType = double;
using StateType = std::vector<DataType>;

namespace inflation {
    const char *sectionName = "cppt_sample";
    const char *paramsSection = "inflation_parameters";
    
    double M_P, V0, eta_R, g_R, lambda_R, alpha, R0, R_init, theta_init, Rdot_init, thetadot_init;

    // function to create a transport::parameters<DataType> object called params
    StateType parameter_creator(double R0, double V0, double eta_R, double g_R, double lambda_R, double alpha) {
        StateType output{R0, V0, eta_R, g_R, lambda_R, alpha};
        return output; 
    }

    // function to create a StateType vector containing intitial conditons
    StateType init_cond_creator(double R_init, double theta_init, double Rdot_init, double thetadot_init) {
        StateType output{R_init, theta_init, Rdot_init, thetadot_init};
        return output;
    }
}

static transport::local_environment env;
static transport::argument_cache arg;
static std::unique_ptr< transport::gelaton_mpi<DataType, StateType> > model;
static std::unique_ptr< transport::twopf_task<DataType> > tk2;

extern "C" {

void * setup(cosmosis::DataBlock * options)
{
    // Read options from the CosmoSIS configuration ini file,
    // passed via the "options" argument
    options->get_val(inflation::sectionName, "M_P", inflation::M_P);

    // Record any configuration information required
    model = std::make_unique< transport::gelaton_mpi<DataType, StateType> > (env, arg);

    // Pass back any object you like
    return options;
}

DATABLOCK_STATUS execute(cosmosis::DataBlock * block, void * config)
{
    // Config is whatever you returned from setup above
    // Block is the collection of parameters and calculations for
    // this set of cosmological parameters
    block->get_val(inflation::paramsSection, "V_0",           inflation::V0);
    block->get_val(inflation::paramsSection, "eta_R",         inflation::eta_R);
    block->get_val(inflation::paramsSection, "g_R",           inflation::g_R);
    block->get_val(inflation::paramsSection, "lambda_R",      inflation::lambda_R);
    block->get_val(inflation::paramsSection, "alpha",         inflation::alpha);
    block->get_val(inflation::paramsSection, "R0",            inflation::R0);
    block->get_val(inflation::paramsSection, "R_init",        inflation::R_init);
    block->get_val(inflation::paramsSection, "theta_init",    inflation::theta_init);
    block->get_val(inflation::paramsSection, "Rdot_init",     inflation::Rdot_init);
    block->get_val(inflation::paramsSection, "thetadot_init", inflation::thetadot_init);

    std::cout << "V0 = " << inflation::V0 << std::endl;
    std::cout << "eta_R = " << inflation::eta_R << std::endl;
    std::cout << "g_R = " << inflation::g_R << std::endl;
    std::cout << "lambda_R = " << inflation::lambda_R << std::endl;
    std::cout << "alpha = " << inflation::alpha << std::endl;
    std::cout << "R0 = " << inflation::R0 << std::endl;
    std::cout << "R_init = " << inflation::R_init << std::endl;
    std::cout << "Rdot_init = " << inflation::Rdot_init << std::endl;
    std::cout << "theta_init = " << inflation::theta_init << std::endl;
    std::cout << "thetadot_init = " << inflation::thetadot_init << std::endl;

    // Set-up initial time for integration (N_init) and N_pre which is used to set the amount of sub-horizon evolution to integrate
    // before the chosen mode crosses the horizon. TODO: maybe move this to the set-up function above to allow it to be set in the
    // cosmosis ini file.
    const double N_init        = 0.0;
    const double N_pre         = 15.0;

    // Create the parameters and initial_conditions objects that CppTransport needs using
    // the two functions defined above.
    using namespace inflation;
    transport::parameters<DataType> params{M_P, inflation::parameter_creator(R0, V0, eta_R, g_R, lambda_R, alpha), model.get()};
    transport::initial_conditions<DataType> ics{"gelaton", params, inflation::init_cond_creator(R_init, theta_init, Rdot_init, thetadot_init), N_init, N_pre};
    
    // Use a silly end value to find nEND and set the time range used by CppT to finish at nEND.
    double Nendhigh = 10000;
    transport::basic_range<DataType> dummy_times{N_init, Nendhigh, 2, transport::spacing::linear};
    transport::background_task<DataType> bkg{ics, dummy_times};

    // From here, we need to enclose the rest of the code in a try-catch statement in order to catch
    // when a particular set of initial conditions fails to integrate or if there is a problem with the
    // physics such as when the end of inflation can't be found or when H is complex etc. these cases are
    // given a specific flag for the type of problem encountered and logged in the cosmosis datablock.
    try {
        // protected code
        double nEND = model->compute_end_of_inflation(&bkg, Nendhigh);
        transport::basic_range<double> times{N_init, nEND, 500, transport::spacing::linear};
//        std::cout << nEND << std::endl;

        // construct a test twopf task to use with the compute_aH function later
        transport::basic_range<double> k_test{exp(0.0), exp(0.0), 1, transport::spacing::log_bottom};

        transport::twopf_task<DataType> tk2_test{"gelaton.twopf_test", ics, times, k_test};
        tk2_test.set_collect_initial_conditions(true).set_adaptive_ics_efolds(5.0);

        // Use the compute aH method to be able to find the values through-out the duration of inflation.
        // These will be used to find appropriate k values exiting at specific e-foldings by setting k=aH
        // at the desired value of N.
        std::vector<double> N;
        std::vector<double> log_aH;
        std::vector<double> log_a2H2M;
        model->compute_aH(&tk2_test, N, log_aH, log_a2H2M);
        // Interpolate the N and log(aH) values.
        transport::spline1d<double> spline(N, log_aH);

        // Construct some (comoving) k values exiting at at nEND-60, nEND-59, nEND-58, ..., nEND-50.
        // std::vector< std::vector<double> > Nk_values;
        std::vector<double> k_values;
        for (int i = 60; i >= 50; --i)
        {
            double N_value = nEND - i;
            if (N_value < 0.0)
            {
                std::cout << "This inflation model gives less than 60 e-folds of inflation!" << std::endl;
                continue;
            } else {
                double k_value = exp( spline(N_value) );
                // std::cout << i << "\t" << N_value << "\t" << k_value << std::endl;
                k_values.push_back(k_value);
            }
        }

        // To get these k numbers to be conventionally normalised, we need the value of aH given at
        // N_pre as defined above. Divide the k numbers by this value to get conventional normalisation.
        for (int i = 0; i < k_values.size(); ++i)
        {
            k_values[i] = k_values[i] / exp( spline(N_pre) );
            // std::cout << k_values[i] << std::endl;
        }

        // use the vector of k values to build a transport::basic_range object to use for integration
        transport::aggregate_range<double> ks;
        transport::aggregate_range<double> kts;
        for (int i = 0; i < k_values.size(); ++i)
        {
            transport::basic_range<double> k_temp{k_values[i], k_values[i], 1, transport::spacing::log_bottom};
            transport::basic_range<double> kt_temp{3.0*k_values[i], 3.0*k_values[i], 1, transport::spacing::log_bottom};
            ks += k_temp;
            kts += kt_temp;
        }

        // Some alphas and betas needed specifically for equilateral and squeezed configs.
        transport::basic_range<double> alpha_equi{0.0, 0.0, 0, transport::spacing::linear};
        transport::basic_range<double> beta_equi{1.0/3.0, 1.0/3.0, 0, transport::spacing::linear};
        transport::basic_range<double> alpha_sqz{0.0, 0.0, 0, transport::spacing::linear};
        transport::basic_range<double> beta_sqz{0.98, 0.9999, 10, transport::spacing::log_bottom};
        
        // construct a twopf task based on the k values generated above
        transport::basic_range<double> times_sample{nEND-11.0, nEND, 12, transport::spacing::linear};
        tk2 = std::make_unique< transport::twopf_task<DataType> > ("gelaton.twopf", ics, times_sample, ks);
        tk2->set_adaptive_ics_efolds(4.5);

        // construct an equilateral threepf task based on the kt values made above
        // transport::threepf_alphabeta_task<DataType> tk3e{"gelaton.threepf-equilateral", ics, times, kts, alpha_equi, beta_equi};
        // tk3e.set_collect_initial_conditions(true).set_adaptive_ics_efolds(2.0);

        // construct a squeezed threepf task based on the kt values made above.
        // transport::threepf_alphabeta_task<DataType> tk3s{"gelaton.threepf-squeezed", ics, times, kts, alpha_sqz, beta_sqz};
        // tk3s.set_collect_initial_conditions(true).set_adaptive_ics_efolds(3.0);
        
        // TODO: remove the paired zeta tasks when I'm sure the unique pointer to tk2 is working properly
        //       our intention is to handle the gauge transformations within the batcher file.
        // construct a pair of zeta task based on the two- and three-point functions above.
        // transport::zeta_twopf_task<DataType> ztk2{"gelaton.twopf-zeta", *tk2};
        // transport::zeta_threepf_task<DataType> ztk3e{"gelaton.threepf-zeta-equilateral", tk3e};
        // transport::zeta_threepf_task<DataType> ztk3s{"gelaton.threepf-zeta-squeezed", tk3s};
        // ztk2.set_paired(true);
        // ztk3e.set_paired(true);
        // ztk3s.set_paired(true);

        // integrate our tasks created above
        // Add a batcher here to collect the data
        std::vector<double> samples;
        boost::filesystem::path lp(boost::filesystem::current_path());
        unsigned int w;
        sampling_integration_batcher batcher(samples, lp, w, model.get(), tk2.get());
        
        // Integrate all of the twopf samples provided above in the tk2 task
        auto db = tk2->get_twopf_database();
        for (auto t = db.record_cbegin(); t != db.record_cend(); ++t)
        {
            model->twopf_kmode(*t, tk2.get(), batcher, 1);
        }

        for (int i = 0; i < samples.size(); ++i)
        {
            std::cout << "Sample no: " << i << " - " << samples[i] << std::endl;
        }

    } catch(transport::end_of_inflation_not_found& xe) {
        // catch 1
        std::cout << "!!! END OF INFLATION NOT FOUND !!!" << std::endl;
    } catch(transport::Hsq_is_negative& xe) {
        // catch 2
        std::cout << "!!! HSQ IS NEGATIVE !!!" << std::endl;
    } catch(transport::integration_produced_nan& xe) {
        // catch 3
        std::cout << "!!! INTEGRATION PRODUCED NAN !!!" << std::endl;
    } catch(transport::no_massless_time& xe) {
        // catch 4
        std::cout << "!!! NO MASSLESS TIME FOR THIS K MODE !!!" << std::endl;
    } catch(transport::eps_is_negative& xe) {
        // catch 5
        std::cout << "!!! EPSILON PARAMATER IS NEGATIVE !!!" << std::endl;
    } catch(transport::eps_too_large& xe) {
        // catch 6
        std::cout << "!!! EPSILON > 3 !!!" << std::endl;
    } catch(transport::V_is_negative& xe) {
        // catch 7
        std::cout << "!!! NEGATIVE POTENTIAL !!!" << std::endl;
    } catch(transport::failed_to_compute_horizon_exit& xe) {
        // catch 8
        std::cout << "!!! FAILED TO COMPUTE HORIZON EXIT FOR ALL K MODES !!!" << std::endl;
    } catch(transport::adaptive_ics_before_Ninit& xe) {
        // catch 9
        std::cout << "!!! THE ADAPTIVE INITIAL CONDITIONS REQUIRE AN INTEGRATION TIME BEFORE N_INITIAL !!!" <<  std::endl;
    }

    DATABLOCK_STATUS status = DBS_SUCCESS;
    return status;
}


int cleanup(void * config)
{
    // Config is whatever you returned from setup above
    // Free it 
    model.release();
}

} // end of extern C
