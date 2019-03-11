// INCLUDES
// Cosmosis includes
#include <cosmosis/datablock/datablock.hh>
#include <cosmosis/datablock/section_names.h>
// CppTransport includes
#include "gelaton_mpi.h"
#include "transport-runtime/utilities/spline1d.h"
#include "transport-runtime/tasks/integration_detail/abstract.h"
#include "transport-runtime/models/advisory_classes.h"
#include "transport-runtime/tasks/integration_detail/twopf_task.h"
#include "transport-runtime/enumerations.h"
#include "transport-runtime/tasks/integration_detail/twopf_db_task.h"
// Other includes
#include <memory>
#include "math.h"
#include <exception>
#include "boost/filesystem.hpp"
#include <boost/math/interpolators/barycentric_rational.hpp>
#include <boost/range/adaptors.hpp>
#include <boost/math/tools/roots.hpp>
// Include batcher file for integration
#include "sampling_integration_batcher.h"

// Conveniences for CppTransport TODO: remove these - in practice we're always using std::vector<double> or doubles.
using DataType = double;
using StateType = std::vector<DataType>;

namespace inflation {
    // Some names for different sections in the cosmosis datablock
    const char *sectionName = "cppt_sample";
    const char *paramsSection = "inflation_parameters";
    const char *twopf_name = "twopf_observables";
    const char *thrpf_name = "thrpf_observables";

    // These are the Lagrangian parameters for our model (including field initial conditions).
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

    // Unsigned ints for capturing failed samples
    unsigned int no_end_inflate, neg_Hsq, integrate_nan, zero_massless, neg_epsilon, large_epsilon, neg_V,
    failed_horizonExit, ics_before_start, inflate60, time_var_pow_spec;
}

// exception for catching when <60 e-folds inflation given as we need at least 60 e-folds for sampling
struct le60inflation : public std::exception {
    const char * what () const throw () {
        return "<60 e-folds!";
    }
};

// exception for catching a power spectrum with too much time-dependence on super-horizon scales
struct time_varying_spectrum : public std::exception {
    const char * what () const throw () {
        return "time varying spectrum";
    }
};

// definition of tolerance for the bisection of physical k values
struct ToleranceCondition {
    bool operator () (double min, double max) {
        return abs(min - max) <= 1E-6;
    }
};

// Set-up a bisection function using a spline to extract a value of N_exit from some desired value of phys_k
double compute_Nexit_for_physical_k (double Phys_k, transport::spline1d<double>& matching_eq, ToleranceCondition tol)
{
    matching_eq.set_offset(std::log(Phys_k));
    std::string task_name = "find N_exit of physical wave-number";
    std::string bracket_error = "extreme values of N didn't bracket the exit value";
    double Nexit;
    Nexit = transport::task_impl::find_zero_of_spline(task_name, bracket_error, matching_eq, tol);
    return Nexit;
}

static transport::local_environment env;
static transport::argument_cache arg;
static std::unique_ptr< transport::gelaton_mpi<DataType, StateType> > model;
static std::unique_ptr< transport::twopf_task<DataType> > tk2;
static std::unique_ptr< transport::threepf_alphabeta_task<DataType> > tk3e;

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
    // Initialise DATABLOCK_STATUS to 0 - this is returned at end of function
    DATABLOCK_STATUS status = (DATABLOCK_STATUS)0;

    // Read in inflation parameters (Lagrangian and field initial values)
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

    // Print out of each of the inflation parameters read-in above
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

    // Set-up initial time for integration (N_init) and N_pre which is used to set the amount of sub-horizon evolution
    // to integrate before the chosen mode crosses the horizon.
    const double N_init        = 0.0;
    const double N_pre         = 15.0;

    // Create the parameters and initial_conditions objects that CppTransport needs using
    // the two functions defined above.
    using namespace inflation;
    transport::parameters<DataType> params{M_P, inflation::parameter_creator(R0, V0, eta_R, g_R, lambda_R, alpha),
                                           model.get()};
    transport::initial_conditions<DataType> ics{"gelaton", params, inflation::init_cond_creator(R_init, theta_init,
            Rdot_init, thetadot_init), N_init, N_pre};
    
    // Use a silly end value to find nEND and set the time range used by CppT to finish at nEND.
    double Nendhigh = 10000;
    transport::basic_range<DataType> dummy_times{N_init, Nendhigh, 2, transport::spacing::linear};
    transport::background_task<DataType> bkg{ics, dummy_times};

    //! Declare the std::vectors needed for storing the observables - here so it's accessible outside of the try-block
    // Twopf observables
    std::vector<double> k_values;
    std::vector<double> A_s;
    std::vector<double> A_t;
    std::vector<double> r;
    std::vector<double> n_s;
    std::vector<double> n_t;
    std::vector<double> n_s_nonLog;
    std::vector<double> n_t_nonLog;
    // Threepf observables
    std::vector<double> B_equi;
    std::vector<double> fNL_equi;

    // From here, we need to enclose the rest of the code in a try-catch statement in order to catch
    // when a particular set of initial conditions fails to integrate or if there is a problem with the
    // physics such as when the end of inflation can't be found or when H is complex etc. these cases are
    // given a specific flag for the type of problem encountered and logged in the cosmosis datablock.
    try {
        // compute nEND-> throw exception struct defined above if we have nEND < 60.0 e-folds
        double nEND = model->compute_end_of_inflation(&bkg, Nendhigh);
        if (nEND < 60.0)
        {
            throw le60inflation();
        }

        // construct a test twopf task to use with the compute_aH function later
        transport::basic_range<double> times{N_init, nEND, 500, transport::spacing::linear};
        transport::basic_range<double> k_test{exp(0.0), exp(0.0), 1, transport::spacing::log_bottom};
        transport::twopf_task<DataType> tk2_test{"gelaton.twopf_test", ics, times, k_test};
        tk2_test.set_collect_initial_conditions(true).set_adaptive_ics_efolds(5.0);

        //! Matching equation for physical wave-numbers
        // Compute the log(H) values across the inflation time range.
        std::vector<double> N_H;
        std::vector<double> log_H;
        model->compute_H(&tk2_test, N_H, log_H);

        // Set-up the different parameters needed for the matching equation
        std::vector<double> N_reversed (N_H.rbegin(), N_H.rend()); // matching equation needs e-folds in reverse order
        std::vector<double> H_reversed (log_H.rbegin(), log_H.rend()); // matching eqn also needs reversed H
        double Hend = log_H.back(); // value of H at the end of inflation
        double norm_const = std::log(243.5363); // dimnless matching equation has constant = M_P / 1E16 GeV
        double k_pivot = 0.05; // pivot scale defined as 0.05 h^-1 Mpc^-1 here
        double e_fold_const = 55.75; // constant defined in the matching eq.

        // Set-up the matching equation solutions across the inflation time range.
        std::vector<double> log_physical_k (H_reversed.size());
        for (int i = 0; i < H_reversed.size(); ++i)
        {
            // log_physical_k[i] = k_pivot * std::exp(e_fold_const - N_reversed[i] - norm_const + log_H[i] - Hend); // exp version
            log_physical_k[i] = std::log(k_pivot) + e_fold_const - N_reversed[i] - norm_const + H_reversed[i] - Hend;
            // std::cout << "log k[" << i << "] = " << log_physical_k[i] << std::endl;
        }

        // Set-up a spline to use with the bisection method defined later.
        transport::spline1d<double> spline_match_eq (N_reversed, log_physical_k);
        // Set-up a tolerance condition for using with the bisection function
        ToleranceCondition tol;
        // Use the bisection method to find the e-fold exit of k pivot.
        double N_pivot_exit = compute_Nexit_for_physical_k(0.05, spline_match_eq, tol);
        std::cout << "e-fold exit for k* is: " << N_pivot_exit << std::endl;

        // TODO: change the construction of k values below to be finding them from the exit times of physical wave numbers in [1e-4, 1] Mpc^(-1).

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
        for (int i = 60; i >= 50; --i)
        {
            double N_value = nEND - i;
            double k_value = exp( spline(N_value) );
            // std::cout << i << "\t" << N_value << "\t" << k_value << std::endl;
            k_values.push_back(k_value);
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
        tk3e = std::make_unique< transport::threepf_alphabeta_task<DataType> > ("gelaton.threepf-equilateral", ics,
                times_sample, kts, alpha_equi, beta_equi);
        tk3e->set_adaptive_ics_efolds(4.5);

        // construct a squeezed threepf task based on the kt values made above.
        // transport::threepf_alphabeta_task<DataType> tk3s{"gelaton.threepf-squeezed", ics, times, kts, alpha_sqz, beta_sqz};
        // tk3s.set_collect_initial_conditions(true).set_adaptive_ics_efolds(3.0);

        //! INTEGRATE OUR TASKS CREATED FOR THE TWO-POINT FUNCTION ABOVE
        // Add a 2pf batcher here to collect the data - this needs a vector to collect the zeta-twopf samples, a boost
        // filesystem path for logging and an unsigned int for logging as well.
        std::vector<double> samples;
        std::vector<double> tens_samples_twpf;
        boost::filesystem::path lp(boost::filesystem::current_path());
        unsigned int w;
        twopf_sampling_batcher batcher(samples, tens_samples_twpf, lp, w, model.get(), tk2.get());
        
        // Integrate all of the twopf samples provided above in the tk2 task - this is working with the new batcher!
        auto db = tk2->get_twopf_database();
        for (auto t = db.record_cbegin(); t != db.record_cend(); ++t)
        {
            model->twopf_kmode(*t, tk2.get(), batcher, 1);
        }

//        for (int i = 0; i < tens_samples_twpf.size(); ++i)
//        {
//            std::cout << "Sample no: " << i << " :-. Zeta 2pf: " << samples[i] << " ; Tensor 2pf: " << tens_samples_twpf[i] << std::endl;
//        }

        // find the std deviation & mean of the power spectrum amplitudes
        std::vector<double> mean(11), std_dev(11);
        for (int i = 0; i < k_values.size(); i++)
        {
            double sum = 0;
            for (int j = 0; j < 13; j++)
            {
                sum += samples[(13*i)+j];
            }
            mean[i] = sum / 13.0;
        }

        for (int i = 0; i < k_values.size(); i++)
        {
            double sum_sq = 0;
            for (int j = 0; j < 13; j++)
            {
                sum_sq += pow(samples[(13*i)+j] - mean[i], 2);
            }
            std_dev[i] = sqrt(sum_sq / (times_sample.size() - 1)); // divide by 12 for N-1 samples
        }

        // find a measure of the dispersion of power spectrum values -> std-dev/mean
        std::vector<double> dispersion(11);
        for (int i = 0; i < mean.size(); ++i)
        {
            dispersion[i] = std_dev[i]/mean[i];
//            std::cout << dispersion[i] << std::endl;
        }

        // throw the time-varying exception defined above if the dispersion is >10%
        for (auto i: dispersion)
        {
            if (i > 0.1)
            {
                throw time_varying_spectrum();
            }
        }

        // find A_s & A_t for each k mode exiting at Nend-10, ..., Nend etc. We take the final time value at Nend to be
        // the amplitude for the scalar and tensor modes. The tensor-to-scalar ratio r is the ratio of these values.
        for (int k = 0; k < k_values.size(); ++k)
        {
            int index = (times_sample.size() * k) + (times_sample.size() -1);
            A_s.push_back(samples[index]);
            A_t.push_back(tens_samples_twpf[index]);
            r.push_back( tens_samples_twpf[index] / samples[index] );
        }

//        for (int i=0; i < A_s.size(); i++)
//        {
//            std::cout << "A_s: " << A_s[i] << std::endl;
//            std::cout << "A_t: " << A_t[i] << std::endl;
//            std::cout << "r: " << r[i] << std::endl;
//        }

        //! Construct log values of A and k for the n_s & n_t splines TODO: compare to non-log diff
        // construct two vectors of log A_s & log A_t values
        std::vector<double> logA_s;
        std::vector<double> logA_t;
        for (int i = 0; i < A_s.size(); ++i)
        {
            logA_s.push_back( log( A_s[i] ) );
            logA_t.push_back( log( A_t[i] ) );
        }
        // construct a vector of ln(k) values
        std::vector<double> logK;
        for (auto& i: k_values)
        {
            logK.push_back(log(i));
        }

        // construct two splines for the log k and log A_s or A_t values
        transport::spline1d<double> ns_spline( logK, logA_s );
        transport::spline1d<double> nt_spline( logK, logA_t );

        // use the eval_diff method in the splines to compute n_s and n_t for each k value
        for (int j = 0; j < logK.size(); ++j)
        {
            double temp = ns_spline.eval_diff(logK[j]) + 1.0; // add 1 for the n_s-1 scalar index convention
            double temp2 = nt_spline.eval_diff(logK[j]);
            n_s.push_back(temp);
            n_t.push_back(temp2);
            std::cout << "n_s for " << logK[j] << " is: " << temp << std::endl;
            std::cout << "n_t for " << logK[j] << " is: " << temp2 << std::endl;
        }

        // Repeat differentiation without using logs
        transport::spline1d<double> ns_spline_nonlog(k_values, A_s);
        transport::spline1d<double> nt_spline_nonlog(k_values, A_t);

        // use the eval_diff method in the splines to compute n_s and n_t for each k value
        for (int i = 0; i < k_values.size(); ++i)
        {
            double temp = ns_spline_nonlog.eval_diff(k_values[i]) * (k_values[i]/A_s[i]) + 1.0;
            double temp2 = nt_spline_nonlog.eval_diff(k_values[i]) * (k_values[i]/A_t[i]);
            n_s_nonLog.emplace_back(temp);
            n_t_nonLog.emplace_back(temp2);
            std::cout << "n_s (non-log) for " << logK[i] << " is: " << temp << std::endl;
            std::cout << "n_t (non-log) for " << logK[i] << " is: " << temp2 << std::endl;
        }

        //! INTEGRATE OUR TASKS CREATED FOR THE EQUILATERAL 3-POINT FUNCTION ABOVE
        // Add a 3pf batcher here to collect the data - this needs 3 vectors for the z2pf, z3pf and redbsp data samples
        // as well as the same boost::filesystem::path and unsigned int variables used in the 2pf batcher.
        std::vector<double> twopf_samples;
        std::vector<double> tens_samples_thpf;
        std::vector<double> threepf_samples;
        std::vector<double> redbsp_samples;
        threepf_sampling_batcher thpf_batcher(twopf_samples, tens_samples_thpf, threepf_samples, redbsp_samples, lp, w, model.get(), tk3e.get());

        // Integrate all of threepf samples provided in the tk3e task
        auto db2 = tk3e->get_threepf_database();
        for (auto t = db2.record_cbegin(); t!= db2.record_cend(); ++t)
        {
            model->threepf_kmode(*t, tk3e.get(), thpf_batcher, 1);
        }

//        for (auto i = 0; i < threepf_samples.size(); i++)
//        {
//            std::cout << "Threepf sample no: " << i << " - " << threepf_samples[i] << " ; Redbsp: " << redbsp_samples[i] << std::endl;
//        }

        // find the bispectrum amplitude and f_NL amplitude at the end of inflation for each k mode
        for (int j = 0; j < kts.size(); j++)
        {
            int index = (times_sample.size() * j) + (times_sample.size() -1);
            B_equi.emplace_back( threepf_samples[index] );
            fNL_equi.emplace_back( redbsp_samples[index] );
        }

    } catch (transport::end_of_inflation_not_found& xe) {
        std::cout << "!!! END OF INFLATION NOT FOUND !!!" << std::endl;
        inflation::no_end_inflate = 1;
    } catch (transport::Hsq_is_negative& xe) {
        std::cout << "!!! HSQ IS NEGATIVE !!!" << std::endl;
        inflation::neg_Hsq = 1;
    } catch (transport::integration_produced_nan& xe) {
        std::cout << "!!! INTEGRATION PRODUCED NAN !!!" << std::endl;
        inflation::integrate_nan = 1;
    } catch (transport::no_massless_time& xe) {
        std::cout << "!!! NO MASSLESS TIME FOR THIS K MODE !!!" << std::endl;
        inflation::zero_massless = 1;
    } catch (transport::eps_is_negative& xe) {
        std::cout << "!!! EPSILON PARAMATER IS NEGATIVE !!!" << std::endl;
        inflation::neg_epsilon = 1;
    } catch (transport::eps_too_large& xe) {
        std::cout << "!!! EPSILON > 3 !!!" << std::endl;
        inflation::large_epsilon = 1;
    } catch (transport::V_is_negative& xe) {
        std::cout << "!!! NEGATIVE POTENTIAL !!!" << std::endl;
        inflation::neg_V = 1;
    } catch (transport::failed_to_compute_horizon_exit& xe) {
        std::cout << "!!! FAILED TO COMPUTE HORIZON EXIT FOR ALL K MODES !!!" << std::endl;
        inflation::failed_horizonExit = 1;
    } catch (transport::adaptive_ics_before_Ninit& xe) {
        std::cout << "!!! THE ADAPTIVE INITIAL CONDITIONS REQUIRE INTEGRATION TIME BEFORE N_INITIAL !!!" <<  std::endl;
        inflation::ics_before_start = 1;
    } catch (le60inflation& xe) {
        std::cout << "!!! WE HAVE LESS THAN 60 E-FOLDS OF INFLATION !!!" << std::endl;
        inflation::inflate60 = 1;
    } catch (time_varying_spectrum& xe) {
        std::cout << "!!! THE POWER SPECTRUM VARIES TOO MUCH !!!" << std::endl;
        inflation::time_var_pow_spec = 1;
    }

    //! Return the calculated observables to the datablock
    // Use the put_val method to add second-order observables (A_s, k_values, A_t, n_s, n_t & r to the datablock
    // currently only adding the values for the first (index 0) item for a k mode exiting at nEND-60.
    status = block->put_val( inflation::twopf_name, "A_s", A_s[0] );
    status = block->put_val( inflation::twopf_name, "k_values", k_values[0] );
    status = block->put_val( inflation::twopf_name, "A_t", A_t[0] );
    status = block->put_val( inflation::twopf_name, "n_s", n_s[0] );
    status = block->put_val( inflation::twopf_name, "n_t", n_t[0] );
    status = block->put_val( inflation::twopf_name, "r", r[0] );
    // Use put_val to put the three-point observables (B_equi, fNL_equi) onto the datablock
    status = block->put_val( inflation::thrpf_name, "B_equi", B_equi[0] );
    status = block->put_val( inflation::thrpf_name, "fNL_equi", fNL_equi[0] );

    // return status variable declared at the start of the function
    return status;
}


int cleanup(void * config)
{
    // Config is whatever you returned from setup above
    // Free it 
    model.release();
}

} // end of extern C
