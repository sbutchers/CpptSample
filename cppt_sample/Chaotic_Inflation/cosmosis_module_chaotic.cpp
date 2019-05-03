// INCLUDES
// Cosmosis includes
#include <cosmosis/datablock/datablock.hh>
#include <cosmosis/datablock/section_names.h>
// CppTransport includes
#include "chaotic_mpi.h"
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
#include <boost/range/adaptors.hpp>
#include <boost/math/tools/roots.hpp>
// Include batcher file for integration
#include "sampling_integration_batcher.h"

namespace inflation {
    // Some names for different sections in the cosmosis datablock
    const char *sectionName = "cppt_sample";
    const char *paramsSection = "inflation_parameters";
    const char *twopf_name = "twopf_observables";
    const char *thrpf_name = "thrpf_observables";
    const char *spec_file = "wavenumber_spectrum";
    const char *fail_names = "failed_samples";

    // These are the Lagrangian parameters for our model (including field initial conditions).
    double M_P, m, phi_init, phiDot_init;

    // no. of k samples for CLASS read in from ini file
    int num_k_samples;

    // function to create a transport::parameters<double> object called params
    std::vector<double> parameter_creator(double m) {
        std::vector<double> output{m};
        return output; 
    }

    // function to create a std::vector<double> vector containing initial conditons
    std::vector<double> init_cond_creator(double phi_init, double phiDot_init) {
        std::vector<double> output{phi_init, phiDot_init};
        return output;
    }

    // ints for capturing failed samples
    int no_end_inflate = 0, neg_Hsq = 0, integrate_nan = 0, zero_massless = 0, neg_epsilon = 0, 
        large_epsilon = 0, neg_V = 0, failed_horizonExit = 0, ics_before_start = 0, inflate60 = 0, 
        time_var_pow_spec = 0;
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
        return abs(min - max) <= 1E-12;
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
    matching_eq.set_offset(0.0);
    return Nexit;
}

// Set-up a function to create a log-spaced std::vector similar to numpy.logspace
class Logspace {
private:
    double curValue, base;

public:
    Logspace(double first, double base) : curValue(first), base(base) {}

    double operator()() {
        double retval = curValue;
        curValue *= base;
        return retval;
    }
};
std::vector<double> pyLogspace (double start, double stop, int num, double base = 10) {
    double logStart = pow(base, start);
    double logBase = pow(base, (stop-start)/num);

    std::vector<double> log_vector;
    log_vector.reserve(num+1);
    std::generate_n(std::back_inserter(log_vector), num+1, Logspace(logStart, logBase));
    return log_vector;
}

// Set-up a function that returns the k-derivative of A_s or A_t as a double.
double spec_derivative(const double k, const double dk, std::vector<double> Amplitude)
{
    // compute d/dk[Amplitude] using a three-point central difference of O(dk^6).
    const double dk1 = dk;
    const double dk2 = dk1*2.0;
    const double dk3 = dk1*3.0;

    const double Amp1 = (Amplitude[4] - Amplitude[2]) / 2;
    const double Amp2 = (Amplitude[5] - Amplitude[1]) / 4;
    const double Amp3 = (Amplitude[6] - Amplitude[0]) / 6;

    const double fifteen_Amp1 = 15 * Amp1;
    const double six_Amp2 = 6 * Amp2;
    const double ten_dk1 = 10 * dk1;

    const double dAdk = ((fifteen_Amp1 - six_Amp2) + Amp3) / ten_dk1;
    const double k_div_A = k / Amplitude[3];
    return (dAdk * k_div_A);
}

// Set-up a dispersion class that has a function which checks the power spectrum values contained in samples with k and
// t samples stored in k_size and time_size for strongly-varying spectrum values. This is when the std deviation of
// the times samples for a given k_sample is 10% of of the mean value.
class dispersion
{
// Constructors etc.
public:
    // Constructor for transport::basic_range k samples
    dispersion(transport::basic_range<double>& k_samples, transport::basic_range<double>& time_samples,
            std::vector<double>& spectrum_samples)
            : k_size(k_samples.size()),
              time_size(time_samples.size()),
              samples(spectrum_samples)
    {
    }
    // Constructor for transport::aggregate_range k samples
    dispersion(transport::aggregate_range<double>& k_samples, transport::basic_range<double>& time_samples,
            std::vector<double>& spectrum_samples)
            : k_size(k_samples.size()),
              time_size(time_samples.size()),
              samples(spectrum_samples)
    {
    }
    // move constructor
    dispersion(dispersion&&) = default;
    // destructor
    virtual ~dispersion() = default;

// Dispersion calculation
public:
    bool dispersion_check()
    {
        // find the mean of the power spectrum amplitudes
        std::vector<double> mean(k_size), std_dev(k_size);
        for (int i = 0; i < k_size; i++)
        {
            double sum = 0;
            for (int j = 0; j < time_size; j++)
            {
                sum += samples[(time_size*i)+j];
            }
            mean[i] = sum / time_size;
        }

        // find sum_square values for standard deviation
        for (int i = 0; i < k_size; i++)
        {
            double sum_sq = 0;
            for (int j = 0; j < time_size; j++)
            {
                sum_sq += pow(samples[(time_size*i)+j] - mean[i], 2);
            }
            std_dev[i] = sqrt(sum_sq / (time_size - 1)); // divide by time_size-1 for N-1 samples
        }

        // find a measure of the dispersion of power spectrum values -> std-dev/mean
        std::vector<double> dispersion(k_size);
        for (int i = 0; i < mean.size(); ++i)
        {
            dispersion[i] = std_dev[i]/mean[i];
        }

        // return true if the dispersion is >10% for any of the k samples
        for (auto i: dispersion)
        {
            if (i > 0.1)
            {
                return true;
            }
        }

        return false;
    }

// Internal data
protected:
    size_t k_size;
    size_t time_size;
    std::vector<double>& samples;
};

// Create instances of the model and separate integration tasks for the two-point function -> a sampling one and a task
// at k_pivot=0.05Mpc^(-1), and three-point function task with k=k_pivot for the equilateral and squeezed configurations
static transport::local_environment env;
static transport::argument_cache arg;
static std::unique_ptr< transport::chaotic_mpi<double, std::vector<double>> > model;
static std::unique_ptr< transport::twopf_task<double> > tk2;
static std::unique_ptr< transport::twopf_task<double> > tk2_piv;
static std::unique_ptr< transport::threepf_alphabeta_task<double> > tk3e;
static std::unique_ptr< transport::threepf_alphabeta_task<double> > tk3s;

extern "C" {

void * setup(cosmosis::DataBlock * options)
{
    // Read options from the CosmoSIS configuration ini file,
    // passed via the "options" argument
    options->get_val(inflation::sectionName, "M_P", inflation::M_P); // TODO: get rid of this option?
    options->get_val(inflation::sectionName, "k_samples", inflation::num_k_samples);

    // Record any configuration information required
    model = std::make_unique< transport::chaotic_mpi<double, std::vector<double>> > (env, arg);

    // Pass back any object you like
    return options;
}

DATABLOCK_STATUS execute(cosmosis::DataBlock * block, void * config)
{
    // Initialise DATABLOCK_STATUS to 0 - this is returned at end of function
    DATABLOCK_STATUS status = (DATABLOCK_STATUS)0;
    // Add failure DATABLOCK_STATUS variable - returned if an integration fails.
    const DATABLOCK_STATUS failure = (DATABLOCK_STATUS)1;

    //! Read in inflation parameters (Lagrangian and field initial values)
    block->get_val(inflation::paramsSection, "m", inflation::m);
    block->get_val(inflation::paramsSection, "phi_init", inflation::phi_init);
    block->get_val(inflation::paramsSection, "phi_dot_init", inflation::phiDot_init);

    // Print out of each of the inflation parameters read-in above
    std::cout << "m = " << inflation::m << std::endl;
    std::cout << "phi_init = " << inflation::phi_init << std::endl;
    std::cout << "phi_dot_init = " << inflation::phiDot_init << std::endl;

    // Set-up initial time for integration (N_init) and N_pre which is used to set the amount of sub-horizon evolution
    // to integrate before the chosen mode crosses the horizon.
    const double N_init        = 0.0;
    const double N_pre         = 15.0;

    // Create the parameters and initial_conditions objects that CppTransport needs using
    // the two functions defined above.
    transport::parameters<double> params{inflation::M_P, inflation::parameter_creator(inflation::m), model.get()};
    transport::initial_conditions<double> ics{"chaotic", params, inflation::init_cond_creator(inflation::phi_init,
            inflation::phiDot_init), N_init, N_pre};
    
    // Use a silly end value to find nEND and set the time range used by CppT to finish at nEND.
    double Nendhigh = 10000;
    transport::basic_range<double> dummy_times{N_init, Nendhigh, 2, transport::spacing::linear};
    transport::background_task<double> bkg{ics, dummy_times};

    //! Declare the doubles needed for storing the observables - here so it's accessible outside of the try-block
    // Pivot task observables
    // Twopf observables
    double k_pivot_cppt;
    double N_pivot_exit;
    double A_s_pivot;
    double A_t_pivot;
    double r_pivot;
    double ns_pivot;
    double nt_pivot;
    std::vector<double> r;
    // Threepf observables (at pivot scale)
    double B_equi_piv;
    double fNL_equi_piv;
    double B_squ_piv;
    double fNL_squ_piv;

    //! Objects needed for creating & storing the big twopf task for passing to a Boltzmann code.
    // Wavenumber k vectors for passing to CLASS, CAMB or another Boltzmann code
    // Use the pyLogspace function to produce log-spaced values between 10^(-6) & 10^(0) Mpc^(-1) with the number of
    // k samples given in 'num_k_samples' read-in above.
    std::vector<double> Phys_waveno_sample = pyLogspace(-6.0, 1.7, inflation::num_k_samples, 10);
    std::vector<double> k_conventional(Phys_waveno_sample.size());
    // Vectors for storing A_s and A_t before writing them to a temporary file
    std::vector<double> A_s;
    std::vector<double> A_t;

    // From here, we need to enclose the rest of the code in a try-catch statement in order to catch
    // when a particular set of initial conditions fails to integrate or if there is a problem with the
    // physics such as when the end of inflation can't be found or when H is complex etc. these cases are
    // given a specific flag for the type of problem encountered and logged in the cosmosis datablock.
    try {
        //! compute nEND-> throw exception struct defined above if we have nEND < 60.0 e-folds
        double nEND = model->compute_end_of_inflation(&bkg, Nendhigh);
        std::cout << "Inflation lasts for: " << nEND << " e-folds." << std::endl;
        if (nEND < 60.0)
        {
            throw le60inflation();
        }

        //! construct a test twopf task to use with the compute_H function later
        transport::basic_range<double> times{N_init, nEND, 500, transport::spacing::linear};
        transport::basic_range<double> k_test{exp(0.0), exp(0.0), 1, transport::spacing::log_bottom};
        transport::twopf_task<double> tk2_test{"chaotic.twopf_test", ics, times, k_test};
        tk2_test.set_collect_initial_conditions(true).set_adaptive_ics_efolds(5.0);

        //! Matching equation for physical wave-numbers
        // Compute the log(H) values across the inflation time range.
        std::vector<double> N_H;
        std::vector<double> log_H;
        model->compute_H(&tk2_test, N_H, log_H);

        // Set-up the different parameters needed for the matching equation
        double Hend = 0.5 * log_H.back(); // value of H at the end of inflation
        double norm_const = std::log(243.5363 * pow(3.0, 0.25)); // dimnless matching eq has const = (3^(1/4)*M_P)/1E16 GeV
        double k_pivot = std::log(0.05); // pivot scale defined as 0.05 Mpc^-1 here
        double e_fold_const = 55.75; // constant defined in the matching eq.
        double constants = e_fold_const + k_pivot + norm_const - Hend; // wrap up constants in a single term

        // Find the matching equation solutions across the inflation time range.
        std::vector<double> log_physical_k (N_H.size());
        for (int i = 0; i < N_H.size(); ++i)
        {
            log_physical_k[i] = log_H[i] - (nEND - N_H[i]) + constants;
        }

        // Set-up a tolerance condition for using with the bisection function
        ToleranceCondition tol;
        // Set-up a spline to use with the bisection method defined later.
        transport::spline1d<double> spline_match_eq (N_H, log_physical_k);

        // Use the bisection method to find the e-fold exit of k pivot.
        N_pivot_exit = compute_Nexit_for_physical_k(0.05, spline_match_eq, tol);
        // std::cout << "e-fold exit for k* is: " << N_pivot_exit << std::endl;
        // std::cout << "k* from spline is:" << spline_match_eq(N_pivot_exit) << std::endl;

        // Construct a vector of exit times (= no. of e-folds BEFORE the end of inflation!)
//        std::vector<double> Phys_k_exits (Phys_waveno_sample.size());
//        for (int i = 0; i < Phys_k_exits.size(); ++i)
//        {
//            Phys_k_exits[i] = compute_Nexit_for_physical_k(Phys_waveno_sample[i], spline_match_eq, tol);
//            std::cout << Phys_waveno_sample[i] << "Mpc^(-1) exits at: " << Phys_k_exits[i] << " e-folds." << std::endl;
//        }

        //! Construct the wave-numbers using a linearity relation.
        // Build CppT normalised wave-numbers by using the linear relation k_phys = gamma * k_cppt and k_cppt[Npre] == 1
        double gamma = spline_match_eq(N_pre);

        for (int i = 0; i < k_conventional.size(); ++i)
        {
            k_conventional[i] = Phys_waveno_sample[i] / exp(gamma);
            //std::cout << Phys_k_exits[i] << "\t" << k_conventional[i] << std::endl;
        }

        // Put these into an aggregate range one-by-one
        transport::aggregate_range<double> ks;
        for (int i = 0; i < k_conventional.size(); ++i)
        {
            transport::basic_range<double> k_temp{k_conventional[i], k_conventional[i], 1, transport::spacing::linear};
            ks += k_temp;
        }

        // Construct a CppT normalised wave-number for the pivot scale (=0.05Mpc^-1) using the linearity constant
        k_pivot_cppt = 0.05 / std::exp(gamma);

        // Use the CppT normalised kpivot value to build a wave-number range for kpivot with some other values to use
        // for finding the spectral indices.
        double dk = 0.0001 * k_pivot_cppt;
        transport::basic_range<double> k_pivot_range{k_pivot_cppt-(3*dk), k_pivot_cppt+(3*dk), 6, transport::spacing::linear};

        // Use the CppT normalised kpivot value to build a range with kt = 3*kpivot only
        transport::basic_range<double> kt_pivot_range{3.0*k_pivot_cppt, 3.0*k_pivot_cppt, 1, transport::spacing::linear};

        //! ################################################################################

        //! OLD WAY OF CONSTRUCTING THE WAVENUMBERS - this found wavenumbers exiting at nEND-60 up to nEND-50 by finding
        //! the value of aH at these times and then rescaling with the value at nPre for CppT normalisation
        // Use the compute aH method to be able to find the values through-out the duration of inflation.
        // These will be used to find appropriate k values exiting at specific e-foldings by setting k=aH
        // at the desired value of N.
//        std::vector<double> N;
//        std::vector<double> log_aH;
//        std::vector<double> log_a2H2M;
//        model->compute_aH(&tk2_test, N, log_aH, log_a2H2M);
//        // Interpolate the N and log(aH) values.
//        transport::spline1d<double> aH_spline(N, log_aH);
//        // Construct some (comoving) k values exiting at at nEND-60, nEND-59, nEND-58, ..., nEND-50.
//        for (int i = 60; i >= 50; --i)
//        {
//            double N_value = nEND - i;
//            double k_value = exp( aH_spline(N_value) );
//            k_values.push_back(k_value);
//        }
//        // To get these k numbers to be conventionally normalised, we need the value of aH given at
//        // N_pre as defined above. Divide the k numbers by this value to get conventional normalisation.
//        for (int i = 0; i < k_values.size(); ++i)
//        {
//            k_values[i] = k_values[i] / exp( aH_spline(N_pre) );
//        }
//        // use the vector of k values to build a transport::basic_range object to use for integration
//        transport::aggregate_range<double> ks;
//        for (int i = 0; i < k_conventional.size(); ++i)
//        {
//            transport::basic_range<double> k_temp{k_conventional[i], k_conventional[i], 1, transport::spacing::log_bottom};
//            ks += k_temp;
//        }
//        // Do the same thing for kt but not as many because 3pf integrations are longer!
//        transport::aggregate_range<double> kts;
//        for (int i = 0; i < k_values.size(); ++i)
//        {
//            transport::basic_range<double> kt_temp{3.0*k_values[i], 3.0*k_values[i], 1, transport::spacing::log_bottom};
//            kts += kt_temp;
//        }

        //! END OF OLD WAY OF CONSTRUCTING WAVENUMBERS

        //! ################################################################################

        //! Construct the integration tasks for the different configurations
        // Some alphas and betas needed specifically for equilateral and squeezed configs.
        transport::basic_range<double> alpha_equi{0.0, 0.0, 0, transport::spacing::linear};
        transport::basic_range<double> beta_equi{1.0/3.0, 1.0/3.0, 0, transport::spacing::linear};
        transport::basic_range<double> alpha_sqz{0.0, 0.0, 0, transport::spacing::linear};
        transport::basic_range<double> beta_sqz{0.98, 0.98, 0, transport::spacing::linear};
        // Set-up a time sample for integrations at the end of inflation so we can extract A_s, A_t etc. while giving
        // a wide enough interval to check the values are stable.
        transport::basic_range<double> times_sample{nEND-11.0, nEND, 12, transport::spacing::linear};
        
        // construct a twopf task based on the k values generated above
        tk2 = std::make_unique< transport::twopf_task<double> > ("chaotic.twopf", ics, times_sample, ks);
        tk2->set_adaptive_ics_efolds(4.5);
        // construct a twopf task for the pivot scale
        tk2_piv = std::make_unique< transport::twopf_task<double> > ("chaotic.twopf-pivot", ics, times_sample, k_pivot_range);
        tk2_piv->set_adaptive_ics_efolds(4.5);
        // construct an equilateral threepf task based on the kt pivot scale made above
        tk3e = std::make_unique< transport::threepf_alphabeta_task<double> > ("chaotic.threepf-equilateral", ics,
                times_sample, kt_pivot_range, alpha_equi, beta_equi);
        tk3e->set_adaptive_ics_efolds(4.5);
        // construct a squeezed threepf task based on the kt pivot scale made above.
        tk3s = std::make_unique< transport::threepf_alphabeta_task<double> > ("chaotic.threepf-squeezed", ics,
                times_sample, kt_pivot_range, alpha_sqz, beta_sqz);
        tk3s->set_adaptive_ics_efolds(4.5);

        //! INTEGRATE OUR TASKS CREATED FOR THE TWO-POINT FUNCTION ABOVE
        // All batchers need the filesystem path and an unsigned int for logging TODO: Double check these are ok to use for every task!
        boost::filesystem::path lp(boost::filesystem::current_path());
        unsigned int w;
        int g = 0;
        bool no_log = true;

        //! Twopf pivot task
        std::vector<double> pivot_twopf_samples;
        std::vector<double> tens_pivot_samples;
        twopf_sampling_batcher pivot_batcher(pivot_twopf_samples, tens_pivot_samples, lp, w, model.get(), tk2_piv.get(), g, no_log);

        // Integrate the pivot task
        auto db_piv = tk2_piv->get_twopf_database();
        for (auto t = db_piv.record_cbegin(); t != db_piv.record_cend(); ++t)
        {
            model->twopf_kmode(*t, tk2_piv.get(), pivot_batcher, 1);
        }

        // Print out of samples for the pivot task
//        for (int i = 0; i < pivot_twopf_samples.size(); ++i)
//        {
//            std::cout << "Sample no: " << i << " :-. Zeta 2pf: " << pivot_twopf_samples[i] << " ; Tensor 2pf: " << tens_pivot_samples[i] << std::endl;
//        }

        std::vector<double> k_pivots;
        for (int i = 0; i < k_pivot_range.size(); ++i)
        {
            k_pivots.push_back(k_pivot_range[i]);
        }

        // Perform a dispersion check on the spectrum values - throw time_varying_spectrum if they're varying.
        dispersion twpf_pivot_dispersion(k_pivot_range, times_sample, pivot_twopf_samples);
        if(twpf_pivot_dispersion.dispersion_check() == true)
        {
            throw time_varying_spectrum();
        }

        // Extract the A_s & a_t values: put the 7 A_s & A_t values into vectors for finding n_s and n_t with, then
        // take the values at index 3 (centre) to get the pivot scale.
        std::vector<double> A_s_spec(k_pivot_range.size());
        std::vector<double> A_t_spec(k_pivot_range.size());
        for (int k = 0; k < k_pivot_range.size(); ++k)
        {
            int index = (times_sample.size() * k) + (times_sample.size() - 1);
            A_s_spec[k] = pivot_twopf_samples[index];
            A_t_spec[k] = tens_pivot_samples[index];
            if (k==3) 
            {
                A_s_pivot = A_s_spec[k];
                A_t_pivot = A_t_spec[k];
                r_pivot = ( A_t_pivot / A_s_pivot );
            }
        }

        // std::cout << "r_pivot is: " << r_pivot << std::endl;
        // std::cout << "A_s (pivot) is: " << A_s_pivot << std::endl;
        // std::cout << "A_t (pivot) is: " << A_t_pivot << std::endl;

        // Use the function defined above to find dA/dk and compute n_s and n_t from those
        ns_pivot = spec_derivative(k_pivot_cppt, dk, A_s_spec) + 1.0;
        nt_pivot = spec_derivative(k_pivot_cppt, dk, A_t_spec);

        transport::spline1d<double> ns_piv_spline(k_pivots, A_s_spec);
        double ns_pivot_spline = ns_piv_spline.eval_diff(k_pivot_cppt) * (k_pivot_cppt / A_s_pivot) + 1.0;

        transport::spline1d<double> nt_piv_spline(k_pivots, A_t_spec);
        double nt_pivot_spline = nt_piv_spline.eval_diff(k_pivot_cppt) * (k_pivot_cppt / A_t_pivot);

        // std::cout << "ns: " << ns_pivot << "\t" << "ns(spline): " << ns_pivot_spline << std::endl;
        // std::cout << "nt: " << nt_pivot << "\t" << "nt(spline): " << nt_pivot_spline << std::endl;

        //! Big twopf task for CLASS or CAMB
        // Add a 2pf batcher here to collect the data - this needs a vector to collect the zeta-twopf samples.
        std::vector<double> samples;
        std::vector<double> tens_samples_twpf;
        twopf_sampling_batcher batcher(samples, tens_samples_twpf, lp, w, model.get(), tk2.get(), g, no_log);

        // Integrate all of the twopf samples provided above in the tk2 task
        auto db = tk2->get_twopf_database();
        for (auto t = db.record_cbegin(); t != db.record_cend(); ++t)
        {
            model->twopf_kmode(*t, tk2.get(), batcher, 1);
        }

//        for (int i = 0; i < tens_samples_twpf.size(); ++i)
//        {
//            std::cout << "Sample no: " << i << " :-. Zeta 2pf: " << samples[i] << " ; Tensor 2pf: " << tens_samples_twpf[i] << std::endl;
//        }

        // Perform a dispersion check on the spectrum values - throw time_varying_spectrum if they're varying.
        dispersion twopf_task_disp(ks, times_sample, samples);
        if ( twopf_task_disp.dispersion_check() == true)
        {
            std::cout << "time-varying spectrum" << std::endl;
            throw time_varying_spectrum();
        }

        // find A_s & A_t for each k mode exiting at Nend-10, ..., Nend etc. We take the final time value at Nend to be
        // the amplitude for the scalar and tensor modes. The tensor-to-scalar ratio r is the ratio of these values.
        for (int k = 0; k < k_conventional.size(); ++k)
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

        //! Integrate the tasks created for the equilateral 3-point function above
        // Add a 3pf batcher here to collect the data - this needs 3 vectors for the z2pf, z3pf and redbsp data samples
        // as well as the same boost::filesystem::path and unsigned int variables used in the 2pf batcher.
       std::vector<double> eq_twopf_samples;
       std::vector<double> eq_tens_samples;
       std::vector<double> eq_threepf_samples;
       std::vector<double> eq_redbsp_samples;
       threepf_sampling_batcher eq_thpf_batcher(eq_twopf_samples, eq_tens_samples, eq_threepf_samples, eq_redbsp_samples,
               lp, w, model.get(), tk3e.get(), g, no_log);

       // Integrate all of the threepf samples provided in the tk3e task
       auto eq_db = tk3e->get_threepf_database();
       for (auto t = eq_db.record_cbegin(); t!= eq_db.record_cend(); ++t)
       {
           model->threepf_kmode(*t, tk3e.get(), eq_thpf_batcher, 1);
       }

       // Print-out of threepf samples
    //    for (auto i = 0; i < eq_threepf_samples.size(); i++)
    //    {
    //        std::cout << "Threepf sample no: " << i << " - " << eq_threepf_samples[i] << " ; Redbsp: " << eq_redbsp_samples[i] << std::endl;
    //    }

       // Perform a dispersion check - throw time_varying_spectrum if spectra aren't stable
       dispersion equi_B_disp_check(kt_pivot_range, times_sample, eq_threepf_samples);
       dispersion equi_fNL_disp_check(kt_pivot_range, times_sample, eq_redbsp_samples);
       if ( (equi_B_disp_check.dispersion_check() == true) or (equi_fNL_disp_check.dispersion_check() == true) ) {
           throw time_varying_spectrum();
       }

       // find the bispectrum amplitude and f_NL amplitude at the end of inflation for the pivot scale
       // do this by taking the value at the end of inflation
       B_equi_piv = eq_threepf_samples.back();
       fNL_equi_piv = eq_redbsp_samples.back();

        //! Integrate the task for the squeezed 3-point function above
//        // Add a 3pf batcher here to collect the data - this needs 3 vectors for the z2pf, z3pf and redbsp data samples
//        // as well as the same boost::filesystem::path and unsigned int variables used in the 2pf batcher.
//        std::vector<double> sq_twopf_samples;
//        std::vector<double> sq_tens_samples;
//        std::vector<double> sq_threepf_samples;
//        std::vector<double> sq_redbsp_samples;
//        threepf_sampling_batcher sq_thpf_batcher(sq_twopf_samples, sq_tens_samples, sq_threepf_samples, sq_redbsp_samples,
//                                                 lp, w, model.get(), tk3s.get());
//
//        // Integrate all of the threepf samples provided in the tk3s task
//        auto sq_db = tk3s->get_threepf_database();
//        for (auto t = sq_db.record_cbegin(); t!= sq_db.record_cend(); ++t)
//        {
//            model->threepf_kmode(*t, tk3s.get(), sq_thpf_batcher, 1);
//        }
//
//        // Print-out of squeezed threepf data
//        for (auto i = 0; i < sq_threepf_samples.size(); i++)
//        {
//            std::cout << "Squeezed threepf sample no: " << i << " - " << sq_threepf_samples[i] << " ; Redbsp: " << sq_redbsp_samples[i] << std::endl;
//        }
//
//        // Perform a dispersion check - throw time_varying_spectrum if spectra aren't stable
//        dispersion sq_B_disp_check(kt_pivot_range, times_sample, sq_threepf_samples);
//        dispersion sq_fNL_disp_check(kt_pivot_range, times_sample, sq_redbsp_samples);
//        if ( (sq_B_disp_check.dispersion_check() == true) or (sq_fNL_disp_check.dispersion_check() == true) ) {
//            throw time_varying_spectrum();
//        }
//
//        // find the bispectrum amplitude and f_NL amplitude at the end of inflation for the pivot scale
//        B_squ_piv = sq_threepf_samples.back();
//        fNL_squ_piv = sq_redbsp_samples.back();

    // Begin catches for different exceptions thrown from a failed integration sample.
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
        std::cout << "!!! EPSILON PARAMETER IS NEGATIVE !!!" << std::endl;
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
    } // end of try-catch block

    // FAILED SAMPLE INFO
    // Use put_val to write the information about any caught exceptions to the datablock.
    status = block->put_val( inflation::fail_names, "no_end_inflation",   inflation::no_end_inflate);
    status = block->put_val( inflation::fail_names, "negative_Hsq",       inflation::neg_Hsq);
    status = block->put_val( inflation::fail_names, "integrate_nan",      inflation::integrate_nan);
    status = block->put_val( inflation::fail_names, "zero_massless_time", inflation::zero_massless);
    status = block->put_val( inflation::fail_names, "negative_epsilon",   inflation::neg_epsilon);
    status = block->put_val( inflation::fail_names, "eps_geq_three",      inflation::large_epsilon);
    status = block->put_val( inflation::fail_names, "negative_pot",       inflation::neg_V);
    status = block->put_val( inflation::fail_names, "noFind_hor_exit",    inflation::failed_horizonExit);
    status = block->put_val( inflation::fail_names, "ICs_before_start",   inflation::ics_before_start);
    status = block->put_val( inflation::fail_names, "leq_60_efolds",      inflation::inflate60);
    status = block->put_val( inflation::fail_names, "varying_Spec",       inflation::time_var_pow_spec);

    // Sum all the failed sample ints and add that to status - if any = 1 then break out of pipeline for this sample.
    int err_sum = inflation::no_end_inflate + inflation::neg_Hsq + inflation::integrate_nan + inflation::zero_massless +
                  inflation::neg_epsilon + inflation::large_epsilon + inflation::neg_V + inflation::failed_horizonExit +
                  inflation::ics_before_start + inflation::inflate60 + inflation::time_var_pow_spec;
    if (err_sum >= 1)
    {
        return failure;
    }

    //! Create a temporary path & file for passing wave-number information to the datablock for class
    boost::filesystem::path temp_path = boost::filesystem::temp_directory_path() / boost::filesystem::unique_path("%%%%-%%%%-%%%%-%%%%.dat");
    std::cout << "Temp. path = " << temp_path.string() << std::endl;
    std::ofstream outf(temp_path.string(), std::ios_base::out | std::ios_base::trunc);
    for (int i = 0; i < Phys_waveno_sample.size(); ++i) {
        outf << Phys_waveno_sample[i] << "\t";
        outf << A_s[i] << "\t";
        outf << A_t[i] << "\n";
    }
    outf.close();

    //! Return the calculated observables to the datablock
    // PIVOT TASKS
    // Use put_val method to add second-order observables (k, A_s, A_t, n_s, n_t, r) at the pivot scale to the datablock
    status = block->put_val( inflation::twopf_name, "k_piv", k_pivot_cppt );
    status = block->put_val( inflation::twopf_name, "N_piv", N_pivot_exit );
    status = block->put_val( inflation::twopf_name, "A_s", A_s_pivot );
    status = block->put_val( inflation::twopf_name, "A_t", A_t_pivot );
    status = block->put_val( inflation::twopf_name, "n_s", ns_pivot );
    status = block->put_val( inflation::twopf_name, "n_t", nt_pivot );
    status = block->put_val( inflation::twopf_name, "r", r_pivot );
    // Use put_val to put the three-point observables (B_equi, fNL_equi) onto the datablock
    status = block->put_val( inflation::thrpf_name, "B_equi", B_equi_piv );
    status = block->put_val( inflation::thrpf_name, "fNL_equi", fNL_equi_piv );
//    status = block->put_val( inflation::thrpf_name, "B_squ", B_squ_piv );
//    status = block->put_val( inflation::thrpf_name, "fNL_squ", fNL_squ_piv );

    // CMB TASK for Boltzmann solver
    // Use put_val to write the temporary file with k, P_s(k) and P_t(k) information for CLASS
    status = block->put_val( inflation::spec_file, "spec_table", temp_path.string() );

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
