#include <vector>
#include <list>
#include <functional>
#include <memory>

#include "boost/filesystem/operations.hpp"
#include "boost/timer/timer.hpp"
#include "boost/log/core.hpp"
#include "boost/log/trivial.hpp"
#include "boost/log/expressions.hpp"
#include "boost/log/attributes.hpp"
#include "boost/log/sources/severity_feature.hpp"
#include "boost/log/sources/severity_logger.hpp"
#include "boost/log/sinks/sync_frontend.hpp"
#include "boost/log/sinks/text_file_backend.hpp"
#include "boost/log/utility/setup/common_attributes.hpp"
#include "boost/log/attributes.hpp"
#include "boost/log/expressions.hpp"
#include "boost/log/support/date_time.hpp"
#include "boost/date_time/posix_time/posix_time.hpp"

#include "transport-runtime/enumerations.h"
#include "transport-runtime/data/batchers/integration_items.h"

#include "transport-runtime/defaults.h"
#include "transport-runtime/data/batchers/postprocess_delegate.h"
#include "transport-runtime/tasks/task_configurations.h"

class sampling_integration_batcher
  {

  // Enum definitions needed for logging
  public:
    //! Logging severity level
    enum class log_severity_level { datapipe_pull, normal, warning, error, critical };

    //! logging sink
    typedef boost::log::sinks::synchronous_sink< boost::log::sinks::text_file_backend > sink_t;
    
    //! logging source
    typedef boost::log::sources::severity_logger<log_severity_level> logger;

  // CONSTRUCTOR, DESTRUCTOR
  public:

    //! constructor
    sampling_integration_batcher(const boost::filesystem::path& lp, unsigned int w,
                                transport::model<double>* m, transport::integration_task<double>* tk,
                                unsigned int g=0, bool no_log=false);

    //! move constructor
    sampling_integration_batcher(sampling_integration_batcher&&) = default;

    ~sampling_integration_batcher() = default;

  // BATCHER INFORMATION - default batcher functions needed in CppT interface
  public:
    bool is_collecting_initial_conditions() {return false;}

    void push_ics(unsigned int k_serial, double t_exit, const std::vector<double>& values)
    {
      return;
    }

  void report_integration_success(boost::timer::nanosecond_type integration, boost::timer::nanosecond_type batching,
                                                         unsigned int kserial, size_t steps, unsigned int refinement)
    {
      return;
    }

  // LOGGING FUNCTION
  public:
    //! Return logger
    logger& get_log() { return(this->log_source); }


  protected:
    //! cache number of fields associated with this integration
//    const unsigned int Nfields;

    //! Log directory path
    boost::filesystem::path logdir_path;

    //! Worker group associated with this batcher;
    //! usually zero unless we are doing parallel batching.
    //! Later, groups identify different integrations which have been chained together
    unsigned int worker_group;
    
    //! Worker number associated with this batcher
    unsigned int worker_number;

    // OTHER INTERNAL DATA
    // LOGGING
    //! Logger source
    boost::log::sources::severity_logger<sampling_integration_batcher::log_severity_level> log_source;
    
    //! Logger sink; note we are forced to use boost::shared_ptr<> because this
    //! is what the Boost.Log API expects
    boost::shared_ptr<sink_t> log_sink;
  };

// DECLARE TWOPF_BATCHER AND METHODS/DATA
class twopf_sampling_batcher: public sampling_integration_batcher
{
public:
  //! Constructor
  twopf_sampling_batcher(std::vector<double>& dt_twopf, std::vector<double>& dt_tenspf,
                         const boost::filesystem::path& lp, unsigned int w,
                         transport::model<double>* m, transport::twopf_task<double>* tk,
                         unsigned int g=0, bool no_log=false);

  //! move constructor
  twopf_sampling_batcher(twopf_sampling_batcher&&) = default;

  //! destructor is default
  virtual ~twopf_sampling_batcher() = default;

  //! BATCHERS
  void push_backg(unsigned int time_serial, unsigned int source_serial, const std::vector<double>& values);

  void push_twopf(unsigned int time_serial, unsigned int k_serial, unsigned int source_serial,
                  const std::vector<double>& twopf_values, const std::vector<double>& backg);

  void push_tensor_twopf(unsigned int time_serial, unsigned int k_serial,
                       unsigned int source_serial, const std::vector<double>& tensor_values);

// INTERNAL DATA
private:
  //! std::vector for collecting the zeta-twopf samples in
  std::vector<double>& zeta_twopf_data;

  //! std::vector for collecting the tensor-twopf samples in
  std::vector<double>& tensor_twopf_data;

protected:
  //! cache number of fields associated with this integration
  const unsigned int Nfields;

  // COMPUTATION AGENT
  // compute delegate
  transport::postprocess_delegate<double> compute_agent;

  //! cache for linear part of gauge transformation
  std::vector<double> gauge_xfm1;
};

// DECLARE THREEPF_BATCHER AND METHODS/DATA
class threepf_sampling_batcher: public sampling_integration_batcher
{
public:
  //! Constructor
  threepf_sampling_batcher(std::vector<double>& dt_twopf, std::vector<double>& dt_tenspf,
          std::vector<double>& dt_thrpf, std::vector<double>& dt_redbsp,
          const boost::filesystem::path& lp, unsigned int w,transport::model<double>* m,
          transport::threepf_task<double>* tk, unsigned int g=0, bool no_log=false);

  //! move constructor
  threepf_sampling_batcher(threepf_sampling_batcher&&) = default;

  //! destructor is default
  virtual ~threepf_sampling_batcher() = default;

  //! BATCHERS
  // Push the background
  void push_backg(unsigned int time_serial, unsigned int source_serial, const std::vector<double>& values);

  // Push a set of initial condition
  void push_kt_ics(unsigned int k_serial, double t_exit, const std::vector<double>& values)
  {
    return;
  }

  // Functions giving batches needed for different n-point functions
  void push_twopf(unsigned int time_serial, unsigned int k_serial, unsigned int source_serial,
                  const std::vector<double>& twopf_values, const std::vector<double>& backg, transport::twopf_type t);

  void push_tensor_twopf(unsigned int time_serial, unsigned int k_serial,
                         unsigned int source_serial, const std::vector<double>& tensor_values);

  void push_threepf(unsigned int time_serial, double t,
                    const transport::threepf_kconfig& kconfig, unsigned int source_serial,
                    const std::vector<double>& threepf,
                    const std::vector<double>& tpf_k1_re, const std::vector<double>& tpf_k1_im,
                    const std::vector<double>& tpf_k2_re, const std::vector<double>& tpf_k2_im,
                    const std::vector<double>& tpf_k3_re, const std::vector<double>& tpf_k3_im, const std::vector<double>& bg);

  // INTERNAL DATA
private:
  //! std::vector for collecting the zeta-twopf samples in
  std::vector<double>& zeta_twopf_data;

  //! std::vector for collecting the tensor-twopf samples in
  std::vector<double>& tensor_twopf_data;

  //! std::vector for collecting the zeta-threepf samples in
  std::vector<double>& zeta_threepf_data;

  //! std::vector for collecting the reduced bispectrum samples in
  std::vector<double>& redbsp_data;

protected:
  //! cache number of fields associated with this integration
  const unsigned int Nfields;

  // COMPUTATION AGENT
  // compute delegate
  transport::postprocess_delegate<double> compute_agent;

  //! cache for linear part of gauge transformation
  std::vector<double> gauge_xfm1;

  //! cache for quadratic part of gauge transformation, 123 permutation
  std::vector<double> gauge_xfm2_123;

  //! cache for quadratic part of gauge transformation, 213 permutation
  std::vector<double> gauge_xfm2_213;

  //! cache for quadratic part of gauge transformation, 312 permutation
  std::vector<double> gauge_xfm2_312;
};

// overload << to push log_severity_level to stream
std::ostream& operator<<(std::ostream& stream, sampling_integration_batcher::log_severity_level level)
  {
    static const std::map< sampling_integration_batcher::log_severity_level, std::string > stringize_map =
      {
        { sampling_integration_batcher::log_severity_level::datapipe_pull, "datapipe" },
        { sampling_integration_batcher::log_severity_level::normal, "normal" },
        { sampling_integration_batcher::log_severity_level::warning, "warning" },
        { sampling_integration_batcher::log_severity_level::error, "error" },
        { sampling_integration_batcher::log_severity_level::critical, "CRITICAL" }
      };
    
    stream << stringize_map.at(level);
    
    return stream;
  }

// overload << to push log_severity_level to Boost.Log log
struct sampling_batcher_severity_tag;
boost::log::formatting_ostream& operator<<(boost::log::formatting_ostream& stream,
                                         const boost::log::to_log_manip<sampling_integration_batcher::log_severity_level, sampling_batcher_severity_tag> manip)
  {
    static const std::map< sampling_integration_batcher::log_severity_level, std::string > stringize_map =
      {
        { sampling_integration_batcher::log_severity_level::datapipe_pull, "data" },
        { sampling_integration_batcher::log_severity_level::normal, "info" },
        { sampling_integration_batcher::log_severity_level::warning, "warn" },
        { sampling_integration_batcher::log_severity_level::error, "err " },
        { sampling_integration_batcher::log_severity_level::critical, "CRIT" }
      };
    
    sampling_integration_batcher::log_severity_level level = manip.get();
    stream << stringize_map.at(level);
    
    return stream;
  }

// GENERIC BATCHER METHODS
// CONSTRUCTOR, DESTRUCTOR
sampling_integration_batcher::sampling_integration_batcher(const boost::filesystem::path& lp,
                            unsigned int w, transport::model<double>* m, transport::integration_task<double>* tk, unsigned int g, bool no_log)
  : logdir_path(lp),
    worker_group(g),
    worker_number(w)
{
  // set up logging
  std::ostringstream log_file;
  log_file << transport::CPPTRANSPORT_LOG_FILENAME_A << worker_number << transport::CPPTRANSPORT_LOG_FILENAME_B;

  boost::filesystem::path log_path = logdir_path / log_file.str();

  if(!no_log)
    {
      boost::shared_ptr<boost::log::core> core = boost::log::core::get();

  //        core->set_filter(boost::log::trivial::severity >= normal);

      boost::shared_ptr<boost::log::sinks::text_file_backend> backend =
        boost::make_shared<boost::log::sinks::text_file_backend>(
          boost::log::keywords::file_name = log_path.string(),
          boost::log::keywords::open_mode = std::ios::app
        );

      // enable auto-flushing of log entries
      // this degrades performance, but we are not writing many entries and they
      // will not be lost in the event of a crash
      backend->auto_flush(true);

      // Wrap it into the frontend and register in the core.
      // The backend requires synchronization in the frontend.
      this->log_sink = boost::make_shared<sink_t>(backend);
      this->log_sink->set_formatter(
        boost::log::expressions::stream
          << boost::log::expressions::format_date_time<boost::posix_time::ptime>("TimeStamp", "%Y-%m-%d %H:%M:%S")
          << " | "
          << boost::log::expressions::attr< sampling_integration_batcher::log_severity_level, sampling_batcher_severity_tag >("Severity")
          << " | "
          << boost::log::expressions::smessage
      );

      core->add_sink(this->log_sink);

      boost::log::add_common_attributes();
    }

}

//! TWOPF_SAMPLING METHODS
// CONSTRUCTOR
twopf_sampling_batcher::twopf_sampling_batcher(std::vector<double> &dt_twopf, std::vector<double> &dt_tenspf,
                                               const boost::filesystem::path &lp,
                                               unsigned int w, transport::model<double> *m,
                                               transport::twopf_task<double> *tk, unsigned int g, bool no_log)
     : sampling_integration_batcher(lp, w, m, tk, g, no_log),
       zeta_twopf_data(dt_twopf),
       tensor_twopf_data(dt_tenspf),
       Nfields(m->get_N_fields()),
       compute_agent(m,tk)
{
  // Ensure the std::vector for the 1st-order gauge transform has the correct dimensions
  gauge_xfm1.resize(2*this->Nfields);
}

// PUSH_BACKG
void twopf_sampling_batcher::push_backg(unsigned int time_serial, unsigned int source_serial,
                                        const std::vector<double> &values)
{
  return;
}

// PUSH_TWOPF
void twopf_sampling_batcher::push_twopf(unsigned int time_serial, unsigned int k_serial, unsigned int source_serial,
                                        const std::vector<double> &twopf_values, const std::vector<double> &backg)
{
  double zeta_twopf = 0.0;
  this->compute_agent.zeta_twopf(twopf_values, backg, zeta_twopf, this->gauge_xfm1);
  this->zeta_twopf_data.push_back(zeta_twopf);
  return;
}

// PUSH_TENSOR_TWOPF
void twopf_sampling_batcher::push_tensor_twopf(unsigned int time_serial, unsigned int k_serial,
                                               unsigned int source_serial, const std::vector<double> &tensor_values)
{
  double tensor_twopf = tensor_values[0];
  this->tensor_twopf_data.emplace_back(tensor_twopf);
}

//! THREEPF_SAMPLING METHODS
// CONSTRUCTOR
threepf_sampling_batcher::threepf_sampling_batcher(std::vector<double> &dt_twopf, std::vector<double> &dt_tenspf,
                                                   std::vector<double> &dt_thrpf, std::vector<double> &dt_redbsp,
                                                   const boost::filesystem::path &lp,
                                                   unsigned int w, transport::model<double> *m,
                                                   transport::threepf_task<double> *tk, unsigned int g, bool no_log)
       :  sampling_integration_batcher(lp, w, m, tk, g, no_log),
          zeta_twopf_data(dt_twopf),
          tensor_twopf_data(dt_tenspf),
          zeta_threepf_data(dt_thrpf),
          redbsp_data(dt_redbsp),
          Nfields(m->get_N_fields()),
          compute_agent(m,tk)
{
  // Ensure the std::vector for the 1st-order gauge transform has the correct dimensions
  gauge_xfm1.resize(2*this->Nfields);

  // Ensure the std::vectors for the quadratic gauge transformations have the correct dimensions
  gauge_xfm2_123.resize(2*this->Nfields * 2*this->Nfields);
  gauge_xfm2_213.resize(2*this->Nfields * 2*this->Nfields);
  gauge_xfm2_312.resize(2*this->Nfields * 2*this->Nfields);
}

// PUSH_BACKG
void threepf_sampling_batcher::push_backg(unsigned int time_serial, unsigned int source_serial,
                                          const std::vector<double> &values)
{
  return;
}

// PUSH_TWOPF
void threepf_sampling_batcher::push_twopf(unsigned int time_serial, unsigned int k_serial, unsigned int source_serial,
                                          const std::vector<double> &twopf_values, const std::vector<double> &backg,
                                          transport::twopf_type t)
{
  double zeta_twopf = 0.0;
  this->compute_agent.zeta_twopf(twopf_values, backg, zeta_twopf, this->gauge_xfm1);
  this->zeta_twopf_data.push_back(zeta_twopf);
  return;
}

// PUSH_TENSOR_TWOPF
void threepf_sampling_batcher::push_tensor_twopf(unsigned int time_serial, unsigned int k_serial,
                                                 unsigned int source_serial,
                                                 const std::vector<double> &tensor_values)
{
  double tensor_twopf = tensor_values[0];
  this->tensor_twopf_data.emplace_back(tensor_twopf);
}

// PUSH_THREEPF
void threepf_sampling_batcher::push_threepf(unsigned int time_serial, double t, const transport::threepf_kconfig &kconfig,
                                            unsigned int source_serial, const std::vector<double> &threepf,
                                            const std::vector<double> &tpf_k1_re, const std::vector<double> &tpf_k1_im,
                                            const std::vector<double> &tpf_k2_re, const std::vector<double> &tpf_k2_im,
                                            const std::vector<double> &tpf_k3_re, const std::vector<double> &tpf_k3_im,
                                            const std::vector<double> &bg)
{
  double zeta_threepf = 0.0;
  double redbsp = 0.0;

  this->compute_agent.zeta_threepf(kconfig, t, threepf, tpf_k1_re, tpf_k1_im, tpf_k2_re, tpf_k2_im, tpf_k3_re, tpf_k3_im, bg, zeta_threepf, redbsp,
                                   this->gauge_xfm1, this->gauge_xfm2_123, this->gauge_xfm2_213, this->gauge_xfm2_312);

  this->zeta_threepf_data.push_back(zeta_threepf);
  this->redbsp_data.push_back(redbsp);

  return;
}