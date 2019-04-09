//
// --@@
// Copyright (c) 2016 University of Sussex. All rights reserved.
//
// This template file is part of the CppTransport platform.
//
// CppTransport is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// CppTransport is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with CppTransport.  If not, see <http://www.gnu.org/licenses/>.
//
// As a special exception, you may create a larger work that contains
// part or all of this template file and distribute that work
// under terms of your choice.  Alternatively, if you modify or redistribute
// this template file itself, you may (at your option) remove this
// special exception, which will cause the template and the resulting
// CppTransport output files to be licensed under the GNU General Public
// License without this special exception.
//
// @license: GPL-2
// @contributor: David Seery <D.Seery@sussex.ac.uk>
// --@@
//
// DO NOT EDIT: GENERATED AUTOMATICALLY BY CppTransport 2018.1
//
// 'quartic_mpi.h' generated from '/home/seanb/cosmosis/modules/cppt_sample/cppt_sample/Quartic_Inflation/quartic.model'
// processed on 2019-Apr-09 14:34:11

// MPI implementation

#ifndef CPPTRANSPORT_QUARTIC_MPI_H   // avoid multiple inclusion
#define CPPTRANSPORT_QUARTIC_MPI_H

#include "transport-runtime/transport.h"

#include "quartic_core.h"

namespace transport
  {

//     phase-space flattener set to 'FLATTEN' 
//     field-space flattener set to 'FIELDS_FLATTEN' 

//     working type set to 'number' 

//     IF fast 

//       set replacement rule 'U2_DECLARE' to "const auto __u2_$M_$N" 

//       set replacement rule 'U2_k1_DECLARE' to "const auto __u2_k1_$M_$N" 
//       set replacement rule 'U2_k2_DECLARE' to "const auto __u2_k2_$M_$N" 
//       set replacement rule 'U2_k3_DECLARE' to "const auto __u2_k3_$M_$N" 

//       set replacement rule 'U3_k1k2k3_DECLARE' to "const auto __u3_k1k2k3_$L_$M_$N" 
//       set replacement rule 'U3_k2k1k3_DECLARE' to "const auto __u3_k2k1k3_$L_$M_$N" 
//       set replacement rule 'U3_k3k1k2_DECLARE' to "const auto __u3_k3k1k2_$L_$M_$N" 

//       set replacement rule 'U2_CONTAINER' to "__u2_$M_$N" 

//       set replacement rule 'U2_k1_CONTAINER' to "__u2_k1_$M_$N" 
//       set replacement rule 'U2_k2_CONTAINER' to "__u2_k2_$M_$N" 
//       set replacement rule 'U2_k3_CONTAINER' to "__u2_k3_$M_$N" 

//       set replacement rule 'U3_k1k2k3_CONTAINER' to "__u3_k1k2k3_$L_$M_$N" 
//       set replacement rule 'U3_k2k1k3_CONTAINER' to "__u3_k2k1k3_$L_$M_$N" 
//       set replacement rule 'U3_k3k1k2_CONTAINER' to "__u3_k3k1k2_$L_$M_$N" 

//     ENDIF fast 

    namespace quartic_pool
      {
        const static std::string backend = "MPI";
        const static std::string pert_stepper = "runge_kutta_dopri5";
        const static std::string back_stepper = "runge_kutta_dopri5";
      }


    // *********************************************************************************************


    // CLASS FOR quartic '*_mpi', ie., an MPI-based implementation
    template <typename number = default_number_type, typename StateType = std::vector<number> >
    class quartic_mpi : public quartic<number>
      {

        // TYPES

      public:

        //! expose floating point value type
        using value_type = number;

        //! expose 2pf/3pf integration state type
        using twopf_state = StateType;
        using threepf_state = StateType;


        // CONSTRUCTOR, DESTRUCTOR

      public:

        //! constructor
        quartic_mpi(local_environment& e, argument_cache& a)
          : quartic<number>(e, a)
          {
#ifdef CPPTRANSPORT_INSTRUMENT
            twopf_setup_timer.stop();
            twopf_u_tensor_timer.stop();
            twopf_transport_eq_timer.stop();

            threepf_setup_timer.stop();
            threepf_u_tensor_timer.stop();
            threepf_transport_eq_timer.stop();

            twopf_items = 0;
            threepf_items = 0;

            twopf_invokations = 0;
            threepf_invokations = 0;
#endif
          }

        // destructor is default unless instrumented
#ifdef CPPTRANSPORT_INSTRUMENT
        //! instrumented destructor
        ~quartic_mpi()
          {
            if(this->twopf_items > 0)
              {
                std::cout << '\n' << "TWOPF INSTRUMENTATION REPORT" << '\n';
                std::cout << "* TOTALS" << '\n';
                std::cout << "  -- setup: " << format_time(twopf_setup_timer.elapsed().user) << " user, " << format_time(twopf_setup_timer.elapsed().system) << " system, " << format_time(twopf_setup_timer.elapsed().wall) << " wall" << '\n';
                std::cout << "  -- U tensors: " << format_time(twopf_u_tensor_timer.elapsed().user) << " user, " << format_time(twopf_u_tensor_timer.elapsed().system) << " system, " << format_time(twopf_u_tensor_timer.elapsed().wall) << " wall" << '\n';
                std::cout << "  -- transport equations: " << format_time(twopf_transport_eq_timer.elapsed().user) << " user, " << format_time(twopf_transport_eq_timer.elapsed().system) << " system, " << format_time(twopf_transport_eq_timer.elapsed().wall) << " wall" << '\n';
                std::cout << "* PER ITEM" << '\n';
                std::cout << "  -- setup: " << format_time(twopf_setup_timer.elapsed().user/this->twopf_items) << " user, " << format_time(twopf_setup_timer.elapsed().system/this->twopf_items) << " system, " << format_time(twopf_setup_timer.elapsed().wall/this->twopf_items) << " wall" << '\n';
                std::cout << "  -- U tensors: " << format_time(twopf_u_tensor_timer.elapsed().user/this->twopf_items) << " user, " << format_time(twopf_u_tensor_timer.elapsed().system/this->twopf_items) << " system, " << format_time(twopf_u_tensor_timer.elapsed().wall/this->twopf_items) << " wall" << '\n';
                std::cout << "  -- transport equations: " << format_time(twopf_transport_eq_timer.elapsed().user/this->twopf_items) << " user, " << format_time(twopf_transport_eq_timer.elapsed().system/this->twopf_items) << " system, " << format_time(twopf_transport_eq_timer.elapsed().wall/this->twopf_items) << " wall" << '\n';
                std::cout << "* PER INVOKATION" << '\n';
                std::cout << "  -- setup: " << format_time(twopf_setup_timer.elapsed().user/this->twopf_invokations) << " user, " << format_time(twopf_setup_timer.elapsed().system/this->twopf_invokations) << " system, " << format_time(twopf_setup_timer.elapsed().wall/this->twopf_invokations) << " wall" << '\n';
                std::cout << "  -- U tensors: " << format_time(twopf_u_tensor_timer.elapsed().user/this->twopf_invokations) << " user, " << format_time(twopf_u_tensor_timer.elapsed().system/this->twopf_invokations) << " system, " << format_time(twopf_u_tensor_timer.elapsed().wall/this->twopf_invokations) << " wall" << '\n';
                std::cout << "  -- transport equations: " << format_time(twopf_transport_eq_timer.elapsed().user/this->twopf_invokations) << " user, " << format_time(twopf_transport_eq_timer.elapsed().system/this->twopf_invokations) << " system, " << format_time(twopf_transport_eq_timer.elapsed().wall/this->twopf_invokations) << " wall" << '\n';
              }

            if(this->threepf_items > 0)
              {
                std::cout << '\n' << "THREEPF INSTRUMENTATION REPORT" << '\n';
                std::cout << "* TOTALS" << '\n';
                std::cout << "  -- setup: " << format_time(threepf_setup_timer.elapsed().user) << " user, " << format_time(threepf_setup_timer.elapsed().system) << " system, " << format_time(threepf_setup_timer.elapsed().wall) << " wall" << '\n';
                std::cout << "  -- U tensors: " << format_time(threepf_u_tensor_timer.elapsed().user) << " user, " << format_time(threepf_u_tensor_timer.elapsed().system) << " system, " << format_time(threepf_u_tensor_timer.elapsed().wall) << " wall" << '\n';
                std::cout << "  -- transport equations: " << format_time(threepf_transport_eq_timer.elapsed().user) << " user, " << format_time(threepf_transport_eq_timer.elapsed().system) << " system, " << format_time(threepf_transport_eq_timer.elapsed().wall) << " wall" << '\n';
                std::cout << "* PER ITEM" << '\n';
                std::cout << "  -- setup: " << format_time(threepf_setup_timer.elapsed().user/this->threepf_items) << " user, " << format_time(threepf_setup_timer.elapsed().system/this->threepf_items) << " system, " << format_time(threepf_setup_timer.elapsed().wall/this->threepf_items) << " wall" << '\n';
                std::cout << "  -- U tensors: " << format_time(threepf_u_tensor_timer.elapsed().user/this->threepf_items) << " user, " << format_time(threepf_u_tensor_timer.elapsed().system/this->threepf_items) << " system, " << format_time(threepf_u_tensor_timer.elapsed().wall/this->threepf_items) << " wall" << '\n';
                std::cout << "  -- transport equations: " << format_time(threepf_transport_eq_timer.elapsed().user/this->threepf_items) << " user, " << format_time(threepf_transport_eq_timer.elapsed().system/this->threepf_items) << " system, " << format_time(threepf_transport_eq_timer.elapsed().wall/this->threepf_items) << " wall" << '\n';
                std::cout << "* PER INVOKATION" << '\n';
                std::cout << "  -- setup: " << format_time(threepf_setup_timer.elapsed().user/this->threepf_invokations) << " user, " << format_time(threepf_setup_timer.elapsed().system/this->threepf_invokations) << " system, " << format_time(threepf_setup_timer.elapsed().wall/this->threepf_invokations) << " wall" << '\n';
                std::cout << "  -- U tensors: " << format_time(threepf_u_tensor_timer.elapsed().user/this->threepf_invokations) << " user, " << format_time(threepf_u_tensor_timer.elapsed().system/this->threepf_invokations) << " system, " << format_time(threepf_u_tensor_timer.elapsed().wall/this->threepf_invokations) << " wall" << '\n';
                std::cout << "  -- transport equations: " << format_time(threepf_transport_eq_timer.elapsed().user/this->threepf_invokations) << " user, " << format_time(threepf_transport_eq_timer.elapsed().system/this->threepf_invokations) << " system, " << format_time(threepf_transport_eq_timer.elapsed().wall/this->threepf_invokations) << " wall" << '\n';
              }
          }
#else
        //! uninstrumented destructor
        virtual ~quartic_mpi() = default;
#endif


        // EXTRACT MODEL INFORMATION -- implements a 'model' interface

      public:

        //! return backend name
        const std::string& get_backend() const override { return(quartic_pool::backend); }

        //! return background stepper name
        const std::string& get_back_stepper() const override { return(quartic_pool::back_stepper); }

        //! return perturbations stepper name
        const std::string& get_pert_stepper() const override { return(quartic_pool::pert_stepper); }

        //! return background tolerances
        std::pair< double, double > get_back_tol() const override { return std::make_pair(9.9999999999999998e-13, 9.9999999999999998e-13); }

        //! return perturbations tolerances
        std::pair< double, double > get_pert_tol() const override { return std::make_pair(1e-14, 1e-14); }


        // BACKEND INTERFACE

      public:

        //! set up a context
        context backend_get_context() override;

        //! report backend type
        worker_type get_backend_type() override;

        //! report backend memory capacity
        unsigned int get_backend_memory() override;

        //! report backend priority
        unsigned int get_backend_priority() override;

        //! integrate background and 2-point function on the CPU
        void backend_process_queue(work_queue<twopf_kconfig_record>& work, const twopf_db_task<number>* tk,
                                   twopf_batcher<number>& batcher, bool silent = false) override;

        //! integrate background, 2-point function and 3-point function on the CPU
        void backend_process_queue(work_queue<threepf_kconfig_record>& work, const threepf_task<number>* tk,
                                   threepf_batcher<number>& batcher, bool silent = false) override;

        //! report 2pf integrator state size
        unsigned int backend_twopf_state_size() const override { return(quartic_pool::twopf_state_size); }

        //! report 3pf integrator state size
        unsigned int backend_threepf_state_size() const override { return(quartic_pool::threepf_state_size); }

        //! query whether backend support collection of per-configuration statistics
        virtual bool supports_per_configuration_statistics() const override { return(true); }


        // INTERNAL API

      public:

        //! integrate a single 2pf k-configuration
        template <typename BatchObject>
        void twopf_kmode(const twopf_kconfig_record& kconfig, const twopf_db_task<number>* tk,
                         BatchObject& batcher, unsigned int refinement_level);

        //! integrate a single 3pf k-configuration
        template <typename BatchObject>
        void threepf_kmode(const threepf_kconfig_record&, const threepf_task<number>* tk,
                           BatchObject& batcher, unsigned int refinement_level);

      protected:

        //! populate initial values for a 2pf configuration
        void populate_twopf_ic(twopf_state& x, unsigned int start, double kmode, double Ninit,
                               const twopf_db_task<number>* tk, const std::vector<number>& ic, double k_normalize=1.0, bool imaginary = false);

        //! populate initial values for a tensor 2pf configuration
        void populate_tensor_ic(twopf_state& x, unsigned int start, double kmode, double Ninit,
                                const twopf_db_task<number>* tk, const std::vector<number>& ic, double k_normalize=1.0);

        //! populate initial values for a 3pf configuration
        void populate_threepf_ic(threepf_state& x, unsigned int start, const threepf_kconfig& kconfig,
                                 double Ninit, const twopf_db_task<number>* tk, const std::vector<number>& ic, double k_normalize=1.0);


        // INTERNAL DATA

      private:

#ifdef CPPTRANSPORT_INSTRUMENT
        boost::timer::cpu_timer twopf_setup_timer;
        boost::timer::cpu_timer twopf_u_tensor_timer;
        boost::timer::cpu_timer twopf_transport_eq_timer;

        unsigned int twopf_items;
        unsigned int twopf_invokations;

        boost::timer::cpu_timer threepf_setup_timer;
        boost::timer::cpu_timer threepf_u_tensor_timer;
        boost::timer::cpu_timer threepf_transport_eq_timer;

        unsigned int threepf_items;
        unsigned int threepf_invokations;
#endif

      };


    // integration - 2pf functor
    template <typename Model>
    class quartic_mpi_twopf_functor
      {

      public:

        //! inherit number type from Model
        using number = typename Model::value_type;

        //! inherit state type from model
        using twopf_state = typename Model::twopf_state;


      public:

        quartic_mpi_twopf_functor(const twopf_db_task<number>* tk, const twopf_kconfig& k
#ifdef CPPTRANSPORT_INSTRUMENT
          ,
            boost::timer::cpu_timer& st,
            boost::timer::cpu_timer& ut,
            boost::timer::cpu_timer& tt,
            unsigned int& in
#endif
        )
          : __params(tk->get_params()),
            __Mp(tk->get_params().get_Mp()),
            __N_horizon_exit(tk->get_N_horizon_crossing()),
            __astar_normalization(tk->get_astar_normalization()),
            __config(k),
            __k(k.k_comoving),

//             ENDIF !fast 

            __raw_params(nullptr)
#ifdef CPPTRANSPORT_INSTRUMENT
            ,
                __setup_timer(st),
                __u_tensor_timer(ut),
                __transport_eq_timer(tt),
                __invokations(in)
#endif
          {
          }

        void set_up_workspace()
          {
//             ENDIF !fast 

            this->__raw_params = new number[1];

            const auto& __pvector = __params.get_vector();
            this->__raw_params[0] = __pvector[0];
          }

        void close_down_workspace()
          {
//             ENDIF !fast 

            delete[] this->__raw_params;
          }

        void operator()(const twopf_state& __x, twopf_state& __dxdt, number __t);

        // adjust horizon exit time, given an initial time N_init which we wish to move to zero
        void rebase_horizon_exit_time(double N_init) { this->__N_horizon_exit -= N_init; }


        // INTERNAL DATA

      private:

        const parameters<number>& __params;

        number __Mp;

        double __N_horizon_exit;

        double __astar_normalization;

        const twopf_kconfig __config;

        const double __k;

        // manage memory ourselves, rather than via an STL container, for maximum performance
        // also avoids copying overheads (the Boost odeint library copies the functor by value)

//         ENDIF !fast 

        number* __raw_params;

#ifdef CPPTRANSPORT_INSTRUMENT
        boost::timer::cpu_timer& __setup_timer;
        boost::timer::cpu_timer& __u_tensor_timer;
        boost::timer::cpu_timer& __transport_eq_timer;
        unsigned int& __invokations;
#endif

      };


    // integration - observer object for 2pf
    template <typename Model, typename BatchObject>
    class quartic_mpi_twopf_observer: public twopf_singleconfig_batch_observer<typename Model::value_type, BatchObject>
      {

      public:

        //! inherit number type from Model
        using number = typename Model::value_type;

        //! inherit state type from model
        using twopf_state = typename Model::twopf_state;


      public:

        quartic_mpi_twopf_observer(BatchObject& b, const twopf_kconfig_record& c,
                                  double t_ics, const time_config_database& t)
          : twopf_singleconfig_batch_observer<number, BatchObject>
              (b, c, t_ics, t,
               quartic_pool::backg_size, quartic_pool::tensor_size, quartic_pool::twopf_size,
               quartic_pool::backg_start, quartic_pool::tensor_start, quartic_pool::twopf_start)
          {
          }

        void operator()(const twopf_state& x, number t);

      };


    // integration - 3pf functor
    template <typename Model>
    class quartic_mpi_threepf_functor
      {

      public:

        //! inherit number type from Model
        using number = typename Model::value_type;

        //! inherit state type from model
        using threepf_state = typename Model::threepf_state;


      public:

        quartic_mpi_threepf_functor(const twopf_db_task<number>* tk, const threepf_kconfig& k
#ifdef CPPTRANSPORT_INSTRUMENT
          ,
          boost::timer::cpu_timer& st,
          boost::timer::cpu_timer& ut,
          boost::timer::cpu_timer& tt,
          unsigned int& in
#endif
        )
          : __params(tk->get_params()),
            __Mp(tk->get_params().get_Mp()),
            __N_horizon_exit(tk->get_N_horizon_crossing()),
            __astar_normalization(tk->get_astar_normalization()),
            __config(k),
            __k1(k.k1_comoving),
            __k2(k.k2_comoving),
            __k3(k.k3_comoving),

//             ENDIF !fast 

            __raw_params(nullptr)
#ifdef CPPTRANSPORT_INSTRUMENT
            ,
            __setup_timer(st),
            __u_tensor_timer(ut),
            __transport_eq_timer(tt),
            __invokations(in)
#endif
          {
          }

        void set_up_workspace()
          {
//             ENDIF !fast 

            this->__raw_params = new number[1];

            const auto& __pvector = __params.get_vector();
            this->__raw_params[0] = __pvector[0];
          }

        void close_down_workspace()
          {
//             ENDIF !fast 

            delete[] this->__raw_params;
          }

        void operator()(const threepf_state& __x, threepf_state& __dxdt, number __dt);

        // adjust horizon exit time, given an initial time N_init which we wish to move to zero
        void rebase_horizon_exit_time(double N_init) { this->__N_horizon_exit -= N_init; }

      private:

        const parameters<number>& __params;

        number __Mp;

        double __N_horizon_exit;

        double __astar_normalization;

        const threepf_kconfig __config;

        const double __k1;
        const double __k2;
        const double __k3;

        // manage memory ourselves, rather than via an STL container, for maximum performance
        // also avoids copying overheads (the Boost odeint library copies the functor by value)

//         ENDIF !fast 

        number* __raw_params;

#ifdef CPPTRANSPORT_INSTRUMENT
        boost::timer::cpu_timer& __setup_timer;
        boost::timer::cpu_timer& __u_tensor_timer;
        boost::timer::cpu_timer& __transport_eq_timer;
        unsigned int& __invokations;
#endif

      };


    // integration - observer object for 3pf
    template <typename Model, typename BatchObject>
    class quartic_mpi_threepf_observer: public threepf_singleconfig_batch_observer<typename Model::value_type, BatchObject>
      {

      public:

        //! inherit number type from Model
        using number = typename Model::value_type;

        //! inherit state type from model
        using threepf_state = typename Model::threepf_state;


      public:
        quartic_mpi_threepf_observer(BatchObject& b, const threepf_kconfig_record& c,
                                    double t_ics, const time_config_database& t)
          : threepf_singleconfig_batch_observer<number, BatchObject>
              (b, c, t_ics, t,
               quartic_pool::backg_size, quartic_pool::tensor_size,
               quartic_pool::twopf_size, quartic_pool::threepf_size,
               quartic_pool::backg_start,
               quartic_pool::tensor_k1_start, quartic_pool::tensor_k2_start, quartic_pool::tensor_k3_start,
               quartic_pool::twopf_re_k1_start, quartic_pool::twopf_im_k1_start,
               quartic_pool::twopf_re_k2_start, quartic_pool::twopf_im_k2_start,
               quartic_pool::twopf_re_k3_start, quartic_pool::twopf_im_k3_start,
               quartic_pool::threepf_start)
          {
          }

        void operator()(const threepf_state& x, number t);

      };


    // BACKEND INTERFACE


    // generate a context
    template <typename number, typename StateType>
    context quartic_mpi<number, StateType>::backend_get_context(void)
      {
        context ctx;

        // set up just one device
        ctx.add_device(quartic_pool::backend);

        return(ctx);
      }


    template <typename number, typename StateType>
    worker_type quartic_mpi<number, StateType>::get_backend_type(void)
      {
        return(worker_type::cpu);
      }


    template <typename number, typename StateType>
    unsigned int quartic_mpi<number, StateType>::get_backend_memory(void)
      {
        return(0);
      }


    template <typename number, typename StateType>
    unsigned int quartic_mpi<number, StateType>::get_backend_priority(void)
      {
        return(1);
      }


    // process work queue for twopf
    template <typename number, typename StateType>
    void quartic_mpi<number, StateType>::backend_process_queue(work_queue<twopf_kconfig_record>& work,
                                                              const twopf_db_task<number>* tk,
                                                              twopf_batcher<number>& batcher, bool silent)
      {
        // set batcher to delayed flushing mode so that we have a chance to unwind failed integrations
        batcher.set_flush_mode(generic_batcher::flush_mode::flush_delayed);

        std::ostringstream work_msg;
        BOOST_LOG_SEV(batcher.get_log(), generic_batcher::log_severity_level::normal)
            << "** MPI compute backend processing twopf task";
        work_msg << work;
        BOOST_LOG_SEV(batcher.get_log(), generic_batcher::log_severity_level::normal) << work_msg.str();
//        std::cerr << work_msg.str();
        if(!silent) this->write_task_data(tk, batcher, 1e-14, 1e-14, 1.0000000000000001e-15, "runge_kutta_dopri5");

        // get work queue for the zeroth device (should be the only device in this backend)
        assert(work.size() == 1);
        const work_queue<twopf_kconfig_record>::device_queue queues = work[0];

        // we expect only one queue on this device
        assert(queues.size() == 1);
        const work_queue<twopf_kconfig_record>::device_work_list list = queues[0];

        for(unsigned int i = 0; i < list.size(); ++i)
          {
            bool success = false;
            unsigned int refinement_level = 0;

            while(!success)
            try
              {
                // write the time history for this k-configuration
                this->twopf_kmode(list[i], tk, batcher, refinement_level);    // logging and report of successful integration are wrapped up in the observer stop_timers() method
                success = true;
               }
            catch(std::overflow_error& xe)
              {
                // unwind any batched results before trying again with a refined mesh
                if(refinement_level == 0) batcher.report_refinement();
                batcher.unbatch(list[i]->serial);
                refinement_level++;

                BOOST_LOG_SEV(batcher.get_log(), generic_batcher::log_severity_level::warning)
                    << "** " << CPPTRANSPORT_RETRY_CONFIG << " " << list[i]->serial << " (" << i+1
                    << " " CPPTRANSPORT_OF << " " << list.size() << "), "
                    << CPPTRANSPORT_REFINEMENT_LEVEL << " = " << refinement_level
                    << " (" << CPPTRANSPORT_REFINEMENT_INTERNAL << xe.what() << ")";
              }
            catch(runtime_exception& xe)
              {
                batcher.report_integration_failure(list[i]->serial);
                batcher.unbatch(list[i]->serial);
                success = true;

                BOOST_LOG_SEV(batcher.get_log(), generic_batcher::log_severity_level::error)
                    << "!! " CPPTRANSPORT_FAILED_CONFIG << " " << list[i]->serial << " (" << i+1
                    << " " CPPTRANSPORT_OF << " " << list.size() << ") | " << list[i];
              }
          }
      }


    template <typename number, typename StateType>
    template <typename BatchObject>
    void quartic_mpi<number, StateType>::twopf_kmode(const twopf_kconfig_record& kconfig,
                                                    const twopf_db_task<number>* tk, BatchObject& batcher,
                                                    unsigned int refinement_level)
      {
        DEFINE_INDEX_TOOLS

        if(refinement_level > tk->get_max_refinements()) throw runtime_exception(exception_type::REFINEMENT_FAILURE, CPPTRANSPORT_REFINEMENT_TOO_DEEP);

        // get time configuration database
        const time_config_database time_db = tk->get_time_config_database(*kconfig);

        // set up a functor to observe the integration
        // this also starts the timers running, so we do it as early as possible
        quartic_mpi_twopf_observer< quartic_mpi<number, StateType>, BatchObject > obs(batcher, kconfig, tk->get_initial_time(*kconfig), time_db);

        // set up a functor to evolve this system
        quartic_mpi_twopf_functor< quartic_mpi<number, StateType> > rhs(tk, *kconfig
#ifdef CPPTRANSPORT_INSTRUMENT
          ,
            this->twopf_setup_timer, this->twopf_u_tensor_timer, this->twopf_transport_eq_timer, this->twopf_invokations
#endif
          );
        rhs.set_up_workspace();

        // set up a state vector
        twopf_state x;
        x.resize(quartic_pool::twopf_state_size);

        // fix initial conditions - background
        const std::vector<number> ics = tk->get_ics_vector(*kconfig);
        x[quartic_pool::backg_start + FLATTEN(0)] = ics[0];
        x[quartic_pool::backg_start + FLATTEN(1)] = ics[1];

        if(batcher.is_collecting_initial_conditions())
          {
            const std::vector<number> ics_1 = tk->get_ics_exit_vector(*kconfig);
            double t_exit = tk->get_ics_exit_time(*kconfig);
            batcher.push_ics(kconfig->serial, t_exit, ics_1);
          }

        // observers expect all correlation functions to be dimensionless and rescaled by the same factors

        // fix initial conditions - tensors (use dimensionless correlation functions)
        this->populate_tensor_ic(x, quartic_pool::tensor_start, kconfig->k_comoving, *(time_db.value_begin()), tk, ics, kconfig->k_comoving);

        // fix initial conditions - 2pf (use dimensionless correlation functions)
        this->populate_twopf_ic(x, quartic_pool::twopf_start, kconfig->k_comoving, *(time_db.value_begin()), tk, ics, kconfig->k_comoving);

        // up to this point the calculation has been done in the user-supplied time variable.
        // However, the integrator apparently performs much better if times are measured from zero (but not yet clear why)
        // TODO: would be nice to remove this in future
        rhs.rebase_horizon_exit_time(tk->get_ics().get_N_initial());
        auto begin_iterator = time_db.value_begin(tk->get_ics().get_N_initial());
        auto end_iterator   = time_db.value_end(tk->get_ics().get_N_initial());

        using boost::numeric::odeint::integrate_times;

        auto stepper = boost::numeric::odeint::make_dense_output< boost::numeric::odeint::runge_kutta_dopri5< twopf_state, number, twopf_state, number, CPPTRANSPORT_ALGEBRA_NAME(twopf_state), CPPTRANSPORT_OPERATIONS_NAME(twopf_state) > >(1e-14, 1e-14);
        size_t steps = integrate_times(stepper, rhs, x, begin_iterator, end_iterator,
                                       static_cast<number>(1.0000000000000001e-15/pow(4.0,refinement_level)), obs);

        obs.stop_timers(steps, refinement_level);
        rhs.close_down_workspace();

#ifdef CPPTRANSPORT_INSTRUMENT
        ++this->twopf_items;
#endif
      }


    // make initial conditions for each component of the 2pf
    // x           - state vector *containing* space for the 2pf (doesn't have to be entirely the 2pf)
    // start       - starting position of twopf components within the state vector
    // kmode       - *comoving normalized* wavenumber for which we will compute the twopf
    // Ninit       - initial time
    // tk          - parent task
    // ics         - initial conditions for the background fields (or fields+momenta)
    // k_normalize - used to adjust ics to be dimensionless, or just 1.0 to get raw correlation function
    // imaginary   - whether to populate using real or imaginary components of the 2pf
    template <typename number, typename StateType>
    void quartic_mpi<number, StateType>::populate_twopf_ic(twopf_state& x, unsigned int start, double kmode,
                                                          double Ninit, const twopf_db_task<number>* tk,
                                                          const std::vector<number>& ics, double k_normalize,
                                                          bool imaginary)
      {
        DEFINE_INDEX_TOOLS

        assert(x.size() >= start);
        assert(x.size() >= start + quartic_pool::twopf_size);

        // populate components of the 2pf
        x[start + FLATTEN(0,0)] = imaginary ? this->make_twopf_im_ic(0, 0, kmode, Ninit, tk, ics, k_normalize) : this->make_twopf_re_ic(0, 0, kmode, Ninit, tk, ics, k_normalize);
        x[start + FLATTEN(0,1)] = imaginary ? this->make_twopf_im_ic(0, 1, kmode, Ninit, tk, ics, k_normalize) : this->make_twopf_re_ic(0, 1, kmode, Ninit, tk, ics, k_normalize);
        x[start + FLATTEN(1,0)] = imaginary ? this->make_twopf_im_ic(1, 0, kmode, Ninit, tk, ics, k_normalize) : this->make_twopf_re_ic(1, 0, kmode, Ninit, tk, ics, k_normalize);
        x[start + FLATTEN(1,1)] = imaginary ? this->make_twopf_im_ic(1, 1, kmode, Ninit, tk, ics, k_normalize) : this->make_twopf_re_ic(1, 1, kmode, Ninit, tk, ics, k_normalize);
      }


    // make initial conditions for the tensor twopf
    template <typename number, typename StateType>
    void quartic_mpi<number, StateType>::populate_tensor_ic(twopf_state& x, unsigned int start, double kmode,
                                                           double Ninit, const twopf_db_task<number>* tk,
                                                           const std::vector<number>& ics, double k_normalize)
      {
        DEFINE_INDEX_TOOLS

        assert(x.size() >= start);
        assert(x.size() >= start + quartic_pool::tensor_size);

        // populate components of the 2pf
        x[start + TENSOR_FLATTEN(0,0)] = this->make_twopf_tensor_ic(0, 0, kmode, Ninit, tk, ics, k_normalize);
        x[start + TENSOR_FLATTEN(0,1)] = this->make_twopf_tensor_ic(0, 1, kmode, Ninit, tk, ics, k_normalize);
        x[start + TENSOR_FLATTEN(1,0)] = this->make_twopf_tensor_ic(1, 0, kmode, Ninit, tk, ics, k_normalize);
        x[start + TENSOR_FLATTEN(1,1)] = this->make_twopf_tensor_ic(1, 1, kmode, Ninit, tk, ics, k_normalize);
      }


    // THREE-POINT FUNCTION INTEGRATION


    template <typename number, typename StateType>
    void quartic_mpi<number, StateType>::backend_process_queue(work_queue<threepf_kconfig_record>& work,
                                                              const threepf_task<number>* tk,
                                                              threepf_batcher<number>& batcher, bool silent)
      {
        // set batcher to delayed flushing mode so that we have a chance to unwind failed integrations
        batcher.set_flush_mode(generic_batcher::flush_mode::flush_delayed);

        std::ostringstream work_msg;
        BOOST_LOG_SEV(batcher.get_log(), generic_batcher::log_severity_level::normal)
          << "** MPI compute backend processing threepf task";
        work_msg << work;
        BOOST_LOG_SEV(batcher.get_log(), generic_batcher::log_severity_level::normal) << work_msg.str();
//        std::cerr << work_msg.str();
        if(!silent) this->write_task_data(tk, batcher, 1e-14, 1e-14, 1.0000000000000001e-15, "runge_kutta_dopri5");

        // get work queue for the zeroth device (should be only one device with this backend)
        assert(work.size() == 1);
        const work_queue<threepf_kconfig_record>::device_queue queues = work[0];

        // we expect only one queue on this device
        assert(queues.size() == 1);
        const work_queue<threepf_kconfig_record>::device_work_list list = queues[0];

        // step through the queue, solving for the three-point functions in each case
        for(unsigned int i = 0; i < list.size(); ++i)
          {
            bool success = false;
            unsigned int refinement_level = 0;

            while(!success)
            try
              {
                // write the time history for this k-configuration
                this->threepf_kmode(list[i], tk, batcher, refinement_level);    // logging and report of successful integration are wrapped up in the observer stop_timers() method
                success = true;
              }
            catch(std::overflow_error& xe)
              {
                // unwind any batched results before trying again with a refined mesh
                if(refinement_level == 0) batcher.report_refinement();
                batcher.unbatch(list[i]->serial);
                refinement_level++;

                BOOST_LOG_SEV(batcher.get_log(), generic_batcher::log_severity_level::warning)
                    << "** " << CPPTRANSPORT_RETRY_CONFIG << " " << list[i]->serial << " (" << i+1
                    << " " << CPPTRANSPORT_OF << " " << list.size() << "), "
                    << CPPTRANSPORT_REFINEMENT_LEVEL << " = " << refinement_level
                    << " (" << CPPTRANSPORT_REFINEMENT_INTERNAL << xe.what() << ")";
              }
            catch(runtime_exception& xe)
              {
                batcher.report_integration_failure(list[i]->serial);
                batcher.unbatch(list[i]->serial);
                success = true;

                BOOST_LOG_SEV(batcher.get_log(), generic_batcher::log_severity_level::normal)
                    << "!! " CPPTRANSPORT_FAILED_CONFIG << " " << list[i]->serial << " (" << i+1
                    << " " << CPPTRANSPORT_OF << " " << list.size() << ") | " << list[i]
                    << " (" << CPPTRANSPORT_FAILED_INTERNAL << xe.what() << ")";
              }
          }
      }


    template <typename number, typename StateType>
    template <typename BatchObject>
    void quartic_mpi<number, StateType>::threepf_kmode(const threepf_kconfig_record& kconfig,
                                                      const threepf_task<number>* tk,
                                                      BatchObject& batcher, unsigned int refinement_level)
      {
        DEFINE_INDEX_TOOLS

        if(refinement_level > tk->get_max_refinements()) throw runtime_exception(exception_type::REFINEMENT_FAILURE, CPPTRANSPORT_REFINEMENT_TOO_DEEP);

        // get list of time steps, and storage list
        const time_config_database time_db = tk->get_time_config_database(*kconfig);

        // set up a functor to observe the integration
        // this also starts the timers running, so we do it as early as possible
        quartic_mpi_threepf_observer< quartic_mpi<number, StateType>, BatchObject > obs(batcher, kconfig, tk->get_initial_time(*kconfig), time_db);

        // set up a functor to evolve this system
        quartic_mpi_threepf_functor< quartic_mpi<number, StateType> >  rhs(tk, *kconfig
#ifdef CPPTRANSPORT_INSTRUMENT
          ,
            this->threepf_setup_timer, this->threepf_u_tensor_timer, this->threepf_transport_eq_timer, this->threepf_invokations
#endif
          );
        rhs.set_up_workspace();

        // set up a state vector
        threepf_state x;
        x.resize(quartic_pool::threepf_state_size);

        // fix initial conditions - background
        // use adaptive ics if enabled
        // (don't need explicit FLATTEN since it would appear on both sides)
        const std::vector<number> ics = tk->get_ics_vector(*kconfig);
        x[quartic_pool::backg_start + FLATTEN(0)] = ics[0];
        x[quartic_pool::backg_start + FLATTEN(1)] = ics[1];

        if(batcher.is_collecting_initial_conditions())
          {
            const std::vector<number> ics_1 = tk->get_ics_exit_vector(*kconfig, threepf_ics_exit_type::smallest_wavenumber_exit);
            const std::vector<number> ics_2 = tk->get_ics_exit_vector(*kconfig, threepf_ics_exit_type::kt_wavenumber_exit);
            double t_exit_1 = tk->get_ics_exit_time(*kconfig, threepf_ics_exit_type::smallest_wavenumber_exit);
            double t_exit_2 = tk->get_ics_exit_time(*kconfig, threepf_ics_exit_type::kt_wavenumber_exit);
            batcher.push_ics(kconfig->serial, t_exit_1, ics_1);
            batcher.push_kt_ics(kconfig->serial, t_exit_2, ics_2);
          }

        // observers expect all correlation functions to be dimensionless and rescaled by the same factors

        // fix initial conditions - tensors (use dimensionless correlation functions, all rescaled by k_t to be consistent)
        this->populate_tensor_ic(x, quartic_pool::tensor_k1_start, kconfig->k1_comoving, *(time_db.value_begin()), tk, ics, kconfig->kt_comoving);
        this->populate_tensor_ic(x, quartic_pool::tensor_k2_start, kconfig->k2_comoving, *(time_db.value_begin()), tk, ics, kconfig->kt_comoving);
        this->populate_tensor_ic(x, quartic_pool::tensor_k3_start, kconfig->k3_comoving, *(time_db.value_begin()), tk, ics, kconfig->kt_comoving);

        // fix initial conditions - real 2pfs (use dimensionless correlation functions, all rescaled by k_t to be consistent)
        this->populate_twopf_ic(x, quartic_pool::twopf_re_k1_start, kconfig->k1_comoving, *(time_db.value_begin()), tk, ics, kconfig->kt_comoving, false);
        this->populate_twopf_ic(x, quartic_pool::twopf_re_k2_start, kconfig->k2_comoving, *(time_db.value_begin()), tk, ics, kconfig->kt_comoving, false);
        this->populate_twopf_ic(x, quartic_pool::twopf_re_k3_start, kconfig->k3_comoving, *(time_db.value_begin()), tk, ics, kconfig->kt_comoving, false);

        // fix initial conditions - imaginary 2pfs (use dimensionless correlation functions, all rescaled by k_t to be consistent)
        this->populate_twopf_ic(x, quartic_pool::twopf_im_k1_start, kconfig->k1_comoving, *(time_db.value_begin()), tk, ics, kconfig->kt_comoving, true);
        this->populate_twopf_ic(x, quartic_pool::twopf_im_k2_start, kconfig->k2_comoving, *(time_db.value_begin()), tk, ics, kconfig->kt_comoving, true);
        this->populate_twopf_ic(x, quartic_pool::twopf_im_k3_start, kconfig->k3_comoving, *(time_db.value_begin()), tk, ics, kconfig->kt_comoving, true);

        // fix initial conditions - threepf (use dimensionless correlation functions)
        this->populate_threepf_ic(x, quartic_pool::threepf_start, *kconfig, *(time_db.value_begin()), tk, ics, kconfig->kt_comoving);

        // up to this point the calculation has been done in the user-supplied time variable.
        // However, the integrator apparently performs much better if times are measured from zero (but not yet clear why)
        // TODO: would be nice to remove this in future
        rhs.rebase_horizon_exit_time(tk->get_ics().get_N_initial());
        auto begin_iterator = time_db.value_begin(tk->get_ics().get_N_initial());
        auto end_iterator   = time_db.value_end(tk->get_ics().get_N_initial());

        using boost::numeric::odeint::integrate_times;

        auto stepper = boost::numeric::odeint::make_dense_output< boost::numeric::odeint::runge_kutta_dopri5< threepf_state, number, threepf_state, number, CPPTRANSPORT_ALGEBRA_NAME(threepf_state), CPPTRANSPORT_OPERATIONS_NAME(threepf_state) > >(1e-14, 1e-14);
        size_t steps = integrate_times(stepper, rhs, x, begin_iterator, end_iterator,
                                       static_cast<number>(1.0000000000000001e-15/pow(4.0,refinement_level)), obs);

        obs.stop_timers(steps, refinement_level);
        rhs.close_down_workspace();

#ifdef CPPTRANSPORT_INSTRUMENT
        ++this->threepf_items;
#endif
      }


    template <typename number, typename StateType>
    void quartic_mpi<number, StateType>::populate_threepf_ic(threepf_state& x, unsigned int start,
                                                            const threepf_kconfig& kconfig, double Ninit,
                                                            const twopf_db_task<number>* tk,
                                                            const std::vector<number>& ics, double k_normalize)
      {
        DEFINE_INDEX_TOOLS

        assert(x.size() >= start);
        assert(x.size() >= start + quartic_pool::threepf_size);

        x[start + FLATTEN(0,0,0)] = this->make_threepf_ic(0, 0, 0, kconfig.k1_comoving, kconfig.k2_comoving, kconfig.k3_comoving, Ninit, tk, ics, k_normalize);
        x[start + FLATTEN(0,0,1)] = this->make_threepf_ic(0, 0, 1, kconfig.k1_comoving, kconfig.k2_comoving, kconfig.k3_comoving, Ninit, tk, ics, k_normalize);
        x[start + FLATTEN(0,1,0)] = this->make_threepf_ic(0, 1, 0, kconfig.k1_comoving, kconfig.k2_comoving, kconfig.k3_comoving, Ninit, tk, ics, k_normalize);
        x[start + FLATTEN(0,1,1)] = this->make_threepf_ic(0, 1, 1, kconfig.k1_comoving, kconfig.k2_comoving, kconfig.k3_comoving, Ninit, tk, ics, k_normalize);
        x[start + FLATTEN(1,0,0)] = this->make_threepf_ic(1, 0, 0, kconfig.k1_comoving, kconfig.k2_comoving, kconfig.k3_comoving, Ninit, tk, ics, k_normalize);
        x[start + FLATTEN(1,0,1)] = this->make_threepf_ic(1, 0, 1, kconfig.k1_comoving, kconfig.k2_comoving, kconfig.k3_comoving, Ninit, tk, ics, k_normalize);
        x[start + FLATTEN(1,1,0)] = this->make_threepf_ic(1, 1, 0, kconfig.k1_comoving, kconfig.k2_comoving, kconfig.k3_comoving, Ninit, tk, ics, k_normalize);
        x[start + FLATTEN(1,1,1)] = this->make_threepf_ic(1, 1, 1, kconfig.k1_comoving, kconfig.k2_comoving, kconfig.k3_comoving, Ninit, tk, ics, k_normalize);
      }


    // IMPLEMENTATION - FUNCTOR FOR 2PF INTEGRATION


    template <typename Model>
    void quartic_mpi_twopf_functor<Model>::operator()(const twopf_state& __x, twopf_state& __dxdt, number __t)
      {
        DEFINE_INDEX_TOOLS
//         release resources 

        const auto __a = std::exp(__t - this->__N_horizon_exit + this->__astar_normalization);

//         parameters resource set to '__raw_params' 
//         coordinates resource set to '__x' 

        // calculation of dV, ddV, dddV has to occur above the temporary pool
//         ENDIF !fast 

#ifdef CPPTRANSPORT_INSTRUMENT
        __setup_timer.resume();
#endif

// BEGIN TEMPORARY POOL (sequence=0) 
const auto _t_0_0 = __x[FLATTEN(1)];
const auto _t_0_1 = 2.0;
const auto _t_0_2 = (_t_0_0*_t_0_0);
const auto _t_0_3 = __Mp;
const auto _t_0_4 = -2.0;
const auto _t_0_5 = 1.0/(_t_0_3*_t_0_3);
const auto _t_0_6 = _t_0_2*_t_0_5;
const auto _t_0_7 = -6.0;
const auto _t_0_8 = _t_0_6+_t_0_7;
const auto _t_0_9 = -1.0;
const auto _t_0_10 = 1.0/(_t_0_8);
const auto _t_0_11 = __raw_params[0];
const auto _t_0_12 = __x[FLATTEN(0)];
const auto _t_0_13 = 4.0;
const auto _t_0_14 = (_t_0_12*_t_0_12*_t_0_12*_t_0_12);
const auto _t_0_15 = -(1.0/2.0);
const auto _InternalHsq = _t_0_10*_t_0_11*_t_0_5*_t_0_14*_t_0_15;
const auto _t_0_16 = (1.0/2.0);
const auto _InternalEps = _t_0_2*_t_0_5*_t_0_16;
const auto _t_0_17 = _InternalEps;
const auto _t_0_18 = -3.0;
const auto _t_0_19 = _t_0_17+_t_0_18;
const auto _t_0_20 = _t_0_0*(_t_0_19);
const auto _t_0_21 = _InternalHsq;
const auto _t_0_22 = 1.0/_t_0_21;
const auto _t_0_23 = 3.0;
const auto _t_0_24 = (_t_0_12*_t_0_12*_t_0_12);
const auto _t_0_25 = _t_0_11*_t_0_22*_t_0_24*_t_0_9;
const auto _t_0_26 = _t_0_20+_t_0_25;
const auto _t_0_27 = 0.0;
const auto _t_0_28 = 1.0;
const auto _t_0_29 = __a;
const auto _t_0_30 = 1.0/(_t_0_29*_t_0_29);
const auto _t_0_31 = __k;
const auto _t_0_32 = (_t_0_31*_t_0_31);
const auto _t_0_33 = _t_0_30*_t_0_22*_t_0_32*_t_0_9;
const auto _t_0_34 = _t_0_0*_t_0_11*_t_0_5*_t_0_22*_t_0_24*_t_0_4;
const auto _t_0_35 = _t_0_2*_t_0_5*(_t_0_19);
const auto _t_0_36 = (_t_0_12*_t_0_12);
const auto _t_0_37 = _t_0_11*_t_0_22*_t_0_36*_t_0_18;
const auto _t_0_38 = _t_0_33+_t_0_34+_t_0_35+_t_0_37;
        // END TEMPORARY POOL (sequence=0) 

        // check FLATTEN functions are being evaluated at compile time
        static_assert(TENSOR_FLATTEN(0,0) == 0, "TENSOR_FLATTEN failure");
        static_assert(FLATTEN(0) == 0, "FLATTEN failure");
        static_assert(FLATTEN(0,0) == 0, "FLATTEN failure");
        static_assert(FLATTEN(0,0,0) == 0, "FLATTEN failure");

        const auto __tensor_twopf_ff = __x[quartic_pool::tensor_start + TENSOR_FLATTEN(0,0)];
        const auto __tensor_twopf_fp = __x[quartic_pool::tensor_start + TENSOR_FLATTEN(0,1)];
        const auto __tensor_twopf_pf = __x[quartic_pool::tensor_start + TENSOR_FLATTEN(1,0)];
        const auto __tensor_twopf_pp = __x[quartic_pool::tensor_start + TENSOR_FLATTEN(1,1)];

#undef __twopf
#define __twopf(a,b) __x[quartic_pool::twopf_start + FLATTEN(a,b)]

#undef __background
#undef __dtwopf
#undef __dtwopf_tensor
#define __background(a)      __dxdt[quartic_pool::backg_start + FLATTEN(a)]
#define __dtwopf_tensor(a,b) __dxdt[quartic_pool::tensor_start + TENSOR_FLATTEN(a,b)]
#define __dtwopf(a,b)        __dxdt[quartic_pool::twopf_start + FLATTEN(a,b)]

#ifdef CPPTRANSPORT_INSTRUMENT
        __setup_timer.stop();
        __u_tensor_timer.resume();
#endif

        // evolve the background
        __background(0) = _t_0_0;
        __background(1) = _t_0_26;

        const auto __Hsq = _t_0_21;
        const auto __eps = _t_0_17;

        // evolve the tensor modes
        const auto __ff = 0.0;
        const auto __fp = 1.0;
        const auto __pf = -__k*__k/(__a*__a*__Hsq);
        const auto __pp = __eps-3.0;
        __dtwopf_tensor(0,0) = __ff*__tensor_twopf_ff + __fp*__tensor_twopf_pf + __ff*__tensor_twopf_ff + __fp*__tensor_twopf_fp;
        __dtwopf_tensor(0,1) = __ff*__tensor_twopf_fp + __fp*__tensor_twopf_pp + __pf*__tensor_twopf_ff + __pp*__tensor_twopf_fp;
        __dtwopf_tensor(1,0) = __pf*__tensor_twopf_ff + __pp*__tensor_twopf_pf + __ff*__tensor_twopf_pf + __fp*__tensor_twopf_pp;
        __dtwopf_tensor(1,1) = __pf*__tensor_twopf_fp + __pp*__tensor_twopf_pp + __pf*__tensor_twopf_pf + __pp*__tensor_twopf_pp;

        // set up components of the u2 tensor
        const auto __u2_0_0 = _t_0_27;
        const auto __u2_0_1 = _t_0_28;
        const auto __u2_1_0 = _t_0_38;
        const auto __u2_1_1 = _t_0_19;

#ifdef CPPTRANSPORT_INSTRUMENT
        __u_tensor_timer.stop();
        __transport_eq_timer.resume();
#endif

        // evolve the 2pf
        // here, we are dealing only with the real part - which is symmetric.
        // so the index placement is not important
        __dtwopf(0, 0)  = 
           + __u2_0_0 * __twopf(0, 0)
           + __u2_0_1 * __twopf(1, 0);
        __dtwopf(0, 1)  = 
           + __u2_0_0 * __twopf(0, 1)
           + __u2_0_1 * __twopf(1, 1);
        __dtwopf(1, 0)  = 
           + __u2_1_0 * __twopf(0, 0)
           + __u2_1_1 * __twopf(1, 0);
        __dtwopf(1, 1)  = 
           + __u2_1_0 * __twopf(0, 1)
           + __u2_1_1 * __twopf(1, 1);
        __dtwopf(0, 0)  += 
           + __u2_0_0 * __twopf(0, 0)
           + __u2_0_1 * __twopf(0, 1);
        __dtwopf(0, 1)  += 
           + __u2_1_0 * __twopf(0, 0)
           + __u2_1_1 * __twopf(0, 1);
        __dtwopf(1, 0)  += 
           + __u2_0_0 * __twopf(1, 0)
           + __u2_0_1 * __twopf(1, 1);
        __dtwopf(1, 1)  += 
           + __u2_1_0 * __twopf(1, 0)
           + __u2_1_1 * __twopf(1, 1);

#ifdef CPPTRANSPORT_STRICT_FP_TEST
        if(std::isnan(__background(0)) || std::isinf(__background(0))) throw runtime_exception(exception_type::INTEGRATION_FAILURE, CPPTRANSPORT_INTEGRATOR_NAN_OR_INF);
        if(std::isnan(__background(1)) || std::isinf(__background(1))) throw runtime_exception(exception_type::INTEGRATION_FAILURE, CPPTRANSPORT_INTEGRATOR_NAN_OR_INF);

        if(std::isnan(__dtwopf_tensor(0,0)) || std::isinf(__dtwopf_tensor(0,0))) throw runtime_exception(exception_type::INTEGRATION_FAILURE, CPPTRANSPORT_INTEGRATOR_NAN_OR_INF);
        if(std::isnan(__dtwopf_tensor(0,1)) || std::isinf(__dtwopf_tensor(0,1))) throw runtime_exception(exception_type::INTEGRATION_FAILURE, CPPTRANSPORT_INTEGRATOR_NAN_OR_INF);
        if(std::isnan(__dtwopf_tensor(1,0)) || std::isinf(__dtwopf_tensor(1,0))) throw runtime_exception(exception_type::INTEGRATION_FAILURE, CPPTRANSPORT_INTEGRATOR_NAN_OR_INF);
        if(std::isnan(__dtwopf_tensor(1,1)) || std::isinf(__dtwopf_tensor(1,1))) throw runtime_exception(exception_type::INTEGRATION_FAILURE, CPPTRANSPORT_INTEGRATOR_NAN_OR_INF);

        if(std::isnan(__dtwopf(0, 0)) || std::isinf(__dtwopf(0, 0))) throw runtime_exception(exception_type::INTEGRATION_FAILURE, CPPTRANSPORT_INTEGRATOR_NAN_OR_INF);
        if(std::isnan(__dtwopf(0, 1)) || std::isinf(__dtwopf(0, 1))) throw runtime_exception(exception_type::INTEGRATION_FAILURE, CPPTRANSPORT_INTEGRATOR_NAN_OR_INF);
        if(std::isnan(__dtwopf(1, 0)) || std::isinf(__dtwopf(1, 0))) throw runtime_exception(exception_type::INTEGRATION_FAILURE, CPPTRANSPORT_INTEGRATOR_NAN_OR_INF);
        if(std::isnan(__dtwopf(1, 1)) || std::isinf(__dtwopf(1, 1))) throw runtime_exception(exception_type::INTEGRATION_FAILURE, CPPTRANSPORT_INTEGRATOR_NAN_OR_INF);
#endif

#ifdef CPPTRANSPORT_INSTRUMENT
        __transport_eq_timer.stop();
        ++__invokations;
#endif
      }


    // IMPLEMENTATION - FUNCTOR FOR 2PF OBSERVATION


    template <typename Model, typename BatchObject>
    void quartic_mpi_twopf_observer<Model, BatchObject>::operator()(const twopf_state& x, number t)
      {
        DEFINE_INDEX_TOOLS

#undef __background
#undef __twopf
#undef __twopf_tensor

#define __background(a)      x[quartic_pool::backg_start + FLATTEN(a)]
#define __twopf(a,b)         x[quartic_pool::twopf_start + FLATTEN(a,b)]
#define __twopf_tensor(a,b) x[quartic_pool::tensor_start + TENSOR_FLATTEN(a,b)]

#ifndef CPPTRANSPORT_NO_STRICT_FP_TEST
        if(std::isnan(__background(0)) || std::isinf(__background(0))) throw runtime_exception(exception_type::INTEGRATION_FAILURE, CPPTRANSPORT_INTEGRATOR_NAN_OR_INF);
        if(std::isnan(__background(1)) || std::isinf(__background(1))) throw runtime_exception(exception_type::INTEGRATION_FAILURE, CPPTRANSPORT_INTEGRATOR_NAN_OR_INF);

        if(std::isnan(__twopf_tensor(0,0)) || std::isinf(__twopf_tensor(0,0))) throw runtime_exception(exception_type::INTEGRATION_FAILURE, CPPTRANSPORT_INTEGRATOR_NAN_OR_INF);
        if(std::isnan(__twopf_tensor(0,1)) || std::isinf(__twopf_tensor(0,1))) throw runtime_exception(exception_type::INTEGRATION_FAILURE, CPPTRANSPORT_INTEGRATOR_NAN_OR_INF);
        if(std::isnan(__twopf_tensor(1,0)) || std::isinf(__twopf_tensor(1,0))) throw runtime_exception(exception_type::INTEGRATION_FAILURE, CPPTRANSPORT_INTEGRATOR_NAN_OR_INF);
        if(std::isnan(__twopf_tensor(1,1)) || std::isinf(__twopf_tensor(1,1))) throw runtime_exception(exception_type::INTEGRATION_FAILURE, CPPTRANSPORT_INTEGRATOR_NAN_OR_INF);

        if(std::isnan(__twopf(0, 0)) || std::isinf(__twopf(0, 0))) throw runtime_exception(exception_type::INTEGRATION_FAILURE, CPPTRANSPORT_INTEGRATOR_NAN_OR_INF);
        if(std::isnan(__twopf(0, 1)) || std::isinf(__twopf(0, 1))) throw runtime_exception(exception_type::INTEGRATION_FAILURE, CPPTRANSPORT_INTEGRATOR_NAN_OR_INF);
        if(std::isnan(__twopf(1, 0)) || std::isinf(__twopf(1, 0))) throw runtime_exception(exception_type::INTEGRATION_FAILURE, CPPTRANSPORT_INTEGRATOR_NAN_OR_INF);
        if(std::isnan(__twopf(1, 1)) || std::isinf(__twopf(1, 1))) throw runtime_exception(exception_type::INTEGRATION_FAILURE, CPPTRANSPORT_INTEGRATOR_NAN_OR_INF);
#endif

        this->start_batching(static_cast<double>(t), this->get_log(), BatchObject::log_severity_level::normal);
        this->push(x);
        this->stop_batching();
      }


    // IMPLEMENTATION - FUNCTOR FOR 3PF INTEGRATION


    template <typename Model>
    void quartic_mpi_threepf_functor<Model>::operator()(const threepf_state& __x, threepf_state& __dxdt, number __t)
      {
        DEFINE_INDEX_TOOLS
//         release resources 

        const auto __a = std::exp(__t - this->__N_horizon_exit + this->__astar_normalization);

//         parameters resource set to '__raw_params' 
//         coordinates resource set to '__x' 

        // calculation of dV, ddV, dddV has to occur above the temporary pool
//         ENDIF !fast 

#ifdef CPPTRANSPORT_INSTRUMENT
        __setup_timer.resume();
#endif

// BEGIN TEMPORARY POOL (sequence=1) 
const auto _t_1_0 = __x[FLATTEN(1)];
const auto _t_1_1 = 2.0;
const auto _t_1_2 = (_t_1_0*_t_1_0);
const auto _t_1_3 = __Mp;
const auto _t_1_4 = -2.0;
const auto _t_1_5 = 1.0/(_t_1_3*_t_1_3);
const auto _t_1_6 = _t_1_2*_t_1_5;
const auto _t_1_7 = -6.0;
const auto _t_1_8 = _t_1_6+_t_1_7;
const auto _t_1_9 = -1.0;
const auto _t_1_10 = 1.0/(_t_1_8);
const auto _t_1_11 = __raw_params[0];
const auto _t_1_12 = __x[FLATTEN(0)];
const auto _t_1_13 = 4.0;
const auto _t_1_14 = (_t_1_12*_t_1_12*_t_1_12*_t_1_12);
const auto _t_1_15 = -(1.0/2.0);
const auto _InternalHsq = _t_1_10*_t_1_11*_t_1_5*_t_1_14*_t_1_15;
const auto _t_1_16 = (1.0/2.0);
const auto _InternalEps = _t_1_2*_t_1_5*_t_1_16;
const auto _t_1_17 = _InternalEps;
const auto _t_1_18 = -3.0;
const auto _t_1_19 = _t_1_17+_t_1_18;
const auto _t_1_20 = _t_1_0*(_t_1_19);
const auto _t_1_21 = _InternalHsq;
const auto _t_1_22 = 1.0/_t_1_21;
const auto _t_1_23 = 3.0;
const auto _t_1_24 = (_t_1_12*_t_1_12*_t_1_12);
const auto _t_1_25 = _t_1_11*_t_1_22*_t_1_24*_t_1_9;
const auto _t_1_26 = _t_1_20+_t_1_25;
const auto _t_1_27 = 0.0;
const auto _t_1_28 = 1.0;
const auto _t_1_29 = __a;
const auto _t_1_30 = 1.0/(_t_1_29*_t_1_29);
const auto _t_1_31 = __k1;
const auto _t_1_32 = (_t_1_31*_t_1_31);
const auto _t_1_33 = _t_1_30*_t_1_32*_t_1_22*_t_1_9;
const auto _t_1_34 = _t_1_0*_t_1_11*_t_1_5*_t_1_22*_t_1_24*_t_1_4;
const auto _t_1_35 = _t_1_2*_t_1_5*(_t_1_19);
const auto _t_1_36 = (_t_1_12*_t_1_12);
const auto _t_1_37 = _t_1_11*_t_1_22*_t_1_36*_t_1_18;
const auto _t_1_38 = _t_1_33+_t_1_34+_t_1_35+_t_1_37;
const auto _t_1_39 = __k2;
const auto _t_1_40 = (_t_1_39*_t_1_39);
const auto _t_1_41 = _t_1_30*_t_1_40*_t_1_22*_t_1_9;
const auto _t_1_42 = _t_1_34+_t_1_35+_t_1_41+_t_1_37;
const auto _t_1_43 = __k3;
const auto _t_1_44 = (_t_1_43*_t_1_43);
const auto _t_1_45 = _t_1_30*_t_1_22*_t_1_44*_t_1_9;
const auto _t_1_46 = _t_1_34+_t_1_35+_t_1_45+_t_1_37;
const auto _t_1_47 = (_t_1_0*_t_1_0*_t_1_0);
const auto _t_1_48 = -4.0;
const auto _t_1_49 = 1.0/(_t_1_3*_t_1_3*_t_1_3*_t_1_3);
const auto _t_1_50 = -(1.0/4.0);
const auto _t_1_51 = _t_1_47*_t_1_49*_t_1_50;
const auto _t_1_52 = _t_1_44*_t_1_9;
const auto _t_1_53 = _t_1_40*_t_1_9;
const auto _t_1_54 = _t_1_32+_t_1_52+_t_1_53;
const auto _t_1_55 = 1.0/(_t_1_43*_t_1_43);
const auto _t_1_56 = (1.0/4.0);
const auto _t_1_57 = _t_1_5*(_t_1_54)*(_t_1_26)*_t_1_55*_t_1_56;
const auto _t_1_58 = 1.0/(_t_1_39*_t_1_39);
const auto _t_1_59 = _t_1_5*_t_1_58*(_t_1_54)*(_t_1_26)*_t_1_56;
const auto _t_1_60 = _t_1_32+_t_1_44+_t_1_53;
const auto _t_1_61 = ((_t_1_60)*(_t_1_60));
const auto _t_1_62 = 1.0/(_t_1_31*_t_1_31);
const auto _t_1_63 = _t_1_61*_t_1_62*_t_1_55;
const auto _t_1_64 = _t_1_63+_t_1_48;
const auto _t_1_65 = -(1.0/32.0);
const auto _t_1_66 = _t_1_2*_t_1_49*(_t_1_26)*(_t_1_64)*_t_1_65;
const auto _t_1_67 = _t_1_32+_t_1_52+_t_1_40;
const auto _t_1_68 = ((_t_1_67)*(_t_1_67));
const auto _t_1_69 = _t_1_62*_t_1_68*_t_1_58;
const auto _t_1_70 = _t_1_69+_t_1_48;
const auto _t_1_71 = _t_1_2*_t_1_49*(_t_1_26)*(_t_1_70)*_t_1_65;
const auto _t_1_72 = _t_1_51+_t_1_57+_t_1_59+_t_1_66+_t_1_71;
const auto _t_1_73 = (1.0/32.0);
const auto _t_1_74 = _t_1_47*_t_1_49*(_t_1_64)*_t_1_73;
const auto _t_1_75 = _t_1_0*_t_1_5*(_t_1_54)*_t_1_55*_t_1_50;
const auto _t_1_76 = _t_1_0*_t_1_5*_t_1_16;
const auto _t_1_77 = _t_1_0*_t_1_62*(_t_1_67)*_t_1_5*_t_1_56;
const auto _t_1_78 = _t_1_74+_t_1_75+_t_1_76+_t_1_77;
const auto _t_1_79 = (_t_1_60)*_t_1_0*_t_1_62*_t_1_5*_t_1_56;
const auto _t_1_80 = _t_1_47*_t_1_49*(_t_1_70)*_t_1_73;
const auto _t_1_81 = _t_1_0*_t_1_5*_t_1_58*(_t_1_54)*_t_1_50;
const auto _t_1_82 = _t_1_79+_t_1_80+_t_1_81+_t_1_76;
const auto _t_1_83 = _t_1_30*_t_1_0*(_t_1_67)*_t_1_5*_t_1_22*_t_1_50;
const auto _t_1_84 = -(9.0/2.0);
const auto _t_1_85 = _t_1_0*_t_1_11*_t_1_5*_t_1_22*_t_1_36*_t_1_84;
const auto _t_1_86 = _t_1_30*_t_1_0*_t_1_5*(_t_1_54)*_t_1_22*_t_1_56;
const auto _t_1_87 = (_t_1_60)*_t_1_30*_t_1_0*_t_1_5*_t_1_22*_t_1_50;
const auto _t_1_88 = ((_t_1_26)*(_t_1_26));
const auto _t_1_89 = ((_t_1_54)*(_t_1_54));
const auto _t_1_90 = _t_1_58*_t_1_89*_t_1_55;
const auto _t_1_91 = _t_1_90+_t_1_48;
const auto _t_1_92 = _t_1_0*_t_1_49*_t_1_88*(_t_1_91)*_t_1_65;
const auto _t_1_93 = _t_1_11*_t_1_22*_t_1_12*_t_1_7;
const auto _t_1_94 = -(3.0/4.0);
const auto _t_1_95 = _t_1_47*_t_1_49*(_t_1_19)*_t_1_94;
const auto _t_1_96 = _t_1_0*_t_1_49*_t_1_88*(_t_1_70)*_t_1_65;
const auto _t_1_97 = _t_1_0*_t_1_49*_t_1_88*(_t_1_64)*_t_1_65;
const auto _t_1_98 = (3.0/4.0);
const auto _t_1_99 = _t_1_2*_t_1_49*(_t_1_26)*_t_1_98;
const auto _t_1_100 = _t_1_83+_t_1_85+_t_1_86+_t_1_87+_t_1_92+_t_1_93+_t_1_95+_t_1_96+_t_1_97+_t_1_99;
const auto _t_1_101 = _t_1_2*_t_1_49*(_t_1_26)*(_t_1_91)*_t_1_73;
const auto _t_1_102 = _t_1_47*_t_1_49*_t_1_56;
const auto _t_1_103 = (_t_1_67)*_t_1_5*_t_1_58*(_t_1_26)*_t_1_56;
const auto _t_1_104 = _t_1_2*_t_1_49*(_t_1_26)*(_t_1_64)*_t_1_73;
const auto _t_1_105 = _t_1_62*(_t_1_67)*_t_1_5*(_t_1_26)*_t_1_56;
const auto _t_1_106 = _t_1_101+_t_1_102+_t_1_103+_t_1_104+_t_1_105;
const auto _t_1_107 = _t_1_2*_t_1_49*(_t_1_26)*(_t_1_70)*_t_1_73;
const auto _t_1_108 = (_t_1_60)*_t_1_5*(_t_1_26)*_t_1_55*_t_1_56;
const auto _t_1_109 = (_t_1_60)*_t_1_62*_t_1_5*(_t_1_26)*_t_1_56;
const auto _t_1_110 = _t_1_101+_t_1_102+_t_1_107+_t_1_108+_t_1_109;
const auto _t_1_111 = _t_1_0*(_t_1_67)*_t_1_5*_t_1_58*_t_1_50;
const auto _t_1_112 = _t_1_47*_t_1_49*(_t_1_91)*_t_1_65;
const auto _t_1_113 = _t_1_0*_t_1_5*_t_1_15;
const auto _t_1_114 = (_t_1_60)*_t_1_0*_t_1_5*_t_1_55*_t_1_50;
const auto _t_1_115 = _t_1_111+_t_1_112+_t_1_113+_t_1_114;
const auto _t_1_116 = _t_1_2*_t_1_49*(_t_1_26)*(_t_1_91)*_t_1_65;
const auto _t_1_117 = (_t_1_60)*_t_1_5*(_t_1_26)*_t_1_55*_t_1_50;
const auto _t_1_118 = (_t_1_60)*_t_1_62*_t_1_5*(_t_1_26)*_t_1_50;
const auto _t_1_119 = _t_1_116+_t_1_51+_t_1_71+_t_1_117+_t_1_118;
const auto _t_1_120 = _t_1_0*(_t_1_67)*_t_1_5*_t_1_58*_t_1_56;
const auto _t_1_121 = _t_1_47*_t_1_49*(_t_1_91)*_t_1_73;
const auto _t_1_122 = (_t_1_60)*_t_1_0*_t_1_5*_t_1_55*_t_1_56;
const auto _t_1_123 = _t_1_120+_t_1_121+_t_1_76+_t_1_122;
const auto _t_1_124 = _t_1_5*(_t_1_54)*(_t_1_26)*_t_1_55*_t_1_50;
const auto _t_1_125 = _t_1_5*_t_1_58*(_t_1_54)*(_t_1_26)*_t_1_50;
const auto _t_1_126 = _t_1_102+_t_1_124+_t_1_125+_t_1_104+_t_1_107;
const auto _t_1_127 = _t_1_47*_t_1_49*(_t_1_64)*_t_1_65;
const auto _t_1_128 = _t_1_0*_t_1_5*(_t_1_54)*_t_1_55*_t_1_56;
const auto _t_1_129 = _t_1_0*_t_1_62*(_t_1_67)*_t_1_5*_t_1_50;
const auto _t_1_130 = _t_1_127+_t_1_128+_t_1_113+_t_1_129;
const auto _t_1_131 = (_t_1_67)*_t_1_5*_t_1_58*(_t_1_26)*_t_1_50;
const auto _t_1_132 = _t_1_62*(_t_1_67)*_t_1_5*(_t_1_26)*_t_1_50;
const auto _t_1_133 = _t_1_116+_t_1_51+_t_1_131+_t_1_66+_t_1_132;
const auto _t_1_134 = (_t_1_60)*_t_1_0*_t_1_62*_t_1_5*_t_1_50;
const auto _t_1_135 = _t_1_47*_t_1_49*(_t_1_70)*_t_1_65;
const auto _t_1_136 = _t_1_0*_t_1_5*_t_1_58*(_t_1_54)*_t_1_56;
const auto _t_1_137 = _t_1_134+_t_1_135+_t_1_136+_t_1_113;
        // END TEMPORARY POOL (sequence=1) 

        // check FLATTEN functions are being evaluated at compile time
        static_assert(TENSOR_FLATTEN(0,0) == 0, "TENSOR_FLATTEN failure");
        static_assert(FLATTEN(0) == 0, "FLATTEN failure");
        static_assert(FLATTEN(0,0) == 0, "FLATTEN failure");
        static_assert(FLATTEN(0,0,0) == 0, "FLATTEN failure");

        const auto __tensor_k1_twopf_ff = __x[quartic_pool::tensor_k1_start + TENSOR_FLATTEN(0,0)];
        const auto __tensor_k1_twopf_fp = __x[quartic_pool::tensor_k1_start + TENSOR_FLATTEN(0,1)];
        const auto __tensor_k1_twopf_pf = __x[quartic_pool::tensor_k1_start + TENSOR_FLATTEN(1,0)];
        const auto __tensor_k1_twopf_pp = __x[quartic_pool::tensor_k1_start + TENSOR_FLATTEN(1,1)];

        const auto __tensor_k2_twopf_ff = __x[quartic_pool::tensor_k2_start + TENSOR_FLATTEN(0,0)];
        const auto __tensor_k2_twopf_fp = __x[quartic_pool::tensor_k2_start + TENSOR_FLATTEN(0,1)];
        const auto __tensor_k2_twopf_pf = __x[quartic_pool::tensor_k2_start + TENSOR_FLATTEN(1,0)];
        const auto __tensor_k2_twopf_pp = __x[quartic_pool::tensor_k2_start + TENSOR_FLATTEN(1,1)];

        const auto __tensor_k3_twopf_ff = __x[quartic_pool::tensor_k3_start + TENSOR_FLATTEN(0,0)];
        const auto __tensor_k3_twopf_fp = __x[quartic_pool::tensor_k3_start + TENSOR_FLATTEN(0,1)];
        const auto __tensor_k3_twopf_pf = __x[quartic_pool::tensor_k3_start + TENSOR_FLATTEN(1,0)];
        const auto __tensor_k3_twopf_pp = __x[quartic_pool::tensor_k3_start + TENSOR_FLATTEN(1,1)];

#undef __twopf_re_k1
#undef __twopf_re_k2
#undef __twopf_re_k3
#undef __twopf_im_k1
#undef __twopf_im_k2
#undef __twopf_im_k3

#undef __threepf

#define __twopf_re_k1(a,b) __x[quartic_pool::twopf_re_k1_start + FLATTEN(a,b)]
#define __twopf_im_k1(a,b) __x[quartic_pool::twopf_im_k1_start + FLATTEN(a,b)]
#define __twopf_re_k2(a,b) __x[quartic_pool::twopf_re_k2_start + FLATTEN(a,b)]
#define __twopf_im_k2(a,b) __x[quartic_pool::twopf_im_k2_start + FLATTEN(a,b)]
#define __twopf_re_k3(a,b) __x[quartic_pool::twopf_re_k3_start + FLATTEN(a,b)]
#define __twopf_im_k3(a,b) __x[quartic_pool::twopf_im_k3_start + FLATTEN(a,b)]

#define __threepf(a,b,c)	 __x[quartic_pool::threepf_start  + FLATTEN(a,b,c)]

#undef __background
#undef __dtwopf_k1_tensor
#undef __dtwopf_k2_tensor
#undef __dtwopf_k3_tensor
#undef __dtwopf_re_k1
#undef __dtwopf_im_k1
#undef __dtwopf_re_k2
#undef __dtwopf_im_k2
#undef __dtwopf_re_k3
#undef __dtwopf_im_k3
#undef __dthreepf
#define __background(a)         __dxdt[quartic_pool::backg_start       + FLATTEN(a)]
#define __dtwopf_k1_tensor(a,b) __dxdt[quartic_pool::tensor_k1_start   + TENSOR_FLATTEN(a,b)]
#define __dtwopf_k2_tensor(a,b) __dxdt[quartic_pool::tensor_k2_start   + TENSOR_FLATTEN(a,b)]
#define __dtwopf_k3_tensor(a,b) __dxdt[quartic_pool::tensor_k3_start   + TENSOR_FLATTEN(a,b)]
#define __dtwopf_re_k1(a,b)     __dxdt[quartic_pool::twopf_re_k1_start + FLATTEN(a,b)]
#define __dtwopf_im_k1(a,b)     __dxdt[quartic_pool::twopf_im_k1_start + FLATTEN(a,b)]
#define __dtwopf_re_k2(a,b)     __dxdt[quartic_pool::twopf_re_k2_start + FLATTEN(a,b)]
#define __dtwopf_im_k2(a,b)     __dxdt[quartic_pool::twopf_im_k2_start + FLATTEN(a,b)]
#define __dtwopf_re_k3(a,b)     __dxdt[quartic_pool::twopf_re_k3_start + FLATTEN(a,b)]
#define __dtwopf_im_k3(a,b)     __dxdt[quartic_pool::twopf_im_k3_start + FLATTEN(a,b)]
#define __dthreepf(a,b,c)       __dxdt[quartic_pool::threepf_start     + FLATTEN(a,b,c)]

#ifdef CPPTRANSPORT_INSTRUMENT
        __setup_timer.stop();
        __u_tensor_timer.resume();
#endif

        // evolve the background
        __background(0) = _t_1_0;
        __background(1) = _t_1_26;

        const auto __Hsq = _t_1_21;
        const auto __eps = _t_1_17;

        // evolve the tensor modes
        const auto __ff = 0.0;
        const auto __fp = 1.0;
        const auto __pp = __eps-3.0;

        auto __pf = -__k1*__k1/(__a*__a*__Hsq);
        __dtwopf_k1_tensor(0,0) = __ff*__tensor_k1_twopf_ff + __fp*__tensor_k1_twopf_pf + __ff*__tensor_k1_twopf_ff + __fp*__tensor_k1_twopf_fp;
        __dtwopf_k1_tensor(0,1) = __ff*__tensor_k1_twopf_fp + __fp*__tensor_k1_twopf_pp + __pf*__tensor_k1_twopf_ff + __pp*__tensor_k1_twopf_fp;
        __dtwopf_k1_tensor(1,0) = __pf*__tensor_k1_twopf_ff + __pp*__tensor_k1_twopf_pf + __ff*__tensor_k1_twopf_pf + __fp*__tensor_k1_twopf_pp;
        __dtwopf_k1_tensor(1,1) = __pf*__tensor_k1_twopf_fp + __pp*__tensor_k1_twopf_pp + __pf*__tensor_k1_twopf_pf + __pp*__tensor_k1_twopf_pp;

        __pf = -__k2*__k2/(__a*__a*__Hsq);
        __dtwopf_k2_tensor(0,0) = __ff*__tensor_k2_twopf_ff + __fp*__tensor_k2_twopf_pf + __ff*__tensor_k2_twopf_ff + __fp*__tensor_k2_twopf_fp;
        __dtwopf_k2_tensor(0,1) = __ff*__tensor_k2_twopf_fp + __fp*__tensor_k2_twopf_pp + __pf*__tensor_k2_twopf_ff + __pp*__tensor_k2_twopf_fp;
        __dtwopf_k2_tensor(1,0) = __pf*__tensor_k2_twopf_ff + __pp*__tensor_k2_twopf_pf + __ff*__tensor_k2_twopf_pf + __fp*__tensor_k2_twopf_pp;
        __dtwopf_k2_tensor(1,1) = __pf*__tensor_k2_twopf_fp + __pp*__tensor_k2_twopf_pp + __pf*__tensor_k2_twopf_pf + __pp*__tensor_k2_twopf_pp;

        __pf = -__k3*__k3/(__a*__a*__Hsq);
        __dtwopf_k3_tensor(0,0) = __ff*__tensor_k3_twopf_ff + __fp*__tensor_k3_twopf_pf + __ff*__tensor_k3_twopf_ff + __fp*__tensor_k3_twopf_fp;
        __dtwopf_k3_tensor(0,1) = __ff*__tensor_k3_twopf_fp + __fp*__tensor_k3_twopf_pp + __pf*__tensor_k3_twopf_ff + __pp*__tensor_k3_twopf_fp;
        __dtwopf_k3_tensor(1,0) = __pf*__tensor_k3_twopf_ff + __pp*__tensor_k3_twopf_pf + __ff*__tensor_k3_twopf_pf + __fp*__tensor_k3_twopf_pp;
        __dtwopf_k3_tensor(1,1) = __pf*__tensor_k3_twopf_fp + __pp*__tensor_k3_twopf_pp + __pf*__tensor_k3_twopf_pf + __pp*__tensor_k3_twopf_pp;

        // set up components of the u2 tensor for k1, k2, k3
        const auto __u2_k1_0_0 = _t_1_27;
        const auto __u2_k1_0_1 = _t_1_28;
        const auto __u2_k1_1_0 = _t_1_38;
        const auto __u2_k1_1_1 = _t_1_19;
        const auto __u2_k2_0_0 = _t_1_27;
        const auto __u2_k2_0_1 = _t_1_28;
        const auto __u2_k2_1_0 = _t_1_42;
        const auto __u2_k2_1_1 = _t_1_19;
        const auto __u2_k3_0_0 = _t_1_27;
        const auto __u2_k3_0_1 = _t_1_28;
        const auto __u2_k3_1_0 = _t_1_46;
        const auto __u2_k3_1_1 = _t_1_19;

        // set up components of the u3 tensor
        const auto __u3_k1k2k3_0_0_0 = _t_1_72;
        const auto __u3_k1k2k3_0_0_1 = _t_1_78;
        const auto __u3_k1k2k3_0_1_0 = _t_1_82;
        const auto __u3_k1k2k3_0_1_1 = _t_1_27;
        const auto __u3_k1k2k3_1_0_0 = _t_1_100;
        const auto __u3_k1k2k3_1_0_1 = _t_1_106;
        const auto __u3_k1k2k3_1_1_0 = _t_1_110;
        const auto __u3_k1k2k3_1_1_1 = _t_1_115;
        const auto __u3_k2k1k3_0_0_0 = _t_1_119;
        const auto __u3_k2k1k3_0_0_1 = _t_1_123;
        const auto __u3_k2k1k3_0_1_0 = _t_1_82;
        const auto __u3_k2k1k3_0_1_1 = _t_1_27;
        const auto __u3_k2k1k3_1_0_0 = _t_1_100;
        const auto __u3_k2k1k3_1_0_1 = _t_1_106;
        const auto __u3_k2k1k3_1_1_0 = _t_1_126;
        const auto __u3_k2k1k3_1_1_1 = _t_1_130;
        const auto __u3_k3k1k2_0_0_0 = _t_1_133;
        const auto __u3_k3k1k2_0_0_1 = _t_1_123;
        const auto __u3_k3k1k2_0_1_0 = _t_1_78;
        const auto __u3_k3k1k2_0_1_1 = _t_1_27;
        const auto __u3_k3k1k2_1_0_0 = _t_1_100;
        const auto __u3_k3k1k2_1_0_1 = _t_1_110;
        const auto __u3_k3k1k2_1_1_0 = _t_1_126;
        const auto __u3_k3k1k2_1_1_1 = _t_1_137;

#ifdef CPPTRANSPORT_INSTRUMENT
        __u_tensor_timer.stop();
        __transport_eq_timer.resume();
#endif

        // evolve the real and imaginary components of the 2pf
        // for the imaginary parts, index placement *does* matter so we must take care
        __dtwopf_re_k1(0, 0)  = 
           + __u2_k1_0_0 * __twopf_re_k1(0, 0)
           + __u2_k1_0_1 * __twopf_re_k1(1, 0);
        __dtwopf_re_k1(0, 1)  = 
           + __u2_k1_0_0 * __twopf_re_k1(0, 1)
           + __u2_k1_0_1 * __twopf_re_k1(1, 1);
        __dtwopf_re_k1(1, 0)  = 
           + __u2_k1_1_0 * __twopf_re_k1(0, 0)
           + __u2_k1_1_1 * __twopf_re_k1(1, 0);
        __dtwopf_re_k1(1, 1)  = 
           + __u2_k1_1_0 * __twopf_re_k1(0, 1)
           + __u2_k1_1_1 * __twopf_re_k1(1, 1);
        __dtwopf_re_k1(0, 0)  += 
           + __u2_k1_0_0 * __twopf_re_k1(0, 0)
           + __u2_k1_0_1 * __twopf_re_k1(0, 1);
        __dtwopf_re_k1(0, 1)  += 
           + __u2_k1_1_0 * __twopf_re_k1(0, 0)
           + __u2_k1_1_1 * __twopf_re_k1(0, 1);
        __dtwopf_re_k1(1, 0)  += 
           + __u2_k1_0_0 * __twopf_re_k1(1, 0)
           + __u2_k1_0_1 * __twopf_re_k1(1, 1);
        __dtwopf_re_k1(1, 1)  += 
           + __u2_k1_1_0 * __twopf_re_k1(1, 0)
           + __u2_k1_1_1 * __twopf_re_k1(1, 1);

        __dtwopf_im_k1(0, 0)  = 
           + __u2_k1_0_0 * __twopf_im_k1(0, 0)
           + __u2_k1_0_1 * __twopf_im_k1(1, 0);
        __dtwopf_im_k1(0, 1)  = 
           + __u2_k1_0_0 * __twopf_im_k1(0, 1)
           + __u2_k1_0_1 * __twopf_im_k1(1, 1);
        __dtwopf_im_k1(1, 0)  = 
           + __u2_k1_1_0 * __twopf_im_k1(0, 0)
           + __u2_k1_1_1 * __twopf_im_k1(1, 0);
        __dtwopf_im_k1(1, 1)  = 
           + __u2_k1_1_0 * __twopf_im_k1(0, 1)
           + __u2_k1_1_1 * __twopf_im_k1(1, 1);
        __dtwopf_im_k1(0, 0)  += 
           + __u2_k1_0_0 * __twopf_im_k1(0, 0)
           + __u2_k1_0_1 * __twopf_im_k1(0, 1);
        __dtwopf_im_k1(0, 1)  += 
           + __u2_k1_1_0 * __twopf_im_k1(0, 0)
           + __u2_k1_1_1 * __twopf_im_k1(0, 1);
        __dtwopf_im_k1(1, 0)  += 
           + __u2_k1_0_0 * __twopf_im_k1(1, 0)
           + __u2_k1_0_1 * __twopf_im_k1(1, 1);
        __dtwopf_im_k1(1, 1)  += 
           + __u2_k1_1_0 * __twopf_im_k1(1, 0)
           + __u2_k1_1_1 * __twopf_im_k1(1, 1);

        __dtwopf_re_k2(0, 0)  = 
           + __u2_k2_0_0 * __twopf_re_k2(0, 0)
           + __u2_k2_0_1 * __twopf_re_k2(1, 0);
        __dtwopf_re_k2(0, 1)  = 
           + __u2_k2_0_0 * __twopf_re_k2(0, 1)
           + __u2_k2_0_1 * __twopf_re_k2(1, 1);
        __dtwopf_re_k2(1, 0)  = 
           + __u2_k2_1_0 * __twopf_re_k2(0, 0)
           + __u2_k2_1_1 * __twopf_re_k2(1, 0);
        __dtwopf_re_k2(1, 1)  = 
           + __u2_k2_1_0 * __twopf_re_k2(0, 1)
           + __u2_k2_1_1 * __twopf_re_k2(1, 1);
        __dtwopf_re_k2(0, 0)  += 
           + __u2_k2_0_0 * __twopf_re_k2(0, 0)
           + __u2_k2_0_1 * __twopf_re_k2(0, 1);
        __dtwopf_re_k2(0, 1)  += 
           + __u2_k2_1_0 * __twopf_re_k2(0, 0)
           + __u2_k2_1_1 * __twopf_re_k2(0, 1);
        __dtwopf_re_k2(1, 0)  += 
           + __u2_k2_0_0 * __twopf_re_k2(1, 0)
           + __u2_k2_0_1 * __twopf_re_k2(1, 1);
        __dtwopf_re_k2(1, 1)  += 
           + __u2_k2_1_0 * __twopf_re_k2(1, 0)
           + __u2_k2_1_1 * __twopf_re_k2(1, 1);

        __dtwopf_im_k2(0, 0)  = 
           + __u2_k2_0_0 * __twopf_im_k2(0, 0)
           + __u2_k2_0_1 * __twopf_im_k2(1, 0);
        __dtwopf_im_k2(0, 1)  = 
           + __u2_k2_0_0 * __twopf_im_k2(0, 1)
           + __u2_k2_0_1 * __twopf_im_k2(1, 1);
        __dtwopf_im_k2(1, 0)  = 
           + __u2_k2_1_0 * __twopf_im_k2(0, 0)
           + __u2_k2_1_1 * __twopf_im_k2(1, 0);
        __dtwopf_im_k2(1, 1)  = 
           + __u2_k2_1_0 * __twopf_im_k2(0, 1)
           + __u2_k2_1_1 * __twopf_im_k2(1, 1);
        __dtwopf_im_k2(0, 0)  += 
           + __u2_k2_0_0 * __twopf_im_k2(0, 0)
           + __u2_k2_0_1 * __twopf_im_k2(0, 1);
        __dtwopf_im_k2(0, 1)  += 
           + __u2_k2_1_0 * __twopf_im_k2(0, 0)
           + __u2_k2_1_1 * __twopf_im_k2(0, 1);
        __dtwopf_im_k2(1, 0)  += 
           + __u2_k2_0_0 * __twopf_im_k2(1, 0)
           + __u2_k2_0_1 * __twopf_im_k2(1, 1);
        __dtwopf_im_k2(1, 1)  += 
           + __u2_k2_1_0 * __twopf_im_k2(1, 0)
           + __u2_k2_1_1 * __twopf_im_k2(1, 1);

        __dtwopf_re_k3(0, 0)  = 
           + __u2_k3_0_0 * __twopf_re_k3(0, 0)
           + __u2_k3_0_1 * __twopf_re_k3(1, 0);
        __dtwopf_re_k3(0, 1)  = 
           + __u2_k3_0_0 * __twopf_re_k3(0, 1)
           + __u2_k3_0_1 * __twopf_re_k3(1, 1);
        __dtwopf_re_k3(1, 0)  = 
           + __u2_k3_1_0 * __twopf_re_k3(0, 0)
           + __u2_k3_1_1 * __twopf_re_k3(1, 0);
        __dtwopf_re_k3(1, 1)  = 
           + __u2_k3_1_0 * __twopf_re_k3(0, 1)
           + __u2_k3_1_1 * __twopf_re_k3(1, 1);
        __dtwopf_re_k3(0, 0)  += 
           + __u2_k3_0_0 * __twopf_re_k3(0, 0)
           + __u2_k3_0_1 * __twopf_re_k3(0, 1);
        __dtwopf_re_k3(0, 1)  += 
           + __u2_k3_1_0 * __twopf_re_k3(0, 0)
           + __u2_k3_1_1 * __twopf_re_k3(0, 1);
        __dtwopf_re_k3(1, 0)  += 
           + __u2_k3_0_0 * __twopf_re_k3(1, 0)
           + __u2_k3_0_1 * __twopf_re_k3(1, 1);
        __dtwopf_re_k3(1, 1)  += 
           + __u2_k3_1_0 * __twopf_re_k3(1, 0)
           + __u2_k3_1_1 * __twopf_re_k3(1, 1);

        __dtwopf_im_k3(0, 0)  = 
           + __u2_k3_0_0 * __twopf_im_k3(0, 0)
           + __u2_k3_0_1 * __twopf_im_k3(1, 0);
        __dtwopf_im_k3(0, 1)  = 
           + __u2_k3_0_0 * __twopf_im_k3(0, 1)
           + __u2_k3_0_1 * __twopf_im_k3(1, 1);
        __dtwopf_im_k3(1, 0)  = 
           + __u2_k3_1_0 * __twopf_im_k3(0, 0)
           + __u2_k3_1_1 * __twopf_im_k3(1, 0);
        __dtwopf_im_k3(1, 1)  = 
           + __u2_k3_1_0 * __twopf_im_k3(0, 1)
           + __u2_k3_1_1 * __twopf_im_k3(1, 1);
        __dtwopf_im_k3(0, 0)  += 
           + __u2_k3_0_0 * __twopf_im_k3(0, 0)
           + __u2_k3_0_1 * __twopf_im_k3(0, 1);
        __dtwopf_im_k3(0, 1)  += 
           + __u2_k3_1_0 * __twopf_im_k3(0, 0)
           + __u2_k3_1_1 * __twopf_im_k3(0, 1);
        __dtwopf_im_k3(1, 0)  += 
           + __u2_k3_0_0 * __twopf_im_k3(1, 0)
           + __u2_k3_0_1 * __twopf_im_k3(1, 1);
        __dtwopf_im_k3(1, 1)  += 
           + __u2_k3_1_0 * __twopf_im_k3(1, 0)
           + __u2_k3_1_1 * __twopf_im_k3(1, 1);

        // evolve the components of the 3pf
        // index placement matters, partly because of the k-dependence
        // but also in the source terms from the imaginary components of the 2pf

        __dthreepf(0, 0, 0)  = 
           + __u2_k1_0_0 * __threepf(0, 0, 0)
           + __u2_k1_0_1 * __threepf(1, 0, 0);
        __dthreepf(0, 0, 1)  = 
           + __u2_k1_0_0 * __threepf(0, 0, 1)
           + __u2_k1_0_1 * __threepf(1, 0, 1);
        __dthreepf(0, 1, 0)  = 
           + __u2_k1_0_0 * __threepf(0, 1, 0)
           + __u2_k1_0_1 * __threepf(1, 1, 0);
        __dthreepf(0, 1, 1)  = 
           + __u2_k1_0_0 * __threepf(0, 1, 1)
           + __u2_k1_0_1 * __threepf(1, 1, 1);
        __dthreepf(1, 0, 0)  = 
           + __u2_k1_1_0 * __threepf(0, 0, 0)
           + __u2_k1_1_1 * __threepf(1, 0, 0);
        __dthreepf(1, 0, 1)  = 
           + __u2_k1_1_0 * __threepf(0, 0, 1)
           + __u2_k1_1_1 * __threepf(1, 0, 1);
        __dthreepf(1, 1, 0)  = 
           + __u2_k1_1_0 * __threepf(0, 1, 0)
           + __u2_k1_1_1 * __threepf(1, 1, 0);
        __dthreepf(1, 1, 1)  = 
           + __u2_k1_1_0 * __threepf(0, 1, 1)
           + __u2_k1_1_1 * __threepf(1, 1, 1);
        __dthreepf(0, 0, 0)  += 
           + __u3_k1k2k3_0_0_0 * __twopf_re_k2(0, 0) * __twopf_re_k3(0, 0)
           + __u3_k1k2k3_0_0_1 * __twopf_re_k2(0, 0) * __twopf_re_k3(1, 0)
           + __u3_k1k2k3_0_1_0 * __twopf_re_k2(1, 0) * __twopf_re_k3(0, 0)
           + __u3_k1k2k3_0_1_1 * __twopf_re_k2(1, 0) * __twopf_re_k3(1, 0);
        __dthreepf(0, 0, 1)  += 
           + __u3_k1k2k3_0_0_0 * __twopf_re_k2(0, 0) * __twopf_re_k3(0, 1)
           + __u3_k1k2k3_0_0_1 * __twopf_re_k2(0, 0) * __twopf_re_k3(1, 1)
           + __u3_k1k2k3_0_1_0 * __twopf_re_k2(1, 0) * __twopf_re_k3(0, 1)
           + __u3_k1k2k3_0_1_1 * __twopf_re_k2(1, 0) * __twopf_re_k3(1, 1);
        __dthreepf(0, 1, 0)  += 
           + __u3_k1k2k3_0_0_0 * __twopf_re_k2(0, 1) * __twopf_re_k3(0, 0)
           + __u3_k1k2k3_0_0_1 * __twopf_re_k2(0, 1) * __twopf_re_k3(1, 0)
           + __u3_k1k2k3_0_1_0 * __twopf_re_k2(1, 1) * __twopf_re_k3(0, 0)
           + __u3_k1k2k3_0_1_1 * __twopf_re_k2(1, 1) * __twopf_re_k3(1, 0);
        __dthreepf(0, 1, 1)  += 
           + __u3_k1k2k3_0_0_0 * __twopf_re_k2(0, 1) * __twopf_re_k3(0, 1)
           + __u3_k1k2k3_0_0_1 * __twopf_re_k2(0, 1) * __twopf_re_k3(1, 1)
           + __u3_k1k2k3_0_1_0 * __twopf_re_k2(1, 1) * __twopf_re_k3(0, 1)
           + __u3_k1k2k3_0_1_1 * __twopf_re_k2(1, 1) * __twopf_re_k3(1, 1);
        __dthreepf(1, 0, 0)  += 
           + __u3_k1k2k3_1_0_0 * __twopf_re_k2(0, 0) * __twopf_re_k3(0, 0)
           + __u3_k1k2k3_1_0_1 * __twopf_re_k2(0, 0) * __twopf_re_k3(1, 0)
           + __u3_k1k2k3_1_1_0 * __twopf_re_k2(1, 0) * __twopf_re_k3(0, 0)
           + __u3_k1k2k3_1_1_1 * __twopf_re_k2(1, 0) * __twopf_re_k3(1, 0);
        __dthreepf(1, 0, 1)  += 
           + __u3_k1k2k3_1_0_0 * __twopf_re_k2(0, 0) * __twopf_re_k3(0, 1)
           + __u3_k1k2k3_1_0_1 * __twopf_re_k2(0, 0) * __twopf_re_k3(1, 1)
           + __u3_k1k2k3_1_1_0 * __twopf_re_k2(1, 0) * __twopf_re_k3(0, 1)
           + __u3_k1k2k3_1_1_1 * __twopf_re_k2(1, 0) * __twopf_re_k3(1, 1);
        __dthreepf(1, 1, 0)  += 
           + __u3_k1k2k3_1_0_0 * __twopf_re_k2(0, 1) * __twopf_re_k3(0, 0)
           + __u3_k1k2k3_1_0_1 * __twopf_re_k2(0, 1) * __twopf_re_k3(1, 0)
           + __u3_k1k2k3_1_1_0 * __twopf_re_k2(1, 1) * __twopf_re_k3(0, 0)
           + __u3_k1k2k3_1_1_1 * __twopf_re_k2(1, 1) * __twopf_re_k3(1, 0);
        __dthreepf(1, 1, 1)  += 
           + __u3_k1k2k3_1_0_0 * __twopf_re_k2(0, 1) * __twopf_re_k3(0, 1)
           + __u3_k1k2k3_1_0_1 * __twopf_re_k2(0, 1) * __twopf_re_k3(1, 1)
           + __u3_k1k2k3_1_1_0 * __twopf_re_k2(1, 1) * __twopf_re_k3(0, 1)
           + __u3_k1k2k3_1_1_1 * __twopf_re_k2(1, 1) * __twopf_re_k3(1, 1);
        __dthreepf(0, 0, 0)  += 
           - __u3_k1k2k3_0_0_0 * __twopf_im_k2(0, 0) * __twopf_im_k3(0, 0)
           - __u3_k1k2k3_0_0_1 * __twopf_im_k2(0, 0) * __twopf_im_k3(1, 0)
           - __u3_k1k2k3_0_1_0 * __twopf_im_k2(1, 0) * __twopf_im_k3(0, 0)
           - __u3_k1k2k3_0_1_1 * __twopf_im_k2(1, 0) * __twopf_im_k3(1, 0);
        __dthreepf(0, 0, 1)  += 
           - __u3_k1k2k3_0_0_0 * __twopf_im_k2(0, 0) * __twopf_im_k3(0, 1)
           - __u3_k1k2k3_0_0_1 * __twopf_im_k2(0, 0) * __twopf_im_k3(1, 1)
           - __u3_k1k2k3_0_1_0 * __twopf_im_k2(1, 0) * __twopf_im_k3(0, 1)
           - __u3_k1k2k3_0_1_1 * __twopf_im_k2(1, 0) * __twopf_im_k3(1, 1);
        __dthreepf(0, 1, 0)  += 
           - __u3_k1k2k3_0_0_0 * __twopf_im_k2(0, 1) * __twopf_im_k3(0, 0)
           - __u3_k1k2k3_0_0_1 * __twopf_im_k2(0, 1) * __twopf_im_k3(1, 0)
           - __u3_k1k2k3_0_1_0 * __twopf_im_k2(1, 1) * __twopf_im_k3(0, 0)
           - __u3_k1k2k3_0_1_1 * __twopf_im_k2(1, 1) * __twopf_im_k3(1, 0);
        __dthreepf(0, 1, 1)  += 
           - __u3_k1k2k3_0_0_0 * __twopf_im_k2(0, 1) * __twopf_im_k3(0, 1)
           - __u3_k1k2k3_0_0_1 * __twopf_im_k2(0, 1) * __twopf_im_k3(1, 1)
           - __u3_k1k2k3_0_1_0 * __twopf_im_k2(1, 1) * __twopf_im_k3(0, 1)
           - __u3_k1k2k3_0_1_1 * __twopf_im_k2(1, 1) * __twopf_im_k3(1, 1);
        __dthreepf(1, 0, 0)  += 
           - __u3_k1k2k3_1_0_0 * __twopf_im_k2(0, 0) * __twopf_im_k3(0, 0)
           - __u3_k1k2k3_1_0_1 * __twopf_im_k2(0, 0) * __twopf_im_k3(1, 0)
           - __u3_k1k2k3_1_1_0 * __twopf_im_k2(1, 0) * __twopf_im_k3(0, 0)
           - __u3_k1k2k3_1_1_1 * __twopf_im_k2(1, 0) * __twopf_im_k3(1, 0);
        __dthreepf(1, 0, 1)  += 
           - __u3_k1k2k3_1_0_0 * __twopf_im_k2(0, 0) * __twopf_im_k3(0, 1)
           - __u3_k1k2k3_1_0_1 * __twopf_im_k2(0, 0) * __twopf_im_k3(1, 1)
           - __u3_k1k2k3_1_1_0 * __twopf_im_k2(1, 0) * __twopf_im_k3(0, 1)
           - __u3_k1k2k3_1_1_1 * __twopf_im_k2(1, 0) * __twopf_im_k3(1, 1);
        __dthreepf(1, 1, 0)  += 
           - __u3_k1k2k3_1_0_0 * __twopf_im_k2(0, 1) * __twopf_im_k3(0, 0)
           - __u3_k1k2k3_1_0_1 * __twopf_im_k2(0, 1) * __twopf_im_k3(1, 0)
           - __u3_k1k2k3_1_1_0 * __twopf_im_k2(1, 1) * __twopf_im_k3(0, 0)
           - __u3_k1k2k3_1_1_1 * __twopf_im_k2(1, 1) * __twopf_im_k3(1, 0);
        __dthreepf(1, 1, 1)  += 
           - __u3_k1k2k3_1_0_0 * __twopf_im_k2(0, 1) * __twopf_im_k3(0, 1)
           - __u3_k1k2k3_1_0_1 * __twopf_im_k2(0, 1) * __twopf_im_k3(1, 1)
           - __u3_k1k2k3_1_1_0 * __twopf_im_k2(1, 1) * __twopf_im_k3(0, 1)
           - __u3_k1k2k3_1_1_1 * __twopf_im_k2(1, 1) * __twopf_im_k3(1, 1);

        __dthreepf(0, 0, 0)  += 
           + __u2_k2_0_0 * __threepf(0, 0, 0)
           + __u2_k2_0_1 * __threepf(0, 1, 0);
        __dthreepf(0, 0, 1)  += 
           + __u2_k2_0_0 * __threepf(0, 0, 1)
           + __u2_k2_0_1 * __threepf(0, 1, 1);
        __dthreepf(0, 1, 0)  += 
           + __u2_k2_1_0 * __threepf(0, 0, 0)
           + __u2_k2_1_1 * __threepf(0, 1, 0);
        __dthreepf(0, 1, 1)  += 
           + __u2_k2_1_0 * __threepf(0, 0, 1)
           + __u2_k2_1_1 * __threepf(0, 1, 1);
        __dthreepf(1, 0, 0)  += 
           + __u2_k2_0_0 * __threepf(1, 0, 0)
           + __u2_k2_0_1 * __threepf(1, 1, 0);
        __dthreepf(1, 0, 1)  += 
           + __u2_k2_0_0 * __threepf(1, 0, 1)
           + __u2_k2_0_1 * __threepf(1, 1, 1);
        __dthreepf(1, 1, 0)  += 
           + __u2_k2_1_0 * __threepf(1, 0, 0)
           + __u2_k2_1_1 * __threepf(1, 1, 0);
        __dthreepf(1, 1, 1)  += 
           + __u2_k2_1_0 * __threepf(1, 0, 1)
           + __u2_k2_1_1 * __threepf(1, 1, 1);
        __dthreepf(0, 0, 0)  += 
           + __u3_k2k1k3_0_0_0 * __twopf_re_k1(0, 0) * __twopf_re_k3(0, 0)
           + __u3_k2k1k3_0_0_1 * __twopf_re_k1(0, 0) * __twopf_re_k3(1, 0)
           + __u3_k2k1k3_0_1_0 * __twopf_re_k1(0, 1) * __twopf_re_k3(0, 0)
           + __u3_k2k1k3_0_1_1 * __twopf_re_k1(0, 1) * __twopf_re_k3(1, 0);
        __dthreepf(0, 0, 1)  += 
           + __u3_k2k1k3_0_0_0 * __twopf_re_k1(0, 0) * __twopf_re_k3(0, 1)
           + __u3_k2k1k3_0_0_1 * __twopf_re_k1(0, 0) * __twopf_re_k3(1, 1)
           + __u3_k2k1k3_0_1_0 * __twopf_re_k1(0, 1) * __twopf_re_k3(0, 1)
           + __u3_k2k1k3_0_1_1 * __twopf_re_k1(0, 1) * __twopf_re_k3(1, 1);
        __dthreepf(0, 1, 0)  += 
           + __u3_k2k1k3_1_0_0 * __twopf_re_k1(0, 0) * __twopf_re_k3(0, 0)
           + __u3_k2k1k3_1_0_1 * __twopf_re_k1(0, 0) * __twopf_re_k3(1, 0)
           + __u3_k2k1k3_1_1_0 * __twopf_re_k1(0, 1) * __twopf_re_k3(0, 0)
           + __u3_k2k1k3_1_1_1 * __twopf_re_k1(0, 1) * __twopf_re_k3(1, 0);
        __dthreepf(0, 1, 1)  += 
           + __u3_k2k1k3_1_0_0 * __twopf_re_k1(0, 0) * __twopf_re_k3(0, 1)
           + __u3_k2k1k3_1_0_1 * __twopf_re_k1(0, 0) * __twopf_re_k3(1, 1)
           + __u3_k2k1k3_1_1_0 * __twopf_re_k1(0, 1) * __twopf_re_k3(0, 1)
           + __u3_k2k1k3_1_1_1 * __twopf_re_k1(0, 1) * __twopf_re_k3(1, 1);
        __dthreepf(1, 0, 0)  += 
           + __u3_k2k1k3_0_0_0 * __twopf_re_k1(1, 0) * __twopf_re_k3(0, 0)
           + __u3_k2k1k3_0_0_1 * __twopf_re_k1(1, 0) * __twopf_re_k3(1, 0)
           + __u3_k2k1k3_0_1_0 * __twopf_re_k1(1, 1) * __twopf_re_k3(0, 0)
           + __u3_k2k1k3_0_1_1 * __twopf_re_k1(1, 1) * __twopf_re_k3(1, 0);
        __dthreepf(1, 0, 1)  += 
           + __u3_k2k1k3_0_0_0 * __twopf_re_k1(1, 0) * __twopf_re_k3(0, 1)
           + __u3_k2k1k3_0_0_1 * __twopf_re_k1(1, 0) * __twopf_re_k3(1, 1)
           + __u3_k2k1k3_0_1_0 * __twopf_re_k1(1, 1) * __twopf_re_k3(0, 1)
           + __u3_k2k1k3_0_1_1 * __twopf_re_k1(1, 1) * __twopf_re_k3(1, 1);
        __dthreepf(1, 1, 0)  += 
           + __u3_k2k1k3_1_0_0 * __twopf_re_k1(1, 0) * __twopf_re_k3(0, 0)
           + __u3_k2k1k3_1_0_1 * __twopf_re_k1(1, 0) * __twopf_re_k3(1, 0)
           + __u3_k2k1k3_1_1_0 * __twopf_re_k1(1, 1) * __twopf_re_k3(0, 0)
           + __u3_k2k1k3_1_1_1 * __twopf_re_k1(1, 1) * __twopf_re_k3(1, 0);
        __dthreepf(1, 1, 1)  += 
           + __u3_k2k1k3_1_0_0 * __twopf_re_k1(1, 0) * __twopf_re_k3(0, 1)
           + __u3_k2k1k3_1_0_1 * __twopf_re_k1(1, 0) * __twopf_re_k3(1, 1)
           + __u3_k2k1k3_1_1_0 * __twopf_re_k1(1, 1) * __twopf_re_k3(0, 1)
           + __u3_k2k1k3_1_1_1 * __twopf_re_k1(1, 1) * __twopf_re_k3(1, 1);
        __dthreepf(0, 0, 0)  += 
           - __u3_k2k1k3_0_0_0 * __twopf_im_k1(0, 0) * __twopf_im_k3(0, 0)
           - __u3_k2k1k3_0_0_1 * __twopf_im_k1(0, 0) * __twopf_im_k3(1, 0)
           - __u3_k2k1k3_0_1_0 * __twopf_im_k1(0, 1) * __twopf_im_k3(0, 0)
           - __u3_k2k1k3_0_1_1 * __twopf_im_k1(0, 1) * __twopf_im_k3(1, 0);
        __dthreepf(0, 0, 1)  += 
           - __u3_k2k1k3_0_0_0 * __twopf_im_k1(0, 0) * __twopf_im_k3(0, 1)
           - __u3_k2k1k3_0_0_1 * __twopf_im_k1(0, 0) * __twopf_im_k3(1, 1)
           - __u3_k2k1k3_0_1_0 * __twopf_im_k1(0, 1) * __twopf_im_k3(0, 1)
           - __u3_k2k1k3_0_1_1 * __twopf_im_k1(0, 1) * __twopf_im_k3(1, 1);
        __dthreepf(0, 1, 0)  += 
           - __u3_k2k1k3_1_0_0 * __twopf_im_k1(0, 0) * __twopf_im_k3(0, 0)
           - __u3_k2k1k3_1_0_1 * __twopf_im_k1(0, 0) * __twopf_im_k3(1, 0)
           - __u3_k2k1k3_1_1_0 * __twopf_im_k1(0, 1) * __twopf_im_k3(0, 0)
           - __u3_k2k1k3_1_1_1 * __twopf_im_k1(0, 1) * __twopf_im_k3(1, 0);
        __dthreepf(0, 1, 1)  += 
           - __u3_k2k1k3_1_0_0 * __twopf_im_k1(0, 0) * __twopf_im_k3(0, 1)
           - __u3_k2k1k3_1_0_1 * __twopf_im_k1(0, 0) * __twopf_im_k3(1, 1)
           - __u3_k2k1k3_1_1_0 * __twopf_im_k1(0, 1) * __twopf_im_k3(0, 1)
           - __u3_k2k1k3_1_1_1 * __twopf_im_k1(0, 1) * __twopf_im_k3(1, 1);
        __dthreepf(1, 0, 0)  += 
           - __u3_k2k1k3_0_0_0 * __twopf_im_k1(1, 0) * __twopf_im_k3(0, 0)
           - __u3_k2k1k3_0_0_1 * __twopf_im_k1(1, 0) * __twopf_im_k3(1, 0)
           - __u3_k2k1k3_0_1_0 * __twopf_im_k1(1, 1) * __twopf_im_k3(0, 0)
           - __u3_k2k1k3_0_1_1 * __twopf_im_k1(1, 1) * __twopf_im_k3(1, 0);
        __dthreepf(1, 0, 1)  += 
           - __u3_k2k1k3_0_0_0 * __twopf_im_k1(1, 0) * __twopf_im_k3(0, 1)
           - __u3_k2k1k3_0_0_1 * __twopf_im_k1(1, 0) * __twopf_im_k3(1, 1)
           - __u3_k2k1k3_0_1_0 * __twopf_im_k1(1, 1) * __twopf_im_k3(0, 1)
           - __u3_k2k1k3_0_1_1 * __twopf_im_k1(1, 1) * __twopf_im_k3(1, 1);
        __dthreepf(1, 1, 0)  += 
           - __u3_k2k1k3_1_0_0 * __twopf_im_k1(1, 0) * __twopf_im_k3(0, 0)
           - __u3_k2k1k3_1_0_1 * __twopf_im_k1(1, 0) * __twopf_im_k3(1, 0)
           - __u3_k2k1k3_1_1_0 * __twopf_im_k1(1, 1) * __twopf_im_k3(0, 0)
           - __u3_k2k1k3_1_1_1 * __twopf_im_k1(1, 1) * __twopf_im_k3(1, 0);
        __dthreepf(1, 1, 1)  += 
           - __u3_k2k1k3_1_0_0 * __twopf_im_k1(1, 0) * __twopf_im_k3(0, 1)
           - __u3_k2k1k3_1_0_1 * __twopf_im_k1(1, 0) * __twopf_im_k3(1, 1)
           - __u3_k2k1k3_1_1_0 * __twopf_im_k1(1, 1) * __twopf_im_k3(0, 1)
           - __u3_k2k1k3_1_1_1 * __twopf_im_k1(1, 1) * __twopf_im_k3(1, 1);

        __dthreepf(0, 0, 0)  += 
           + __u2_k3_0_0 * __threepf(0, 0, 0)
           + __u2_k3_0_1 * __threepf(0, 0, 1);
        __dthreepf(0, 0, 1)  += 
           + __u2_k3_1_0 * __threepf(0, 0, 0)
           + __u2_k3_1_1 * __threepf(0, 0, 1);
        __dthreepf(0, 1, 0)  += 
           + __u2_k3_0_0 * __threepf(0, 1, 0)
           + __u2_k3_0_1 * __threepf(0, 1, 1);
        __dthreepf(0, 1, 1)  += 
           + __u2_k3_1_0 * __threepf(0, 1, 0)
           + __u2_k3_1_1 * __threepf(0, 1, 1);
        __dthreepf(1, 0, 0)  += 
           + __u2_k3_0_0 * __threepf(1, 0, 0)
           + __u2_k3_0_1 * __threepf(1, 0, 1);
        __dthreepf(1, 0, 1)  += 
           + __u2_k3_1_0 * __threepf(1, 0, 0)
           + __u2_k3_1_1 * __threepf(1, 0, 1);
        __dthreepf(1, 1, 0)  += 
           + __u2_k3_0_0 * __threepf(1, 1, 0)
           + __u2_k3_0_1 * __threepf(1, 1, 1);
        __dthreepf(1, 1, 1)  += 
           + __u2_k3_1_0 * __threepf(1, 1, 0)
           + __u2_k3_1_1 * __threepf(1, 1, 1);
        __dthreepf(0, 0, 0)  += 
           + __u3_k3k1k2_0_0_0 * __twopf_re_k1(0, 0) * __twopf_re_k2(0, 0)
           + __u3_k3k1k2_0_0_1 * __twopf_re_k1(0, 0) * __twopf_re_k2(0, 1)
           + __u3_k3k1k2_0_1_0 * __twopf_re_k1(0, 1) * __twopf_re_k2(0, 0)
           + __u3_k3k1k2_0_1_1 * __twopf_re_k1(0, 1) * __twopf_re_k2(0, 1);
        __dthreepf(0, 0, 1)  += 
           + __u3_k3k1k2_1_0_0 * __twopf_re_k1(0, 0) * __twopf_re_k2(0, 0)
           + __u3_k3k1k2_1_0_1 * __twopf_re_k1(0, 0) * __twopf_re_k2(0, 1)
           + __u3_k3k1k2_1_1_0 * __twopf_re_k1(0, 1) * __twopf_re_k2(0, 0)
           + __u3_k3k1k2_1_1_1 * __twopf_re_k1(0, 1) * __twopf_re_k2(0, 1);
        __dthreepf(0, 1, 0)  += 
           + __u3_k3k1k2_0_0_0 * __twopf_re_k1(0, 0) * __twopf_re_k2(1, 0)
           + __u3_k3k1k2_0_0_1 * __twopf_re_k1(0, 0) * __twopf_re_k2(1, 1)
           + __u3_k3k1k2_0_1_0 * __twopf_re_k1(0, 1) * __twopf_re_k2(1, 0)
           + __u3_k3k1k2_0_1_1 * __twopf_re_k1(0, 1) * __twopf_re_k2(1, 1);
        __dthreepf(0, 1, 1)  += 
           + __u3_k3k1k2_1_0_0 * __twopf_re_k1(0, 0) * __twopf_re_k2(1, 0)
           + __u3_k3k1k2_1_0_1 * __twopf_re_k1(0, 0) * __twopf_re_k2(1, 1)
           + __u3_k3k1k2_1_1_0 * __twopf_re_k1(0, 1) * __twopf_re_k2(1, 0)
           + __u3_k3k1k2_1_1_1 * __twopf_re_k1(0, 1) * __twopf_re_k2(1, 1);
        __dthreepf(1, 0, 0)  += 
           + __u3_k3k1k2_0_0_0 * __twopf_re_k1(1, 0) * __twopf_re_k2(0, 0)
           + __u3_k3k1k2_0_0_1 * __twopf_re_k1(1, 0) * __twopf_re_k2(0, 1)
           + __u3_k3k1k2_0_1_0 * __twopf_re_k1(1, 1) * __twopf_re_k2(0, 0)
           + __u3_k3k1k2_0_1_1 * __twopf_re_k1(1, 1) * __twopf_re_k2(0, 1);
        __dthreepf(1, 0, 1)  += 
           + __u3_k3k1k2_1_0_0 * __twopf_re_k1(1, 0) * __twopf_re_k2(0, 0)
           + __u3_k3k1k2_1_0_1 * __twopf_re_k1(1, 0) * __twopf_re_k2(0, 1)
           + __u3_k3k1k2_1_1_0 * __twopf_re_k1(1, 1) * __twopf_re_k2(0, 0)
           + __u3_k3k1k2_1_1_1 * __twopf_re_k1(1, 1) * __twopf_re_k2(0, 1);
        __dthreepf(1, 1, 0)  += 
           + __u3_k3k1k2_0_0_0 * __twopf_re_k1(1, 0) * __twopf_re_k2(1, 0)
           + __u3_k3k1k2_0_0_1 * __twopf_re_k1(1, 0) * __twopf_re_k2(1, 1)
           + __u3_k3k1k2_0_1_0 * __twopf_re_k1(1, 1) * __twopf_re_k2(1, 0)
           + __u3_k3k1k2_0_1_1 * __twopf_re_k1(1, 1) * __twopf_re_k2(1, 1);
        __dthreepf(1, 1, 1)  += 
           + __u3_k3k1k2_1_0_0 * __twopf_re_k1(1, 0) * __twopf_re_k2(1, 0)
           + __u3_k3k1k2_1_0_1 * __twopf_re_k1(1, 0) * __twopf_re_k2(1, 1)
           + __u3_k3k1k2_1_1_0 * __twopf_re_k1(1, 1) * __twopf_re_k2(1, 0)
           + __u3_k3k1k2_1_1_1 * __twopf_re_k1(1, 1) * __twopf_re_k2(1, 1);
        __dthreepf(0, 0, 0)  += 
           - __u3_k3k1k2_0_0_0 * __twopf_im_k1(0, 0) * __twopf_im_k2(0, 0)
           - __u3_k3k1k2_0_0_1 * __twopf_im_k1(0, 0) * __twopf_im_k2(0, 1)
           - __u3_k3k1k2_0_1_0 * __twopf_im_k1(0, 1) * __twopf_im_k2(0, 0)
           - __u3_k3k1k2_0_1_1 * __twopf_im_k1(0, 1) * __twopf_im_k2(0, 1);
        __dthreepf(0, 0, 1)  += 
           - __u3_k3k1k2_1_0_0 * __twopf_im_k1(0, 0) * __twopf_im_k2(0, 0)
           - __u3_k3k1k2_1_0_1 * __twopf_im_k1(0, 0) * __twopf_im_k2(0, 1)
           - __u3_k3k1k2_1_1_0 * __twopf_im_k1(0, 1) * __twopf_im_k2(0, 0)
           - __u3_k3k1k2_1_1_1 * __twopf_im_k1(0, 1) * __twopf_im_k2(0, 1);
        __dthreepf(0, 1, 0)  += 
           - __u3_k3k1k2_0_0_0 * __twopf_im_k1(0, 0) * __twopf_im_k2(1, 0)
           - __u3_k3k1k2_0_0_1 * __twopf_im_k1(0, 0) * __twopf_im_k2(1, 1)
           - __u3_k3k1k2_0_1_0 * __twopf_im_k1(0, 1) * __twopf_im_k2(1, 0)
           - __u3_k3k1k2_0_1_1 * __twopf_im_k1(0, 1) * __twopf_im_k2(1, 1);
        __dthreepf(0, 1, 1)  += 
           - __u3_k3k1k2_1_0_0 * __twopf_im_k1(0, 0) * __twopf_im_k2(1, 0)
           - __u3_k3k1k2_1_0_1 * __twopf_im_k1(0, 0) * __twopf_im_k2(1, 1)
           - __u3_k3k1k2_1_1_0 * __twopf_im_k1(0, 1) * __twopf_im_k2(1, 0)
           - __u3_k3k1k2_1_1_1 * __twopf_im_k1(0, 1) * __twopf_im_k2(1, 1);
        __dthreepf(1, 0, 0)  += 
           - __u3_k3k1k2_0_0_0 * __twopf_im_k1(1, 0) * __twopf_im_k2(0, 0)
           - __u3_k3k1k2_0_0_1 * __twopf_im_k1(1, 0) * __twopf_im_k2(0, 1)
           - __u3_k3k1k2_0_1_0 * __twopf_im_k1(1, 1) * __twopf_im_k2(0, 0)
           - __u3_k3k1k2_0_1_1 * __twopf_im_k1(1, 1) * __twopf_im_k2(0, 1);
        __dthreepf(1, 0, 1)  += 
           - __u3_k3k1k2_1_0_0 * __twopf_im_k1(1, 0) * __twopf_im_k2(0, 0)
           - __u3_k3k1k2_1_0_1 * __twopf_im_k1(1, 0) * __twopf_im_k2(0, 1)
           - __u3_k3k1k2_1_1_0 * __twopf_im_k1(1, 1) * __twopf_im_k2(0, 0)
           - __u3_k3k1k2_1_1_1 * __twopf_im_k1(1, 1) * __twopf_im_k2(0, 1);
        __dthreepf(1, 1, 0)  += 
           - __u3_k3k1k2_0_0_0 * __twopf_im_k1(1, 0) * __twopf_im_k2(1, 0)
           - __u3_k3k1k2_0_0_1 * __twopf_im_k1(1, 0) * __twopf_im_k2(1, 1)
           - __u3_k3k1k2_0_1_0 * __twopf_im_k1(1, 1) * __twopf_im_k2(1, 0)
           - __u3_k3k1k2_0_1_1 * __twopf_im_k1(1, 1) * __twopf_im_k2(1, 1);
        __dthreepf(1, 1, 1)  += 
           - __u3_k3k1k2_1_0_0 * __twopf_im_k1(1, 0) * __twopf_im_k2(1, 0)
           - __u3_k3k1k2_1_0_1 * __twopf_im_k1(1, 0) * __twopf_im_k2(1, 1)
           - __u3_k3k1k2_1_1_0 * __twopf_im_k1(1, 1) * __twopf_im_k2(1, 0)
           - __u3_k3k1k2_1_1_1 * __twopf_im_k1(1, 1) * __twopf_im_k2(1, 1);

#ifdef CPPTRANSPORT_INSTRUMENT
        __transport_eq_timer.stop();
        ++__invokations;
#endif
      }


    // IMPLEMENTATION - FUNCTOR FOR 3PF OBSERVATION


    template <typename Model, typename BatchObject>
    void quartic_mpi_threepf_observer<Model, BatchObject>::operator()(const threepf_state& x, number t)
      {
        DEFINE_INDEX_TOOLS

#undef __background
#undef __twopf_k1_tensor
#undef __twopf_k2_tensor
#undef __twopf_k3_tensor
#undef __twopf_re_k1
#undef __twopf_im_k1
#undef __twopf_re_k2
#undef __twopf_im_k2
#undef __twopf_re_k3
#undef __twopf_im_k3
#undef __threepf

#define __background(a)        x[quartic_pool::backg_start       + FLATTEN(a)]
#define __twopf_k1_tensor(a,b) x[quartic_pool::tensor_k1_start   + TENSOR_FLATTEN(a,b)]
#define __twopf_k2_tensor(a,b) x[quartic_pool::tensor_k2_start   + TENSOR_FLATTEN(a,b)]
#define __twopf_k3_tensor(a,b) x[quartic_pool::tensor_k3_start   + TENSOR_FLATTEN(a,b)]
#define __twopf_re_k1(a,b)     x[quartic_pool::twopf_re_k1_start + FLATTEN(a,b)]
#define __twopf_im_k1(a,b)     x[quartic_pool::twopf_im_k1_start + FLATTEN(a,b)]
#define __twopf_re_k2(a,b)     x[quartic_pool::twopf_re_k2_start + FLATTEN(a,b)]
#define __twopf_im_k2(a,b)     x[quartic_pool::twopf_im_k2_start + FLATTEN(a,b)]
#define __twopf_re_k3(a,b)     x[quartic_pool::twopf_re_k3_start + FLATTEN(a,b)]
#define __twopf_im_k3(a,b)     x[quartic_pool::twopf_im_k3_start + FLATTEN(a,b)]
#define __threepf(a,b,c)       x[quartic_pool::threepf_start     + FLATTEN(a,b,c)]

#ifndef CPPTRANSPORT_NO_STRICT_FP_TEST
        if(std::isnan(__background(0)) || std::isinf(__background(0))) throw runtime_exception(exception_type::INTEGRATION_FAILURE, CPPTRANSPORT_INTEGRATOR_NAN_OR_INF);
        if(std::isnan(__background(1)) || std::isinf(__background(1))) throw runtime_exception(exception_type::INTEGRATION_FAILURE, CPPTRANSPORT_INTEGRATOR_NAN_OR_INF);

        if(std::isnan(__twopf_k1_tensor(0,0)) || std::isinf(__twopf_k1_tensor(0,0))) throw runtime_exception(exception_type::INTEGRATION_FAILURE, CPPTRANSPORT_INTEGRATOR_NAN_OR_INF);
        if(std::isnan(__twopf_k1_tensor(0,1)) || std::isinf(__twopf_k1_tensor(0,1))) throw runtime_exception(exception_type::INTEGRATION_FAILURE, CPPTRANSPORT_INTEGRATOR_NAN_OR_INF);
        if(std::isnan(__twopf_k1_tensor(1,0)) || std::isinf(__twopf_k1_tensor(1,0))) throw runtime_exception(exception_type::INTEGRATION_FAILURE, CPPTRANSPORT_INTEGRATOR_NAN_OR_INF);
        if(std::isnan(__twopf_k1_tensor(1,1)) || std::isinf(__twopf_k1_tensor(1,1))) throw runtime_exception(exception_type::INTEGRATION_FAILURE, CPPTRANSPORT_INTEGRATOR_NAN_OR_INF);

        if(std::isnan(__twopf_k2_tensor(0,0)) || std::isinf(__twopf_k2_tensor(0,0))) throw runtime_exception(exception_type::INTEGRATION_FAILURE, CPPTRANSPORT_INTEGRATOR_NAN_OR_INF);
        if(std::isnan(__twopf_k2_tensor(0,1)) || std::isinf(__twopf_k2_tensor(0,1))) throw runtime_exception(exception_type::INTEGRATION_FAILURE, CPPTRANSPORT_INTEGRATOR_NAN_OR_INF);
        if(std::isnan(__twopf_k2_tensor(1,0)) || std::isinf(__twopf_k2_tensor(1,0))) throw runtime_exception(exception_type::INTEGRATION_FAILURE, CPPTRANSPORT_INTEGRATOR_NAN_OR_INF);
        if(std::isnan(__twopf_k2_tensor(1,1)) || std::isinf(__twopf_k2_tensor(1,1))) throw runtime_exception(exception_type::INTEGRATION_FAILURE, CPPTRANSPORT_INTEGRATOR_NAN_OR_INF);

        if(std::isnan(__twopf_k3_tensor(0,0)) || std::isinf(__twopf_k3_tensor(0,0))) throw runtime_exception(exception_type::INTEGRATION_FAILURE, CPPTRANSPORT_INTEGRATOR_NAN_OR_INF);
        if(std::isnan(__twopf_k3_tensor(0,1)) || std::isinf(__twopf_k3_tensor(0,1))) throw runtime_exception(exception_type::INTEGRATION_FAILURE, CPPTRANSPORT_INTEGRATOR_NAN_OR_INF);
        if(std::isnan(__twopf_k3_tensor(1,0)) || std::isinf(__twopf_k3_tensor(1,0))) throw runtime_exception(exception_type::INTEGRATION_FAILURE, CPPTRANSPORT_INTEGRATOR_NAN_OR_INF);
        if(std::isnan(__twopf_k3_tensor(1,1)) || std::isinf(__twopf_k3_tensor(1,1))) throw runtime_exception(exception_type::INTEGRATION_FAILURE, CPPTRANSPORT_INTEGRATOR_NAN_OR_INF);

        if(std::isnan(__twopf_re_k1(0, 0)) || std::isinf(__twopf_re_k1(0, 0))) throw runtime_exception(exception_type::INTEGRATION_FAILURE, CPPTRANSPORT_INTEGRATOR_NAN_OR_INF);
        if(std::isnan(__twopf_re_k1(0, 1)) || std::isinf(__twopf_re_k1(0, 1))) throw runtime_exception(exception_type::INTEGRATION_FAILURE, CPPTRANSPORT_INTEGRATOR_NAN_OR_INF);
        if(std::isnan(__twopf_re_k1(1, 0)) || std::isinf(__twopf_re_k1(1, 0))) throw runtime_exception(exception_type::INTEGRATION_FAILURE, CPPTRANSPORT_INTEGRATOR_NAN_OR_INF);
        if(std::isnan(__twopf_re_k1(1, 1)) || std::isinf(__twopf_re_k1(1, 1))) throw runtime_exception(exception_type::INTEGRATION_FAILURE, CPPTRANSPORT_INTEGRATOR_NAN_OR_INF);
        if(std::isnan(__twopf_im_k1(0, 0)) || std::isinf(__twopf_im_k1(0, 0))) throw runtime_exception(exception_type::INTEGRATION_FAILURE, CPPTRANSPORT_INTEGRATOR_NAN_OR_INF);
        if(std::isnan(__twopf_im_k1(0, 1)) || std::isinf(__twopf_im_k1(0, 1))) throw runtime_exception(exception_type::INTEGRATION_FAILURE, CPPTRANSPORT_INTEGRATOR_NAN_OR_INF);
        if(std::isnan(__twopf_im_k1(1, 0)) || std::isinf(__twopf_im_k1(1, 0))) throw runtime_exception(exception_type::INTEGRATION_FAILURE, CPPTRANSPORT_INTEGRATOR_NAN_OR_INF);
        if(std::isnan(__twopf_im_k1(1, 1)) || std::isinf(__twopf_im_k1(1, 1))) throw runtime_exception(exception_type::INTEGRATION_FAILURE, CPPTRANSPORT_INTEGRATOR_NAN_OR_INF);

        if(std::isnan(__twopf_re_k2(0, 0)) || std::isinf(__twopf_re_k2(0, 0))) throw runtime_exception(exception_type::INTEGRATION_FAILURE, CPPTRANSPORT_INTEGRATOR_NAN_OR_INF);
        if(std::isnan(__twopf_re_k2(0, 1)) || std::isinf(__twopf_re_k2(0, 1))) throw runtime_exception(exception_type::INTEGRATION_FAILURE, CPPTRANSPORT_INTEGRATOR_NAN_OR_INF);
        if(std::isnan(__twopf_re_k2(1, 0)) || std::isinf(__twopf_re_k2(1, 0))) throw runtime_exception(exception_type::INTEGRATION_FAILURE, CPPTRANSPORT_INTEGRATOR_NAN_OR_INF);
        if(std::isnan(__twopf_re_k2(1, 1)) || std::isinf(__twopf_re_k2(1, 1))) throw runtime_exception(exception_type::INTEGRATION_FAILURE, CPPTRANSPORT_INTEGRATOR_NAN_OR_INF);
        if(std::isnan(__twopf_im_k2(0, 0)) || std::isinf(__twopf_im_k2(0, 0))) throw runtime_exception(exception_type::INTEGRATION_FAILURE, CPPTRANSPORT_INTEGRATOR_NAN_OR_INF);
        if(std::isnan(__twopf_im_k2(0, 1)) || std::isinf(__twopf_im_k2(0, 1))) throw runtime_exception(exception_type::INTEGRATION_FAILURE, CPPTRANSPORT_INTEGRATOR_NAN_OR_INF);
        if(std::isnan(__twopf_im_k2(1, 0)) || std::isinf(__twopf_im_k2(1, 0))) throw runtime_exception(exception_type::INTEGRATION_FAILURE, CPPTRANSPORT_INTEGRATOR_NAN_OR_INF);
        if(std::isnan(__twopf_im_k2(1, 1)) || std::isinf(__twopf_im_k2(1, 1))) throw runtime_exception(exception_type::INTEGRATION_FAILURE, CPPTRANSPORT_INTEGRATOR_NAN_OR_INF);

        if(std::isnan(__twopf_re_k3(0, 0)) || std::isinf(__twopf_re_k3(0, 0))) throw runtime_exception(exception_type::INTEGRATION_FAILURE, CPPTRANSPORT_INTEGRATOR_NAN_OR_INF);
        if(std::isnan(__twopf_re_k3(0, 1)) || std::isinf(__twopf_re_k3(0, 1))) throw runtime_exception(exception_type::INTEGRATION_FAILURE, CPPTRANSPORT_INTEGRATOR_NAN_OR_INF);
        if(std::isnan(__twopf_re_k3(1, 0)) || std::isinf(__twopf_re_k3(1, 0))) throw runtime_exception(exception_type::INTEGRATION_FAILURE, CPPTRANSPORT_INTEGRATOR_NAN_OR_INF);
        if(std::isnan(__twopf_re_k3(1, 1)) || std::isinf(__twopf_re_k3(1, 1))) throw runtime_exception(exception_type::INTEGRATION_FAILURE, CPPTRANSPORT_INTEGRATOR_NAN_OR_INF);
        if(std::isnan(__twopf_im_k3(0, 0)) || std::isinf(__twopf_im_k3(0, 0))) throw runtime_exception(exception_type::INTEGRATION_FAILURE, CPPTRANSPORT_INTEGRATOR_NAN_OR_INF);
        if(std::isnan(__twopf_im_k3(0, 1)) || std::isinf(__twopf_im_k3(0, 1))) throw runtime_exception(exception_type::INTEGRATION_FAILURE, CPPTRANSPORT_INTEGRATOR_NAN_OR_INF);
        if(std::isnan(__twopf_im_k3(1, 0)) || std::isinf(__twopf_im_k3(1, 0))) throw runtime_exception(exception_type::INTEGRATION_FAILURE, CPPTRANSPORT_INTEGRATOR_NAN_OR_INF);
        if(std::isnan(__twopf_im_k3(1, 1)) || std::isinf(__twopf_im_k3(1, 1))) throw runtime_exception(exception_type::INTEGRATION_FAILURE, CPPTRANSPORT_INTEGRATOR_NAN_OR_INF);

        if(std::isnan(__threepf(0, 0, 0)) || std::isinf(__threepf(0, 0, 0))) throw runtime_exception(exception_type::INTEGRATION_FAILURE, CPPTRANSPORT_INTEGRATOR_NAN_OR_INF);
        if(std::isnan(__threepf(0, 0, 1)) || std::isinf(__threepf(0, 0, 1))) throw runtime_exception(exception_type::INTEGRATION_FAILURE, CPPTRANSPORT_INTEGRATOR_NAN_OR_INF);
        if(std::isnan(__threepf(0, 1, 0)) || std::isinf(__threepf(0, 1, 0))) throw runtime_exception(exception_type::INTEGRATION_FAILURE, CPPTRANSPORT_INTEGRATOR_NAN_OR_INF);
        if(std::isnan(__threepf(0, 1, 1)) || std::isinf(__threepf(0, 1, 1))) throw runtime_exception(exception_type::INTEGRATION_FAILURE, CPPTRANSPORT_INTEGRATOR_NAN_OR_INF);
        if(std::isnan(__threepf(1, 0, 0)) || std::isinf(__threepf(1, 0, 0))) throw runtime_exception(exception_type::INTEGRATION_FAILURE, CPPTRANSPORT_INTEGRATOR_NAN_OR_INF);
        if(std::isnan(__threepf(1, 0, 1)) || std::isinf(__threepf(1, 0, 1))) throw runtime_exception(exception_type::INTEGRATION_FAILURE, CPPTRANSPORT_INTEGRATOR_NAN_OR_INF);
        if(std::isnan(__threepf(1, 1, 0)) || std::isinf(__threepf(1, 1, 0))) throw runtime_exception(exception_type::INTEGRATION_FAILURE, CPPTRANSPORT_INTEGRATOR_NAN_OR_INF);
        if(std::isnan(__threepf(1, 1, 1)) || std::isinf(__threepf(1, 1, 1))) throw runtime_exception(exception_type::INTEGRATION_FAILURE, CPPTRANSPORT_INTEGRATOR_NAN_OR_INF);
#endif

        this->start_batching(static_cast<double>(t), this->get_log(), BatchObject::log_severity_level::normal);
        this->push(x);
        this->stop_batching();
      }


    }   // namespace transport


#endif  // CPPTRANSPORT_QUARTIC_MPI_H

