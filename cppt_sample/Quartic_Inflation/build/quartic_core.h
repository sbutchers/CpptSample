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
// processed on 2019-Apr-09 14:34:10

#ifndef CPPTRANSPORT_QUARTIC_CORE_H   // avoid multiple inclusion
#define CPPTRANSPORT_QUARTIC_CORE_H

#include <assert.h>
#include <cmath>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <memory>

#include "boost/numeric/odeint.hpp"
#include "boost/range/algorithm.hpp"
#include "boost/optional.hpp"

#include "Eigen/Core"

#include "transport-runtime/transport.h"
#include "transport-runtime/models/canonical_model.h"
#include "transport-runtime/models/odeint_defaults.h"


// #define CPPTRANSPORT_INSTRUMENT
// #define CPPTRANSPORT_NO_STRICT_FP_TEST


namespace transport
  {

    // Literal data pool
    namespace quartic_pool
      {
        static std::vector<std::string> field_names = { "phi" };
        static std::vector<std::string> latex_names = { "\\phi" };
        static std::vector<std::string> param_names = { "lambda" };
        static std::vector<std::string> platx_names = { "\\lambda" };
        static std::vector<std::string> state_names = { "phi", "__dphi" };

        static std::string              name        = "Quartic inflation";
        static std::string              citeguide   = "CppTransport user guide arXiv:16xx.yyyy";
        static std::string              description = "Quartic model with V=m^2phi^2/2";
        static std::string              license     = "CC BY";
        static unsigned int             revision    = 1;

        static std::vector<std::string> references  = {};
        static std::vector<std::string> urls        = { "http://transportmethod.com" };

        static author_db::value_type    auth_arr[]  = { { "Sean Butchers", std::make_unique<author_record>("Sean Butchers", "smlb20@sussex.ac.uk", "Astronomy Centre, University of Sussex") } };
        static author_db                author(std::make_move_iterator(std::begin(auth_arr)), std::make_move_iterator(std::end(auth_arr)));

        static std::string              unique_id   = "b5776652-458d-8ed9-f6fb-7efdec3514a3";

        constexpr unsigned int backg_size         = (2*1);
        constexpr unsigned int twopf_size         = ((2*1)*(2*1));
        constexpr unsigned int tensor_size        = (4);
        constexpr unsigned int threepf_size       = ((2*1)*(2*1)*(2*1));

        constexpr unsigned int backg_start        = 0;
        constexpr unsigned int tensor_start       = backg_start + backg_size;         // for twopf state vector
        constexpr unsigned int tensor_k1_start    = tensor_start;                     // for threepf state vector
        constexpr unsigned int tensor_k2_start    = tensor_k1_start + tensor_size;
        constexpr unsigned int tensor_k3_start    = tensor_k2_start + tensor_size;
        constexpr unsigned int twopf_start        = tensor_k1_start + tensor_size;    // for twopf state vector
        constexpr unsigned int twopf_re_k1_start  = tensor_k3_start + tensor_size;    // for threepf state vector
        constexpr unsigned int twopf_im_k1_start  = twopf_re_k1_start + twopf_size;
        constexpr unsigned int twopf_re_k2_start  = twopf_im_k1_start + twopf_size;
        constexpr unsigned int twopf_im_k2_start  = twopf_re_k2_start + twopf_size;
        constexpr unsigned int twopf_re_k3_start  = twopf_im_k2_start + twopf_size;
        constexpr unsigned int twopf_im_k3_start  = twopf_re_k3_start + twopf_size;
        constexpr unsigned int threepf_start      = twopf_im_k3_start + twopf_size;

        constexpr unsigned int backg_state_size   = backg_size;
        constexpr unsigned int twopf_state_size   = backg_size + tensor_size + twopf_size;
        constexpr unsigned int threepf_state_size = backg_size + 3*tensor_size + 6*twopf_size + threepf_size;

        constexpr unsigned int u2_size            = ((2*1)*(2*1));
        constexpr unsigned int u3_size            = ((2*1)*(2*1)*(2*1));


        // FLATTENING FUNCTIONS
        constexpr unsigned int SPECIES       (unsigned int a)                                 { return flatten_impl::species(a, 1); }
        constexpr unsigned int MOMENTUM      (unsigned int a)                                 { return flatten_impl::momentum(a, 1); }
        constexpr unsigned int IS_FIELD      (unsigned int a)                                 { return flatten_impl::is_field(a, 1); }
        constexpr unsigned int IS_MOMENTUM   (unsigned int a)                                 { return flatten_impl::is_momentum(a, 1); }

        constexpr unsigned int FLATTEN       (unsigned int a)                                 { return flatten_impl::flatten(a, 1); }
        constexpr unsigned int FLATTEN       (unsigned int a, unsigned int b)                 { return flatten_impl::flatten(a, b, 1); }
        constexpr unsigned int FLATTEN       (unsigned int a, unsigned int b, unsigned int c) { return flatten_impl::flatten(a, b, c, 1); }

        constexpr unsigned int FIELDS_FLATTEN(unsigned int a)                                 { return flatten_impl::fields_flatten(a, 1); }
        constexpr unsigned int FIELDS_FLATTEN(unsigned int a, unsigned int b)                 { return flatten_impl::fields_flatten(a, b, 1); }
        constexpr unsigned int FIELDS_FLATTEN(unsigned int a, unsigned int b, unsigned int c) { return flatten_impl::fields_flatten(a, b, c, 1); }

        constexpr unsigned int TENSOR_FLATTEN(unsigned int a, unsigned int b)                 { return flatten_impl::tensor_flatten(a, b); }

      }


#undef DEFINE_INDEX_TOOLS
#define DEFINE_INDEX_TOOLS \
    using quartic_pool::SPECIES; \
    using quartic_pool::MOMENTUM; \
    using quartic_pool::IS_FIELD; \
    using quartic_pool::IS_MOMENTUM; \
    using quartic_pool::FLATTEN; \
    using quartic_pool::FIELDS_FLATTEN; \
    using quartic_pool::TENSOR_FLATTEN; \


//     phase-space flattener set to 'FLATTEN' 
//     field-space flattener set to 'FIELDS_FLATTEN' 

//     working type set to 'number' 


    // *********************************************************************************************


    // CLASS FOR quartic CORE
    // contains code and functionality shared by all the compute backends (OpenMP, MPI, OpenCL, CUDA, ...)
    // these backends are implemented by classes which inherit from this common core
    template <typename number>
    class quartic : public canonical_model<number>
      {

        // CONSTRUCTOR, DESTRUCTOR

      public:

        //! constructor
        quartic(local_environment& e, argument_cache& a);

        //! destructor
		    virtual ~quartic();


        // EXTRACT MODEL INFORMATION -- implements a 'model' interface

      public:

        const std::string& get_name() const override { return(quartic_pool::name); }

        const author_db& get_authors() const override { return(quartic_pool::author); }

        const std::string& get_citeguide() const override { return(quartic_pool::citeguide); }

        const std::string& get_description() const override { return(quartic_pool::description); }

        const std::string& get_license() const override { return(quartic_pool::license); }

        unsigned int get_revision() const override { return(quartic_pool::revision); }

        const std::vector<std::string>& get_references() const override { return(quartic_pool::references); }

        const std::vector<std::string>& get_urls() const override { return(quartic_pool::urls); }

        unsigned int get_N_fields() const override { return(1); }

        unsigned int get_N_params() const override { return(1); }

        const std::vector< std::string >& get_field_names() const override { return(quartic_pool::field_names); }

        const std::vector< std::string >& get_f_latex_names() const override { return(quartic_pool::latex_names); }

        const std::vector< std::string >& get_param_names() const override { return(quartic_pool::param_names); }

        const std::vector< std::string >& get_p_latex_names() const override { return(quartic_pool::platx_names); }

        const std::vector< std::string >& get_state_names() const override { return(quartic_pool::state_names); }


        // INDEX FLATTENING FUNCTIONS -- implements an 'abstract_flattener' interface

      public:

        unsigned int flatten(unsigned int a)                                         const override { return(a); };
        unsigned int flatten(unsigned int a, unsigned int b)                         const override { return(2*1*a + b); };
        unsigned int flatten(unsigned int a, unsigned int b, unsigned int c)         const override { return(2*1*2*1*a + 2*1*b + c); };

        unsigned int fields_flatten(unsigned int a)                                  const override { return(a); };
        unsigned int fields_flatten(unsigned int a, unsigned int b)                  const override { return(1*a + b); };
        unsigned int fields_flatten(unsigned int a, unsigned int b, unsigned int c)  const override { return(1*1*a + 1*b + c); };

        unsigned int tensor_flatten(unsigned int a, unsigned int b)                  const override { return(2*a + b); }


        // INDEX TRAITS -- implements an 'abstract_flattener' interface

      public:

        unsigned int species(unsigned int a)                                         const override { return((a >= 1) ? a-1 : a); };
        unsigned int momentum(unsigned int a)                                        const override { return((a >= 1) ? a : a+1); };
        unsigned int is_field(unsigned int a)                                        const override { return(a < 1); }
        unsigned int is_momentum(unsigned int a)                                     const override { return(a >= 1 && a <= 2*1); }


        // COMPUTE BASIC PHYSICAL QUANTITIES -- implements a 'model'/'canonical_model' interface

      public:

        // Over-ride functions inherited from 'model'
        number H(const parameters<number>& __params, const flattened_tensor<number>& __coords) const override;
        number epsilon(const parameters<number>& __params, const flattened_tensor<number>& __coords) const override;
        number eta(const parameters<number>& __params, const flattened_tensor<number>& __coords) const override;

        // Over-ride functions inherited from 'canonical_model'
        number V(const parameters<number>& __params, const flattened_tensor<number>& __coords) const override;


        // INITIAL CONDITIONS HANDLING -- implements a 'model' interface

      public:

        void validate_ics(const parameters<number>& p, const flattened_tensor<number>& input, flattened_tensor<number>& output) override;


        // PARAMETER HANDLING -- implements a 'model' interface

      public:

        void validate_params(const flattened_tensor<number>& input, flattened_tensor<number>& output) override;


        // CALCULATE MODEL-SPECIFIC QUANTITIES -- implements a 'model' interface

      public:

        // calculate gauge transformations to zeta
        void compute_gauge_xfm_1(const twopf_db_task<number>* __task, const flattened_tensor<number>& __state, flattened_tensor<number>& __dN) override;
        void compute_gauge_xfm_2(const twopf_db_task<number>* __task, const flattened_tensor<number>& __state, double __k, double __k1, double __k2, double __N, flattened_tensor<number>& __ddN) override;

        // calculate tensor quantities, including the 'flow' tensors u2, u3 and the basic tensors A, B, C from which u3 is built
        void u2(const twopf_db_task<number>* __task, const flattened_tensor<number>& __fields, double __k, double __N, flattened_tensor<number>& __u2) override;
        void u3(const twopf_db_task<number>* __task, const flattened_tensor<number>& __fields, double __km, double __kn, double __kr, double __N, flattened_tensor<number>& __u3) override;

        void A(const twopf_db_task<number>* __task, const flattened_tensor<number>& __fields, double __km, double __kn, double __kr, double __N, flattened_tensor<number>& __A) override;
        void B(const twopf_db_task<number>* __task, const flattened_tensor<number>& __fields, double __km, double __kn, double __kr, double __N, flattened_tensor<number>& __B) override;
        void C(const twopf_db_task<number>* __task, const flattened_tensor<number>& __fields, double __km, double __kn, double __kr, double __N, flattened_tensor<number>& __C) override;

        // calculate mass matrix
        void M(const integration_task<number>* __task, const flattened_tensor<number>& __fields, double __N, flattened_tensor<number>& __M) override;

        // calculate raw mass spectrum
        void mass_spectrum(const integration_task<number>* __task, const flattened_tensor<number>& __fields, double __N, flattened_tensor<number>& __M, flattened_tensor<number>& __E) override;

        // calculate the sorted mass spectrum, normalized to H^2 if desired
        void sorted_mass_spectrum(const integration_task<number>* __task, const flattened_tensor<number>& __fields, double __N, bool __norm, flattened_tensor<number>& __M, flattened_tensor<number>& __E) override;

        // BACKEND INTERFACE (PARTIAL IMPLEMENTATION -- WE PROVIDE A COMMON BACKGROUND INTEGRATOR)

      public:

        void backend_process_backg(const background_task<number>* tk, backg_history<number>& solution, bool silent=false) override;

        double compute_end_of_inflation(const integration_task<number>* tk, double search_time=CPPTRANSPORT_DEFAULT_END_OF_INFLATION_SEARCH) override;

		    void compute_aH(const integration_task<number>* tk, std::vector<double>& N,
		                    flattened_tensor<number>& log_aH, flattened_tensor<number>& log_a2H2M,
		                    boost::optional<double> largest_k = boost::none) override;

        void compute_H(const integration_task<number>* tk, std::vector<double>& N,
                       flattened_tensor<number>& log_H, boost::optional<double> largest_k = boost::none) override;


        // CALCULATE INITIAL CONDITIONS FOR N-POINT FUNCTIONS

      protected:

        number make_twopf_re_ic(unsigned int __i, unsigned int __j, double __k, double __Ninit,
                                const twopf_db_task<number>* __task, const flattened_tensor<number>& __fields, double __k_norm);

        number make_twopf_im_ic(unsigned int __i, unsigned int __j, double __k, double __Ninit,
                                const twopf_db_task<number>* __task, const flattened_tensor<number>& __fields, double __k_norm);

        number make_twopf_tensor_ic(unsigned int __i, unsigned int __j, double __k, double __Ninit,
                                    const twopf_db_task<number>* __task, const flattened_tensor<number>& __fields, double __k_norm);

        number make_threepf_ic(unsigned int __i, unsigned int __j, unsigned int __k,
                               double kmode_1, double kmode_2, double kmode_3, double __Ninit,
                               const twopf_db_task<number>* __task, const flattened_tensor<number>& __fields, double __k_norm);


        // INTERNAL DATA

      private:

        number* __A_k1k2k3;
        number* __A_k1k3k2;
        number* __B_k1k2k3;
        number* __B_k1k3k2;
        number* __C_k1k2k3;
        number* __C_k1k3k2;

        number* __A_k2k1k3;
        number* __A_k2k3k1;
        number* __B_k2k1k3;
        number* __B_k2k3k1;
        number* __C_k2k1k3;
        number* __C_k2k3k1;

        number* __A_k3k1k2;
        number* __A_k3k2k1;
        number* __B_k3k1k2;
        number* __B_k3k2k1;
        number* __C_k3k1k2;
        number* __C_k3k2k1;

//         ENDIF !fast 

        number* __raw_params;

        //! workspace: Eigen matrix representing mass matrix
        Eigen::Matrix<number, 1, 1> __mass_matrix;

      };


    // integration - background functor
    template <typename number>
    class quartic_background_functor
      {

      public:

        quartic_background_functor(const parameters<number>& p)
          : __params(p)
          {
            __Mp = p.get_Mp();
          }

        void set_up_workspace()
          {
//             ENDIF !fast 

            __raw_params = new number[1];

            const auto& __pvector = __params.get_vector();
            __raw_params[0] = __pvector[0];
          }

        void close_down_workspace()
          {
//             ENDIF !fast 

            delete[] __raw_params;
          }

        void operator ()(const backg_state<number>& __x, backg_state<number>& __dxdt, number __t);

      protected:

        const parameters<number>& __params;

//         ENDIF !fast 

        number* __raw_params;

        number __Mp;

      };


    // integration - observer object for background only
    template <typename number>
    class quartic_background_observer
      {

      public:

        quartic_background_observer(backg_history<number>& h, const time_config_database& t)
          : history(h),
            time_db(t)
          {
            current_step = time_db.record_begin();
          }

        void operator ()(const backg_state<number>& x, number t);

      private:

        backg_history<number>& history;

        const time_config_database& time_db;

        time_config_database::const_record_iterator current_step;

      };


    // CLASS quartic -- CONSTRUCTORS, DESTRUCTORS


    template <typename number>
    quartic<number>::quartic(local_environment& e, argument_cache& a)
      : canonical_model<number>("b5776652-458d-8ed9-f6fb-7efdec3514a3", 201801, e, a)
      {
        __A_k1k2k3 = new number[1 * 1 * 1];
        __A_k1k3k2 = new number[1 * 1 * 1];
        __B_k1k2k3 = new number[1 * 1 * 1];
        __B_k1k3k2 = new number[1 * 1 * 1];
        __C_k1k2k3 = new number[1 * 1 * 1];
        __C_k1k3k2 = new number[1 * 1 * 1];

        __A_k2k1k3 = new number[1 * 1 * 1];
        __A_k2k3k1 = new number[1 * 1 * 1];
        __B_k2k1k3 = new number[1 * 1 * 1];
        __B_k2k3k1 = new number[1 * 1 * 1];
        __C_k2k1k3 = new number[1 * 1 * 1];
        __C_k2k3k1 = new number[1 * 1 * 1];

        __A_k3k1k2 = new number[1 * 1 * 1];
        __A_k3k2k1 = new number[1 * 1 * 1];
        __B_k3k1k2 = new number[1 * 1 * 1];
        __B_k3k2k1 = new number[1 * 1 * 1];
        __C_k3k1k2 = new number[1 * 1 * 1];
        __C_k3k2k1 = new number[1 * 1 * 1];

//         ENDIF !fast 

        __raw_params = new number[1];
      }


    template <typename number>
    quartic<number>::~quartic()
      {
        delete[] __A_k1k2k3;
        delete[] __A_k1k3k2;
        delete[] __B_k1k2k3;
        delete[] __B_k1k3k2;
        delete[] __C_k1k2k3;
        delete[] __C_k1k3k2;

        delete[] __A_k2k1k3;
        delete[] __A_k2k3k1;
        delete[] __B_k2k1k3;
        delete[] __B_k2k3k1;
        delete[] __C_k2k1k3;
        delete[] __C_k2k3k1;

        delete[] __A_k3k1k2;
        delete[] __A_k3k2k1;
        delete[] __B_k3k1k2;
        delete[] __B_k3k2k1;
        delete[] __C_k3k1k2;
        delete[] __C_k3k2k1;

//         ENDIF !fast 

        delete[] __raw_params;
      }


    // INTERFACE: COMPUTE BASIC PHYSICAL QUANTITIES


    template <typename number>
    number quartic<number>::H(const parameters<number>& __params, const flattened_tensor<number>& __coords) const
      {
        assert(__coords.size() == 2*1);
        if(__coords.size() != 2*1)
          {
            std::ostringstream msg;
            msg << CPPTRANSPORT_WRONG_COORDS_A << __coords.size() << CPPTRANSPORT_WRONG_ICS_B << 2*1 << "]";
            throw std::out_of_range(msg.str());
          }

        DEFINE_INDEX_TOOLS
//         release resources 
        const auto __Mp = __params.get_Mp();
        const auto& __param_vector = __params.get_vector();

//         parameters resource set to '__param_vector' 
//         coordinates resource set to '__coords' 

// BEGIN TEMPORARY POOL (sequence=0) 
const auto _t_0_0 = __Mp;
const auto _t_0_1 = -2.0;
const auto _t_0_2 = 1.0/(_t_0_0*_t_0_0);
const auto _t_0_3 = __param_vector[0];
const auto _t_0_4 = __coords[FLATTEN(1)];
const auto _t_0_5 = 2.0;
const auto _t_0_6 = (_t_0_4*_t_0_4);
const auto _t_0_7 = _t_0_6*_t_0_2;
const auto _t_0_8 = -6.0;
const auto _t_0_9 = _t_0_7+_t_0_8;
const auto _t_0_10 = -1.0;
const auto _t_0_11 = 1.0/(_t_0_9);
const auto _t_0_12 = __coords[FLATTEN(0)];
const auto _t_0_13 = 4.0;
const auto _t_0_14 = (_t_0_12*_t_0_12*_t_0_12*_t_0_12);
const auto _t_0_15 = -(1.0/2.0);
const auto _InternalHsq = _t_0_2*_t_0_3*_t_0_11*_t_0_14*_t_0_15;
const auto _t_0_16 = _InternalHsq;
        // END TEMPORARY POOL (sequence=0) 

        return std::sqrt(_t_0_16);
      }


    template <typename number>
    number quartic<number>::epsilon(const parameters<number>& __params, const flattened_tensor<number>& __coords) const
      {
        assert(__coords.size() == 2*1);
        if(__coords.size() != 2*1)
          {
            std::ostringstream msg;
            msg << CPPTRANSPORT_WRONG_COORDS_A << __coords.size() << CPPTRANSPORT_WRONG_ICS_B << 2*1 << "]";
            throw std::out_of_range(msg.str());
          }

        DEFINE_INDEX_TOOLS
//         release resources 
        const auto __Mp = __params.get_Mp();
        // note canonical epsilon doesn't depend on any Lagrangian parameters, only the field derivatives

//         coordinates resource set to '__coords' 

// BEGIN TEMPORARY POOL (sequence=1) 
const auto _t_1_0 = __coords[FLATTEN(1)];
const auto _t_1_1 = 2.0;
const auto _t_1_2 = (_t_1_0*_t_1_0);
const auto _t_1_3 = __Mp;
const auto _t_1_4 = -2.0;
const auto _t_1_5 = 1.0/(_t_1_3*_t_1_3);
const auto _t_1_6 = (1.0/2.0);
const auto _InternalEps = _t_1_2*_t_1_5*_t_1_6;
const auto _t_1_7 = _InternalEps;
        // END TEMPORARY POOL (sequence=1) 

        return _t_1_7;
      }


    template <typename number>
    number quartic<number>::eta(const parameters<number>& __params, const flattened_tensor<number>& __coords) const
      {
        assert(__coords.size() == 2*1);
        if(__coords.size() != 2*1)
          {
            std::ostringstream msg;
            msg << CPPTRANSPORT_WRONG_COORDS_A << __coords.size() << CPPTRANSPORT_WRONG_ICS_B << 2*1 << "]";
            throw std::out_of_range(msg.str());
          }

        DEFINE_INDEX_TOOLS
//         release resources 
        const auto __Mp = __params.get_Mp();
        const auto& __param_vector = __params.get_vector();

//         parameters resource set to '__param_vector' 
//         coordinates resource set to '__coords' 

// BEGIN TEMPORARY POOL (sequence=2) 
const auto _t_2_0 = __coords[FLATTEN(1)];
const auto _t_2_1 = 2.0;
const auto _t_2_2 = (_t_2_0*_t_2_0);
const auto _t_2_3 = __Mp;
const auto _t_2_4 = -2.0;
const auto _t_2_5 = 1.0/(_t_2_3*_t_2_3);
const auto _t_2_6 = (1.0/2.0);
const auto _InternalEps = _t_2_2*_t_2_5*_t_2_6;
const auto _t_2_7 = __param_vector[0];
const auto _t_2_8 = _t_2_2*_t_2_5;
const auto _t_2_9 = -6.0;
const auto _t_2_10 = _t_2_8+_t_2_9;
const auto _t_2_11 = -1.0;
const auto _t_2_12 = 1.0/(_t_2_10);
const auto _t_2_13 = __coords[FLATTEN(0)];
const auto _t_2_14 = 4.0;
const auto _t_2_15 = (_t_2_13*_t_2_13*_t_2_13*_t_2_13);
const auto _t_2_16 = -(1.0/2.0);
const auto _InternalHsq = _t_2_5*_t_2_7*_t_2_12*_t_2_15*_t_2_16;
const auto _t_2_17 = _InternalEps;
const auto _t_2_18 = -3.0;
const auto _t_2_19 = _t_2_17+_t_2_18;
const auto _t_2_20 = _t_2_17*(_t_2_19)*_t_2_1;
const auto _t_2_21 = _InternalHsq;
const auto _t_2_22 = 1.0/_t_2_21;
const auto _t_2_23 = 3.0;
const auto _t_2_24 = (_t_2_13*_t_2_13*_t_2_13);
const auto _t_2_25 = _t_2_0*_t_2_5*_t_2_7*_t_2_22*_t_2_24*_t_2_11;
const auto _t_2_26 = _t_2_20+_t_2_25;
const auto _t_2_27 = 1.0/_t_2_17;
const auto _InternalEta = (_t_2_26)*_t_2_27;
const auto _t_2_28 = _InternalEta;
        // END TEMPORARY POOL (sequence=2) 

        return _t_2_28;
      }


    template <typename number>
    number quartic<number>::V(const parameters<number>& __params, const flattened_tensor<number>& __coords) const
      {
        assert(__coords.size() == 2*1);
        if(__coords.size() != 2*1)
          {
            std::ostringstream msg;
            msg << CPPTRANSPORT_WRONG_COORDS_A << __coords.size() << CPPTRANSPORT_WRONG_ICS_B << 2*1 << "]";
            throw std::out_of_range(msg.str());
          }

        DEFINE_INDEX_TOOLS
//         release resources 
        const auto __Mp = __params.get_Mp();
        const auto& __param_vector = __params.get_vector();

//         parameters resource set to '__param_vector' 
//         coordinates resource set to '__coords' 

// BEGIN TEMPORARY POOL (sequence=3) 
const auto _t_3_0 = __param_vector[0];
const auto _t_3_1 = __coords[FLATTEN(0)];
const auto _t_3_2 = 4.0;
const auto _t_3_3 = (_t_3_1*_t_3_1*_t_3_1*_t_3_1);
const auto _t_3_4 = (1.0/4.0);
const auto _InternalV = _t_3_0*_t_3_3*_t_3_4;
const auto _t_3_5 = _InternalV;
        // END TEMPORARY POOL (sequence=3) 

        return _t_3_5;
      }


//     ENDIF !fast 


    // Handle initial conditions


    template <typename number>
    void quartic<number>::validate_ics(const parameters<number>& __params, const flattened_tensor<number>& __input, flattened_tensor<number>& __output)
      {
        __output.clear();
        __output.reserve(2*1);
        __output.insert(__output.end(), __input.begin(), __input.end());

        if(__input.size() == 1)  // initial conditions for momenta *were not* supplied -- need to compute them
          {
            DEFINE_INDEX_TOOLS
//             release resources 

            // supply the missing initial conditions using a slow-roll approximation
            const auto __Mp = __params.get_Mp();

            const auto& __pvector = __params.get_vector();
//             parameters resource set to '__pvector' 
//             coordinates resource set to '__input' 

// BEGIN TEMPORARY POOL (sequence=4) 
const auto _t_4_0 = __input[FLATTEN(0)];
const auto _t_4_1 = 4.0;
const auto _t_4_2 = (_t_4_0*_t_4_0*_t_4_0*_t_4_0);
const auto _t_4_3 = __pvector[0];
const auto _t_4_4 = (1.0/4.0);
const auto _InternalV = _t_4_2*_t_4_3*_t_4_4;
const auto _t_4_6 = -1.0;
const auto _t_4_5 = _InternalV;
const auto _t_4_7 = 1.0/_t_4_5;
const auto _t_4_8 = 3.0;
const auto _t_4_9 = (_t_4_0*_t_4_0*_t_4_0);
const auto _t_4_11 = 2.0;
const auto _t_4_10 = __Mp;
const auto _t_4_12 = (_t_4_10*_t_4_10);
const auto _t_4_13 = _t_4_7*_t_4_9*_t_4_12*_t_4_3*_t_4_6;
            // END TEMPORARY POOL (sequence=4) 

            // force unroll to make explicit that we wish to populate array elements
            __output.push_back(_t_4_13);
          }
        else if(__input.size() == 2*1)  // initial conditions for momenta *were* supplied
          {
            // need do nothing
          }
        else
          {
            std::ostringstream msg;

            msg << CPPTRANSPORT_WRONG_ICS_A << __input.size()
                << CPPTRANSPORT_WRONG_ICS_B << 1
                << CPPTRANSPORT_WRONG_ICS_C << 2*1 << "]";

            throw std::out_of_range(msg.str());
          }
      }


    // Handle parameters


    template <typename number>
    void quartic<number>::validate_params(const flattened_tensor<number>& input, flattened_tensor<number>& output)
      {
        output.clear();

        if(input.size() == 1)
          {
            output.assign(input.begin(), input.end());
          }
        else
          {
            std::ostringstream msg;

            msg << CPPTRANSPORT_WRONG_PARAMS_A << input.size() << CPPTRANSPORT_WRONG_PARAMS_B << 1 << "]";

            throw std::out_of_range(msg.str());
          }
      }


    // set up initial conditions for the real part of the equal-time two-point function
    // __i,__j  -- label component of the twopf for which we wish to compute initial conditions
    // __k      -- *comoving normalized* wavenumber for which we wish to assign initial conditions
    // __Ninit  -- initial time
    // __fields -- vector of initial conditions for the background fields (or fields+momenta)
    template <typename number>
    number quartic<number>::make_twopf_re_ic(unsigned int __i, unsigned int __j, double __k, double __Ninit,
                                            const twopf_db_task<number>* __task, const flattened_tensor<number>& __fields,
                                            double __k_norm)
      {
        DEFINE_INDEX_TOOLS
//         release resources 
        const auto& __pvector = __task->get_params().get_vector();
        __raw_params[0] = __pvector[0];

        const auto __Mp = __task->get_params().get_Mp();
        const auto __a = std::exp(__Ninit - __task->get_N_horizon_crossing() + __task->get_astar_normalization());

//         parameters resource set to '__raw_params' 
//         coordinates resource set to '__fields' 

// BEGIN TEMPORARY POOL (sequence=5) 
const auto _t_5_0 = __Mp;
const auto _t_5_1 = -2.0;
const auto _t_5_2 = 1.0/(_t_5_0*_t_5_0);
const auto _t_5_3 = __fields[FLATTEN(1)];
const auto _t_5_4 = 2.0;
const auto _t_5_5 = (_t_5_3*_t_5_3);
const auto _t_5_6 = _t_5_2*_t_5_5;
const auto _t_5_7 = -6.0;
const auto _t_5_8 = _t_5_6+_t_5_7;
const auto _t_5_9 = -1.0;
const auto _t_5_10 = 1.0/(_t_5_8);
const auto _t_5_11 = __raw_params[0];
const auto _t_5_12 = __fields[FLATTEN(0)];
const auto _t_5_13 = 4.0;
const auto _t_5_14 = (_t_5_12*_t_5_12*_t_5_12*_t_5_12);
const auto _t_5_15 = -(1.0/2.0);
const auto _InternalHsq = _t_5_10*_t_5_11*_t_5_2*_t_5_14*_t_5_15;
const auto _t_5_16 = _InternalHsq;
        // END TEMPORARY POOL (sequence=5) 

        const auto __Hsq = _t_5_16;
        const auto __N = std::log(__k / (__a * std::sqrt(__Hsq)));

        number __tpf = 0.0;

        // NOTE - conventions for the scale factor are
        //   a = exp(t), where t is the user-defined time (usually = 0 at the start of the integration)
        //   so usually this is zero

        if(IS_FIELD(__i) && IS_FIELD(__j))              // field-field correlation function
          {
            // LEADING-ORDER INITIAL CONDITION
            auto __leading = (SPECIES(__i) == SPECIES(__j) ? 1.0 : 0.0);
            auto __subl    = 0.0;
            auto __subsubl = 0.0;

            // NEXT-ORDER INITIAL CONDITION - induces rapid onset of subhorizon oscillations
//              auto __leading = (SPECIES(__i) == SPECIES(__j) ? 1.0 : 0.0) * (1.0 - 2.0*__eps*(1.0-__N));
//              auto __subl    = (SPECIES(__i) == SPECIES(__j) ? 1.0 : 0.0) * (1.0 - 2.0*__eps*(1.0-__N))
//                               + (3.0/2.0)*__M[SPECIES(__i)][SPECIES(__j)];
//              auto __subsubl = (9.0/4.0)*__M[SPECIES(__i)][SPECIES(__j)];

            __tpf = + __leading                             / (2.0*__k*__a*__a)
                    + __subl*__Hsq                          / (2.0*__k*__k*__k)
                    + __subsubl*__Hsq*__Hsq*__a*__a / (2.0*__k*__k*__k*__k*__k);
          }
        else if((IS_FIELD(__i) && IS_MOMENTUM(__j))     // field-momentum or momentum-field correlation function
                || (IS_MOMENTUM(__i) && IS_FIELD(__j)))
          {
            // LEADING-ORDER INITIAL CONDITION
            auto __leading = (SPECIES(__i) == SPECIES(__j) ? 1.0 : 0.0) * (-1.0);
            auto __subl    = 0.0;
            auto __subsubl = 0.0;

            // NEXT-ORDER INITIAL CONDITION - induces slow onset of subhorizon oscillations
//              auto __leading = (SPECIES(__i) == SPECIES(__j) ? 1.0 : 0.0) * (-1.0 + __eps*(1.0-2.0*__N));
//              auto __subl    = (SPECIES(__i) == SPECIES(__j) ? 1.0 : 0.0) * (- __eps);
//              auto __subsubl = (9.0/4.0)*__M[SPECIES(__i)][SPECIES(__j)];

            __tpf = + __leading                             / (2.0*__k*__a*__a)
                    + __subl*__Hsq                          / (2.0*__k*__k*__k)
                    + __subsubl*__Hsq*__Hsq*__a*__a / (2.0*__k*__k*__k*__k*__k);
          }
        else if(IS_MOMENTUM(__i) && IS_MOMENTUM(__j))   // momentum-momentum correlation function
          {
            // LEADING-ORDER INITIAL CONDITION
            auto __leading = (SPECIES(__i) == SPECIES(__j) ? 1.0 : 0.0);
            auto __subl    = 0.0;
            auto __subsubl = 0.0;

            // NEXT-ORDER INITIAL CONDITION - induces rapid onset of subhorizon oscillations
//              auto __leading = (SPECIES(__i) == SPECIES(__j) ? 1.0 : 0.0) * (1.0 - 2.0*__eps*(1.0-__N));
//              auto __subl    = (SPECIES(__i) == SPECIES(__j) ? 1.0 : 0.0) * 2.0*__eps
//                               - (3.0/2.0)*__M[SPECIES(__i)][SPECIES(__j)];
//              auto __subsubl = - (3.0/4.0)*__M[SPECIES(__i)][SPECIES(__j)];

            __tpf = + __k*__leading   / (2.0*__Hsq*__a*__a*__a*__a)
                    + __subl          / (2.0*__k*__a*__a)
                    + __subsubl*__Hsq / (2.0*__k*__k*__k);
          }
        else
          {
            assert(false);
          }

        // return value, rescaled to give dimensionless correlation function
        return(__tpf * __k_norm*__k_norm*__k_norm);
      }


  // set up initial conditions for the imaginary part of the equal-time two-point function
  template <typename number>
  number quartic<number>::make_twopf_im_ic(unsigned int __i, unsigned int __j, double __k, double __Ninit,
                                          const twopf_db_task<number>* __task, const flattened_tensor<number>& __fields,
                                          double __k_norm)
    {
      DEFINE_INDEX_TOOLS
//       release resources 
      const auto& __pvector = __task->get_params().get_vector();
      __raw_params[0] = __pvector[0];

      const auto __Mp = __task->get_params().get_Mp();
      const auto __a = std::exp(__Ninit - __task->get_N_horizon_crossing() + __task->get_astar_normalization());

//       parameters resource set to '__raw_params' 
//       coordinates resource set to '__fields' 

// BEGIN TEMPORARY POOL (sequence=6) 
const auto _t_6_0 = __Mp;
const auto _t_6_1 = -2.0;
const auto _t_6_2 = 1.0/(_t_6_0*_t_6_0);
const auto _t_6_3 = __fields[FLATTEN(1)];
const auto _t_6_4 = 2.0;
const auto _t_6_5 = (_t_6_3*_t_6_3);
const auto _t_6_6 = _t_6_2*_t_6_5;
const auto _t_6_7 = -6.0;
const auto _t_6_8 = _t_6_6+_t_6_7;
const auto _t_6_9 = -1.0;
const auto _t_6_10 = 1.0/(_t_6_8);
const auto _t_6_11 = __raw_params[0];
const auto _t_6_12 = __fields[FLATTEN(0)];
const auto _t_6_13 = 4.0;
const auto _t_6_14 = (_t_6_12*_t_6_12*_t_6_12*_t_6_12);
const auto _t_6_15 = -(1.0/2.0);
const auto _InternalHsq = _t_6_10*_t_6_11*_t_6_2*_t_6_14*_t_6_15;
const auto _t_6_16 = _InternalHsq;
      // END TEMPORARY POOL (sequence=6) 

      const auto __Hsq = _t_6_16;
      const auto __N = std::log(__k / (__a * std::sqrt(__Hsq)));

      number __tpf = 0.0;

      // only the field-momentum and momentum-field correlation functions have imaginary parts
      if(IS_FIELD(__i) && IS_MOMENTUM(__j))
        {
          // LEADING-ORDER INITIAL CONDITION
          auto __leading = (SPECIES(__i) == SPECIES(__j) ? 1.0 : 0.0);

          // NEXT-ORDER INITIAL CONDITION
//            auto __leading = (SPECIES(__i) == SPECIES(__j) ? 1.0 : 0.0) * (1.0 - 2.0*__eps*(1.0-__N));

          __tpf = + __leading / (2.0*std::sqrt(__Hsq)*__a*__a*__a);
        }
      else if(IS_MOMENTUM(__i) && IS_FIELD(__j))
        {
          // LEADING-ORDER INITIAL CONDITION
          auto __leading = (SPECIES(__i) == SPECIES(__j) ? 1.0 : 0.0);

          // NEXT-ORDER INITIAL CONDITION
//            auto __leading = (SPECIES(__i) == SPECIES(__j) ? 1.0 : 0.0) * (1.0 - 2.0*__eps*(1.0-__N));

          __tpf = - __leading / (2.0*std::sqrt(__Hsq)*__a*__a*__a);
        }

      // return value, rescaled to give dimensionless correlation function
      return(__tpf * __k_norm*__k_norm*__k_norm);
    }


    // set up initial conditions for the real part of the equal-time tensor two-point function
    template <typename number>
    number quartic<number>::make_twopf_tensor_ic(unsigned int __i, unsigned int __j, double __k, double __Ninit,
                                                const twopf_db_task<number>* __task, const flattened_tensor<number>& __fields,
                                                double __k_norm)
      {
        DEFINE_INDEX_TOOLS
//         release resources 
        const auto& __pvector = __task->get_params().get_vector();
        __raw_params[0] = __pvector[0];

        const auto __Mp = __task->get_params().get_Mp();
        const auto __a = std::exp(__Ninit - __task->get_N_horizon_crossing() + __task->get_astar_normalization());

// BEGIN TEMPORARY POOL (sequence=7) 
const auto _t_7_0 = __Mp;
const auto _t_7_1 = -2.0;
const auto _t_7_2 = 1.0/(_t_7_0*_t_7_0);
const auto _t_7_3 = __fields[FLATTEN(1)];
const auto _t_7_4 = 2.0;
const auto _t_7_5 = (_t_7_3*_t_7_3);
const auto _t_7_6 = _t_7_2*_t_7_5;
const auto _t_7_7 = -6.0;
const auto _t_7_8 = _t_7_6+_t_7_7;
const auto _t_7_9 = -1.0;
const auto _t_7_10 = 1.0/(_t_7_8);
const auto _t_7_11 = __raw_params[0];
const auto _t_7_12 = __fields[FLATTEN(0)];
const auto _t_7_13 = 4.0;
const auto _t_7_14 = (_t_7_12*_t_7_12*_t_7_12*_t_7_12);
const auto _t_7_15 = -(1.0/2.0);
const auto _InternalHsq = _t_7_10*_t_7_11*_t_7_2*_t_7_14*_t_7_15;
const auto _t_7_16 = _InternalHsq;
        // END TEMPORARY POOL (sequence=7) 

//         parameters resource set to '__raw_params' 
//         coordinates resource set to '__fields' 

        const auto __Hsq = _t_7_16;
        const auto __N = std::log(__k / (__a * std::sqrt(__Hsq)));

        number __tpf = 0.0;

        if(__i == 0 && __j == 0)                                      // h-h correlation function
          {
            // LEADING-ORDER INITIAL CONDITION
            __tpf = 1.0 / (__Mp*__Mp*__k*__a*__a);
//            __tpf = 1.0 / (2.0*__k*__a*__a);
          }
        else if((__i == 0 && __j == 1) || (__i == 1 && __j == 0))     // h-dh or dh-h correlation function
          {
            // LEADING ORDER INITIAL CONDITION
            __tpf = -1.0 / (__Mp*__Mp*__k*__a*__a);
//            __tpf = -1.0 / (2.0*__k*__a*__a);
          }
        else if(__i == 1 && __j == 1)                                 // dh-dh correlation function
          {
            // LEADING ORDER INITIAL CONDITION
            __tpf = __k / (__Mp*__Mp*__Hsq*__a*__a*__a*__a);
//            __tpf = __k / (2.0*__Hsq*__a*__a*__a*__a);
          }
        else
          {
            assert(false);
          }

        // return value, rescaled to give dimensionless correlation function
        return(__tpf * __k_norm*__k_norm*__k_norm);
      }


    // set up initial conditions for the real part of the equal-time three-point function
    template <typename number>
    number quartic<number>::make_threepf_ic(unsigned int __i, unsigned int __j, unsigned int __k,
                                           double __k1, double __k2, double __k3, double __Ninit,
                                           const twopf_db_task<number>* __task, const flattened_tensor<number>& __fields,
                                           double __k_norm)
      {
        DEFINE_INDEX_TOOLS
//         release resources 
        const auto& __pvector = __task->get_params().get_vector();
        __raw_params[0] = __pvector[0];

        const auto __Mp = __task->get_params().get_Mp();
        const auto __a = std::exp(__Ninit - __task->get_N_horizon_crossing() + __task->get_astar_normalization());

//         ENDIF !fast 

// BEGIN TEMPORARY POOL (sequence=8) 
const auto _t_8_0 = __Mp;
const auto _t_8_1 = -2.0;
const auto _t_8_2 = 1.0/(_t_8_0*_t_8_0);
const auto _t_8_3 = __fields[FLATTEN(1)];
const auto _t_8_4 = 2.0;
const auto _t_8_5 = (_t_8_3*_t_8_3);
const auto _t_8_6 = _t_8_2*_t_8_5;
const auto _t_8_7 = -6.0;
const auto _t_8_8 = _t_8_6+_t_8_7;
const auto _t_8_9 = -1.0;
const auto _t_8_10 = 1.0/(_t_8_8);
const auto _t_8_11 = __raw_params[0];
const auto _t_8_12 = __fields[FLATTEN(0)];
const auto _t_8_13 = 4.0;
const auto _t_8_14 = (_t_8_12*_t_8_12*_t_8_12*_t_8_12);
const auto _t_8_15 = -(1.0/2.0);
const auto _InternalHsq = _t_8_10*_t_8_11*_t_8_2*_t_8_14*_t_8_15;
const auto _t_8_16 = _InternalHsq;
const auto _t_8_17 = (1.0/2.0);
const auto _InternalEps = _t_8_2*_t_8_5*_t_8_17;
const auto _t_8_18 = -4.0;
const auto _t_8_19 = 1.0/(_t_8_0*_t_8_0*_t_8_0*_t_8_0);
const auto _t_8_20 = 3.0;
const auto _t_8_21 = (_t_8_3*_t_8_3*_t_8_3);
const auto _t_8_22 = (1.0/4.0);
const auto _t_8_23 = _t_8_19*_t_8_21*_t_8_22;
const auto _t_8_24 = __k1;
const auto _t_8_25 = (_t_8_24*_t_8_24);
const auto _t_8_26 = __k3;
const auto _t_8_27 = (_t_8_26*_t_8_26);
const auto _t_8_28 = _t_8_27*_t_8_9;
const auto _t_8_29 = __k2;
const auto _t_8_30 = (_t_8_29*_t_8_29);
const auto _t_8_31 = _t_8_25+_t_8_28+_t_8_30;
const auto _t_8_32 = 1.0/(_t_8_29*_t_8_29);
const auto _t_8_33 = _InternalEps;
const auto _t_8_34 = -3.0;
const auto _t_8_35 = _t_8_33+_t_8_34;
const auto _t_8_36 = _t_8_3*(_t_8_35);
const auto _t_8_37 = (_t_8_12*_t_8_12*_t_8_12);
const auto _t_8_38 = 1.0/_t_8_16;
const auto _t_8_39 = _t_8_11*_t_8_37*_t_8_38*_t_8_9;
const auto _t_8_40 = _t_8_36+_t_8_39;
const auto _t_8_41 = (_t_8_31)*_t_8_2*_t_8_32*(_t_8_40)*_t_8_22;
const auto _t_8_42 = _t_8_30*_t_8_9;
const auto _t_8_43 = _t_8_25+_t_8_28+_t_8_42;
const auto _t_8_44 = ((_t_8_43)*(_t_8_43));
const auto _t_8_45 = 1.0/(_t_8_26*_t_8_26);
const auto _t_8_46 = _t_8_32*_t_8_44*_t_8_45;
const auto _t_8_47 = _t_8_46+_t_8_18;
const auto _t_8_48 = (1.0/32.0);
const auto _t_8_49 = _t_8_19*(_t_8_40)*_t_8_5*(_t_8_47)*_t_8_48;
const auto _t_8_50 = _t_8_25+_t_8_27+_t_8_42;
const auto _t_8_51 = ((_t_8_50)*(_t_8_50));
const auto _t_8_52 = 1.0/(_t_8_24*_t_8_24);
const auto _t_8_53 = _t_8_51*_t_8_52*_t_8_45;
const auto _t_8_54 = _t_8_53+_t_8_18;
const auto _t_8_55 = _t_8_19*(_t_8_40)*(_t_8_54)*_t_8_5*_t_8_48;
const auto _t_8_56 = _t_8_52*(_t_8_31)*_t_8_2*(_t_8_40)*_t_8_22;
const auto _t_8_57 = _t_8_23+_t_8_41+_t_8_49+_t_8_55+_t_8_56;
const auto _t_8_58 = ((_t_8_31)*(_t_8_31));
const auto _t_8_59 = _t_8_52*_t_8_58*_t_8_32;
const auto _t_8_60 = _t_8_59+_t_8_18;
const auto _t_8_61 = _t_8_19*(_t_8_40)*(_t_8_60)*_t_8_5*_t_8_48;
const auto _t_8_62 = (_t_8_50)*_t_8_2*(_t_8_40)*_t_8_45*_t_8_22;
const auto _t_8_63 = (_t_8_50)*_t_8_52*_t_8_2*(_t_8_40)*_t_8_22;
const auto _t_8_64 = _t_8_23+_t_8_61+_t_8_49+_t_8_62+_t_8_63;
const auto _t_8_65 = _t_8_2*_t_8_32*(_t_8_43)*_t_8_3*_t_8_22;
const auto _t_8_66 = -(1.0/4.0);
const auto _t_8_67 = (_t_8_50)*_t_8_52*_t_8_2*_t_8_3*_t_8_66;
const auto _t_8_68 = -(1.0/32.0);
const auto _t_8_69 = _t_8_19*(_t_8_60)*_t_8_21*_t_8_68;
const auto _t_8_70 = _t_8_2*_t_8_3*_t_8_15;
const auto _t_8_71 = _t_8_65+_t_8_67+_t_8_69+_t_8_70;
const auto _t_8_72 = _t_8_2*(_t_8_43)*_t_8_45*_t_8_3*_t_8_22;
const auto _t_8_73 = _t_8_52*(_t_8_31)*_t_8_2*_t_8_3*_t_8_66;
const auto _t_8_74 = _t_8_19*(_t_8_54)*_t_8_21*_t_8_68;
const auto _t_8_75 = _t_8_72+_t_8_73+_t_8_70+_t_8_74;
const auto _t_8_76 = ((_t_8_40)*(_t_8_40));
const auto _t_8_77 = -(1.0/96.0);
const auto _t_8_78 = _t_8_19*_t_8_76*(_t_8_54)*_t_8_3*_t_8_77;
const auto _t_8_79 = _t_8_19*_t_8_76*(_t_8_60)*_t_8_3*_t_8_77;
const auto _t_8_80 = _t_8_19*(_t_8_40)*_t_8_5*_t_8_22;
const auto _t_8_81 = _t_8_11*_t_8_12*_t_8_38*_t_8_1;
const auto _t_8_82 = _t_8_19*_t_8_76*_t_8_3*(_t_8_47)*_t_8_77;
const auto _t_8_83 = _t_8_19*_t_8_21*(_t_8_35)*_t_8_66;
const auto _t_8_84 = (_t_8_12*_t_8_12);
const auto _t_8_85 = -(3.0/2.0);
const auto _t_8_86 = _t_8_11*_t_8_2*_t_8_84*_t_8_38*_t_8_3*_t_8_85;
const auto _t_8_87 = _t_8_78+_t_8_79+_t_8_80+_t_8_81+_t_8_82+_t_8_83+_t_8_86;
const auto _t_8_88 = _t_8_2*(_t_8_43)*(_t_8_40)*_t_8_45*_t_8_66;
const auto _t_8_89 = _t_8_2*_t_8_32*(_t_8_43)*(_t_8_40)*_t_8_66;
const auto _t_8_90 = _t_8_23+_t_8_88+_t_8_89+_t_8_61+_t_8_55;
const auto _t_8_91 = (_t_8_50)*_t_8_2*_t_8_45*_t_8_3*_t_8_66;
const auto _t_8_92 = _t_8_19*_t_8_21*(_t_8_47)*_t_8_68;
const auto _t_8_93 = (_t_8_31)*_t_8_2*_t_8_32*_t_8_3*_t_8_66;
const auto _t_8_94 = _t_8_91+_t_8_92+_t_8_70+_t_8_93;
        // END TEMPORARY POOL (sequence=8) 

//         parameters resource set to '__raw_params' 
//         coordinates resource set to '__fields' 
//         ENDIF !fast 

        const auto __Hsq = _t_8_16;

        const auto __kt      = __k1 + __k2 + __k3;
        const auto __Ksq     = __k1*__k2 + __k1*__k3 + __k2*__k3;
        const auto __kprod3  = 4.0 * __k1*__k1*__k1 * __k2*__k2*__k2 * __k3*__k3*__k3;

        const auto __k2dotk3 = (__k1*__k1 - __k2*__k2 - __k3*__k3)/2.0;
        const auto __k1dotk3 = (__k2*__k2 - __k1*__k1 - __k3*__k3)/2.0;
        const auto __k1dotk2 = (__k3*__k3 - __k1*__k1 - __k2*__k2)/2.0;

        number __tpf = 0.0;

//         set replacement rule 'B_k1k2k3' to "__B_k1k2k3[FIELDS_FLATTEN($a,$b,$c)]" 
//         set replacement rule 'B_k1k3k2' to "__B_k1k3k2[FIELDS_FLATTEN($a,$b,$c)]" 
//         set replacement rule 'C_k1k2k3' to "__C_k1k2k3[FIELDS_FLATTEN($a,$b,$c)]" 
//         set replacement rule 'C_k1k3k2' to "__C_k1k3k2[FIELDS_FLATTEN($a,$b,$c)]" 
//         set replacement rule 'Atilde_k1k2k3' to "__A_k1k2k3[FIELDS_FLATTEN($a,$b,$c)]" 
//         set replacement rule 'Atilde_k1k3k2' to "__A_k1k3k2[FIELDS_FLATTEN($a,$b,$c)]" 

//         set replacement rule 'B_k2k3k1' to "__B_k2k3k1[FIELDS_FLATTEN($a,$b,$c)]" 
//         set replacement rule 'B_k2k1k3' to "__B_k2k1k3[FIELDS_FLATTEN($a,$b,$c)]" 
//         set replacement rule 'C_k2k3k1' to "__C_k2k3k1[FIELDS_FLATTEN($a,$b,$c)]" 
//         set replacement rule 'C_k2k1k3' to "__C_k2k1k3[FIELDS_FLATTEN($a,$b,$c)]" 
//         set replacement rule 'Atilde_k2k3k1' to "__A_k2k3k1[FIELDS_FLATTEN($a,$b,$c)]" 
//         set replacement rule 'Atilde_k2k1k3' to "__A_k2k1k3[FIELDS_FLATTEN($a,$b,$c)]" 

//         set replacement rule 'B_k3k1k2' to "__B_k3k1k2[FIELDS_FLATTEN($a,$b,$c)]" 
//         set replacement rule 'B_k3k2k1' to "__B_k3k2k1[FIELDS_FLATTEN($a,$b,$c)]" 
//         set replacement rule 'C_k3k1k2' to "__C_k3k1k2[FIELDS_FLATTEN($a,$b,$c)]" 
//         set replacement rule 'C_k3k2k1' to "__C_k3k2k1[FIELDS_FLATTEN($a,$b,$c)]" 
//         set replacement rule 'Atilde_k3k1k2' to "__A_k3k1k2[FIELDS_FLATTEN($a,$b,$c)]" 
//         set replacement rule 'Atilde_k3k2k1' to "__A_k3k2k1[FIELDS_FLATTEN($a,$b,$c)]" 

        __B_k1k2k3[FIELDS_FLATTEN(0,0,0)] = _t_8_57;
        __B_k1k3k2[FIELDS_FLATTEN(0,0,0)] = _t_8_64;
        __C_k1k2k3[FIELDS_FLATTEN(0,0,0)] = _t_8_71;
        __C_k1k3k2[FIELDS_FLATTEN(0,0,0)] = _t_8_75;
        __A_k1k2k3[FIELDS_FLATTEN(0,0,0)] = _t_8_87;
        __A_k1k3k2[FIELDS_FLATTEN(0,0,0)] = _t_8_87;

        __B_k2k3k1[FIELDS_FLATTEN(0,0,0)] = _t_8_90;
        __B_k2k1k3[FIELDS_FLATTEN(0,0,0)] = _t_8_57;
        __C_k2k3k1[FIELDS_FLATTEN(0,0,0)] = _t_8_94;
        __C_k2k1k3[FIELDS_FLATTEN(0,0,0)] = _t_8_71;
        __A_k2k3k1[FIELDS_FLATTEN(0,0,0)] = _t_8_87;
        __A_k2k1k3[FIELDS_FLATTEN(0,0,0)] = _t_8_87;

        __B_k3k1k2[FIELDS_FLATTEN(0,0,0)] = _t_8_64;
        __B_k3k2k1[FIELDS_FLATTEN(0,0,0)] = _t_8_90;
        __C_k3k1k2[FIELDS_FLATTEN(0,0,0)] = _t_8_75;
        __C_k3k2k1[FIELDS_FLATTEN(0,0,0)] = _t_8_94;
        __A_k3k1k2[FIELDS_FLATTEN(0,0,0)] = _t_8_87;
        __A_k3k2k1[FIELDS_FLATTEN(0,0,0)] = _t_8_87;

        unsigned int total_fields  = (IS_FIELD(__i)    ? 1 : 0) + (IS_FIELD(__j)    ? 1 : 0) + (IS_FIELD(__k)    ? 1 : 0);
        unsigned int total_momenta = (IS_MOMENTUM(__i) ? 1 : 0) + (IS_MOMENTUM(__j) ? 1 : 0) + (IS_MOMENTUM(__k) ? 1 : 0);

        assert(total_fields + total_momenta == 3);

        switch(total_fields)
          {
            case 3:   // field-field-field correlation function
              {
                assert(total_fields == 3);

                // prefactor here is dimension 5
                auto __prefactor = (__k1*__k1) * (__k2*__k2) * (__k3*__k3) / (__kt * __a*__a*__a*__a);

                // these components are dimension 1
                // note factor of 2 compared to analytic calculation, from symmetrization over beta, gamma
                __tpf  = (SPECIES(__j) == SPECIES(__k) ? __fields[MOMENTUM(__i)] : 0.0) * __k2dotk3 / (2.0*__Mp*__Mp);
                __tpf += (SPECIES(__i) == SPECIES(__k) ? __fields[MOMENTUM(__j)] : 0.0) * __k1dotk3 / (2.0*__Mp*__Mp);
                __tpf += (SPECIES(__i) == SPECIES(__j) ? __fields[MOMENTUM(__k)] : 0.0) * __k1dotk2 / (2.0*__Mp*__Mp);

                // these components are dimension 1
                __tpf += - (__C_k1k2k3[FIELDS_FLATTEN(SPECIES(__i), SPECIES(__j), SPECIES(__k))] + __C_k2k1k3[FIELDS_FLATTEN(SPECIES(__j), SPECIES(__i), SPECIES(__k))])*__k1*__k2/2.0;
                __tpf += - (__C_k1k3k2[FIELDS_FLATTEN(SPECIES(__i), SPECIES(__k), SPECIES(__j))] + __C_k3k1k2[FIELDS_FLATTEN(SPECIES(__k), SPECIES(__i), SPECIES(__j))])*__k1*__k3/2.0;
                __tpf += - (__C_k2k3k1[FIELDS_FLATTEN(SPECIES(__j), SPECIES(__k), SPECIES(__i))] + __C_k3k2k1[FIELDS_FLATTEN(SPECIES(__k), SPECIES(__j), SPECIES(__i))])*__k2*__k3/2.0;

                // these components are dimension 1
                __tpf += __a*__a * __Hsq * (__A_k1k2k3[FIELDS_FLATTEN(SPECIES(__i), SPECIES(__j), SPECIES(__k))] + __A_k1k3k2[FIELDS_FLATTEN(SPECIES(__i), SPECIES(__k), SPECIES(__j))]) / 2.0;
                __tpf += __a*__a * __Hsq * (__A_k2k1k3[FIELDS_FLATTEN(SPECIES(__j), SPECIES(__i), SPECIES(__k))] + __A_k2k3k1[FIELDS_FLATTEN(SPECIES(__j), SPECIES(__k), SPECIES(__i))]) / 2.0;
                __tpf += __a*__a * __Hsq * (__A_k3k1k2[FIELDS_FLATTEN(SPECIES(__k), SPECIES(__i), SPECIES(__j))] + __A_k3k2k1[FIELDS_FLATTEN(SPECIES(__k), SPECIES(__j), SPECIES(__i))]) / 2.0;

                // these components are dimension 1
                __tpf += __a*__a * __Hsq * (__B_k1k2k3[FIELDS_FLATTEN(SPECIES(__i), SPECIES(__j), SPECIES(__k))] + __B_k2k1k3[FIELDS_FLATTEN(SPECIES(__j), SPECIES(__i), SPECIES(__k))]) * ( (__k1+__k2)*__k3 / (__k1*__k2) - __Ksq / (__k1*__k2)) / 2.0;
                __tpf += __a*__a * __Hsq * (__B_k1k3k2[FIELDS_FLATTEN(SPECIES(__i), SPECIES(__k), SPECIES(__j))] + __B_k3k1k2[FIELDS_FLATTEN(SPECIES(__k), SPECIES(__i), SPECIES(__j))]) * ( (__k1+__k3)*__k2 / (__k1*__k3) - __Ksq / (__k1*__k3)) / 2.0;
                __tpf += __a*__a * __Hsq * (__B_k2k3k1[FIELDS_FLATTEN(SPECIES(__j), SPECIES(__k), SPECIES(__i))] + __B_k3k2k1[FIELDS_FLATTEN(SPECIES(__k), SPECIES(__j), SPECIES(__i))]) * ( (__k2+__k3)*__k1 / (__k2*__k3) - __Ksq / (__k2*__k3)) / 2.0;

                __tpf *= __prefactor / __kprod3;

                break;
              }

            case 2:   // field-field-momentum correlation function
              {
                assert(total_fields == 2);

                auto __momentum_k = (IS_MOMENTUM(__i) ? __k1 : 0.0) + (IS_MOMENTUM(__j) ? __k2 : 0.0) + (IS_MOMENTUM(__k) ? __k3 : 0.0);

                // prefactor has dimension 3
                auto __prefactor = __k1*__k2*__k3 / (__kt * __a*__a*__a*__a);

                // component of prefactor that should not be symmetrized; has dimension 2
                auto __mom_factor_1 = __momentum_k*__momentum_k*(__kt-__momentum_k);

                // these components are dimension 1
                // note factor of 2 compared to analytic calculation, from symmetrization over beta, gamma
                auto __tpf_1  = (SPECIES(__j) == SPECIES(__k) ? __fields[MOMENTUM(__i)] : 0.0) * __k2dotk3 / (2.0*__Mp*__Mp);
                     __tpf_1 += (SPECIES(__i) == SPECIES(__k) ? __fields[MOMENTUM(__j)] : 0.0) * __k1dotk3 / (2.0*__Mp*__Mp);
                     __tpf_1 += (SPECIES(__i) == SPECIES(__j) ? __fields[MOMENTUM(__k)] : 0.0) * __k1dotk2 / (2.0*__Mp*__Mp);

                // these components are dimension 1
                     __tpf_1 += - (__C_k1k2k3[FIELDS_FLATTEN(SPECIES(__i), SPECIES(__j), SPECIES(__k))] + __C_k2k1k3[FIELDS_FLATTEN(SPECIES(__j), SPECIES(__i), SPECIES(__k))])*__k1*__k2 / 2.0;
                     __tpf_1 += - (__C_k1k3k2[FIELDS_FLATTEN(SPECIES(__i), SPECIES(__k), SPECIES(__j))] + __C_k3k1k2[FIELDS_FLATTEN(SPECIES(__k), SPECIES(__i), SPECIES(__j))])*__k1*__k3 / 2.0;
                     __tpf_1 += - (__C_k2k3k1[FIELDS_FLATTEN(SPECIES(__j), SPECIES(__k), SPECIES(__i))] + __C_k3k2k1[FIELDS_FLATTEN(SPECIES(__k), SPECIES(__j), SPECIES(__i))])*__k2*__k3 / 2.0;

                // these components are dimension 1
                     __tpf_1 += __a*__a * __Hsq * (__A_k1k2k3[FIELDS_FLATTEN(SPECIES(__i), SPECIES(__j), SPECIES(__k))] + __A_k1k3k2[FIELDS_FLATTEN(SPECIES(__i), SPECIES(__k), SPECIES(__j))]) / 2.0;
                     __tpf_1 += __a*__a * __Hsq * (__A_k2k1k3[FIELDS_FLATTEN(SPECIES(__j), SPECIES(__i), SPECIES(__k))] + __A_k2k3k1[FIELDS_FLATTEN(SPECIES(__j), SPECIES(__k), SPECIES(__i))]) / 2.0;
                     __tpf_1 += __a*__a * __Hsq * (__A_k3k1k2[FIELDS_FLATTEN(SPECIES(__k), SPECIES(__i), SPECIES(__j))] + __A_k3k2k1[FIELDS_FLATTEN(SPECIES(__k), SPECIES(__j), SPECIES(__i))]) / 2.0;

                __tpf = __prefactor * __mom_factor_1 * __tpf_1 / __kprod3;

                // this prefactor has dimension 3
                auto __mom_factor_2 = __momentum_k;

                // these components are dimension 3
                // note factor of 2 compared to analytic calculation, from symmetrization over beta, gamma
                auto __tpf_2  = - (SPECIES(__j) == SPECIES(__k) ? __fields[MOMENTUM(__i)] : 0.0) * (__Ksq + __k1*__k2*__k3/__kt) * __k2dotk3 / (2.0*__Mp*__Mp);
                     __tpf_2 += - (SPECIES(__i) == SPECIES(__k) ? __fields[MOMENTUM(__j)] : 0.0) * (__Ksq + __k1*__k2*__k3/__kt) * __k1dotk3 / (2.0*__Mp*__Mp);
                     __tpf_2 += - (SPECIES(__i) == SPECIES(__j) ? __fields[MOMENTUM(__k)] : 0.0) * (__Ksq + __k1*__k2*__k3/__kt) * __k1dotk2 / (2.0*__Mp*__Mp);

                // these components are dimension 3
                     __tpf_2 += (__C_k1k2k3[FIELDS_FLATTEN(SPECIES(__i), SPECIES(__j), SPECIES(__k))] + __C_k2k1k3[FIELDS_FLATTEN(SPECIES(__j), SPECIES(__i), SPECIES(__k))])*__k1*__k1*__k2*__k2*(1.0+__k3/__kt) / 2.0;
                     __tpf_2 += (__C_k1k3k2[FIELDS_FLATTEN(SPECIES(__i), SPECIES(__k), SPECIES(__j))] + __C_k3k1k2[FIELDS_FLATTEN(SPECIES(__k), SPECIES(__i), SPECIES(__j))])*__k1*__k1*__k3*__k3*(1.0+__k2/__kt) / 2.0;
                     __tpf_2 += (__C_k2k3k1[FIELDS_FLATTEN(SPECIES(__j), SPECIES(__k), SPECIES(__i))] + __C_k3k2k1[FIELDS_FLATTEN(SPECIES(__k), SPECIES(__j), SPECIES(__i))])*__k2*__k2*__k3*__k3*(1.0+__k1/__kt) / 2.0;

                // these components are dimension 3
                     __tpf_2 += (__B_k2k3k1[FIELDS_FLATTEN(SPECIES(__j), SPECIES(__k), SPECIES(__i))] + __B_k3k2k1[FIELDS_FLATTEN(SPECIES(__k), SPECIES(__j), SPECIES(__i))])*__k1*__k1*__k2*__k3 / 2.0;
                     __tpf_2 += (__B_k1k3k2[FIELDS_FLATTEN(SPECIES(__i), SPECIES(__k), SPECIES(__j))] + __B_k3k1k2[FIELDS_FLATTEN(SPECIES(__k), SPECIES(__i), SPECIES(__j))])*__k2*__k2*__k1*__k3 / 2.0;
                     __tpf_2 += (__B_k1k2k3[FIELDS_FLATTEN(SPECIES(__i), SPECIES(__j), SPECIES(__k))] + __B_k2k1k3[FIELDS_FLATTEN(SPECIES(__j), SPECIES(__i), SPECIES(__k))])*__k3*__k3*__k1*__k2 / 2.0;

                // these components are dimension 3
                     __tpf_2 += - __a*__a * __Hsq * (__A_k1k2k3[FIELDS_FLATTEN(SPECIES(__i), SPECIES(__j), SPECIES(__k))] + __A_k1k3k2[FIELDS_FLATTEN(SPECIES(__i), SPECIES(__k), SPECIES(__j))]) * (__Ksq - __k1*__k2*__k3/__kt) / 2.0;
                     __tpf_2 += - __a*__a * __Hsq * (__A_k2k1k3[FIELDS_FLATTEN(SPECIES(__j), SPECIES(__i), SPECIES(__k))] + __A_k2k3k1[FIELDS_FLATTEN(SPECIES(__j), SPECIES(__k), SPECIES(__i))]) * (__Ksq - __k1*__k2*__k3/__kt) / 2.0;
                     __tpf_2 += - __a*__a * __Hsq * (__A_k3k1k2[FIELDS_FLATTEN(SPECIES(__k), SPECIES(__i), SPECIES(__j))] + __A_k3k2k1[FIELDS_FLATTEN(SPECIES(__k), SPECIES(__j), SPECIES(__i))]) * (__Ksq - __k1*__k2*__k3/__kt) / 2.0;

                __tpf += __prefactor * __mom_factor_2 * __tpf_2 / __kprod3;

                break;
              }

            case 1:   // field-momentum-momentum correlation function
              {
                assert(total_fields == 1);

                // this prefactor has dimension 3
                auto __prefactor = (__k1*__k2*__k3) / (__kt * __Hsq * __a*__a*__a*__a*__a*__a);

                // component of prefactor that should not be symmetrized
                auto __mom_factor1 = (IS_FIELD(__i) ? __k2*__k2*__k3*__k3*__k1 : 0.0) + (IS_FIELD(__j) ? __k1*__k1*__k3*__k3*__k2 : 0.0) + (IS_FIELD(__k) ? __k1*__k1*__k2*__k2*__k3 : 0.0);

                // these components have dimension 1
                // note factor of 2 compared to analytic calculation, from symmetrization over beta, gamma
                auto __tpf_1  = - (SPECIES(__j) == SPECIES(__k) ? __fields[MOMENTUM(__i)] : 0.0) * __k2dotk3 / (2.0*__Mp*__Mp);
                     __tpf_1 += - (SPECIES(__i) == SPECIES(__k) ? __fields[MOMENTUM(__j)] : 0.0) * __k1dotk3 / (2.0*__Mp*__Mp);
                     __tpf_1 += - (SPECIES(__i) == SPECIES(__j) ? __fields[MOMENTUM(__k)] : 0.0) * __k1dotk2 / (2.0*__Mp*__Mp);

                // these components have dimension 1
                     __tpf_1 += (__C_k1k2k3[FIELDS_FLATTEN(SPECIES(__i), SPECIES(__j), SPECIES(__k))] + __C_k2k1k3[FIELDS_FLATTEN(SPECIES(__j), SPECIES(__i), SPECIES(__k))]) * __k1*__k2 / 2.0;
                     __tpf_1 += (__C_k1k3k2[FIELDS_FLATTEN(SPECIES(__i), SPECIES(__k), SPECIES(__j))] + __C_k3k1k2[FIELDS_FLATTEN(SPECIES(__k), SPECIES(__i), SPECIES(__j))]) * __k1*__k3 / 2.0;
                     __tpf_1 += (__C_k2k3k1[FIELDS_FLATTEN(SPECIES(__j), SPECIES(__k), SPECIES(__i))] + __C_k3k2k1[FIELDS_FLATTEN(SPECIES(__k), SPECIES(__j), SPECIES(__i))]) * __k2*__k3 / 2.0;

                // these components have dimension 1
                     __tpf_1 += - __a*__a * __Hsq * (__A_k1k2k3[FIELDS_FLATTEN(SPECIES(__i), SPECIES(__j), SPECIES(__k))] + __A_k1k3k2[FIELDS_FLATTEN(SPECIES(__i), SPECIES(__k), SPECIES(__j))]) / 2.0;
                     __tpf_1 += - __a*__a * __Hsq * (__A_k2k1k3[FIELDS_FLATTEN(SPECIES(__j), SPECIES(__i), SPECIES(__k))] + __A_k2k3k1[FIELDS_FLATTEN(SPECIES(__j), SPECIES(__k), SPECIES(__i))]) / 2.0;
                     __tpf_1 += - __a*__a * __Hsq * (__A_k3k1k2[FIELDS_FLATTEN(SPECIES(__k), SPECIES(__i), SPECIES(__j))] + __A_k3k2k1[FIELDS_FLATTEN(SPECIES(__k), SPECIES(__j), SPECIES(__i))]) / 2.0;

                // these components have dimension 1
                     __tpf_1 += - __a*__a * __Hsq * (__B_k1k2k3[FIELDS_FLATTEN(SPECIES(__i), SPECIES(__j), SPECIES(__k))] + __B_k2k1k3[FIELDS_FLATTEN(SPECIES(__j), SPECIES(__i), SPECIES(__k))]) * (__k1+__k2)*__k3 / (__k1*__k2) / 2.0;
                     __tpf_1 += - __a*__a * __Hsq * (__B_k1k3k2[FIELDS_FLATTEN(SPECIES(__i), SPECIES(__k), SPECIES(__j))] + __B_k3k1k2[FIELDS_FLATTEN(SPECIES(__k), SPECIES(__i), SPECIES(__j))]) * (__k1+__k3)*__k2 / (__k1*__k3) / 2.0;
                     __tpf_1 += - __a*__a * __Hsq * (__B_k2k3k1[FIELDS_FLATTEN(SPECIES(__j), SPECIES(__k), SPECIES(__i))] + __B_k3k2k1[FIELDS_FLATTEN(SPECIES(__k), SPECIES(__j), SPECIES(__i))]) * (__k2+__k3)*__k1 / (__k2*__k3) / 2.0;

                __tpf = __prefactor * __mom_factor1 * __tpf_1 / __kprod3;

                // prefactor has dimension 4
                auto __mom_factor2 = (IS_FIELD(__i) ? __k2*__k2*__k3*__k3 : 0.0) + (IS_FIELD(__j) ? __k1*__k1*__k3*__k3 : 0.0) + (IS_FIELD(__k) ? __k1*__k1*__k2*__k2 : 0.0);

                // these components have dimension 2
                auto __tpf_2  = __a*__a * __Hsq * (__B_k1k2k3[FIELDS_FLATTEN(SPECIES(__i), SPECIES(__j), SPECIES(__k))] + __B_k2k1k3[FIELDS_FLATTEN(SPECIES(__j), SPECIES(__i), SPECIES(__k))]) * __k3 / 2.0;
                     __tpf_2 += __a*__a * __Hsq * (__B_k1k3k2[FIELDS_FLATTEN(SPECIES(__i), SPECIES(__k), SPECIES(__j))] + __B_k3k1k2[FIELDS_FLATTEN(SPECIES(__k), SPECIES(__i), SPECIES(__j))]) * __k2 / 2.0;
                     __tpf_2 += __a*__a * __Hsq * (__B_k2k3k1[FIELDS_FLATTEN(SPECIES(__j), SPECIES(__k), SPECIES(__i))] + __B_k3k2k1[FIELDS_FLATTEN(SPECIES(__k), SPECIES(__j), SPECIES(__i))]) * __k1 / 2.0;

                __tpf += __prefactor * __mom_factor2 * __tpf_2 / __kprod3;

                break;
              }

            case 0:   // momentum-momentum-momentum correlation function
              {
                assert(total_fields == 0);

                // prefactor is dimension 3
                auto __prefactor = (__k1*__k1) * (__k2*__k2) * (__k3*__k3) / (__kt * __Hsq * __a*__a*__a*__a*__a*__a);

                // these components are dimension 3
                // note factor of 2 compared to analytic calculation, from symmetrization over beta, gamma
                __tpf  = (SPECIES(__j) == SPECIES(__k) ? __fields[MOMENTUM(__i)] : 0.0) * (__Ksq + __k1*__k2*__k3/__kt) * __k2dotk3 / (2.0*__Mp*__Mp);
                __tpf += (SPECIES(__i) == SPECIES(__k) ? __fields[MOMENTUM(__j)] : 0.0) * (__Ksq + __k1*__k2*__k3/__kt) * __k1dotk3 / (2.0*__Mp*__Mp);
                __tpf += (SPECIES(__i) == SPECIES(__j) ? __fields[MOMENTUM(__k)] : 0.0) * (__Ksq + __k1*__k2*__k3/__kt) * __k1dotk2 / (2.0*__Mp*__Mp);

                // these components are dimension 3
                __tpf += - (__C_k1k2k3[FIELDS_FLATTEN(SPECIES(__i), SPECIES(__j), SPECIES(__k))] + __C_k2k1k3[FIELDS_FLATTEN(SPECIES(__j), SPECIES(__i), SPECIES(__k))])*__k1*__k1*__k2*__k2*(1.0+__k3/__kt) / 2.0;
                __tpf += - (__C_k1k3k2[FIELDS_FLATTEN(SPECIES(__i), SPECIES(__k), SPECIES(__j))] + __C_k3k1k2[FIELDS_FLATTEN(SPECIES(__k), SPECIES(__i), SPECIES(__j))])*__k1*__k1*__k3*__k3*(1.0+__k2/__kt) / 2.0;
                __tpf += - (__C_k2k3k1[FIELDS_FLATTEN(SPECIES(__j), SPECIES(__k), SPECIES(__i))] + __C_k3k2k1[FIELDS_FLATTEN(SPECIES(__k), SPECIES(__j), SPECIES(__i))])*__k2*__k2*__k3*__k3*(1.0+__k1/__kt) / 2.0;

                // these components are dimension 3
                __tpf += - (__B_k2k3k1[FIELDS_FLATTEN(SPECIES(__j), SPECIES(__k), SPECIES(__i))] + __B_k3k2k1[FIELDS_FLATTEN(SPECIES(__k), SPECIES(__j), SPECIES(__i))])*__k1*__k1*__k2*__k3/2.0;
                __tpf += - (__B_k1k3k2[FIELDS_FLATTEN(SPECIES(__i), SPECIES(__k), SPECIES(__j))] + __B_k3k1k2[FIELDS_FLATTEN(SPECIES(__k), SPECIES(__i), SPECIES(__j))])*__k2*__k2*__k1*__k3/2.0;
                __tpf += - (__B_k1k2k3[FIELDS_FLATTEN(SPECIES(__i), SPECIES(__j), SPECIES(__k))] + __B_k2k1k3[FIELDS_FLATTEN(SPECIES(__j), SPECIES(__i), SPECIES(__k))])*__k3*__k3*__k1*__k2/2.0;

                // these components are dimension 3
                __tpf += __a*__a * __Hsq * (__A_k1k2k3[FIELDS_FLATTEN(SPECIES(__i), SPECIES(__j), SPECIES(__k))] + __A_k1k2k3[FIELDS_FLATTEN(SPECIES(__i), SPECIES(__k), SPECIES(__j))]) * (__Ksq - __k1*__k2*__k3/__kt) / 2.0;
                __tpf += __a*__a * __Hsq * (__A_k2k1k3[FIELDS_FLATTEN(SPECIES(__j), SPECIES(__i), SPECIES(__k))] + __A_k2k3k1[FIELDS_FLATTEN(SPECIES(__j), SPECIES(__k), SPECIES(__i))]) * (__Ksq - __k1*__k2*__k3/__kt) / 2.0;
                __tpf += __a*__a * __Hsq * (__A_k3k1k2[FIELDS_FLATTEN(SPECIES(__k), SPECIES(__i), SPECIES(__j))] + __A_k3k2k1[FIELDS_FLATTEN(SPECIES(__k), SPECIES(__j), SPECIES(__i))]) * (__Ksq - __k1*__k2*__k3/__kt) / 2.0;

                __tpf *= __prefactor / __kprod3;

                break;
              }

            default:
              assert(false);
          }

        // return value, rescaled to give dimensionless correlation function
        return(__tpf * __k_norm*__k_norm*__k_norm*__k_norm*__k_norm*__k_norm);
      }


    // CALCULATE GAUGE TRANSFORMATIONS


    template <typename number>
    void quartic<number>::compute_gauge_xfm_1(const twopf_db_task<number>* __task,
                                             const flattened_tensor<number>& __state,
                                             flattened_tensor<number>& __dN)
      {
        DEFINE_INDEX_TOOLS
//         release resources 
        const auto& __pvector = __task->get_params().get_vector();
        __raw_params[0] = __pvector[0];

        const auto __Mp = __task->get_params().get_Mp();

//         parameters resource set to '__raw_params' 
//         coordinates resource set to '__state' 

// BEGIN TEMPORARY POOL (sequence=9) 
const auto _t_9_0 = __Mp;
const auto _t_9_1 = -2.0;
const auto _t_9_2 = 1.0/(_t_9_0*_t_9_0);
const auto _t_9_3 = __state[FLATTEN(1)];
const auto _t_9_4 = 2.0;
const auto _t_9_5 = (_t_9_3*_t_9_3);
const auto _t_9_6 = (1.0/2.0);
const auto _InternalEps = _t_9_2*_t_9_5*_t_9_6;
const auto _t_9_8 = -1.0;
const auto _t_9_7 = _InternalEps;
const auto _t_9_9 = 1.0/_t_9_7;
const auto _t_9_10 = -(1.0/2.0);
const auto _t_9_11 = _t_9_2*_t_9_9*_t_9_3*_t_9_10;
const auto _t_9_12 = 0.0;
        // END TEMPORARY POOL (sequence=9) 

        __dN[FLATTEN(0)] = _t_9_11;
        __dN[FLATTEN(1)] = _t_9_12;
      }


    template <typename number>
    void quartic<number>::compute_gauge_xfm_2(const twopf_db_task<number>* __task,
                                             const flattened_tensor<number>& __state,
                                             double __k, double __k1, double __k2, double __N,
                                             flattened_tensor<number>& __ddN)
      {
        DEFINE_INDEX_TOOLS
//         release resources 
        const auto& __pvector = __task->get_params().get_vector();
        __raw_params[0] = __pvector[0];

        const auto __Mp = __task->get_params().get_Mp();
        const auto __a = std::exp(__N - __task->get_N_horizon_crossing() + __task->get_astar_normalization());

//         parameters resource set to '__raw_params' 
//         coordinates resource set to '__state' 
//         ENDIF !fast 

// BEGIN TEMPORARY POOL (sequence=10) 
const auto _t_10_0 = __raw_params[0];
const auto _t_10_1 = __state[FLATTEN(0)];
const auto _t_10_2 = 4.0;
const auto _t_10_3 = (_t_10_1*_t_10_1*_t_10_1*_t_10_1);
const auto _t_10_4 = __Mp;
const auto _t_10_5 = -2.0;
const auto _t_10_6 = 1.0/(_t_10_4*_t_10_4);
const auto _t_10_7 = __state[FLATTEN(1)];
const auto _t_10_8 = 2.0;
const auto _t_10_9 = (_t_10_7*_t_10_7);
const auto _t_10_10 = _t_10_6*_t_10_9;
const auto _t_10_11 = -6.0;
const auto _t_10_12 = _t_10_10+_t_10_11;
const auto _t_10_13 = -1.0;
const auto _t_10_14 = 1.0/(_t_10_12);
const auto _t_10_15 = -(1.0/2.0);
const auto _InternalHsq = _t_10_0*_t_10_3*_t_10_6*_t_10_14*_t_10_15;
const auto _t_10_16 = (1.0/2.0);
const auto _InternalEps = _t_10_6*_t_10_9*_t_10_16;
const auto _t_10_17 = 3.0;
const auto _t_10_18 = (_t_10_1*_t_10_1*_t_10_1);
const auto _t_10_19 = _InternalHsq;
const auto _t_10_20 = 1.0/_t_10_19;
const auto _InternalZ2p = _t_10_0*_t_10_18*_t_10_6*_t_10_20*_t_10_7;
const auto _t_10_21 = -4.0;
const auto _t_10_22 = 1.0/(_t_10_4*_t_10_4*_t_10_4*_t_10_4);
const auto _t_10_23 = _InternalEps;
const auto _t_10_24 = 1.0/_t_10_23;
const auto _t_10_25 = 1.0/(_t_10_23*_t_10_23);
const auto _t_10_26 = _t_10_0*_t_10_18*_t_10_6*_t_10_25*_t_10_20*_t_10_7;
const auto _t_10_27 = 6.0;
const auto _t_10_28 = _t_10_24*_t_10_27;
const auto _t_10_29 = _t_10_26+_t_10_28+_t_10_5;
const auto _t_10_30 = (1.0/4.0);
const auto _t_10_31 = _t_10_22*_t_10_24*(_t_10_29)*_t_10_9*_t_10_30;
const auto _t_10_32 = _t_10_22*_t_10_25*_t_10_9*_t_10_16;
const auto _t_10_33 = __k1;
const auto _t_10_34 = (_t_10_33*_t_10_33);
const auto _t_10_35 = __k;
const auto _t_10_36 = 1.0/(_t_10_35*_t_10_35);
const auto _t_10_37 = _t_10_34*_t_10_6*_t_10_24*_t_10_36*_t_10_15;
const auto _t_10_38 = __k2;
const auto _t_10_39 = (_t_10_38*_t_10_38);
const auto _t_10_40 = (_t_10_35*_t_10_35);
const auto _t_10_41 = _t_10_40*_t_10_13;
const auto _t_10_42 = _t_10_34+_t_10_39+_t_10_41;
const auto _t_10_43 = _t_10_6*_t_10_24*(_t_10_42)*_t_10_36*_t_10_30;
const auto _t_10_44 = _t_10_32+_t_10_37+_t_10_43;
const auto _t_10_45 = _t_10_6*_t_10_24*_t_10_39*_t_10_36*_t_10_15;
const auto _t_10_46 = _t_10_32+_t_10_43+_t_10_45;
const auto _t_10_47 = 0.0;
        // END TEMPORARY POOL (sequence=10) 

        __ddN[FLATTEN(0,0)] = _t_10_31;
        __ddN[FLATTEN(0,1)] = _t_10_44;
        __ddN[FLATTEN(1,0)] = _t_10_46;
        __ddN[FLATTEN(1,1)] = _t_10_47;
      }


    // CALCULATE TENSOR QUANTITIES


    template <typename number>
    void quartic<number>::u2(const twopf_db_task<number>* __task,
                            const flattened_tensor<number>& __fields, double __k, double __N,
                            flattened_tensor<number>& __u2)
      {
        DEFINE_INDEX_TOOLS
//         release resources 
        const auto& __pvector = __task->get_params().get_vector();
        __raw_params[0] = __pvector[0];

        const auto __Mp = __task->get_params().get_Mp();
        const auto __a = std::exp(__N - __task->get_N_horizon_crossing() + __task->get_astar_normalization());

//         parameters resource set to '__raw_params' 
//         coordinates resource set to '__fields' 
//         ENDIF !fast 

// BEGIN TEMPORARY POOL (sequence=11) 
const auto _t_11_0 = __Mp;
const auto _t_11_1 = -2.0;
const auto _t_11_2 = 1.0/(_t_11_0*_t_11_0);
const auto _t_11_3 = __fields[FLATTEN(1)];
const auto _t_11_4 = 2.0;
const auto _t_11_5 = (_t_11_3*_t_11_3);
const auto _t_11_6 = _t_11_2*_t_11_5;
const auto _t_11_7 = -6.0;
const auto _t_11_8 = _t_11_6+_t_11_7;
const auto _t_11_9 = -1.0;
const auto _t_11_10 = 1.0/(_t_11_8);
const auto _t_11_11 = __raw_params[0];
const auto _t_11_12 = __fields[FLATTEN(0)];
const auto _t_11_13 = 4.0;
const auto _t_11_14 = (_t_11_12*_t_11_12*_t_11_12*_t_11_12);
const auto _t_11_15 = -(1.0/2.0);
const auto _InternalHsq = _t_11_10*_t_11_11*_t_11_2*_t_11_14*_t_11_15;
const auto _t_11_16 = (1.0/2.0);
const auto _InternalEps = _t_11_2*_t_11_5*_t_11_16;
const auto _t_11_17 = 0.0;
const auto _t_11_18 = 1.0;
const auto _t_11_19 = __a;
const auto _t_11_20 = 1.0/(_t_11_19*_t_11_19);
const auto _t_11_21 = _InternalHsq;
const auto _t_11_22 = 1.0/_t_11_21;
const auto _t_11_23 = __k;
const auto _t_11_24 = (_t_11_23*_t_11_23);
const auto _t_11_25 = _t_11_20*_t_11_22*_t_11_24*_t_11_9;
const auto _t_11_26 = 3.0;
const auto _t_11_27 = (_t_11_12*_t_11_12*_t_11_12);
const auto _t_11_28 = _t_11_11*_t_11_2*_t_11_27*_t_11_22*_t_11_3*_t_11_1;
const auto _t_11_29 = (_t_11_12*_t_11_12);
const auto _t_11_30 = -3.0;
const auto _t_11_31 = _t_11_11*_t_11_29*_t_11_22*_t_11_30;
const auto _t_11_32 = _InternalEps;
const auto _t_11_33 = _t_11_32+_t_11_30;
const auto _t_11_34 = _t_11_2*_t_11_5*(_t_11_33);
const auto _t_11_35 = _t_11_25+_t_11_28+_t_11_31+_t_11_34;
        // END TEMPORARY POOL (sequence=11) 

        __u2[FLATTEN(0,0)] = _t_11_17;
        __u2[FLATTEN(0,1)] = _t_11_18;
        __u2[FLATTEN(1,0)] = _t_11_35;
        __u2[FLATTEN(1,1)] = _t_11_33;
      }


    template <typename number>
    void quartic<number>::u3(const twopf_db_task<number>* __task,
                            const flattened_tensor<number>& __fields, double __k1, double __k2, double __k3, double __N,
                            flattened_tensor<number>& __u3)
      {
        DEFINE_INDEX_TOOLS
//         release resources 
        const auto& __pvector = __task->get_params().get_vector();
        __raw_params[0] = __pvector[0];

        const auto __Mp = __task->get_params().get_Mp();
        const auto __a = std::exp(__N - __task->get_N_horizon_crossing() + __task->get_astar_normalization());

//         parameters resource set to '__raw_params' 
//         coordinates resource set to '__fields' 
//         ENDIF !fast 

// BEGIN TEMPORARY POOL (sequence=12) 
const auto _t_12_0 = __Mp;
const auto _t_12_1 = -2.0;
const auto _t_12_2 = 1.0/(_t_12_0*_t_12_0);
const auto _t_12_3 = __fields[FLATTEN(1)];
const auto _t_12_4 = 2.0;
const auto _t_12_5 = (_t_12_3*_t_12_3);
const auto _t_12_6 = _t_12_2*_t_12_5;
const auto _t_12_7 = -6.0;
const auto _t_12_8 = _t_12_6+_t_12_7;
const auto _t_12_9 = -1.0;
const auto _t_12_10 = 1.0/(_t_12_8);
const auto _t_12_11 = __raw_params[0];
const auto _t_12_12 = __fields[FLATTEN(0)];
const auto _t_12_13 = 4.0;
const auto _t_12_14 = (_t_12_12*_t_12_12*_t_12_12*_t_12_12);
const auto _t_12_15 = -(1.0/2.0);
const auto _InternalHsq = _t_12_10*_t_12_11*_t_12_2*_t_12_14*_t_12_15;
const auto _t_12_16 = (1.0/2.0);
const auto _InternalEps = _t_12_2*_t_12_5*_t_12_16;
const auto _t_12_17 = -4.0;
const auto _t_12_18 = 1.0/(_t_12_0*_t_12_0*_t_12_0*_t_12_0);
const auto _t_12_19 = 3.0;
const auto _t_12_20 = (_t_12_3*_t_12_3*_t_12_3);
const auto _t_12_21 = -(1.0/4.0);
const auto _t_12_22 = _t_12_18*_t_12_20*_t_12_21;
const auto _t_12_23 = __k1;
const auto _t_12_24 = (_t_12_23*_t_12_23);
const auto _t_12_25 = __k3;
const auto _t_12_26 = (_t_12_25*_t_12_25);
const auto _t_12_27 = _t_12_26*_t_12_9;
const auto _t_12_28 = __k2;
const auto _t_12_29 = (_t_12_28*_t_12_28);
const auto _t_12_30 = _t_12_29*_t_12_9;
const auto _t_12_31 = _t_12_24+_t_12_27+_t_12_30;
const auto _t_12_32 = _InternalEps;
const auto _t_12_33 = -3.0;
const auto _t_12_34 = _t_12_32+_t_12_33;
const auto _t_12_35 = _t_12_3*(_t_12_34);
const auto _t_12_36 = (_t_12_12*_t_12_12*_t_12_12);
const auto _t_12_37 = _InternalHsq;
const auto _t_12_38 = 1.0/_t_12_37;
const auto _t_12_39 = _t_12_11*_t_12_36*_t_12_38*_t_12_9;
const auto _t_12_40 = _t_12_35+_t_12_39;
const auto _t_12_41 = 1.0/(_t_12_25*_t_12_25);
const auto _t_12_42 = (1.0/4.0);
const auto _t_12_43 = _t_12_2*(_t_12_31)*(_t_12_40)*_t_12_41*_t_12_42;
const auto _t_12_44 = 1.0/(_t_12_28*_t_12_28);
const auto _t_12_45 = _t_12_2*_t_12_44*(_t_12_31)*(_t_12_40)*_t_12_42;
const auto _t_12_46 = 1.0/(_t_12_23*_t_12_23);
const auto _t_12_47 = _t_12_24+_t_12_27+_t_12_29;
const auto _t_12_48 = ((_t_12_47)*(_t_12_47));
const auto _t_12_49 = _t_12_46*_t_12_48*_t_12_44;
const auto _t_12_50 = _t_12_49+_t_12_17;
const auto _t_12_51 = -(1.0/32.0);
const auto _t_12_52 = _t_12_18*(_t_12_40)*(_t_12_50)*_t_12_5*_t_12_51;
const auto _t_12_53 = _t_12_24+_t_12_26+_t_12_30;
const auto _t_12_54 = ((_t_12_53)*(_t_12_53));
const auto _t_12_55 = _t_12_54*_t_12_46*_t_12_41;
const auto _t_12_56 = _t_12_55+_t_12_17;
const auto _t_12_57 = _t_12_18*(_t_12_40)*(_t_12_56)*_t_12_5*_t_12_51;
const auto _t_12_58 = _t_12_22+_t_12_43+_t_12_45+_t_12_52+_t_12_57;
const auto _t_12_59 = _t_12_2*(_t_12_31)*_t_12_41*_t_12_3*_t_12_21;
const auto _t_12_60 = _t_12_46*(_t_12_47)*_t_12_2*_t_12_3*_t_12_42;
const auto _t_12_61 = _t_12_2*_t_12_3*_t_12_16;
const auto _t_12_62 = (1.0/32.0);
const auto _t_12_63 = _t_12_18*(_t_12_56)*_t_12_20*_t_12_62;
const auto _t_12_64 = _t_12_59+_t_12_60+_t_12_61+_t_12_63;
const auto _t_12_65 = _t_12_2*_t_12_44*(_t_12_31)*_t_12_3*_t_12_21;
const auto _t_12_66 = (_t_12_53)*_t_12_46*_t_12_2*_t_12_3*_t_12_42;
const auto _t_12_67 = _t_12_18*(_t_12_50)*_t_12_20*_t_12_62;
const auto _t_12_68 = _t_12_65+_t_12_66+_t_12_67+_t_12_61;
const auto _t_12_69 = 0.0;
const auto _t_12_70 = ((_t_12_40)*(_t_12_40));
const auto _t_12_71 = _t_12_18*_t_12_70*(_t_12_56)*_t_12_3*_t_12_51;
const auto _t_12_72 = __a;
const auto _t_12_73 = 1.0/(_t_12_72*_t_12_72);
const auto _t_12_74 = _t_12_73*(_t_12_47)*_t_12_2*_t_12_38*_t_12_3*_t_12_21;
const auto _t_12_75 = _t_12_18*_t_12_70*(_t_12_50)*_t_12_3*_t_12_51;
const auto _t_12_76 = (3.0/4.0);
const auto _t_12_77 = _t_12_18*(_t_12_40)*_t_12_5*_t_12_76;
const auto _t_12_78 = _t_12_11*_t_12_12*_t_12_38*_t_12_7;
const auto _t_12_79 = ((_t_12_31)*(_t_12_31));
const auto _t_12_80 = _t_12_44*_t_12_79*_t_12_41;
const auto _t_12_81 = _t_12_80+_t_12_17;
const auto _t_12_82 = _t_12_18*_t_12_70*_t_12_3*(_t_12_81)*_t_12_51;
const auto _t_12_83 = (_t_12_53)*_t_12_73*_t_12_2*_t_12_38*_t_12_3*_t_12_21;
const auto _t_12_84 = -(3.0/4.0);
const auto _t_12_85 = _t_12_18*_t_12_20*(_t_12_34)*_t_12_84;
const auto _t_12_86 = _t_12_73*_t_12_2*(_t_12_31)*_t_12_38*_t_12_3*_t_12_42;
const auto _t_12_87 = (_t_12_12*_t_12_12);
const auto _t_12_88 = -(9.0/2.0);
const auto _t_12_89 = _t_12_11*_t_12_2*_t_12_87*_t_12_38*_t_12_3*_t_12_88;
const auto _t_12_90 = _t_12_71+_t_12_74+_t_12_75+_t_12_77+_t_12_78+_t_12_82+_t_12_83+_t_12_85+_t_12_86+_t_12_89;
const auto _t_12_91 = _t_12_18*_t_12_20*_t_12_42;
const auto _t_12_92 = (_t_12_47)*_t_12_2*_t_12_44*(_t_12_40)*_t_12_42;
const auto _t_12_93 = _t_12_18*(_t_12_40)*_t_12_5*(_t_12_81)*_t_12_62;
const auto _t_12_94 = _t_12_18*(_t_12_40)*(_t_12_56)*_t_12_5*_t_12_62;
const auto _t_12_95 = _t_12_46*(_t_12_47)*_t_12_2*(_t_12_40)*_t_12_42;
const auto _t_12_96 = _t_12_91+_t_12_92+_t_12_93+_t_12_94+_t_12_95;
const auto _t_12_97 = _t_12_18*(_t_12_40)*(_t_12_50)*_t_12_5*_t_12_62;
const auto _t_12_98 = (_t_12_53)*_t_12_2*(_t_12_40)*_t_12_41*_t_12_42;
const auto _t_12_99 = (_t_12_53)*_t_12_46*_t_12_2*(_t_12_40)*_t_12_42;
const auto _t_12_100 = _t_12_91+_t_12_97+_t_12_93+_t_12_98+_t_12_99;
const auto _t_12_101 = (_t_12_53)*_t_12_2*_t_12_41*_t_12_3*_t_12_21;
const auto _t_12_102 = _t_12_18*_t_12_20*(_t_12_81)*_t_12_51;
const auto _t_12_103 = _t_12_2*_t_12_3*_t_12_15;
const auto _t_12_104 = (_t_12_47)*_t_12_2*_t_12_44*_t_12_3*_t_12_21;
const auto _t_12_105 = _t_12_101+_t_12_102+_t_12_103+_t_12_104;
        // END TEMPORARY POOL (sequence=12) 

        __u3[FLATTEN(0,0,0)] = _t_12_58;
        __u3[FLATTEN(0,0,1)] = _t_12_64;
        __u3[FLATTEN(0,1,0)] = _t_12_68;
        __u3[FLATTEN(0,1,1)] = _t_12_69;
        __u3[FLATTEN(1,0,0)] = _t_12_90;
        __u3[FLATTEN(1,0,1)] = _t_12_96;
        __u3[FLATTEN(1,1,0)] = _t_12_100;
        __u3[FLATTEN(1,1,1)] = _t_12_105;
      }


    template <typename number>
    void quartic<number>::A(const twopf_db_task<number>* __task,
                           const flattened_tensor<number>& __fields, double __k1, double __k2, double __k3, double __N,
                           flattened_tensor<number>& __A)
      {
        DEFINE_INDEX_TOOLS
//         release resources 
        const auto& __pvector = __task->get_params().get_vector();
        __raw_params[0] = __pvector[0];

        const auto __Mp = __task->get_params().get_Mp();
        const auto __a = std::exp(__N - __task->get_N_horizon_crossing() + __task->get_astar_normalization());

//         parameters resource set to '__raw_params' 
//         coordinates resource set to '__fields' 
//         ENDIF !fast 

// BEGIN TEMPORARY POOL (sequence=13) 
const auto _t_13_0 = __Mp;
const auto _t_13_1 = -2.0;
const auto _t_13_2 = 1.0/(_t_13_0*_t_13_0);
const auto _t_13_3 = __fields[FLATTEN(1)];
const auto _t_13_4 = 2.0;
const auto _t_13_5 = (_t_13_3*_t_13_3);
const auto _t_13_6 = _t_13_2*_t_13_5;
const auto _t_13_7 = -6.0;
const auto _t_13_8 = _t_13_6+_t_13_7;
const auto _t_13_9 = -1.0;
const auto _t_13_10 = 1.0/(_t_13_8);
const auto _t_13_11 = __raw_params[0];
const auto _t_13_12 = __fields[FLATTEN(0)];
const auto _t_13_13 = 4.0;
const auto _t_13_14 = (_t_13_12*_t_13_12*_t_13_12*_t_13_12);
const auto _t_13_15 = -(1.0/2.0);
const auto _InternalHsq = _t_13_10*_t_13_11*_t_13_2*_t_13_14*_t_13_15;
const auto _t_13_16 = (1.0/2.0);
const auto _InternalEps = _t_13_2*_t_13_5*_t_13_16;
const auto _t_13_17 = -4.0;
const auto _t_13_18 = 1.0/(_t_13_0*_t_13_0*_t_13_0*_t_13_0);
const auto _t_13_19 = _InternalEps;
const auto _t_13_20 = -3.0;
const auto _t_13_21 = _t_13_19+_t_13_20;
const auto _t_13_22 = _t_13_3*(_t_13_21);
const auto _t_13_23 = 3.0;
const auto _t_13_24 = (_t_13_12*_t_13_12*_t_13_12);
const auto _t_13_25 = _InternalHsq;
const auto _t_13_26 = 1.0/_t_13_25;
const auto _t_13_27 = _t_13_11*_t_13_24*_t_13_26*_t_13_9;
const auto _t_13_28 = _t_13_22+_t_13_27;
const auto _t_13_29 = ((_t_13_28)*(_t_13_28));
const auto _t_13_30 = __k1;
const auto _t_13_31 = (_t_13_30*_t_13_30);
const auto _t_13_32 = __k3;
const auto _t_13_33 = (_t_13_32*_t_13_32);
const auto _t_13_34 = __k2;
const auto _t_13_35 = (_t_13_34*_t_13_34);
const auto _t_13_36 = _t_13_35*_t_13_9;
const auto _t_13_37 = _t_13_31+_t_13_33+_t_13_36;
const auto _t_13_38 = ((_t_13_37)*(_t_13_37));
const auto _t_13_39 = 1.0/(_t_13_30*_t_13_30);
const auto _t_13_40 = 1.0/(_t_13_32*_t_13_32);
const auto _t_13_41 = _t_13_38*_t_13_39*_t_13_40;
const auto _t_13_42 = _t_13_41+_t_13_17;
const auto _t_13_43 = -(1.0/96.0);
const auto _t_13_44 = _t_13_18*_t_13_29*(_t_13_42)*_t_13_3*_t_13_43;
const auto _t_13_45 = __a;
const auto _t_13_46 = 1.0/(_t_13_45*_t_13_45);
const auto _t_13_47 = _t_13_33*_t_13_9;
const auto _t_13_48 = _t_13_31+_t_13_47+_t_13_35;
const auto _t_13_49 = -(1.0/12.0);
const auto _t_13_50 = _t_13_46*(_t_13_48)*_t_13_2*_t_13_26*_t_13_3*_t_13_49;
const auto _t_13_51 = ((_t_13_48)*(_t_13_48));
const auto _t_13_52 = 1.0/(_t_13_34*_t_13_34);
const auto _t_13_53 = _t_13_39*_t_13_51*_t_13_52;
const auto _t_13_54 = _t_13_53+_t_13_17;
const auto _t_13_55 = _t_13_18*_t_13_29*(_t_13_54)*_t_13_3*_t_13_43;
const auto _t_13_56 = (1.0/4.0);
const auto _t_13_57 = _t_13_18*(_t_13_28)*_t_13_5*_t_13_56;
const auto _t_13_58 = _t_13_11*_t_13_12*_t_13_26*_t_13_1;
const auto _t_13_59 = _t_13_31+_t_13_47+_t_13_36;
const auto _t_13_60 = ((_t_13_59)*(_t_13_59));
const auto _t_13_61 = _t_13_52*_t_13_60*_t_13_40;
const auto _t_13_62 = _t_13_61+_t_13_17;
const auto _t_13_63 = _t_13_18*_t_13_29*_t_13_3*(_t_13_62)*_t_13_43;
const auto _t_13_64 = (_t_13_37)*_t_13_46*_t_13_2*_t_13_26*_t_13_3*_t_13_49;
const auto _t_13_65 = (_t_13_3*_t_13_3*_t_13_3);
const auto _t_13_66 = -(1.0/4.0);
const auto _t_13_67 = _t_13_18*_t_13_65*(_t_13_21)*_t_13_66;
const auto _t_13_68 = (1.0/12.0);
const auto _t_13_69 = _t_13_46*_t_13_2*(_t_13_59)*_t_13_26*_t_13_3*_t_13_68;
const auto _t_13_70 = (_t_13_12*_t_13_12);
const auto _t_13_71 = -(3.0/2.0);
const auto _t_13_72 = _t_13_11*_t_13_2*_t_13_70*_t_13_26*_t_13_3*_t_13_71;
const auto _t_13_73 = _t_13_44+_t_13_50+_t_13_55+_t_13_57+_t_13_58+_t_13_63+_t_13_64+_t_13_67+_t_13_69+_t_13_72;
        // END TEMPORARY POOL (sequence=13) 

        __A[FIELDS_FLATTEN(0,0,0)] = _t_13_73;
      }


    template <typename number>
    void quartic<number>::B(const twopf_db_task<number>* __task,
                           const flattened_tensor<number>& __fields, double __k1, double __k2, double __k3, double __N,
                           flattened_tensor<number>& __B)
      {
        DEFINE_INDEX_TOOLS
//         release resources 
        const auto& __pvector = __task->get_params().get_vector();
        __raw_params[0] = __pvector[0];

        const auto __Mp = __task->get_params().get_Mp();
        const auto __a = std::exp(__N - __task->get_N_horizon_crossing() + __task->get_astar_normalization());

//         parameters resource set to '__raw_params' 
//         coordinates resource set to '__fields' 
//         ENDIF !fast 

// BEGIN TEMPORARY POOL (sequence=14) 
const auto _t_14_0 = __Mp;
const auto _t_14_1 = -2.0;
const auto _t_14_2 = 1.0/(_t_14_0*_t_14_0);
const auto _t_14_3 = __fields[FLATTEN(1)];
const auto _t_14_4 = 2.0;
const auto _t_14_5 = (_t_14_3*_t_14_3);
const auto _t_14_6 = _t_14_2*_t_14_5;
const auto _t_14_7 = -6.0;
const auto _t_14_8 = _t_14_6+_t_14_7;
const auto _t_14_9 = -1.0;
const auto _t_14_10 = 1.0/(_t_14_8);
const auto _t_14_11 = __raw_params[0];
const auto _t_14_12 = __fields[FLATTEN(0)];
const auto _t_14_13 = 4.0;
const auto _t_14_14 = (_t_14_12*_t_14_12*_t_14_12*_t_14_12);
const auto _t_14_15 = -(1.0/2.0);
const auto _InternalHsq = _t_14_10*_t_14_11*_t_14_2*_t_14_14*_t_14_15;
const auto _t_14_16 = (1.0/2.0);
const auto _InternalEps = _t_14_2*_t_14_5*_t_14_16;
const auto _t_14_17 = -4.0;
const auto _t_14_18 = 1.0/(_t_14_0*_t_14_0*_t_14_0*_t_14_0);
const auto _t_14_19 = 3.0;
const auto _t_14_20 = (_t_14_3*_t_14_3*_t_14_3);
const auto _t_14_21 = (1.0/4.0);
const auto _t_14_22 = _t_14_18*_t_14_20*_t_14_21;
const auto _t_14_23 = __k1;
const auto _t_14_24 = (_t_14_23*_t_14_23);
const auto _t_14_25 = __k3;
const auto _t_14_26 = (_t_14_25*_t_14_25);
const auto _t_14_27 = _t_14_26*_t_14_9;
const auto _t_14_28 = __k2;
const auto _t_14_29 = (_t_14_28*_t_14_28);
const auto _t_14_30 = _t_14_24+_t_14_27+_t_14_29;
const auto _t_14_31 = 1.0/(_t_14_28*_t_14_28);
const auto _t_14_32 = _InternalEps;
const auto _t_14_33 = -3.0;
const auto _t_14_34 = _t_14_32+_t_14_33;
const auto _t_14_35 = _t_14_3*(_t_14_34);
const auto _t_14_36 = (_t_14_12*_t_14_12*_t_14_12);
const auto _t_14_37 = _InternalHsq;
const auto _t_14_38 = 1.0/_t_14_37;
const auto _t_14_39 = _t_14_11*_t_14_36*_t_14_38*_t_14_9;
const auto _t_14_40 = _t_14_35+_t_14_39;
const auto _t_14_41 = (_t_14_30)*_t_14_2*_t_14_31*(_t_14_40)*_t_14_21;
const auto _t_14_42 = _t_14_29*_t_14_9;
const auto _t_14_43 = _t_14_24+_t_14_27+_t_14_42;
const auto _t_14_44 = ((_t_14_43)*(_t_14_43));
const auto _t_14_45 = 1.0/(_t_14_25*_t_14_25);
const auto _t_14_46 = _t_14_31*_t_14_44*_t_14_45;
const auto _t_14_47 = _t_14_46+_t_14_17;
const auto _t_14_48 = (1.0/32.0);
const auto _t_14_49 = _t_14_18*(_t_14_40)*_t_14_5*(_t_14_47)*_t_14_48;
const auto _t_14_50 = _t_14_24+_t_14_26+_t_14_42;
const auto _t_14_51 = ((_t_14_50)*(_t_14_50));
const auto _t_14_52 = 1.0/(_t_14_23*_t_14_23);
const auto _t_14_53 = _t_14_51*_t_14_52*_t_14_45;
const auto _t_14_54 = _t_14_53+_t_14_17;
const auto _t_14_55 = _t_14_18*(_t_14_40)*(_t_14_54)*_t_14_5*_t_14_48;
const auto _t_14_56 = _t_14_52*(_t_14_30)*_t_14_2*(_t_14_40)*_t_14_21;
const auto _t_14_57 = _t_14_22+_t_14_41+_t_14_49+_t_14_55+_t_14_56;
        // END TEMPORARY POOL (sequence=14) 

        __B[FIELDS_FLATTEN(0,0,0)] = _t_14_57;
      }


    template <typename number>
    void quartic<number>::C(const twopf_db_task<number>* __task,
                           const flattened_tensor<number>& __fields, double __k1, double __k2, double __k3, double __N,
                           flattened_tensor<number>& __C)
      {
        DEFINE_INDEX_TOOLS
//         release resources 
        const auto& __pvector = __task->get_params().get_vector();
        __raw_params[0] = __pvector[0];

        const auto __Mp = __task->get_params().get_Mp();
        const auto __a = std::exp(__N - __task->get_N_horizon_crossing() + __task->get_astar_normalization());

//         parameters resource set to '__raw_params' 
//         coordinates resource set to '__fields' 

// BEGIN TEMPORARY POOL (sequence=15) 
const auto _t_15_1 = -2.0;
const auto _t_15_0 = __Mp;
const auto _t_15_2 = 1.0/(_t_15_0*_t_15_0);
const auto _t_15_3 = __k2;
const auto _t_15_4 = 1.0/(_t_15_3*_t_15_3);
const auto _t_15_6 = 2.0;
const auto _t_15_5 = __k1;
const auto _t_15_7 = (_t_15_5*_t_15_5);
const auto _t_15_8 = __k3;
const auto _t_15_9 = (_t_15_8*_t_15_8);
const auto _t_15_10 = -1.0;
const auto _t_15_11 = _t_15_9*_t_15_10;
const auto _t_15_12 = (_t_15_3*_t_15_3);
const auto _t_15_13 = _t_15_12*_t_15_10;
const auto _t_15_14 = _t_15_7+_t_15_11+_t_15_13;
const auto _t_15_15 = __fields[FLATTEN(1)];
const auto _t_15_16 = (1.0/4.0);
const auto _t_15_17 = _t_15_2*_t_15_4*(_t_15_14)*_t_15_15*_t_15_16;
const auto _t_15_18 = _t_15_7+_t_15_9+_t_15_13;
const auto _t_15_19 = 1.0/(_t_15_5*_t_15_5);
const auto _t_15_20 = -(1.0/4.0);
const auto _t_15_21 = (_t_15_18)*_t_15_19*_t_15_2*_t_15_15*_t_15_20;
const auto _t_15_22 = -4.0;
const auto _t_15_23 = 1.0/(_t_15_0*_t_15_0*_t_15_0*_t_15_0);
const auto _t_15_24 = _t_15_7+_t_15_11+_t_15_12;
const auto _t_15_25 = ((_t_15_24)*(_t_15_24));
const auto _t_15_26 = _t_15_19*_t_15_25*_t_15_4;
const auto _t_15_27 = _t_15_26+_t_15_22;
const auto _t_15_28 = 3.0;
const auto _t_15_29 = (_t_15_15*_t_15_15*_t_15_15);
const auto _t_15_30 = -(1.0/32.0);
const auto _t_15_31 = _t_15_23*(_t_15_27)*_t_15_29*_t_15_30;
const auto _t_15_32 = -(1.0/2.0);
const auto _t_15_33 = _t_15_2*_t_15_15*_t_15_32;
const auto _t_15_34 = _t_15_17+_t_15_21+_t_15_31+_t_15_33;
        // END TEMPORARY POOL (sequence=15) 

        __C[FIELDS_FLATTEN(0,0,0)] = _t_15_34;
      }


    template <typename number>
    void quartic<number>::M(const integration_task<number>* __task, const flattened_tensor<number>& __fields, double __N,
                           flattened_tensor<number>& __M)
      {
        DEFINE_INDEX_TOOLS
//         release resources 
        const auto& __pvector = __task->get_params().get_vector();
        __raw_params[0] = __pvector[0];

        const auto __Mp = __task->get_params().get_Mp();

//         parameters resource set to '__raw_params' 
//         coordinates resource set to '__fields' 
//         ENDIF !fast 

// BEGIN TEMPORARY POOL (sequence=16) 
const auto _t_16_0 = __Mp;
const auto _t_16_1 = -2.0;
const auto _t_16_2 = 1.0/(_t_16_0*_t_16_0);
const auto _t_16_3 = __fields[FLATTEN(1)];
const auto _t_16_4 = 2.0;
const auto _t_16_5 = (_t_16_3*_t_16_3);
const auto _t_16_6 = _t_16_2*_t_16_5;
const auto _t_16_7 = -6.0;
const auto _t_16_8 = _t_16_6+_t_16_7;
const auto _t_16_9 = -1.0;
const auto _t_16_10 = 1.0/(_t_16_8);
const auto _t_16_11 = __raw_params[0];
const auto _t_16_12 = __fields[FLATTEN(0)];
const auto _t_16_13 = 4.0;
const auto _t_16_14 = (_t_16_12*_t_16_12*_t_16_12*_t_16_12);
const auto _t_16_15 = -(1.0/2.0);
const auto _InternalHsq = _t_16_10*_t_16_11*_t_16_2*_t_16_14*_t_16_15;
const auto _t_16_16 = (1.0/2.0);
const auto _InternalEps = _t_16_2*_t_16_5*_t_16_16;
const auto _t_16_17 = 3.0;
const auto _t_16_18 = (_t_16_12*_t_16_12*_t_16_12);
const auto _t_16_19 = _InternalHsq;
const auto _t_16_20 = 1.0/_t_16_19;
const auto _t_16_21 = _t_16_11*_t_16_2*_t_16_18*_t_16_20*_t_16_3*_t_16_4;
const auto _t_16_22 = (_t_16_12*_t_16_12);
const auto _t_16_23 = _t_16_11*_t_16_22*_t_16_20*_t_16_17;
const auto _t_16_24 = _InternalEps;
const auto _t_16_25 = -3.0;
const auto _t_16_26 = _t_16_24+_t_16_25;
const auto _t_16_27 = _t_16_2*_t_16_5*(_t_16_26)*_t_16_9;
const auto _t_16_28 = _t_16_21+_t_16_23+_t_16_27;
        // END TEMPORARY POOL (sequence=16) 

        __M[FIELDS_FLATTEN(0,0)] = _t_16_28;
      }


    template <typename number>
    void quartic<number>::sorted_mass_spectrum(const integration_task<number>* __task, const flattened_tensor<number>& __fields,
                                              double __N, bool __norm, flattened_tensor<number>& __M, flattened_tensor<number>& __E)
      {
        // get raw, unsorted mass spectrum in __E
        this->mass_spectrum(__task, __fields, __N, __M, __E);

        // sort mass spectrum into order
        std::sort(__E.begin(), __E.end());

        // if normalized values requested, divide through by 1/H^2
        if(__norm)
          {
            DEFINE_INDEX_TOOLS

//             release resources 
            const auto& __pvector = __task->get_params().get_vector();
            __raw_params[0] = __pvector[0];

            const auto __Mp = __task->get_params().get_Mp();

// BEGIN TEMPORARY POOL (sequence=17) 
const auto _t_17_0 = __Mp;
const auto _t_17_1 = -2.0;
const auto _t_17_2 = 1.0/(_t_17_0*_t_17_0);
const auto _t_17_3 = __fields[FLATTEN(1)];
const auto _t_17_4 = 2.0;
const auto _t_17_5 = (_t_17_3*_t_17_3);
const auto _t_17_6 = _t_17_2*_t_17_5;
const auto _t_17_7 = -6.0;
const auto _t_17_8 = _t_17_6+_t_17_7;
const auto _t_17_9 = -1.0;
const auto _t_17_10 = 1.0/(_t_17_8);
const auto _t_17_11 = __raw_params[0];
const auto _t_17_12 = __fields[FLATTEN(0)];
const auto _t_17_13 = 4.0;
const auto _t_17_14 = (_t_17_12*_t_17_12*_t_17_12*_t_17_12);
const auto _t_17_15 = -(1.0/2.0);
const auto _InternalHsq = _t_17_10*_t_17_11*_t_17_2*_t_17_14*_t_17_15;
const auto _t_17_16 = _InternalHsq;
            // END TEMPORARY POOL (sequence=17) 

//             parameters resource set to '__raw_params' 
//             coordinates resource set to '__fields' 

            const auto __Hsq = _t_17_16;

            __E[0] = __E[0] / __Hsq;
          };
      }


    template <typename number>
    void quartic<number>::mass_spectrum(const integration_task<number>* __task, const flattened_tensor<number>& __fields,
                                       double __N, flattened_tensor<number>& __M, flattened_tensor<number>& __E)
      {
        DEFINE_INDEX_TOOLS

        // write mass matrix (in canonical format) into __M
        this->M(__task, __fields, __N, __M);

        // copy elements of the mass matrix into an Eigen matrix
        __mass_matrix(0,0) = __M[FIELDS_FLATTEN(0,0)];

        // extract eigenvalues from this matrix
        // In general Eigen would give us complex results, which we'd like to avoid. That can be done by
        // forcing Eigen to use a self-adjoint matrix, which has guaranteed real eigenvalues
        auto __evalues = __mass_matrix.template selfadjointView<Eigen::Upper>().eigenvalues();

        // copy eigenvalues into output matrix
        __E[FIELDS_FLATTEN(0)] = __evalues(0);
      }


    template <typename number>
    void quartic<number>::backend_process_backg(const background_task<number>* tk, backg_history<number>& solution, bool silent)
      {
        DEFINE_INDEX_TOOLS

        assert(tk != nullptr);

        const time_config_database time_db = tk->get_time_config_database();

        solution.clear();
        solution.reserve(time_db.size());

//        if(!silent) this->write_task_data(tk, std::cout, 9.9999999999999998e-13, 9.9999999999999998e-13, 1e-10, "runge_kutta_dopri5");

        // set up an observer which writes to this history vector
        quartic_background_observer<number> obs(solution, time_db);

        // set up a functor to evolve this system
        quartic_background_functor<number> system(tk->get_params());
        system.set_up_workspace();

        auto ics = tk->get_ics_vector();

        backg_state<number> x(quartic_pool::backg_state_size);
        x[FLATTEN(0)] = ics[0];
        x[FLATTEN(1)] = ics[1];

        using boost::numeric::odeint::integrate_times;

        auto stepper = boost::numeric::odeint::make_dense_output< boost::numeric::odeint::runge_kutta_dopri5< backg_state<number>, number, backg_state<number>, number, CPPTRANSPORT_ALGEBRA_NAME(backg_state<number>), CPPTRANSPORT_OPERATIONS_NAME(backg_state<number>) > >(1e-12, 1e-12);
        integrate_times(stepper, system, x, time_db.value_begin(), time_db.value_end(),
                        static_cast<number>(1e-10), obs);

        system.close_down_workspace();
      }


    namespace quartic_impl
      {

        template <typename number>
        class EpsilonUnityPredicate
          {
          public:
            EpsilonUnityPredicate(const parameters<number>& p)
              : __Mp(p.get_Mp())
              {
              }

            bool operator()(const std::pair< backg_state<number>, number >& __x)
              {
                DEFINE_INDEX_TOOLS
//                 release resources 
// BEGIN TEMPORARY POOL (sequence=18) 
const auto _t_18_0 = __x.first[FLATTEN(1)];
const auto _t_18_1 = 2.0;
const auto _t_18_2 = (_t_18_0*_t_18_0);
const auto _t_18_3 = __Mp;
const auto _t_18_4 = -2.0;
const auto _t_18_5 = 1.0/(_t_18_3*_t_18_3);
const auto _t_18_6 = (1.0/2.0);
const auto _InternalEps = _t_18_2*_t_18_5*_t_18_6;
const auto _t_18_7 = _InternalEps;
                // END TEMPORARY POOL (sequence=18) 

//                 coordinates resource set to '__x.first' 

                const auto __eps = _t_18_7;

                return (__eps >= 1.0 || __eps < 0.0);
              }

          private:

            //! cache Planck mass
            const number __Mp;

          };

      }   // namespace quartic_impl


    template <typename number>
    double quartic<number>::compute_end_of_inflation(const integration_task<number>* tk, double search_time)
      {
        DEFINE_INDEX_TOOLS

        assert(tk != nullptr);

        // set up a functor to evolve this system
        quartic_background_functor<number> system(tk->get_params());
        system.set_up_workspace();

        auto ics = tk->get_ics_vector();

        backg_state<number> x(quartic_pool::backg_state_size);
        x[FLATTEN(0)] = ics[0];
        x[FLATTEN(1)] = ics[1];

		    // find point where epsilon = 1
        using boost::numeric::odeint::make_adaptive_time_range;

        auto stepper = boost::numeric::odeint::make_dense_output< boost::numeric::odeint::runge_kutta_dopri5< backg_state<number>, number, backg_state<number>, number, CPPTRANSPORT_ALGEBRA_NAME(backg_state<number>), CPPTRANSPORT_OPERATIONS_NAME(backg_state<number>) > >(1e-12, 1e-12);
        auto range = make_adaptive_time_range(stepper, system, x, tk->get_N_initial(), tk->get_N_initial()+search_time, 1e-10);

        double Nend = 0.0;
        try
          {
            // returns the first iterator in 'range' for which the predicate EpsilonUnityPredicate() is satisfied
            auto iter = boost::find_if(range, quartic_impl::EpsilonUnityPredicate<number>(tk->get_params()));
            if(iter == boost::end(range)) throw end_of_inflation_not_found{};

            // may require explicit narrowing cast to double if working type is wider
            Nend = static_cast<double>(iter->second);
          }
        catch(integration_produced_nan& xe)
          {
            throw end_of_inflation_not_found{};
          }
        catch(Hsq_is_negative& xe)
          {
            throw end_of_inflation_not_found{};
          }

        system.close_down_workspace();

        return Nend;
      };


    namespace quartic_impl
      {

        template <typename number>
        class aHAggregatorPredicate
          {
          public:
            aHAggregatorPredicate(const integration_task<number>* tk, model<number>* m, std::vector<double>& N,
                                  flattened_tensor<number>& log_aH, flattened_tensor<number>& log_a2H2M,
                                  boost::optional<double>& lk)
              : params(tk->get_params()),
                task(tk),
                mdl(m),
                N_vector(N),
                log_aH_vector(log_aH),
                log_a2H2M_vector(log_a2H2M),
                largest_k(lk),
                flat_M(1*1),
                flat_E(1),
                N_horizon_crossing(tk->get_N_horizon_crossing()),
                astar_normalization(tk->get_astar_normalization()),
                __Mp(params.get_Mp()),
                __params(params.get_vector())
              {
              }

            number largest_evalue(const backg_state<number>& fields, number N)
              {
                this->mdl->mass_spectrum(this->task, fields, static_cast<double>(N), this->flat_M, this->flat_E);

                // step through eigenvalue vector, extracting largest absolute value
                number largest_eigenvalue = -std::numeric_limits<number>().max();

                if(std::abs(this->flat_E[0]) > largest_eigenvalue) { largest_eigenvalue = std::abs(this->flat_E[0]); }

                return largest_eigenvalue;
              }

            bool operator()(const std::pair< backg_state<number>, number >& __x)
              {
                DEFINE_INDEX_TOOLS
//                 release resources 
// BEGIN TEMPORARY POOL (sequence=19) 
const auto _t_19_0 = __x.first[FLATTEN(0)];
const auto _t_19_1 = 4.0;
const auto _t_19_2 = (_t_19_0*_t_19_0*_t_19_0*_t_19_0);
const auto _t_19_3 = __Mp;
const auto _t_19_4 = -2.0;
const auto _t_19_5 = 1.0/(_t_19_3*_t_19_3);
const auto _t_19_6 = __params[0];
const auto _t_19_7 = __x.first[FLATTEN(1)];
const auto _t_19_8 = 2.0;
const auto _t_19_9 = (_t_19_7*_t_19_7);
const auto _t_19_10 = _t_19_9*_t_19_5;
const auto _t_19_11 = -6.0;
const auto _t_19_12 = _t_19_10+_t_19_11;
const auto _t_19_13 = -1.0;
const auto _t_19_14 = 1.0/(_t_19_12);
const auto _t_19_15 = -(1.0/2.0);
const auto _InternalHsq = _t_19_2*_t_19_5*_t_19_6*_t_19_14*_t_19_15;
const auto _t_19_16 = _InternalHsq;
                // END TEMPORARY POOL (sequence=19) 

//                 parameters resource set to '__params' 
//                 coordinates resource set to '__x.first' 

                const auto __Hsq = _t_19_16;
                const auto __H   = std::sqrt(__Hsq);

                const auto __N   = __x.second - this->N_horizon_crossing + this->astar_normalization;

                this->N_vector.push_back(static_cast<double>(__x.second));
                this->log_aH_vector.push_back(__N + std::log(__H)); // = log(aH)
                this->log_a2H2M_vector.push_back(2.0*__N + 2.0*std::log(__H)
                                                 + std::log(this->largest_evalue(__x.first, __x.second))); // = log(a^2 H^2 * largest eigenvalue)

                // if a largest k-mode was provided,
                // are we now at a point where we have comfortably covered the horizon crossing time for it?
                if(!this->largest_k) return false;

                if(std::log(*largest_k) - __N - std::log(__H) < -0.5) return true;
                return false;
              }

          private:

            //! pointer to model object
            model<number>* mdl;

            //! point to task object
            const integration_task<number>* task;

            //! parameters for the model in use
            const parameters<number>& params;

            //! cache parameters vectors
            const flattened_tensor<number>& __params;

            //! cache Planck mass
            const number __Mp;

            //! output vector for times N
            std::vector<double>& N_vector;

            //! output vector for values log(aH)
            flattened_tensor<number>& log_aH_vector;

            //! output vector for field values
            flattened_tensor<number>& log_a2H2M_vector;

            //! working space for calculation of mass matrix
            flattened_tensor<number> flat_M;

            //! working space for calculation of mass eigenvalues
            flattened_tensor<number> flat_E;

            //! largest k-mode for which we are trying to find a horizon-exit time
            const boost::optional<double>& largest_k;

            //! time of horizon crossing
            const double N_horizon_crossing;

            //! normalization of ln(a) at horizon crossing time
            const double astar_normalization;

          };

        template <typename number>
        class HAggregatorPredicate
        {
        public:
          HAggregatorPredicate(const integration_task<number>* tk, model<number>* m, std::vector<double>& N,
                               flattened_tensor<number>& log_H, boost::optional<double>& lk)
            : params(tk->get_params()),
              task(tk),
              mdl(m),
              N_vector(N),
              log_H_vector(log_H),
              largest_k(lk),
              __Mp(params.get_Mp()),
              __params(params.get_vector())
          {
          }

          bool operator()(const std::pair< backg_state<number>, number >& __x)
          {
            DEFINE_INDEX_TOOLS
//             release resources 
// BEGIN TEMPORARY POOL (sequence=20) 
const auto _t_20_0 = __x.first[FLATTEN(0)];
const auto _t_20_1 = 4.0;
const auto _t_20_2 = (_t_20_0*_t_20_0*_t_20_0*_t_20_0);
const auto _t_20_3 = __Mp;
const auto _t_20_4 = -2.0;
const auto _t_20_5 = 1.0/(_t_20_3*_t_20_3);
const auto _t_20_6 = __params[0];
const auto _t_20_7 = __x.first[FLATTEN(1)];
const auto _t_20_8 = 2.0;
const auto _t_20_9 = (_t_20_7*_t_20_7);
const auto _t_20_10 = _t_20_9*_t_20_5;
const auto _t_20_11 = -6.0;
const auto _t_20_12 = _t_20_10+_t_20_11;
const auto _t_20_13 = -1.0;
const auto _t_20_14 = 1.0/(_t_20_12);
const auto _t_20_15 = -(1.0/2.0);
const auto _InternalHsq = _t_20_2*_t_20_5*_t_20_6*_t_20_14*_t_20_15;
const auto _t_20_16 = _InternalHsq;
            // END TEMPORARY POOL (sequence=20) 

//             parameters resource set to '__params' 
//             coordinates resource set to '__x.first' 

            const auto __Hsq = _t_20_16;
            const auto __H   = std::sqrt(__Hsq);

            this->N_vector.push_back(static_cast<double>(__x.second));
            this->log_H_vector.push_back(std::log(__H)); // = log(H)

            // if a largest k-mode was provided,
            // are we now at a point where we have comfortably covered the horizon crossing time for it?
            if(!this->largest_k) return false;

            // if(std::log(*largest_k) - __N - std::log(__H) < -0.5) return true;
            // return false;
          }

        private:
          //! pointer to model object
          model<number>* mdl;

          //! point to task object
          const integration_task<number>* task;

          //! parameters for the model in use
          const parameters<number>& params;

          //! cache parameters vectors
          const flattened_tensor<number>& __params;

          //! cache Planck mass
          const number __Mp;

          //! output vector for times N
          std::vector<double>& N_vector;

          //! output vector for values log(aH)
          flattened_tensor<number>& log_H_vector;

          //! largest k-mode for which we are trying to find a horizon-exit time
          const boost::optional<double>& largest_k;
        };

      }   // namespace quartic_impl


		template <typename number>
		void quartic<number>::compute_aH(const integration_task<number>* tk, std::vector<double>& N,
		                                flattened_tensor<number>& log_aH, flattened_tensor<number>& log_a2H2M,
		                                boost::optional<double> largest_k)
			{
        DEFINE_INDEX_TOOLS

				N.clear();
				log_aH.clear();
        log_a2H2M.clear();

				// set up a functor to evolve the system
				quartic_background_functor<number> system(tk->get_params());
        system.set_up_workspace();

				auto ics = tk->get_ics_vector();

				backg_state<number> x(quartic_pool::backg_state_size);
				x[FLATTEN(0)] = ics[0];
				x[FLATTEN(1)] = ics[1];

        double N_range = 0.0;
        bool found_end = false;
        try
          {
            N_range   = tk->get_N_end_of_inflation();
            found_end = true;
          }
        catch(end_of_inflation_not_found& xe)
          {
            // try to fall back on a sensible default
            N_range = tk->get_N_initial() + CPPTRANSPORT_DEFAULT_END_OF_INFLATION_SEARCH;
          }

        using boost::numeric::odeint::make_adaptive_time_range;

        auto stepper = boost::numeric::odeint::make_dense_output< boost::numeric::odeint::runge_kutta_dopri5< backg_state<number>, number, backg_state<number>, number, CPPTRANSPORT_ALGEBRA_NAME(backg_state<number>), CPPTRANSPORT_OPERATIONS_NAME(backg_state<number>) > >(1e-12, 1e-12);
        auto range = make_adaptive_time_range(stepper, system, x, tk->get_N_initial(), N_range, 1e-10);

        quartic_impl::aHAggregatorPredicate<number> aggregator(tk, this, N, log_aH, log_a2H2M, largest_k);

				// step through iterators, finding first point which is comfortably after time when largest_k has left
				// the horizon
				// aggregator writes N, log_aH and the field values into the output vectors at each iteration
				auto iter = boost::find_if(range, aggregator);

				// if we got to the end of the range, then we didn't cover all exit times up to largest_k
				// so something has gone wrong
				if(iter == boost::end(range) && largest_k)
					{
            throw failed_to_compute_horizon_exit(tk->get_N_initial(), N_range, found_end, log_aH.size(),
                                                 (N.size() > 0 ? N.back() : 0.0),
                                                 (log_aH.size() > 0 ? static_cast<double>(log_aH.back()) : 0.0),
                                                 largest_k.get());
					}

        system.close_down_workspace();
			}

    template <typename number>
    void quartic<number>::compute_H(const integration_task<number>* tk, std::vector<double>& N,
                                   flattened_tensor<number>& log_H, boost::optional<double> largest_k)
    {
      DEFINE_INDEX_TOOLS

      N.clear();
      log_H.clear();

      // set up a functor to evolve the system
      quartic_background_functor<number> system(tk->get_params());
      system.set_up_workspace();

      auto ics = tk->get_ics_vector();

      backg_state<number> x(quartic_pool::backg_state_size);
      x[FLATTEN(0)] = ics[0];
      x[FLATTEN(1)] = ics[1];

      double N_range = 0.0;
      bool found_end = false;
      try
      {
        N_range   = tk->get_N_end_of_inflation();
        found_end = true;
      }
      catch(end_of_inflation_not_found& xe)
      {
        // try to fall back on a sensible default
        N_range = tk->get_N_initial() + CPPTRANSPORT_DEFAULT_END_OF_INFLATION_SEARCH;
      }

      using boost::numeric::odeint::make_adaptive_time_range;

      auto stepper = boost::numeric::odeint::make_dense_output< boost::numeric::odeint::runge_kutta_dopri5< backg_state<number>, number, backg_state<number>, number, CPPTRANSPORT_ALGEBRA_NAME(backg_state<number>), CPPTRANSPORT_OPERATIONS_NAME(backg_state<number>) > >(1e-12, 1e-12);
      auto range = make_adaptive_time_range(stepper, system, x, tk->get_N_initial(), N_range, 1e-10);

      quartic_impl::HAggregatorPredicate<number> aggregator(tk, this, N, log_H, largest_k);

      // step through iterators, finding first point which is comfortably after time when largest_k has left
      // the horizon
      // aggregator writes N, log_aH and the field values into the output vectors at each iteration
      auto iter = boost::find_if(range, aggregator);

      // if we got to the end of the range, then we didn't cover all exit times up to largest_k
      // so something has gone wrong
      if(iter == boost::end(range) && largest_k)
      {
        throw failed_to_compute_horizon_exit(tk->get_N_initial(), N_range, found_end, log_H.size(),
                                             (N.size() > 0 ? N.back() : 0.0),
                                             (log_H.size() > 0 ? static_cast<double>(log_H.back()) : 0.0),
                                             largest_k.get());
      }

      system.close_down_workspace();
    }

    // IMPLEMENTATION - FUNCTOR FOR BACKGROUND INTEGRATION


    template <typename number>
    void quartic_background_functor<number>::operator()(const backg_state<number>& __x, backg_state<number>& __dxdt, number __t)
      {
        DEFINE_INDEX_TOOLS
//         release resources 

//         parameters resource set to '__raw_params' 
//         coordinates resource set to '__x' 

//         ENDIF !fast 

// BEGIN TEMPORARY POOL (sequence=21) 
const auto _t_21_0 = __x[FLATTEN(1)];
const auto _t_21_1 = 2.0;
const auto _t_21_2 = (_t_21_0*_t_21_0);
const auto _t_21_3 = __Mp;
const auto _t_21_4 = -2.0;
const auto _t_21_5 = 1.0/(_t_21_3*_t_21_3);
const auto _t_21_6 = _t_21_2*_t_21_5;
const auto _t_21_7 = -6.0;
const auto _t_21_8 = _t_21_6+_t_21_7;
const auto _t_21_9 = -1.0;
const auto _t_21_10 = 1.0/(_t_21_8);
const auto _t_21_11 = __raw_params[0];
const auto _t_21_12 = __x[FLATTEN(0)];
const auto _t_21_13 = 4.0;
const auto _t_21_14 = (_t_21_12*_t_21_12*_t_21_12*_t_21_12);
const auto _t_21_15 = -(1.0/2.0);
const auto _InternalHsq = _t_21_10*_t_21_11*_t_21_5*_t_21_14*_t_21_15;
const auto _t_21_16 = _InternalHsq;
const auto _t_21_17 = (1.0/2.0);
const auto _InternalEps = _t_21_2*_t_21_5*_t_21_17;
const auto _t_21_18 = _InternalEps;
const auto _t_21_19 = (1.0/4.0);
const auto _InternalV = _t_21_11*_t_21_14*_t_21_19;
const auto _t_21_20 = _InternalV;
const auto _t_21_21 = -3.0;
const auto _t_21_22 = _t_21_18+_t_21_21;
const auto _t_21_23 = _t_21_0*(_t_21_22);
const auto _t_21_24 = 1.0/_t_21_16;
const auto _t_21_25 = 3.0;
const auto _t_21_26 = (_t_21_12*_t_21_12*_t_21_12);
const auto _t_21_27 = _t_21_11*_t_21_24*_t_21_26*_t_21_9;
const auto _t_21_28 = _t_21_23+_t_21_27;
        // END TEMPORARY POOL (sequence=21) 

        const auto __Hsq = _t_21_16;
        const auto __eps = _t_21_18;
        const auto __V = _t_21_20;

        // check whether 0 < epsilon < 3
        if(__eps < 0.0) throw eps_is_negative(static_cast<double>(__t), static_cast<double>(__Hsq), static_cast<double>(__eps), static_cast<double>(__V), __x, quartic_pool::state_names);
        if(__eps > 3.0) throw eps_too_large(static_cast<double>(__t), static_cast<double>(__Hsq), static_cast<double>(__eps), static_cast<double>(__V), __x, quartic_pool::state_names);

        // check whether potential is +ve definite
        if(__V < 0.0) throw V_is_negative(static_cast<double>(__t), static_cast<double>(__Hsq), static_cast<double>(__eps), static_cast<double>(__V), __x, quartic_pool::state_names);

        // check whether Hsq is positive
        if(__Hsq < 0.0) throw Hsq_is_negative(static_cast<double>(__t), static_cast<double>(__Hsq), static_cast<double>(__eps), static_cast<double>(__V), __x, quartic_pool::state_names);

        // check for nan being produced
        if(std::isnan(__x[0])) throw integration_produced_nan(static_cast<double>(__t), static_cast<double>(__Hsq), static_cast<double>(__eps), static_cast<double>(__V), __x, quartic_pool::state_names);
        if(std::isnan(__x[1])) throw integration_produced_nan(static_cast<double>(__t), static_cast<double>(__Hsq), static_cast<double>(__eps), static_cast<double>(__V), __x, quartic_pool::state_names);

        __dxdt[FLATTEN(0)] = _t_21_0;
        __dxdt[FLATTEN(1)] = _t_21_28;
      }


    // IMPLEMENTATION - FUNCTOR FOR BACKGROUND OBSERVATION


    template <typename number>
    void quartic_background_observer<number>::operator()(const backg_state<number>& x, number t)
      {
        if(this->current_step != this->time_db.record_end() && this->current_step->is_stored())
          {
            this->history.push_back(x);
          }
        this->current_step++;
      }


  };   // namespace transport


#endif  // CPPTRANSPORT_QUARTIC_CORE_H

