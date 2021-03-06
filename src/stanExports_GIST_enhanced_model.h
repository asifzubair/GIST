// Generated by rstantools.  Do not edit by hand.

/*
    GIST is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    GIST is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with GIST.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef MODELS_HPP
#define MODELS_HPP
#define STAN__SERVICES__COMMAND_HPP
#include <rstan/rstaninc.hpp>
// Code generated by Stan version 2.21.0
#include <stan/model/model_header.hpp>
namespace model_GIST_enhanced_model_namespace {
using std::istream;
using std::string;
using std::stringstream;
using std::vector;
using stan::io::dump;
using stan::math::lgamma;
using stan::model::prob_grad;
using namespace stan::math;
static int current_statement_begin__;
stan::io::program_reader prog_reader__() {
    stan::io::program_reader reader;
    reader.add_event(0, 0, "start", "model_GIST_enhanced_model");
    reader.add_event(68, 66, "end", "model_GIST_enhanced_model");
    return reader;
}
#include <stan_meta_header.hpp>
class model_GIST_enhanced_model
  : public stan::model::model_base_crtp<model_GIST_enhanced_model> {
private:
        int numGenes;
        int numCellTypes;
        vector_d exprMixVec;
        double prior_phi;
        double prior_lambda;
        int prior_index;
        matrix_d sigMat;
        vector_d alpha;
public:
    model_GIST_enhanced_model(stan::io::var_context& context__,
        std::ostream* pstream__ = 0)
        : model_base_crtp(0) {
        ctor_body(context__, 0, pstream__);
    }
    model_GIST_enhanced_model(stan::io::var_context& context__,
        unsigned int random_seed__,
        std::ostream* pstream__ = 0)
        : model_base_crtp(0) {
        ctor_body(context__, random_seed__, pstream__);
    }
    void ctor_body(stan::io::var_context& context__,
                   unsigned int random_seed__,
                   std::ostream* pstream__) {
        typedef double local_scalar_t__;
        boost::ecuyer1988 base_rng__ =
          stan::services::util::create_rng(random_seed__, 0);
        (void) base_rng__;  // suppress unused var warning
        current_statement_begin__ = -1;
        static const char* function__ = "model_GIST_enhanced_model_namespace::model_GIST_enhanced_model";
        (void) function__;  // dummy to suppress unused var warning
        size_t pos__;
        (void) pos__;  // dummy to suppress unused var warning
        std::vector<int> vals_i__;
        std::vector<double> vals_r__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning
        try {
            // initialize data block variables from context__
            current_statement_begin__ = 11;
            context__.validate_dims("data initialization", "numGenes", "int", context__.to_vec());
            numGenes = int(0);
            vals_i__ = context__.vals_i("numGenes");
            pos__ = 0;
            numGenes = vals_i__[pos__++];
            check_greater_or_equal(function__, "numGenes", numGenes, 0);
            current_statement_begin__ = 12;
            context__.validate_dims("data initialization", "numCellTypes", "int", context__.to_vec());
            numCellTypes = int(0);
            vals_i__ = context__.vals_i("numCellTypes");
            pos__ = 0;
            numCellTypes = vals_i__[pos__++];
            check_greater_or_equal(function__, "numCellTypes", numCellTypes, 0);
            current_statement_begin__ = 13;
            validate_non_negative_index("exprMixVec", "numGenes", numGenes);
            context__.validate_dims("data initialization", "exprMixVec", "vector_d", context__.to_vec(numGenes));
            exprMixVec = Eigen::Matrix<double, Eigen::Dynamic, 1>(numGenes);
            vals_r__ = context__.vals_r("exprMixVec");
            pos__ = 0;
            size_t exprMixVec_j_1_max__ = numGenes;
            for (size_t j_1__ = 0; j_1__ < exprMixVec_j_1_max__; ++j_1__) {
                exprMixVec(j_1__) = vals_r__[pos__++];
            }
            current_statement_begin__ = 15;
            context__.validate_dims("data initialization", "prior_phi", "double", context__.to_vec());
            prior_phi = double(0);
            vals_r__ = context__.vals_r("prior_phi");
            pos__ = 0;
            prior_phi = vals_r__[pos__++];
            check_greater_or_equal(function__, "prior_phi", prior_phi, 0);
            current_statement_begin__ = 16;
            context__.validate_dims("data initialization", "prior_lambda", "double", context__.to_vec());
            prior_lambda = double(0);
            vals_r__ = context__.vals_r("prior_lambda");
            pos__ = 0;
            prior_lambda = vals_r__[pos__++];
            check_greater_or_equal(function__, "prior_lambda", prior_lambda, 1);
            current_statement_begin__ = 18;
            context__.validate_dims("data initialization", "prior_index", "int", context__.to_vec());
            prior_index = int(0);
            vals_i__ = context__.vals_i("prior_index");
            pos__ = 0;
            prior_index = vals_i__[pos__++];
            check_greater_or_equal(function__, "prior_index", prior_index, 0);
            current_statement_begin__ = 21;
            validate_non_negative_index("sigMat", "numGenes", numGenes);
            validate_non_negative_index("sigMat", "numCellTypes", numCellTypes);
            context__.validate_dims("data initialization", "sigMat", "matrix_d", context__.to_vec(numGenes,numCellTypes));
            sigMat = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>(numGenes, numCellTypes);
            vals_r__ = context__.vals_r("sigMat");
            pos__ = 0;
            size_t sigMat_j_2_max__ = numCellTypes;
            size_t sigMat_j_1_max__ = numGenes;
            for (size_t j_2__ = 0; j_2__ < sigMat_j_2_max__; ++j_2__) {
                for (size_t j_1__ = 0; j_1__ < sigMat_j_1_max__; ++j_1__) {
                    sigMat(j_1__, j_2__) = vals_r__[pos__++];
                }
            }
            // initialize transformed data variables
            current_statement_begin__ = 26;
            validate_non_negative_index("alpha", "numCellTypes", numCellTypes);
            alpha = Eigen::Matrix<double, Eigen::Dynamic, 1>(numCellTypes);
            stan::math::fill(alpha, DUMMY_VAR__);
            // execute transformed data statements
            current_statement_begin__ = 27;
            for (int i = 1; i <= numCellTypes; ++i) {
                current_statement_begin__ = 31;
                stan::model::assign(alpha, 
                            stan::model::cons_list(stan::model::index_uni(i), stan::model::nil_index_list()), 
                            1, 
                            "assigning variable alpha");
            }
            // validate transformed data
            current_statement_begin__ = 26;
            check_greater_or_equal(function__, "alpha", alpha, 0);
            // validate, set parameter ranges
            num_params_r__ = 0U;
            param_ranges_i__.clear();
            current_statement_begin__ = 39;
            validate_non_negative_index("estimatedProportionsVecSimp", "numCellTypes", numCellTypes);
            num_params_r__ += (numCellTypes - 1);
            current_statement_begin__ = 41;
            num_params_r__ += 1;
            current_statement_begin__ = 44;
            num_params_r__ += 1;
            current_statement_begin__ = 45;
            num_params_r__ += 1;
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }
    }
    ~model_GIST_enhanced_model() { }
    void transform_inits(const stan::io::var_context& context__,
                         std::vector<int>& params_i__,
                         std::vector<double>& params_r__,
                         std::ostream* pstream__) const {
        typedef double local_scalar_t__;
        stan::io::writer<double> writer__(params_r__, params_i__);
        size_t pos__;
        (void) pos__; // dummy call to supress warning
        std::vector<double> vals_r__;
        std::vector<int> vals_i__;
        current_statement_begin__ = 39;
        if (!(context__.contains_r("estimatedProportionsVecSimp")))
            stan::lang::rethrow_located(std::runtime_error(std::string("Variable estimatedProportionsVecSimp missing")), current_statement_begin__, prog_reader__());
        vals_r__ = context__.vals_r("estimatedProportionsVecSimp");
        pos__ = 0U;
        validate_non_negative_index("estimatedProportionsVecSimp", "numCellTypes", numCellTypes);
        context__.validate_dims("parameter initialization", "estimatedProportionsVecSimp", "vector_d", context__.to_vec(numCellTypes));
        Eigen::Matrix<double, Eigen::Dynamic, 1> estimatedProportionsVecSimp(numCellTypes);
        size_t estimatedProportionsVecSimp_j_1_max__ = numCellTypes;
        for (size_t j_1__ = 0; j_1__ < estimatedProportionsVecSimp_j_1_max__; ++j_1__) {
            estimatedProportionsVecSimp(j_1__) = vals_r__[pos__++];
        }
        try {
            writer__.simplex_unconstrain(estimatedProportionsVecSimp);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(std::runtime_error(std::string("Error transforming variable estimatedProportionsVecSimp: ") + e.what()), current_statement_begin__, prog_reader__());
        }
        current_statement_begin__ = 41;
        if (!(context__.contains_r("nu")))
            stan::lang::rethrow_located(std::runtime_error(std::string("Variable nu missing")), current_statement_begin__, prog_reader__());
        vals_r__ = context__.vals_r("nu");
        pos__ = 0U;
        context__.validate_dims("parameter initialization", "nu", "double", context__.to_vec());
        double nu(0);
        nu = vals_r__[pos__++];
        try {
            writer__.scalar_lb_unconstrain(1, nu);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(std::runtime_error(std::string("Error transforming variable nu: ") + e.what()), current_statement_begin__, prog_reader__());
        }
        current_statement_begin__ = 44;
        if (!(context__.contains_r("sigma")))
            stan::lang::rethrow_located(std::runtime_error(std::string("Variable sigma missing")), current_statement_begin__, prog_reader__());
        vals_r__ = context__.vals_r("sigma");
        pos__ = 0U;
        context__.validate_dims("parameter initialization", "sigma", "double", context__.to_vec());
        double sigma(0);
        sigma = vals_r__[pos__++];
        try {
            writer__.scalar_lb_unconstrain(0, sigma);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(std::runtime_error(std::string("Error transforming variable sigma: ") + e.what()), current_statement_begin__, prog_reader__());
        }
        current_statement_begin__ = 45;
        if (!(context__.contains_r("beta0")))
            stan::lang::rethrow_located(std::runtime_error(std::string("Variable beta0 missing")), current_statement_begin__, prog_reader__());
        vals_r__ = context__.vals_r("beta0");
        pos__ = 0U;
        context__.validate_dims("parameter initialization", "beta0", "double", context__.to_vec());
        double beta0(0);
        beta0 = vals_r__[pos__++];
        try {
            writer__.scalar_unconstrain(beta0);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(std::runtime_error(std::string("Error transforming variable beta0: ") + e.what()), current_statement_begin__, prog_reader__());
        }
        params_r__ = writer__.data_r();
        params_i__ = writer__.data_i();
    }
    void transform_inits(const stan::io::var_context& context,
                         Eigen::Matrix<double, Eigen::Dynamic, 1>& params_r,
                         std::ostream* pstream__) const {
      std::vector<double> params_r_vec;
      std::vector<int> params_i_vec;
      transform_inits(context, params_i_vec, params_r_vec, pstream__);
      params_r.resize(params_r_vec.size());
      for (int i = 0; i < params_r.size(); ++i)
        params_r(i) = params_r_vec[i];
    }
    template <bool propto__, bool jacobian__, typename T__>
    T__ log_prob(std::vector<T__>& params_r__,
                 std::vector<int>& params_i__,
                 std::ostream* pstream__ = 0) const {
        typedef T__ local_scalar_t__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // dummy to suppress unused var warning
        T__ lp__(0.0);
        stan::math::accumulator<T__> lp_accum__;
        try {
            stan::io::reader<local_scalar_t__> in__(params_r__, params_i__);
            // model parameters
            current_statement_begin__ = 39;
            Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> estimatedProportionsVecSimp;
            (void) estimatedProportionsVecSimp;  // dummy to suppress unused var warning
            if (jacobian__)
                estimatedProportionsVecSimp = in__.simplex_constrain(numCellTypes, lp__);
            else
                estimatedProportionsVecSimp = in__.simplex_constrain(numCellTypes);
            current_statement_begin__ = 41;
            local_scalar_t__ nu;
            (void) nu;  // dummy to suppress unused var warning
            if (jacobian__)
                nu = in__.scalar_lb_constrain(1, lp__);
            else
                nu = in__.scalar_lb_constrain(1);
            current_statement_begin__ = 44;
            local_scalar_t__ sigma;
            (void) sigma;  // dummy to suppress unused var warning
            if (jacobian__)
                sigma = in__.scalar_lb_constrain(0, lp__);
            else
                sigma = in__.scalar_lb_constrain(0);
            current_statement_begin__ = 45;
            local_scalar_t__ beta0;
            (void) beta0;  // dummy to suppress unused var warning
            if (jacobian__)
                beta0 = in__.scalar_constrain(lp__);
            else
                beta0 = in__.scalar_constrain();
            // model body
            {
            current_statement_begin__ = 52;
            local_scalar_t__ b_alpha(DUMMY_VAR__);
            (void) b_alpha;  // dummy to suppress unused var warning
            stan::math::initialize(b_alpha, DUMMY_VAR__);
            stan::math::fill(b_alpha, DUMMY_VAR__);
            stan::math::assign(b_alpha,(prior_lambda * prior_phi));
            current_statement_begin__ = 53;
            local_scalar_t__ b_beta(DUMMY_VAR__);
            (void) b_beta;  // dummy to suppress unused var warning
            stan::math::initialize(b_beta, DUMMY_VAR__);
            stan::math::fill(b_beta, DUMMY_VAR__);
            stan::math::assign(b_beta,(prior_lambda * (1 - prior_phi)));
            current_statement_begin__ = 56;
            lp_accum__.add(dirichlet_log<propto__>(estimatedProportionsVecSimp, alpha));
            current_statement_begin__ = 57;
            lp_accum__.add(beta_log<propto__>(get_base1(estimatedProportionsVecSimp, prior_index, "estimatedProportionsVecSimp", 1), b_alpha, b_beta));
            current_statement_begin__ = 62;
            lp_accum__.add(gamma_log<propto__>(nu, 2, 0.1));
            current_statement_begin__ = 64;
            lp_accum__.add(student_t_log<propto__>(exprMixVec, nu, add(beta0, multiply(sigMat, estimatedProportionsVecSimp)), sigma));
            }
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }
        lp_accum__.add(lp__);
        return lp_accum__.sum();
    } // log_prob()
    template <bool propto, bool jacobian, typename T_>
    T_ log_prob(Eigen::Matrix<T_,Eigen::Dynamic,1>& params_r,
               std::ostream* pstream = 0) const {
      std::vector<T_> vec_params_r;
      vec_params_r.reserve(params_r.size());
      for (int i = 0; i < params_r.size(); ++i)
        vec_params_r.push_back(params_r(i));
      std::vector<int> vec_params_i;
      return log_prob<propto,jacobian,T_>(vec_params_r, vec_params_i, pstream);
    }
    void get_param_names(std::vector<std::string>& names__) const {
        names__.resize(0);
        names__.push_back("estimatedProportionsVecSimp");
        names__.push_back("nu");
        names__.push_back("sigma");
        names__.push_back("beta0");
    }
    void get_dims(std::vector<std::vector<size_t> >& dimss__) const {
        dimss__.resize(0);
        std::vector<size_t> dims__;
        dims__.resize(0);
        dims__.push_back(numCellTypes);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dimss__.push_back(dims__);
    }
    template <typename RNG>
    void write_array(RNG& base_rng__,
                     std::vector<double>& params_r__,
                     std::vector<int>& params_i__,
                     std::vector<double>& vars__,
                     bool include_tparams__ = true,
                     bool include_gqs__ = true,
                     std::ostream* pstream__ = 0) const {
        typedef double local_scalar_t__;
        vars__.resize(0);
        stan::io::reader<local_scalar_t__> in__(params_r__, params_i__);
        static const char* function__ = "model_GIST_enhanced_model_namespace::write_array";
        (void) function__;  // dummy to suppress unused var warning
        // read-transform, write parameters
        Eigen::Matrix<double, Eigen::Dynamic, 1> estimatedProportionsVecSimp = in__.simplex_constrain(numCellTypes);
        size_t estimatedProportionsVecSimp_j_1_max__ = numCellTypes;
        for (size_t j_1__ = 0; j_1__ < estimatedProportionsVecSimp_j_1_max__; ++j_1__) {
            vars__.push_back(estimatedProportionsVecSimp(j_1__));
        }
        double nu = in__.scalar_lb_constrain(1);
        vars__.push_back(nu);
        double sigma = in__.scalar_lb_constrain(0);
        vars__.push_back(sigma);
        double beta0 = in__.scalar_constrain();
        vars__.push_back(beta0);
        double lp__ = 0.0;
        (void) lp__;  // dummy to suppress unused var warning
        stan::math::accumulator<double> lp_accum__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning
        if (!include_tparams__ && !include_gqs__) return;
        try {
            if (!include_gqs__ && !include_tparams__) return;
            if (!include_gqs__) return;
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }
    }
    template <typename RNG>
    void write_array(RNG& base_rng,
                     Eigen::Matrix<double,Eigen::Dynamic,1>& params_r,
                     Eigen::Matrix<double,Eigen::Dynamic,1>& vars,
                     bool include_tparams = true,
                     bool include_gqs = true,
                     std::ostream* pstream = 0) const {
      std::vector<double> params_r_vec(params_r.size());
      for (int i = 0; i < params_r.size(); ++i)
        params_r_vec[i] = params_r(i);
      std::vector<double> vars_vec;
      std::vector<int> params_i_vec;
      write_array(base_rng, params_r_vec, params_i_vec, vars_vec, include_tparams, include_gqs, pstream);
      vars.resize(vars_vec.size());
      for (int i = 0; i < vars.size(); ++i)
        vars(i) = vars_vec[i];
    }
    std::string model_name() const {
        return "model_GIST_enhanced_model";
    }
    void constrained_param_names(std::vector<std::string>& param_names__,
                                 bool include_tparams__ = true,
                                 bool include_gqs__ = true) const {
        std::stringstream param_name_stream__;
        size_t estimatedProportionsVecSimp_j_1_max__ = numCellTypes;
        for (size_t j_1__ = 0; j_1__ < estimatedProportionsVecSimp_j_1_max__; ++j_1__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "estimatedProportionsVecSimp" << '.' << j_1__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        param_name_stream__.str(std::string());
        param_name_stream__ << "nu";
        param_names__.push_back(param_name_stream__.str());
        param_name_stream__.str(std::string());
        param_name_stream__ << "sigma";
        param_names__.push_back(param_name_stream__.str());
        param_name_stream__.str(std::string());
        param_name_stream__ << "beta0";
        param_names__.push_back(param_name_stream__.str());
        if (!include_gqs__ && !include_tparams__) return;
        if (include_tparams__) {
        }
        if (!include_gqs__) return;
    }
    void unconstrained_param_names(std::vector<std::string>& param_names__,
                                   bool include_tparams__ = true,
                                   bool include_gqs__ = true) const {
        std::stringstream param_name_stream__;
        size_t estimatedProportionsVecSimp_j_1_max__ = (numCellTypes - 1);
        for (size_t j_1__ = 0; j_1__ < estimatedProportionsVecSimp_j_1_max__; ++j_1__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "estimatedProportionsVecSimp" << '.' << j_1__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        param_name_stream__.str(std::string());
        param_name_stream__ << "nu";
        param_names__.push_back(param_name_stream__.str());
        param_name_stream__.str(std::string());
        param_name_stream__ << "sigma";
        param_names__.push_back(param_name_stream__.str());
        param_name_stream__.str(std::string());
        param_name_stream__ << "beta0";
        param_names__.push_back(param_name_stream__.str());
        if (!include_gqs__ && !include_tparams__) return;
        if (include_tparams__) {
        }
        if (!include_gqs__) return;
    }
}; // model
}  // namespace
typedef model_GIST_enhanced_model_namespace::model_GIST_enhanced_model stan_model;
#ifndef USING_R
stan::model::model_base& new_model(
        stan::io::var_context& data_context,
        unsigned int seed,
        std::ostream* msg_stream) {
  stan_model* m = new stan_model(data_context, seed, msg_stream);
  return *m;
}
#endif
#endif
