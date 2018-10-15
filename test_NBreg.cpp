// We can now use the BH package
// [[Rcpp::depends(BH)]]

#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <Rmath.h>
#include <cmath>
#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/split.hpp>

// [[Rcpp::depends("RcppArmadillo")]]

template<typename T>
void printR_obj(const T& obj){
	Rcpp::Rcout << obj << std::endl;	
}

// const double LOG2PI = log(2*arma::datum::pi);

// --------------------
// Intermediate Functions
// --------------------

// [[Rcpp::export]]
double Rcpp_norm(const arma::vec& a){
	return arma::norm(a);
}

// [[Rcpp::export]]
double Rcpp_logSumExp(const arma::vec& log_x){
	if( log_x.n_elem == 1 ){
		return log_x.at(0);
	} else {
		double max_val = max(log_x);
		arma::vec log_x_2 = log_x - max_val;
		return log(arma::sum(arma::exp(log_x_2))) + max_val;
	}
}

// [[Rcpp::export]]
double Rcpp_min_diff(const arma::vec& x){
	double min_diff = 100.0;
	arma::uword i,j;
	arma::uword len_x = x.n_elem;
	
	if( len_x == 1 ){
		min_diff = 0.0;
	} else {
		double abs_diff;
		for(i = 0; i < len_x - 1; i++){
		for(j = i + 1; j < len_x; j++){
			abs_diff = std::abs(x.at(i) - x.at(j));
			if(abs_diff < min_diff) min_diff = abs_diff;
		}
		}
	}
	
	return min_diff;
}

// [[Rcpp::export]]
double Rcpp_log_binom_coef(const double& n, const double& x){
	double out = std::lgamma(n + 1.0);
	out -= std::lgamma(x + 1.0) + std::lgamma(n - x + 1.0);
	return out;
}

// --------------------
// Negative Binomial related functions
// --------------------
// printR_obj();

// ----------------
// Log Likelihood

// [[Rcpp::export]]
double Rcpp_log_NB(const arma::uword& x,
	const double& mu,const double& phi){
	
	double loglike = 0.0;
	double n = 1.0 / phi;
	// lg = std::lgamma(x + 1.0)
	
	if( x > 0 ){
		loglike += std::lgamma(x + n) -
			std::lgamma(n) -
			std::lgamma(x + 1.0);
		loglike += x * std::log(mu);
	}
	
	loglike += n * std::log(n) -
		(n + x) * std::log(n + mu);
	
	return loglike;
}

// [[Rcpp::export]]
double Rcpp_log_NB_2(const arma::uword& x,
	const double& mu,const double& phi,
	const double& lgx1){
	
	double loglike = 0.0;
	double n = 1.0 / phi;
	// lgx1 = std::lgamma(x + 1.0)
	
	if( x > 0 ){
		loglike += std::lgamma(x + n) -
			std::lgamma(n) - lgx1;
		loglike += x * std::log(mu);
	}
	
	loglike += n * std::log(n) -
		(n + x) * std::log(n + mu);
	
	return loglike;
}

// ----------------
// Gradient

// [[Rcpp::export]]
arma::vec Rcpp_grad_log_NB(const arma::uword& x,
	const double& mu,const double& phi){

	arma::vec grad = arma::zeros<arma::vec>(2);
	double n = 1.0 / phi;
	// part mu
	grad.at(0) = x / mu - (n + x) / (n + mu);
	
	// part n
	if(x > 0){
		grad.at(1) += R::digamma(x + n) - R::digamma(n);
	}
	grad.at(1) += 1.0 + std::log(n) -
		( (n + x) / (n + mu) + 
			std::log(n + mu) );
	grad.at(1) *= -1.0 * n * n;
	return grad;
}


// ----------------
// NB Regression (with offset)
// ----------------

// ----------------
// Log Likelihood

// [[Rcpp::export]]
double Rcpp_NB_reg_LL(const arma::vec& Y,
	const arma::mat& X,const arma::vec& O,
	const arma::vec& PARAMS){
	
	arma::uword ii;
	double mu;
	double LL = 0.0;
	arma::uword pp = X.n_cols;
	arma::vec BETA = PARAMS.subvec(0,pp-1);
	double PHI = std::exp(PARAMS.at(pp));
	
	for(ii = 0; ii < Y.n_elem; ii++){
		mu = std::exp( arma::dot(X.row(ii).t(),BETA) + O.at(ii) );
		LL += Rcpp_log_NB(Y.at(ii),mu,PHI);
	}
	
	return LL;
}

// ----------------
// Gradient

// [[Rcpp::export]]
arma::vec Rcpp_NB_reg_GRAD(const arma::vec& Y,
	const arma::mat& X,const arma::vec& O,
	const arma::vec& PARAMS){
	
	arma::uword ii;
	arma::uword pp = X.n_cols;
	arma::vec BETA = PARAMS.subvec(0,pp-1);
	double PHI = std::exp(PARAMS.at(pp));
	double mu;
	double nn = 1.0 / PHI;
	arma::vec GRAD = arma::zeros<arma::vec>(pp + 1);
	double digam_nn = R::digamma(nn);
	double log_nn = std::log(nn);
	
	for(ii = 0; ii < Y.n_elem; ii++){
		mu = std::exp( arma::dot(X.row(ii).t(),BETA) + O.at(ii) );
		// Part BETA
		GRAD.subvec(0,pp - 1) += (Y.at(ii) / mu - (nn + Y.at(ii)) / (nn + mu)) * mu * X.row(ii).t();
		// Part PHI*
		GRAD.at(pp) += (R::digamma(Y.at(ii) + nn) - digam_nn + 
			1.0 + log_nn - (nn + Y.at(ii)) / (nn + mu) - 
			std::log(nn + mu)) * (-1.0 * nn);
	}

	return GRAD;
}

// ----------------
// Log Likelihood

// [[Rcpp::export]]
double Rcpp_NB_reg_LL_logscale_PHI(const arma::vec& Y,
	const arma::mat& X,const arma::vec& O,
	const arma::vec& PARAMS){
	
	arma::uword ii;
	double mu;
	double LL = 0.0;
	arma::uword pp = X.n_cols;
	arma::vec BETA = PARAMS.subvec(0,pp-1);
	double PHI = PARAMS.at(pp);
	
	for(ii = 0; ii < Y.n_elem; ii++){
		mu = std::exp( arma::dot(X.row(ii).t(),BETA) + O.at(ii) );
		LL += Rcpp_log_NB(Y.at(ii),mu,PHI);
	}
	
	return LL;
}

// ----------------
// Gradient

// [[Rcpp::export]]
arma::vec Rcpp_NB_reg_GRAD_logscale_PHI(const arma::vec& Y,
	const arma::mat& X,const arma::vec& O,
	const arma::vec& PARAMS){
	
	arma::uword ii;
	arma::uword pp = X.n_cols;
	arma::vec BETA = PARAMS.subvec(0,pp-1);
	double PHI = PARAMS.at(pp);
	double mu;
	double nn = 1.0 / PHI;
	arma::vec GRAD = arma::zeros<arma::vec>(pp + 1);
	double digam_nn = R::digamma(nn);
	double log_nn = std::log(nn);
	
	for(ii = 0; ii < Y.n_elem; ii++){
		mu = std::exp( arma::dot(X.row(ii).t(),BETA) + O.at(ii) );
		// Part BETA
		GRAD.subvec(0,pp - 1) += (Y.at(ii) / mu - (nn + Y.at(ii)) / (nn + mu)) * mu * X.row(ii).t();
		// Part PHI*
		GRAD.at(pp) += (R::digamma(Y.at(ii) + nn) - digam_nn + 
			1.0 + log_nn - (nn + Y.at(ii)) / (nn + mu) - 
			std::log(nn + mu)) * (-1.0 * nn);
	}

	return GRAD;
}



// ----------------
// Hessian

// [[Rcpp::export]]
arma::mat Rcpp_NB_reg_HESS(const arma::vec& Y,
	const arma::mat& X,const arma::vec& O,
	const arma::vec& PARAMS){
	
	arma::uword ii;
	arma::uword pp = X.n_cols;
	arma::vec BETA = PARAMS.subvec(0,pp-1);
	double PHI = std::exp(PARAMS.at(pp));
	double mu;
	double nn = 1.0 / PHI;
	arma::mat HESS = arma::zeros<arma::mat>(pp + 1,pp + 1);
	arma::vec hess_beta_phi = arma::zeros<arma::vec>(pp);
	double trigam_nn = R::trigamma(nn);
	double digam_nn = R::digamma(nn);
	double log_nn = std::log(nn);
	
	for(ii = 0; ii < Y.n_elem; ii++){
		mu = std::exp( arma::dot(X.row(ii).t(),BETA) + O.at(ii) );
		// Part2 BETA
		HESS.submat(0,0,pp - 1,pp - 1) += 
			-1.0 * (nn + Y.at(ii)) / std::pow(nn + mu,2.0) * 
			nn * mu * X.row(ii).t() * X.row(ii);
		
		// Part2 BETA*PHI*
		hess_beta_phi = -1.0 * (Y.at(ii) - mu) / std::pow(nn + mu,2.0) * nn * mu * X.row(ii).t();
		HESS(arma::span(0,pp - 1),pp) += hess_beta_phi;
		HESS(pp,arma::span(0,pp - 1)) += hess_beta_phi.t();
		
		// Part2 PHI*
		HESS.at(pp,pp) += nn * ( nn * (R::trigamma(Y.at(ii) + nn) - trigam_nn + PHI +
			(Y.at(ii) - mu)/std::pow(nn + mu,2.0) - 1.0/(nn + mu) ) + 
			R::digamma(Y.at(ii) + nn) - digam_nn + 1.0 + 
			log_nn - (nn + Y.at(ii)) / (nn + mu) - std::log(nn + mu)
			);
	}

	return HESS;
}

// ----------------
// Optimize

// [[Rcpp::export]]
Rcpp::List Rcpp_NB_reg(const arma::vec& Y,
	const arma::mat& X,const arma::vec& O,
	const arma::vec& params0,
	const arma::uword& max_iter = 4e3,
	const double& eps = 1e-7){
	
	// arma::uword N = X.n_rows;
	arma::uword pp = X.n_cols;
	arma::uword jj,uu;
	arma::uword converge = 0;
	arma::uword iter = 0;
	double curr_LL = 0.0;
	double old_LL,new_LL,norm_GRAD,norm_FS;
	double new_LL_1,new_LL_2;
	
	// Initialize params, remember to unconstrain PHI
	arma::vec old_PARAMS = params0;
	arma::vec new_PARAMS = arma::zeros<arma::vec>(pp + 1);
	arma::vec curr_PARAMS = new_PARAMS;
	arma::vec curr_GRAD = new_PARAMS;
	arma::vec up_GRAD = new_PARAMS;
	arma::vec curr_FS = new_PARAMS;
	arma::vec up_FS = new_PARAMS;
	arma::mat curr_HESS = arma::zeros<arma::mat>(pp + 1,pp + 1);
	arma::vec new_params1 = arma::zeros<arma::vec>(pp + 1);
	arma::vec new_params2 = new_params1;
	
	while(iter < max_iter){
		old_LL = Rcpp_NB_reg_LL(Y,X,O,old_PARAMS);
		curr_GRAD = Rcpp_NB_reg_GRAD(Y,X,O,old_PARAMS);
		norm_GRAD = Rcpp_norm(curr_GRAD);
		up_GRAD = curr_GRAD / std::max(1.0,norm_GRAD);
		curr_HESS = Rcpp_NB_reg_HESS(Y,X,O,old_PARAMS);
		curr_FS = arma::inv_sympd(-1.0 * curr_HESS) * curr_GRAD;
		norm_FS = Rcpp_norm(curr_FS);
		up_FS = curr_FS / std::max(1.0,norm_FS);
		uu = 0;
		
		for(jj = 0; jj <= 20; jj++){
			new_params1 = old_PARAMS + up_GRAD / std::pow(4,jj);
			new_LL_1 = Rcpp_NB_reg_LL(Y,X,O,new_params1);
			new_params2 = old_PARAMS + up_FS / std::pow(4,jj);
			new_LL_2 = Rcpp_NB_reg_LL(Y,X,O,new_params2);
			new_LL = std::max(new_LL_1,new_LL_2);
			if( new_LL > old_LL ){
				if( new_LL == new_LL_1 ){
					old_PARAMS = new_params1;
				} else {
					old_PARAMS = new_params2;
				}
				old_LL = new_LL;
				uu = 1;
				break;
			}
		}
		
		if( uu == 0 ){
			break;
		}
		
		if(iter > 0){
			if( std::abs(curr_LL - old_LL) < eps && Rcpp_norm(curr_PARAMS - old_PARAMS) < eps ){
				curr_GRAD = Rcpp_NB_reg_GRAD(Y,X,O,old_PARAMS);
				norm_GRAD = Rcpp_norm(curr_GRAD);
				curr_HESS = Rcpp_NB_reg_HESS(Y,X,O,old_PARAMS);
				norm_FS = Rcpp_norm(arma::inv_sympd(-1.0 * curr_HESS) * curr_GRAD);
				if( norm_GRAD < eps && norm_FS < eps ){
					converge = 1;
					break;
				}
			}
		}
		
		curr_PARAMS = old_PARAMS;
		curr_LL = old_LL;
		iter++;
	}
	
	// return R_NilValue;
	return Rcpp::List::create(
		Rcpp::Named("converge",converge),
		Rcpp::Named("LL",old_LL),
		Rcpp::Named("iter",iter),
		Rcpp::Named("norm_FS",norm_FS),
		Rcpp::Named("norm_GRAD",norm_GRAD),
		Rcpp::Named("PARAMS",Rcpp::NumericVector(old_PARAMS.begin(),old_PARAMS.end()))
		);
}

// [[Rcpp::export]]
void Rcpp_NB_reg2(const arma::vec& Y,
	const arma::mat& X,const arma::vec& O,
	const arma::vec& params0,
	arma::vec& PARAMS,
	const arma::uword& max_iter = 4e3,
	const double& eps = 1e-7){
	
	// arma::uword N = X.n_rows;
	arma::uword pp = X.n_cols;
	arma::uword jj,uu;
	// arma::uword converge = 0;
	arma::uword iter = 0;
	double curr_LL = 0.0;
	double old_LL,new_LL,norm_GRAD,norm_FS;
	double new_LL_1,new_LL_2;
	
	// Initialize params, remember to unconstrain PHI
	arma::vec old_PARAMS = params0;
	arma::vec new_PARAMS = arma::zeros<arma::vec>(pp + 1);
	arma::vec curr_PARAMS = new_PARAMS;
	arma::vec curr_GRAD = new_PARAMS;
	arma::vec up_GRAD = new_PARAMS;
	arma::vec curr_FS = new_PARAMS;
	arma::vec up_FS = new_PARAMS;
	arma::mat curr_HESS = arma::zeros<arma::mat>(pp + 1,pp + 1);
	arma::vec new_params1 = arma::zeros<arma::vec>(pp + 1);
	arma::vec new_params2 = new_params1;
	
	while(iter < max_iter){
		old_LL = Rcpp_NB_reg_LL(Y,X,O,old_PARAMS);
		curr_GRAD = Rcpp_NB_reg_GRAD(Y,X,O,old_PARAMS);
		norm_GRAD = Rcpp_norm(curr_GRAD);
		up_GRAD = curr_GRAD / std::max(1.0,norm_GRAD);
		curr_HESS = Rcpp_NB_reg_HESS(Y,X,O,old_PARAMS);
		curr_FS = arma::inv_sympd(-1.0 * curr_HESS) * curr_GRAD;
		norm_FS = Rcpp_norm(curr_FS);
		up_FS = curr_FS / std::max(1.0,norm_FS);
		uu = 0;
		
		for(jj = 0; jj <= 20; jj++){
			new_params1 = old_PARAMS + up_GRAD / std::pow(4,jj);
			new_LL_1 = Rcpp_NB_reg_LL(Y,X,O,new_params1);
			new_params2 = old_PARAMS + up_FS / std::pow(4,jj);
			new_LL_2 = Rcpp_NB_reg_LL(Y,X,O,new_params2);
			new_LL = std::max(new_LL_1,new_LL_2);
			if( new_LL > old_LL ){
				if( new_LL == new_LL_1 ){
					old_PARAMS = new_params1;
				} else {
					old_PARAMS = new_params2;
				}
				old_LL = new_LL;
				uu = 1;
				break;
			}
		}
		
		if( uu == 0 ){
			break;
		}
		
		if(iter > 0){
			if( std::abs(curr_LL - old_LL) < eps && Rcpp_norm(curr_PARAMS - old_PARAMS) < eps ){
				curr_GRAD = Rcpp_NB_reg_GRAD(Y,X,O,old_PARAMS);
				norm_GRAD = Rcpp_norm(curr_GRAD);
				curr_HESS = Rcpp_NB_reg_HESS(Y,X,O,old_PARAMS);
				norm_FS = Rcpp_norm(arma::inv_sympd(-1.0 * curr_HESS) * curr_GRAD);
				if( norm_GRAD < eps && norm_FS < eps ){
					// converge = 1;
					break;
				}
			}
		}
		
		curr_PARAMS = old_PARAMS;
		curr_LL = old_LL;
		iter++;
	}
	
	PARAMS = old_PARAMS;
	
	// return R_NilValue;
	/*
	return Rcpp::List::create(
		Rcpp::Named("converge",converge),
		Rcpp::Named("LL",old_LL),
		Rcpp::Named("iter",iter),
		Rcpp::Named("norm_FS",norm_FS),
		Rcpp::Named("norm_GRAD",norm_GRAD),
		Rcpp::Named("PARAMS",Rcpp::NumericVector(old_PARAMS.begin(),old_PARAMS.end()))
		);
	*/
}

// [[Rcpp::export]]
Rcpp::List Rcpp_NB_reg_BFGS(const arma::vec& Y,
	const arma::mat& X,const arma::vec& O,
	const arma::vec& params0,
	const arma::uword& max_iter = 4e3,
	const double& eps = 1e-7,const bool& show = true){
	
	arma::uword num_params = params0.n_elem;
	arma::uword iter = 0;
	arma::uword jj,uu;
	arma::uword converge = 0;
	
	arma::vec xk = params0;
	arma::mat inv_Bk = arma::eye<arma::mat>(num_params,num_params);
	arma::vec curr_xk = arma::zeros<arma::vec>(num_params);
	arma::mat I_num_params = arma::eye<arma::mat>(num_params,num_params);
	arma::vec new_xk = arma::zeros<arma::vec>(num_params);
	arma::vec gr_k = arma::zeros<arma::vec>(num_params);
	arma::vec p_k = arma::zeros<arma::vec>(num_params);
	arma::vec s_k = arma::zeros<arma::vec>(num_params);
	arma::vec y_k = arma::zeros<arma::vec>(num_params);
	arma::mat ISYT = arma::zeros<arma::mat>(num_params,num_params);
	
	double old_LL,new_LL,inv_norm_p_k,tmp_alpha,ys;
	double fnscale = -1.0; // For maximization
	double curr_LL = 0.0;
	
	while(iter < max_iter){
		// Calculate Direction p_k
		gr_k = fnscale * Rcpp_NB_reg_GRAD(Y,X,O,xk);
		p_k = -1.0 * inv_Bk * gr_k;
		inv_norm_p_k = 1.0 / std::max(1.0,Rcpp_norm(p_k));
		
		// Line search for new xk
		uu = 0;
		old_LL = fnscale * Rcpp_NB_reg_LL(Y,X,O,xk);
		for(jj = 0; jj <= 30; jj++){
			tmp_alpha = inv_norm_p_k / std::pow(4,jj);
			new_xk = xk + tmp_alpha * p_k;
			new_LL = fnscale * Rcpp_NB_reg_LL(Y,X,O,new_xk);
			if(new_LL < old_LL){ // minimizing
				s_k = tmp_alpha * p_k;
				y_k = fnscale * Rcpp_NB_reg_GRAD(Y,X,O,new_xk) - gr_k;
				ys = arma::dot(y_k,s_k);
				if( ys > 0.0 ){
					if(show) printR_obj("Update inv_Bk");
					ISYT = I_num_params - (s_k * y_k.t()) / ys;
					inv_Bk = ISYT * inv_Bk * ISYT.t() + s_k * s_k.t() / ys;
					xk = new_xk;
					old_LL = new_LL;
					uu = 1;
					break;
				}
			}
		}
		
		if( uu == 0 ) { // aka no update
			if( Rcpp_norm(gr_k) > 1.0 ){
				if(show) printR_obj("Reset inv_Bk");
				inv_Bk = I_num_params;
			} else {
				if(show) printR_obj("Failed line search");
				break;
			}
		}
		
		// Check Convergence
                // Chong: the optimization will only stop if both conditions are met.
		if( iter > 0 ){
			if( std::abs(curr_LL - old_LL) < eps &&
				Rcpp_norm(curr_xk - xk) < eps ){
				gr_k = Rcpp_NB_reg_GRAD(Y,X,O,xk);
				if( Rcpp_norm(gr_k) < eps ){
					converge = 1;
					break;
				}
			}
		}
		
		curr_xk = xk;
		curr_LL = old_LL;
		iter++;
	}
	
	old_LL = Rcpp_NB_reg_LL(Y,X,O,xk);
	return Rcpp::List::create(
		Rcpp::Named("converge",converge),
		Rcpp::Named("LL",old_LL),
		Rcpp::Named("iter",iter),
		Rcpp::Named("norm_GRAD",Rcpp_norm(gr_k)),
		Rcpp::Named("PAR",Rcpp::NumericVector(xk.begin(),xk.end()))
		);
}











