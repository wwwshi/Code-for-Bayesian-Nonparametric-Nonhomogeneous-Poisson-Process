#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

//calculate VN
// [[Rcpp::export]]
NumericVector fn_calc_VN(int n_grid, double gamma, double lambda) {
	NumericVector VN(n_grid + 10L, 0);
	double b, r, m;
	for (int t = 1; t < n_grid + 11; t++) {
		r = R_NegInf;
		for (int k = t; k < n_grid + 101; k++) {
			NumericVector l1 = log(Range(k - t + 1, k));
			NumericVector l2 = log(Range(k * gamma, k * gamma + n_grid - 1));
			b = sum(l1) - sum(l2) + dpois(NumericVector::create(k - 1), lambda, true)[0];
			m = max(NumericVector::create(b, r));
			r = log(exp(r - m) + exp(b - m)) + m;
		}
		VN[t - 1] = r;
	}
	return(VN);
}

// [[Rcpp::export]]
void fn_remove_grid(int id_grid, IntegerVector z, IntegerVector K, IntegerVector N_k, 
	NumericVector lambda, int n_grid){
	int cur_z, i;

	cur_z = z[id_grid];
	z[id_grid] = (-1);
	--N_k[cur_z - 1]; //decrease the number of grids by 1 in the corresponding cluster.

	//If a grid is clustered by itself
	if (N_k[cur_z - 1] == 0) {
		for (i = 0; i < n_grid; i++) {
			if (z[i] > cur_z) {
				--z[i];
			}
		}

		for (i = cur_z; i < n_grid; i++) {
			N_k[i - 1] = N_k[i];
			lambda[i - 1] = lambda[i];
		}
		N_k[(n_grid - 1)] = 0;
		lambda[(n_grid - 1)] = 0.0;
		--K[0];
	}	
}

//update z
// [[Rcpp::export]]
void fn_z_update(NumericVector data, double g, double a, double b, NumericVector VN,
	IntegerVector z, IntegerVector K, IntegerVector N_k, NumericVector lambda, int n_grid){
	NumericVector prob(n_grid + 0L, 0);
	double lambda_k;
	NumericVector data_i;
	int i, k;

	for (i = 0; i < n_grid; i++)
	{
		fn_remove_grid(i, z, K, N_k, lambda, n_grid);
		std::fill(prob.begin(), prob.end(), 0);
		data_i = data[i];

		//print(z);
		//print(N_k);
		//print(lambda);
		//print(K);

		for (k = 0; k < K[0]; k++) {
			lambda_k = lambda[k];
			prob[k] = log(g + N_k[k]) + dpois(data_i, lambda_k, true)[0];
		}

		//Prob of starting a new cluster
		prob[K[0]] = log(g) - lfactorial(data_i)[0] + lgamma(data[i] + a) - lgamma(a) + 
			a * log(b) - (data[i] + a) * log(b + 1) + VN[K[0]] - VN[K[0] - 1];

		// Sample the cluster
		IntegerVector range = seq(0, K[0]);
		NumericVector tprob = prob[range];
		tprob = exp(tprob);
		//printf("i = %d \n", i);
		//print(tprob);
		IntegerVector ind = Rcpp::RcppArmadillo::sample(range, 1, true, tprob);
		k = ind[0];
		z[i] = k + 1;
		//printf("\n selected k= %d", k);
		
		if (k < K[0]) { //Join one existing cluster
			++N_k[k];
		} else { //Start a new cluster
			N_k[k] = 1;
			lambda[k] = rgamma(1, a, 1 / b)[0];
			++K[0];
		}

		//print(z);
		//print(N_k);
		//print(lambda);
		//print(K);
	}
}