/*
 *	Author: Kan Ouivirach
 *	File: backward.cc
 *	Purpose: Implement HMM class: Backward algorithm
 *	Credit: Tapas Kanungo, kanungo@cfar.umd.edu
 */

#include "hmm.h"

void HMM::backward( int T, int *O, double **beta ) {
	double sum;
/*
	cout << "Observation: ";
	for( int i = 1; i <= T; i++ ) {
		cout << O[i] << " ";
	}
	cout << endl;
*/
	// 1. Initialization 
	for( int i = 1; i <= this->N; i++ ) {
		beta[T][i] = 1.0;
	}
	// 2. Induction 
	for( int t = T - 1; t >= 1; t-- ) {
		for( int i = 1; i <= this->N; i++ ) {
			sum = 0.0;
			for( int j = 1; j <= this->N; j++ ) {
				sum += this->A[i][j] * this->B[j][O[t+1]] * beta[t+1][j];
			}
			beta[t][i] = sum;
		}
	}
	// 3. Termination 
	this->loglik = 0.0;
	for( int i = 1; i <= this->N; i++ ) {
		this->loglik += beta[1][i];
	}
	this->loglik = log( this->loglik );
}

void HMM::backwardWithScale( int T, int *O, double **beta, double *scale, int type = 0 ) {
	double sum;
/*
	cout << "Observation: ";
	for( int i = 1; i <= T; i++ ) {
		cout << O[i] << " ";
	}
	cout << endl;
*/
	// 1. Initialization 
//	if( scale[T] == 0 ) {
//		scale[T] = 1e-80;
//	}
	for( int i = 1; i <= this->N; i++ ) {
		beta[T][i] = 1.0/scale[T];
	}
	// 2. Induction 
	for( int t = T - 1; t >= 1; t-- ) {
		for( int i = 1; i <= this->N; i++ ) {
			sum = 0.0;
			for( int j = 1; j <= this->N; j++ ) {
				if( beta[t+1][j] == 0 ) {
					beta[t+1][j] = 1e-80;
				}
				sum += this->A[i][j] * this->B[j][O[t+1]] * beta[t+1][j];
			//	cout << this->A[i][j] << ", " << this->B[j][O[t+1]] << ", " << beta[t+1][j] << endl;
			//	cout << this->A[i][j] * this->B[j][O[t+1]] * beta[t+1][j] << endl;
			}
		//	if( sum == 0 ) {
		//		sum = 1e-80;
		//	}
			beta[t][i] = sum/scale[t];
		//	cout << "sum: " << sum << endl;
		}
	}
}

void HMM::backwardContWithScale( int T, double **O, double **beta, double *scale ) {
	double sum;
/*
	cout << "Observation: ";
	for( int i = 1; i <= T; i++ ) {
		cout << O[i] << " ";
	}
	cout << endl;
*/
	// 1. Initialization 
	if( scale[T] == 0 ) {
		scale[T] = 1e-80;
	}
	for( int i = 1; i <= this->N; i++ ) {
		beta[T][i] = 1.0/scale[T];
	}
	// 2. Induction 
	for( int t = T - 1; t >= 1; t-- ) {
		for( int i = 1; i <= this->N; i++ ) {
			sum = 0.0;
			for( int j = 1; j <= this->N; j++ ) {
//				sum += this->A[i][j] * this->B[j][O[t+1]] * beta[t+1][j];
				sum += this->A[i][j] * computeGaussian( j, O[t+1] ) * beta[t+1][j];
			}
			beta[t][i] = sum/scale[t];
		}
	}
}



