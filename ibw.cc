/*
 *	Author: Kan Ouivirach
 *	File: baumwelch.cc
 *	Purpose: Implement HMM class: incremental Baum-Welch algorithm
 *	Credit: Tapas Kanungo, kanungo@cfar.umd.edu for the conventional Baum-Welch algorithm
 */

#include <cstring>
#include <sstream>
#include <fstream>
#include "hmm.h"
#include "utility.h"

// only used in the IBW algorithm
// 
void HMM::computeGamma( int T, double **alpha, double **gamma ) {
	double denominator;
	for( int t = 1; t <= T; t++ ) {
		denominator = 0.0;
		for( int i = 1; i <= this->N; i++ ) {
			gamma[t][i] = alpha[t][i];
			denominator += gamma[t][i];
		//	cout << i << "-denom: " << denominator << endl;
		}
	//	cout << endl;
		for( int i = 1; i <= this->N; i++ ) {
			gamma[t][i] = gamma[t][i]/denominator;
		}
	}
/*
	cout << "GAMMA" << endl;
	for( int t = 1; t <= T; t++ ) {
		for( int i = 1; i <= this->N; i++ ) {
			cout << gamma[t][i] << " ";
		}
		cout << endl;
	}
//*/
}

// only used in the IBW algorithm
// 
void HMM::computeXi( int T, int *O, double **alpha, double ***xi ) {
//	cout << "inside Xi" << endl;
	double sum;

	for( int t = 1; t <= T - 1; t++ ) {
		sum = 0.0;
		for( int i = 1; i <= this->N; i++ ) {
			for( int j = 1; j <= this->N; j++ ) {
				xi[t][i][j] = alpha[t][i] * this->A[i][j] * this->B[j][O[t+1]];
				sum += xi[t][i][j];
			}
		}
		//cout << "N: " << this->N << endl;
		for( int i = 1; i <= this->N; i++ ) {
			for( int j = 1; j <= this->N; j++ ) {
				xi[t][i][j] /= sum;
			//	cout << "xi: " << xi[t][i][j] << " ";
			}
		//	cout << endl;
		}
	}
}

double ***HMM::allocXi( int T, int N ) {
//	cout << "inside allocXi" << endl;
	double ***xi;
//cout << "??" << endl;
	xi = (double ***)malloc( ( T + 5 ) * sizeof(double **) );
	xi--;
//cout << "??" << endl;
//	cout << "T: " << T << endl;
	for( int t = 1; t <= T; t++ ) {
		xi[t] = (double **) calloc( (unsigned) ( N + 5 ), sizeof(double*) );
		for( int i = 1; i <= N; i++ ) {
			xi[t][i] = (double *) calloc( (unsigned) ( N + 5 ), sizeof(double) );
//			cout << t << ", " << i << endl;
		}
	}

	for( int t = 1; t <= T; t++ ) {
		for( int i = 1; i <= N; i++ ) {
			for( int j = 1; j <= N; j++ ) {
			//	cout << "j: " << j << endl;
				xi[t][i][j] = 0;
//				cout << xi[t][i][j] << " ";
			}
//			cout << endl;
		}
//		cout << endl;
	}
//	cout << endl;
//	cout << "before return allocXi" << endl;
	return xi;
}


void HMM::train_ibw( int T, int *O, double **alpha, double **beta, double **gamma, int *pniter, double *plogprobinit, double *plogprobfinal ) {
	const double DELTA = 0.001;

	int l = 0;	// the number of iterations

	double logprobf;

	double numeratorA, denominatorA;
	double numeratorB, denominatorB;

	double ***xi, *scale;
	double delta, logprobprev;

	xi = allocXi( T, this->N );
	scale = (double *) calloc ( (unsigned) ( T + 1 ), sizeof( double ) );

	// For each observation, estimate the parameters.
	// 
	for( int t = 2; t <= T; t++ ) {
		cout << "T: " << t << endl;

		this->forwardWithScale( t, O, alpha, scale, 0 );
		logprobf = this->getLoglik();
		*plogprobinit = logprobf;

		// For this algorithm, the beta cannot be computed incrementally since there are no observations 
		// after the current time are available.
		// this->backwardWithScale( T, O, beta, scale );
		// Set the beta value approximately to be 1 as suggested by Stenger et al. (2001)
		//
		for( int tmp = 1; tmp <= t; tmp++ ) {
			for( int i = 1; i <= this->N; i++ ) {
				beta[tmp][i] = 1;
			}
		}
		this->computeGamma( t, alpha, gamma );
		this->computeXi( t, O, alpha, xi );
		logprobprev = logprobf;

		do {
			// No need to reestimate the pi
			//
		//	for( int i = 1; i <= this->N; i++ ) {
			//	this->pi[i] = 0.001 + 0.999 * gamma[1][i];
		//		this->pi[i] = gamma[1][i];
		//	}
			// reestimate transition matrix and symbol prob in each state 
			for( int i = 1; i <= this->N; i++ ) {
				denominatorA = 0.0;
				for( int tmp = 1; tmp <= t - 1; tmp++ ) {
					denominatorA += gamma[tmp][i];
				//	cout << "denoA: " << denominatorA << endl;
				}
				for( int j = 1; j <= this->N; j++ ) {
					numeratorA = 0.0;
					for( int tmp = 1; tmp <= t - 2; tmp++ ) {
						numeratorA += gamma[tmp][i];
					}
					numeratorA = ( this->A[i][j] * numeratorA ) + xi[t-1][i][j];
					this->A[i][j] = numeratorA/denominatorA;
				//	cout << numeratorA << " ";
				}
			//	cout << endl;

				denominatorB = denominatorA + gamma[t][i];
				numeratorB = 0.0;
				for( int tmp = 1; tmp <= t - 1; tmp++ ) {
					numeratorB += gamma[tmp][i];
				}
				for( int k = 1; k <= this->M; k++ ) {
					numeratorB = this->B[i][k] * numeratorB;
					if( O[t] == k ) {
						numeratorB += gamma[t][i];
					}
				//	this->B[i][k] = 0.001 + 0.999 * ( numeratorB/denominatorB );
					this->B[i][k] = numeratorB/denominatorB;
				}
			}
			this->forwardWithScale( t, O, alpha, scale, 0 );
		//	this->backwardWithScale( T, O, beta, scale );
			for( int tmp = 1; tmp <= t; tmp++ ) {
				for( int i = 1; i <= this->N; i++ ) {
					beta[tmp][i] = 1;
				}
			}
			this->computeGamma( t, alpha, gamma );
			this->computeXi( t, O, alpha, xi );
			logprobf = this->getLoglik();

		//	cout << "logprobf (after): " << logprobf << endl;

			// compute the difference between the previous log probability and the current one
			delta = logprobf - logprobprev;
	//		cout << "logprobf: " << logprobf << ", logprobprev: " << logprobprev << endl;
			cout << "delta: " << delta << endl;
			logprobprev = logprobf;
			l++;	// count the number of iterations

			
			// Debugging: Print the model for each time step to see how the model changes
			// 
			char c;
			cin >> c;
			printHMM();
		}
		while( delta > DELTA );	// if log likelihood does not change much, then exit
		*pniter = l;		
	}
	*plogprobfinal = logprobf;
}



