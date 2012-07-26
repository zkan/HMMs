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

void HMM::train_iml( int T, int *O, double **alpha, double **beta, double **gamma, int niter, double *plogprobinit, 
				double *plogprobfinal, double *pi_t, double **numeratorA_t, double *denominatorA_t, double **numeratorB_t, double *denominatorB_t, 
				int NS, ostringstream &os_seq_name, double forgetting_factor, int actual_obs_length, int num_of_seqs_trained, double logprob_trained ) {
	int l = 0;
	double logprobf = 0;
	
//	cout << "num_of_seqs_trained: " << num_of_seqs_trained << endl;
//	cout << "logprob_trained: " << logprob_trained << endl;
	
//	exit( 1 );
	
//	cout << os_seq_name.str() << endl; 
//	if( this->getNS() == 1 ) exit( 1 );

  const double DELTA = 1e-05;

	double **numeratorA, *denominatorA;
	double **numeratorB, *denominatorB;

	double ***xi, *scale;
	double delta, logprobprev;

	xi = allocXi( T + 1, this->N + 1 );
	scale = (double *) calloc ( (unsigned) ( T + 1 ), sizeof( double ) );
	for( int i = 1; i <= T; i++ ) {
		scale[i] = 0;
	}

	numeratorA = ( double** ) calloc( ( unsigned ) ( this->N + 1 ), sizeof( double* ) );
	for( int i = 1; i <= this->N + 1; i++ ) {
		numeratorA[i] = ( double* ) calloc( ( unsigned ) ( this->N + 1 ), sizeof( double ) );
	}
	denominatorA = ( double* ) calloc( ( unsigned ) ( this->N + 1 ), sizeof( double ) );
	numeratorB = ( double** ) calloc( ( unsigned ) ( this->N + 1 ), sizeof( double* ) );
	for( int i = 1; i <= this->N + 1; i++ ) {
		numeratorB[i] = ( double* ) calloc( ( unsigned ) ( this->M + 1 ), sizeof( double ) );
	}
	denominatorB = ( double* ) calloc( ( unsigned ) ( this->N + 1 ), sizeof( double ) );
	// initialize numeratorA, denominatorA, numeratorB, and denominatorB
	for( int i = 1; i <= this->N; i++ ) {
		denominatorA[i] = 0;
		for( int j = 1; j <= this->N; j++ ) {
			numeratorA[i][j] = 0;
		}
	}
	for( int i = 1; i <= this->N; i++ ) {
		denominatorB[i] = 0;
		for( int j = 1; j <= this->M; j++ ) {
			numeratorB[i][j] = 0;
		}
	}
	// initialize temporary Pi
	double *tmp_pi = (double*) calloc( (unsigned) ( this->N + 1 ), sizeof( double ) );
	for( int i = 1; i <= this->N; i++ ) {
		tmp_pi[i] = 0;
	}
	// --------------------------------------------------------------

	// set T again
	T = actual_obs_length;
//	cout << "T: " << T << endl; exit( 1 );

	// ignore the sequence length of zero
	if( T == 0 ) {
		exit( 1 );
	}

	this->forwardWithScale( T, O, alpha, scale, 0 );
	this->backwardWithScale( T, O, beta, scale, 0 );
	this->computeGamma( T, alpha, beta, gamma );
	this->computeXi( T, O, alpha, beta, xi );

	// E-step: calculate the parameters, and add up with the sufficient statistics
	// add the sufficient statistic to Pi
	for( int i = 1; i <= this->N; i++ ) {
//		if( num_of_seqs_trained != 0 ) {
//			tmp_pi[i] = gamma[1][i] + pi_t[i]/num_of_seqs_trained;
//		}
//		else {
			tmp_pi[i] = gamma[1][i] + pi_t[i];
//		}
	}
	double sumGamma = 0;
	double sumXi = 0;
	// add the sufficient statistics to A and B 
	for( int i = 1; i <= this->N; i++ ) {
		sumGamma = 0;
		for( int t = 1; t <= T - 1; t++ ) {
			sumGamma += gamma[t][i];
		}
		denominatorA[i] = sumGamma + ( forgetting_factor * denominatorA_t[i] );
		for( int j = 1; j <= this->N; j++ ) {
			sumXi = 0;
			for( int t = 1; t <= T - 1; t++ ) {
				sumXi += xi[t][i][j];
			}
			numeratorA[i][j] = sumXi + ( forgetting_factor * numeratorA_t[i][j] );
		}
		sumGamma = 0;
		for( int t = 1; t <= T; t++ ) {
			sumGamma += gamma[t][i];
		}
		denominatorB[i] = sumGamma + ( forgetting_factor * denominatorB_t[i] );

		for( int k = 1; k <= this->M; k++ ) {
			sumGamma = 0;
			for( int t = 1; t <= T; t++ ) {
				if( O[t] == k ) {
					sumGamma += gamma[t][i];
				}
			}
			numeratorB[i][k] = sumGamma + ( forgetting_factor * numeratorB_t[i][k] );
		}
	}

	logprobf = this->getLoglik()/T;
//	cout << logprobf << endl; exit( 1 );
//	cout << this->getLogProbFinal() << endl; cout << this->getNS() << endl; if( this->getNS() == 1 ) exit( 1 );
	if( num_of_seqs_trained != 0 ) { 
		logprobf += logprob_trained/num_of_seqs_trained;
	//	cout << this->getLoglik() << " + " << this->getLogProbFinal()/this->getNS() << endl; if( this->getNS() == 1 ) exit( 1 );
	}
	*plogprobinit = logprobf;
	logprobprev = logprobf;

	do {
		// M-step: update the HMM parameters using the sufficient statistics
		// normalize Pi
	  double dSumGamma = 0;
	  for( int i = 1; i <= this->N; i++ ) {
	    dSumGamma += tmp_pi[i];
	  }
		for( int i = 1; i <= this->N; i++ ) {
//			if( num_of_seqs_trained != 0 ) {
		  //this->pi[i] = tmp_pi[i]/( num_of_seqs_trained + 1 );
		  this->pi[i] = tmp_pi[i]/dSumGamma;
//			}
//			else {
//				this->pi[i] = tmp_pi[i];
//			}
		}
		for( int i = 1; i <= this->N; i++ ) {
			for( int j = 1; j <= this->N; j++ ) {
				this->A[i][j] = numeratorA[i][j]/denominatorA[i];
			}
			for( int k = 1; k <= this->M; k++ ) {
				this->B[i][k] = numeratorB[i][k]/denominatorB[i];
			}
		}
		// to avoid some parameters to receive 0, add a constant e in pi, A, and B.
//*
		const double e = 1e-80;
		for( int i = 1; i <= this->N; i++ ) {
			this->pi[i] += e;
			for( int j = 1; j <= this->N; j++ ) {
				this->A[i][j] += e;
			}
		}
		for( int i = 1; i <= this->N; i++ ) {
			for( int j = 1; j <= this->M; j++ ) {
				this->B[i][j] += e;
			}
		}
//*/
		// ignore the sequence length of zero
		if( T == 0 ) {
			exit( 1 );
		}
	
		this->forwardWithScale( T, O, alpha, scale, 0 );
		this->backwardWithScale( T, O, beta, scale, 0 );
		this->computeGamma( T, alpha, beta, gamma );
		this->computeXi( T, O, alpha, beta, xi );

		// E-step: calculate the parameters, and add up with the sufficient statistics
		// add the sufficient statistic to Pi
		for( int i = 1; i <= this->N; i++ ) {
//			if( num_of_seqs_trained != 0 ) {
//				tmp_pi[i] = gamma[1][i] + pi_t[i]/num_of_seqs_trained;
//			}
//			else {
				tmp_pi[i] = gamma[1][i] + pi_t[i];
//			}
		}
		double sumGamma = 0;
		double sumXi = 0;
		// add the sufficient statistics to A and B 
		for( int i = 1; i <= this->N; i++ ) {
			sumGamma = 0;
			for( int t = 1; t <= T - 1; t++ ) {
				sumGamma += gamma[t][i];
			}
			denominatorA[i] = sumGamma + ( forgetting_factor * denominatorA_t[i] );
			for( int j = 1; j <= this->N; j++ ) {
				sumXi = 0;
				for( int t = 1; t <= T - 1; t++ ) {
					sumXi += xi[t][i][j];
				}
				numeratorA[i][j] = sumXi + ( forgetting_factor * numeratorA_t[i][j] );
			}
			sumGamma = 0;
			for( int t = 1; t <= T; t++ ) {
				sumGamma += gamma[t][i];
			}
			denominatorB[i] = sumGamma + ( forgetting_factor * denominatorB_t[i] );

			for( int k = 1; k <= this->M; k++ ) {
				sumGamma = 0;
				for( int t = 1; t <= T; t++ ) {
					if( O[t] == k ) {
						sumGamma += gamma[t][i];
					}
				}
				numeratorB[i][k] = sumGamma + ( forgetting_factor * numeratorB_t[i][k] );
			}
		}

		logprobf = this->getLoglik()/T;
		if( num_of_seqs_trained != 0 ) { 
			logprobf += logprob_trained/num_of_seqs_trained;
		}

		delta = logprobf - logprobprev;
	//	cout << "delta: " << delta << endl;
	//	cout << "logprof: " << logprobf << endl;
	//	cout << logprobf/( NS + num_of_seqs_trained ) << ", ";
		logprobprev = logprobf;
		l++;	// count the number of iterations
		
		if(delta < DELTA) {
		  break;
		}

	}
//	while( delta > DELTA );
	while( l < niter );

	*plogprobfinal = logprobf/( num_of_seqs_trained + 1 );

	// store the sufficient statistics
	for( int i = 1; i <= this->N; i++ ) {
		this->pi_t[i] = tmp_pi[i];
	}
	for( int i = 1; i <= this->N; i++ ) {
		this->denominatorA_t[i] = denominatorA[i];
		for( int j = 1; j <= this->N; j++ ) {
			this->numeratorA_t[i][j] = numeratorA[i][j];
		}
	}
	for( int i = 1; i <= this->N; i++ ) {
		this->denominatorB_t[i] = denominatorB[i];
		for( int j = 1; j <= this->M; j++ ) {
			this->numeratorB_t[i][j] = numeratorB[i][j];
		}
	}
	this->logprobfinal += logprobf;
	this->NS++;

//	free( numeratorA );
//	free( denominatorA );
//	free( numeratorB );
//	free( denominatorB );
}
