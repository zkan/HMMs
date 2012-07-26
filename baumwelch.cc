/*
 *	Author: Kan Ouivirach
 *	File: baumwelch.cc
 *	Purpose: Implement HMM class: conventional Baum-Welch algorithm, incremental Baum-Welch algorithm, incremental Maximum-Likelihood algorithm
 *	Credit: Tapas Kanungo, kanungo@cfar.umd.edu for the conventional Baum-Welch algorithm
 */

#include <cstring>
#include <sstream>
#include <fstream>
#include "hmm.h"
#include "utility.h"

string itos( int i ) {
	stringstream s;
	s << i;
	return s.str();
}

void HMM::baumwelch( int T, int *O, double **alpha, double **beta, double **gamma, int *pniter, double *plogprobinit, double *plogprobfinal ) {
	const double DELTA = 1e-05;

//	int i, j, k;
	int l = 0;

	double logprobf;
	double numeratorA, denominatorA;
	double numeratorB, denominatorB;

	double ***xi, *scale;
	double delta, deltaprev, logprobprev;

	deltaprev = 10e-70;

	xi = allocXi( T, this->N );
	scale = (double *) calloc ( (unsigned) ( T + 1 ), sizeof( double ) );

	this->forwardWithScale( T, O, alpha, scale, 0 );
/*
	cout << "Alpha" << endl;
	for( int q = 1; q <= this->N; q++ ) {
		for( int w = 1; w <= T; w++ ) {
			cout << alpha[q][w] << " ";
		}
		cout << endl;
	}
//*/
	logprobf = this->getLoglik()/T;
//	cout << this->getLoglik() << endl;
//	exit( 1 );
	*plogprobinit = logprobf;
	this->backwardWithScale( T, O, beta, scale, 0 );
/*
	cout << "Beta:" << endl;
	for( int q = 1; q <= this->N; q++ ) {
		for( int w = 1; w <= T; w++ ) {
			cout << beta[q][w] << " ";
		}
		cout << endl;
	}
*/
	this->computeGamma( T, alpha, beta, gamma );
	this->computeXi( T, O, alpha, beta, xi );
	logprobprev = logprobf;
	
	do {
		// reestimate frequency of state i in time t = 1
		for( int i = 1; i <= this->N; i++ ) {
		//	this->pi[i] = 0.001 + 0.999 * gamma[1][i];
			this->pi[i] = gamma[1][i];
		}
		// reestimate transition matrix and symbol prob in each state 
		for( int i = 1; i <= this->N; i++ ) {
			denominatorA = 0.0;
			for( int t = 1; t <= T - 1; t++ ) {
				denominatorA += gamma[t][i];
			//	cout << "denoA: " << denominatorA << endl;
			}
		//	cout << "denoA: " << denominatorA << endl;
			for( int j = 1; j <= this->N; j++ ) {
				numeratorA = 0.0;
				for( int t = 1; t <= T - 1; t++ ) {
			//		cout << "xi: " << xi[t][i][j] << endl;
					numeratorA += xi[t][i][j];
				}
			//	cout << "numeratorA: " << numeratorA << endl;
			//	cout << "denominatorA: " << denominatorA << endl;
			//	this->A[i][j] = 0.001 + 0.999 * ( numeratorA/denominatorA );
		//		cout << 0.999*numeratorA/denominatorA << endl;
				this->A[i][j] = numeratorA/denominatorA;
			//	cout << numeratorA << " ";
			}
		//	cout << endl;

			denominatorB = denominatorA + gamma[T][i];
			for( int k = 1; k <= this->M; k++ ) {
				numeratorB = 0.0;
				for( int t = 1; t <= T; t++ ) {
					if( O[t] == k ) {
						numeratorB += gamma[t][i];
					}
				}
			//	this->B[i][k] = 0.001 + 0.999 * ( numeratorB/denominatorB );
				this->B[i][k] = numeratorB/denominatorB;
			}
		}
		this->forwardWithScale( T, O, alpha, scale, 0 );
		this->backwardWithScale( T, O, beta, scale, 0 );
		this->computeGamma( T, alpha, beta, gamma );
		this->computeXi( T, O, alpha, beta, xi );
		logprobf = this->getLoglik()/T;
	//	cout << "logprobf (after): " << logprobf << endl;

		// compute the difference between the previous log probability and the current one
		delta = logprobf - logprobprev;
//		cout << "logprobf: " << logprobf << ", logprobprev: " << logprobprev << endl;
		cout << "delta: " << delta << endl;
		logprobprev = logprobf;
		l++;	// count the number of iterations
	}
	while( delta > DELTA );	// if log likelihood does not change much, then exit 
//	while( l < 150 );

	*pniter = l;
	*plogprobfinal = logprobf;
}

void HMM::train_bw( char *dir_name, int T, int *O, double **alpha, double **beta, double **gamma, int *pniter, double *plogprobinit, double *plogprobfinal ) {
	const double DELTA = 1e-05;
	int l = 0;

	double logprobf = 0;

	double **numeratorA, *denominatorA;
	double **numeratorB, *denominatorB;

	double ***xi, *scale;
	double delta, logprobprev;

	xi = allocXi( T + 1, this->N + 1 );
	scale = ( double* ) calloc ( ( unsigned ) ( T + 1 ), sizeof( double ) );

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
	double *tmp_pi = ( double* ) calloc( ( unsigned ) ( this->N + 1 ), sizeof( double ) );
	for( int i = 1; i <= this->N; i++ ) {
		tmp_pi[i] = 0;
	}
	
	char **files = get_all_files( dir_name );
	if( !files ) {
	//	fprintf( stderr, "%s: %s: something went wrong\n", argv[0], argv[2] );
	//  return EXIT_FAILURE;
		fprintf( stderr, "Error: cannot get all files under the directory!\n" );
		exit( 1 );
	}
//*
	ostringstream os_seq_name_tmp;
	int NS = 0;	// number of sequences
	cout << "### Sequences ###" << endl;
	for( int i = 0; files[i]; ++i ) {
		os_seq_name_tmp.clear();
		os_seq_name_tmp.str( "" );
		// ignore . and ..
		if( files[i][0] != '.' ) {
			os_seq_name_tmp << dir_name << files[i];
			cout << os_seq_name_tmp.str() << endl;
			os_seq_name_tmp.clear();
			os_seq_name_tmp.str( "" );
			NS++;
		}
	}
	
	os_seq_name_tmp.clear();
	os_seq_name_tmp.str( "" );
	cout << "### Training process ###" << endl;
	for( int i = 0; files[i]; ++i ) {
		os_seq_name_tmp.str( "" );
		// ignore . and ..
		if( files[i][0] != '.' ) {
//			printf( "%s\n", files[i] );
			os_seq_name_tmp << dir_name << files[i];
			cout << os_seq_name_tmp.str() << endl;
	//		cout << dir_name << files[i] << endl;

			FILE *fp;
			fp = fopen( os_seq_name_tmp.str().c_str(), "r" );
			if( fp == NULL ) {
				fprintf( stderr, "Error: File %s not found\n", os_seq_name_tmp.str().c_str() );
				exit( 1 );
			}
			readSequence( fp, &T, &O );
			if( fp ) {
				fclose( fp );
				fp = NULL;
			}
			os_seq_name_tmp.str( "" );
			
			// ignore the sequence length of zero
			if( T == 0 ) {
				continue;
			}
		
			this->forwardWithScale( T, O, alpha, scale, 0 );
			this->backwardWithScale( T, O, beta, scale, 0 );
			this->computeGamma( T, alpha, beta, gamma );
			this->computeXi( T, O, alpha, beta, xi );

			logprobf += this->getLoglik()/T;	// using the log likelihood per observation

			for( int i = 1; i <= this->N; i++ ) {
				tmp_pi[i] += gamma[1][i];
			}
			double sumGamma = 0;
			double sumXi = 0;
			for( int i = 1; i <= this->N; i++ ) {
				sumGamma = 0;
				for( int t = 1; t <= T - 1; t++ ) {
					sumGamma += gamma[t][i];
				}
				denominatorA[i] += sumGamma;
				for( int j = 1; j <= this->N; j++ ) {
					sumXi = 0;
					for( int t = 1; t <= T - 1; t++ ) {
						sumXi += xi[t][i][j];
					}
					numeratorA[i][j] += sumXi;
				}
				sumGamma = 0;
				for( int t = 1; t <= T; t++ ) {
					sumGamma += gamma[t][i];
				}
				denominatorB[i] += sumGamma;

				for( int k = 1; k <= this->M; k++ ) {
					sumGamma = 0;
					for( int t = 1; t <= T; t++ ) {
						if( O[t] == k ) {
							sumGamma += gamma[t][i];
						}
					}
					numeratorB[i][k] += sumGamma;
				}
			}
		}
	}
	
//*/
	// read the entire training data (for Baum-Welch only) 
/*
	for( int n = 1; n <= NS; n++ ) {
	//	if( n == 4 ) continue;
		FILE *fpx;

		ostringstream os_seq_name_tmp;	// store the sequence name temporarily
		os_seq_name_tmp << os_seq_name.str() << n << ".seq";
//		cout << os_seq_name_tmp.str() << endl;

		fpx = fopen( os_seq_name_tmp.str().c_str(), "r" );
		if( fpx == NULL ) {
			fprintf( stderr, "Error: File %s not found\n", os_seq_name_tmp.str().c_str() );
			exit( 1 );
		}
		readSequence( fpx, &T, &O );
		if( fpx ) {
			fclose( fpx );
			fpx = NULL;
		}
		os_seq_name_tmp.seekp( 0 );
		
		this->forwardWithScale( T, O, alpha, scale, 0 );
		this->backwardWithScale( T, O, beta, scale, 0 );
		this->computeGamma( T, alpha, beta, gamma );
		this->computeXi( T, O, alpha, beta, xi );

		logprobf += this->getLoglik()/T;	// using the log likelihood per observation

		for( int i = 1; i <= this->N; i++ ) {
			tmp_pi[i] += gamma[1][i];
		}
		double sumGamma = 0;
		double sumXi = 0;
		for( int i = 1; i <= this->N; i++ ) {
			sumGamma = 0;
			for( int t = 1; t <= T - 1; t++ ) {
				sumGamma += gamma[t][i];
			}
			denominatorA[i] += sumGamma;
			for( int j = 1; j <= this->N; j++ ) {
				sumXi = 0;
				for( int t = 1; t <= T - 1; t++ ) {
					sumXi += xi[t][i][j];
				}
				numeratorA[i][j] += sumXi;
			//	cout << "x" << numeratorA[i][j] << " ";
			}
			sumGamma = 0;
			for( int t = 1; t <= T; t++ ) {
				sumGamma += gamma[t][i];
			}
			denominatorB[i] += sumGamma;
//			cout << denominatorB[i] << endl;

			for( int k = 1; k <= this->M; k++ ) {
				sumGamma = 0;
				for( int t = 1; t <= T; t++ ) {
					if( O[t] == k ) {
						sumGamma += gamma[t][i];
					}
				}
				numeratorB[i][k] += sumGamma;
			}
		}
	}
//*/
	*plogprobinit = logprobf/NS;
	logprobprev = logprobf/NS;

	do {
		// normalize Pi
		for( int i = 1; i <= this->N; i++ ) {
			this->pi[i] = tmp_pi[i]/NS;
		}
		for( int i = 1; i <= this->N; i++ ) {
			for( int j = 1; j <= this->N; j++ ) {
				this->A[i][j] = numeratorA[i][j]/denominatorA[i];
			}
			for( int k = 1; k <= this->M; k++ ) {
				this->B[i][k] = numeratorB[i][k]/denominatorB[i];
			}
		}

		// reset the log likelihood from the forward algorithm
		logprobf = 0;
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
		for( int i = 1; i <= this->N; i++ ) {
			tmp_pi[i] = 0;
		}
//*
		for( int i = 0; files[i]; ++i ) {
			os_seq_name_tmp.str( "" );
			// ignore . and ..
			if( files[i][0] != '.' ) {
				os_seq_name_tmp << dir_name << files[i];

				FILE *fp;
				fp = fopen( os_seq_name_tmp.str().c_str(), "r" );
				if( fp == NULL ) {
					fprintf( stderr, "Error: File %s not found\n", os_seq_name_tmp.str().c_str() );
					exit( 1 );
				}
				readSequence( fp, &T, &O );
				if( fp ) {
					fclose( fp );
					fp = NULL;
				}
				os_seq_name_tmp.str( "" );
				
				// ignore the sequence length of zero
				if( T == 0 ) {
					continue;
				}
		
				this->forwardWithScale( T, O, alpha, scale, 0 );
				this->backwardWithScale( T, O, beta, scale, 0 );
				this->computeGamma( T, alpha, beta, gamma );
				this->computeXi( T, O, alpha, beta, xi );

				logprobf += this->getLoglik()/T;	// using the log likelihood per observation

				for( int i = 1; i <= this->N; i++ ) {
					tmp_pi[i] += gamma[1][i];
				}
				double sumGamma = 0;
				double sumXi = 0;
				for( int i = 1; i <= this->N; i++ ) {
					sumGamma = 0;
					for( int t = 1; t <= T - 1; t++ ) {
						sumGamma += gamma[t][i];
					}
					denominatorA[i] += sumGamma;
					for( int j = 1; j <= this->N; j++ ) {
						sumXi = 0;
						for( int t = 1; t <= T - 1; t++ ) {
							sumXi += xi[t][i][j];
						}
						numeratorA[i][j] += sumXi;
					//	cout << "x" << numeratorA[i][j] << " ";
					}
					sumGamma = 0;
					for( int t = 1; t <= T; t++ ) {
						sumGamma += gamma[t][i];
					}
					denominatorB[i] += sumGamma;
		//			cout << denominatorB[i] << endl;

					for( int k = 1; k <= this->M; k++ ) {
						sumGamma = 0;
						for( int t = 1; t <= T; t++ ) {
							if( O[t] == k ) {
								sumGamma += gamma[t][i];
							}
						}
						numeratorB[i][k] += sumGamma;
					}
				}
			}
		}
//*/
/*
		// read the entire training data (for Baum-Welch only)
		for( int n = 1; n <= NS; n++ ) {

		//	if( n == 4 ) continue;
			FILE *fpxx;
	//		string ss;
	//		char tmp[1000] = "";
	//		char str[100] = "/home/zkan/Desktop/ECTI-workspace/data/seq/a.";
	
	//		ss = itos( n );
	//		strcpy( tmp, ss.c_str() );
	//		strcat( str, strcat( tmp, ".seq" ) );

			ostringstream os_seq_name_tmp;	// for storing the sequence name temporarily
			os_seq_name_tmp << os_seq_name.str() << n << ".seq";
	//		cout << os_seq_name_tmp.str() << endl;

			fpxx = fopen( os_seq_name_tmp.str().c_str(), "r" );
			if( fpxx == NULL ) {
				fprintf( stderr, "Error: File %s not found\n", os_seq_name_tmp.str().c_str() );
				exit( 1 );
			}
			os_seq_name_tmp.seekp( 0 );
			
			readSequence( fpxx, &T, &O );
			this->forwardWithScale( T, O, alpha, scale, 0 );
			this->backwardWithScale( T, O, beta, scale, 0 );
//			cout << "while reestimating" << endl;
			this->computeGamma( T, alpha, beta, gamma );
			this->computeXi( T, O, alpha, beta, xi );
		
			if( fpxx ) {
				fclose( fpxx );
				fpxx = NULL;
			}

			logprobf += this->getLoglik()/T;
//			cout << "logprobf: " << logprobf << endl;

			for( int i = 1; i <= this->N; i++ ) {
				tmp_pi[i] += gamma[1][i];
			}

			double sumGamma;
			double sumXi;
			for( int i = 1; i <= this->N; i++ ) {
				sumGamma = 0;
				for( int t = 1; t <= T - 1; t++ ) {
					sumGamma += gamma[t][i];
				}
				denominatorA[i] += sumGamma;
				for( int j = 1; j <= this->N; j++ ) {
					sumXi = 0;
					for( int t = 1; t <= T - 1; t++ ) {
						sumXi += xi[t][i][j];
					}
					numeratorA[i][j] += sumXi;
				}
				//denominatorB[i] = denominatorA[i] + gamma[T][i];
				sumGamma = 0;
				for( int t = 1; t <= T; t++ ) {
					sumGamma += gamma[t][i];
				}
				denominatorB[i] += sumGamma;

				for( int k = 1; k <= this->M; k++ ) {
					sumGamma = 0;
					for( int t = 1; t <= T; t++ ) {
						if( O[t] == k ) {
							sumGamma += gamma[t][i];
						}
					}
					numeratorB[i][k] += sumGamma;
				}
			}
		}
//*/
		delta = logprobf/NS - logprobprev;
	//	cout << "logprobf: " << logprobf/NS << ", logprobprev: " << logprobprev << endl;
		cout << "delta: " << delta << endl;
		logprobprev = logprobf/NS;
	//	cout << "logprobf: " << logprobf/NS << endl;
	//	cout << logprobf/NS << ", ";
		l++;	// count the number of iterations

	}
	while( delta > DELTA );
//	while( l < 150 );

	for( int i = 0; files[i]; ++i ) {
		free( files[i] );
	}
	free( files );

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
	this->logprobfinal = logprobf;
	this->NS = NS;

	*pniter = l;
	*plogprobfinal = logprobf/NS;
}



//*

//*/

/*
void HMM::baumWelchWithSS( int T, int *O, double **alpha, double **beta, double **gamma, int *pniter, double *plogprobinit, double *plogprobfinal, 
				double *pi_t, double **numeratorA_t, double *denominatorA_t, double **numeratorB_t, double *denominatorB_t, 
				int NS, double forgetting_factor, int actual_obs_length ) {

	const double DELTA = 0.0025;
	int l = 0;
	double logprobf = 0;

	double **numeratorA, *denominatorA;
	double **numeratorB, *denominatorB;

	double ***xi, *scale;
	double delta, logprobprev;

	xi = allocXi( T + 1, this->N + 5 );
	scale = (double *) calloc ( (unsigned) ( T + 1 ), sizeof( double ) );

	numeratorA = (double **) calloc( (unsigned) ( this->N + 100 ), sizeof( double* ) );
	denominatorA = (double *) calloc( (unsigned) ( this->N + 100 ), sizeof( double ) );
	for( int i = 1; i <= this->N + 100; i++ ) {
		numeratorA[i] = (double *) calloc( (unsigned) ( this->N + 100 ), sizeof( double* ) );
	}
	numeratorB = (double **) calloc( (unsigned) ( this->N + 100 ), sizeof( double* ) );
	denominatorB = (double *) calloc( (unsigned) ( this->N + 100 ), sizeof( double ) );
	for( int i = 1; i <= this->N + 100; i++ ) {
		numeratorB[i] = (double *) calloc( (unsigned) ( this->M + 100 ), sizeof( double* ) );
	}
	// Initialize numeratorA, denominatorA, numeratorB, denominatorB
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
	// Initialize temporary Pi
	double *tmp_pi = (double*) calloc( (unsigned) ( this->N + 1 ), sizeof( double ) );
	for( int i = 1; i <= this->N; i++ ) {
		tmp_pi[i] = 0;
	}
	// Read entire training data (for testing Baum-Welch only)
	for( int n = 1; n <= NS; n++ ) {

		FILE *fpx;
		string ss;
		char tmp[1000] = "";
		char str[100] = "/home/zkan/PhD/svn/src/HMMs/sequence/tmp/";

		ss = itos( n );
		strcpy( tmp, ss.c_str() );
		strcat( str, strcat( tmp, ".seq" ) );
		fpx = fopen( str, "r" );
		if( fpx == NULL ) {
			fprintf( stderr, "Error: File %s not found\n", str );
			exit( 1 );
		}

		readSequence( fpx, &T, &O );
		this->forwardWithScale( T, O, alpha, scale );
		this->backwardWithScale( T, O, beta, scale );
		this->computeGamma( T, alpha, beta, gamma );
		this->computeXi( T, O, alpha, beta, xi );
		
		fclose( fpx );

		logprobf += this->getLoglik();

		for( int i = 1; i <= this->N; i++ ) {
			tmp_pi[i] += gamma[1][i];
		//	tmp_pi[i] += gamma[1][i] + pi_t[i];
		}
		double sumGamma = 0;
		double sumXi = 0;
		for( int i = 1; i <= this->N; i++ ) {
			sumGamma = 0;
			for( int t = 1; t <= T - 1; t++ ) {
				sumGamma += gamma[t][i];
			}
			denominatorA[i] += sumGamma;
			for( int j = 1; j <= this->N; j++ ) {
				sumXi = 0;
				for( int t = 1; t <= T - 1; t++ ) {
					sumXi += xi[t][i][j];
				}
				numeratorA[i][j] += sumXi;
			//	cout << "x" << numeratorA[i][j] << " ";
			}
			//denominatorB[i] = denominatorA[i] + gamma[T][i];
			sumGamma = 0;
			for( int t = 1; t <= T; t++ ) {
				sumGamma += gamma[t][i];
			}
			denominatorB[i] += sumGamma;
//			cout << denominatorB[i] << endl;

			for( int k = 1; k <= this->M; k++ ) {
				sumGamma = 0;
				for( int t = 1; t <= T; t++ ) {
					if( O[t] == k ) {
						sumGamma += gamma[t][i];
					}
				}
				numeratorB[i][k] += sumGamma;
			}
		}
		// calculate only the new suff stat of the new sequence and add to the old one
		// 
	//	for( int i = 1; i <= this->N; i++ ) {
	//		denominatorA[i] += denominatorA_t[i];
	//		for( int j = 1; j <= this->N; j++ ) {
	//			numeratorA[i][j] += numeratorA_t[i][j];
	//		}
	//		denominatorB[i] += denominatorB_t[i];
	//		for( int k = 1; k <= this->M; k++ ) {
	//			numeratorB[i][k] += numeratorB_t[i][k];
	//		}
	//	}
		
	}
	*plogprobinit = logprobf/NS;
	logprobprev = logprobf/NS;

	do {
		// update the HMM parameters
		for( int i = 1; i <= this->N; i++ ) {
			this->pi[i] = tmp_pi[i]/NS;
		}
		for( int i = 1; i <= this->N; i++ ) {
			for( int j = 1; j <= this->N; j++ ) {
				this->A[i][j] = numeratorA[i][j]/denominatorA[i];
			}
			for( int k = 1; k <= this->M; k++ ) {
				this->B[i][k] = numeratorB[i][k]/denominatorB[i];
			}
		}

		// Reset Log likelihood from Forward algorithm
		logprobf = 0;
		// Initialize numeratorA, denominatorA, numeratorB, denominatorB
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
		// Initialize temporary Pi
		for( int i = 1; i <= this->N; i++ ) {
			tmp_pi[i] = 0;
		}
		// Read entire training data (for testing Baum-Welch only)
		for( int n = 1; n <= NS; n++ ) {
			FILE *fpxx;
			string ss;
			char tmp[1000] = "";
			char str[100] = "/home/zkan/PhD/svn/src/HMMs/sequence/tmp/";
	
			ss = itos( n );
			strcpy( tmp, ss.c_str() );
			strcat( str, strcat( tmp, ".seq" ) );

			fpxx = fopen( str, "r" );
			if( fpxx == NULL ) {
				fprintf( stderr, "Error: File %s not found\n", str );
				exit( 1 );
			}
			readSequence( fpxx, &T, &O );
			this->forwardWithScale( T, O, alpha, scale );
			this->backwardWithScale( T, O, beta, scale );
//			cout << "while reestimating" << endl;
			this->computeGamma( T, alpha, beta, gamma );
			this->computeXi( T, O, alpha, beta, xi );
		
//			cout << "after first acc" << endl;	
			fclose( fpxx );

			logprobf += this->getLoglik();
//			cout << "logprobf: " << logprobf << endl;

			for( int i = 1; i <= this->N; i++ ) {
				tmp_pi[i] += gamma[1][i];
			}
			double sumGamma;
			double sumXi;
			for( int i = 1; i <= this->N; i++ ) {
				sumGamma = 0;
				for( int t = 1; t <= T - 1; t++ ) {
					sumGamma += gamma[t][i];
				}
				denominatorA[i] += sumGamma;
				for( int j = 1; j <= this->N; j++ ) {
					sumXi = 0;
					for( int t = 1; t <= T - 1; t++ ) {
						sumXi += xi[t][i][j];
					}
					numeratorA[i][j] += sumXi;
				}
				//denominatorB[i] = denominatorA[i] + gamma[T][i];
				sumGamma = 0;
				for( int t = 1; t <= T; t++ ) {
					sumGamma += gamma[t][i];
				}
				denominatorB[i] += sumGamma;

				for( int k = 1; k <= this->M; k++ ) {
					sumGamma = 0;
					for( int t = 1; t <= T; t++ ) {
						if( O[t] == k ) {
							sumGamma += gamma[t][i];
						}
					}
					numeratorB[i][k] += sumGamma;
				}
			}
			// calculate only the new suff stat of the new sequence and add to the old one
			// 
			//	for( int i = 1; i <= this->N; i++ ) {
			//		denominatorA[i] += denominatorA_t[i];
			//		for( int j = 1; j <= this->N; j++ ) {
			//			numeratorA[i][j] += numeratorA_t[i][j];
			//		}
			//		denominatorB[i] += denominatorB_t[i];
			//		for( int k = 1; k <= this->M; k++ ) {
			//			numeratorB[i][k] += numeratorB_t[i][k];
			//		}
			//	}
		}
		delta = logprobf/NS - logprobprev;
	//	cout << "logprobf: " << logprobf/NS << ", logprobprev: " << logprobprev << endl;
		cout << "delta: " << delta << endl;
		logprobprev = logprobf/NS;
		cout << "logprobf: " << logprobf/NS << endl;
	//	cout << logprobf/NS << endl;
		l++;	// count the number of iterations

	}
	while( delta > DELTA );
//	while(true);

	*pniter = l;
	*plogprobfinal = logprobf/NS;

	// store the sufficient statistics
	// 
	for( int i = 1; i <= this->N; i++ ) {
	//	pi_t[i] = ( forgetting_factor * pi_t[i] ) + this->pi[i];
	//	pi_t[i] = tmp_pi[i];
		pi_t[i] += this->pi[i];
	}
	for( int i = 1; i <= this->N; i++ ) {
	//	denominatorA_t[i] = ( forgetting_factor * denominatorA_t[i] ) + denominatorA[i];
		denominatorA_t[i] = denominatorA[i];
		for( int j = 1; j <= this->N; j++ ) {
		//	numeratorA_t[i][j] = ( forgetting_factor * numeratorA_t[i][j] ) + numeratorA[i][j];
			numeratorA_t[i][j] = numeratorA[i][j];
		}
	}
	for( int i = 1; i <= this->N; i++ ) {
	//	denominatorB_t[i] = ( forgetting_factor * denominatorB_t[i] ) + denominatorB[i];
		denominatorB_t[i] = denominatorB[i];
		for( int j = 1; j <= this->M; j++ ) {
		//	numeratorB_t[i][j] = ( forgetting_factor * numeratorB_t[i][j] ) + numeratorB[i][j];
			numeratorB_t[i][j] = numeratorB[i][j];
		}
	}
}
//*/










