/*
 *	Author: Kan Ouivirach
 *	File: baumwelch.cc
 *	Purpose: Implement HMM class: conventional continuous Baum-Welch algorithm
 *	Credit: Tapas Kanungo, kanungo@cfar.umd.edu for the conventional Baum-Welch algorithm
 */

#include <cstring>
#include <sstream>
#include <fstream>
#include "hmm.h"
#include "utility.h"

void HMM::computeContXi( int T, double **O, double **alpha, double **beta, double ***xi ) {
	double sum;
	for( int t = 1; t <= T - 1; t++ ) {
		sum = 0.0;
		for( int i = 1; i <= this->N; i++ ) {
			for( int j = 1; j <= this->N; j++ ) {
				xi[t][i][j] = alpha[t][i] * beta[t+1][j] * (this->A[i][j]) * computeGaussian( j, O[t+1] );
				sum += xi[t][i][j];
			}
		}
		for( int i = 1; i <= this->N; i++ ) {
			for( int j = 1; j <= this->N; j++ ) {
				xi[t][i][j] /= sum;
			}
		}
	}
}

double ***HMM::allocXi( int T, int N ) {
	double ***xi;
	xi = (double ***)malloc( ( T + 5 ) * sizeof(double **) );
	xi--;
	for( int t = 1; t <= T; t++ ) {
		xi[t] = (double **) calloc( (unsigned) ( N + 5 ), sizeof(double*) );
		for( int i = 1; i <= N; i++ ) {
			xi[t][i] = (double *) calloc( (unsigned) ( N + 5 ), sizeof(double) );
		}
	}
	for( int t = 1; t <= T; t++ ) {
		for( int i = 1; i <= N; i++ ) {
			for( int j = 1; j <= N; j++ ) {
				xi[t][i][j] = 0;
			}
		}
	}
	return xi;
}

void HMM::train_bw_cont( int NS, int T, double **O, double **alpha, double **beta, double **gamma, int *pniter, double *plogprobinit, double *plogprobfinal ) {
//	cout << "training.." << endl;

	const double DELTA = 1e-20;
	int l = 0;

	double logprobf = 0;

	double **numeratorA, *denominatorA;
	double *denominatorB;
//	double *numeratorB;

	double ***xi, *scale;
	double delta, logprobprev;

	double **numerator_mean;
	double ***numerator_cov;
	Matrix m_tmp_mean( 1, this->dim );
	Matrix m_tmp_cov( this->dim, this->dim );
	/* numerator mean and cov */
	numerator_mean = (double **) calloc( (unsigned) ( this->N + 100 ), sizeof(double*) );
	for( int i = 1; i <= this->N + 100; i++ ) {
		numerator_mean[i] = (double *) calloc( (unsigned) ( this->dim + 100 ), sizeof(double) );
	}
	numerator_cov = (double ***) calloc( (unsigned) ( this->N + 100 ), sizeof(double**) );
	for( int i = 1; i <= this->N + 100; i++ ) {
		numerator_cov[i] = (double **)calloc( (unsigned) ( this->dim + 100 ), sizeof(double*) );
		for( int j = 1; j <= this->dim + 100; j++ ) {
			numerator_cov[i][j] = (double *)calloc( (unsigned) ( this->dim + 100 ), sizeof(double) );
		}
	}

	xi = allocXi( T + 1, this->N + 5 );
	scale = (double *) calloc ( (unsigned) ( T + 1 ), sizeof( double ) );

	numeratorA = (double **) calloc( (unsigned) ( this->N + 100 ), sizeof( double* ) );
	denominatorA = (double *) calloc( (unsigned) ( this->N + 100 ), sizeof( double ) );
	for( int i = 1; i <= this->N + 100; i++ ) {
		numeratorA[i] = (double *) calloc( (unsigned) ( this->N + 100 ), sizeof( double* ) );
	}
//	numeratorB = (double *) calloc( (unsigned) ( this->N + 100 ), sizeof( double ) );
	denominatorB = (double *) calloc( (unsigned) ( this->N + 100 ), sizeof( double ) );
	/* Initialize numeratorA, denominatorA, numeratorB, denominatorB */
	for( int i = 1; i <= this->N; i++ ) {
		denominatorA[i] = 0;
		for( int j = 1; j <= this->N; j++ ) {
			numeratorA[i][j] = 0;
		}
	}
	for( int i = 1; i <= this->N; i++ ) {
//		numeratorB[i] = 0;
		denominatorB[i] = 0;
	}
	/* Initialize temporary Pi */
	double *tmp_pi = (double*) calloc( (unsigned) (this->N + 5), sizeof( double ) );
	for( int i = 1; i <= this->N; i++ ) {
		tmp_pi[i] = 0;
	}
	/* Initialize numerator mean and cov */
	for( int i = 1; i <= this->N; i++ ) {
		for( int j = 1; j <= this->dim; j++ ) {
			numerator_mean[i][j] = 0;
		}
	}
	for( int i = 1; i <= this->N; i++ ) {
		for( int j = 1; j <= this->dim; j++ ) {
			for( int k = 1; k <= this->dim; k++ ) {
				numerator_cov[i][j][k] = 0;
			}
		}
	}
	/* Read entire training data (for testing Baum-Welch only) */
	for( int n = 1; n <= NS; n++ ) {
		FILE *fpx;
		string ss;
		char tmp[1000] = "";
	//	char str[100] = "sequence/x";
	//	char str[100] = "sequence/test/test";
	//	char str[100] = "sequence/test";
	//	char str[100] = "sequence/data_3s5d/data";
		char str[100] = "../data/model/all_sequences/cont/";
		//char str[100] = "sequence/all_cont/";

	//	cout << "n: " << n << endl;
//n = 7;
		ss = itos( n );
		strcpy( tmp, ss.c_str() );
		strcat( str, strcat( tmp, ".seq" ) );

		fpx = fopen( str, "r" );
		if( fpx == NULL ) {
			fprintf( stderr, "Error: File %s not found\n", str );
			exit( 1 );
		}

		readContSequence( fpx, &T, &O );
		this->forwardContWithScale( T, O, alpha, scale );
		this->backwardContWithScale( T, O, beta, scale );
		this->computeGamma( T, alpha, beta, gamma );
		this->computeContXi( T, O, alpha, beta, xi );
		fclose( fpx );

		logprobf += this->getLoglik();
	//	cout << "logprobf: " << logprobf << endl;
	
		for( int i = 1; i <= this->N; i++ ) {
			tmp_pi[i] += gamma[1][i];
		}
		double sumGamma = 0;
		double sumXi = 0;
		for( int i = 1; i <= this->N; i++ ) {
			sumGamma = 0;
			for( int t = 1; t <= T - 1; t++ ) {
				sumGamma += gamma[t][i];
			//	cout << gamma[t][i] << " ";
			}
			//cout << endl;
			denominatorA[i] += sumGamma;
			for( int j = 1; j <= this->N; j++ ) {
				sumXi = 0;
				for( int t = 1; t <= T - 1; t++ ) {
					sumXi += xi[t][i][j];
				}
				numeratorA[i][j] += sumXi;
			}
//			denominatorB[i] = denominatorA[i] + gamma[T][i];
			sumGamma = 0;
			for( int t = 1; t <= T; t++ ) {
				sumGamma += gamma[t][i];
			}
			denominatorB[i] += sumGamma;

			/* Reset tmp mean matrix */
			for( int d = 0; d < this->dim; d++ ) {
				m_tmp_mean( 0, d ) = 0;
			}
			for( int t = 1; t <= T; t++ ) {
				for( int d = 0; d < this->dim; d++ ) {
					m_tmp_mean( 0, d ) += gamma[t][i] * O[t][d+1];
				}
			}
			/* Store in numerator mean and cov */
			for( int d = 0; d < this->dim; d++ ) {
				numerator_mean[i][d+1] += m_tmp_mean( 0, d );
			}
		}
		for( int i = 1; i <= this->N; i++ ) {
			/* Reset tmp cov matrix */
			for( int row = 0; row < this->dim; row++ ) {
				for( int col = 0; col < this->dim; col++ ) {
					m_tmp_cov( row, col ) = 0;
				}
			}
			for( int t = 1; t <= T; t++ ) {
				Matrix m_T( 1, this->dim );
				for( int col = 0; col < this->dim; col++ ) {
					m_T( 0, col ) = O[t][col + 1] - this->meanMat[i][col + 1];
				}
				Matrix m = ~m_T;
				Matrix multi_result( this->dim, this->dim );
				multi_result = m * m_T;
				for( int row = 0; row < this->dim; row++ ) {
					for( int col = 0; col < this->dim; col++ ) {
						m_tmp_cov( row, col ) += gamma[t][i] * multi_result( row, col );
					}
				}
			}
			for( int row = 0; row < this->dim; row++ ) {
				for( int col = 0; col < this->dim; col++ ) {
					numerator_cov[i][row + 1][col + 1] += m_tmp_cov( row, col );
				}
			}
			
		}
	}
	*plogprobinit = logprobf/NS;
	logprobprev = logprobf/NS;
//	cout << "logprobf: " << logprobf/NS << endl;

//	printContHMM();
//	exit( 1 );

	do {
		/* Normalize Pi */
		for( int i = 1; i <= this->N; i++ ) {
			this->pi[i] = tmp_pi[i]/NS;
		}
		for( int i = 1; i <= this->N; i++ ) {
			for( int j = 1; j <= this->N; j++ ) {
				this->A[i][j] = numeratorA[i][j]/denominatorA[i];
			//	cout << denominatorA[i] << " ";
			}
			//cout << endl;
			/* Update mean first and then update cov */
			for( int d = 1; d <= this->dim; d++ ) {
				this->meanMat[i][d] = numerator_mean[i][d]/denominatorB[i];
			//	cout << denominatorB[i] << " ";
			}
			//cout << endl;
			/* Update cov */
			for( int row = 1; row <= this->dim; row++ ) {
				for( int col = 1; col <= this->dim; col++ ) {
					this->covMat[i][row][col] = numerator_cov[i][row][col]/denominatorB[i];
				}
			}
		}

//		printContHMM();
//		exit( 1 );

		/* Reset Log likelihood from Forward algorithm */
		logprobf = 0;
		/* Reset numeratorA, denominatorA, numeratorB, denominatorB */
		for( int i = 1; i <= this->N; i++ ) {
			denominatorA[i] = 0;
			for( int j = 1; j <= this->N; j++ ) {
				numeratorA[i][j] = 0;		
			}
		}
		for( int i = 1; i <= this->N; i++ ) {
		//	numeratorB[i] = 0;
			denominatorB[i] = 0;
		}
		/* Reset temporary Pi */
		for( int i = 1; i <= this->N; i++ ) {
			tmp_pi[i] = 0;
		}
		/* Reset numerator mean and cov */
		for( int i = 1; i <= this->N; i++ ) {
			for( int j = 1; j <= this->dim; j++ ) {
				numerator_mean[i][j] = 0;
			}
		}
		for( int i = 1; i <= this->N; i++ ) {
			for( int j = 1; j <= this->dim; j++ ) {
				for( int k = 1; k <= this->dim; k++ ) {
					numerator_cov[i][j][k] = 0;
				}
			}
		}
		/* Read entire training data (for testing Baum-Welch only) */
		for( int n = 1; n <= NS; n++ ) {
			FILE *fpxx;
			string ss;
			char tmp[1000] = "";
		//	char str[100] = "sequence/x";
		//	char str[100] = "sequence/test/test";
		//	char str[100] = "sequence/test";
			char str[100] = "../data/model/all_sequences/cont/";
			//char str[100] = "sequence/all_cont/";
			
			//cout << "n: " << n << endl;
//n = 7;
			ss = itos( n );
			strcpy( tmp, ss.c_str() );
			strcat( str, strcat( tmp, ".seq" ) );

			fpxx = fopen( str, "r" );
			if( fpxx == NULL ) {
				fprintf( stderr, "Error: File %s not found\n", str );
				exit( 1 );
			}
			readContSequence( fpxx, &T, &O );
			this->forwardContWithScale( T, O, alpha, scale );
			this->backwardContWithScale( T, O, beta, scale );
			this->computeGamma( T, alpha, beta, gamma );
			this->computeContXi( T, O, alpha, beta, xi );
			fclose( fpxx );

			logprobf += this->getLoglik();
			//cout << "logprobf: " << logprobf << endl;

			for( int i = 1; i <= this->N; i++ ) {
				tmp_pi[i] += gamma[1][i];
			//	cout << "gamma: " << gamma[1][i] << endl;
			}
/*
			for( int i = 1; i <= this->N; i++ ) {
				// Reset tmp mean matrix 
				for( int d = 0; d < this->dim; d++ ) {
					m_tmp_mean( 0, d ) = 0;
				}
				for( int j = 1; j <= this->N; j++ ) {
					for( int t = 1; t <= T - 1; t++ ) {
						numeratorA[i][j] += xi[t][i][j];
						denominatorA[i][j] += gamma[t][i];
					}
				}
				for( int t = 1; t <= T; t++ ) {
					for( int d = 0; d < this->dim; d++ ) {
						m_tmp_mean( 0, d ) += gamma[t][i] * O[t][d+1];
					//	cout << m_tmp_mean << endl;
					}
					denominatorB[i] += gamma[t][i];
				}
				// Store in numerator mean 
				for( int d = 0; d < this->dim; d++ ) {
					numerator_mean[i][d+1] += m_tmp_mean( 0, d );
				//	cout << numerator_mean[i][d+1] << endl;
				}
				//cout << endl;
			}
*/
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

				/* Reset tmp mean matrix */
				for( int d = 0; d < this->dim; d++ ) {
					m_tmp_mean( 0, d ) = 0;	
				}
				for( int t = 1; t <= T; t++ ) {
					for( int d = 0; d < this->dim; d++ ) {
						m_tmp_mean( 0, d ) += gamma[t][i] * O[t][d+1];
					}
				}
				/* Store in numerator mean and cov */
				for( int d = 0; d < this->dim; d++ ) {
					numerator_mean[i][d+1] += m_tmp_mean( 0, d );
				}
			}
			for( int i = 1; i <= this->N; i++ ) {
				/* Reset tmp cov matrix */
				for( int row = 0; row < this->dim; row++ ) {
					for( int col = 0; col < this->dim; col++ ) {
						m_tmp_cov( row, col ) = 0;
					}
				}
				for( int t = 1; t <= T; t++ ) {
					Matrix m_T( 1, this->dim );
					for( int col = 0; col < this->dim; col++ ) {
						m_T( 0, col ) = O[t][col + 1] - this->meanMat[i][col + 1];
					}
					Matrix m = ~m_T;
					Matrix multi_result( this->dim, this->dim );
					multi_result = m * m_T;
					for( int row = 0; row < this->dim; row++ ) {
						for( int col = 0; col < this->dim; col++ ) {
							m_tmp_cov( row, col ) += gamma[t][i] * multi_result( row, col );
						}
					}
				}
				for( int row = 0; row < this->dim; row++ ) {
					for( int col = 0; col < this->dim; col++ ) {
						numerator_cov[i][row + 1][col + 1] += m_tmp_cov( row, col );
					}
				}
				
			}
		}
		delta = logprobf/NS - logprobprev;
	//	printf( "%.30lf - %.30lf\n", logprobf/NS, logprobprev );
		cout << "delta: " << delta << endl;
		//cout << logprobf/NS << ", ";

		logprobprev = logprobf/NS;
		l++;	// count the number of iterations

		if( l > 150 ) break;
	}
	while( delta > DELTA );

	*pniter = l;
	*plogprobfinal = logprobf/NS;
}


