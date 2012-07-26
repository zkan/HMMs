#include "hmm.h"

string func_itos( int i ) {
	stringstream s;
	s << i;
	return s.str();
}

int main( int argc, char **argv ) {
	FILE *fp;
	int T;
	int NS = 1;
	HMM *hmm;
	int *O;		/* observations */
	double **alpha, **beta, **gamma;
	double logprobinit, logprobfinal;
	int niter;	/* number of iterations */

	// forgetting factor for sufficient statistics
	double forgetting_factor = 0.8;

	int c;
	extern char *optarg;
	extern int optind, opterr, optopt;

	while( ( c = getopt( argc, argv, "N:f:" ) ) != EOF ) {
		switch( c ) {
			case 'N':
				sscanf( optarg, "%d", &NS );
				break;
			case 'f':
				sscanf( optarg, "%lf", &forgetting_factor );
				break;
		}
	}

	if( ( ( argc - optind ) != 1 ) || ( argc == 1 ) ) {
		cout << "Usage error" << endl;
		cout << "Usage: xxx <model.hmm>" << endl;
		cout << "Option: " << endl;
		cout << "-N <# of sequences>" << endl;
		exit( 1 );
	}
	
	fp = fopen( argv[argc - 1], "r" );
	if( fp == NULL ) {
		fprintf( stderr, "Error: File %s not found\n", argv[1] );
		exit( 1 );
	}

	hmm = new HMM();
	hmm->readHMMWithSuffStats( fp );

	double *pi_t = (double *) calloc( (unsigned) ( hmm->getNumStates() + 5 ), sizeof( double ) );
	double **numeratorA_t, *denominatorA_t;
	double **numeratorB_t, *denominatorB_t;
	numeratorA_t = (double **) calloc( (unsigned) ( hmm->getNumStates() + 5 ), sizeof( double* ) );
	denominatorA_t = (double *) calloc( (unsigned) ( hmm->getNumStates() + 5 ), sizeof( double ) );
	for( int i = 1; i <= hmm->getNumStates() + 5; i++ ) {
		numeratorA_t[i] = (double *) calloc( (unsigned) ( hmm->getNumStates() + 5 ), sizeof( double* ) );
	}
	numeratorB_t = (double **) calloc( (unsigned) ( hmm->getNumStates() + 5 ), sizeof( double* ) );
	denominatorB_t = (double *) calloc( (unsigned) ( hmm->getNumStates()+ 5 ), sizeof( double ) );
	for( int i = 1; i <= hmm->getNumStates() + 5; i++ ) {
		numeratorB_t[i] = (double *) calloc( (unsigned) ( hmm->getNumSymbols() + 5 ), sizeof( double* ) );
	}
	// Initialize numeratorA, denominatorA, numeratorB, denominatorB
	for( int i = 1; i <= hmm->getNumStates(); i++ ) {
		denominatorA_t[i] = 0;
		for( int j = 1; j <= hmm->getNumStates(); j++ ) {
			numeratorA_t[i][j] = 0;
		}
	}
	for( int i = 1; i <= hmm->getNumStates(); i++ ) {
		denominatorB_t[i] = 0;
		for( int j = 1; j <= hmm->getNumSymbols(); j++ ) {
			numeratorB_t[i][j] = 0;
		}
	}
	hmm->printSuffStats();
//return 1;
//*
	pi_t = hmm->getPi_t();
	numeratorA_t = hmm->getNumeratorA_t();
	denominatorA_t = hmm->getDenominatorA_t();
	numeratorB_t = hmm->getNumeratorB_t();
	denominatorB_t = hmm->getDenominatorB_t();
//*/
	fp = fopen( argv[argc - 1], "r" );
	if( fp == NULL ) {
		fprintf( stderr, "Error: File %s not found\n", argv[1] );
		exit( 1 );
	}
	hmm->readHMM( fp );
//	hmm->printHMM();
	
	// print the sufficient statistics
	cout << endl << "## SUFFICEINT STATISTICS ##" << endl;
	cout << "T= " << NS << endl;
	cout << "Pi: " << endl;
	for( int i = 1; i <= hmm->getNumStates(); i++ ) {
		cout << pi_t[i] << " ";
	}
	cout << endl;
	cout << "Numerator A: " << endl;
	for( int i = 1; i <= hmm->getNumStates(); i++ ) {
		for( int j = 1; j <= hmm->getNumStates(); j++ ) {
			cout << numeratorA_t[i][j] << " ";
		}
		cout << endl;
	}
	cout << "Denominator A: " << endl;
	for( int i = 1; i <= hmm->getNumStates(); i++ ) {
		cout << denominatorA_t[i] << " ";
	}
	cout << endl;
	cout << "Numerator B: " << endl;
	for( int i = 1; i <= hmm->getNumStates(); i++ ) {
		for( int j = 1; j <= hmm->getNumSymbols(); j++ ) {
			cout << numeratorB_t[i][j] << " ";
		}
		cout << endl;
	}
	cout << "Denominator B: " << endl;
	for( int i = 1; i <= hmm->getNumStates(); i++ ) {
		cout << denominatorB_t[i] << " ";
	}
	cout << endl;

//	return 1;

	// the sufficient statistics
	double ***xi_t;
	double **gamma_t;

//	scale = (double *) calloc ( (unsigned) ( T + 1 ), sizeof( double ) );

	T = 1000;

	xi_t = hmm->allocXi( T + 1, hmm->getNumStates() + 1 );
	gamma_t = (double **) calloc( (unsigned) ( T + 1 ), sizeof(double*) );
	for( int i = 1; i <= T; i++ ) {
		gamma_t[i] = (double *) calloc( (unsigned) ( hmm->getNumStates() + 1 ), sizeof(double) );
	}
	// Set the sufficient statistics to 0
	// 
	for( int t = 1; t <= T; t++ ) {
		for( int i = 1; i <= hmm->getNumStates(); i++ ) {
			gamma_t[t][i] = 0;
		}
	}
	for( int t = 1; t <= T - 1; t++ ) {
		for( int i = 1; i <= hmm->getNumStates(); i++ ) {
			for( int j = 1; j <= hmm->getNumStates(); j++ ) {
				xi_t[t][i][j] = 0;
			}
		}
	}

	int num_of_seqs_trained = 41;
	double logprob_trained = -6.010170E-01 * num_of_seqs_trained;

	int actual_obs_length = 0;
	for( int n = 1; n <= NS; n++ ) {
//		scanf( "%s", filename );
//		if( !strcmp( filename, "q" ) ) {
//			break;
//		}

//		fp = fopen( argv[argc - 1], "r" );

		// Read the sequence name
		// 
	//	cout << func_itos( n ) << endl;

		string ss;
		char tmp[1000] = "";
		char filename[100] = "/home/zkan/Desktop/accv-video/final/activity-group-labeling/4/";
		ss = func_itos( n );
		strcpy( tmp, ss.c_str() );
		strcat( filename, strcat( tmp, ".seq" ) );

//		fp = fopen( filename, "r" );
		fp = fopen( filename, "r" );
		if( fp == NULL ) {
			fprintf( stderr, "Error: File %s not found\n", filename );
			exit( 1 );
		}		
		hmm->readSequence( fp, &T, &O );
	//	cout << T << endl;
	
		ostringstream os_seq_name;
		os_seq_name << argv[optind + 1];

		// avoid the situation when the next sequence is shorter
		actual_obs_length = T;
		T = 1000;

		alpha = (double **) calloc( (unsigned) ( T + 1 ), sizeof(double*) );
		for( int i = 1; i <= T; i++ ) {
			alpha[i] = (double *) calloc( (unsigned) ( hmm->getNumStates() + 1 ), sizeof(double) );
		}
		beta = (double **) calloc( (unsigned) ( T + 1 ), sizeof(double*) );
		for( int i = 1; i <= T; i++ ) {
			beta[i] = (double *) calloc( (unsigned) ( hmm->getNumStates() + 1 ), sizeof(double) );
		}
		gamma = (double **) calloc( (unsigned) ( T + 1 ), sizeof(double*) );
		for( int i = 1; i <= T; i++ ) {
			gamma[i] = (double *) calloc( (unsigned) ( hmm->getNumStates() + 1 ), sizeof(double) );
		}

		hmm->baumWelchWithSS( T, O, alpha, beta, gamma, &niter, &logprobinit, &logprobfinal, 
					pi_t, numeratorA_t, denominatorA_t, numeratorB_t, denominatorB_t, n, os_seq_name, 
					forgetting_factor, actual_obs_length, num_of_seqs_trained, logprob_trained );

	//	cout << "Gamma: " << endl;
	//	for( int t = 1; t <= T; t++ ) {
	//		for( int s = 1; s <= hmm->getNumStates(); s++ ) {
	//			cout << gamma[t][s] << " ";
	//		}
	//		cout << endl;
	//	}
				
		fprintf( stdout, "Number of iterations: %d\n", niter );
		fprintf( stdout, "Log prob(O|init model): %E\n", logprobinit );
		fprintf( stdout, "Log prob(O|estimated model): %E\n", logprobfinal );

		hmm->printHMM();

		// -------------------------------------------------------------------------------------------------
		cout << endl << "## Checksum ##" << endl;
		double **tmp_A;		/* Temporal Current Transition probability */
		double **tmp_B;		/* Temporal Current Emission probability */
		double *tmp_pi;		/* Temporal Current Initial state distribution */
		tmp_A = (double **) calloc( (unsigned) ( hmm->getNumStates() + 1 ), sizeof(double*) );
		for( int i = 1; i <= hmm->getNumStates(); i++ ) {
			tmp_A[i] = (double *) calloc( (unsigned) ( hmm->getNumStates() + 1 ), sizeof(double) );
		}
		tmp_B = (double **) calloc( (unsigned) ( hmm->getNumStates() + 1 ), sizeof(double*) );
		for( int i = 1; i <= hmm->getNumStates(); i++ ) {
			tmp_B[i] = (double *) calloc( (unsigned) ( hmm->getNumSymbols() + 1 ), sizeof(double) );
		}
		tmp_pi = (double *) calloc( (unsigned) ( hmm->getNumStates() + 1 ), sizeof( double ) );

		tmp_A = hmm->getA();
		tmp_B = hmm->getB();
		tmp_pi = hmm->getPi();
//*
		float checkSum = 0;
		cout << "A: " << endl;
		for( int i = 1; i <= hmm->getNumStates(); i++ ) {
			cout << "sum of row " << i << ": ";
			for( int j = 1; j <= hmm->getNumStates(); j++ ) {
				checkSum += tmp_A[i][j];
			}
			cout << checkSum << endl;
			checkSum = 0;
		}
		cout << "B: " << endl;
		for( int i = 1; i <= hmm->getNumStates(); i++ ) {
			cout << "sum of row " << i << ": ";
			for( int j = 1; j <= hmm->getNumSymbols(); j++ ) {
				checkSum += tmp_B[i][j];
			}
			cout << checkSum << endl;
			checkSum = 0;
		}
		cout << "Pi: " << endl;
		cout << "sum of Pi: ";
		for( int i = 1; i <= hmm->getNumStates(); i++ ) {
			checkSum += tmp_pi[i];
		}
		cout << checkSum << endl;
		checkSum = 0;
//*/
		// -------------------------------------------------------------------------------------------------
/*
		checkSum = 0;
		cout << "Gamma:" << endl;
		for( int t = 1; t <= actual_obs_length; t++ ) {
			cout << "T= " << t << endl;
			for( int i = 1; i <= hmm->getNumStates(); i++ ) {
				checkSum += gamma_t[t][i];
			}
			cout << checkSum << endl;
			checkSum = 0;
		}
//*/
/*
		cout << "Xi:" << endl;
		for( int t = 1; t <= actual_obs_length - 1; t++ ) {
			for( int i = 1; i <= hmm->getNumStates(); i++ ) {
				for( int j = 1; j <= hmm->getNumSymbols(); j++ ) {
					cout << xi_t[t][i][j] << " ";
				}
				cout << endl;
			}
			cout << endl;
		}
		cout << endl;
//*/
//*
		// print the sufficient statistics
		cout << endl << "## SUFFICEINT STATISTICS ##" << endl;
		cout << "T= " << NS << endl;
		cout << "Pi: " << endl;
	//	pi_t = hmm->getPi();
		for( int i = 1; i <= hmm->getNumStates(); i++ ) {
			cout << pi_t[i] << " ";
		}
		cout << endl;
		cout << "Numerator A: " << endl;
		for( int i = 1; i <= hmm->getNumStates(); i++ ) {
			for( int j = 1; j <= hmm->getNumStates(); j++ ) {
				cout << numeratorA_t[i][j] << " ";
			}
			cout << endl;
		}
		cout << "Denominator A: " << endl;
		for( int i = 1; i <= hmm->getNumStates(); i++ ) {
			cout << denominatorA_t[i] << " ";
		}
		cout << endl;
		cout << "Numerator B: " << endl;
		for( int i = 1; i <= hmm->getNumStates(); i++ ) {
			for( int j = 1; j <= hmm->getNumSymbols(); j++ ) {
				cout << numeratorB_t[i][j] << " ";
			}
			cout << endl;
		}
		cout << "Denominator B: " << endl;
		for( int i = 1; i <= hmm->getNumStates(); i++ ) {
			cout << denominatorB_t[i] << " ";
		}
		cout << endl;

	//	hmm->printSuffStats();
//*/

//		free( alpha );
//		free( beta );
//		free( gamma );
//		free( gamma_t );
//		free( xi_t );

		fprintf( stdout, "Log prob(O|init model): %E\n", logprobinit );
		fprintf( stdout, "Log prob(O|estimated model): %E\n", logprobfinal );
	}

//	fclose( fp );

	return 1;
}



