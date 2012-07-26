/*
 *
 * HMMs with Sufficient Statistics by Kan Ouivirach ( zkan at cs dot ait dot ac dot th )
 *
 */

#include "hmm.h"

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
	double forgetting_factor = 1;

	// mean and variance of the log likelihood for setting the alarm threshold
	double sqrMean = 0, sum_of_sqrX = 0;
	double mean = 0, sd = 0;
	// Alarm threshold
	double Th = 1000000;

	// initialize the log likelihood
	logprobinit = 0;
	logprobfinal = 0;

	int c;
	extern char *optarg;
	extern int optind, opterr, optopt;
//*
	while( ( c = getopt( argc, argv, "f:" ) ) != EOF ) {
		switch( c ) {
//			case 'N':
//				sscanf( optarg, "%d", &NS );
//				break;
			case 'f':
				sscanf( optarg, "%lf", &forgetting_factor );
				break;
		}
	}

	if( ( ( argc - optind ) != 1 ) || ( argc == 1 ) ) {
		cout << "Usage error" << endl;
		cout << "Usage: xxx <model.hmm>" << endl;
		cout << "Option: " << endl;
		cout << "-f <forgetting factor>" << endl;
		exit( 1 );
	}
//*/	
	fp = fopen( argv[argc - 1], "r" );
	if( fp == NULL ) {
		fprintf( stderr, "Error: File %s not found\n", argv[1] );
		exit( 1 );
	}

	hmm = new HMM();
	hmm->readHMM( fp );
	hmm->printHMM();

//	double ***xi_t;
//	xi = allocXi( T, this->N );
//	scale = (double *) calloc ( (unsigned) ( T + 1 ), sizeof( double ) );

//	double **gamma_t;

	double *pi_t = (double *) calloc( (unsigned) ( hmm->getNumStates() + 1 ), sizeof( double ) );
	for( int i = 1; i <= hmm->getNumStates(); i++ ) {
		pi_t[i] = 0;
	}
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
	// Initialize numeratorA_t, denominatorA_t, numeratorB_t, denominatorB_t
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

	char filename[100];
	while( true ) {
		scanf( "%s", filename );
		if( !strcmp( filename, "q" ) ) {
			break;
		}
//		fp = fopen( argv[argc - 1], "r" );
		fp = fopen( filename, "r" );
		if( fp == NULL ) {
			fprintf( stderr, "Error: File %s not found\n", filename );
			exit( 1 );
		}
		
		int actual_obs_length;
		hmm->readSequence( fp, &actual_obs_length, &O );

		T = 200; // prevent the null pointer since each observation sequence has the different length
	//	T = actual_obs_length;

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

		for( int i = 1; i <= T; i++ ) {
			for( int j = 1; j <= hmm->getNumStates(); j++ ) {
				alpha[i][j] = 0;
				beta[i][j] = 0;
				gamma[i][j] = 0;
			}
		}
		
//		hmm->baumWelchWithSS( T, O, alpha, beta, gamma, &niter, &logprobinit, &logprobfinal, 
//					pi_t, numeratorA_t, denominatorA_t, numeratorB_t, denominatorB_t, NS, forgetting_factor, actual_obs_length );

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

		// calculate mean and variance of the log likelihood (per observation) before training
		double loglik_per_obs = logprobinit/actual_obs_length;

		sum_of_sqrX += pow( loglik_per_obs, 2 );
		mean += loglik_per_obs;
		sqrMean = pow( mean/NS, 2 );
		sd = sqrt( sum_of_sqrX/NS - sqrMean );

		cout << endl << "## Mean and SD of the log likelihood per observation ##" << endl;
		// ref: http://en.wikipedia.org/wiki/Standard_deviation
		cout << "Mean: " << mean/NS << endl;
		cout << "SD: " << sd << endl;

		cout << "Threshold: " << Th << endl; 
		fprintf( stdout, "Log prob per observation(O|init model): %lf\n", logprobinit/actual_obs_length );
		if( logprobinit/actual_obs_length < Th ) {
			cout << "ALARM!!!" << endl;
		}
		// Update the alarm threshold
		Th = mean/NS - sd;

		// print the sufficient statistics

//		cout << endl << "## SUFFICEINT STATISTICS ##" << endl;
		cout << "T= " << NS << endl;
/*		cout << "Forgetting factor= " << forgetting_factor << endl;
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
//*/
		free( alpha );
		free( beta );
		free( gamma );
	//	free( O );

		NS++;
	}

	fclose( fp );

	return 1;
}



